#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle

#include <boost/test/unit_test.hpp>

#include "conversion.hpp"
#include "softsusy.h"
#include "CMSSM_two_scale_model.hpp"
#include "test_CMSSM.hpp"
#include "stopwatch.hpp"

#define private public

#include "weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

BOOST_AUTO_TEST_CASE( test_rho_2 )
{
   Weinberg_angle wein;
   double r;

   r = 0.1;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 1.8;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 1.9;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 2.0;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);

   r = 2.1;
   BOOST_CHECK_CLOSE_FRACTION(rho2(r), wein.rho_2(r), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(fs, ss, input);
   fs.set_scale(173.34);
   ss.setMu(173.34);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);

   double outrho = 1.0, outsin = 0.48;
   const double alphaMsbar = ss.displayDataSet().displayAlpha(legacy::ALPHA);
   const double alphaDrbar = ss.qedSusythresh(alphaMsbar, scale);

   const double ss_delta_vb =
      ss.deltaVb(outrho, outsin, alphaDrbar, 0. /* ignored */);

   const double gY          = ss.displayGaugeCoupling(1) * sqrt(0.6);
   const double g2          = ss.displayGaugeCoupling(2);
   const double hmu         = ss.displayYukawaElement(YE, 2, 2);
   const double mselL       = ss.displayDrBarPars().me(1, 1);
   const double msmuL       = ss.displayDrBarPars().me(1, 2);
   const double msnue       = ss.displayDrBarPars().msnu(1);
   const double msnumu      = ss.displayDrBarPars().msnu(2);
   const Eigen::ArrayXd mneut = ToEigenArray(ss.displayDrBarPars().mnBpmz);
   const Eigen::MatrixXcd n   = ToEigenMatrix(ss.displayDrBarPars().nBpmz);
   const Eigen::ArrayXd mch   = ToEigenArray(ss.displayDrBarPars().mchBpmz);
   const Eigen::MatrixXcd u   = ToEigenMatrix(ss.displayDrBarPars().uBpmz);
   const Eigen::MatrixXcd v   = ToEigenMatrix(ss.displayDrBarPars().vBpmz);

   const auto MSe(fs.get_MSe());
   const auto ZE(fs.get_ZE());
   const auto MSv(fs.get_MSv());
   const auto ZV(fs.get_ZV());

   const double fs_gY          = fs.get_g1() * sqrt(0.6);
   const double fs_g2          = fs.get_g2();
   const double fs_hmu         = fs.get_Ye(1,1);
   double fs_mselL             = 0.;
   double fs_msmuL             = 0.;
   double fs_msnue             = 0.;
   double fs_msnumu            = 0.;
   const Eigen::ArrayXd fs_mneut = fs.get_MChi();
   const Eigen::MatrixXcd fs_n   = fs.get_ZN();
   const Eigen::ArrayXd fs_mch   = fs.get_MCha();
   const Eigen::MatrixXcd fs_u   = fs.get_UM();
   const Eigen::MatrixXcd fs_v   = fs.get_UP();

   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      fs_mselL += AbsSqr(ZE(i,0))*MSe(i);
      fs_msmuL += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      fs_msnue  += AbsSqr(ZV(i,0))*MSv(i);
      fs_msnumu += AbsSqr(ZV(i,1))*MSv(i);
   }

   BOOST_CHECK_CLOSE_FRACTION(fs_mselL , mselL , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msmuL , msmuL , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msnue , msnue , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msnumu, msnumu, 1.0e-10);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.msel_drbar = mselL;
   data.msmul_drbar = msmuL;
   data.msve_drbar = msnue;
   data.msvm_drbar = msnumu;
   data.mn_drbar = mneut;
   data.mc_drbar = mch;
   data.zn = n;
   data.um = u;
   data.up = v;
   data.gY = gY;
   data.g2 = g2;

   // test with SoftSusy parameters
   double fs_delta_vb =
      Weinberg_angle::calculate_delta_vb(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_vb, fs_delta_vb, 1.0e-10);

   // test with FlexibleSUSY CMSSM parameters

   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.msel_drbar = fs_mselL;
   data.msmul_drbar = fs_msmuL;
   data.msve_drbar = fs_msnue;
   data.msvm_drbar = fs_msnumu;
   data.mn_drbar = fs_mneut;
   data.mc_drbar = fs_mch;
   data.zn = fs_n;
   data.um = fs_u;
   data.up = fs_v;
   data.gY = fs_gY;
   data.g2 = fs_g2;

   fs_delta_vb =
      Weinberg_angle::calculate_delta_vb(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_vb, fs_delta_vb, 3.0e-9);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(fs, ss, input);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaMsbar = ss.displayDataSet().displayAlpha(legacy::ALPHA);
   const double alphaDrbar = ss.qedSusythresh(alphaMsbar, scale);
   const double pizztMZ    = ss.piZZT(mz_pole, scale, true);
   const double piwwt0     = ss.piWWT(0., scale, true);
   const double gfermi     = Electroweak_constants::gfermi;
   softsusy::GMU = gfermi;

   const double ss_delta_r =
      ss.dR(outrho, outsin, alphaDrbar, pizztMZ, piwwt0);

   const double gY          = ss.displayGaugeCoupling(1) * sqrt(0.6);
   const double g2          = ss.displayGaugeCoupling(2);
   const double hmu         = ss.displayYukawaElement(YE, 2, 2);
   const double mselL       = ss.displayDrBarPars().me(1, 1);
   const double msmuL       = ss.displayDrBarPars().me(1, 2);
   const double msnue       = ss.displayDrBarPars().msnu(1);
   const double msnumu      = ss.displayDrBarPars().msnu(2);
   const Eigen::ArrayXd mneut = ToEigenArray(ss.displayDrBarPars().mnBpmz);
   const Eigen::MatrixXcd n   = ToEigenMatrix(ss.displayDrBarPars().nBpmz);
   const Eigen::ArrayXd mch   = ToEigenArray(ss.displayDrBarPars().mchBpmz);
   const Eigen::MatrixXcd u   = ToEigenMatrix(ss.displayDrBarPars().uBpmz);
   const Eigen::MatrixXcd v   = ToEigenMatrix(ss.displayDrBarPars().vBpmz);
   const double mt          = ss.displayDataSet().displayPoleMt();
   const double g3          = ss.displayGaugeCoupling(3);
   const double tanBeta     = ss.displayTanb();
   const double mh          = ss.displayDrBarPars().mh0(1);
   const double alpha       = ss.h1s2Mix();

   const auto MSe(fs.get_MSe());
   const auto ZE(fs.get_ZE());
   const auto MSv(fs.get_MSv());
   const auto ZV(fs.get_ZV());

   const double fs_gY          = fs.get_g1() * sqrt(0.6);
   const double fs_g2          = fs.get_g2();
   const double fs_hmu         = fs.get_Ye(1,1);
   double fs_mselL             = 0.;
   double fs_msmuL             = 0.;
   double fs_msnue             = 0.;
   double fs_msnumu            = 0.;
   const Eigen::ArrayXd fs_mneut = fs.get_MChi();
   const Eigen::MatrixXcd fs_n   = fs.get_ZN();
   const Eigen::ArrayXd fs_mch   = fs.get_MCha();
   const Eigen::MatrixXcd fs_u   = fs.get_UM();
   const Eigen::MatrixXcd fs_v   = fs.get_UP();
   const double fs_mt          = ss.displayDataSet().displayPoleMt();
   const double fs_g3          = fs.get_g3();
   const double fs_tanBeta     = fs.get_vu() / fs.get_vd();
   const double fs_mh          = fs.get_Mhh(0);
   const double fs_alpha       = fs.get_ZH(0,1);

   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      fs_mselL += AbsSqr(ZE(i,0))*MSe(i);
      fs_msmuL += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      fs_msnue  += AbsSqr(ZV(i,0))*MSv(i);
      fs_msnumu += AbsSqr(ZV(i,1))*MSv(i);
   }

   BOOST_CHECK_CLOSE_FRACTION(fs_mselL , mselL , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msmuL , msmuL , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msnue , msnue , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_msnumu, msnumu, 1.0e-10);

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, mt, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_g3, g3, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_tanBeta, tanBeta, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_mh, mh, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_alpha, alpha, 1.0e-10);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = gfermi;
   data.self_energy_z_at_mz = pizztMZ;
   data.self_energy_w_at_0 = piwwt0;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = mt;
   data.mh_drbar = mh;
   data.hmix_12 = alpha;
   data.msel_drbar = mselL;
   data.msmul_drbar = msmuL;
   data.msve_drbar = msnue;
   data.msvm_drbar = msnumu;
   data.mn_drbar = mneut;
   data.mc_drbar = mch;
   data.zn = n;
   data.um = u;
   data.up = v;
   data.gY = gY;
   data.g2 = g2;
   data.g3 = g3;
   data.tan_beta = tanBeta;

   // test with SoftSusy parameters
   double fs_delta_r =
      Weinberg_angle::calculate_delta_r(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_r, fs_delta_r, 1.0e-10);

   // test with FlexibleSUSY CMSSM parameters

   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = gfermi;
   data.self_energy_z_at_mz = pizztMZ;
   data.self_energy_w_at_0 = piwwt0;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = fs_mt;
   data.mh_drbar = fs_mh;
   data.hmix_12 = fs_alpha;
   data.msel_drbar = fs_mselL;
   data.msmul_drbar = fs_msmuL;
   data.msve_drbar = fs_msnue;
   data.msvm_drbar = fs_msnumu;
   data.mn_drbar = fs_mneut;
   data.mc_drbar = fs_mch;
   data.zn = fs_n;
   data.um = fs_u;
   data.up = fs_v;
   data.gY = fs_gY;
   data.g2 = fs_g2;
   data.g3 = fs_g3;
   data.tan_beta = fs_tanBeta;

   fs_delta_r =
      Weinberg_angle::calculate_delta_r(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_r, fs_delta_r, 4.0e-9);
}

BOOST_AUTO_TEST_CASE( test_delta_rho )
{
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(fs, ss, input);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   double outrho = 1.0, outsin = 0.48;
   const double alphaMsbar = ss.displayDataSet().displayAlpha(legacy::ALPHA);
   const double alphaDrbar = ss.qedSusythresh(alphaMsbar, scale);
   const double pizztMZ    = ss.piZZT(mz_pole, scale, true);
   const double piwwtMW    = ss.piWWT(mw_pole, scale, true);
   const double gfermi     = Electroweak_constants::gfermi;
   softsusy::GMU = gfermi;

   const double ss_delta_rho =
      ss.dRho(outrho, outsin, alphaDrbar, pizztMZ, piwwtMW);

   const double mt          = ss.displayDataSet().displayPoleMt();
   const double g3          = ss.displayGaugeCoupling(3);
   const double tanBeta     = ss.displayTanb();
   const double mh          = ss.displayDrBarPars().mh0(1);
   const double alpha       = ss.h1s2Mix();

   const double fs_mt          = ss.displayDataSet().displayPoleMt();
   const double fs_g3          = fs.get_g3();
   const double fs_tanBeta     = fs.get_vu() / fs.get_vd();
   const double fs_mh          = fs.get_Mhh(0);
   const double fs_alpha       = fs.get_ZH(0,1);

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, mt, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_g3, g3, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_tanBeta, tanBeta, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_mh, mh, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_alpha, alpha, 1.0e-10);

   Weinberg_angle::Data data;
   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = gfermi;
   data.self_energy_z_at_mz = pizztMZ;
   data.self_energy_w_at_mw = piwwtMW;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = mt;
   data.mh_drbar = mh;
   data.hmix_12 = alpha;
   data.g3 = g3;
   data.tan_beta = tanBeta;

   // test with SoftSusy parameters
   double fs_delta_rho =
      Weinberg_angle::calculate_delta_rho(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_rho, fs_delta_rho, 1.0e-10);

   // test with FlexibleSUSY CMSSM parameters

   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = gfermi;
   data.self_energy_z_at_mz = pizztMZ;
   data.self_energy_w_at_mw = piwwtMW;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = fs_mt;
   data.mh_drbar = fs_mh;
   data.hmix_12 = fs_alpha;
   data.g3 = fs_g3;
   data.tan_beta = fs_tanBeta;

   fs_delta_rho =
      Weinberg_angle::calculate_delta_rho(outrho, outsin, data);

   BOOST_CHECK_CLOSE_FRACTION(ss_delta_rho, fs_delta_rho, 1.0e-10);
}

void setup_data(const CMSSM_input_parameters& input,
                CMSSM<Two_scale>& fs,
                MssmSoftsusy& ss,
                Weinberg_angle::Data& data)
{
   setup_CMSSM_const_non_3rd_gen(fs, ss, input);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);
   BOOST_CHECK_EQUAL(fs.get_scale(), scale);

   const double gfermi = Electroweak_constants::gfermi;
   const double pizztMZ = ss.piZZT(mz_pole, scale, true);
   const double piwwt0  = ss.piWWT(0., scale, true);
   const double piwwtMW = ss.piWWT(mw_pole, scale, true);
   const double alphaMsbar = ss.displayDataSet().displayAlpha(legacy::ALPHA);
   const double alphaDrbar = ss.qedSusythresh(alphaMsbar, scale);
   const double gY         = fs.get_g1() * sqrt(0.6);
   const double g2         = fs.get_g2();
   const double hmu        = fs.get_Ye(1,1);
   const double g3         = fs.get_g3();
   const double mt         = ss.displayDataSet().displayPoleMt();
   const double mh         = fs.get_Mhh(0);
   const double alpha      = fs.get_ZH(0,1);
   const double tanBeta    = fs.get_vu() / fs.get_vd();
   double mselL            = 0.;
   double msmuL            = 0.;
   double msnue            = 0.;
   double msnumu           = 0.;
   const Eigen::ArrayXd mneut = fs.get_MChi();
   const Eigen::MatrixXcd n   = fs.get_ZN();
   const Eigen::ArrayXd mch   = fs.get_MCha();
   const Eigen::MatrixXcd u   = fs.get_UM();
   const Eigen::MatrixXcd v   = fs.get_UP();

   const double ss_pizztMZ = ss.piZZT(mz_pole, scale, false);
   const double ss_piwwt0  = ss.piWWT(0., scale, false);
   const double ss_piwwtMW = ss.piWWT(mw_pole, scale, false);
   const double fs_pizztMZ = Re(fs.self_energy_VZ_1loop(mz_pole));
   const double fs_piwwt0  = Re(fs.self_energy_VWm_1loop(0));
   const double fs_piwwtMW = Re(fs.self_energy_VWm_1loop(mw_pole));

   BOOST_CHECK_CLOSE_FRACTION(ss_pizztMZ, fs_pizztMZ, 5.0e-08);
   BOOST_CHECK_CLOSE_FRACTION(ss_piwwtMW, fs_piwwtMW, 3.0e-05);
   BOOST_CHECK_CLOSE_FRACTION(ss_piwwt0 , fs_piwwt0 , 5.0e-04);

   const auto MSe(fs.get_MSe());
   const auto ZE(fs.get_ZE());
   const auto MSv(fs.get_MSv());
   const auto ZV(fs.get_ZV());

   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      mselL += AbsSqr(ZE(i,0))*MSe(i);
      msmuL += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      msnue  += AbsSqr(ZV(i,0))*MSv(i);
      msnumu += AbsSqr(ZV(i,1))*MSv(i);
   }

   softsusy::GMU = gfermi;
   softsusy::PRINTOUT = 10;

   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = gfermi;
   data.self_energy_z_at_mz = pizztMZ;
   data.self_energy_w_at_mw = piwwtMW;
   data.self_energy_w_at_0 = piwwt0;
   data.mw_pole = mw_pole;
   data.mz_pole = mz_pole;
   data.mt_pole = mt;
   data.mh_drbar = mh;
   data.hmix_12 = alpha;
   data.msel_drbar = mselL;
   data.msmul_drbar = msmuL;
   data.msve_drbar = msnue;
   data.msvm_drbar = msnumu;
   data.mn_drbar = mneut;
   data.mc_drbar = mch;
   data.zn = n;
   data.um = u;
   data.up = v;
   data.gY = gY;
   data.g2 = g2;
   data.g3 = g3;
   data.tan_beta = tanBeta;
}

BOOST_AUTO_TEST_CASE( test_rho_sinTheta )
{
   Stopwatch stopwatch;
   Weinberg_angle::Data data;
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_data(input, fs, ss, data);

   const double tol = 1.0e-8;
   const int maxTries = 20;

   BOOST_TEST_MESSAGE("running MssmSoftsusy::rhohat() ...");
   const double rho_start = 1.0, sin_start = 0.48;
   double outrho = rho_start, outsin = sin_start;
   stopwatch.start();
   ss.rhohat(outrho, outsin,
             data.alpha_em_drbar,
             data.self_energy_z_at_mz,
             data.self_energy_w_at_0,
             data.self_energy_w_at_mw,
             tol, maxTries);
   stopwatch.stop();
   const double ss_time = stopwatch.get_time_in_seconds();
   BOOST_TEST_MESSAGE("MssmSoftsusy::rhohat() finished!");
   BOOST_TEST_MESSAGE("MssmSoftsusy::rhohat() takes " << ss_time << " seconds");

   Weinberg_angle weinberg;
   weinberg.set_number_of_iterations(maxTries);
   weinberg.set_precision_goal(tol);
   weinberg.set_data(data);

   BOOST_TEST_MESSAGE("running Weinberg_angle::calculate() ...");

   stopwatch.start();
   const int error = weinberg.calculate(rho_start, sin_start);
   stopwatch.stop();
   const double fs_time = stopwatch.get_time_in_seconds();

   BOOST_REQUIRE(error == 0);
   const double fs_sintheta = weinberg.get_sin_theta();

   BOOST_TEST_MESSAGE("running Weinberg_angle::calculate() finished!");
   BOOST_TEST_MESSAGE("Weinberg_angle::calculate() takes " << fs_time
                 << " seconds");

   const double fs_rhohat = weinberg.get_rho_hat();

   BOOST_CHECK_CLOSE_FRACTION(outsin, fs_sintheta, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(outrho, fs_rhohat  , 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_self_energy_top_correction )
{
   Weinberg_angle::Data data;
   CMSSM<Two_scale> fs;
   MssmSoftsusy ss;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_data(input, fs, ss, data);

   const double scale = ss.displayMu();
   const double mw_pole = ss.displayMw();
   const double mz_pole = ss.displayMz();
   const double mt_drbar = fs.get_MFu(2);
   const double mb_drbar = fs.get_MFd(2);

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);
   BOOST_CHECK_EQUAL(fs.get_scale(), mz_pole);
   BOOST_CHECK_CLOSE_FRACTION(mt_drbar, ss.displayDrBarPars().mt, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(mb_drbar, ss.displayDrBarPars().mb, 1.0e-10);

   // reference values
   const double pizztMZ = ss.piZZT(mz_pole, scale, true);
   const double piwwtMW = ss.piWWT(mw_pole, scale, true);
   const double piwwt0  = ss.piWWT(0.     , scale, true);

   const double ss_pizztMZ = ss.piZZT(mz_pole, scale, false);
   const double ss_piwwtMW = ss.piWWT(mw_pole, scale, false);
   const double ss_piwwt0  = ss.piWWT(0.     , scale, false);
   const double fs_pizztMZ = Re(fs.self_energy_VZ_1loop(mz_pole));
   const double fs_piwwtMW = Re(fs.self_energy_VWm_1loop(mw_pole));
   const double fs_piwwt0  = Re(fs.self_energy_VWm_1loop(0));

   BOOST_CHECK_CLOSE_FRACTION(ss_pizztMZ, fs_pizztMZ, 5.0e-08);
   BOOST_CHECK_CLOSE_FRACTION(ss_piwwtMW, fs_piwwtMW, 3.0e-05);
   BOOST_CHECK_CLOSE_FRACTION(ss_piwwt0 , fs_piwwt0 , 5.0e-04);

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale = scale;
   se_data.mt_pole = data.mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY = data.gY;
   se_data.g2 = data.g2;

   const double fs_pizztMZ_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_z(fs_pizztMZ, mz_pole, se_data);

   const double fs_piwwtMW_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(fs_piwwtMW, mw_pole, se_data);

   const double fs_piwwt0_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(fs_piwwt0, 0., se_data);

   BOOST_CHECK_CLOSE_FRACTION(fs_pizztMZ_corrected, pizztMZ, 5.0e-08);
   BOOST_CHECK_CLOSE_FRACTION(fs_piwwtMW_corrected, piwwtMW, 3.0e-05);
   BOOST_CHECK_CLOSE_FRACTION(fs_piwwt0_corrected , piwwt0 , 4.0e-04);
}
