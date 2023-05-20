
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_weinberg_angle_meta

#include <boost/test/unit_test.hpp>

#include "CMSSM_mass_eigenstates.hpp"
#include "test_CMSSMNoFV.hpp"
#include "stopwatch.hpp"

#define private public

#include "weinberg_angle.hpp"
#include "CMSSM_weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

void setup_CMSSM_const_non_3rd_gen(CMSSM_mass_eigenstates& model,
                                   const CMSSM_input_parameters& input)
{
   setup_CMSSM_const(model, input);
   model.set_thresholds(2);

   const double ymu = 0.001;

   Eigen::Matrix<double,3,3> Ye = model.get_Ye();
   Ye(1,1) = ymu;
   model.set_Ye(Ye);
}

void setup_data(const CMSSM_mass_eigenstates& model,
                Weinberg_angle::Data& data)
{
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;
   // use mt_drbar instead of mt_pole such that tests work up to high precision
   const double mt_pole = 165.0;
   const double scale   = model.get_scale();
   const double gY      = model.get_g1() * Sqrt(0.6);
   const double g2      = model.get_g2();
   const double e_drbar = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));

   BOOST_CHECK_EQUAL(scale, mz_pole);

   data.fermi_contant  = Electroweak_constants::gfermi;
   data.mw_pole        = mw_pole;
   data.mz_pole        = mz_pole;
   data.mt_pole        = mt_pole;
   data.scale          = scale;
   data.gY             = gY;
   data.g2             = g2;
   data.g3             = model.get_g3();
   data.alpha_em_drbar = Sqr(e_drbar) / (4.0 * Pi);
   data.tan_beta       = model.get_vu() / model.get_vd();
   data.hmix_12        = model.get_ZH(0,1);
   data.mh_drbar       = model.get_Mhh(0);
   data.mn_drbar       = model.get_MChi();
   data.mc_drbar       = model.get_MCha();
   data.zn             = model.get_ZN();
   data.um             = model.get_UM();
   data.up             = model.get_UP();

   double msel  = 0.;
   double msmul = 0.;
   double msve  = 0.;
   double msvm  = 0.;
   const auto MSe = model.get_MSe();
   const auto ZE  = model.get_ZE();
   const auto MSv = model.get_MSv();
   const auto ZV  = model.get_ZV();

   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      msel  += AbsSqr(ZE(i,0))*MSe(i);
      msmul += AbsSqr(ZE(i,1))*MSe(i);
   }

   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      msve += AbsSqr(ZV(i,0))*MSv(i);
      msvm += AbsSqr(ZV(i,1))*MSv(i);
   }

   data.msel_drbar  = msel;
   data.msmul_drbar = msmul;
   data.msve_drbar  = msve;
   data.msvm_drbar  = msvm;

   const double pizztMZ = Re(model.self_energy_VZ_1loop(mz_pole));
   const double piwwtMW = Re(model.self_energy_VWm_1loop(mw_pole));
   const double piwwt0  = Re(model.self_energy_VWm_1loop(0.));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = model.get_MFu(2);
   se_data.mb_drbar = model.get_MFd(2);
   se_data.gY       = gY;
   se_data.g2       = g2;

   const double pizztMZ_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole, se_data);

   const double piwwtMW_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwtMW, mw_pole, se_data);

   const double piwwt0_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0., se_data);

   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
}

BOOST_AUTO_TEST_CASE( test_rho_2 )
{
   double r;

   r = 0.1;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                              CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 1.8;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                              CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 1.9;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                              CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 2.0;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                              CMSSM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 2.1;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                              CMSSM_weinberg_angle::rho_2(r), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   CMSSM_mass_eigenstates model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   const double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   const double delta_vb_1 = Weinberg_angle::calculate_delta_vb(outrho, outsin, data);

   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   CMSSM_weinberg_angle wein(&model, sm_parameters);
   const double delta_vb_2 = wein.calculate_delta_vb(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_vb_1, delta_vb_2, 6.0e-9);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   CMSSM_mass_eigenstates model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   const double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   const double delta_r_1 = Weinberg_angle::calculate_delta_r(outrho, outsin, data);

   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   CMSSM_weinberg_angle wein(&model, sm_parameters);
   // initialize self-energies
   wein.pizzt_MZ = wein.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein.piwwt_MW = wein.calculate_self_energy_VWm(Electroweak_constants::MW);
   wein.piwwt_0  = wein.calculate_self_energy_VWm(0.);
   const double delta_r_2 = wein.calculate_delta_r_hat(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_r_1, delta_r_2, 4.0e-6);
}

BOOST_AUTO_TEST_CASE( test_sin_theta )
{
   CMSSM_mass_eigenstates model;
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const_non_3rd_gen(model, input);

   Stopwatch stopwatch;
   const double tol = 1.0e-10;
   const int maxTries = 20;
   const double rho_start = 1.0, sin_start = 0.48;

   stopwatch.start();
   Weinberg_angle::Data data;
   setup_data(model, data);
   Weinberg_angle wein;
   wein.set_data(data);
   wein.set_number_of_iterations(maxTries);
   wein.set_precision_goal(tol);
   const int error = wein.calculate(rho_start, sin_start);
   const double sin_theta_1 = wein.get_sin_theta();
   stopwatch.stop();
   const double time_1 = stopwatch.get_time_in_seconds();
   BOOST_REQUIRE(error == 0);
   BOOST_TEST_MESSAGE("calculation with Weinberg_angle took "
                      << time_1 << " seconds" << '\n');

   stopwatch.start();
   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   CMSSM_weinberg_angle pwein(&model, sm_parameters);
   pwein.set_number_of_iterations(maxTries);
   pwein.set_precision_goal(tol);
   double sin_theta_2;
   BOOST_REQUIRE_NO_THROW(sin_theta_2 = pwein.calculate(sin_start).first);
   stopwatch.stop();
   const double time_2 = stopwatch.get_time_in_seconds();
   BOOST_TEST_MESSAGE("calculation with CMSSM_weinberg_angle took "
                      << time_2 << " seconds");

   BOOST_CHECK_CLOSE_FRACTION(sin_theta_1, sin_theta_2, 2.0e-6);
}
