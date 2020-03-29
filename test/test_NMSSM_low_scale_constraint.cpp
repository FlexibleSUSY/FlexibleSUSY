
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "test_NMSSM.hpp"
#include "NMSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_low_scale_constraint.hpp"
#include "nmssmsoftsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"

softsusy::QedQcd convert(const softsusy::QedQcd_legacy& ql)
{
   softsusy::QedQcd qn;

   qn.setAlphas(flexiblesusy::ToEigenArray(ql.displayAlphas()));
   qn.setMasses(flexiblesusy::ToEigenArray(ql.displayMass()));
   qn.set_input(ql.display_input());
   qn.setPoleMb(ql.displayPoleMb());
   qn.setCKM(ql.displayCKM());
   qn.setPMNS(ql.displayPMNS());
   qn.set_number_of_parameters(ql.howMany());
   qn.set_scale(ql.displayMu());
   qn.set_loops(ql.displayLoops());
   qn.set_thresholds(ql.displayThresholds());

   return qn;
}

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   QedQcd_legacy qedqcd;
   setup_NMSSM_const(m, s, input);
   s.setData(qedqcd);

   m.calculate_DRbar_masses();
   s.calcDrBarPars();

   NMSSM_low_scale_constraint<Two_scale> constraint(&m, convert(qedqcd));

   const double alpha_em = qedqcd.displayAlpha(legacy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(legacy::ALPHAS);
   const double scale = m.get_scale();

   const double delta_alpha_em_fs = constraint.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_fs  = constraint.calculate_delta_alpha_s(alpha_s);

   const double delta_alpha_em_ss = 1.0 - alpha_em / s.qedSusythresh(alpha_em, scale);
   const double delta_alpha_s_ss  = 1.0 - alpha_s  / s.qcdSusythresh(alpha_s , scale);

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_fs, delta_alpha_em_ss, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_fs , delta_alpha_s_ss , 1.0e-12);
}

BOOST_AUTO_TEST_CASE( test_low_energy_constraint )
{
   NMSSM<Two_scale> m;
   m.set_thresholds(2);
   NmssmSoftsusy s;
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   QedQcd_legacy qedqcd;
   qedqcd.setPoleMt(175.);       // non-default
   qedqcd.setMass(legacy::mBottom, 4.3); // non-default
   setup_NMSSM_const(m, s, input);
   s.setData(qedqcd);

   m.calculate_DRbar_masses();
   s.calcDrBarPars();

   NMSSM_low_scale_constraint<Two_scale> constraint(&m, convert(qedqcd));

   const double TanBeta = input.TanBeta;
   const double g1 = m.get_g1();
   const double g2 = m.get_g2();

   const double ss_mt = s.calcRunningMt();
   const double ss_mb = s.calcRunningMb();
   const double ss_me = s.calcRunningMtau();
   const double MZ    = s.displayMz();
   const double pizzt = s.piZZT(MZ, s.displayMu());
   const double ss_MZ = Sqrt(Sqr(MZ) + pizzt);
   const double ss_new_vev = s.getVev();

   const double fs_mt = m.calculate_MFu_DRbar(qedqcd.displayPoleMt(), 2);
   const double fs_mb = m.calculate_MFd_DRbar(qedqcd.displayMass(legacy::mBottom), 2);
   const double fs_me = m.calculate_MFe_DRbar(qedqcd.displayMass(legacy::mTau), 2);
   const double fs_MZ = m.calculate_MVZ_DRbar(Electroweak_constants::MZ);
   const double fs_old_vd = m.get_vd();
   const double fs_old_vu = m.get_vu();
   // const double fs_old_vev = Sqrt(Sqr(fs_old_vu) + Sqr(fs_old_vd));
   const double fs_new_vd = (2*fs_MZ)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vu = (2*fs_MZ*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vev = Sqrt(Sqr(fs_new_vu) + Sqr(fs_new_vd));

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, ss_mt, 9.5e-05);
   BOOST_CHECK_CLOSE_FRACTION(fs_mb, ss_mb, 2.6e-14);
   BOOST_CHECK_CLOSE_FRACTION(fs_me, ss_me, 2.1e-04);
   BOOST_CHECK_CLOSE_FRACTION(fs_MZ, ss_MZ, 5.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_new_vev, ss_new_vev, 5.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_old_vu / fs_old_vd, s.displayTanb(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_new_vu / fs_new_vd, s.displayTanb(), 1.0e-10);

   // apply constraints
   constraint.apply();
   s.sparticleThresholdCorrections(input.TanBeta);

   BOOST_CHECK_CLOSE_FRACTION(m.get_g1(), s.displayGaugeCoupling(1), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g2(), s.displayGaugeCoupling(2), 0.007);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g3(), s.displayGaugeCoupling(3), 1.0e-12);

   // test off-diagonal elements
   BOOST_TEST_MESSAGE("testing off-diagonal yukawa elements");
   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         if (i == k)
            continue;
         BOOST_TEST_MESSAGE("testing yukawa elements " << i << ", " << k);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(i-1,k-1), s.displayYukawaMatrix(YU)(i,k), 1.0e-5);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(i-1,k-1), s.displayYukawaMatrix(YD)(i,k), 1.0e-5);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(i-1,k-1), s.displayYukawaMatrix(YE)(i,k), 1.0e-5);
      }
   }

   // The following Yukawa couplings differ a lot from Softsusy,
   // because Ben uses the new vev (= the value of the vev after
   // sparticleThresholdCorrections() was called) to calculate the
   // Yukawa couplings.  We use the old vev (= combination of vu, vd
   // from the last run) to calculate the Yukawa couplings.

   BOOST_TEST_MESSAGE("testing diagonal yukawa elements");
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(0,0), s.displayYukawaMatrix(YU)(1,1), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(0,0), s.displayYukawaMatrix(YD)(1,1), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(0,0), s.displayYukawaMatrix(YE)(1,1), 0.005);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(1,1), s.displayYukawaMatrix(YU)(2,2), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(1,1), s.displayYukawaMatrix(YD)(2,2), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(1,1), s.displayYukawaMatrix(YE)(2,2), 0.005);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(2,2), s.displayYukawaMatrix(YU)(3,3), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(2,2), s.displayYukawaMatrix(YD)(3,3), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(2,2), s.displayYukawaMatrix(YE)(3,3), 0.005);

   BOOST_TEST_MESSAGE("testing running VEV");
   const double running_vev = Sqrt(Sqr(m.get_vu()) +  Sqr(m.get_vd()));
   BOOST_CHECK_CLOSE_FRACTION(running_vev, s.displayHvev(), 1.0e-9);
}
