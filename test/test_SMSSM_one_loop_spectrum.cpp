
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void ensure_tree_level_ewsb(SMSSM<Two_scale>& m, NmssmSoftsusy& s,
                            const SMSSM_input_parameters& input)
{
   const double precision = 1.0e-5;

   // initial guess
   m.set_Kappa(0.1);
   m.set_vS(5000.);
   m.set_ms2(-Sqr(input.m0));
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   s.setKappa(m.get_Kappa());
   s.setSvev(m.get_vS());
   s.setMsSquared(m.get_ms2());
   s.setMh1Squared(m.get_mHd2());
   s.setMh2Squared(m.get_mHu2());

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_CHECK_EQUAL(error, 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3(), precision);

   softsusy::Z3 = false;
   s.rewsbTreeLevel(1);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Mu() , s.displaySusyMu()   , precision * 15.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_BMu(), s.displayM3Squared(), precision * 10.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_LL1(), s.displayXiS()      , precision * 5.);

   m.set_Mu(s.displaySusyMu());
   m.set_BMu(s.displayM3Squared());
   m.set_LL1(s.displayXiS());
}

void ensure_one_loop_ewsb(SMSSM<Two_scale>& m, NmssmSoftsusy& s)
{
   const double precision = 1.0e-5;

   s.calcDrBarPars();
   m.calculate_DRbar_masses();

   const double mt = s.displayDrBarPars().mt;
   const int signMu = 1;

   m.set_ewsb_iteration_precision(precision);
   m.solve_ewsb_one_loop();

   softsusy::Z3 = false;
   softsusy::numRewsbLoops = 1;
   s.rewsb(signMu, mt);

   const double mu_ss  = s.displaySusyMu();
   const double bmu_ss = s.displayM3Squared();
   const double xis_ss = s.displayXiS();

   const double mu_fs  = m.get_Mu();
   const double bmu_fs = m.get_BMu();
   const double xis_fs = m.get_LL1();

   BOOST_CHECK_CLOSE_FRACTION(mu_ss , mu_fs , 1.6e-8);
   BOOST_CHECK_CLOSE_FRACTION(bmu_ss, bmu_fs, 6.4e-8);
   BOOST_CHECK_CLOSE_FRACTION(xis_ss, xis_fs, 2.1e-9);

   m.set_Mu(s.displaySusyMu());
   m.set_BMu(s.displayM3Squared());
   m.set_LL1(s.displayXiS());
}

BOOST_AUTO_TEST_CASE( test_SMSSM_pole_masses )
{
   SMSSM_input_parameters input;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaSInput = 100.;
   input.m0       = 540.;
   input.Azero    = -350.;
   input.MSInput  = 290.;
   input.BMSInput = 400.;
   input.L1Input  = 300.;
   SMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   setup_SMSSM_const(m, s, input);

   softsusy::Z3 = false;
   ensure_tree_level_ewsb(m, s, input);

   softsusy::numHiggsMassLoops = 1;
   s.physical(1);
   m.set_pole_mass_loop_order(1);
   m.calculate_DRbar_masses();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   BOOST_CHECK_EQUAL(s.displayMu(), m.get_scale());

   // neutralinos
   const DoubleVector MChi(ToDoubleVector(m.get_physical().MChi));
   const DoubleVector mneut = s.displayPhys().mneut.apply(fabs);
   BOOST_CHECK_CLOSE(mneut(1), MChi(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(2), MChi(2), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(3), MChi(3), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(4), MChi(4), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(5), MChi(5), 1.0e-12);

   // charginos
   const DoubleVector MCha(ToDoubleVector(m.get_physical().MCha));
   const DoubleVector mch = s.displayPhys().mch.apply(fabs);
   BOOST_CHECK_CLOSE(mch(1), MCha(1), 5.0e-9);
   BOOST_CHECK_CLOSE(mch(2), MCha(2), 2.0e-9);

   // photon, gluon mass
   const double vp = m.get_physical().MVP;
   const double vg = m.get_physical().MVG;
   BOOST_CHECK_EQUAL(vp, 0.0);
   BOOST_CHECK_EQUAL(vg, 0.0);

   // gluinos
   const double MGlu = m.get_physical().MGlu;
   const double mGluino = s.displayPhys().mGluino;
   BOOST_CHECK_CLOSE(MGlu, mGluino, 4.0e-5);

   // down-type squarks
   const DoubleVector Sd(ToDoubleVector(m.get_physical().MSd));
   const DoubleVector md(s.displayPhys().md.flatten().sort());
   BOOST_CHECK_CLOSE(Sd(1), md(1), 4.0e-8);
   BOOST_CHECK_CLOSE(Sd(2), md(2), 4.0e-8);
   BOOST_CHECK_CLOSE(Sd(3), md(3), 4.0e-8);
   BOOST_CHECK_CLOSE(Sd(4), md(4), 4.0e-8);
   BOOST_CHECK_CLOSE(Sd(5), md(5), 4.0e-8);
   BOOST_CHECK_CLOSE(Sd(6), md(6), 4.0e-8);

   // up-type squarks
   const DoubleVector Su(ToDoubleVector(m.get_physical().MSu));
   const DoubleVector mu(s.displayPhys().mu.flatten().sort());
   BOOST_CHECK_CLOSE(Su(1), mu(1), 4.0e-8);
   BOOST_CHECK_CLOSE(Su(2), mu(2), 4.0e-8);
   BOOST_CHECK_CLOSE(Su(3), mu(3), 4.0e-8);
   BOOST_CHECK_CLOSE(Su(4), mu(4), 4.0e-8);
   BOOST_CHECK_CLOSE(Su(5), mu(5), 4.0e-8);
   BOOST_CHECK_CLOSE(Su(6), mu(6), 3.0e-8);

   // down-type sleptons
   const DoubleVector Se(ToDoubleVector(m.get_physical().MSe));
   const DoubleVector me(s.displayPhys().me.flatten().sort());
   BOOST_CHECK_CLOSE(Se(1), me(1), 2.0e-9);
   BOOST_CHECK_CLOSE(Se(2), me(2), 2.0e-9);
   BOOST_CHECK_CLOSE(Se(3), me(3), 2.0e-9);
   BOOST_CHECK_CLOSE(Se(4), me(4), 2.0e-9);
   BOOST_CHECK_CLOSE(Se(5), me(5), 2.0e-9);
   BOOST_CHECK_CLOSE(Se(6), me(6), 2.0e-9);

   // up-type sleptons
   const DoubleVector Sv(ToDoubleVector(m.get_physical().MSv));
   const DoubleVector msnu(s.displayPhys().msnu.sort());
   BOOST_CHECK_CLOSE(Sv(1), msnu(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Sv(2), msnu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sv(3), msnu(3), 1.0e-12);

   // neutrinos
   const DoubleVector MFv(ToDoubleVector(m.get_physical().MFv));
   BOOST_CHECK_EQUAL(MFv(1), 0.0);
   BOOST_CHECK_EQUAL(MFv(2), 0.0);
   BOOST_CHECK_EQUAL(MFv(3), 0.0);

   // leptons
   const DoubleVector MFe(ToDoubleVector(m.get_physical().MFe));
   BOOST_CHECK_EQUAL(MFe(1), 0.0);
   BOOST_CHECK_EQUAL(MFe(2), 0.0);
   // BOOST_CHECK_CLOSE(MFe(3), s.displayPhys().mtau, 1.0e-12);

   // ups
   const DoubleVector MFu(ToDoubleVector(m.get_physical().MFu));
   BOOST_CHECK_EQUAL(MFu(1), 0.0);
   BOOST_CHECK_EQUAL(MFu(2), 0.0);
   // BOOST_CHECK_CLOSE(MFu(3), s.displayPhys().mt, 1.0e-12);

   // downs
   const DoubleVector MFd(ToDoubleVector(m.get_physical().MFd));
   BOOST_CHECK_EQUAL(MFd(1), 0.0);
   BOOST_CHECK_EQUAL(MFd(2), 0.0);
   // BOOST_CHECK_CLOSE(MFd(3), s.displayPhys().mb, 1.0e-12);

   ensure_one_loop_ewsb(m, s);
   m.calculate_pole_masses();
   s.physical(1);

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   // neutral CP even Higgs
   const DoubleVector hh(ToDoubleVector(m.get_physical().Mhh));
   const DoubleVector mh0(s.displayPhys().mh0);
   BOOST_CHECK_CLOSE(hh(1), mh0(1), 0.0015);
   BOOST_CHECK_CLOSE(hh(2), mh0(2), 2.0e-5);
   BOOST_CHECK_CLOSE(hh(3), mh0(3), 0.008);

   // neutral CP odd Higgs
   const DoubleVector Ah(ToDoubleVector(m.get_physical().MAh));
   const DoubleVector mA0(s.displayPhys().mA0);
   BOOST_CHECK_CLOSE(Ah(2), mA0(1), 3.0e-5);
   BOOST_CHECK_CLOSE(Ah(3), mA0(2), 7.0e-4);
   // charged Higgs
   const DoubleVector Hpm(ToDoubleVector(m.get_physical().MHpm));
   const double mHpm = s.displayPhys().mHpm;
   BOOST_CHECK_CLOSE(Hpm(2), mHpm, 0.0009);
}

BOOST_AUTO_TEST_CASE( test_self_energies )
{
   SMSSM_input_parameters input;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaSInput = 100.;
   input.m0       = 540.;
   input.Azero    = -350.;
   input.MSInput  = 290.;
   input.BMSInput = 400.;
   input.L1Input  = 300.;
   SMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   setup_SMSSM_const(m, s, input);

   softsusy::Z3 = false;
   softsusy::numHiggsMassLoops = 1;
   softsusy::MIXING = 1;
   ensure_tree_level_ewsb(m, s, input);

   m.calculate_DRbar_masses();
   s.calcDrBarPars();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   BOOST_CHECK_EQUAL(s.displayMu(), m.get_scale());

   const double gp = s.displayGaugeCoupling(1) * sqrt(0.6);
   const double g2 =  s.displayGaugeCoupling(2);
   const double sw = gp / sqrt(sqr(gp) + sqr(g2));

   const double nmT1 = s.doCalcTadpole1oneLoop(m.get_MFu()(2), sw);
   const double nmT2 = s.doCalcTadpole2oneLoop(m.get_MFu()(2), sw);
   const double nmTS = s.doCalcTadpoleSoneLoop(m.get_MFu()(2), sw);

   const double fsT1 = m.tadpole_hh_1loop(0).real() / m.get_vd();
   const double fsT2 = m.tadpole_hh_1loop(1).real() / m.get_vu();
   const double fsTS = m.tadpole_hh_1loop(2).real() / m.get_vS();

   BOOST_CHECK_CLOSE(nmT1, fsT1, 1.1e-12);
   BOOST_CHECK_CLOSE(nmT2, fsT2, 1.0e-12);
   BOOST_CHECK_CLOSE(nmTS, fsTS, 2.0e-12);

   /// CP-even Higgs self-energies
   const double momentum = 100.;
   const double nmHSE11 = s.pis1s1(momentum,s.displayMu());
   const double fsHSE11 = m.self_energy_hh_1loop(momentum, 0, 0).real();
   const double nmHSE12 = s.pis1s2(momentum,s.displayMu());
   const double fsHSE12 = m.self_energy_hh_1loop(momentum, 0, 1).real();
   const double nmHSE22 = s.pis2s2(momentum,s.displayMu());
   const double fsHSE22 = m.self_energy_hh_1loop(momentum, 1, 1).real();
   const double nmHSE13 = s.pis1s3(momentum,s.displayMu());
   const double fsHSE13 = m.self_energy_hh_1loop(momentum, 0, 2).real();
   const double nmHSE23 = s.pis2s3(momentum,s.displayMu());
   const double fsHSE23 = m.self_energy_hh_1loop(momentum, 1, 2).real();
   const double nmHSE33 = s.pis3s3(momentum,s.displayMu());
   const double fsHSE33 = m.self_energy_hh_1loop(momentum, 2, 2).real();

   BOOST_CHECK_CLOSE(nmHSE11, fsHSE11, 1.0e-12);
   BOOST_CHECK_CLOSE(nmHSE12, fsHSE12, 1.0e-12);
   BOOST_CHECK_CLOSE(nmHSE13, fsHSE13, 2.0e-12);
   BOOST_CHECK_CLOSE(nmHSE22, fsHSE22, 1.0e-12);
   BOOST_CHECK_CLOSE(nmHSE23, fsHSE23, 2.0e-12);
   BOOST_CHECK_CLOSE(nmHSE33, fsHSE33, 5.0e-11);

   /// CP-odd self-energies
   const double nmPHSE11 = s.pip1p1(momentum,s.displayMu());
   const double fsPHSE11 = m.self_energy_Ah_1loop(momentum, 0, 0).real();
   const double nmPHSE12 = s.pip1p2(momentum,s.displayMu());
   const double fsPHSE12 = m.self_energy_Ah_1loop(momentum, 0, 1).real();
   const double nmPHSE22 = s.pip2p2(momentum,s.displayMu());
   const double fsPHSE22 = m.self_energy_Ah_1loop(momentum, 1, 1).real();
   const double nmPHSE13 = s.pip1p3(momentum,s.displayMu());
   const double fsPHSE13 = m.self_energy_Ah_1loop(momentum, 0, 2).real();
   const double nmPHSE23 = s.pip2p3(momentum,s.displayMu());
   const double fsPHSE23 = m.self_energy_Ah_1loop(momentum, 1, 2).real();
   const double nmPHSE33 = s.pip3p3(momentum,s.displayMu());
   const double fsPHSE33 = m.self_energy_Ah_1loop(momentum, 2, 2).real();

   BOOST_CHECK_CLOSE(nmPHSE11, fsPHSE11, 9.0e-12);
   BOOST_CHECK_CLOSE(nmPHSE12, fsPHSE12, 7.0e-12);
   BOOST_CHECK_CLOSE(nmPHSE13, fsPHSE13, 7.0e-12);
   BOOST_CHECK_CLOSE(nmPHSE22, fsPHSE22, 2.0e-12);
   BOOST_CHECK_CLOSE(nmPHSE23, fsPHSE23, 2.0e-11);
   BOOST_CHECK_CLOSE(nmPHSE33, fsPHSE33, 1.0e-12);
}
