
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NUTNMSSM_spectrum

#include <boost/test/unit_test.hpp>

#define private public

#include "nmssmsoftsusy.h"
#include "error.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "NUTNMSSM_two_scale_model.hpp"
#include "NUTNMSSM_input_parameters.hpp"
#include "NUTNMSSM_two_scale_high_scale_constraint.hpp"
#include "NUTNMSSM_two_scale_susy_scale_constraint.hpp"
#include "NUTNMSSM_two_scale_low_scale_constraint.hpp"
#include "NUTNMSSM_two_scale_convergence_tester.hpp"
#include "NUTNMSSM_two_scale_initial_guesser.hpp"
#include "test_NUTNMSSM.hpp"

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

class SoftSusy_error : public Error {
public:
   SoftSusy_error(const std::string& msg)
      : Error(msg) {}
   virtual ~SoftSusy_error() {}
};

class SoftSusy_NoConvergence_error : public SoftSusy_error {
public:
   SoftSusy_NoConvergence_error(const std::string& msg_)
      : SoftSusy_error(msg_) {}
   virtual ~SoftSusy_NoConvergence_error() {}
};

class SoftSusy_NonPerturbative_error : public SoftSusy_error {
public:
   SoftSusy_NonPerturbative_error(const std::string& msg_)
      : SoftSusy_error(msg_) {}
   virtual ~SoftSusy_NonPerturbative_error() {}
};

class SoftSusy_tester {
public:
   SoftSusy_tester()
      : mx(0.0), msusy(0.0), softSusy(), gaugeUnification(true), loops(1) {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   sPhysical get_physical() const { return softSusy.displayPhys(); }
   NmssmSoftsusy get_model() const { return softSusy; }
   void reset() {
      mx = 0.0;
      msusy = 0.0;
      gaugeUnification = true;
      loops = 1;
   }
   void set_loops(int l) { loops = l; }
   void test(const NUTNMSSM_input_parameters& pp, double mxGuess, const QedQcd_legacy& qedqcd) {
      // run softsusy
      softsusy::numRewsbLoops = loops;
      softsusy::numHiggsMassLoops = loops;
      softsusy::TOLERANCE = 1.0e-4;
      softsusy::Z3 = true;
      softsusy::GUTlambda = false;
      softsusy::GUTkappa = false;
      softsusy::GUTsVev = false;
      softsusy::SoftHiggsOut = true;
#ifdef ENABLE_VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      DoubleVector pars(8);
      pars(1) = pp.m0;
      pars(2) = pp.m12;
      pars(3) = pp.Azero;
      pars(4) = 0.; // Mu
      pars(5) = 0.; // BMu / (Cos(Beta) Sin(Beta))
      pars(6) = 0.; // xiS
      pars(7) = pp.ALambdaInput;
      pars(8) = pp.AKappaInput;
      DoubleVector nmpars(5);
      nmpars(1) = pp.LambdaInput;
      nmpars(2) = pp.KappaInput;
      nmpars(3) = pp.MuEff * Sqrt(2.) / pp.LambdaInput;
      nmpars(4) = 0.;
      nmpars(5) = 0.;

      softSusy.setAlternativeMs(false);

      try {
         softSusy.lowOrg(NmssmSemiMsugraNoSoftHiggsMassBcs, mxGuess, pars,
                         nmpars, 1, pp.TanBeta, qedqcd, gaugeUnification);
      } catch (const std::string& str) {
         BOOST_TEST_MESSAGE("SoftSusy problem: " << str);
         throw SoftSusy_error(str);
      }

      mx = softSusy.displayMxBC();
      msusy = softSusy.displayMsusy();
      softsusy::PRINTOUT = 0;

      if (softSusy.displayProblem().test()) {
         std::stringstream ss;
         ss << "SoftSusy problem: " << softSusy.displayProblem();
         BOOST_TEST_MESSAGE(ss.str());
         if (softSusy.displayProblem().noConvergence)
            throw SoftSusy_NoConvergence_error(ss.str());
         else if (softSusy.displayProblem().nonperturbative)
            throw SoftSusy_NonPerturbative_error(ss.str());
         else
            throw SoftSusy_error(ss.str());
      }
   }
private:
   double mx, msusy;
   NmssmSoftsusy softSusy;
   bool gaugeUnification;
   int loops;
};

class NUTNMSSM_tester {
public:
   NUTNMSSM_tester()
      : mx(0.0), msusy(0.0), mssm()
      , high_constraint(NULL), susy_constraint(NULL), low_constraint(NULL)
      , loops(1) {}
   ~NUTNMSSM_tester() {
      delete high_constraint;
      delete susy_constraint;
      delete low_constraint;
   }
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   NUTNMSSM_physical get_physical() const { return mssm.get_physical(); }
   NUTNMSSM<Two_scale> get_model() const { return mssm; }
   void reset() {
      mx = 0.0;
      msusy = 0.0;
      mssm.clear();
      high_constraint = NULL;
      susy_constraint = NULL;
      low_constraint = NULL;
      loops = 1;
   }
   void set_loops(int l) { loops = l; }
   void set_low_scale_constraint(NUTNMSSM_low_scale_constraint<Two_scale>* c) { low_constraint = c; }
   void set_susy_scale_constraint(NUTNMSSM_susy_scale_constraint<Two_scale>* c) { susy_constraint = c; }
   void set_high_scale_constraint(NUTNMSSM_high_scale_constraint<Two_scale>* c) { high_constraint = c; }
   void setup_default_constaints(const NUTNMSSM_input_parameters& pp, const QedQcd_legacy& qedqcd) {
      if (!high_constraint)
         high_constraint = new NUTNMSSM_high_scale_constraint<Two_scale>(&mssm);
      if (!susy_constraint)
         susy_constraint = new NUTNMSSM_susy_scale_constraint<Two_scale>(&mssm, convert(qedqcd));
      if (!low_constraint)
         low_constraint = new NUTNMSSM_low_scale_constraint<Two_scale>(&mssm, convert(qedqcd));
   }
   void test(const NUTNMSSM_input_parameters& pp, const QedQcd_legacy& qedqcd) {
      setup_default_constaints(pp, qedqcd);

      mssm.clear();
      mssm.set_loops(2);
      mssm.set_thresholds(2);
      mssm.set_ewsb_loop_order(loops);
      mssm.set_pole_mass_loop_order(loops);
      mssm.set_input_parameters(pp);
      mssm.set_precision(1.0e-4); // == softsusy::TOLERANCE
      mssm.do_force_output(true);

      high_constraint->clear();
      susy_constraint->clear();
      low_constraint ->clear();
      high_constraint->set_model(&mssm);
      susy_constraint->set_model(&mssm);
      low_constraint ->set_model(&mssm);
      low_constraint ->set_sm_parameters(convert(qedqcd));
      high_constraint->initialize();
      susy_constraint->initialize();
      low_constraint ->initialize();

      NUTNMSSM_convergence_tester<Two_scale> convergence_tester(&mssm, 1.0e-4);
      convergence_tester.set_max_iterations(100);
      NUTNMSSM_initial_guesser<Two_scale> initial_guesser(&mssm, convert(qedqcd),
                                                      *low_constraint,
                                                      *susy_constraint,
                                                      *high_constraint);
      Two_scale_increasing_precision precision(10.0, 1.0e-6);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&convergence_tester);
      solver.set_running_precision(&precision);
      solver.set_initial_guesser(&initial_guesser);
      solver.add(low_constraint, &mssm);
      solver.add(high_constraint, &mssm);
      solver.add(susy_constraint, &mssm);

      try {
         solver.solve();
         mssm.run_to(susy_constraint->get_scale());
         mssm.solve_ewsb();
         mssm.calculate_spectrum();
         mssm.run_to(Electroweak_constants::MZ);
      } catch (const Error& error) {
         mssm.get_problems().flag_thrown(error.what());
         BOOST_TEST_MESSAGE("FlexibleSUSY error: " << error.what());
      } catch (const std::string& str) {
         mssm.get_problems().flag_thrown(str);
         BOOST_TEST_MESSAGE("FlexibleSUSY error: " << str);
      }

      mx = high_constraint->get_scale();
      msusy = susy_constraint->get_scale();
   }
private:
   double mx, msusy;
   NUTNMSSM<Two_scale> mssm;
   NUTNMSSM_high_scale_constraint<Two_scale>* high_constraint;
   NUTNMSSM_susy_scale_constraint<Two_scale>* susy_constraint;
   NUTNMSSM_low_scale_constraint<Two_scale>*  low_constraint;
   int loops;
};

void set_S1(NUTNMSSM_input_parameters& pp, softsusy::QedQcd_legacy& qedqcd)
{
   pp.m0 = 500.;
   pp.m12 = 500.;
   pp.TanBeta = 10.;
   pp.Azero = -1500.;
   pp.LambdaInput = 0.1;
   pp.KappaInput = 0.11;
   pp.ALambdaInput = -1500.;
   pp.AKappaInput = -36.;
   pp.MuEff = 965;

   qedqcd.setAlpha(legacy::ALPHA , 1./127.944);
   qedqcd.setAlpha(legacy::ALPHAS, 1.185e-01);
   softsusy::GMU = 1.1663787e-5;
   softsusy::MZ = 91.1876;
   qedqcd.setPoleMZ(softsusy::MZ);
   qedqcd.setMass(legacy::mBottom, 4.18000000E+00);
   qedqcd.setMbMb(4.18000000E+00);
   qedqcd.setPoleMt(1.73070000E+02);
   qedqcd.setMass(legacy::mTau, 1.77682);
   qedqcd.setPoleMtau(1.77682);

   qedqcd.toMz();
}

void set_BP1(NUTNMSSM_input_parameters& pp, softsusy::QedQcd_legacy& qedqcd)
{
   pp.m0 = 2400;
   pp.m12 = 550;
   pp.TanBeta = 2.665;
   pp.Azero = -972.1;
   pp.LambdaInput = 0.646;
   pp.KappaInput = 0.377;
   pp.ALambdaInput = -511.0;
   pp.AKappaInput = -845.7;
   pp.MuEff = 120.5;

   qedqcd.setAlpha(legacy::ALPHA , 1./127.944);
   qedqcd.setAlpha(legacy::ALPHAS, 1.185e-01);
   softsusy::GMU = 1.1663787e-5;
   softsusy::MZ = 91.1876;
   qedqcd.setPoleMZ(softsusy::MZ);
   qedqcd.setMass(legacy::mBottom, 4.18000000E+00);
   qedqcd.setMbMb(4.18000000E+00);
   qedqcd.setPoleMt(1.73070000E+02);
   qedqcd.setMass(legacy::mTau, 1.77682);
   qedqcd.setPoleMtau(1.77682);

   qedqcd.toMz();
}

void set_BP2(NUTNMSSM_input_parameters& pp, softsusy::QedQcd_legacy& qedqcd)
{
   pp.m0 = 2450;
   pp.m12 = 550;
   pp.TanBeta = 4.229;
   pp.Azero = -1923.9;
   pp.LambdaInput = 0.683;
   pp.KappaInput = 0.093;
   pp.ALambdaInput = 1774.9;
   pp.AKappaInput = 2533.4;
   pp.MuEff = 229.2;

   qedqcd.setAlpha(legacy::ALPHA , 1./127.944);
   qedqcd.setAlpha(legacy::ALPHAS, 1.185e-01);
   softsusy::GMU = 1.1663787e-5;
   softsusy::MZ = 91.1876;
   qedqcd.setPoleMZ(softsusy::MZ);
   qedqcd.setMass(legacy::mBottom, 4.18000000E+00);
   qedqcd.setMbMb(4.18000000E+00);
   qedqcd.setPoleMt(1.73070000E+02);
   qedqcd.setMass(legacy::mTau, 1.77682);
   qedqcd.setPoleMtau(1.77682);

   qedqcd.toMz();
}

void set_BP3(NUTNMSSM_input_parameters& pp, softsusy::QedQcd_legacy& qedqcd)
{
   pp.m0 = 2400;
   pp.m12 = 600;
   pp.TanBeta = 3.042;
   pp.Azero = -1956.4;
   pp.LambdaInput = 0.650;
   pp.KappaInput = 0.164;
   pp.ALambdaInput = 763.8;
   pp.AKappaInput = 1268.2;
   pp.MuEff = 265.2;

   qedqcd.setAlpha(legacy::ALPHA , 1./127.944);
   qedqcd.setAlpha(legacy::ALPHAS, 1.185e-01);
   softsusy::GMU = 1.1663787e-5;
   softsusy::MZ = 91.1876;
   qedqcd.setPoleMZ(softsusy::MZ);
   qedqcd.setMass(legacy::mBottom, 4.18000000E+00);
   qedqcd.setMbMb(4.18000000E+00);
   qedqcd.setPoleMt(1.73070000E+02);
   qedqcd.setMass(legacy::mTau, 1.77682);
   qedqcd.setPoleMtau(1.77682);

   qedqcd.toMz();
}

void compare_tadpoles_0loop(NUTNMSSM<Two_scale> fs, NmssmSoftsusy ss)
{
   copy_parameters(fs, ss);

   ss.setTadpole1Ms(0.);
   ss.setTadpole2Ms(0.);
   ss.setTadpoleSMs(0.);

   softsusy::SoftHiggsOut = true;
   ss.rewsbTreeLevel(1);
   fs.solve_ewsb_tree_level();

   BOOST_CHECK_CLOSE_FRACTION(ss.displayMh1Squared(), fs.get_mHd2(), 1.e-10);
   BOOST_CHECK_CLOSE_FRACTION(ss.displayMh2Squared(), fs.get_mHu2(), 1.e-10);
   BOOST_CHECK_CLOSE_FRACTION(ss.displayMsSquared() , fs.get_ms2() , 1.e-10);
}

void compare_tadpoles_1loop(NUTNMSSM<Two_scale> fs, NmssmSoftsusy ss)
{
   copy_parameters(fs, ss);

   ss.setTadpole1Ms(0.);
   ss.setTadpole2Ms(0.);
   ss.setTadpoleSMs(0.);

   softsusy::SoftHiggsOut = true;
   ss.calcDrBarPars();
   fs.set_ewsb_loop_order(1);
   fs.calculate_DRbar_masses();

   const double mt = ss.displayDrBarPars().mt;
   const double sinthDRbar = ss.calcSinthdrbar();

   ss.calcTadpole1Ms1loop(mt, sinthDRbar);
   ss.calcTadpole2Ms1loop(mt, sinthDRbar);
   ss.calcTadpoleSMs1loop(mt, sinthDRbar);

   const double td = Re(fs.tadpole_hh_1loop(0));
   const double tu = Re(fs.tadpole_hh_1loop(1));
   const double ts = Re(fs.tadpole_hh_1loop(2));

   const double vd = fs.get_vd();
   const double vu = fs.get_vu();
   const double vs = fs.get_vS();

   // check equality of tadpoles
   BOOST_CHECK_CLOSE_FRACTION(td / vd, ss.displayTadpole1Ms1loop(), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(tu / vu, ss.displayTadpole2Ms1loop(), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(ts / vs, ss.displayTadpoleSMs1loop(), 1.0e-6);
}

void compare_tadpoles_2loop(NUTNMSSM<Two_scale> fs, NmssmSoftsusy ss)
{
   copy_parameters(fs, ss);

   ss.setTadpole1Ms(0.);
   ss.setTadpole2Ms(0.);
   ss.setTadpoleSMs(0.);

   softsusy::SoftHiggsOut = true;

   ss.calcDrBarPars();
   fs.calculate_DRbar_masses();

   const double mt = ss.displayDrBarPars().mt;
   const double sinthDRbar = ss.calcSinthdrbar();
   const double vd = fs.get_vd();
   const double vu = fs.get_vu();
   const double vs = fs.get_vS();

   const double td_fs = fs.tadpole_hh_1loop(0).real();
   const double tu_fs = fs.tadpole_hh_1loop(1).real();
   const double ts_fs = fs.tadpole_hh_1loop(2).real();

   double td_ss = ss.doCalcTadpole1oneLoop(mt, sinthDRbar);
   double tu_ss = ss.doCalcTadpole2oneLoop(mt, sinthDRbar);
   double ts_ss = ss.doCalcTadpoleSoneLoop(mt, sinthDRbar);

   // check equality of 1-loop tadpoles
   BOOST_CHECK_CLOSE_FRACTION(td_fs / vd, td_ss, 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(tu_fs / vu, tu_ss, 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(ts_fs / vs, ts_ss, 1.0e-6);

   // make sure the one-loop tadpoles are calculated correctly
   softsusy::numRewsbLoops = 1;
   ss.doTadpoles(mt, sinthDRbar);

   td_ss = ss.displayTadpole1Ms();
   tu_ss = ss.displayTadpole2Ms();
   ts_ss = ss.displayTadpoleSMs();

   // check equality of 1-loop tadpoles again
   BOOST_CHECK_CLOSE_FRACTION(td_fs / vd, td_ss, 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(tu_fs / vu, tu_ss, 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(ts_fs / vs, ts_ss, 1.0e-6);

   // calculate 2-loop tadpoles
   softsusy::numRewsbLoops = 2;
   ss.doTadpoles(mt, sinthDRbar);

   const double td_1_and_2loop_ss = ss.displayTadpole1Ms();
   const double tu_1_and_2loop_ss = ss.displayTadpole2Ms();
   const double ts_1_and_2loop_ss = ss.displayTadpoleSMs();

   const auto two_loop_tadpole(fs.tadpole_hh_2loop());

   // check equality of 1-loop tadpoles again
   // works only if amu = the lightest CP-even Higgs,
   // but not the goldstone boson
   BOOST_CHECK_CLOSE_FRACTION(two_loop_tadpole[0] / vd, td_1_and_2loop_ss - td_ss, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(two_loop_tadpole[1] / vu, tu_1_and_2loop_ss - tu_ss, 1.0e-11);
   BOOST_CHECK_CLOSE_FRACTION(two_loop_tadpole[2] / vs, ts_1_and_2loop_ss - ts_ss, 1.0e-11);
}

BOOST_AUTO_TEST_CASE( test_NUTNMSSM_spectrum )
{
   NUTNMSSM_input_parameters pp;
   softsusy::QedQcd_legacy qedqcd;
   set_S1(pp, qedqcd);

   NUTNMSSM<Two_scale> _model(pp);
   const NUTNMSSM_high_scale_constraint<Two_scale> high_constraint(&_model);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   NUTNMSSM_tester nmssm_tester;
   BOOST_REQUIRE_NO_THROW(nmssm_tester.test(pp, qedqcd));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess, qedqcd));

   BOOST_CHECK_CLOSE_FRACTION(nmssm_tester.get_mx(), softSusy_tester.get_mx(), 0.18);
   BOOST_CHECK_CLOSE_FRACTION(nmssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.006);

   // compare model parameters
   const NmssmSoftsusy ss(softSusy_tester.get_model());
   const NUTNMSSM<Two_scale> fs(nmssm_tester.get_model());

   compare_tadpoles_0loop(fs, ss);
   compare_tadpoles_1loop(fs, ss);
   compare_tadpoles_2loop(fs, ss);

   BOOST_CHECK_EQUAL(ss.displayLoops()     , fs.get_loops());
   BOOST_CHECK_EQUAL(ss.displayMu()        , fs.get_scale());
   BOOST_CHECK(fs.get_thresholds());
   BOOST_CHECK(ss.displayThresholds());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00076);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0015);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.00013);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(0,0), ss.displayYukawaMatrix(YU)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(1,1), ss.displayYukawaMatrix(YU)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(2,2), ss.displayYukawaMatrix(YU)(3,3), 0.0012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(0,0), ss.displayYukawaMatrix(YD)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(1,1), ss.displayYukawaMatrix(YD)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(2,2), ss.displayYukawaMatrix(YD)(3,3), 0.0075);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(0,0), ss.displayYukawaMatrix(YE)(1,1), 0.0093);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(1,1), ss.displayYukawaMatrix(YE)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(2,2), ss.displayYukawaMatrix(YE)(3,3), 0.0096);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Kappa()  , ss.displayKappa(), 0.00001);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Lambdax(), ss.displayLambda(), 0.00003);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_vS()     , ss.displaySvev(), 0.00001);

   BOOST_CHECK_EQUAL(ss.displaySusyMu()    , 0.);
   BOOST_CHECK_EQUAL(ss.displayM3Squared() , 0.);
   BOOST_CHECK_EQUAL(ss.displayMspSquared(), 0.);
   BOOST_CHECK_EQUAL(ss.displayXiS()       , 0.);
   BOOST_CHECK_EQUAL(ss.displayXiF()       , 0.);
   BOOST_CHECK_EQUAL(ss.displayMupr()      , 0.);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassB() , ss.displayGaugino(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassWB(), ss.displayGaugino(2), 0.0046);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassG() , ss.displayGaugino(3), 0.0051);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.07);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.009);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ms2() , ss.displayMsSquared() , 0.015);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(0,0), ss.displaySoftMassSquared(mQl)(1,1), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(1,1), ss.displaySoftMassSquared(mQl)(2,2), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(2,2), ss.displaySoftMassSquared(mQl)(3,3), 0.012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(0,0), ss.displaySoftMassSquared(mUr)(1,1), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(1,1), ss.displaySoftMassSquared(mUr)(2,2), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(2,2), ss.displaySoftMassSquared(mUr)(3,3), 0.012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(0,0), ss.displaySoftMassSquared(mDr)(1,1), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(1,1), ss.displaySoftMassSquared(mDr)(2,2), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(2,2), ss.displaySoftMassSquared(mDr)(3,3), 0.012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(0,0), ss.displaySoftMassSquared(mLl)(1,1), 0.01);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(1,1), ss.displaySoftMassSquared(mLl)(2,2), 0.01);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(2,2), ss.displaySoftMassSquared(mLl)(3,3), 0.01);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(0,0), ss.displaySoftMassSquared(mEr)(1,1), 0.0025);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(1,1), ss.displaySoftMassSquared(mEr)(2,2), 0.0025);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(2,2), ss.displaySoftMassSquared(mEr)(3,3), 0.0013);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(0,0), ss.displayTrilinear(UA)(1,1), 0.018);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(1,1), ss.displayTrilinear(UA)(2,2), 0.018);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(2,2), ss.displayTrilinear(UA)(3,3), 0.008);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(0,0), ss.displayTrilinear(DA)(1,1), 0.019);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(1,1), ss.displayTrilinear(DA)(2,2), 0.019);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(2,2), ss.displayTrilinear(DA)(3,3), 0.016);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(0,0), ss.displayTrilinear(EA)(1,1), 0.021);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(1,1), ss.displayTrilinear(EA)(2,2), 0.021);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(2,2), ss.displayTrilinear(EA)(3,3), 0.021);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_TLambdax(), ss.displayTrialambda(), 0.001);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TKappa()  , ss.displayTriakappa() , 0.023);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double vs = fs.get_vS();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 3.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(vs     , ss.displaySvev(), 0.004);

   // comparing tree-level masses

   const DoubleVector MCha(ToDoubleVector(fs.get_MCha())),
      MChi(ToDoubleVector(fs.get_MChi())),
      MHpm(ToDoubleVector(fs.get_MHpm())),
      MAh(ToDoubleVector(fs.get_MAh())),
      Mhh(ToDoubleVector(fs.get_Mhh()));
   const DoubleVector mch(ss.displayDrBarPars().mchBpmz),
      mn(ss.displayDrBarPars().mnBpmz);
   const double MwRun = fs.get_MVWm();
   const double MzRun = fs.get_MVZ();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const DoubleVector mA0(ss.displayDrBarPars().mA0);
   const DoubleVector mh0(ss.displayDrBarPars().mh0);

   BOOST_TEST_MESSAGE("SoftSUSY    :\n mh_tree = " << mh0 << " mA_tree = " << mA0);
   BOOST_TEST_MESSAGE("FlexibleSUSY:\n mh_tree = " << Mhh << " mA_tree = " << MAh);

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha(1), mch(1), 0.002);
   BOOST_CHECK_CLOSE_FRACTION(MCha(2), mch(2), 0.008);

   BOOST_CHECK_CLOSE_FRACTION(MChi(1), mn(1), 0.0060);
   BOOST_CHECK_CLOSE_FRACTION(MChi(2), mn(2), 0.0041);
   BOOST_CHECK_CLOSE_FRACTION(MChi(3), mn(3), 0.0045);
   BOOST_CHECK_CLOSE_FRACTION(MChi(4), mn(4), 0.0040);
   BOOST_CHECK_CLOSE_FRACTION(MChi(5), mn(5), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun , 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0(1), 0.008);
   BOOST_CHECK_CLOSE_FRACTION(MAh(3), mA0(2), 0.0032);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0(1), 0.0007);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mh0(2), 0.0032);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(3), mh0(3), 0.004);

   // down-type squarks
   const DoubleVector Sd(ToDoubleVector(fs.get_MSd()));
   const DoubleVector md(ss.displayDrBarPars().md.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Sd(1), md(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(2), md(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(3), md(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(4), md(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(5), md(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(6), md(6), 0.006);

   // up-type squarks
   const DoubleVector Su(ToDoubleVector(fs.get_MSu()));
   const DoubleVector mu(ss.displayDrBarPars().mu.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Su(1), mu(1), 0.0061);
   BOOST_CHECK_CLOSE_FRACTION(Su(2), mu(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(3), mu(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(4), mu(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(5), mu(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(6), mu(6), 0.006);

   // sleptons
   const DoubleVector Se(ToDoubleVector(fs.get_MSe()));
   const DoubleVector me(ss.displayDrBarPars().me.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Se(1), me(1), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(2), me(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(3), me(3), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(4), me(4), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(Se(5), me(5), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(Se(6), me(6), 0.005);

   // sneutrinos
   const DoubleVector msnu(ss.displayDrBarPars().msnu.sort());
   const DoubleVector Snu(ToDoubleVector(fs.get_MSv()));
   BOOST_CHECK_CLOSE_FRACTION(Snu(1), msnu(1), 0.0055);
   BOOST_CHECK_CLOSE_FRACTION(Snu(2), msnu(2), 0.0055);
   BOOST_CHECK_CLOSE_FRACTION(Snu(3), msnu(3), 0.0085);

   BOOST_CHECK_EQUAL(fs.get_MVP(), 0.0);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MVZ() , MzRun, 0.022);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MVWm(), MwRun, 0.024);

   BOOST_CHECK_EQUAL(fs.get_MVG(), 0.0);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MGlu(), ss.displayDrBarPars().mGluino, 0.005);

   BOOST_CHECK_EQUAL(fs.get_MFv()(0), 0.0);
   BOOST_CHECK_EQUAL(fs.get_MFv()(1), 0.0);
   BOOST_CHECK_EQUAL(fs.get_MFv()(2), 0.0);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFe()(2), ss.displayDrBarPars().mtau, 0.001);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFu()(2), ss.displayDrBarPars().mt  , 0.0097);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFd()(2), ss.displayDrBarPars().mb  , 0.0027);


   // comparing pole masses

   const DoubleVector
      MCha_1l(ToDoubleVector(fs.get_physical().MCha)),
      MChi_1l(ToDoubleVector(fs.get_physical().MChi)),
      MHpm_1l(ToDoubleVector(fs.get_physical().MHpm)),
      MAh_1l(ToDoubleVector(fs.get_physical().MAh)),
      Mhh_1l(ToDoubleVector(fs.get_physical().Mhh));
   const DoubleVector mch_1l(ss.displayPhys().mch),
      mn_1l(ss.displayPhys().mneut.apply(fabs));
   const double mHpm_1l = ss.displayPhys().mHpm;
   const DoubleVector mA0_1l(ss.displayPhys().mA0);
   const DoubleVector mh0_1l(ss.displayPhys().mh0);

   BOOST_TEST_MESSAGE("SoftSUSY    :\n mh_1l = " << mh0_1l << " mA_1l = " << mA0_1l);
   BOOST_TEST_MESSAGE("FlexibleSUSY:\n mh_1l = " << Mhh_1l << " mA_1l = " << MAh_1l);

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha_1l(1), mch_1l(1), 0.0016);
   BOOST_CHECK_CLOSE_FRACTION(MCha_1l(2), mch_1l(2), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(1), mn_1l(1), 0.0058);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(2), mn_1l(2), 0.0017);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(3), mn_1l(3), 0.0044);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(4), mn_1l(4), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MHpm_1l(2), mHpm_1l  , 0.0032);
   BOOST_CHECK_CLOSE_FRACTION(MAh_1l(2) , mA0_1l(1), 0.008);
   BOOST_CHECK_CLOSE_FRACTION(MAh_1l(3) , mA0_1l(2), 0.0032);

   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(1), mh0_1l(1), 0.0002);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(2), mh0_1l(2), 0.0033);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(3), mh0_1l(3), 0.001);

   // down-type squarks
   const DoubleVector Sd_1l(ToDoubleVector(fs.get_physical().MSd));
   const DoubleVector md_1l(ss.displayPhys().md.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(1), md_1l(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(2), md_1l(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(3), md_1l(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(4), md_1l(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(5), md_1l(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(6), md_1l(6), 0.006);

   // up-type squarks
   const DoubleVector Su_1l(ToDoubleVector(fs.get_physical().MSu));
   const DoubleVector mu_1l(ss.displayPhys().mu.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(1), mu_1l(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(2), mu_1l(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(3), mu_1l(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(4), mu_1l(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(5), mu_1l(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(6), mu_1l(6), 0.006);

   // sleptons
   const DoubleVector Se_1l(ToDoubleVector(fs.get_physical().MSe));
   const DoubleVector me_1l(ss.displayPhys().me.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(1), me_1l(1), 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(2), me_1l(2), 0.0008);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(3), me_1l(3), 0.0008);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(4), me_1l(4), 0.0050);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(5), me_1l(5), 0.0050);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(6), me_1l(6), 0.0050);

   // sneutrinos
   const DoubleVector msnu_1l(ss.displayPhys().msnu.sort());
   const DoubleVector Snu_1l(ToDoubleVector(fs.get_physical().MSv));
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(1), msnu_1l(1), 0.0056);
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(2), msnu_1l(2), 0.0056);
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(3), msnu_1l(3), 0.0090);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_physical().MGlu, ss.displayPhys().mGluino, 0.005);

   // comparing 2-loop pole masses

   nmssm_tester.reset();
   nmssm_tester.set_loops(2);
   softSusy_tester.reset();
   softSusy_tester.set_loops(2);
   BOOST_REQUIRE_NO_THROW(nmssm_tester.test(pp, qedqcd));
   BOOST_CHECK_NO_THROW(softSusy_tester.test(pp, mxGuess, qedqcd));

   // compare model parameters
   const NmssmSoftsusy ss_2l(softSusy_tester.get_model());
   const NUTNMSSM<Two_scale> fs_2l(nmssm_tester.get_model());

   const DoubleVector
      MHpm_2l(ToDoubleVector(fs_2l.get_physical().MHpm)),
      MAh_2l(ToDoubleVector(fs_2l.get_physical().MAh)),
      Mhh_2l(ToDoubleVector(fs_2l.get_physical().Mhh));
   const double mHpm_2l = ss_2l.displayPhys().mHpm;
   const DoubleVector mA_2l(ss_2l.displayPhys().mA0);
   const DoubleVector mh_2l(ss_2l.displayPhys().mh0);

   BOOST_CHECK_EQUAL(fs_2l.get_loops(), 2);
   BOOST_CHECK_EQUAL(fs_2l.get_loops(), ss_2l.displayLoops());

   BOOST_CHECK_CLOSE_FRACTION(MHpm_2l(2), mHpm_2l, 0.004);

   BOOST_CHECK_CLOSE_FRACTION(MAh_2l(2), mA_2l(1), 0.008);
   BOOST_CHECK_CLOSE_FRACTION(MAh_2l(3), mA_2l(2), 0.004);

   BOOST_CHECK_CLOSE_FRACTION(Mhh_2l(1), mh_2l(1), 8.0-05);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_2l(2), mh_2l(2), 0.004);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_2l(3), mh_2l(3), 0.0005);

   BOOST_TEST_MESSAGE("SoftSUSY    :\n mh_2l = " << mh_2l  << " mA_2l = " << mA_2l);
   BOOST_TEST_MESSAGE("FlexibleSUSY:\n mh_2l = " << Mhh_2l << " mA_2l = " << MAh_2l);
}

void test_NUTNMSSM_spectrum_with_fermi_constant_input_for_point(
   const NUTNMSSM_input_parameters& pp,
   const softsusy::QedQcd_legacy& qedqcd)
{
   NUTNMSSM<Two_scale> _model(pp);
   const NUTNMSSM_high_scale_constraint<Two_scale> high_constraint(&_model);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   NUTNMSSM_tester nmssm_tester;
   BOOST_REQUIRE_NO_THROW(nmssm_tester.test(pp, qedqcd));

   SoftSusy_tester softSusy_tester;
   BOOST_CHECK_NO_THROW(softSusy_tester.test(pp, mxGuess, qedqcd));

   BOOST_CHECK_CLOSE_FRACTION(nmssm_tester.get_mx(), softSusy_tester.get_mx(), 0.04);
   BOOST_CHECK_CLOSE_FRACTION(nmssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.0004);

   // compare model parameters
   const NmssmSoftsusy ss(softSusy_tester.get_model());
   const NUTNMSSM<Two_scale> fs(nmssm_tester.get_model());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 1.0e-04);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0003);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 2.0e-05);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Kappa() , ss.displayKappa(), 3e-07);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_vS() , ss.displaySvev(), 1.0e-7);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.07);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 1.0e-04);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ms2(), ss.displayMsSquared(), 2e-05);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 8.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0002);

   // comparing tree-level masses

   const DoubleVector MHpm(ToDoubleVector(fs.get_MHpm())),
      MAh(ToDoubleVector(fs.get_MAh())),
      Mhh(ToDoubleVector(fs.get_Mhh()));
   const double MwRun = fs.get_MVWm();
   const double MzRun = fs.get_MVZ();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const DoubleVector mA(ss.displayDrBarPars().mA0);
   const DoubleVector mh(ss.displayDrBarPars().mh0);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm, 0.003);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA(1), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(MAh(3), mA(2), 0.003);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh(1), 0.00015);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mh(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(3), mh(3), 0.0005);

   BOOST_TEST_MESSAGE("SoftSUSY    :\n mh_tree = " << mh  << " mA_tree = " << mA);
   BOOST_TEST_MESSAGE("FlexibleSUSY:\n mh_tree = " << Mhh << " mA_tree = " << MAh);

   // comparing 1-loop pole masses

   const DoubleVector
      MHpm_1l(ToDoubleVector(fs.get_physical().MHpm)),
      MAh_1l(ToDoubleVector(fs.get_physical().MAh)),
      Mhh_1l(ToDoubleVector(fs.get_physical().Mhh));
   const double mHpm_1l = ss.displayPhys().mHpm;
   const DoubleVector mA_1l(ss.displayPhys().mA0);
   const DoubleVector mh_1l(ss.displayPhys().mh0);

   BOOST_CHECK_CLOSE_FRACTION(MHpm_1l(2), mHpm_1l, 0.003);

   BOOST_CHECK_CLOSE_FRACTION(MAh_1l(2), mA_1l(1), 0.004);
   BOOST_CHECK_CLOSE_FRACTION(MAh_1l(3), mA_1l(2), 0.003);

   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(1), mh_1l(1), 0.0002);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(2), mh_1l(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(3), mh_1l(3), 0.0005);

   BOOST_TEST_MESSAGE("SoftSUSY    :\n mh_1l = " << mh_1l  << " mA_1l = " << mA_1l);
   BOOST_TEST_MESSAGE("FlexibleSUSY:\n mh_1l = " << Mhh_1l << " mA_1l = " << MAh_1l);

   // comparing 2-loop pole masses

   nmssm_tester.reset();
   nmssm_tester.set_loops(2);
   softSusy_tester.reset();
   softSusy_tester.set_loops(2);
   BOOST_REQUIRE_NO_THROW(nmssm_tester.test(pp, qedqcd));
   BOOST_CHECK_NO_THROW(softSusy_tester.test(pp, mxGuess, qedqcd));

   // compare model parameters
   const NmssmSoftsusy ss_2l(softSusy_tester.get_model());
   const NUTNMSSM<Two_scale> fs_2l(nmssm_tester.get_model());

   const DoubleVector
      MHpm_2l(ToDoubleVector(fs_2l.get_physical().MHpm)),
      MAh_2l(ToDoubleVector(fs_2l.get_physical().MAh)),
      Mhh_2l(ToDoubleVector(fs_2l.get_physical().Mhh));
   const double mHpm_2l = ss_2l.displayPhys().mHpm;
   const DoubleVector mA_2l(ss_2l.displayPhys().mA0);
   const DoubleVector mh_2l(ss_2l.displayPhys().mh0);

   BOOST_CHECK_EQUAL(fs_2l.get_loops(), 2);
   BOOST_CHECK_EQUAL(fs_2l.get_loops(), ss_2l.displayLoops());

   BOOST_CHECK_CLOSE_FRACTION(MHpm_2l(2), mHpm_2l, 0.003);

   BOOST_CHECK_CLOSE_FRACTION(MAh_2l(2), mA_2l(1), 0.004);
   BOOST_CHECK_CLOSE_FRACTION(MAh_2l(3), mA_2l(2), 0.004);

   BOOST_CHECK_CLOSE_FRACTION(Mhh_2l(1), mh_2l(1), 0.0002);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_2l(2), mh_2l(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_2l(3), mh_2l(3), 0.0005);

   BOOST_TEST_MESSAGE("SoftSUSY    :\n mh_2l = " << mh_2l  << " mA_2l = " << mA_2l);
   BOOST_TEST_MESSAGE("FlexibleSUSY:\n mh_2l = " << Mhh_2l << " mA_2l = " << MAh_2l);
}

BOOST_AUTO_TEST_CASE( test_NUTNMSSM_spectrum_with_fermi_constant_input )
{
   // standard NUTNMSSM testing point S1
   {
      BOOST_TEST_MESSAGE("testing S1 ...");
      softsusy::QedQcd_legacy qedqcd;
      NUTNMSSM_input_parameters pp;
      set_S1(pp, qedqcd);
      test_NUTNMSSM_spectrum_with_fermi_constant_input_for_point(pp, qedqcd);
   }

   // // NUTNMSSM point BP1
   // {
   //    BOOST_TEST_MESSAGE("testing BP1 ...");
   //    softsusy::QedQcd_legacy qedqcd;
   //    NUTNMSSM_input_parameters pp;
   //    set_BP1(pp, qedqcd);
   //    test_NUTNMSSM_spectrum_with_fermi_constant_input_for_point(pp, qedqcd);
   // }

   // // NUTNMSSM point BP2
   // {
   //    BOOST_TEST_MESSAGE("testing BP2 ...");
   //    softsusy::QedQcd_legacy qedqcd;
   //    NUTNMSSM_input_parameters pp;
   //    set_BP2(pp, qedqcd);
   //    test_NUTNMSSM_spectrum_with_fermi_constant_input_for_point(pp, qedqcd);
   // }

   // // NUTNMSSM point BP3
   // {
   //    BOOST_TEST_MESSAGE("testing BP3 ...");
   //    softsusy::QedQcd_legacy qedqcd;
   //    NUTNMSSM_input_parameters pp;
   //    set_BP3(pp, qedqcd);
   //    test_NUTNMSSM_spectrum_with_fermi_constant_input_for_point(pp, qedqcd);
   // }
}
