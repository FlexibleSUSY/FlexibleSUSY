
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMMassWInput_spectrum

#include <boost/test/unit_test.hpp>

#define private public

#include "softsusy.h"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "weinberg_angle.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "CMSSMMassWInput_two_scale_model.hpp"
#include "CMSSMMassWInput_input_parameters.hpp"
#include "CMSSMMassWInput_two_scale_high_scale_constraint.hpp"
#include "CMSSMMassWInput_two_scale_susy_scale_constraint.hpp"
#include "CMSSMMassWInput_two_scale_low_scale_constraint.hpp"
#include "CMSSMMassWInput_two_scale_convergence_tester.hpp"
#include "CMSSMMassWInput_two_scale_initial_guesser.hpp"
#include "test_CMSSM_like.hpp"

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

softsusy::QedQcd_legacy convert(const softsusy::QedQcd& ql)
{
   softsusy::QedQcd_legacy qn;

   qn.setAlphas(flexiblesusy::ToDoubleVector(ql.displayAlphas()));
   qn.setMasses(flexiblesusy::ToDoubleVector(ql.displayMass()));
   qn.set_input(ql.display_input());
   qn.setPoleMb(ql.displayPoleMb());
   qn.setCKM(ql.displayCKM());
   qn.setPMNS(ql.displayPMNS());
   qn.setPars(ql.get_number_of_parameters());
   qn.setMu(ql.get_scale());
   qn.setLoops(ql.get_loops());
   qn.setThresholds(ql.get_thresholds());

   return qn;
}

using namespace weinberg_angle;

/**
 * @class CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint
 *
 * Replacement class for CMSSMMassWInput_low_scale_constraint, which calculates
 * the gauge couplings at the low scale as Softsusy does it.
 */
class CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint
   : public CMSSMMassWInput_low_scale_constraint<Two_scale> {
public:
   CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint()
      : CMSSMMassWInput_low_scale_constraint<Two_scale>() {}
   CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint(CMSSMMassWInput<Two_scale>* model_, const CMSSMMassWInput_input_parameters&, const QedQcd_legacy& qedqcd_)
      : CMSSMMassWInput_low_scale_constraint<Two_scale>(model_,convert(qedqcd_)) {}
   virtual ~CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint() {}

   virtual void apply();
};

void CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint::apply()
{
   assert(model && "Error: CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint:"
          " model pointer must not be zero");

   // save old model parmeters
   const CMSSMMassWInput<Two_scale> mssm(*model);

   // run CMSSMMassWInput_low_scale_constraint::apply(), without the gauge
   // couplings
   model->calculate_DRbar_masses();

   if (model->get_problems().have_problem()) {
      BOOST_TEST_MESSAGE("Problem in CMSSMMassWInput_precise_gauge_couplings_"
                    "low_scale_constraint::apply(): "
                    << model->get_problems());
   }

   update_scale();
   calculate_DRbar_gauge_couplings();
   calculate_running_SM_masses();

   const double TanBeta = model->get_input().TanBeta;
   const double g1 = model->get_g1();
   const double g2 = model->get_g2();

   model->set_vd((2*mZ_run)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
      )));
   model->set_vu((2*mZ_run*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta))));

   calculate_DRbar_yukawa_couplings();

   const Eigen::Matrix<double,3,3> new_Yu(model->get_Yu()),
      new_Yd(model->get_Yd()), new_Ye(model->get_Ye());

   // Now calculate the gauge couplings using
   // MssmSoftsusy::sparticleThresholdCorrections
   MssmSoftsusy softsusy;
   copy_parameters(mssm, softsusy);
   softsusy.setData(convert(qedqcd));
   softsusy.setMw(qedqcd.displayPoleMW());

   // prevent tan(beta) from being reset
   softsusy.setSetTbAtMX(true);
   softsusy.sparticleThresholdCorrections(model->get_input().TanBeta);

   if (softsusy.displayProblem().test()) {
      std::ostringstream ss;
      ss << "Softsusy problem in CMSSMMassWInput_precise_gauge_couplings_"
           "low_scale_constraint::apply(): "
           "error while calculating the sparticle thresholds: "
         << softsusy.displayProblem();
      BOOST_TEST_MESSAGE(ss.str());
   }

   BOOST_TEST_MESSAGE("Difference (g1_FlexibleSUSY - g1_softsusy)(MZ) = "
                 << new_g1 - softsusy.displayGaugeCoupling(1));
   BOOST_TEST_MESSAGE("Difference (g2_FlexibleSUSY - g2_softsusy)(MZ) = "
                 << new_g2 - softsusy.displayGaugeCoupling(2));
   BOOST_TEST_MESSAGE("Difference (g3_FlexibleSUSY - g3_softsusy)(MZ) = "
                 << new_g3 - softsusy.displayGaugeCoupling(3));

   BOOST_TEST_MESSAGE("Difference (Yu_FlexibleSUSY - Yu_softsusy)(MZ) = "
                 << ToDoubleMatrix(new_Yu) - softsusy.displayYukawaMatrix(YU));
   BOOST_TEST_MESSAGE("Difference (Yd_FlexibleSUSY - Yd_softsusy)(MZ) = "
                 << ToDoubleMatrix(new_Yd) - softsusy.displayYukawaMatrix(YD));
   BOOST_TEST_MESSAGE("Difference (Ye_FlexibleSUSY - Ye_softsusy)(MZ) = "
                 << ToDoubleMatrix(new_Ye) - softsusy.displayYukawaMatrix(YE));

   model->set_g1(softsusy.displayGaugeCoupling(1));
   model->set_g2(softsusy.displayGaugeCoupling(2));
   model->set_g3(softsusy.displayGaugeCoupling(3));
}

/**
 * @class CMSSMMassWInput_weinberg_angle_low_scale_constraint
 *
 * Replacement class for CMSSMMassWInput_low_scale_constraint, which calculates
 * the gauge couplings at the low scale as Softsusy does it.
 */
class CMSSMMassWInput_weinberg_angle_low_scale_constraint
   : public CMSSMMassWInput_low_scale_constraint<Two_scale> {
public:
   CMSSMMassWInput_weinberg_angle_low_scale_constraint()
      : CMSSMMassWInput_low_scale_constraint<Two_scale>() {}
   CMSSMMassWInput_weinberg_angle_low_scale_constraint(
      CMSSMMassWInput<Two_scale>* model_,
      const CMSSMMassWInput_input_parameters&,
      const QedQcd_legacy& qedqcd_)
      : CMSSMMassWInput_low_scale_constraint<Two_scale>(model_,convert(qedqcd_)) {}
   virtual ~CMSSMMassWInput_weinberg_angle_low_scale_constraint() {}

   virtual void apply();

private:
   void calculate_DRbar_gauge_couplings_gmu();
   void fill_data(Weinberg_angle::Data&, double);
};

void CMSSMMassWInput_weinberg_angle_low_scale_constraint::apply()
{
   assert(model && "Error: CMSSMMassWInput_weinberg_angle_low_scale_constraint:"
          " model pointer must not be zero");

   // save old model parmeters
   const CMSSMMassWInput<Two_scale> mssm(*model);

   model->calculate_DRbar_masses();
   update_scale();
   calculate_DRbar_gauge_couplings();
   calculate_running_SM_masses();

   // save gauge couplings calculated from MW, MZ
   const double g1_mw_mz = new_g1;
   const double g2_mw_mz = new_g2;

   calculate_DRbar_gauge_couplings_gmu();

   // save gauge couplings calculated from GF, MZ
   const double g1_gf_mz = new_g1;
   const double g2_gf_mz = new_g2;

   BOOST_TEST_MESSAGE("Difference (g1_mw_mz - g1_gf_mz)(MZ) = "
                 << g1_mw_mz - g1_gf_mz);
   BOOST_TEST_MESSAGE("Difference (g2_mw_mz - g2_gf_mz)(MZ) = "
                 << g2_mw_mz - g2_gf_mz);

   const auto TanBeta = model->get_input().TanBeta;
   const auto g1 = model->get_g1();
   const auto g2 = model->get_g2();

   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   model->set_vd((2*mZ_run)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
      )));
   model->set_vu((2*mZ_run*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta))));

   // Now calculate the gauge couplings using
   // MssmSoftsusy::sparticleThresholdCorrections
   MssmSoftsusy softsusy;
   copy_parameters(mssm, softsusy);
   softsusy.setData(convert(qedqcd));
   softsusy.setMw(qedqcd.displayPoleMW());

   // prevent tan(beta) from being reset
   softsusy.setSetTbAtMX(true);
   softsusy::PRINTOUT = 10;
   softsusy.sparticleThresholdCorrections(model->get_input().TanBeta);

   if (softsusy.displayProblem().test()) {
      std::ostringstream ss;
      ss << "Softsusy problem in CMSSMMassWInput_precise_gauge_couplings_"
           "low_scale_constraint::apply(): "
           "error while calculating the sparticle thresholds: "
         << softsusy.displayProblem();
      BOOST_TEST_MESSAGE(ss.str());
   }

   // save gauge couplings calculated from GF, MZ with Softsusy
   const double g1_ss = softsusy.displayGaugeCoupling(1);
   const double g2_ss = softsusy.displayGaugeCoupling(2);

   BOOST_TEST_MESSAGE("Difference (g1_gf_mz - g1_ss)(MZ) = "
                 << g1_gf_mz - g1_ss);
   BOOST_TEST_MESSAGE("Difference (g2_gf_mz - g2_ss)(MZ) = "
                 << g2_gf_mz - g2_ss);

   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);
}

void CMSSMMassWInput_weinberg_angle_low_scale_constraint::calculate_DRbar_gauge_couplings_gmu()
{
   assert(qedqcd.get_scale() == get_scale() && "Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");

   const double alpha_em = qedqcd.displayAlpha(ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(ALPHAS);

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds()) {
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);
      delta_alpha_s  = calculate_delta_alpha_s(alpha_s);
   }

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   mZ_run = qedqcd.displayPoleMZ();
   mW_run = qedqcd.displayPoleMW();

   if (model->get_thresholds()) {
      mZ_run = model->calculate_MVZ_DRbar(qedqcd.displayPoleMZ());
      mW_run = model->calculate_MVWm_DRbar(qedqcd.displayPoleMW());
   }

   const double AlphaS = alpha_s_drbar;
   const double EDRbar = e_drbar;

   Weinberg_angle::Data data;
   fill_data(data, alpha_em_drbar);

   Weinberg_angle weinberg;
   weinberg.set_data(data);

   const int error = weinberg.calculate();
   BOOST_REQUIRE(error == 0);
   const double ThetaW = ArcSin(weinberg.get_sin_theta());

   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaW);
   new_g2 = EDRbar*Csc(ThetaW);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

   BOOST_TEST_MESSAGE("ThetaW = " << ThetaW);
   BOOST_TEST_MESSAGE("EDRbar = " << EDRbar);
   BOOST_TEST_MESSAGE("=> g1 = " << new_g1);
   BOOST_TEST_MESSAGE("   g2 = " << new_g2);
}

void CMSSMMassWInput_weinberg_angle_low_scale_constraint::fill_data(
   Weinberg_angle::Data& data, double alpha_em_drbar)
{
   const double scale   = model->get_scale();
   const double mw_pole = qedqcd.displayPoleMW();
   const double mz_pole = qedqcd.displayPoleMZ();
   const double mt_pole = qedqcd.displayPoleMt();
   const double gfermi     = softsusy::GMU;
   const double mt_drbar   = model->get_MFu(2);
   const double mb_drbar   = model->get_MFd(2);
   const double alphaDrbar = alpha_em_drbar;
   const double gY         = model->get_g1() * sqrt(0.6);
   const double g2         = model->get_g2();
   const double hmu        = model->get_Ye(1,1);
   const double g3         = model->get_g3();
   const double mt         = qedqcd.displayPoleMt();
   const double mh         = model->get_Mhh(0);
   const double alpha      = model->get_ZH(0,1);
   const double tanBeta    = model->get_vu() / model->get_vd();
   double mselL            = 0.;
   double msmuL            = 0.;
   double msnue            = 0.;
   double msnumu           = 0.;
   const Eigen::ArrayXd mneut = model->get_MChi();
   const Eigen::MatrixXcd n   = model->get_ZN();
   const Eigen::ArrayXd mch   = model->get_MCha();
   const Eigen::MatrixXcd u   = model->get_UM();
   const Eigen::MatrixXcd v   = model->get_UP();

   const auto MSe(model->get_MSe());
   const auto ZE(model->get_ZE());
   const auto MSv(model->get_MSv());
   const auto ZV(model->get_ZV());

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

   BOOST_REQUIRE(mw_pole > 0.);
   BOOST_REQUIRE(mz_pole > 0.);
   BOOST_CHECK_EQUAL(scale, mz_pole);

   const double pizztMZ = Re(model->self_energy_VZ_1loop(mz_pole));
   const double piwwt0  = Re(model->self_energy_VWm_1loop(0));
   const double piwwtMW = Re(model->self_energy_VWm_1loop(mw_pole));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale = scale;
   se_data.mt_pole = mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY = gY;
   se_data.g2 = g2;

   const double pizztMZ_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole, se_data);

   const double piwwtMW_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwtMW, mw_pole, se_data);

   const double piwwt0_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0., se_data);

   // check corrected self-energies
   {
      model->MFu(2) = mt_pole;
      const double pizztMZ_check = Re(model->self_energy_VZ_1loop(mz_pole));
      const double piwwt0_check  = Re(model->self_energy_VWm_1loop(0));
      const double piwwtMW_check = Re(model->self_energy_VWm_1loop(mw_pole));
      model->MFu(2) = mt_drbar;

      BOOST_CHECK_CLOSE_FRACTION(pizztMZ_check, pizztMZ_corrected, 1.0e-10);
      BOOST_CHECK_CLOSE_FRACTION(piwwtMW_check, piwwtMW_corrected, 1.0e-10);
      BOOST_CHECK_CLOSE_FRACTION(piwwt0_check , piwwt0_corrected , 1.0e-10);
   }

   data.scale = scale;
   data.alpha_em_drbar = alphaDrbar;
   data.fermi_contant = gfermi;
   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0 = piwwt0_corrected;
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

/**
 * @class CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint
 *
 * Replacement class for CMSSMMassWInput_susy_scale_constraint, which does the
 * one-loop ewsb at the susy scale as Softsusy does it.
 */
class CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint
   : public CMSSMMassWInput_susy_scale_constraint<Two_scale> {
public:
   CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint()
      : CMSSMMassWInput_susy_scale_constraint<Two_scale>() {}
   CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint(CMSSMMassWInput<Two_scale>* model_, const CMSSMMassWInput_input_parameters& inputPars_)
      : CMSSMMassWInput_susy_scale_constraint<Two_scale>(model_, QedQcd()) {}
   virtual ~CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint() {}

   virtual void apply();
};

void CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint::apply()
{
   assert(model && "Error: CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint:"
          " model pointer must not be zero");

   // save old model parmeters
   model->calculate_DRbar_masses();
   const CMSSMMassWInput<Two_scale> mssm(*model);

   CMSSMMassWInput_susy_scale_constraint<Two_scale>::apply();

   // Now do the one-loop EWSB using MssmSoftsusy::rewsb
   MssmSoftsusy softsusy;
   softsusy.setAlternativeMs(false);
   copy_parameters(mssm, softsusy);
   softsusy.calcDrBarPars();
   const double new_Msusy = softsusy.calcMs();

   const int signMu = model->get_input().SignMu;
   const double mt = softsusy.displayDrBarPars().mt;
   DoubleVector highScaleSoftPars(3);
   highScaleSoftPars(1) = model->get_input().m0;
   highScaleSoftPars(2) = model->get_input().m12;
   highScaleSoftPars(3) = model->get_input().Azero;

   softsusy.rewsb(signMu, mt, highScaleSoftPars);

   const double new_Mu  = softsusy.displaySusyMu();
   const double new_BMu = softsusy.displayM3Squared();

   BOOST_TEST_MESSAGE("Difference (Mu_FlexibleSUSY - Mu_softsusy)(Msusy) = "
                 << model->get_Mu() - new_Mu);
   BOOST_TEST_MESSAGE("Difference (BMu_FlexibleSUSY - BMu_softsusy)(Msusy) = "
                 << model->get_BMu() - new_BMu);
   BOOST_TEST_MESSAGE("Difference (mt_FlexibleSUSY - mt_softsusy)(Msusy) = "
                 << model->get_MFu()(2) - mt);
   BOOST_TEST_MESSAGE("Difference (Msusy_FlexibleSUSY - Msusy_softsusy)(Msusy) = "
                 << get_scale() - new_Msusy);

   model->set_Mu(new_Mu);
   model->set_BMu(new_BMu);
   scale = new_Msusy;
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
      : mx(0.0), msusy(0.0), softSusy(), gaugeUnification(true) {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   sPhysical get_physical() const { return softSusy.displayPhys(); }
   MssmSoftsusy get_model() const { return softSusy; }
   void test(const CMSSMMassWInput_input_parameters& pp, double mxGuess, const QedQcd_legacy& qedqcd = QedQcd_legacy()) {
      // run softsusy
      softsusy::numRewsbLoops = 1;
      softsusy::numHiggsMassLoops = 1;
      softsusy::TOLERANCE = 1.0e-4;
#ifdef ENABLE_VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      DoubleVector pars(3);
      pars(1) = pp.m0;
      pars(2) = pp.m12;
      pars(3) = pp.Azero;
      softSusy.setAlternativeMs(false);
      softSusy.lowOrg(sugraBcs, mxGuess, pars, pp.SignMu, pp.TanBeta,
                      qedqcd, gaugeUnification);
      mx = softSusy.displayMxBC();
      msusy = softSusy.displayMsusy();
      softsusy::PRINTOUT = 0;

      if (softSusy.displayProblem().test()) {
         std::stringstream ss;
         ss << "SoftSusy problem: " << softSusy.displayProblem();
         VERBOSE_MSG(ss.str());
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
   MssmSoftsusy softSusy;
   bool gaugeUnification;
};

class CMSSMMassWInput_tester {
public:
   CMSSMMassWInput_tester()
      : mx(0.0), msusy(0.0), mssm()
      , ewsb_loop_order(1), pole_mass_loop_order(1)
      , high_constraint(NULL), susy_constraint(NULL), low_constraint(NULL) {}
   ~CMSSMMassWInput_tester() {
      delete high_constraint;
      delete susy_constraint;
      delete low_constraint;
   }
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   CMSSMMassWInput_physical get_physical() const { return mssm.get_physical(); }
   const Problems& get_problems() const { return mssm.get_problems(); }
   CMSSMMassWInput<Two_scale> get_model() const { return mssm; }
   void set_ewsb_loop_order(int l) { ewsb_loop_order = l; }
   void set_pole_mass_loop_order(int l) { pole_mass_loop_order = l; }
   void set_low_scale_constraint(CMSSMMassWInput_low_scale_constraint<Two_scale>* c) { low_constraint = c; }
   void set_susy_scale_constraint(CMSSMMassWInput_susy_scale_constraint<Two_scale>* c) { susy_constraint = c; }
   void set_high_scale_constraint(CMSSMMassWInput_high_scale_constraint<Two_scale>* c) { high_constraint = c; }
   void setup_default_constaints(const CMSSMMassWInput_input_parameters& pp, const QedQcd_legacy& qedqcd) {
      if (!high_constraint)
         high_constraint = new CMSSMMassWInput_high_scale_constraint<Two_scale>(&mssm);
      if (!susy_constraint)
         susy_constraint = new CMSSMMassWInput_susy_scale_constraint<Two_scale>(&mssm, convert(qedqcd));
      if (!low_constraint)
         low_constraint = new CMSSMMassWInput_low_scale_constraint<Two_scale>(&mssm, convert(qedqcd));
   }
   void test(const CMSSMMassWInput_input_parameters& pp, const QedQcd_legacy& qedqcd = QedQcd_legacy()) {
      setup_default_constaints(pp, qedqcd);

      const double precision_goal = softsusy::TOLERANCE;

      mssm.clear();
      mssm.set_loops(2);
      mssm.set_thresholds(2);
      mssm.set_ewsb_loop_order(ewsb_loop_order);
      mssm.set_pole_mass_loop_order(pole_mass_loop_order);
      mssm.set_input_parameters(pp);
      mssm.set_precision(precision_goal);

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

      CMSSMMassWInput_convergence_tester<Two_scale> convergence_tester(&mssm, precision_goal);
      CMSSMMassWInput_initial_guesser<Two_scale> initial_guesser(&mssm, convert(qedqcd),
                                                      *low_constraint,
                                                      *susy_constraint,
                                                      *high_constraint);
      Two_scale_increasing_precision precision(10.0, precision_goal);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&convergence_tester);
      solver.set_running_precision(&precision);
      solver.set_initial_guesser(&initial_guesser);
      solver.add(low_constraint, &mssm);
      solver.add(high_constraint, &mssm);
      solver.add(susy_constraint, &mssm);
      solver.solve();
      mssm.run_to(susy_constraint->get_scale());
      mssm.solve_ewsb();
      mssm.calculate_spectrum();
      mssm.run_to(Electroweak_constants::MZ);

      mx = high_constraint->get_scale();
      msusy = susy_constraint->get_scale();
   }
private:
   double mx, msusy;
   CMSSMMassWInput<Two_scale> mssm;
   CMSSMMassWInput_high_scale_constraint<Two_scale>* high_constraint;
   CMSSMMassWInput_susy_scale_constraint<Two_scale>* susy_constraint;
   CMSSMMassWInput_low_scale_constraint<Two_scale>*  low_constraint;
   int ewsb_loop_order, pole_mass_loop_order;
};

BOOST_AUTO_TEST_CASE( test_CMSSMMassWInput_spectrum )
{
   CMSSMMassWInput_input_parameters pp;
   pp.m0 = 125.;
   pp.m12 = 500.;
   pp.TanBeta = 10.;
   pp.SignMu = 1;
   pp.Azero = 0.;

   CMSSMMassWInput<Two_scale> _model(pp);
   const CMSSMMassWInput_high_scale_constraint<Two_scale> high_constraint(&_model);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   CMSSMMassWInput_tester mssm_tester;
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 0.13);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.006);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const CMSSMMassWInput<Two_scale> fs(mssm_tester.get_model());

   BOOST_CHECK_EQUAL(ss.displayLoops()     , fs.get_loops());
   BOOST_CHECK_EQUAL(ss.displayMu()        , fs.get_scale());
   // BOOST_CHECK_EQUAL(ss.displayThresholds(), fs.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00076);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0011);
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

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassB() , ss.displayGaugino(1), 0.0044);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassWB(), ss.displayGaugino(2), 0.0046);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassG() , ss.displayGaugino(3), 0.0051);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.0043);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.009);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.0017);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.0077);

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

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-7);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0093);

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
   const double mA0 = ss.displayDrBarPars().mA0(1);
   const double mh0 = ss.displayDrBarPars().mh0(1);
   const double mH0 = ss.displayDrBarPars().mh0(2);

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha(1), mch(1), 0.0013);
   BOOST_CHECK_CLOSE_FRACTION(MCha(2), mch(2), 0.008);

   BOOST_CHECK_CLOSE_FRACTION(MChi(1), mn(1), 0.0042);
   BOOST_CHECK_CLOSE_FRACTION(MChi(2), mn(2), 0.0041);
   BOOST_CHECK_CLOSE_FRACTION(MChi(3), mn(3), 0.0043);
   BOOST_CHECK_CLOSE_FRACTION(MChi(4), mn(4), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0  , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0, 0.00007);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mH0, 0.004);

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
   BOOST_CHECK_CLOSE_FRACTION(Su(1), mu(1), 0.006);
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
   const DoubleVector msnu(ss.displayDrBarPars().msnu);
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

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFe()(2), ss.displayDrBarPars().mtau, 0.0011);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFu()(2), ss.displayDrBarPars().mt  , 0.0097);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFd()(2), ss.displayDrBarPars().mb  , 0.0027);


   // comparing one-loop masses

   const DoubleVector
      MCha_1l(ToDoubleVector(fs.get_physical().MCha)),
      MChi_1l(ToDoubleVector(fs.get_physical().MChi)),
      MHpm_1l(ToDoubleVector(fs.get_physical().MHpm)),
      MAh_1l(ToDoubleVector(fs.get_physical().MAh)),
      Mhh_1l(ToDoubleVector(fs.get_physical().Mhh));
   const DoubleVector mch_1l(ss.displayPhys().mch),
      mn_1l(ss.displayPhys().mneut.apply(fabs));
   const double mHpm_1l = ss.displayPhys().mHpm;
   const double mA0_1l  = ss.displayPhys().mA0(1);
   const double mh0_1l  = ss.displayPhys().mh0(1);
   const double mH0_1l  = ss.displayPhys().mh0(2);

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha_1l(1), mch_1l(1), 0.0011);
   BOOST_CHECK_CLOSE_FRACTION(MCha_1l(2), mch_1l(2), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(1), mn_1l(1), 0.0042);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(2), mn_1l(2), 0.0011);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(3), mn_1l(3), 0.0043);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(4), mn_1l(4), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MHpm_1l(2), mHpm_1l , 0.004);
   BOOST_CHECK_CLOSE_FRACTION(MAh_1l(2) , mA0_1l  , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(1), mh0_1l, 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(2), mH0_1l, 0.004);

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
   const DoubleVector msnu_1l(ss.displayPhys().msnu);
   const DoubleVector Snu_1l(ToDoubleVector(fs.get_physical().MSv));
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(1), msnu_1l(1), 0.0056);
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(2), msnu_1l(2), 0.0056);
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(3), msnu_1l(3), 0.0090);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_physical().MGlu, ss.displayPhys().mGluino, 0.005);
}

// ===== test with gauge couplings determined from the Rho parameter (SoftSUSY variant) =====

BOOST_AUTO_TEST_CASE( test_CMSSMMassWInput_spectrum_with_Softsusy_gauge_couplings )
{
   CMSSMMassWInput_input_parameters pp;
   pp.m0 = 125.;
   pp.m12 = 500.;
   pp.TanBeta = 10.;
   pp.SignMu = 1;
   pp.Azero = 0.;
   softsusy::QedQcd_legacy qedqcd;

   CMSSMMassWInput<Two_scale> _model(pp);
   const CMSSMMassWInput_high_scale_constraint<Two_scale> high_constraint(&_model);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   CMSSMMassWInput_tester mssm_tester;
   mssm_tester.set_low_scale_constraint(new CMSSMMassWInput_precise_gauge_couplings_low_scale_constraint(&_model, pp, qedqcd));
   mssm_tester.set_susy_scale_constraint(new CMSSMMassWInput_softsusy_ewsb_susy_scale_constraint(&_model, pp));
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 0.0007);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 2.0e-5);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const CMSSMMassWInput<Two_scale> fs(mssm_tester.get_model());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.0000023);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0000066);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.0000010);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.0012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.0024);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.0022);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-7);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0068);

   // comparing tree-level masses

   const DoubleVector MHpm(ToDoubleVector(fs.get_MHpm())),
      MAh(ToDoubleVector(fs.get_MAh())),
      Mhh(ToDoubleVector(fs.get_Mhh()));
   const double MwRun = fs.get_MVWm();
   const double MzRun = fs.get_MVZ();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const double mA0 = ss.displayDrBarPars().mA0(1);
   const double mh0 = ss.displayDrBarPars().mh0(1);
   const double mH0 = ss.displayDrBarPars().mh0(2);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm, 8.e-5);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 1.0e-10); // for RXi(VZ) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0, 8.e-5);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0, 4.e-6);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mH0, 8.e-5);
}

// ===== test with gauge couplings determined from the Rho parameter (FlexibleSUSY variant) =====

BOOST_AUTO_TEST_CASE( test_CMSSMMassWInput_spectrum_with_weinberg_angle )
{
   CMSSMMassWInput_input_parameters pp;
   pp.m0 = 125.;
   pp.m12 = 500.;
   pp.TanBeta = 10.;
   pp.SignMu = 1;
   pp.Azero = 0.;
   softsusy::QedQcd_legacy qedqcd;

   CMSSMMassWInput<Two_scale> _model(pp);
   const CMSSMMassWInput_high_scale_constraint<Two_scale> high_constraint(&_model);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   CMSSMMassWInput_tester mssm_tester;
   mssm_tester.set_low_scale_constraint(new CMSSMMassWInput_weinberg_angle_low_scale_constraint(&_model, pp, qedqcd));
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 6.2e-4);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 1.5e-5);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const CMSSMMassWInput<Two_scale> fs(mssm_tester.get_model());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.0000023);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0000066);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.0000010);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.0012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.0024);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.0022);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-7);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0068);

   // comparing tree-level masses

   const DoubleVector MHpm(ToDoubleVector(fs.get_MHpm())),
      MAh(ToDoubleVector(fs.get_MAh())),
      Mhh(ToDoubleVector(fs.get_Mhh()));
   const double MwRun = fs.get_MVWm();
   const double MzRun = fs.get_MVZ();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const double mA0 = ss.displayDrBarPars().mA0(1);
   const double mh0 = ss.displayDrBarPars().mh0(1);
   const double mH0 = ss.displayDrBarPars().mh0(2);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm, 8.e-5);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 1.0e-10); // for RXi(VZ) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0, 8.e-5);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0, 3.e-6);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mH0, 8.e-5);
}


/**
 * @class CMSSMMassWInput_iterative_low_scale_constraint
 *
 * Replacement class for CMSSMMassWInput_low_scale_constraint, which inputs the
 * Higgs mass and eliminates TanBeta.
 */
class CMSSMMassWInput_iterative_low_scale_constraint
   : public CMSSMMassWInput_low_scale_constraint<Two_scale> {
public:
   CMSSMMassWInput_iterative_low_scale_constraint()
      : CMSSMMassWInput_low_scale_constraint<Two_scale>() {}
   CMSSMMassWInput_iterative_low_scale_constraint(CMSSMMassWInput<Two_scale>* model_, const CMSSMMassWInput_input_parameters&, const QedQcd_legacy& qedqcd_)
      : CMSSMMassWInput_low_scale_constraint<Two_scale>(model_,convert(qedqcd_)) {}
   virtual ~CMSSMMassWInput_iterative_low_scale_constraint() {}

   virtual void apply();
};

void CMSSMMassWInput_iterative_low_scale_constraint::apply()
{
   assert(model && "Error: CMSSMMassWInput_low_scale_constraint:"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();
   calculate_DRbar_gauge_couplings();
   calculate_DRbar_yukawa_couplings();

   auto func = [this](const Eigen::Matrix<double,2,1>& x) {
      const double vd = x(0);
      const double vu = x(1);

      if (vd < std::numeric_limits<double>::epsilon() ||
          vu < std::numeric_limits<double>::epsilon())
         return std::numeric_limits<double>::max();

      model->set_vd(vd);
      model->set_vu(vu);

      model->calculate_DRbar_masses();
      model->calculate_Mhh_pole();
      model->calculate_MVZ_pole();

      const double mH = model->get_physical().Mhh(0);
      const double mZ = model->get_physical().MVZ;

      #define LowEnergyConstant(p) Electroweak_constants::p
      #define STANDARD_DEVIATION(p) Electroweak_constants::Error_##p

      return Sqr(LowEnergyConstant(MZ) - mZ)/Sqr(STANDARD_DEVIATION(MZ))
         + Sqr(LowEnergyConstant(MH) - mH)/Sqr(STANDARD_DEVIATION(MH)*10);
   };

   Minimizer<2> minimizer(func, 100, 1.0e-2);
   Eigen::Matrix<double,2,1> start;
   start << model->get_vd(), model->get_vu();

   const int status = minimizer.minimize(start);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_TEST_MESSAGE("chi^2 = " << minimizer.get_minimum_value());
   BOOST_TEST_MESSAGE("New vd = " << model->get_vd() << ", vu = " << model->get_vu());
   BOOST_TEST_MESSAGE("Predicted tan(beta) = " << model->get_vu() / model->get_vd());

   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);
}

BOOST_AUTO_TEST_CASE( test_CMSSMMassWInput_spectrum_higgs_iteration )
{
   CMSSMMassWInput_input_parameters pp;
   pp.m0 = 500.;
   pp.m12 = 500.;
   pp.Azero = 1000.;
   pp.SignMu = 1;
   pp.TanBeta = 30.;
   softsusy::QedQcd_legacy qedqcd;

   CMSSMMassWInput<Two_scale> _model(pp);
   CMSSMMassWInput_tester mssm_tester;
   mssm_tester.set_low_scale_constraint(new CMSSMMassWInput_iterative_low_scale_constraint(&_model, pp, qedqcd));
   // BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   try {
      mssm_tester.test(pp);
   } catch (Error& error) {
      ERROR(error.what());
   }
}

BOOST_AUTO_TEST_CASE( test_CMSSMMassWInput_EWSB_problems )
{
   CMSSMMassWInput_input_parameters pp;
   pp.m0 = 5000.;
   pp.m12 = 5000.;
   pp.Azero = 5000.;
   pp.SignMu = 1;
   pp.TanBeta = 22.45;

   CMSSMMassWInput_tester mssm_tester;
   mssm_tester.set_ewsb_loop_order(2);
   mssm_tester.set_pole_mass_loop_order(2);
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   if (mssm_tester.get_problems().no_ewsb()) {
      BOOST_ERROR("no ewsb");
      BOOST_TEST_MESSAGE("\tProblems: " << mssm_tester.get_problems());
   }
}
