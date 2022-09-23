
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_higgs_iteration

#include <boost/test/unit_test.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include "minimizer.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_input_parameters.hpp"
#include "CMSSM_two_scale_high_scale_constraint.hpp"
#include "CMSSM_two_scale_susy_scale_constraint.hpp"
#include "CMSSM_two_scale_low_scale_constraint.hpp"
#include "CMSSM_two_scale_convergence_tester.hpp"
#include "CMSSM_two_scale_initial_guesser.hpp"
#include "test_CMSSM.hpp"

#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARD_DEVIATION(p) Electroweak_constants::Error_##p

BOOST_AUTO_TEST_CASE( test_copy_Minimizer )
{
   Minimizer<2> minimizer1(nullptr, 100, 1.0e-2);
   Minimizer<2> minimizer2(minimizer1);
}

BOOST_AUTO_TEST_CASE( test_CMSSM_higgs_iteration )
{
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;
   CMSSM<Two_scale> model(input);

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   model.set_scale(scale);
   model.set_loops(1);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   auto func = [&model](const Eigen::Matrix<double,2,1>& x) -> double {
      const double vd = x(0);
      const double vu = x(1);

      model.set_vd(vd);
      model.set_vu(vu);

      model.calculate_DRbar_masses();
      model.calculate_Mhh_pole();
      model.calculate_MVZ_pole();

      const double mH = model.get_physical().Mhh(1);
      const double mZ = model.get_physical().MVZ;

      return Sqr(LowEnergyConstant(MZ) - mZ)/Sqr(STANDARD_DEVIATION(MZ))
         + Sqr(LowEnergyConstant(MH) - mH)/Sqr(STANDARD_DEVIATION(MH));
   };

   Minimizer<2> minimizer(func, 100, 1.0e-2);
   Eigen::Matrix<double,2,1> start;
   start << model.get_vd(), model.get_vu();

   const int status = minimizer.minimize(start);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_TEST_MESSAGE("New vd = " << model.get_vd() << ", vu = " << model.get_vu());
   BOOST_TEST_MESSAGE("Predicted tan(beta) = " << model.get_vu() / model.get_vd());

   // check how close we got
   model.calculate_DRbar_masses();
   model.calculate_Mhh_pole();
   model.calculate_MVZ_pole();

   // const double mH = model.get_physical().Mhh(1);
   const double mZ = model.get_physical().MVZ;

   // BOOST_CHECK_CLOSE_FRACTION(mH, 125., 0.400);
   BOOST_CHECK_CLOSE_FRACTION(mZ, 91. , 0.003);
}
