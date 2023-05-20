
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_two_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "wrappers.hpp"
#include "ew_input.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "lowe.h"

using namespace flexiblesusy;
using namespace softsusy;

void ensure_tree_level_ewsb(CMSSM<Two_scale>& m)
{
   // ensure that the EWSB eqs. are satisfied (Drees p.222)
   const double vu = m.get_vu();
   const double vd = m.get_vd();
   const double gY = m.get_g1() * Sqrt(0.6);
   const double g2 = m.get_g2();
   const double Mu = m.get_Mu();
   const double BMu = m.get_BMu();
   const double mHd2 = BMu*vu/vd - (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   const double mHu2 = BMu*vd/vu + (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   m.set_mHd2(mHd2);
   m.set_mHu2(mHu2);
}

template <class T1, class TInput>
void setup_CMSSM_const(T1& m, const TInput& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4 * Pi * alpha1);
   const double g2 = Sqrt(4 * Pi * alpha2);
   const double g3 = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = Sqrt(2.0);
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

   m.set_scale(scale);
   m.set_loops(2);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
}

BOOST_AUTO_TEST_CASE( test_CMSSM_two_loop_top_pole_mass )
{
   const QedQcd qedqcd;
   CMSSM_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.SignMu = 1;
   input.Azero = 0.;
   CMSSM<Two_scale> m(input);
   m.set_thresholds(2);
   m.do_calculate_sm_pole_masses(true);
   setup_CMSSM_const(m, input);

   const double mt_pole_input = qedqcd.displayPoleMt();
   const double vu = m.get_vu();

   BOOST_TEST_MESSAGE("mt_pole(input) = " << mt_pole_input);

   // calculate DR-bar masses
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.solve_ewsb_tree_level();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   BOOST_TEST_MESSAGE("mt_drbar(guess) = " << m.get_MFu(2));

   int iterations = 100;

   // calculate top DR-bar mass from top pole mass using two-loop
   // corrections
   do {
      Eigen::Matrix<double,3,3> mt_drbar_2loop(Eigen::Matrix<double,3,3>::Zero());
      mt_drbar_2loop(0,0) = qedqcd.displayMass(mUp);
      mt_drbar_2loop(1,1) = qedqcd.displayMass(mCharm);
      mt_drbar_2loop(2,2) = m.calculate_MFu_DRbar(mt_pole_input, 2);
      m.set_Yu(((1.4142135623730951*mt_drbar_2loop)/vu).transpose());

      m.calculate_DRbar_masses();
      m.solve_ewsb_tree_level();

      if (m.get_problems().have_problem()) {
         std::ostringstream ostr;
         m.get_problems().print_problems(ostr);
         BOOST_FAIL(ostr.str());
      }
   } while (--iterations);

   BOOST_TEST_MESSAGE("mt_drbar(2-loop) = " << m.get_MFu(2));

   m.set_pole_mass_loop_order(2);
   m.calculate_MFu_pole();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double mt_pole_2loop  = m.get_physical().MFu(2);

   BOOST_TEST_MESSAGE("mt_pole(2-loop) = " << mt_pole_2loop);

   BOOST_CHECK_CLOSE_FRACTION(mt_pole_input, mt_pole_2loop, 1.4e-3);
}
