
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "wrappers.hpp"
#include "pv.hpp"
#include "SM_two_scale_model.hpp"
#include "lowe.h"

using namespace flexiblesusy;
using namespace softsusy;
using namespace passarino_veltman;


BOOST_AUTO_TEST_CASE( test_SM_four_loop_top_pole_mass )
{
   const QedQcd qedqcd;
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   Loop_corrections tlc;
   tlc.top_qcd = 3;               // add 4L QCD correction

   SM<Two_scale> m;
   m.do_calculate_sm_pole_masses(true);
   m.set_thresholds(4);           // add 4L QCD correction
   m.set_ewsb_loop_order(4);      // add 4L QCD correction
   m.set_pole_mass_loop_order(4); // add 4L QCD correction
   m.set_loop_corrections(tlc);

   setup_SM_const(m, input);

   const double mt_pole_input = qedqcd.displayPoleMt();
   const double v = m.get_v();

   BOOST_TEST_MESSAGE("mt_pole(input) = " << mt_pole_input);

   // calculate DR-bar masses
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.solve_ewsb();

   BOOST_TEST_MESSAGE("mt_MSbar(guess) = " << m.get_MFu(2));

   int iterations = 100;

   // calculate top DR-bar mass from top pole mass using three-loop
   // corrections
   do {
      Eigen::Matrix<double,3,3> mt_MSbar_3loop(Eigen::Matrix<double,3,3>::Zero());
      mt_MSbar_3loop(2,2) = m.calculate_MFu_DRbar(mt_pole_input, 2);
      m.set_Yu(((1.4142135623730951*mt_MSbar_3loop)/v).transpose());
      m.calculate_DRbar_masses();
      m.solve_ewsb();
   } while (--iterations);

   BOOST_TEST_MESSAGE("mt_MSbar(4-loop) = " << m.get_MFu(2));

   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double mt_pole_3loop  = m.get_physical().MFu(2);

   BOOST_TEST_MESSAGE("mt_pole(4-loop) = " << mt_pole_3loop);

   BOOST_CHECK_CLOSE_FRACTION(mt_pole_input, mt_pole_3loop, 2.0e-4);
}
