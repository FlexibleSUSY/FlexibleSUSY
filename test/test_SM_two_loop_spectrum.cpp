
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "loop_libraries/loop_library.hpp"
#include "test_SM.hpp"
#include "wrappers.hpp"
#include "SM_two_scale_model.hpp"
#include "lowe.h"

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_SM_two_loop_top_pole_mass )
{
   Loop_library::set(0);
   const QedQcd qedqcd;
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   m.do_calculate_sm_pole_masses(true);
   m.set_thresholds(2);
   setup_SM_const(m, input);

   const double mt_pole_input = qedqcd.displayPoleMt();
   const double v = m.get_v();

   BOOST_TEST_MESSAGE("mt_pole(input) = " << mt_pole_input);

   // calculate DR-bar masses
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.solve_ewsb();

   BOOST_TEST_MESSAGE("mt_drbar(guess) = " << m.get_MFu(2));

   int iterations = 100;

   // calculate top DR-bar mass from top pole mass using two-loop
   // corrections
   do {
      Eigen::Matrix<double,3,3> mt_drbar_2loop(Eigen::Matrix<double,3,3>::Zero());
      mt_drbar_2loop(2,2) = m.calculate_MFu_DRbar(mt_pole_input, 2);
      m.set_Yu(((1.4142135623730951*mt_drbar_2loop)/v).transpose());
      m.calculate_DRbar_masses();
      m.solve_ewsb();
   } while (--iterations);

   BOOST_TEST_MESSAGE("mt_drbar(2-loop) = " << m.get_MFu(2));

   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double mt_pole_2loop  = m.get_physical().MFu(2);

   BOOST_TEST_MESSAGE("mt_pole(2-loop) = " << mt_pole_2loop);

   BOOST_CHECK_CLOSE_FRACTION(mt_pole_input, mt_pole_2loop, 3.0e-4);
}
