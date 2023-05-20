
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_susy_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_CMSSM.hpp"

#include <cmath>

#define private public

#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_two_scale_susy_scale_constraint.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"

BOOST_AUTO_TEST_CASE( test_susy_scale_constraint )
{
   QedQcd qedqcd;
   CMSSM_input_parameters input;
   CMSSM<Two_scale> m; MssmSoftsusy s;
   setup_CMSSM(m, s, input);

   CMSSM_susy_scale_constraint<Two_scale> constraint(&m, qedqcd);
   constraint.apply();

   const auto tadpole = m.tadpole_equations();

   // check that EWSB eqs. are fulfilled
   BOOST_CHECK_SMALL(std::fabs(tadpole[0]), 1.0);
   BOOST_CHECK_SMALL(std::fabs(tadpole[1]), 2.3);
}
