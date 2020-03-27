
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMNoFV_onshell

#include <boost/test/unit_test.hpp>

#include "test_CMSSMNoFV.hpp"
#include "wrappers.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2_1loop_helpers.hpp"
#include "gm2_2loop_helpers.hpp"
#include "test_MSSMNoFV_onshell.hpp"

using namespace flexiblesusy;
using namespace gm2calc;

BOOST_AUTO_TEST_CASE( test_gm2_standard_point )
{
   gm2calc::MSSMNoFV_onshell osmodel(setup());

   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MM()    ,  0.103936   , 3e-06);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MT()    ,  173.3      , 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MW()    ,  80.404     , 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MZ()    ,  91.1876    , 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_TB()    ,  10         , 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Yu(0,0) , 1.41242e-05 , 4e-06);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Yu(1,1) , 0.0074682   , 7e-7);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Yu(2,2) , 1.01762     , 3e-6);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Yd(0,0) , 0.000279538 , 5e-7);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Yd(1,1) , 0.00612041  , 5e-7);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Yd(2,2) , 0.156175706 , 5e-7);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_BMu()   , 49493.5     , 5e-7);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_Mu()    , 619.858     , 6e-7);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MassB() , 211.722     , 3e-6);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MassWB(), 401.057     , 5e-7);
   BOOST_CHECK_CLOSE_FRACTION(std::sqrt(osmodel.get_ml2(1,1)), 356.042, 2e-6);
   BOOST_CHECK_CLOSE_FRACTION(std::sqrt(osmodel.get_me2(1,1)), 225.075, 3e-6);

   const double gm2_1l = gm2calc::calculate_amu_1loop_non_tan_beta_resummed(osmodel);
   const double gm2_2l_tanb_approx =  + (gm2calc::tan_beta_cor(osmodel) - 1.) * gm2_1l;

   BOOST_CHECK_CLOSE_FRACTION(gm2_1l               ,  8.91837e-10, 4e-7);
   // BOOST_CHECK_CLOSE_FRACTION(amu1Lapprox(osmodel),  9.04621e-10, 3e-7);
   BOOST_CHECK_CLOSE_FRACTION(amu1LWHnu(osmodel)   ,  7.59385e-10, 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(amu1LBmuLmuR(osmodel),  4.32919e-10, 3e-6);
   BOOST_CHECK_CLOSE_FRACTION(amu1LBHmuL(osmodel)  ,  4.51653e-11, 3e-6);
   BOOST_CHECK_CLOSE_FRACTION(amu1LWHmuL(osmodel)  , -1.58811e-10, 3e-6);
   BOOST_CHECK_CLOSE_FRACTION(amu1LBHmuR(osmodel)  , -1.73907e-10, 3e-6);
   BOOST_CHECK_CLOSE_FRACTION(gm2_2l_tanb_approx   ,  1.80958e-11, 3e-6);

}
