
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_loopfunctions

#include <boost/test/unit_test.hpp>
#include "wrappers.hpp"

#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_two_scale_model.hpp"

using namespace flexiblesusy;

namespace
{

/// gives access to the protected loop functions of the model
struct Test_model : CMSSM<Two_scale> {
   using CMSSM_mass_eigenstates::A0;
   using CMSSM_mass_eigenstates::B0;
   using CMSSM_mass_eigenstates::B00;
   using CMSSM_mass_eigenstates::B1;
   using CMSSM_mass_eigenstates::B22;
   using CMSSM_mass_eigenstates::F0;
   using CMSSM_mass_eigenstates::G0;
   using CMSSM_mass_eigenstates::H0;
};

} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_A0 )
{
   Test_model model;
   model.set_scale(100.);

   BOOST_CHECK_EQUAL(model.A0(0.), 0.);
}

BOOST_AUTO_TEST_CASE( test_B0 )
{
   Test_model model;
   const double scale = 100.;
   const double scale2 = Sqr(scale);
   const double p = 91.0;
   const double p2 = Sqr(p);
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.B0(0., 0., 0.), 0.);

   BOOST_CHECK_CLOSE(model.B0(p2 , 0., 0.), 2. - log(p2/scale2), 1.0e-10);
   BOOST_CHECK_CLOSE(model.B0(0. , p2, 0.), 1. - log(p2/scale2), 1.0e-10);
   BOOST_CHECK_EQUAL(model.B0(0. , 0., p2), model.B0(0., p2, 0.));

   BOOST_CHECK_CLOSE(model.B0(p2, p2, 0.), 2. - log(p2/scale2), 0.005);
   BOOST_CHECK_CLOSE(model.B0(0., p2, p2), 0. - log(p2/scale2), 0.005);
   BOOST_CHECK_EQUAL(model.B0(p2, 0., p2), model.B0(p2 , p2 , 0.));
}

BOOST_AUTO_TEST_CASE( test_B1 )
{
   Test_model model;
   const double scale = 100.;
   const double p = 91.0;
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.B1(0., 0., 0.), 0.);

   BOOST_CHECK_CLOSE(model.B1(p , 0., 0.), -0.5 * model.B0(p, 0., 0.), 1.0e-10);
   BOOST_CHECK_CLOSE(model.B1(0., p , p ), -0.5 * model.B0(0., p, p ), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_B00 )
{
   Test_model model;
   const double scale = 100.;
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.B00(0., 0., 0.), 0.);
}

BOOST_AUTO_TEST_CASE( test_B22 )
{
   Test_model model;
   const double scale = 100.;
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.B22(0., 0., 0.), 0.);
}

BOOST_AUTO_TEST_CASE( test_H0 )
{
   Test_model model;
   const double scale = 100.;
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.H0(0., 0., 0.), 0.);
}

BOOST_AUTO_TEST_CASE( test_F0 )
{
   Test_model model;
   const double scale = 100.;
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.F0(0., 0., 0.), 0.);
}

BOOST_AUTO_TEST_CASE( test_G0 )
{
   Test_model model;
   const double scale = 100.;
   model.set_scale(scale);

   BOOST_CHECK_EQUAL(model.G0(0., 0., 0.), 0.);
}
