#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_error

#include <boost/test/unit_test.hpp>
#include <type_traits>
#include "error.hpp"

#define THROW(Exception_t,...)                  \
   do {                                         \
      throw Exception_t(__VA_ARGS__);           \
   } while (false)

BOOST_AUTO_TEST_CASE( test_catch_by_std_exception )
{
   using namespace flexiblesusy;

   using ET = std::exception;

   BOOST_REQUIRE_THROW(THROW(FatalError)    , ET);
   BOOST_REQUIRE_THROW(THROW(SetupError, ""), ET);
   BOOST_REQUIRE_THROW(THROW(NoConvergenceError, 1, ""), ET);
   BOOST_REQUIRE_THROW(THROW(NoSinThetaWConvergenceError, 1, 0.2), ET);
   BOOST_REQUIRE_THROW(THROW(NonPerturbativeSinThetaW), ET);
   BOOST_REQUIRE_THROW(THROW(NonPerturbativeRunningError, 0.), ET);
   BOOST_REQUIRE_THROW(THROW(RGIntegrationError,
                             RGIntegrationError::Reason::too_many_steps, 0., 0,
                             0.),
                       ET);
   BOOST_REQUIRE_THROW(THROW(NonPerturbativeRunningQedQcdError, ""), ET);
   BOOST_REQUIRE_THROW(THROW(OutOfMemoryError, ""), ET);
   BOOST_REQUIRE_THROW(THROW(OutOfBoundsError, ""), ET);
   BOOST_REQUIRE_THROW(THROW(ReadError, ""), ET);
   BOOST_REQUIRE_THROW(THROW(PhysicalError, ""), ET);
   BOOST_REQUIRE_THROW(THROW(HimalayaError, ""), ET);
}

/*
BOOST_AUTO_TEST_CASE( test_is_nothrow_constructible )
{
   using namespace flexiblesusy;

   BOOST_CHECK(std::is_nothrow_constructible_v<FatalError>);
}
*/

BOOST_AUTO_TEST_CASE(test_non_perturbative_running_error_message)
{
   using namespace flexiblesusy;

   const NonPerturbativeRunningError e(123.);
   BOOST_TEST_MESSAGE(e.what_detailed());
   BOOST_CHECK(e.what_detailed().find("123") != std::string::npos);

   const RGIntegrationError r(RGIntegrationError::Reason::too_many_steps, 117.,
                              3, 1e-14);
   BOOST_TEST_MESSAGE(r.what_detailed());
   BOOST_CHECK(r.what_detailed().find("too many steps") != std::string::npos);
   BOOST_CHECK(r.get_reason() == RGIntegrationError::Reason::too_many_steps);
}
