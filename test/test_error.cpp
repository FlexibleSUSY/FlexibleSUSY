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

   BOOST_CHECK(std::is_nothrow_constructible<FatalError>::value);
}
*/
