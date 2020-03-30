
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_composite_root_finder

#include <boost/test/unit_test.hpp>

#include "composite_root_finder.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

typedef Eigen::Matrix<double,2,1> EV2_t;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   auto parabola = [](const EV2_t& x) -> EV2_t {
      const double y = x(0);
      const double z = x(1);
      EV2_t f;
      f << y*(y - 5.0), z*(z - 1.0);
      return f;
   };

   Eigen::Matrix<double,2,1> start;
   start << 10, 10;

   Composite_root_finder<2> ss;
   ss.set_precision(1e-5);

   BOOST_REQUIRE_NO_THROW(ss.solve(parabola, start));

   const auto sol = ss.get_solution();

   BOOST_CHECK_CLOSE_FRACTION(sol(0), 5.0, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(sol(1), 1.0, 1.0e-5);
}
