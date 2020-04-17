#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_twoloop_as

#include <boost/test/unit_test.hpp>
#include "mssm_twoloop_as.hpp"
#include <cmath>

using namespace flexiblesusy;
using namespace flexiblesusy::mssm_twoloop_as;

struct Results {
   Parameters pars;
   Real das_as_as{};
   Real das_at_as{};
   Real das_ab_as{};
};

Parameters make_point(Real MS, Real tb, Real xt)
{
   Parameters pars;

   pars.g3   = 1.2l;
   pars.yt   = 0.9l;
   pars.yb   = 0.1l;
   pars.mt   = 173.34l;
   pars.mb   = 2.4l;
   pars.mg   = MS * 0.9;
   pars.mst1 = MS;
   pars.mst2 = MS * 1.1;
   pars.msb1 = MS * 1.2;
   pars.msb2 = MS * 1.3;
   pars.msd1 = MS * 1.4;
   pars.msd2 = MS * 1.5;
   pars.xt   = xt;
   pars.xb   = xt * 1.1;
   pars.mw   = 83.0l;
   pars.mz   = 91.0l;
   pars.mh   = 125.0l;
   pars.mH   = MS * 1.6;
   pars.mC   = MS * 1.7;
   pars.mA   = MS * 1.8;
   pars.mu   = MS * 1.9;
   pars.tb   = tb;
   pars.Q    = MS * 2.0;

   return pars;
}

BOOST_AUTO_TEST_CASE(test_points)
{
   std::vector<Results> res = {
      { make_point(1000.0l, 10.0l, 2.0l), 0.000715816532730145103587l, 0.000108116526832918386358l, -2.30035474310885696247e-05l }
   };

   for (const auto r: res) {
      const Real das_as_as = delta_alpha_s_2loop_as_as(r.pars);
      const Real das_at_as = delta_alpha_s_2loop_at_as(r.pars);
      const Real das_ab_as = delta_alpha_s_2loop_ab_as(r.pars);

      BOOST_CHECK_CLOSE_FRACTION(das_as_as, r.das_as_as, 1.0e-10l);
      BOOST_CHECK_CLOSE_FRACTION(das_at_as, r.das_at_as, 1.0e-10l);
      BOOST_CHECK_CLOSE_FRACTION(das_ab_as, r.das_ab_as, 1.0e-10l);
   }
}
