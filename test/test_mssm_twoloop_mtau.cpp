#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_twoloop_mtau

#include <boost/test/unit_test.hpp>
#include "mssm_twoloop_mtau.hpp"
#include <cmath>

using namespace flexiblesusy;
using namespace flexiblesusy::mssm_twoloop_mtau;

struct Results {
   Parameters pars;
   Real dmtau_atau_atau{};
   Real dmtau_atau_at{};
   Real dmtau_atau_ab{};
};

Parameters make_point(Real MS, Real tb, Real xt)
{
   Parameters pars;

   pars.yt   = 0.9l;
   pars.yb   = 0.1l;
   pars.ytau = 0.05l;
   pars.mt   = 173.34l;
   pars.mb   = 2.4l;
   pars.mtau = 1.7l;
   pars.mst1 = MS;
   pars.mst2 = MS * 1.1;
   pars.msb1 = MS * 1.2;
   pars.msb2 = MS * 1.3;
   pars.mstau1 = MS * 1.4;
   pars.mstau2 = MS * 1.5;
   pars.msntau = MS * 1.6;
   pars.xt   = xt;
   pars.xb   = xt * 1.1;
   pars.xtau = xt * 1.2;
   pars.mw   = 83.0l;
   pars.mz   = 91.0l;
   pars.mh   = 125.0l;
   pars.mH   = MS * 1.7;
   pars.mC   = MS * 1.8;
   pars.mA   = MS * 1.9;
   pars.mu   = MS * 2.0;
   pars.tb   = tb;
   pars.Q    = MS * 2.1;

   return pars;
}

BOOST_AUTO_TEST_CASE(test_points)
{
   std::vector<Results> res = {
      { make_point(1000.0l, 10.0l, 2.0l), -7.2740436727506185726e-09, -1.47266541588342011481e-07, -1.71438422230067884779e-07 }
   };

   for (const auto& r: res) {
      const Real dmtau_atau_atau = delta_mtau_2loop_atau_atau(r.pars);
      const Real dmtau_atau_at   = delta_mtau_2loop_atau_at(r.pars);
      const Real dmtau_atau_ab   = delta_mtau_2loop_atau_ab(r.pars);

      BOOST_CHECK_CLOSE_FRACTION(dmtau_atau_atau, r.dmtau_atau_atau, 1.0e-10l);
      BOOST_CHECK_CLOSE_FRACTION(dmtau_atau_at  , r.dmtau_atau_at  , 1.0e-10l);
      BOOST_CHECK_CLOSE_FRACTION(dmtau_atau_ab  , r.dmtau_atau_ab  , 1.0e-10l);
   }
}
