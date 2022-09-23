#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_twoloop_mt

#include <boost/test/unit_test.hpp>
#include "benchmark.hpp"
#include "mssm_twoloop_mt.hpp"
#include "wrappers.hpp"
#include <cmath>

using namespace flexiblesusy;
using namespace flexiblesusy::mssm_twoloop_mt;

const double g3 = std::sqrt(4*Pi* 0.1184);
const double mt = 173.34;

BOOST_AUTO_TEST_CASE(test_universal)
{
   std::vector<Parameters> pars = {
      // g3 mt mg mst1 mst2 msusy Xt Q
      { g3, mt, 1000., 1000., 1000., 1000.,                  0., 1000. },
      { g3, mt, 1000., 1000., 1000., 1000.,                  0., 1100. },
      { g3, mt, 1000.,  900., 1100., 1000.,  std::sqrt(6.)*1000, 1000. },
      { g3, mt, 1000.,  900., 1100., 1000.,  std::sqrt(6.)*1000, 1100. },
      { g3, mt, 1000.,  900., 1100., 1000., -std::sqrt(6.)*1000, 1000. },
      { g3, mt, 1000.,  900., 1100., 1000., -std::sqrt(6.)*1000, 1100. }
   };

   std::cout << "# universal SUSY masses" << std::endl;
   for (const auto p: pars) {
      const auto dmt = dMt_over_mt_2loop(p);
      BOOST_CHECK(std::isfinite(dmt));
      std::cout << dmt << std::endl;
   }
}

BOOST_AUTO_TEST_CASE(test_non_universal)
{
   std::vector<Parameters> pars = {
      // g3 mt mg mst1 mst2 msusy Xt Q
      { g3, mt, 1000., 1200., 1300., 1400.,                  0., 1000. },
      { g3, mt, 1000., 1200., 1300., 1400.,                  0., 1100. },
      { g3, mt, 1000., 1200., 1300., 1400.,  std::sqrt(6.)*1300, 1000. },
      { g3, mt, 1000., 1200., 1300., 1400.,  std::sqrt(6.)*1300, 1100. },
      { g3, mt, 1000., 1200., 1300., 1400., -std::sqrt(6.)*1300, 1000. },
      { g3, mt, 1000., 1200., 1300., 1400., -std::sqrt(6.)*1300, 1100. }
   };

   std::cout << "# non-universal SUSY masses" << std::endl;
   for (const auto p: pars) {
      const auto dmt_1L_qcd  = dMt_over_mt_1loop_qcd(p);
      const auto dmt_1L_sqcd = dMt_over_mt_1loop_susy(p);
      const auto dmt_2L_qcd  = dMt_over_mt_2loop_qcd(p);
      const auto dmt_2L_sqcd = dMt_over_mt_2loop_susy(p);

      BOOST_CHECK(std::isfinite(dmt_1L_qcd));
      BOOST_CHECK(std::isfinite(dmt_1L_sqcd));
      BOOST_CHECK(std::isfinite(dmt_2L_qcd));
      BOOST_CHECK(std::isfinite(dmt_2L_sqcd));

      std::cout << dmt_1L_qcd << ' '
                << dmt_1L_sqcd << ' '
                << dmt_2L_qcd << ' '
                << dmt_2L_sqcd << std::endl;
   }
}

struct Results {
   Parameters pars;
   double qcd_1l{};
   double sqcd_1l{};
   double qcd_2l{};
   double sqcd_2l{};
};

Parameters make_point()
{
   Parameters pars;

   const double mt = 173.0;
   const double mst1 = 200.0;
   const double mst2 = 300.0;
   const double SX = 0.5;
   const double xt = SX / (2.0 * mt);

   pars.g3    = 2.0;
   pars.mt    = mt;
   pars.mg    = 400.0;
   pars.mst1  = mst1;
   pars.mst2  = mst2;
   pars.msusy = 500.0;
   pars.xt    = xt;
   pars.Q     = 100.0;

   return pars;
}

BOOST_AUTO_TEST_CASE(test_points)
{
   std::vector<Results> res = {
      { make_point(), 0.057796019624082492786, 0.074892665606551035487, 0.042346025594701857529, 0.089577801683073006380 }
   };

   for (const auto& r: res) {
      const double qcd_1l  = dMt_over_mt_1loop_qcd(r.pars);
      const double sqcd_1l = dMt_over_mt_1loop_susy(r.pars);
      const double qcd_2l  = dMt_over_mt_2loop_qcd(r.pars);
      const double sqcd_2l = dMt_over_mt_2loop_susy(r.pars);

      BOOST_CHECK_CLOSE_FRACTION(qcd_1l , r.qcd_1l , 1.0e-10);
      BOOST_CHECK_CLOSE_FRACTION(sqcd_1l, r.sqcd_1l, 1.0e-10);
      BOOST_CHECK_CLOSE_FRACTION(qcd_2l , r.qcd_2l , 1.0e-10);
      BOOST_CHECK_CLOSE_FRACTION(sqcd_2l, r.sqcd_2l, 1.0e-10);
   }
}

BOOST_AUTO_TEST_CASE(test_bench)
{
   const int N = 100'000;
   const Parameters point = make_point();

   Stopwatch sw;
   sw.start();

   for (int i = 0; i < N; ++i) {
      do_not_optimize(dMt_over_mt_1loop_qcd(point));
      do_not_optimize(dMt_over_mt_1loop_susy(point));
      do_not_optimize(dMt_over_mt_2loop_qcd(point));
      do_not_optimize(dMt_over_mt_2loop_susy(point));
   }

   sw.stop();

   const double t = sw.get_time_in_seconds();

   BOOST_TEST_MESSAGE("Evaluation of loop corrections " << N << " times took " << t << " s");
}
