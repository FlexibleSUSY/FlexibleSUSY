#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_higher_loop_higgs

#include <boost/test/unit_test.hpp>

#include "sm_fourloophiggs.hpp"
#include "sm_threeloophiggs.hpp"

#include <array>
#include <cmath>
#include <numbers>

using namespace flexiblesusy;

namespace
{

constexpr double zeta3 = 1.2020569031595942853997381615114;
constexpr double li4_half = 0.51747906167389938633075816189887; // Li_4(1/2)
constexpr double pi = std::numbers::pi;
constexpr double ln2 = std::numbers::ln2;

double loop_factor(int loops) { return std::pow(16 * pi * pi, -loops); }

// Benchmark point of arXiv:1508.00912, eqs. (5.6)-(5.9), at Q = M_t
struct Point {
   double Q = 173.34;
   double v = 246.647;
   double lambda = 0.12597;
   double yt = 0.93690;
   double g3 = 1.1666;
   [[nodiscard]] double mt() const { return yt * v / std::numbers::sqrt2; }
   [[nodiscard]] double mh() const { return std::sqrt(2 * lambda) * v; }
};

} // anonymous namespace

// arXiv:1407.4336v3, eqs. (3.3), (3.5)-(3.8)
BOOST_AUTO_TEST_CASE(test_3loop_qcd)
{
   const double c0 = -3776. / 9 + 320 * zeta3 + 704 * std::pow(pi, 4) / 135 +
                     256. / 9 * ln2 * ln2 * (pi * pi - ln2 * ln2) -
                     2048. / 3 * li4_half;
   const double c1 = 128 * zeta3 + 2056. / 3;
   const double d0 =
      760. / 3 - 16 * pi * pi / 3 + 576 * zeta3 + 496 * std::pow(pi, 4) / 15 +
      512. / 3 * ln2 * ln2 * (pi * pi - ln2 * ln2) - 4096 * li4_half;
   const double d1 = -344 + 48 * pi * pi + 960 * zeta3;

   for (const double Q : {86.67, 173.34, 346.68}) {
      Point p;
      p.Q = Q;
      const double t = p.mt() * p.mt();
      const double L = std::log(t / (Q * Q));

      const double as_as = std::pow(p.g3, 4) * p.yt * p.yt * t *
                           (c0 + c1 * L + 160 * L * L - 736 * L * L * L) *
                           loop_factor(3);
      const double at_as = p.g3 * p.g3 * std::pow(p.yt, 4) * t *
                           (d0 + d1 * L - 360 * L * L + 240 * L * L * L) *
                           loop_factor(3);

      BOOST_TEST_CONTEXT("Q = " << Q)
      {
         BOOST_CHECK_CLOSE_FRACTION(
            sm_threeloophiggs::delta_mh_3loop_at_as_as_sm(Q, p.mt(), p.yt,
                                                          p.g3),
            as_as, 1e-12);
         BOOST_CHECK_CLOSE_FRACTION(
            sm_threeloophiggs::delta_mh_3loop_at_at_as_sm(Q, p.mt(), p.yt,
                                                          p.g3),
            at_as, 1e-12);
      }
   }
}

// arXiv:1407.4336v3, eqs. (3.4), (3.9), (3.10) (corrected in v3)
BOOST_AUTO_TEST_CASE(test_3loop_non_qcd)
{
   const double e0 = -673 - 17 * pi * pi / 2 - 1962 * zeta3 -
                     88 * std::pow(pi, 4) / 15 -
                     32 * ln2 * ln2 * (pi * pi - ln2 * ln2) + 768 * li4_half;
   const double e1 = -3759. / 2 - 39 * pi * pi - 144 * zeta3;

   for (const double Q : {86.67, 173.34, 346.68}) {
      Point p;
      p.Q = Q;
      const double t = p.mt() * p.mt();
      const double h = p.mh() * p.mh();
      const double L = std::log(t / (Q * Q));
      const double LH = std::log(h / (Q * Q));

      const double expected =
         std::pow(p.yt, 6) * t *
         (e0 + 36 * LH + e1 * L + 756 * LH * L - 513. / 2 * L * L +
          324 * LH * L * L - 333 * L * L * L) *
         loop_factor(3);

      BOOST_TEST_CONTEXT("Q = " << Q)
      {
         BOOST_CHECK_CLOSE_FRACTION(
            sm_threeloophiggs::delta_mh_3loop_at_at_at_sm(Q, p.mt(), p.yt,
                                                          p.mh()),
            expected, 1e-12);
      }
   }
}

// Delta M_h^2 = (d^2/dv^2 - 1/v d/dv) V, eq. (3.1) of arXiv:1407.4336, with
// the 4-loop effective potential V = g3^6 t^2 P(ln t) of eq. (5.3) of
// arXiv:1508.00912, gives g3^6 yt^2 t 2 (2 P + 3 P' + P''). The decimal
// coefficients of eq. (5.3) are given with 7 significant digits, which
// limits the agreement to about 1e-5 after cancellations.
BOOST_AUTO_TEST_CASE(test_4loop_qcd)
{
   const std::array<double, 5> c = {59366.97, -54056.36, 27699.06, -9144, 1380};

   for (const double Q : {86.67, 173.34, 346.68}) {
      Point p;
      p.Q = Q;
      const double t = p.mt() * p.mt();
      const double L = std::log(t / (Q * Q));

      double P = 0, dP = 0, ddP = 0;
      for (int i = 0; i < 5; ++i) {
         P += c[i] * std::pow(L, i);
         if (i >= 1) {
            dP += i * c[i] * std::pow(L, i - 1);
         }
         if (i >= 2) {
            ddP += i * (i - 1) * c[i] * std::pow(L, i - 2);
         }
      }

      const double expected = std::pow(p.g3, 6) * p.yt * p.yt * t * 2 *
                              (2 * P + 3 * dP + ddP) * loop_factor(4);

      BOOST_TEST_CONTEXT("Q = " << Q)
      {
         BOOST_CHECK_CLOSE_FRACTION(
            sm_fourloophiggs::delta_mh_4loop_at_as_as_as_sm(Q, p.mt(), p.yt,
                                                            p.g3),
            expected, 1e-4);
      }
   }
}
