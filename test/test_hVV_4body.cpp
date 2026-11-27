#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_hVV_4body

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"
#include "decays/decay_functions.hpp"
#include "wrappers.hpp"

#include <cmath>

using namespace flexiblesusy;

namespace {

constexpr double mW = 80.385;
constexpr double GammaW = 2.085;
constexpr double mZ = 91.1876;
constexpr double GammaZ = 2.4952;

void check(double mH, double mV, double GammaV, double expected)
{
   const auto result = hVV_4body(mH, mV, GammaV);
   BOOST_CHECK_CLOSE_FRACTION(result.value, expected, 1e-6);
   // error estimate is below the threshold for the warning in the decay code
   BOOST_CHECK_LE(result.error, 1e-3*result.value);
}

} // anonymous namespace

// The reference values were obtained with a tolerance of 1e-9.  They agree
// with 1000000 point MISER integrations within the statistical errors.

BOOST_AUTO_TEST_CASE(test_below_single_vector_threshold)
{
   check(70., mW, GammaW, 3.330360444e-05);
}

BOOST_AUTO_TEST_CASE(test_between_thresholds)
{
   check(125.09, mW, GammaW, 0.01304931174);
   check(125.09, mZ, GammaZ, 0.003260272064);
   check(155., mW, GammaW, 0.1830731849);
}

BOOST_AUTO_TEST_CASE(test_close_to_pair_threshold)
{
   check(160.7, mW, GammaW, 0.6289885477);
   check(182., mZ, GammaZ, 0.5955580784);
}

BOOST_AUTO_TEST_CASE(test_deterministic)
{
   const auto first = hVV_4body(125.09, mW, GammaW);
   const auto second = hVV_4body(125.09, mW, GammaW);
   BOOST_CHECK_EQUAL(first.value, second.value);
}

namespace
{

/// SM widths of H -> VV, with the normalizations of
/// src/decays/specializations/H/decay_H_to_{WpWm,ZZ}.inc
struct SM_HVV {
   double mV{};
   double GammaV{};
   bool is_Z{};
   int offshell_VV_decays{};
   int four_body_calls{};

   static constexpr double v = 246.22;
   double g2() const { return 2. * mW / v; }
   double cw2() const { return Sqr(mW / mZ); }
   double ghVV() const
   {
      return is_Z ? g2() * mZ / std::sqrt(cw2()) : g2() * mW;
   }
   double symmetry() const { return is_Z ? 0.5 : 1.; }

   double width(double mH)
   {
      const double x = Sqr(mV / mH);

      const auto four_body = [&] {
         four_body_calls++;
         const auto integral = hVV_4body(mH, mV, GammaV);
         const double norm = 1. / (64. * Cube(Pi)) * Sqr(ghVV()) * Cube(mH) /
                             Power4(mV) * symmetry();
         return HVV_4body_result{integral.value * norm, integral.error * norm};
      };
      const auto three_body = [&] {
         if (is_Z) {
            const double sw2 = 1. - cw2();
            const double deltaV =
               7. / 12. - 10. / 9. * sw2 + 40. / 27. * Sqr(sw2);
            return 3. / (512. * Power3(Pi)) / mH * deltaV * RT(x) / x *
                   Sqr(ghVV() * g2()) / (1. - sw2);
         }
         return 1. / (768. * Power3(Pi) * mH) * RT(x) / x * Sqr(ghVV() * g2()) /
                2. * 9.;
      };
      const auto two_body = [&] {
         return Sqr(ghVV()) * Cube(mH) / (64. * Pi * Power4(mV)) *
                std::sqrt(KallenLambda(1., x, x)) *
                (1. - 4. * x + 12. * x * x) * symmetry();
      };

      return hVV_width(x, offshell_VV_decays, four_body, three_body, two_body)
         .value;
   }
};

/// checks that the width changes by less than 1% between neighbouring
/// points in [mH_min, mH_max]
void check_continuous(SM_HVV& hvv, double mH_min, double mH_max, double step)
{
   const int n = static_cast<int>((mH_max - mH_min) / step);
   double previous = hvv.width(mH_min);
   for (int i = 1; i <= n; ++i) {
      const double mH = mH_min + i * step;
      const double current = hvv.width(mH);
      BOOST_TEST_CONTEXT("mH = " << mH << ", offshell_VV_decays = "
                                 << hvv.offshell_VV_decays
                                 << ", mV = " << hvv.mV)
      {
         BOOST_TEST(std::isfinite(current));
         BOOST_TEST(current > 0.);
         BOOST_TEST(std::abs(current - previous) < 1e-2 * previous);
      }
      previous = current;
   }
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_width_continuous_at_VV_threshold)
{
   for (const int flag : {1, 2}) {
      SM_HVV W{
         .mV = mW, .GammaV = GammaW, .is_Z = false, .offshell_VV_decays = flag};
      check_continuous(W, 2. * mW - 2., 2. * mW + 4., 0.005);
      SM_HVV Z{
         .mV = mZ, .GammaV = GammaZ, .is_Z = true, .offshell_VV_decays = flag};
      check_continuous(Z, 2. * mZ - 2., 2. * mZ + 4., 0.005);
   }
}

BOOST_AUTO_TEST_CASE(test_width_continuous_below_VV_threshold)
{
   // single and double off-shell widths cross at about 0.96*2*mV
   SM_HVV W{.mV = mW, .GammaV = GammaW, .is_Z = false, .offshell_VV_decays = 1};
   check_continuous(W, 146., 2. * mW - 1., 0.01);
   SM_HVV Z{.mV = mZ, .GammaV = GammaZ, .is_Z = true, .offshell_VV_decays = 1};
   check_continuous(Z, 166., 2. * mZ - 1., 0.01);
}

BOOST_AUTO_TEST_CASE(test_width_exactly_at_VV_threshold)
{
   for (const int flag : {1, 2}) {
      SM_HVV W{
         .mV = mW, .GammaV = GammaW, .is_Z = false, .offshell_VV_decays = flag};
      const double width = W.width(2. * mW);
      BOOST_TEST(width == hVV_4body(2. * mW, mW, GammaW).value /
                             (64. * Cube(Pi)) * Sqr(W.ghVV()) * Cube(2. * mW) /
                             Power4(mW),
                 boost::test_tools::tolerance(1e-12));
   }
}

BOOST_AUTO_TEST_CASE(test_width_formulas_away_from_thresholds)
{
   // on-shell width only, without calculating the double off-shell one
   for (const int flag : {0, 1, 2}) {
      SM_HVV W{
         .mV = mW, .GammaV = GammaW, .is_Z = false, .offshell_VV_decays = flag};
      const double mH = 200.;
      const double x = Sqr(mW / mH);
      BOOST_TEST(W.width(mH) == Sqr(W.ghVV()) * Cube(mH) /
                                   (64. * Pi * Power4(mW)) *
                                   std::sqrt(KallenLambda(1., x, x)) *
                                   (1. - 4. * x + 12. * x * x),
                 boost::test_tools::tolerance(1e-14));
      BOOST_TEST(W.four_body_calls == 0);
   }

   // no off-shell decays
   SM_HVV W0{
      .mV = mW, .GammaV = GammaW, .is_Z = false, .offshell_VV_decays = 0};
   BOOST_TEST(W0.width(125.) == 0.);
   BOOST_TEST(W0.width(2. * mW) == 0.);
   BOOST_TEST(W0.width(70.) == 0.);
   BOOST_TEST(W0.four_body_calls == 0);

   // single off-shell width only
   SM_HVV W1{
      .mV = mW, .GammaV = GammaW, .is_Z = false, .offshell_VV_decays = 1};
   const double mH = 125.;
   const double x = Sqr(mW / mH);
   BOOST_TEST(W1.width(mH) == 1. / (768. * Power3(Pi) * mH) * RT(x) / x *
                                 Sqr(W1.ghVV() * W1.g2()) / 2. * 9.,
              boost::test_tools::tolerance(1e-14));
   BOOST_TEST(W1.four_body_calls == 0);

   // double off-shell width below the single vector threshold
   BOOST_TEST(W1.width(70.) > 0.);
   BOOST_TEST(W1.four_body_calls == 1);
}
