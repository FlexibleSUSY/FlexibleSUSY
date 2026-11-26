#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_twoloophiggs

#include <boost/test/unit_test.hpp>
#include "sm_twoloophiggs.hpp"

#include <array>
#include <cmath>

using namespace flexiblesusy::sm_twoloophiggs;

namespace {

constexpr double Pi = 3.141592653589793;
constexpr double twoLoop = 4.010149318236068752e-05; // 1/(4 Pi)^4

double sqr(double x) noexcept { return x*x; }
double pow3(double x) noexcept { return x*x*x; }
double pow4(double x) noexcept { return sqr(sqr(x)); }

/// point of the generated 2-loop diagrams, without the factor 1/(16 Pi^2)^2
struct Generated_point {
   double Q2, mt, mb, mtau, mh, yt, ytau, g3;
   double tadpole, tadpole_no_g3, self_energy, self_energy_no_g3;
};

/**
 * tadpole_hh_2loop()/v and self_energy_hh_2loop() of the generated
 * Landau-gauge SM with lambda -> 0 and m_h set by hand, from
 * 2L_HiggsMass/theory_calc/checks/diagrams_ref.csv.  Points 3 and 4
 * have m_h = 200 GeV and enlarged y_b = 0.1, y_tau = 0.05,
 * respectively.
 */
const std::array<Generated_point, 4> generated_points{{
   {30046.7556, 163.63218384676523, 2.7227637839558421, 1.7497734255695216,
    124.1313667278574, 0.93516212702513601, 0.01, 1.1648463470930048,
    -1883238.2433658631, 971803.69602276909, -1645179.4676236892, 1133209.86625979},
   {250000, 163.63218384676523, 2.7227637839558421, 1.7497734255695216,
    124.1313667278574, 0.93516212702513601, 0.01, 1.1648463470930048,
    -11562642.208930062, 4268509.0947625097, -20695636.442419447, 8086378.2567607174},
   {30046.7556, 163.63218384676523, 2.7227637839558421, 1.7497734255695216,
    200, 0.93516212702513601, 0.01, 1.1648463470930048,
    -2043917.9336217134, 811124.00576691865, -2336913.4474476278, 441475.88643585122},
   {30046.7556, 163.63218384676523, 17.497734255695217, 8.7488671278476087,
    124.1313667278574, 0.93516212702513601, 0.050000000000000003, 1.1648463470930048,
    -1899835.7624272695, 961247.4002247582, -1682451.241175825, 1109743.1146792651},
}};

/// previous built-in O((alpha_t + alpha_b)^2) self-energy (m_h = 0, expanded in r)
double self_energy_at_at_expanded(double scale, double mt, double yt, double mb)
{
   const double Pi2 = sqr(Pi);
   const double yt4 = pow4(yt);
   const double mt2 = sqr(mt);
   const double r = sqr(mb)/mt2;
   const double LogT = std::log(mt2/sqr(scale));
   const double LogT2 = sqr(LogT);
   const double LogTB = 0.5*std::log(r);
   const double LogTB2 = sqr(LogTB);

   return twoLoop*(
      3*mt2*yt4*(19 + Pi2 - 27*LogT + 9*LogT2)
      - 3*mt2*yt4*(5 + 3*Pi2 - 3*LogT + 9*LogT2)*r
      + 1.5*mt2*yt4*(35 + 6*Pi2 + 6*LogT - 18*LogT2 - 24*LogTB*(1 + 3*LogT))*sqr(r)
      - 0.5*mt2*yt4*(-29 + 6*Pi2 - 144*LogTB2 + 162*LogT - 54*LogT2
                     - 24*LogTB*(-8 + 9*LogT))*pow3(r));
}

/// previous built-in O((alpha_t + alpha_b)^2) tadpole (m_h = 0, expanded in r)
double tadpole_at_at_expanded(double scale, double mt, double yt, double mb)
{
   const double Pi2 = sqr(Pi);
   const double yt4 = pow4(yt);
   const double mt2 = sqr(mt);
   const double r = sqr(mb)/mt2;
   const double LogT = std::log(mt2/sqr(scale));
   const double LogT2 = sqr(LogT);
   const double LogTB = 0.5*std::log(r);
   const double LogTB2 = sqr(LogTB);

   return twoLoop*(
      mt2*yt4*(45 + Pi2 - 39*LogT + 9*LogT2)
      - 3*mt2*yt4*(5 + Pi2 - 5*LogT + 3*LogT2)*r
      + 1.5*mt2*yt4*(5 + 2*Pi2 + LogTB*(8 - 24*LogT) + 10*LogT - 6*LogT2)*sqr(r)
      - 1./6.*mt2*yt4*(-185 + 6*Pi2 - 144*LogTB2 + 234*LogT - 54*LogT2
                       - 24*LogTB*(-14 + 9*LogT))*pow3(r));
}

/// previous built-in O(alpha_tau^2) self-energy (m_h = 0)
double self_energy_atau_atau_massless_h(double scale, double mtau, double ytau)
{
   const double L = std::log(sqr(mtau/scale));
   return twoLoop*sqr(mtau)*pow4(ytau)*(19 + sqr(Pi) - 27*L + 9*sqr(L));
}

/// previous built-in O(alpha_tau^2) tadpole (m_h = 0)
double tadpole_atau_atau_massless_h(double scale, double mtau, double ytau)
{
   const double L = std::log(sqr(mtau/scale));
   return twoLoop/3*sqr(mtau)*pow4(ytau)*(45 + sqr(Pi) - 39*L + 9*sqr(L));
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_generated_diagrams )
{
   const double eps = 1e-9;

   for (const auto& p : generated_points) {
      // Arrange
      const double scale = std::sqrt(p.Q2);
      const double yb = p.yt*p.mb/p.mt;

      // Act
      const double tadpole_no_g3 =
         (tadpole_higgs_2loop_at_at_sm(scale, p.mt, p.yt, p.mb, p.mh)
          + tadpole_higgs_2loop_atau_atau_sm(scale, p.mtau, p.ytau, p.mh))/twoLoop;
      const double self_energy_no_g3 =
         (self_energy_higgs_2loop_at_at_sm(0, scale, p.mt, p.yt, p.mb, p.mh)
          + self_energy_higgs_2loop_atau_atau_sm(0, scale, p.mtau, p.ytau, p.mh))
         /twoLoop;
      const double tadpole_g3 =
         (tadpole_higgs_2loop_at_as_sm(scale, p.mt, p.yt, p.g3)
          + tadpole_higgs_2loop_ab_as_sm(scale, p.mb, yb, p.g3))/twoLoop;
      const double self_energy_g3 =
         (self_energy_higgs_2loop_at_as_sm(0, scale, p.mt, p.yt, p.g3)
          + self_energy_higgs_2loop_ab_as_sm(0, scale, p.mb, yb, p.g3))/twoLoop;

      // Assert
      BOOST_CHECK_CLOSE_FRACTION(tadpole_no_g3, p.tadpole_no_g3, eps);
      BOOST_CHECK_CLOSE_FRACTION(self_energy_no_g3, p.self_energy_no_g3, eps);
      BOOST_CHECK_CLOSE_FRACTION(tadpole_g3, p.tadpole - p.tadpole_no_g3, eps);
      BOOST_CHECK_CLOSE_FRACTION(
         self_energy_g3, p.self_energy - p.self_energy_no_g3, eps);
   }
}

BOOST_AUTO_TEST_CASE( test_massless_higgs_limit )
{
   const double eps = 1e-9;
   const double mt = 163.63218384676523, yt = 0.93516212702513601;
   const double mb = 1e-2*mt; // r = mb^2/mt^2 = 1e-4
   const double mtau = 1.7497734255695216, ytau = 0.01;

   // m_h -> 0 numerically and exactly (m_h = 0, used by the model code)
   for (const auto [scale, ratio] : {std::array{mt, 1e-6}, std::array{500., 1e-6},
                                     std::array{mt, 0.}, std::array{500., 0.}}) {
      // Arrange
      const double mh_t = ratio*mt;
      const double mh_tau = ratio*mtau;

      // Act
      const double se =
         self_energy_higgs_2loop_at_at_sm(0, scale, mt, yt, mb, mh_t);
      const double ta = tadpole_higgs_2loop_at_at_sm(scale, mt, yt, mb, mh_t);
      const double se_tau =
         self_energy_higgs_2loop_atau_atau_sm(0, scale, mtau, ytau, mh_tau);
      const double ta_tau =
         tadpole_higgs_2loop_atau_atau_sm(scale, mtau, ytau, mh_tau);

      // Assert
      BOOST_CHECK_CLOSE_FRACTION(
         se, self_energy_at_at_expanded(scale, mt, yt, mb), eps);
      BOOST_CHECK_CLOSE_FRACTION(
         ta, tadpole_at_at_expanded(scale, mt, yt, mb), eps);
      BOOST_CHECK_CLOSE_FRACTION(
         se_tau, self_energy_atau_atau_massless_h(scale, mtau, ytau), eps);
      BOOST_CHECK_CLOSE_FRACTION(
         ta_tau, tadpole_atau_atau_massless_h(scale, mtau, ytau), eps);
   }
}

BOOST_AUTO_TEST_CASE( test_lepton_limit_of_quark_result )
{
   const double eps = 1e-9;
   const double mtau = 1.7497734255695216, ytau = 0.01;
   // suppresses the y_nu terms below eps; phi_xyz overflows for mnu^2/mh^2 < 1e-16
   const double mnu = 1e-5*mtau;
   const double mh = 124.1313667278574;

   for (const double scale : {173.34, 500.}) {
      // Arrange: quark result with N_c = 3 and (t, b) -> (tau, nu)
      const double se_quark =
         self_energy_higgs_2loop_at_at_sm(0, scale, mtau, ytau, mnu, mh);
      const double ta_quark =
         tadpole_higgs_2loop_at_at_sm(scale, mtau, ytau, mnu, mh);

      // Act
      const double se_tau =
         self_energy_higgs_2loop_atau_atau_sm(0, scale, mtau, ytau, mh);
      const double ta_tau = tadpole_higgs_2loop_atau_atau_sm(scale, mtau, ytau, mh);

      // Assert
      BOOST_CHECK_CLOSE_FRACTION(se_tau, se_quark/3, eps);
      BOOST_CHECK_CLOSE_FRACTION(ta_tau, ta_quark/3, eps);
   }
}
