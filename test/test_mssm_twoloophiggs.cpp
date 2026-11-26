#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_twoloophiggs

#include <boost/test/unit_test.hpp>
#include "mssm_twoloophiggs.hpp"

#include <Eigen/Core>
#include <array>
#include <cmath>

using namespace flexiblesusy::mssm_twoloophiggs;

namespace {

/// input of the built-in functions, masses squared in GeV^2, mu in the DSZ sign
struct Input {
   double scale2, tanb, vev2, mt2, mb2, mtau2, mA2, mu, mg, mst12, mst22, sxt,
      msb12, msb22, sxb, mstau12, mstau22, sintau, msv2, gs;
};

/// tadpoles over the VEVs, self-energy with tadpoles (11, 12, 22), dMA
using Correction = std::array<double, 6>;

enum Order : int { at_as, ab_as, at_at, atau_atau };

struct Reference_point {
   Input input;
   std::array<Correction, 4> correction; ///< indexed by Order
};

/**
 * Two-loop tadpoles and self-energies of the generated (SARAH) Landau-gauge
 * MSSM in the exact gaugeless limit, per built-in order, from
 * 2L_HiggsMass/theory_calc/MSSM2L_test_reference.wls.  Points: generic;
 * s_t = 0; s_t = 0 and mst1 = mst2; tanb = 1860/61 with mu > 0 (DSZ
 * sign); generic at a low scale.
 */
const std::array<Reference_point, 5> reference_points{{
   {{2.2500000000000000e6, 10.476190476190476, 60516.000000000000, 25600.000000000000, 9.0000000000000000, 3.2400000000000000, 640000.00000000000, -900.00000000000000, 2000.0000000000000, 2.8900000000000000e6, 5.2900000000000000e6, 0.60000000000000000, 3.6100000000000000e6, 4.4104448400000000e6, 0.38461538461538462, 202500.00000000000, 360000.00000000000, 0.28000000000000000, 214844.76000000000, 1.0500000000000000},
    {{{-13149.582211578270, -3294.0460336232851, 5.6747583922929109, 11.766010539420002, -2727.0970155501333, 13100.180893105617},
      {-62146.599562911706, -63.299362785693867, 139.22917727382962, 1.4795009721196040, 0.015671869422135787, 6308.6257606658403},
      {-2.9420560864103663e6, -197445.80521050473, 586.87479966865700, 954.70694094999501, -261.26412333299123, 921944.88911545146},
      {2387.4312611011902, 14.711380970129796, 31.816439265255273, 2.0405414564548305, 0.13086989488413719, -973.54036413335013}}}},
   {{2.2500000000000000e6, 10.476190476190476, 60516.000000000000, 25600.000000000000, 9.0000000000000000, 3.2400000000000000, 640000.00000000000, -900.00000000000000, 2000.0000000000000, 2.8900000000000000e6, 5.2900000000000000e6, 0., 2.5600000000000000e6, 4.6178048400000000e6, 0.38461538461538462, 202500.00000000000, 360000.00000000000, 0.28000000000000000, 214844.76000000000, 1.0500000000000000},
    {{{904.16619273178563, 6083.4927294158110, 0., -34.028520842614628, -2123.6349376023352, -1077.2671323318615},
      {-273765.55416645851, -100.69546800302044, 12919.307896911734, 49.694477601543311, 0.19113836474237096, 10719.790434439696},
      {-1.0479552564842855e8, 54519.772855652897, 1.1873876886177301e6, -31741.542560413654, -1970.0940865949872, 3.9936891548645536e6},
      {2387.4312611011902, 14.711380970129796, 31.816439265255273, 2.0405414564548305, 0.13086989488413719, -973.54036413335013}}}},
   {{2.2500000000000000e6, 10.476190476190476, 60516.000000000000, 25600.000000000000, 9.0000000000000000, 3.2400000000000000, 640000.00000000000, -900.00000000000000, 2000.0000000000000, 2.8900000000000000e6, 2.8900000000000000e6, 0., 2.5600000000000000e6, 4.6178048400000000e6, 0.38461538461538462, 202500.00000000000, 360000.00000000000, 0.28000000000000000, 214844.76000000000, 1.0500000000000000},
    {{{-1122.9656386348767, 5158.3055549926977, 0., -42.209324519698560, -2062.1402300796290, 1159.5929058533363},
      {-273765.55416645851, -100.69546800302044, 12919.307896911734, 49.694477601543311, 0.19113836474237096, 10719.790434439696},
      {-1.0474750327846115e8, 54941.531222248441, 1.1887825049662552e6, -31631.412813974666, -1951.5115138787760, 3.9781605112373805e6},
      {2387.4312611011902, 14.711380970129796, 31.816439265255273, 2.0405414564548305, 0.13086989488413719, -973.54036413335013}}}},
   {{2.2500000000000000e6, 30.491803278688525, 60516.000000000000, 25600.000000000000, 9.0000000000000000, 3.2400000000000000, 640000.00000000000, 700.00000000000000, 2000.0000000000000, 2.8900000000000000e6, 5.2900000000000000e6, 0.60000000000000000, 3.6100000000000000e6, 4.4104448400000000e6, 0.38461538461538462, 202500.00000000000, 360000.00000000000, 0.28000000000000000, 214844.76000000000, 1.0500000000000000},
    {{{29530.777433075150, -3418.4353276172948, 3.4055409673694681, -9.6081051054671307, -2702.4737161784221, -29663.280060433722},
      {-712833.89065872929, 142.15504066995345, 2174.3335869961844, -13.173282827845884, 0.079674433446685436, -162296.06314552081},
      {-3.5484092543243351e8, -1.7824260817397643e6, 270937.11343632821, 97966.570861712938, -11938.049605810188, -1.0507741791984730e8},
      {5.4415518316632803e6, -2789.6801623495510, 228386.49314301546, -3566.7880994138819, 55.703717836971368, 4.9557530483693141e6}}}},
   {{25600.000000000000, 10.476190476190476, 60516.000000000000, 25600.000000000000, 9.0000000000000000, 3.2400000000000000, 640000.00000000000, -900.00000000000000, 2000.0000000000000, 2.8900000000000000e6, 5.2900000000000000e6, 0.60000000000000000, 3.6100000000000000e6, 4.4104448400000000e6, 0.38461538461538462, 202500.00000000000, 360000.00000000000, 0.28000000000000000, 214844.76000000000, 1.0500000000000000},
    {{{595959.21565988477, 588037.39107871218, -28.371582925893745, -104.62242240451363, 1351.0159317563325, -595342.22311789763},
      {2.3799951784129958e6, 2357.9952156269885, 38.148588915822849, 0.68544514445667192, 0.010199136610213090, -235518.55282400398},
      {-1.9398421644397755e7, -2.4022555841286992e6, -40996.700888004167, -4398.9621734699069, -31077.526766755676, 4.6384460385789362e6},
      {-6423.7106095323649, -39.211054021982550, -51.320519409909216, -3.2901120073052260, -0.21092708197507879, 2600.9212581769191}}}}
}};

/// results of the built-in functions
struct Builtin {
   Eigen::Matrix<double, 2, 1> tadpole;
   Eigen::Matrix<double, 2, 2> self_energy;
   Eigen::Matrix<double, 2, 2> self_energy_with_tadpoles;
   Eigen::Matrix<double, 2, 2> pseudoscalar;
   double dma;
};

Builtin builtin(const Input& p, Order order)
{
   const double cxt = std::sqrt(1 - p.sxt*p.sxt);
   const double cxb = std::sqrt(1 - p.sxb*p.sxb);
   const double costau = std::sqrt(1 - p.sintau*p.sintau);
   const double cotb = 1/p.tanb;

   switch (order) {
   case at_as:
      return {
         tadpole_higgs_2loop_at_as_mssm(p.mt2, p.mg, p.mst12, p.mst22, p.sxt, cxt,
                                        p.scale2, p.mu, p.tanb, p.vev2, p.gs),
         self_energy_higgs_2loop_at_as_mssm(p.mt2, p.mg, p.mst12, p.mst22, p.sxt,
                                            cxt, p.scale2, p.mu, p.tanb, p.vev2,
                                            p.gs),
         self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
            p.mt2, p.mg, p.mst12, p.mst22, p.sxt, cxt, p.scale2, p.mu, p.tanb,
            p.vev2, p.gs),
         self_energy_pseudoscalar_2loop_at_as_mssm(
            p.mt2, p.mg, p.mst12, p.mst22, p.sxt, cxt, p.scale2, p.mu, p.tanb,
            p.vev2, p.gs),
         self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
            p.mt2, p.mg, p.mst12, p.mst22, p.sxt, cxt, p.scale2, p.mu, p.tanb,
            p.vev2, p.gs)};
   case ab_as:
      return {
         tadpole_higgs_2loop_ab_as_mssm(p.mb2, p.mg, p.msb12, p.msb22, p.sxb, cxb,
                                        p.scale2, p.mu, cotb, p.vev2, p.gs),
         self_energy_higgs_2loop_ab_as_mssm(p.mb2, p.mg, p.msb12, p.msb22, p.sxb,
                                            cxb, p.scale2, p.mu, cotb, p.vev2,
                                            p.gs),
         self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
            p.mb2, p.mg, p.msb12, p.msb22, p.sxb, cxb, p.scale2, p.mu, cotb,
            p.vev2, p.gs),
         self_energy_pseudoscalar_2loop_ab_as_mssm(
            p.mb2, p.mg, p.msb12, p.msb22, p.sxb, cxb, p.scale2, p.mu, cotb,
            p.vev2, p.gs),
         self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
            p.mb2, p.mg, p.msb12, p.msb22, p.sxb, cxb, p.scale2, p.mu, cotb,
            p.vev2, p.gs)};
   case at_at:
      return {
         tadpole_higgs_2loop_at_at_mssm(p.mt2, p.mb2, p.mA2, p.mst12, p.mst22,
                                        p.msb12, p.msb22, p.sxt, cxt, p.sxb, cxb,
                                        p.scale2, p.mu, p.tanb, p.vev2),
         self_energy_higgs_2loop_at_at_mssm(p.mt2, p.mb2, p.mA2, p.mst12, p.mst22,
                                            p.msb12, p.msb22, p.sxt, cxt, p.sxb,
                                            cxb, p.scale2, p.mu, p.tanb, p.vev2),
         self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
            p.mt2, p.mb2, p.mA2, p.mst12, p.mst22, p.msb12, p.msb22, p.sxt, cxt,
            p.sxb, cxb, p.scale2, p.mu, p.tanb, p.vev2),
         self_energy_pseudoscalar_2loop_at_at_mssm(
            p.mt2, p.mb2, p.mA2, p.mst12, p.mst22, p.msb12, p.msb22, p.sxt, cxt,
            p.sxb, cxb, p.scale2, p.mu, p.tanb, p.vev2),
         self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
            p.mt2, p.mb2, p.mA2, p.mst12, p.mst22, p.msb12, p.msb22, p.sxt, cxt,
            p.sxb, cxb, p.scale2, p.mu, p.tanb, p.vev2)};
   case atau_atau:
      break;
   }

   return {
      tadpole_higgs_2loop_atau_atau_mssm(p.mtau2, p.mA2, p.msv2, p.mstau12,
                                         p.mstau22, p.sintau, costau, p.scale2,
                                         p.mu, p.tanb, p.vev2),
      self_energy_higgs_2loop_atau_atau_mssm(p.mtau2, p.mA2, p.msv2, p.mstau12,
                                             p.mstau22, p.sintau, costau,
                                             p.scale2, p.mu, p.tanb, p.vev2),
      self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
         p.mtau2, p.mA2, p.msv2, p.mstau12, p.mstau22, p.sintau, costau, p.scale2,
         p.mu, p.tanb, p.vev2),
      self_energy_pseudoscalar_2loop_atau_atau_mssm(
         p.mtau2, p.mA2, p.msv2, p.mstau12, p.mstau22, p.sintau, costau, p.scale2,
         p.mu, p.tanb, p.vev2),
      self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
         p.mtau2, p.mA2, p.msv2, p.mstau12, p.mstau22, p.sintau, costau, p.scale2,
         p.mu, p.tanb, p.vev2)};
}

/// reference self-energy with tadpoles as a matrix
Eigen::Matrix<double, 2, 2> self_energy_with_tadpoles(const Correction& c)
{
   Eigen::Matrix<double, 2, 2> result;
   result << c[2], c[3], c[3], c[4];
   return result;
}

/// dMA rotated into the CP-even (sign = -1) or CP-odd (sign = +1) basis
Eigen::Matrix<double, 2, 2> rotate(double dma, double tanb, double sign)
{
   Eigen::Matrix<double, 2, 2> result;
   result << tanb*tanb, sign*tanb, sign*tanb, 1;
   return dma/(1 + tanb*tanb)*result;
}

/**
 * Relative tolerance per order and point index.  The O((at+ab)^2) routines
 * (ddstad, ddshiggs, ddsodd) replace s_t -> 1e-5 (eps_st) and T2 -> T1 (1 + 1e-5)
 * (eps_t1) internally, which limits their accuracy at the points 2 and 3 with
 * s_t = 0 and mst1 = mst2 to O(1e-4).
 */
double tolerance(Order order, std::size_t point)
{
   if (order == at_at && (point == 1 || point == 2)) {
      return 2e-4;
   }
   return 1e-9;
}

/// |a - b| <= eps*scale element-wise
template <typename T>
void check_close(const T& a, const T& b, double scale, double eps)
{
   BOOST_CHECK_SMALL(((a - b)/scale).cwiseAbs().maxCoeff(), eps);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_with_tadpoles_vs_generated )
{
   for (std::size_t n = 0; n < reference_points.size(); ++n) {
      for (int o = at_as; o <= atau_atau; ++o) {
         // Arrange
         const auto order = static_cast<Order>(o);
         const Reference_point& ref = reference_points[n];
         const Correction& c = ref.correction[order];
         const Eigen::Matrix<double, 2, 1> tadpole(c[0], c[1]);
         const Eigen::Matrix<double, 2, 2> self_energy = self_energy_with_tadpoles(c);
         const double eps = tolerance(order, n);
         BOOST_TEST_MESSAGE("point " << n + 1 << ", order " << o);

         // Act
         const Builtin b = builtin(ref.input, order);

         // Assert
         check_close(b.tadpole, tadpole, tadpole.cwiseAbs().maxCoeff(), eps);
         check_close(b.self_energy_with_tadpoles, self_energy,
                     self_energy.cwiseAbs().maxCoeff(), eps);
         BOOST_CHECK_SMALL(std::abs((b.dma - c[5])/c[5]), eps);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_without_tadpoles_vs_generated )
{
   for (std::size_t n = 0; n < reference_points.size(); ++n) {
      for (int o = at_as; o <= atau_atau; ++o) {
         // Arrange
         const auto order = static_cast<Order>(o);
         const Reference_point& ref = reference_points[n];
         const Correction& c = ref.correction[order];
         const double tanb = ref.input.tanb;
         const Eigen::Matrix<double, 2, 2> tadpoles =
            Eigen::Matrix<double, 2, 1>(c[0], c[1]).asDiagonal();
         const Eigen::Matrix<double, 2, 2> self_energy =
            self_energy_with_tadpoles(c) + rotate(c[5], tanb, -1) + tadpoles;
         const Eigen::Matrix<double, 2, 2> pseudoscalar =
            rotate(c[5], tanb, 1) + tadpoles;
         const double eps = tolerance(order, n);
         BOOST_TEST_MESSAGE("point " << n + 1 << ", order " << o);

         // Act
         const Builtin b = builtin(ref.input, order);

         // Assert
         check_close(b.self_energy, self_energy, self_energy.cwiseAbs().maxCoeff(),
                     eps);
         check_close(b.pseudoscalar, pseudoscalar,
                     pseudoscalar.cwiseAbs().maxCoeff(), eps);
      }
   }
}