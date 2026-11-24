// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

/**
 * Physical masses depend only on rephasing-invariant combinations of
 * the complex MSSM parameters.  An R rotation (Higgs superfields with
 * R charge 0, matter superfields with R charge 1) maps
 *
 *    M_i -> M_i exp(-i phi),  T_f -> T_f exp(-i phi),  mu -> mu exp(i phi)
 *
 * and leaves B mu, the VEVs, the Yukawa couplings and the soft scalar
 * masses unchanged.  Physical observables must therefore be the same
 * before and after the rotation.  This tests the lepton EDMs, the
 * anomalous magnetic moment of the muon and all decay widths (which also
 * involve the running of the complex parameters to the decay scales).
 * The point has lepton flavour violation (off-diagonal ml2, me2 and
 * complex off-diagonal A_e), so that l -> l' gamma, l -> 3 l' and
 * mu -> e conversion are non-zero and the complex flavour mixing of the
 * sleptons is tested too.
 *
 * A second test rephases the lepton fields, which must leave the same
 * lepton observables unchanged (the Wilson coefficients of l -> 3 l' and
 * mu -> e conversion then change by a common phase).
 *
 * When a new observable is added to FlexibleSUSY, add it here too.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMCPV_rephasing_observables

#include <boost/test/unit_test.hpp>

#include "MSSMCPV_amm.hpp"
#include "MSSMCPV_edm.hpp"
#include "MSSMCPV_input_parameters.hpp"
#include "MSSMCPV_two_scale_spectrum_generator.hpp"
#include "cxx_qft/MSSMCPV_fields.hpp"
#include "decays/MSSMCPV_decays.hpp"
#include "decays/flexibledecay_settings.hpp"
#include "lowe.h"
#include "observables/MSSMCPV_br_l_to_3l.hpp"
#include "observables/MSSMCPV_br_l_to_l_gamma.hpp"
#include "observables/MSSMCPV_l_to_l_conversion.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace flexiblesusy;

namespace
{

MSSMCPV_input_parameters make_input(double phi)
{
   const std::complex<double> e = std::polar(1., phi);
   const std::complex<double> M1 = 300. * e, M2 = 400. * e, M3 = 2000. * e;
   const std::complex<double> A = 1500. * e;
   const double msf2 = 1500. * 1500.;

   MSSMCPV_input_parameters input;
   input.TanBeta = 10.;
   input.MSUSY = 1500.;
   input.ReM1Input = M1.real();
   input.ImM1Input = M1.imag();
   input.ReM2Input = M2.real();
   input.ImM2Input = M2.imag();
   input.ReM3Input = M3.real();
   input.ImM3Input = M3.imag();
   input.ReMuInput = 450.;
   input.ImMuInput = 0.;
   input.mA2Input = 1000. * 1000.;
   input.etaInput = 0.;
   input.mq2Input.setIdentity();
   input.mq2Input *= msf2;
   input.mu2Input = input.mq2Input;
   input.md2Input = input.mq2Input;
   input.ml2Input = input.mq2Input;
   input.me2Input = input.mq2Input;
   input.ReAuInput.setIdentity();
   input.ReAuInput *= A.real();
   input.ReAdInput = input.ReAuInput;
   input.ReAeInput = input.ReAuInput;
   input.ImAuInput.setIdentity();
   input.ImAuInput *= A.imag();
   input.ImAdInput = input.ImAuInput;
   input.ImAeInput = input.ImAuInput;

   // lepton flavour violation, with non-degenerate sleptons (near-degenerate
   // ones would make the mixing, and hence the widths, sensitive to the
   // numerical precision of the RG running to the decay scales)
   input.ml2Input.diagonal() << 1300. * 1300., 1500. * 1500., 1700. * 1700.;
   input.me2Input.diagonal() << 1200. * 1200., 1400. * 1400., 1600. * 1600.;
   input.ml2Input(0, 1) = input.ml2Input(1, 0) = 0.1 * msf2;
   input.me2Input(1, 2) = input.me2Input(2, 1) = 0.1 * msf2;
   const std::complex<double> Ae01 = std::polar(300., 0.7);
   const std::complex<double> Ae12 = std::polar(200., -1.1);
   input.ReAeInput(0, 1) = Ae01.real();
   input.ImAeInput(0, 1) = Ae01.imag();
   input.ReAeInput(1, 2) = Ae12.real();
   input.ImAeInput(1, 2) = Ae12.imag();

   return input;
}

/// applies the R rotation to the parameters of the model
template <class Model>
void rotate(Model& m, double phi)
{
   const std::complex<double> e = std::polar(1., phi);
   m.set_MassB(m.get_MassB() / e);
   m.set_MassWB(m.get_MassWB() / e);
   m.set_MassG(m.get_MassG() / e);
   m.set_TYu(m.get_TYu() / e);
   m.set_TYd(m.get_TYd() / e);
   m.set_TYe(m.get_TYe() / e);
   m.set_Mu(m.get_Mu() * e);
}

/// rephases the lepton (super)fields, l_i -> exp(i alpha_i) l_i, e_i ->
/// exp(i alpha_i) e_i, by redefining the parameters
template <class Model>
void rephase_leptons(Model& m, const std::array<double, 3>& alpha)
{
   Eigen::Matrix<std::complex<double>, 3, 3> D;
   D.setZero();
   for (int i = 0; i < 3; i++) {
      D(i, i) = std::polar(1., alpha[i]);
   }
   m.set_Ye((D * m.get_Ye() * D).eval());
   m.set_TYe((D * m.get_TYe() * D).eval());
   m.set_ml2((D.adjoint() * m.get_ml2() * D).eval());
   // me2 is the soft mass of e_R = e^*, so it transforms the other way
   m.set_me2((D * m.get_me2() * D.adjoint()).eval());
}

/// spectrum at the test point
MSSMCPV_mass_eigenstates run_point(const softsusy::QedQcd& qedqcd,
                                   const Spectrum_generator_settings& settings)
{
   MSSMCPV_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, make_input(std::numbers::pi / 3.));

   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   return spectrum_generator.get_model();
}

Spectrum_generator_settings make_settings()
{
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1e-5);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 1);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, 1);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 1);
   return settings;
}

/// What must agree between two points
enum class Compare : bool {
   all,       ///< everything (R rotation)
   abs_values ///< Wilson coefficients only up to a common phase (rephasing)
};

/// compares the result of an l -> 3 l' or conversion calculation: the
/// branching ratio or conversion rate (element 0) and the Wilson
/// coefficients (elements 1-12), relative to the largest coefficient
void check_lfv_result(const Eigen::Array<std::complex<double>, 13, 1>& a,
                      const Eigen::Array<std::complex<double>, 13, 1>& b,
                      Compare compare, const std::string& name)
{
   BOOST_TEST_MESSAGE(name << ": " << a(0).real() << " vs " << b(0).real());
   BOOST_CHECK(a(0).real() > 0.);
   BOOST_CHECK_CLOSE_FRACTION(a(0).real(), b(0).real(), 1e-8);

   const auto ca = a.tail<12>();
   const auto cb = b.tail<12>();
   const double scale = ca.abs().maxCoeff();
   BOOST_CHECK(scale > 0.);
   const double diff = compare == Compare::all
                          ? (ca - cb).abs().maxCoeff()
                          : (ca.abs() - cb.abs()).abs().maxCoeff();
   BOOST_TEST_MESSAGE(name << ": max. difference of the coefficients "
                           << diff / scale << " (relative to the largest)");
   BOOST_CHECK_SMALL(diff / scale, 1e-8);
}

/// compares the lepton EDMs, a_mu, BR(l -> l' gamma), BR(l -> 3 l') and the
/// mu -> e conversion rates of two models
void check_lepton_observables(const MSSMCPV_mass_eigenstates& mA,
                              const MSSMCPV_mass_eigenstates& mB,
                              const softsusy::QedQcd& qedqcd,
                              const Spectrum_generator_settings& settings,
                              Compare compare)
{
   using MSSMCPV_cxx_diagrams::fields::Fe;

   for (int i = 0; i < 3; i++) {
      const double edmA = MSSMCPV_edm::calculate_edm<Fe>(mA, qedqcd, i);
      const double edmB = MSSMCPV_edm::calculate_edm<Fe>(mB, qedqcd, i);
      BOOST_TEST_MESSAGE("EDM(" << i << "): " << edmA << " vs " << edmB);
      BOOST_CHECK(std::abs(edmA) > 0.);
      BOOST_CHECK_CLOSE_FRACTION(edmA, edmB, 1e-8);
   }

   // the model instantiates the AMM only for the muon
   const double ammA = MSSMCPV_amm::calculate_amm<Fe>(mA, qedqcd, settings, 1);
   const double ammB = MSSMCPV_amm::calculate_amm<Fe>(mB, qedqcd, settings, 1);
   BOOST_TEST_MESSAGE("a_mu: " << ammA << " vs " << ammB);
   BOOST_CHECK_CLOSE_FRACTION(ammA, ammB, 1e-8);

   const Physical_input physical_input;
   const std::array<std::array<int, 2>, 3> transitions{
      {{1, 0}, {2, 0}, {2, 1}}};
   for (const auto& [i, j] : transitions) {
      const double brA = MSSMCPV_br_l_to_l_gamma::calculate_Fe_to_Fe_VP(
         i, j, mA, qedqcd, physical_input);
      const double brB = MSSMCPV_br_l_to_l_gamma::calculate_Fe_to_Fe_VP(
         i, j, mB, qedqcd, physical_input);
      BOOST_TEST_MESSAGE("BR(Fe(" << i << ") -> Fe(" << j << ") gamma): " << brA
                                  << " vs " << brB);
      BOOST_CHECK(brA > 0.);
      BOOST_CHECK_CLOSE_FRACTION(brA, brB, 1e-8);
   }

   // l_I -> l_O l_A l_A-bar
   const std::array<std::array<int, 3>, 5> decays_3l{
      {{1, 0, 0}, {2, 0, 0}, {2, 1, 1}, {2, 1, 0}, {2, 0, 1}}};
   for (const auto& [nI, nO, nA] : decays_3l) {
      std::ostringstream name;
      name << "BR(Fe(" << nI << ") -> Fe(" << nO << ") Fe(" << nA << ") Fe("
           << nA << ")bar)";
      check_lfv_result(
         MSSMCPV_br_l_to_3l::calculate_Fe_to_FeFebarFe_for_All_1loop(
            nI, nO, nA, mA, qedqcd),
         MSSMCPV_br_l_to_3l::calculate_Fe_to_FeFebarFe_for_All_1loop(
            nI, nO, nA, mB, qedqcd),
         compare, name.str());
   }

   using MSSMCPV_l_to_l_conversion::Nucleus;
   const LToLConversion_settings conversion_settings;
   // Ti is listed in Nucleus, but has no data (throws)
   const std::array<std::pair<Nucleus, const char*>, 2> nuclei{
      {{Nucleus::Al, "Al"}, {Nucleus::Au, "Au"}}};
   for (const auto& [nucleus, nucleus_name] : nuclei) {
      check_lfv_result(MSSMCPV_l_to_l_conversion::calculate_FeFe_forAll_1loop(
                          1, 0, nucleus, mA, conversion_settings, qedqcd),
                       MSSMCPV_l_to_l_conversion::calculate_FeFe_forAll_1loop(
                          1, 0, nucleus, mB, conversion_settings, qedqcd),
                       compare,
                       std::string("CR(mu -> e, ") + nucleus_name + ")");
   }
}

/// all partial widths, keyed by (initial PDG, sorted final PDGs)
template <class Table>
std::map<std::vector<int>, double> collect_widths(const Table& table)
{
   std::map<std::vector<int>, double> widths;
   for (const auto& list : table) {
      const std::vector<int> key{list.get_particle_id(), 0};
      widths[key] = list.get_total_width();
      for (const auto& entry : list) {
         std::vector<int> k{list.get_particle_id()};
         auto fs = entry.second.get_final_state_particle_ids();
         std::sort(fs.begin(), fs.end());
         k.insert(k.end(), fs.begin(), fs.end());
         widths[k] = entry.second.get_width();
      }
   }
   return widths;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_MSSMCPV_observables_rephasing)
{
   const double phi = std::numbers::pi / 3.;
   const softsusy::QedQcd qedqcd;
   const auto settings = make_settings();

   auto mA = run_point(qedqcd, settings);
   auto mB = mA;
   rotate(mB, phi);

   for (auto* m : {&mA, &mB}) {
      m->calculate_DRbar_masses();
      m->calculate_pole_masses();
   }

   check_lepton_observables(mA, mB, qedqcd, settings, Compare::all);

   const Physical_input physical_input;
   const FlexibleDecay_settings fd_settings;

   MSSMCPV_decays decA(mA, qedqcd, physical_input, fd_settings);
   MSSMCPV_decays decB(mB, qedqcd, physical_input, fd_settings);
   decA.calculate_decays();
   decB.calculate_decays();

   const auto wA = collect_widths(decA.get_decay_table());
   const auto wB = collect_widths(decB.get_decay_table());

   BOOST_CHECK_EQUAL(wA.size(), wB.size());

   // the widths are calculated after the running to the decay scales, whose
   // numerical error differs between the two points; suppressed channels
   // (e.g. flavour-changing slepton -> slepton Z) are sensitive to it
   const double tol = 1e-7;

   int n_checked = 0, n_failed = 0;
   for (const auto& [key, widthA] : wA) {
      const auto it = wB.find(key);
      BOOST_REQUIRE(it != wB.end());
      const double widthB = it->second;
      if (std::abs(widthA) < 1e-15 && std::abs(widthB) < 1e-15) {
         continue;
      }
      n_checked++;
      const double rel = std::abs(widthA - widthB) /
                         std::max(std::abs(widthA), std::abs(widthB));
      if (rel > tol) {
         n_failed++;
         std::ostringstream ks;
         for (const auto k : key)
            ks << k << ' ';
         BOOST_TEST_MESSAGE("width {" << ks.str() << "}: " << widthA << " vs "
                                      << widthB << " rel. diff. " << rel);
      }
      BOOST_CHECK_CLOSE_FRACTION(widthA, widthB, tol);
   }
   BOOST_TEST_MESSAGE("widths compared: " << n_checked
                                          << ", differing: " << n_failed);
}

/// The decays aren't checked here: the decoupling scheme sets Ye to the
/// real diagonal matrix of the lepton masses (the basis of the low-scale
/// constraint), which undoes the rephasing of Ye but not of TYe, ml2 and
/// me2.
BOOST_AUTO_TEST_CASE(test_MSSMCPV_lepton_rephasing)
{
   const softsusy::QedQcd qedqcd;
   const auto settings = make_settings();

   auto mA = run_point(qedqcd, settings);
   auto mB = mA;
   rephase_leptons(mB, {0.4, -1.3, 2.1});

   for (auto* m : {&mA, &mB}) {
      m->calculate_DRbar_masses();
      m->calculate_pole_masses();
   }

   check_lepton_observables(mA, mB, qedqcd, settings, Compare::abs_values);
}
