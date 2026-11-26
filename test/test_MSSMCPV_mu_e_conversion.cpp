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
 * Compares the mu -> e conversion rates in aluminium, titanium and gold
 * with the rates from the Wilson coefficients of SARAH 4.15.5/SPheno 4.0.7
 * (FlavorKit), put into eq. (14) of Kitano, Koike, Okada [hep-ph/0203110]
 * with the overlap integrals and capture rates used by FlexibleSUSY
 * (tables 2, 4 for Ti, and 8).  At these
 * points of a low-scale MSSM with mixing of the right-handed smuon and
 * selectron, the photon dipole and monopole (charge radius) contributions
 * are of similar size and add up.  Until 2026 the dipole contribution had
 * the wrong sign relative to the others, which gave a conversion rate
 * smaller by a factor 27 (mu = 500 GeV) or 7 (mu = 800 GeV) in Al, and
 * titanium had no overlap integrals (an exception was thrown).  The
 * remaining differences of a few percent come from the different masses
 * and couplings of the two codes (BR(mu -> e gamma) differs by 10%).
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMCPV_mu_e_conversion

#include <boost/test/unit_test.hpp>

#include "MSSMCPV_input_parameters.hpp"
#include "MSSMCPV_two_scale_spectrum_generator.hpp"
#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "observables/MSSMCPV_l_to_l_conversion.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "spectrum_generator_settings.hpp"

#include <array>
#include <cstddef>
#include <utility>

using namespace flexiblesusy;

namespace
{

/// low-scale MSSM: tan(beta) = 3, M1, M2, M3 = 300, 600, 2000 GeV,
/// mA = 1 TeV, sleptons 800 GeV, squarks 2 TeV, A = 0, (me2)_12 = 0.1 *
/// (800 GeV)^2, all at Q = 1 TeV
MSSMCPV_input_parameters make_input(double mu)
{
   MSSMCPV_input_parameters input;
   input.TanBeta = 3.;
   input.MSUSY = 1000.;
   input.ReM1Input = 300.;
   input.ReM2Input = 600.;
   input.ReM3Input = 2000.;
   input.ReMuInput = mu;
   input.mA2Input = 1000. * 1000.;
   input.etaInput = 0.;
   input.mq2Input.setIdentity();
   input.mq2Input *= 2000. * 2000.;
   input.mu2Input = input.mq2Input;
   input.md2Input = input.mq2Input;
   input.ml2Input.setIdentity();
   input.ml2Input *= 800. * 800.;
   input.me2Input = input.ml2Input;
   input.me2Input(0, 1) = input.me2Input(1, 0) = 0.1 * 800. * 800.;
   input.ReAuInput.setZero();
   input.ReAdInput.setZero();
   input.ReAeInput.setZero();
   input.ImAuInput.setZero();
   input.ImAdInput.setZero();
   input.ImAeInput.setZero();
   return input;
}

using MSSMCPV_l_to_l_conversion::Nucleus;

/// CR(mu -> e) in the nuclei Al, Ti, Au
std::array<double, 3> conversion_rates(double mu)
{
   Loop_library::set(1);

   const softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1e-5);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 1);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, 1);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

   MSSMCPV_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, make_input(mu));
   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   const auto model = spectrum_generator.get_model();
   std::array<double, 3> res{};
   const std::array<Nucleus, 3> nuclei{Nucleus::Al, Nucleus::Ti, Nucleus::Au};
   for (std::size_t i = 0; i < nuclei.size(); i++) {
      res.at(i) =
         MSSMCPV_l_to_l_conversion::calculate_FeFe_forAll_1loop(
            1, 0, nuclei.at(i), model, LToLConversion_settings{}, qedqcd)(0)
            .real();
   }

   return res;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_mu_e_conversion)
{
   const std::array<const char*, 3> names{"Al", "Ti", "Au"};

   // mu and CR(mu -> e) in Al, Ti, Au from the SPheno/FlavorKit Wilson
   // coefficients
   const std::array<std::pair<double, std::array<double, 3>>, 2> points{{
      {500., {1.92e-13, 2.95e-13, 2.34e-13}},
      {800., {1.11e-13, 1.78e-13, 1.52e-13}},
   }};

   for (const auto& [mu, cr_ref] : points) {
      const auto cr = conversion_rates(mu);
      for (std::size_t i = 0; i < cr.size(); i++) {
         BOOST_TEST_MESSAGE("mu = " << mu << ": CR(mu -> e, " << names.at(i)
                                    << ") = " << cr.at(i) << ", reference "
                                    << cr_ref.at(i));
         BOOST_CHECK_CLOSE_FRACTION(cr.at(i), cr_ref.at(i), 0.15);
      }
   }
}
