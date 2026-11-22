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
 * In U(1)xU(1) models with kinetic mixing the photon is a mixture of
 * B, W3 and B' as soon as the gauge couplings are not in the
 * triangular basis (gX1 != 0), e.g. after running away from the scale
 * where gX1 = 0 is imposed.  The vertices must then be written with
 * the photon-Z-Z' mixing matrix, not with two mixing angles that
 * assume gX1 = 0.  Checks: the photon doesn't couple to neutral
 * fermions, and its coupling to charged leptons is vectorial and the
 * same for all of them (including the vector-like ones).
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_U1xMSSM3G_photon_couplings

#include <boost/test/unit_test.hpp>

#include "U1xMSSM3G_input_parameters.hpp"
#include "U1xMSSM3G_slha_io.hpp"
#include "U1xMSSM3G_two_scale_spectrum_generator.hpp"
#include "cxx_qft/U1xMSSM3G_context_base.hpp"
#include "cxx_qft/U1xMSSM3G_fields.hpp"
#include "cxx_qft/U1xMSSM3G_vertices.hpp"

#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "spectrum_generator_settings.hpp"

#include <array>
#include <cmath>

using namespace flexiblesusy;
using namespace flexiblesusy::U1xMSSM3G_cxx_diagrams;

BOOST_AUTO_TEST_CASE(test_photon_couplings_away_from_input_scale)
{
   Loop_library::set(0);

   U1xMSSM3G_slha_io slha_io;
   slha_io.read_from_file("model_files/U1xMSSM3G/LesHouches.in.U1xMSSM3G");

   softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;
   U1xMSSM3G_input_parameters input;
   slha_io.fill(settings);
   slha_io.fill(qedqcd);
   slha_io.fill(input);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 0);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 0);

   U1xMSSM3G_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   // gX1 = 0 is imposed at the SUSY scale; run away from it
   auto m = spectrum_generator.get_model();
   m.run_to(500.);
   m.calculate_DRbar_masses();
   BOOST_REQUIRE_GT(std::abs(m.get_gX1()), 1e-3);

   const context_base context{&m};

   // reference: photon coupling to the electron
   const auto vee =
      Vertex<fields::bar<fields::Fe>::type, fields::Fe, fields::VP>::evaluate(
         std::array<int, 2>{0, 0}, context);
   const double e = std::abs(vee.left());
   BOOST_REQUIRE_GT(e, 0.2);
   BOOST_CHECK_CLOSE_FRACTION(std::abs(vee.right()), e, 1e-10);

   // all charged leptons, including the vector-like ones Fe(4..6)
   for (int i = 0; i < 6; i++) {
      const auto v =
         Vertex<fields::bar<fields::Fe>::type, fields::Fe,
                fields::VP>::evaluate(std::array<int, 2>{i, i}, context);
      BOOST_CHECK_CLOSE_FRACTION(std::abs(v.left()), e, 1e-10);
      BOOST_CHECK_CLOSE_FRACTION(std::abs(v.right()), e, 1e-10);
   }

   // no coupling of the photon to neutral fermions
   for (int i = 0; i < 8; i++) {
      for (int k = 0; k < 8; k++) {
         const auto v = Vertex<fields::Chi, fields::Chi, fields::VP>::evaluate(
            std::array<int, 2>{i, k}, context);
         BOOST_CHECK_SMALL(std::abs(v.left()), 1e-10 * e);
         BOOST_CHECK_SMALL(std::abs(v.right()), 1e-10 * e);
      }
   }
   for (int i = 0; i < 6; i++) {
      const auto v =
         Vertex<fields::bar<fields::Fv>::type, fields::Fv,
                fields::VP>::evaluate(std::array<int, 2>{i, i}, context);
      BOOST_CHECK_SMALL(std::abs(v.left()), 1e-10 * e);
      BOOST_CHECK_SMALL(std::abs(v.right()), 1e-10 * e);
   }
}
