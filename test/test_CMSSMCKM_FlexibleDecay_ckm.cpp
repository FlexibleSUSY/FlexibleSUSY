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
 * The decays must use the CKM matrix of the input in models that set
 * the quark Yukawa couplings with it (CMSSMCKM).
 *
 * The decoupling scheme matches the running quark masses of the model to
 * those of the SM.  The CKM matrix changes the Yukawa matrices (Yu =
 * sqrt(2) m_u V_CKM / vu), but not their singular values, so the masses
 * must stay the same; the mixing matrices must reproduce the CKM matrix.
 *
 * In the decays this is checked with the charged Higgs decays into
 * quarks.  At tan(beta) = 10
 * the coupling of H- -> b ubar_i is dominated by m_b tan(beta) V_ib, so
 *
 *    Gamma(H- -> b ubar) / Gamma(H- -> b cbar) = |V_ub|^2 / |V_cb|^2
 *
 * up to corrections of order (m_c/(m_b tan^2(beta)))^2 and phase space
 * effects of order m_c^2/m_H^2, both below 1e-4.  With a unit CKM matrix
 * H- -> b ubar is zero.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCKM_FlexibleDecay_ckm

#include <boost/test/unit_test.hpp>

#include "CMSSMCKM_input_parameters.hpp"
#include "CMSSMCKM_mass_eigenstates_decoupling_scheme.hpp"
#include "CMSSMCKM_slha_io.hpp"
#include "CMSSMCKM_two_scale_spectrum_generator.hpp"
#include "decays/CMSSMCKM_decays.hpp"

#include "decays/flexibledecay_settings.hpp"
#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "standard_model.hpp"

#include <Eigen/Core>

#include <cmath>
#include <complex>

using namespace flexiblesusy;

namespace
{

struct Point {
   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   CMSSMCKM_input_parameters input;
   CMSSMCKM_mass_eigenstates model;
};

/// calculates the spectrum of the shipped CMSSMCKM input
Point calculate_point()
{
   Loop_library::set(-1);

   CMSSMCKM_slha_io slha_io;
   slha_io.read_from_file("model_files/CMSSMCKM/LesHouches.in.CMSSMCKM");

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   Spectrum_generator_settings settings;
   CMSSMCKM_input_parameters input;
   slha_io.fill(settings);
   slha_io.fill(qedqcd);
   slha_io.fill(physical_input);
   slha_io.fill(input);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

   CMSSMCKM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   return {qedqcd, physical_input, input, spectrum_generator.get_model()};
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_decoupling_scheme_with_ckm)
{
   const auto point = calculate_point();
   const auto ckm = point.qedqcd.get_complex_ckm();

   standard_model::Standard_model sm;
   sm.initialise_from_input(point.qedqcd);
   sm.solve_ewsb_tree_level();
   sm.calculate_DRbar_masses();

   CMSSMCKM_mass_eigenstates_decoupling_scheme dm_ckm(point.input);
   dm_ckm.fill_from(point.model);
   dm_ckm.fill_from(sm, ckm);

   CMSSMCKM_mass_eigenstates_decoupling_scheme dm_unit(point.input);
   dm_unit.fill_from(point.model);
   dm_unit.fill_from(sm);

   // running quark masses are matched to the SM, with and without CKM
   for (int i = 0; i < 3; i++) {
      BOOST_CHECK_CLOSE_FRACTION(dm_ckm.get_MFu(i), sm.get_MFu(i), 1e-10);
      BOOST_CHECK_CLOSE_FRACTION(dm_ckm.get_MFd(i), sm.get_MFd(i), 1e-10);
      BOOST_CHECK_CLOSE_FRACTION(dm_ckm.get_MFu(i), dm_unit.get_MFu(i), 1e-10);
      BOOST_CHECK_CLOSE_FRACTION(dm_ckm.get_MFd(i), dm_unit.get_MFd(i), 1e-10);
   }

   // the mixing matrices reproduce the CKM matrix (up to the phases of
   // the quark fields)
   const Eigen::Matrix<std::complex<double>, 3, 3> ckm_model =
      dm_ckm.get_ZUL() * dm_ckm.get_ZDL().adjoint();
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_SMALL(std::abs(ckm_model(i, k)) - std::abs(ckm(i, k)),
                           1e-10);
      }
   }

   // the CP-violating phase: the Jarlskog invariant doesn't depend on
   // the phases of the quark fields, but changes sign under V -> V^*,
   // which the moduli don't see
   const auto jarlskog =
      [](const Eigen::Matrix<std::complex<double>, 3, 3>& V) {
         return std::imag(V(0, 1) * V(1, 2) * std::conj(V(0, 2)) *
                          std::conj(V(1, 1)));
      };
   BOOST_TEST_MESSAGE("Jarlskog invariant: model "
                      << jarlskog(ckm_model) << ", input " << jarlskog(ckm));
   BOOST_REQUIRE_GT(std::abs(jarlskog(ckm)), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(jarlskog(ckm_model), jarlskog(ckm), 1e-8);
}

BOOST_AUTO_TEST_CASE(test_charged_Higgs_decays_with_ckm)
{
   const auto point = calculate_point();

   const FlexibleDecay_settings flexibledecay_settings;
   CMSSMCKM_decays decays(point.model, point.qedqcd, point.physical_input,
                          flexibledecay_settings);
   decays.calculate_decays();
   BOOST_REQUIRE(!decays.get_problems().have_problem());

   // physical charged Higgs boson (index 0 is the Goldstone boson)
   const auto& hpm_decays = decays.get_Hpm_decays(1);
   const double width_bu = hpm_decays.get_decay({-2, 5}).get_width();
   const double width_bc = hpm_decays.get_decay({-4, 5}).get_width();

   const auto ckm = point.qedqcd.get_complex_ckm();
   const double expected = std::norm(ckm(0, 2)) / std::norm(ckm(1, 2));

   BOOST_TEST_MESSAGE("Gamma(H- -> b ubar)/Gamma(H- -> b cbar) = "
                      << width_bu / width_bc
                      << ", |V_ub/V_cb|^2 = " << expected);

   BOOST_REQUIRE_GT(width_bc, 0.);
   BOOST_CHECK_CLOSE_FRACTION(width_bu / width_bc, expected, 1e-3);
}
