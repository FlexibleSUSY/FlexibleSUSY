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
 * The anomalous magnetic moment must be finite (not NaN) for complex
 * parameters.  With complex couplings the imaginary part of the form
 * factor F2 cancels only up to rounding, which must not be mistaken
 * for a non-zero imaginary part.
 *
 * The test calculates a_mu at two physically equivalent points: A with
 * a phase on M_i and the A-terms, and B obtained from A by an R
 * rotation (M_i, T_f -> M_i, T_f exp(-i phi), mu -> mu exp(i phi)),
 * which moves the phase to mu.  Both must give the same a_mu.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMCPV_amm_complex

#include <boost/test/unit_test.hpp>

#include "MSSMCPV_amm.hpp"
#include "MSSMCPV_input_parameters.hpp"
#include "MSSMCPV_two_scale_spectrum_generator.hpp"
#include "cxx_qft/MSSMCPV_fields.hpp"

#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "spectrum_generator_settings.hpp"

#include <cmath>
#include <complex>
#include <numbers>

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

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_MSSMCPV_amm_complex_parameters)
{
   Loop_library::set(0);

   const double phi = std::numbers::pi / 3.;

   const softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1e-5);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 1);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, 1);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 1);

   MSSMCPV_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, make_input(phi));

   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   // same Lagrangian parameters, two phase conventions
   auto mA = spectrum_generator.get_model();
   auto mB = mA;
   rotate(mB, phi);

   for (auto* m : {&mA, &mB}) {
      m->calculate_DRbar_masses();
      m->calculate_pole_masses();
   }

   using MSSMCPV_cxx_diagrams::fields::Fe;

   const double amuA = MSSMCPV_amm::calculate_amm<Fe>(mA, qedqcd, settings, 1);
   const double amuB = MSSMCPV_amm::calculate_amm<Fe>(mB, qedqcd, settings, 1);

   BOOST_TEST_MESSAGE("a_mu: " << amuA << " vs " << amuB);

   BOOST_CHECK(std::isfinite(amuA));
   BOOST_CHECK(std::isfinite(amuB));
   BOOST_CHECK(std::abs(amuA) > 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(amuA, amuB, 1e-10);
}
