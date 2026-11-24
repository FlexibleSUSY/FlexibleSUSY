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
 * masses unchanged.  The chargino and neutralino pole masses must
 * therefore be the same before and after the rotation.  This tests the
 * complex conjugations of the fermion self-energies in the 1-loop mass
 * matrices.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMCPV_rephasing_invariance

#include <boost/test/unit_test.hpp>

#include "MSSMCPV_input_parameters.hpp"
#include "MSSMCPV_two_scale_spectrum_generator.hpp"

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

BOOST_AUTO_TEST_CASE(test_MSSMCPV_chargino_neutralino_rephasing)
{
   Loop_library::set(0);

   const double phi = std::numbers::pi / 3.;

   const softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1e-5);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 1);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, 1);
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 0);
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

   const auto& pA = mA.get_physical();
   const auto& pB = mB.get_physical();

   // the phases must matter, otherwise the test is void
   BOOST_CHECK_GT(std::abs(mA.get_MassWB().imag()), 100.);
   BOOST_CHECK_GT(std::abs(mB.get_Mu().imag()), 100.);

   for (int i = 0; i < 2; i++) {
      BOOST_TEST_MESSAGE("MCha(" << i << "): " << pA.MCha(i) << " vs "
                                 << pB.MCha(i));
      BOOST_CHECK_CLOSE_FRACTION(pA.MCha(i), pB.MCha(i), 1e-9);
   }

   for (int i = 0; i < 4; i++) {
      BOOST_TEST_MESSAGE("MChi(" << i << "): " << pA.MChi(i) << " vs "
                                 << pB.MChi(i));
      BOOST_CHECK_CLOSE_FRACTION(pA.MChi(i), pB.MChi(i), 1e-9);
   }
}
