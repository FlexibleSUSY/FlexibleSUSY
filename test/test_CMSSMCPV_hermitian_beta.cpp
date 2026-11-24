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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCPV_hermitian_beta

#include <boost/test/unit_test.hpp>

#include "CMSSMCPV_input_parameters.hpp"
#include "CMSSMCPV_soft_parameters.hpp"
#include "CMSSMCPV_two_scale_spectrum_generator.hpp"
#include "lowe.h"
#include "spectrum_generator_settings.hpp"

#include <Eigen/Core>

#include <complex>
#include <string>

using namespace flexiblesusy;

namespace
{

using Matrix3c = Eigen::Matrix<std::complex<double>, 3, 3>;

Matrix3c make_complex_matrix(double scale, double phase)
{
   Matrix3c m;

   for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
         m(i, k) = std::polar(scale * (1. + i + 0.3 * k) / (1. + 2 * (i + k)),
                              phase * (1 + i - k));
      }
   }

   return m;
}

Matrix3c make_hermitian_matrix(double phase)
{
   const Matrix3c a = make_complex_matrix(1., phase);
   return a * a.adjoint();
}

/// checks that m is exactly hermitian
void check_exactly_hermitian(const Matrix3c& m, const std::string& name)
{
   BOOST_TEST_MESSAGE("checking " << name);

   for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_EQUAL(std::imag(m(i, i)), 0.);
      for (int k = i + 1; k < 3; ++k) {
         BOOST_CHECK_EQUAL(std::real(m(i, k)), std::real(m(k, i)));
         BOOST_CHECK_EQUAL(std::imag(m(i, k)), -std::imag(m(k, i)));
      }
   }
}

} // anonymous namespace

/// The beta functions of the hermitian soft scalar masses must be
/// exactly hermitian. Otherwise rounding errors make the imaginary
/// parts of their diagonal elements non-zero, and the RG integrator
/// can fail to integrate this noise.
BOOST_AUTO_TEST_CASE(test_CMSSMCPV_soft_mass_beta_functions_hermitian)
{
   CMSSMCPV_soft_parameters p;
   p.set_scale(1000.);
   p.set_g1(0.47);
   p.set_g2(0.63);
   p.set_g3(1.05);
   p.set_vd(24.);
   p.set_vu(243.);
   p.set_Mu(std::polar(500., 0.3));
   p.set_Yu(make_complex_matrix(0.9, 0.1));
   p.set_Yd(make_complex_matrix(0.2, 0.4));
   p.set_Ye(make_complex_matrix(0.1, 0.7));
   p.set_TYu(make_complex_matrix(800., 0.2));
   p.set_TYd(make_complex_matrix(200., 0.5));
   p.set_TYe(make_complex_matrix(100., 0.9));
   p.set_BMu(std::polar(1e5, 0.2));
   p.set_mq2(1e6 * make_hermitian_matrix(0.3));
   p.set_ml2(4e5 * make_hermitian_matrix(0.6));
   p.set_mHd2(1e6);
   p.set_mHu2(-1e6);
   p.set_md2(9e5 * make_hermitian_matrix(0.8));
   p.set_mu2(8e5 * make_hermitian_matrix(1.1));
   p.set_me2(3e5 * make_hermitian_matrix(1.3));
   p.set_MassB(std::polar(300., 0.1));
   p.set_MassWB(std::polar(600., 0.2));
   p.set_MassG(std::polar(1800., 0.3));

   for (int loops = 1; loops <= 2; ++loops) {
      p.set_loops(loops);
      const CMSSMCPV_soft_parameters beta = p.calc_beta();

      check_exactly_hermitian(beta.get_mq2(), "beta_mq2");
      check_exactly_hermitian(beta.get_ml2(), "beta_ml2");
      check_exactly_hermitian(beta.get_md2(), "beta_md2");
      check_exactly_hermitian(beta.get_mu2(), "beta_mu2");
      check_exactly_hermitian(beta.get_me2(), "beta_me2");
   }
}

/// This point failed with "non-perturbative running" in Im(mu2(2,2))
/// before the beta functions of the soft scalar masses were made
/// exactly hermitian.
BOOST_AUTO_TEST_CASE(test_CMSSMCPV_running_with_vev_phase)
{
   CMSSMCPV_input_parameters input;
   input.m0 = 1400.;
   input.m12 = 3500.;
   input.TanBeta = 10.;
   input.CosPhiMu = 1.;
   input.Azero = -1400.;
   input.etaInput = 0.1;

   const softsusy::QedQcd qedqcd;

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 0.);

   CMSSMCPV_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   BOOST_TEST_MESSAGE(spectrum_generator.get_problems());
   BOOST_CHECK(!spectrum_generator.get_problems().have_problem());

   const auto& model = spectrum_generator.get_model();

   BOOST_CHECK(!model.get_problems().no_perturbative());
   BOOST_CHECK(!model.get_problems().have_non_perturbative_parameter());

   for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_EQUAL(std::imag(model.get_mu2(i, i)), 0.);
   }
}
