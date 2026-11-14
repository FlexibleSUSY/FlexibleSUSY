#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_decay_amplitudes

#include <boost/test/unit_test.hpp>

#include "decays/decay_amplitudes.hpp"

#include <complex>

using namespace flexiblesusy;

namespace {

double square_FFS(double m_decay, double m_fermion, double m_scalar,
                  std::complex<double> left, std::complex<double> right)
{
   return Decay_amplitude_FFS{m_decay, m_fermion, m_scalar, left, right}.square();
}

} // anonymous namespace

// The spin-averaged |M|^2 of F(M) -> F(m) S(ms) with the vertex L P_L + R P_R is
//
//    (M^2 + m^2 - ms^2)/2 (|L|^2 + |R|^2) + 2 M m Re(L R^*)
//
// Same as SquareAmp_FtoFS * helfactor in SARAH 4.15.5, and consistent with
// FermionToFermionScalar in SPheno 4.0.7.

BOOST_AUTO_TEST_CASE(test_FFS_interference_term)
{
   BOOST_CHECK_CLOSE_FRACTION(square_FFS(500., 100., 50., 1., 1.), 357500., 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(square_FFS(500., 100., 50., 0.3, -0.7), 53675., 1e-14);
}

BOOST_AUTO_TEST_CASE(test_FFS_massless_scalar_keeps_interference_term)
{
   BOOST_CHECK_CLOSE_FRACTION(square_FFS(500., 100., 0., 1., 1.), 360000., 1e-14);
}

BOOST_AUTO_TEST_CASE(test_FFS_no_interference_with_relative_phase)
{
   const std::complex<double> i(0., 1.);

   BOOST_CHECK_CLOSE_FRACTION(square_FFS(500., 100., 50., 1., i), 257500., 1e-14);
}

// Decays of a vector boson V(M), averaged over its 3 polarizations.  The
// reference values follow from VectorBosonToTwoScalars,
// VectorBosonToTwoFermionsC, VectorBosonToScalarAndVectorBoson and
// VectorBosonToTwoVectorBosons in SPheno 4.0.7, whose widths are
// sqrt(lambda)/(48 Pi M^3) times the polarization-summed |M|^2.

BOOST_AUTO_TEST_CASE(test_VSS_square)
{
   Decay_amplitude_VSS amp;
   amp.m_decay = 3000.;
   amp.m_scalar_1 = 500.;
   amp.m_scalar_2 = 300.;
   amp.form_factor = 0.4;

   // |g|^2 lambda(M^2, m1^2, m2^2) / (3 M^2)
   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 443885.037037037, 1e-14);
}

BOOST_AUTO_TEST_CASE(test_VFF_square)
{
   Decay_amplitude_VFF amp;
   amp.m_decay = 3000.;
   amp.m_fermion_1 = 173.;
   amp.m_fermion_2 = 4.;
   amp.form_factor_left = 0.3;
   amp.form_factor_right = std::complex<double>(-0.2, 0.1);

   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 838431.84702742, 1e-14);

   // no interference with a relative phase of i
   amp.form_factor_right = std::complex<double>(0., 0.3);
   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 1078197.3347495398, 1e-14);

   // massless fermions: 2 M^2 (|L|^2 + |R|^2) / 3
   amp.m_fermion_1 = 0.;
   amp.m_fermion_2 = 0.;
   amp.form_factor_right = -0.2;
   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 780000., 1e-14);
}

BOOST_AUTO_TEST_CASE(test_VVS_square)
{
   Decay_amplitude_VVS amp;
   amp.m_decay = 3000.;
   amp.m_vector = 91.;
   amp.m_scalar = 125.;
   amp.form_factor = 50.;

   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 227719.22656039931, 1e-14);

   amp.m_vector = 0.;
   BOOST_CHECK_EQUAL(amp.square(), 0.);
}

BOOST_AUTO_TEST_CASE(test_VVV_square)
{
   Decay_amplitude_VVV amp;
   amp.m_decay = 3000.;
   amp.m_vector_1 = 80.;
   amp.m_vector_2 = 90.;
   amp.form_factor = 0.01;

   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 118692738.56017244, 1e-13);

   amp.m_vector_2 = 0.;
   BOOST_CHECK_EQUAL(amp.square(), 0.);
}
