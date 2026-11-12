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
