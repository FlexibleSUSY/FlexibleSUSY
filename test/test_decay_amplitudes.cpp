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

// F(M) -> F(m) V(mV) with the vertex gamma^mu (L P_L + R P_R): the
// spin-averaged |M|^2 is the bracket of FermionToFermionVectorBoson in
// SPheno 4.0.5 (width = sqrt(lambda)/(16 Pi M^3) |M|^2),
//
//    (|L|^2 + |R|^2) ((M^2 - m^2)^2/mV^2 + M^2 + m^2 - 2 mV^2)/2
//    - 6 M m Re(L^* R)
BOOST_AUTO_TEST_CASE(test_FFV_square_massive)
{
   Decay_amplitude_FFV amp;
   amp.m_decay = 500.;
   amp.m_fermion = 100.;
   amp.m_vector = 80.;
   amp.form_factor_gam_left = {0.3, 0.1};
   amp.form_factor_gam_right = {-0.7, 0.2};

   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 2969868., 1e-14);
}

// F -> F gamma with an amplitude that fulfils the Ward identity: the
// p-form factors are fixed by k_mu Gamma^mu = 0.  The value follows from the
// trace with the polarization sum -g_{mu nu}, and is the mV -> 0 limit of
// the massive result.
BOOST_AUTO_TEST_CASE(test_FFV_square_massless_gauge_invariant)
{
   const double m_decay = 500.;
   const double m_fermion = 100.;
   const std::complex<double> left{0.3, 0.1};
   const std::complex<double> right{-0.2, 0.05};
   const double p_dot_k = 0.5 * (m_decay * m_decay - m_fermion * m_fermion);

   Decay_amplitude_FFV amp;
   amp.m_decay = m_decay;
   amp.m_fermion = m_fermion;
   amp.m_vector = 0.;
   amp.form_factor_gam_left = left;
   amp.form_factor_gam_right = right;
   amp.form_factor_p_1 = (left * m_fermion - right * m_decay) / p_dot_k;
   amp.form_factor_p_2 = (right * m_fermion - left * m_decay) / p_dot_k;

   BOOST_CHECK_CLOSE_FRACTION(amp.square(), 48050., 1e-13);

   // massive result for a tiny vector mass, with the p-form factors
   // fixed for that mass
   const double m_vector = 1e-2;
   const double p_dot_k_massive = p_dot_k + 0.5 * m_vector * m_vector;
   Decay_amplitude_FFV amp_massive = amp;
   amp_massive.m_vector = m_vector;
   amp_massive.form_factor_p_1 =
      (left * m_fermion - right * m_decay) / p_dot_k_massive;
   amp_massive.form_factor_p_2 =
      (right * m_fermion - left * m_decay) / p_dot_k_massive;

   BOOST_CHECK_CLOSE_FRACTION(amp_massive.square(), amp.square(), 1e-6);
}

// multiplying an amplitude by a factor multiplies |M|^2 by |factor|^2 and
// keeps the masses and the threshold for massless vectors
BOOST_AUTO_TEST_CASE(test_amplitude_times_factor)
{
   const std::complex<double> factor{0.6, -0.8}; // |factor| = 1
   const std::complex<double> factor2{1.5, 2.};  // |factor2|^2 = 6.25

   Decay_amplitude_SSV ssv;
   ssv.m_decay = 500.;
   ssv.m_scalar = 100.;
   ssv.m_vector = 80.;
   ssv.massless_vector_threshold = 1.;
   ssv.form_factor = 0.3;
   BOOST_CHECK_EQUAL((factor * ssv).massless_vector_threshold, 1.);
   BOOST_CHECK_EQUAL((ssv * factor).massless_vector_threshold, 1.);
   BOOST_CHECK_CLOSE_FRACTION((factor2 * ssv).square(), 6.25 * ssv.square(),
                              1e-14);

   const Decay_amplitude_FFS ffs{.m_decay = 500.,
                                 .m_fermion = 100.,
                                 .m_scalar = 50.,
                                 .form_factor_left = {0.3, 0.1},
                                 .form_factor_right = {-0.7, 0.2}};
   BOOST_CHECK_CLOSE_FRACTION((factor2 * ffs).square(), 6.25 * ffs.square(),
                              1e-14);
   BOOST_CHECK_CLOSE_FRACTION((ffs * factor2).square(), 6.25 * ffs.square(),
                              1e-14);

   Decay_amplitude_FFV ffv;
   ffv.m_decay = 500.;
   ffv.m_fermion = 100.;
   ffv.m_vector = 80.;
   ffv.massless_vector_threshold = 1.;
   ffv.form_factor_gam_left = {0.3, 0.1};
   ffv.form_factor_gam_right = {-0.7, 0.2};
   ffv.form_factor_p_1 = {1e-3, 2e-4};
   ffv.form_factor_p_2 = {-5e-4, 1e-3};
   BOOST_CHECK_EQUAL((factor * ffv).massless_vector_threshold, 1.);
   BOOST_CHECK_CLOSE_FRACTION((factor2 * ffv).square(), 6.25 * ffv.square(),
                              1e-14);
   BOOST_CHECK_CLOSE_FRACTION((ffv * factor2).square(), 6.25 * ffv.square(),
                              1e-14);
}
