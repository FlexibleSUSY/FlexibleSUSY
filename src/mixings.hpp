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

#ifndef MIXINGS_H
#define MIXINGS_H

#include "numerics2.hpp"
#include <Eigen/Core>
#include <complex>
#include <type_traits>

namespace flexiblesusy {

/**
 * Imaginary parts of mixing matrix elements up to this size are taken
 * as rounding errors.  The SLHA output prints the O(1) elements with 8
 * significant digits, so smaller imaginary parts are not visible.
 */
constexpr double mixing_imag_tolerance = 1e-10;

/**
 * @param z mixing matrix
 * @return true if an element of @p z has an imaginary part larger than
 * rounding errors (always false for a real matrix)
 */
template <class Derived>
bool has_imaginary_part(const Eigen::DenseBase<Derived>& z)
{
   using Scalar = typename Derived::Scalar;
   if constexpr (std::is_same_v<Scalar, std::complex<double>>) {
      return z.size() > 0 &&
             z.imag().cwiseAbs().maxCoeff() > mixing_imag_tolerance;
   } else {
      return false;
   }
}

/**
 * Fixes the phases of the rows of a scalar mixing matrix for the SLHA
 * output: each row is multiplied by a phase such that its element with
 * the largest absolute value is real and positive.
 *
 * The phase of each row of a complex scalar mixing matrix is arbitrary
 * (a phase of the mass eigenstate field) and comes out of the
 * diagonalization at random.  Without CP violation in the mass matrix,
 * every row is then real (up to rounding errors), and the real part
 * written to the SLHA output is the whole matrix.  A real matrix is
 * left unchanged.
 *
 * @param z mixing matrix
 */
template <class Derived>
void convert_scalar_mixings_to_slha(Eigen::MatrixBase<Derived>& z)
{
   using Scalar = typename Derived::Scalar;
   if constexpr (std::is_same_v<Scalar, std::complex<double>>) {
      for (Eigen::Index i = 0; i < z.rows(); i++) {
         Eigen::Index k = 0;
         const double max_abs = z.row(i).cwiseAbs().maxCoeff(&k);
         if (max_abs > 0.) {
            z.row(i) *= std::conj(z(i, k)) / max_abs;
         }
      }
   }
}

namespace detail
{

/// assigns @p src to @p dst, taking the real part if @p dst is real
template <class DerivedDst, class DerivedSrc>
void assign_mixing(Eigen::MatrixBase<DerivedDst>& dst,
                   const Eigen::MatrixBase<DerivedSrc>& src)
{
   using Scalar = typename DerivedDst::Scalar;
   if constexpr (std::is_same_v<Scalar, std::complex<double>>) {
      dst = src;
   } else {
      dst = src.real();
   }
}

} // namespace detail

/**
 * Rotates the mixing matrix of 2n sfermions from the interaction basis
 * (f_L, f_R) of the model to the super-CKM (super-PMNS) basis of the
 * SLHA output, in which the Yukawa couplings are diagonal.
 *
 * The rotation matrices are those that convert the Yukawa couplings,
 * trilinear couplings and soft-breaking squared masses to SLHA
 * convention, e.g. m_Q^2 -> vL m_Q^2 vL^dagger and
 * m_d^2 -> vR^* m_d^2 vR^T, so the sfermion fields transform as
 * f_L -> vL f_L and f_R -> vR^* f_R:
 *
 *    z -> z * diag(vL^dagger, vR^T)
 *
 * A real @p z (in models with real parameters, where the rotations are
 * real) stays real.
 *
 * @param z sfermion mixing matrix (2n x 2n)
 * @param vL left-handed fermion mixing matrix in SLHA convention (n x n)
 * @param vR right-handed fermion mixing matrix in SLHA convention (n x n)
 */
template <class DerivedZ, class DerivedL, class DerivedR>
void convert_sfermion_mixings_to_slha(Eigen::MatrixBase<DerivedZ>& z,
                                      const Eigen::MatrixBase<DerivedL>& vL,
                                      const Eigen::MatrixBase<DerivedR>& vR)
{
   const Eigen::Index n = vL.rows();
   Eigen::MatrixXcd zc = z.template cast<std::complex<double>>();
   zc.leftCols(n) = (zc.leftCols(n) * vL.adjoint()).eval();
   zc.rightCols(n) = (zc.rightCols(n) * vR.transpose()).eval();
   detail::assign_mixing(z, zc);
}

/**
 * Rotates the mixing matrix of n left-handed sfermions (e.g. sneutrinos)
 * from the interaction basis to the super-CKM (super-PMNS) basis of the
 * SLHA output: z -> z * vL^dagger, see
 * convert_sfermion_mixings_to_slha(z, vL, vR).
 *
 * @param z sfermion mixing matrix (n x n)
 * @param vL left-handed fermion mixing matrix in SLHA convention (n x n)
 */
template <class DerivedZ, class DerivedL>
void convert_sfermion_mixings_to_slha(Eigen::MatrixBase<DerivedZ>& z,
                                      const Eigen::MatrixBase<DerivedL>& vL)
{
   const Eigen::MatrixXcd zc =
      z.template cast<std::complex<double>>() * vL.adjoint();
   detail::assign_mixing(z, zc);
}

void convert_symmetric_fermion_mixings_to_slha(double&,
                                               Eigen::Matrix<double, 1, 1>&);

template<int N>
void convert_symmetric_fermion_mixings_to_slha(Eigen::Array<double, N, 1>&,
                                               Eigen::Matrix<double, N, N>&)
{
}

void convert_symmetric_fermion_mixings_to_slha(double&,
                                               Eigen::Matrix<std::complex<double>, 1, 1>&);

/**
 * Converts the given vector of masses and the corresponding (complex)
 * mixing matrix to SLHA convention: Matrix rows whose imaginary part
 * is larger than their real part are multiplied by i and the
 * corresponding mass eigenvalue is multiplied by -1.  As a result the
 * mixing matrix will be real (up to rounding errors) and the mass
 * eigenvalues might be positive or negative.  It is assumed that these
 * mixings result from diagonalizing a symmetric fermion mass matrix in
 * the convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263.
 *
 * A row is compared as a whole, so that a real row with imaginary
 * parts from rounding errors is not rotated (which would make it
 * purely imaginary).
 *
 * If the mass matrix violates CP, the mixing matrix can't be made
 * real.  Then the masses and the mixing matrix are left unchanged:
 * the masses stay positive and the mixing matrix complex, as SLHA2
 * (arXiv:0801.0045) asks for the CP-violating case.
 *
 * @param m vector of masses
 * @param z mixing matrix
 */
template<int N>
void convert_symmetric_fermion_mixings_to_slha(Eigen::Array<double, N, 1>& m,
                                               Eigen::Matrix<std::complex<double>, N, N>& z)
{
   Eigen::Array<double, N, 1> m_slha(m);
   Eigen::Matrix<std::complex<double>, N, N> z_slha(z);

   for (int i = 0; i < N; i++) {
      if (z_slha.row(i).imag().squaredNorm() >
          z_slha.row(i).real().squaredNorm()) {
         z_slha.row(i) *= std::complex<double>(0.0, 1.0);
         m_slha(i) *= -1;
      }
   }

   if (!has_imaginary_part(z_slha)) {
      m = m_slha;
      z = z_slha;
   }
}

void convert_symmetric_fermion_mixings_to_hk(double&,
                                             Eigen::Matrix<double, 1, 1>&);

template<int N>
void convert_symmetric_fermion_mixings_to_hk(Eigen::Array<double, N, 1>&,
                                             Eigen::Matrix<double, N, N>&)
{
}

void convert_symmetric_fermion_mixings_to_hk(double&,
                                             Eigen::Matrix<std::complex<double>, 1, 1>&);

/**
 * Converts the given vector of masses and the corresponding (real)
 * mixing matrix to Haber-Kane convention (Phys. Rept. 117 (1985)
 * 75-263): Masses are positive and mixing matrices can be complex.
 *
 * @param m vector of masses
 * @param z mixing matrix
 */
template<int N>
void convert_symmetric_fermion_mixings_to_hk(Eigen::Array<double, N, 1>& m,
                                             Eigen::Matrix<std::complex<double>, N, N>& z)
{
   for (int i = 0; i < N; i++) {
      if (m(i) < 0.) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
      }
   }
}

} // namespace flexiblesusy

#endif
