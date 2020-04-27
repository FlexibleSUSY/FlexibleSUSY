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

#ifndef LINALG2_H
#define LINALG2_H

#include <cstdint>
#include <limits>
#include <cctype>
#include <cmath>
#include <complex>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>

namespace flexiblesusy {

#define MAX_(i, j) (((i) > (j)) ? (i) : (j))
#define MIN_(i, j) (((i) < (j)) ? (i) : (j))

template<class Real, class Scalar, int M, int N>
void svd_eigen
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M> *u,
 Eigen::Matrix<Scalar, N, N> *vh)
{
    Eigen::JacobiSVD<Eigen::Matrix<Scalar, M, N> >
	svd(m, (u ? Eigen::ComputeFullU : 0) | (vh ? Eigen::ComputeFullV : 0));
    s = svd.singularValues();
    if (u)  { *u  = svd.matrixU(); }
    if (vh) { *vh = svd.matrixV().adjoint(); }
}

template<class Real, class Scalar, int N>
void hermitian_eigen
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N> *z)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar,N,N> >
	es(m, z ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
    w = es.eigenvalues();
    if (z) { *z = es.eigenvectors(); }
}

/**
 * Template version of DDISNA from LAPACK.
 */
template<int M, int N, class Real>
void disna(const char& JOB, const Eigen::Array<Real, MIN_(M, N), 1>& D,
	   Eigen::Array<Real, MIN_(M, N), 1>& SEP, int& INFO)
{
//  -- LAPACK computational routine (version 3.4.0) --
//  -- LAPACK is a software package provided by Univ. of Tennessee,    --
//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
//     November 2011
//
//  =====================================================================

//     .. Parameters ..
      const Real         ZERO = 0;
//     .. Local Scalars ..
      bool               DECR, EIGEN, INCR, LEFT, RIGHT, SINGUL;
      int                I, K;
      Real               ANORM, EPS, NEWGAP, OLDGAP, SAFMIN, THRESH;
//
//     Test the input arguments
//
      INFO = 0;
      EIGEN = std::toupper(JOB) == 'E';
      LEFT  = std::toupper(JOB) == 'L';
      RIGHT = std::toupper(JOB) == 'R';
      SINGUL = LEFT || RIGHT;
      if (EIGEN) {
	 K = M;
      } else if (SINGUL) {
         K = MIN_(M, N);
      }
      if (!EIGEN && !SINGUL) {
         INFO = -1;
      } else if (M < 0) {
         INFO = -2;
      } else if (K < 0) {
         INFO = -3;
      } else {
         INCR = true;
         DECR = true;
         for (I = 0; I < K - 1; I++) {
            if (INCR) {
               INCR = INCR && D(I) <= D(I+1);
            }
            if (DECR) {
	       DECR = DECR && D(I) >= D(I+1);
            }
	 }
         if (SINGUL && K > 0) {
            if (INCR) {
               INCR = INCR && ZERO <= D(0);
            }
            if (DECR) {
               DECR = DECR && D(K-1) >= ZERO;
            }
         }
         if (!(INCR || DECR)) {
            INFO = -4;
         }
      }
      if (INFO != 0) {
         // CALL XERBLA( 'DDISNA', -INFO )
         return;
      }
//
//     Quick return if possible
//
      if (K == 0) {
         return;
      }
//
//     Compute reciprocal condition numbers
//
      if (K == 1) {
         SEP(0) = std::numeric_limits<Real>::max();
      } else {
         OLDGAP = std::fabs(D(1) - D(0));
         SEP(0) = OLDGAP;
         for (I = 1; I < K - 1; I++) {
            NEWGAP = std::fabs(D(I+1) - D(I));
            SEP(I) = std::min(OLDGAP, NEWGAP);
            OLDGAP = NEWGAP;
	 }
         SEP(K-1) = OLDGAP;
      }
      if (SINGUL) {
         if ((LEFT && M > N) || (RIGHT && M < N)) {
            if (INCR) {
               SEP( 0 ) = std::min(SEP( 0 ), D( 0 ));
            }
            if (DECR) {
               SEP(K-1) = std::min(SEP(K-1), D(K-1));
            }
         }
      }
//
//     Ensure that reciprocal condition numbers are not less than
//     threshold, in order to limit the size of the error bound
//
      // Note std::numeric_limits<double>::epsilon() == 2 * DLAMCH('E')
      // since  the former is the smallest eps such that 1.0 + eps > 1.0
      // while DLAMCH('E') is the smallest eps such that 1.0 - eps < 1.0
      EPS = std::numeric_limits<Real>::epsilon();
      SAFMIN = std::numeric_limits<Real>::min();
      ANORM = std::max(std::fabs(D(0)), std::fabs(D(K-1)));
      if (ANORM == ZERO) {
         THRESH = EPS;
      } else {
         THRESH = std::max(EPS*ANORM, SAFMIN);
      }
      for (I = 0; I < K; I++) {
	 SEP(I) = std::max(SEP(I), THRESH);
      }
}


template<class Real, class Scalar, int M, int N>
void svd_internal
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M> *u,
 Eigen::Matrix<Scalar, N, N> *vh)
{
    svd_eigen(m, s, u, vh);
}

template<class Real, class Scalar, int M, int N>
void svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M> *u  = 0,
 Eigen::Matrix<Scalar, N, N> *vh = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, MIN_(M, N), 1> *u_errbd = 0,
 Eigen::Array<Real, MIN_(M, N), 1> *v_errbd = 0)
{
    svd_internal(m, s, u, vh);

    // see http://www.netlib.org/lapack/lug/node96.html
    if (!s_errbd) { return; }
    const Real EPSMCH = std::numeric_limits<Real>::epsilon();
    *s_errbd = EPSMCH * s[0];

    Eigen::Array<Real, MIN_(M, N), 1> RCOND;
    int INFO;
    if (u_errbd) {
	disna<M, N>('L', s, RCOND, INFO);
	u_errbd->fill(*s_errbd);
	*u_errbd /= RCOND;
    }
    if (v_errbd) {
	disna<M, N>('R', s, RCOND, INFO);
	v_errbd->fill(*s_errbd);
	*v_errbd /= RCOND;
    }
}

/**
 * Singular value decomposition of M-by-N matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u * sigma * vh    // LAPACK convention
 *
 * and `(s >= 0).all()`.  Elements of s are in descending order.  The
 * above decomposition can be put in the form
 *
 *     m == u * s.matrix().asDiagonal() * vh
 *
 * if `M == N`.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m, u, and vh
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M unitary matrix
 * @param[out] vh     N-by-N unitary matrix
 */
template<class Real, class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh)
{
    svd_errbd(m, s, &u, &vh);
}

/**
 * Same as svd(m, s, u, vh) except that an approximate error bound for
 * the singular values is returned as well.  The error bound is
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of svd(m, s, u, vh) for the other parameters.
 */
template<class Real, class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 Real& s_errbd)
{
    svd_errbd(m, s, &u, &vh, &s_errbd);
}

/**
 * Same as svd(m, s, u, vh, s_errbd) except that approximate error
 * bounds for the singular vectors are returned as well.  The error
 * bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of svd(m, s, u, vh, s_errbd) for the other
 * parameters.
 */
template<class Real, class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 Real& s_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& v_errbd)
{
    svd_errbd(m, s, &u, &vh, &s_errbd, &u_errbd, &v_errbd);
}

/**
 * Returns singular values of M-by-N matrix m via s such that
 * `(s >= 0).all()`.  Elements of s are in descending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m, u, and vh
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 */
template<class Real, class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s)
{
    svd_errbd(m, s);
}

/**
 * Same as svd(m, s) except that an approximate error bound for the
 * singular values is returned as well.  The error bound is estimated
 * following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of svd(m, s) for the other parameters.
 */
template<class Real, class Scalar, int M, int N>
void svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Real& s_errbd)
{
    svd_errbd(m, s, 0, 0, &s_errbd);
}

// Eigen::SelfAdjointEigenSolver seems to be faster than ZHEEV of ATLAS

template<class Real, class Scalar, int N>
void diagonalize_hermitian_internal
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N> *z)
{
    hermitian_eigen(m, w, z);
}

template<class Real, class Scalar, int N>
void diagonalize_hermitian_errbd
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N> *z = 0,
 Real *w_errbd = 0,
 Eigen::Array<Real, N, 1> *z_errbd = 0)
{
    diagonalize_hermitian_internal(m, w, z);

    // see http://www.netlib.org/lapack/lug/node89.html
    if (!w_errbd) { return; }
    const Real EPSMCH = std::numeric_limits<Real>::epsilon();
    Real mnorm = std::max(std::abs(w[0]), std::abs(w[N-1]));
    *w_errbd = EPSMCH * mnorm;

    if (!z_errbd) { return; }
    Eigen::Array<Real, N, 1> RCONDZ;
    int INFO;
    disna<N, N>('E', w, RCONDZ, INFO);
    z_errbd->fill(*w_errbd);
    *z_errbd /= RCONDZ;
}

/**
 * Diagonalizes N-by-N hermitian matrix m so that
 *
 *     m == z * w.matrix().asDiagonal() * z.adjoint()
 *
 * Elements of w are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m and z
 * @tparam     N      number of rows and columns in m and z
 * @param[in]  m      N-by-N matrix to be diagonalized
 * @param[out] w      array of length N to contain eigenvalues
 * @param[out] z      N-by-N unitary matrix
 */
template<class Real, class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z)
{
    diagonalize_hermitian_errbd(m, w, &z);
}

/**
 * Same as diagonalize_hermitian(m, w, z) except that an approximate
 * error bound for the eigenvalues is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] w_errbd approximate error bound for the elements of w
 *
 * See the documentation of diagonalize_hermitian(m, w, z) for the
 * other parameters.
 */
template<class Real, class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 Real& w_errbd)
{
    diagonalize_hermitian_errbd(m, w, &z, &w_errbd);
}

/**
 * Same as diagonalize_hermitian(m, w, z, w_errbd) except that
 * approximate error bounds for the eigenvectors are returned as well.
 * The error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] z_errbd array of approximate error bounds for z
 *
 * See the documentation of diagonalize_hermitian(m, w, z, w_errbd)
 * for the other parameters.
 */
template<class Real, class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 Real& w_errbd,
 Eigen::Array<Real, N, 1>& z_errbd)
{
    diagonalize_hermitian_errbd(m, w, &z, &w_errbd, &z_errbd);
}

/**
 * Returns eigenvalues of N-by-N hermitian matrix m via w.
 * Elements of w are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m and z
 * @tparam     N      number of rows and columns in m and z
 * @param[in]  m      N-by-N matrix to be diagonalized
 * @param[out] w      array of length N to contain eigenvalues
 */
template<class Real, class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w)
{
    diagonalize_hermitian_errbd(m, w);
}

/**
 * Same as diagonalize_hermitian(m, w) except that an approximate
 * error bound for the eigenvalues is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] w_errbd approximate error bound for the elements of w
 *
 * See the documentation of diagonalize_hermitian(m, w) for the other
 * parameters.
 */
template<class Real, class Scalar, int N>
void diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Real& w_errbd)
{
    diagonalize_hermitian_errbd(m, w, 0, &w_errbd);
}

template<class Real, int N>
void diagonalize_symmetric_errbd
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N> *u = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, N, 1> *u_errbd = 0)
{
    if (!u) {
	svd_errbd(m, s, u, u, s_errbd, u_errbd);
	return;
    }
    Eigen::Matrix<std::complex<Real>, N, N> vh;
    svd_errbd(m, s, u, &vh, s_errbd, u_errbd);
    // see Eq. (5) of https://doi.org/10.1016/j.amc.2014.01.170
    *u *= (u->adjoint() * vh.transpose()).sqrt().eval();
}

/**
 * Diagonalizes N-by-N complex symmetric matrix m so that
 *
 *     m == u * s.matrix().asDiagonal() * u.transpose()
 *
 * and `(s >= 0).all()`.  Elements of s are in descending order.
 *
 * @tparam     Real type of real and imaginary parts
 * @tparam     N    number of rows and columns in m and u
 * @param[in]  m    N-by-N complex symmetric matrix to be decomposed
 * @param[out] s    array of length N to contain singular values
 * @param[out] u    N-by-N complex unitary matrix
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u)
{
    diagonalize_symmetric_errbd(m, s, &u);
}

/**
 * Same as diagonalize_symmetric(m, s, u) except that an approximate
 * error bound for the singular values is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of diagonalize_symmetric(m, s, u) for the
 * other parameters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd)
{
    diagonalize_symmetric_errbd(m, s, &u, &s_errbd);
}

/**
 * Same as diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned as
 * well.  The error bounds are estimated following the method
 * presented at http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of diagonalize_symmetric(m, s, u, s_errbd)
 * for the other parameters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd,
 Eigen::Array<Real, N, 1>& u_errbd)
{
    diagonalize_symmetric_errbd(m, s, &u, &s_errbd, &u_errbd);
}

/**
 * Returns singular values of N-by-N complex symmetric matrix m via s
 * such that `(s >= 0).all()`.  Elements of s are in descending order.
 *
 * @tparam     Real type of real and imaginary parts
 * @tparam     N    number of rows and columns in m and u
 * @param[in]  m    N-by-N complex symmetric matrix to be decomposed
 * @param[out] s    array of length N to contain singular values
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s)
{
    diagonalize_symmetric_errbd(m, s);
}

/**
 * Same as diagonalize_symmetric(m, s) except that an approximate
 * error bound for the singular values is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of diagonalize_symmetric(m, s) for the other
 * parameters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Real& s_errbd)
{
    diagonalize_symmetric_errbd(m, s, 0, &s_errbd);
}

template<class Real>
struct FlipSignOp {
    std::complex<Real> operator() (const std::complex<Real>& z) const {
	return z.real() < 0 ? std::complex<Real>(0,1) :
	    std::complex<Real>(1,0);
    }
};

template<class Real, int N>
void diagonalize_symmetric_errbd
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N> *u = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, N, 1> *u_errbd = 0)
{
    Eigen::Matrix<Real, N, N> z;
    diagonalize_hermitian_errbd(m, s, u ? &z : 0, s_errbd, u_errbd);
    // see http://forum.kde.org/viewtopic.php?f=74&t=62606
    if (u) {
        *u = z * s.template cast<std::complex<Real>>()
                    .unaryExpr(FlipSignOp<Real>())
                    .matrix()
                    .asDiagonal();
    }
    s = s.abs();
}

/**
 * Diagonalizes N-by-N real symmetric matrix m so that
 *
 *     m == u * s.matrix().asDiagonal() * u.transpose()
 *
 * and `(s >= 0).all()`.  Order of elements of s is *unspecified*.
 *
 * @tparam     Real type of real and imaginary parts
 * @tparam     N    number of rows and columns of m
 * @param[in]  m    N-by-N real symmetric matrix to be decomposed
 * @param[out] s    array of length N to contain singular values
 * @param[out] u    N-by-N complex unitary matrix
 *
 * @note Use diagonalize_hermitian() unless sign of `s[i]` matters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u)
{
    diagonalize_symmetric_errbd(m, s, &u);
}

/**
 * Same as diagonalize_symmetric(m, s, u) except that an approximate
 * error bound for the singular values is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of diagonalize_symmetric(m, s, u) for the
 * other parameters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd)
{
    diagonalize_symmetric_errbd(m, s, &u, &s_errbd);
}

/**
 * Same as diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned as
 * well.  The error bounds are estimated following the method
 * presented at http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of diagonalize_symmetric(m, s, u, s_errbd)
 * for the other parameters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd,
 Eigen::Array<Real, N, 1>& u_errbd)
{
    diagonalize_symmetric_errbd(m, s, &u, &s_errbd, &u_errbd);
}

/**
 * Returns singular values of N-by-N real symmetric matrix m via s
 * such that `(s >= 0).all()`.  Order of elements of s is
 * *unspecified*.
 *
 * @tparam     Real type of elements of m and s
 * @tparam     N    number of rows and columns of m
 * @param[in]  m    N-by-N real symmetric matrix to be decomposed
 * @param[out] s    array of length N to contain singular values
 *
 * @note Use diagonalize_hermitian() unless sign of `s[i]` matters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s)
{
    diagonalize_symmetric_errbd(m, s);
}

/**
 * Same as diagonalize_symmetric(m, s) except that an approximate
 * error bound for the singular values is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of diagonalize_symmetric(m, s) for the other
 * parameters.
 */
template<class Real, int N>
void diagonalize_symmetric
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Real& s_errbd)
{
    diagonalize_symmetric_errbd(m, s, 0, &s_errbd);
}

template<class Real, class Scalar, int M, int N>
void reorder_svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M> *u  = 0,
 Eigen::Matrix<Scalar, N, N> *vh = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, MIN_(M, N), 1> *u_errbd = 0,
 Eigen::Array<Real, MIN_(M, N), 1> *v_errbd = 0)
{
    svd_errbd(m, s, u, vh, s_errbd, u_errbd, v_errbd);
    s.reverseInPlace();
    if (u) {
	Eigen::PermutationMatrix<M> p;
	p.setIdentity();
	p.indices().template segment<MIN_(M, N)>(0).reverseInPlace();
	*u *= p;
    }
    if (vh) {
	Eigen::PermutationMatrix<N> p;
	p.setIdentity();
	p.indices().template segment<MIN_(M, N)>(0).reverseInPlace();
	vh->transpose() *= p;
    }
    if (u_errbd) { u_errbd->reverseInPlace(); }
    if (v_errbd) { v_errbd->reverseInPlace(); }
}

/**
 * Singular value decomposition of M-by-N matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u * sigma * vh    // LAPACK convention
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.  The
 * above decomposition can be put in the form
 *
 *     m == u * s.matrix().asDiagonal() * vh
 *
 * if `M == N`.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m, u, and vh
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M unitary matrix
 * @param[out] vh     N-by-N unitary matrix
 */
template<class Real, class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh)
{
    reorder_svd_errbd(m, s, &u, &vh);
}

/**
 * Same as reorder_svd(m, s, u, vh) except that an approximate error
 * bound for the singular values is returned as well.  The error bound
 * is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of reorder_svd(m, s, u, vh) for the other
 * parameters.
 */
template<class Real, class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 Real& s_errbd)
{
    reorder_svd_errbd(m, s, &u, &vh, &s_errbd);
}

/**
 * Same as reorder_svd(m, s, u, vh, s_errbd) except that approximate
 * error bounds for the singular vectors are returned as well.  The
 * error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of reorder_svd(m, s, u, vh, s_errbd) for the
 * other parameters.
 */
template<class Real, class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& vh,
 Real& s_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& v_errbd)
{
    reorder_svd_errbd(m, s, &u, &vh, &s_errbd, &u_errbd, &v_errbd);
}

/**
 * Returns singular values of M-by-N matrix m via s such that
 * `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m, u, and vh
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 */
template<class Real, class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s)
{
    reorder_svd_errbd(m, s);
}

/**
 * Same as reorder_svd(m, s) except that an approximate error bound
 * for the singular values is returned as well.  The error bound is
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of reorder_svd(m, s) for the other
 * parameters.
 */
template<class Real, class Scalar, int M, int N>
void reorder_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Real& s_errbd)
{
    reorder_svd_errbd(m, s, 0, 0, &s_errbd);
}

template<class Real, int N>
void reorder_diagonalize_symmetric_errbd
(const Eigen::Matrix<std::complex<Real>, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N> *u = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, N, 1> *u_errbd = 0)
{
    diagonalize_symmetric_errbd(m, s, u, s_errbd, u_errbd);
    s.reverseInPlace();
    if (u) { *u = u->rowwise().reverse().eval(); }
    if (u_errbd) { u_errbd->reverseInPlace(); }
}

template<class Real, int N>
void reorder_diagonalize_symmetric_errbd
(const Eigen::Matrix<Real, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N> *u = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, N, 1> *u_errbd = 0)
{
    diagonalize_symmetric_errbd(m, s, u, s_errbd, u_errbd);
    Eigen::PermutationMatrix<N> p;
    p.setIdentity();
    std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
              [&s] (int i, int j) { return s[i] < s[j]; });
#if EIGEN_VERSION_AT_LEAST(3,1,4)
    s.matrix().transpose() *= p;
    if (u_errbd) { u_errbd->matrix().transpose() *= p; }
#else
    Eigen::Map<Eigen::Matrix<Real, N, 1> >(s.data()).transpose() *= p;
    if (u_errbd) {
        Eigen::Map<Eigen::Matrix<Real, N, 1>>(u_errbd->data()).transpose() *= p;
    }
#endif
    if (u) { *u *= p; }
}

/**
 * Diagonalizes N-by-N symmetric matrix m so that
 *
 *     m == u * s.matrix().asDiagonal() * u.transpose()
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m
 * @tparam     N      number of rows and columns in m and u
 * @param[in]  m      N-by-N symmetric matrix to be decomposed
 * @param[out] s      array of length N to contain singular values
 * @param[out] u      N-by-N complex unitary matrix
 */
template<class Real, class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u)
{
    reorder_diagonalize_symmetric_errbd(m, s, &u);
}

/**
 * Same as reorder_diagonalize_symmetric(m, s, u) except that an
 * approximate error bound for the singular values is returned as
 * well.  The error bound is estimated following the method presented
 * at http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of reorder_diagonalize_symmetric(m, s, u) for
 * the other parameters.
 */
template<class Real, class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd)
{
    reorder_diagonalize_symmetric_errbd(m, s, &u, &s_errbd);
}

/**
 * Same as reorder_diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned as
 * well.  The error bounds are estimated following the method
 * presented at http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of reorder_diagonalize_symmetric(m, s, u,
 * s_errbd) for the other parameters.
 */
template<class Real, class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd,
 Eigen::Array<Real, N, 1>& u_errbd)
{
    reorder_diagonalize_symmetric_errbd(m, s, &u, &s_errbd, &u_errbd);
}

/**
 * Returns singular values of N-by-N symmetric matrix m via s such
 * that `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m
 * @tparam     N      number of rows and columns in m and u
 * @param[in]  m      N-by-N symmetric matrix to be decomposed
 * @param[out] s      array of length N to contain singular values
 */
template<class Real, class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s)
{
    reorder_diagonalize_symmetric_errbd(m, s);
}

/**
 * Same as reorder_diagonalize_symmetric(m, s) except that an
 * approximate error bound for the singular values is returned as
 * well.  The error bound is estimated following the method presented
 * at http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of reorder_diagonalize_symmetric(m, s) for
 * the other parameters.
 */
template<class Real, class Scalar, int N>
void reorder_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Real& s_errbd)
{
    reorder_diagonalize_symmetric_errbd(m, s, 0, &s_errbd);
}

template<class Real, class Scalar, int M, int N>
void fs_svd_errbd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M> *u = 0,
 Eigen::Matrix<Scalar, N, N> *v = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, MIN_(M, N), 1> *u_errbd = 0,
 Eigen::Array<Real, MIN_(M, N), 1> *v_errbd = 0)
{
    reorder_svd_errbd(m, s, u, v, s_errbd, u_errbd, v_errbd);
    if (u) { u->transposeInPlace(); }
}

/**
 * Singular value decomposition of M-by-N matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u.transpose() * sigma * v
 *     // convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.  The
 * above decomposition can be put in the form
 *
 *     m == u.transpose() * s.matrix().asDiagonal() * v
 *
 * if `M == N`.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m, u, and v
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M unitary matrix
 * @param[out] v      N-by-N unitary matrix
 */
template<class Real, class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v)
{
    fs_svd_errbd(m, s, &u, &v);
}

/**
 * Same as fs_svd(m, s, u, v) except that an approximate error bound
 * for the singular values is returned as well.  The error bound is
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_svd(m, s, u, v) for the other
 * parameters.
 */
template<class Real, class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v,
 Real& s_errbd)
{
    fs_svd_errbd(m, s, &u, &v, &s_errbd);
}

/**
 * Same as fs_svd(m, s, u, v, s_errbd) except that approximate error
 * bounds for the singular vectors are returned as well.  The error
 * bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of fs_svd(m, s, u, v, s_errbd) for the other
 * parameters.
 */
template<class Real, class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<Scalar, M, M>& u,
 Eigen::Matrix<Scalar, N, N>& v,
 Real& s_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& v_errbd)
{
    fs_svd_errbd(m, s, &u, &v, &s_errbd, &u_errbd, &v_errbd);
}

/**
 * Returns singular values of M-by-N matrix m via s such that
 * `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m, u, and v
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 */
template<class Real, class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s)
{
    fs_svd_errbd(m, s);
}

/**
 * Same as fs_svd(m, s) except that an approximate error bound for the
 * singular values is returned as well.  The error bound is estimated
 * following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_svd(m, s) for the other parameters.
 */
template<class Real, class Scalar, int M, int N>
void fs_svd
(const Eigen::Matrix<Scalar, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Real& s_errbd)
{
    fs_svd_errbd(m, s, 0, 0, &s_errbd);
}

/**
 * Singular value decomposition of M-by-N *real* matrix m such that
 *
 *     sigma.setZero(); sigma.diagonal() = s;
 *     m == u.transpose() * sigma * v
 *     // convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.  The
 * above decomposition can be put in the form
 *
 *     m == u.transpose() * s.matrix().asDiagonal() * v
 *
 * if `M == N`.
 *
 * @tparam     Real   type of real and imaginary parts
 * @tparam     M      number of rows in m
 * @tparam     N      number of columns in m
 * @param[in]  m      M-by-N *real* matrix to be decomposed
 * @param[out] s      array of length min(M,N) to contain singular values
 * @param[out] u      M-by-M *complex* unitary matrix
 * @param[out] v      N-by-N *complex* unitary matrix
 *
 * @note This is a convenience overload for the case where the type of
 * u and v (complex) differs from that of m (real).  Mathematically,
 * real u and v are enough to accommodate SVD of any real m.
 */
template<class Real, int M, int N>
void fs_svd
(const Eigen::Matrix<Real, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<std::complex<Real>, M, M>& u,
 Eigen::Matrix<std::complex<Real>, N, N>& v)
{
    fs_svd(m.template cast<std::complex<Real> >().eval(), s, u, v);
}

/**
 * Same as fs_svd(m, s, u, v) except that an approximate error bound
 * for the singular values is returned as well.  The error bound is
 * estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_svd(m, s, u, v) for the other
 * parameters.
 */
template<class Real, int M, int N>
void fs_svd
(const Eigen::Matrix<Real, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<std::complex<Real>, M, M>& u,
 Eigen::Matrix<std::complex<Real>, N, N>& v,
 Real& s_errbd)
{
    fs_svd(m.template cast<std::complex<Real> >().eval(), s, u, v, s_errbd);
}

/**
 * Same as fs_svd(m, s, u, v, s_errbd) except that approximate error
 * bounds for the singular vectors are returned as well.  The error
 * bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 * @param[out] v_errbd array of approximate error bounds for vh
 *
 * See the documentation of fs_svd(m, s, u, v, s_errbd) for the other
 * parameters.
 */
template<class Real, int M, int N>
void fs_svd
(const Eigen::Matrix<Real, M, N>& m,
 Eigen::Array<Real, MIN_(M, N), 1>& s,
 Eigen::Matrix<std::complex<Real>, M, M>& u,
 Eigen::Matrix<std::complex<Real>, N, N>& v,
 Real& s_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& u_errbd,
 Eigen::Array<Real, MIN_(M, N), 1>& v_errbd)
{
    fs_svd(m.template cast<std::complex<Real> >().eval(), s, u, v,
	   s_errbd, u_errbd, v_errbd);
}

template<class Real, class Scalar, int N>
void fs_diagonalize_symmetric_errbd
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N> *u = 0,
 Real *s_errbd = 0,
 Eigen::Array<Real, N, 1> *u_errbd = 0)
{
    reorder_diagonalize_symmetric_errbd(m, s, u, s_errbd, u_errbd);
    if (u) { u->transposeInPlace(); }
}

/**
 * Diagonalizes N-by-N symmetric matrix m so that
 *
 *     m == u.transpose() * s.matrix().asDiagonal() * u
 *     // convention of Haber and Kane, Phys. Rept. 117 (1985) 75-263
 *
 * and `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m
 * @tparam     N      number of rows and columns in m and u
 * @param[in]  m      N-by-N symmetric matrix to be decomposed
 * @param[out] s      array of length N to contain singular values
 * @param[out] u      N-by-N complex unitary matrix
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u)
{
    fs_diagonalize_symmetric_errbd(m, s, &u);
}

/**
 * Same as fs_diagonalize_symmetric(m, s, u) except that an
 * approximate error bound for the singular values is returned as
 * well.  The error bound is estimated following the method presented
 * at http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_diagonalize_symmetric(m, s, u) for the
 * other parameters.
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd)
{
    fs_diagonalize_symmetric_errbd(m, s, &u, &s_errbd);
}

/**
 * Same as fs_diagonalize_symmetric(m, s, u, s_errbd) except that
 * approximate error bounds for the singular vectors are returned as
 * well.  The error bounds are estimated following the method
 * presented at http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] u_errbd array of approximate error bounds for u
 *
 * See the documentation of fs_diagonalize_symmetric(m, s, u, s_errbd)
 * for the other parameters.
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Eigen::Matrix<std::complex<Real>, N, N>& u,
 Real& s_errbd,
 Eigen::Array<Real, N, 1>& u_errbd)
{
    fs_diagonalize_symmetric_errbd(m, s, &u, &s_errbd, &u_errbd);
}

/**
 * Returns singular values of N-by-N symmetric matrix m via s such
 * that `(s >= 0).all()`.  Elements of s are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m
 * @tparam     N      number of rows and columns in m and u
 * @param[in]  m      N-by-N symmetric matrix to be decomposed
 * @param[out] s      array of length N to contain singular values
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s)
{
    fs_diagonalize_symmetric_errbd(m, s);
}

/**
 * Same as fs_diagonalize_symmetric(m, s) except that an approximate
 * error bound for the singular values is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node96.html.
 *
 * @param[out] s_errbd approximate error bound for the elements of s
 *
 * See the documentation of fs_diagonalize_symmetric(m, s) for the
 * other parameters.
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_symmetric
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& s,
 Real& s_errbd)
{
    fs_diagonalize_symmetric_errbd(m, s, 0, &s_errbd);
}

template<class Real, class Scalar, int N>
void fs_diagonalize_hermitian_errbd
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N> *z = 0,
 Real *w_errbd = 0,
 Eigen::Array<Real, N, 1> *z_errbd = 0)
{
    diagonalize_hermitian_errbd(m, w, z, w_errbd, z_errbd);
    Eigen::PermutationMatrix<N> p;
    p.setIdentity();
    std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
              [&w] (int i, int j) { return std::abs(w[i]) < std::abs(w[j]); });
#if EIGEN_VERSION_AT_LEAST(3,1,4)
    w.matrix().transpose() *= p;
    if (z_errbd) { z_errbd->matrix().transpose() *= p; }
#else
    Eigen::Map<Eigen::Matrix<Real, N, 1> >(w.data()).transpose() *= p;
    if (z_errbd) {
        Eigen::Map<Eigen::Matrix<Real, N, 1>>(z_errbd->data()).transpose() *= p;
    }
#endif
    if (z) { *z = (*z * p).adjoint().eval(); }
}

/**
 * Diagonalizes N-by-N hermitian matrix m so that
 *
 *     m == z.adjoint() * w.matrix().asDiagonal() * z    // convention of SARAH
 *
 * w is arranged so that `abs(w[i])` are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m and z
 * @tparam     N      number of rows and columns in m and z
 * @param[in]  m      N-by-N matrix to be diagonalized
 * @param[out] w      array of length N to contain eigenvalues
 * @param[out] z      N-by-N unitary matrix
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z)
{
    fs_diagonalize_hermitian_errbd(m, w, &z);
}

/**
 * Same as fs_diagonalize_hermitian(m, w, z) except that an
 * approximate error bound for the eigenvalues is returned as well.
 * The error bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] w_errbd approximate error bound for the elements of w
 *
 * See the documentation of fs_diagonalize_hermitian(m, w, z) for the
 * other parameters.
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 Real& w_errbd)
{
    fs_diagonalize_hermitian_errbd(m, w, &z, &w_errbd);
}

/**
 * Same as fs_diagonalize_hermitian(m, w, z, w_errbd) except that
 * approximate error bounds for the eigenvectors are returned as well.
 * The error bounds are estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] z_errbd array of approximate error bounds for z
 *
 * See the documentation of fs_diagonalize_hermitian(m, w, z, w_errbd)
 * for the other parameters.
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Eigen::Matrix<Scalar, N, N>& z,
 Real& w_errbd,
 Eigen::Array<Real, N, 1>& z_errbd)
{
    fs_diagonalize_hermitian_errbd(m, w, &z, &w_errbd, &z_errbd);
}

/**
 * Returns eigenvalues of N-by-N hermitian matrix m via w.
 * w is arranged so that `abs(w[i])` are in ascending order.
 *
 * @tparam     Real   type of real and imaginary parts of Scalar
 * @tparam     Scalar type of elements of m and z
 * @tparam     N      number of rows and columns in m and z
 * @param[in]  m      N-by-N matrix to be diagonalized
 * @param[out] w      array of length N to contain eigenvalues
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w)
{
    fs_diagonalize_hermitian_errbd(m, w);
}

/**
 * Same as fs_diagonalize_hermitian(m, w) except that an approximate
 * error bound for the eigenvalues is returned as well.  The error
 * bound is estimated following the method presented at
 * http://www.netlib.org/lapack/lug/node89.html.
 *
 * @param[out] w_errbd approximate error bound for the elements of w
 *
 * See the documentation of fs_diagonalize_hermitian(m, w) for the
 * other parameters.
 */
template<class Real, class Scalar, int N>
void fs_diagonalize_hermitian
(const Eigen::Matrix<Scalar, N, N>& m,
 Eigen::Array<Real, N, 1>& w,
 Real& w_errbd)
{
    fs_diagonalize_hermitian_errbd(m, w, 0, &w_errbd);
}

/**
 * Calculates the mass of a singlet from a (possibly complex)
 * numerical value by taking the magnitude of the value.
 *
 * @param value numerical value
 * @return mass
 */
template <typename T>
double calculate_singlet_mass(T value) noexcept
{
   return std::abs(value);
}

/**
 * Calculates the mass of a Majoran fermion singlet from a (possibly
 * complex) numerical value by taking the magnitude of the value.
 *
 * The phase is set to exp(i theta/2), where theta is the phase angle
 * of the complex value.  If the value is pure real, then the phase
 * will be set to 1.  If the value is purely imaginary, then the phase
 * will be set to \f$e^{i \pi/2}\f$.
 *
 * @param value numerical value
 * @param[out] phase phase
 * @return mass
 */
template <typename T>
double calculate_majorana_singlet_mass(T value, std::complex<double>& phase)
{
   phase = std::polar(1., 0.5 * std::arg(std::complex<double>(value)));
   return std::abs(value);
}

/**
 * Calculates the mass of a Dirac fermion singlet from a (possibly
 * complex) numerical value by taking the magnitude of the value.
 *
 * The phase is set to exp(i theta), where theta is the phase angle of
 * the complex value.  If the value is pure real, then the phase will
 * be set to 1.  If the value is purely imaginary, then the phase will
 * be set to \f$e^{i \pi}\f$.
 *
 * @param value numerical value
 * @param[out] phase phase
 * @return mass
 */
template <typename T>
double calculate_dirac_singlet_mass(T value, std::complex<double>& phase)
{
   phase = std::polar(1., std::arg(std::complex<double>(value)));
   return std::abs(value);
}

} // namespace flexiblesusy

#endif // LINALG2_H
