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

#include "diagonalization.hpp"
#include "error.hpp"
#include "wrappers.hpp"
#include <sstream>

#ifndef ENABLE_DEBUG
#define ENABLE_DEBUG
#endif

namespace flexiblesusy {

namespace {
   double AbsSqrt_d(double x)
   {
      return AbsSqrt(x);
   }
}

softsusy::DoubleVector AbsSqrt(const softsusy::DoubleVector& x)
{
   return x.apply(AbsSqrt_d);
}

/**
 * diag = u m u^T
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const softsusy::DoubleMatrix& m, softsusy::DoubleMatrix& u,
                 softsusy::DoubleVector& eigenvalues)
{
   const int c = m.displayCols();
#ifdef ENABLE_DEBUG
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       u.displayCols() != c || u.displayRows() !=c || c > 10) {
      std::ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m;
      throw DiagonalizationError(ii.str());
   }
#endif
   // Numerical recipes routine replaces argument, so make a copy of elements
   softsusy::DoubleMatrix a(m);
   int nrot;
   diagonaliseJac(a, c, eigenvalues, u, &nrot);
   u.associateOrderAbs(eigenvalues);
   u = u.transpose();
}

/**
 * Returns a phase rotation matrix which will transform a mixing
 * matrix such that the corresponding eigenvalues are positive.
 *
 * @param v eigen values
 *
 * @return phase rotation matrix
 */
softsusy::ComplexMatrix phase_rotation(const softsusy::DoubleVector& v)
{
   const int len = v.size();
   softsusy::ComplexMatrix rot(len, len);

   for (int i = 1; i <= len; i++) {
      if (v(i) < 0.0)
         rot(i, i) = softsusy::Complex(0.0, -1.0);
      else
         rot(i, i) = softsusy::Complex(1.0, 0.0);
   }

   return rot;
}

/**
 * diag = u^* m u^+
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const softsusy::DoubleMatrix& m, softsusy::ComplexMatrix& u, softsusy::DoubleVector& eigenvalues)
{
   const int c = m.displayCols();
#ifdef ENABLE_DEBUG
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       u.displayCols() != c || u.displayRows() !=c || c > 10) {
      std::ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m;
      throw DiagonalizationError(ii.str());
   }
#endif
   // Numerical recipes routine replaces argument, so make a copy of elements
   softsusy::DoubleMatrix a(m), k(c, c);
   int nrot;
   diagonaliseJac(a, c, eigenvalues, k, &nrot);

   // So Far, V is such that the eigenvalues are in some random order.
   // We must now re-order the rows of V to get the eigenvalues in the
   // correct order (taken to be increasing abs(eigenvalues)).
   k.associateOrderAbs(eigenvalues);

   // We want to change the PHASES of the neutralino mixing matrix in order to
   // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
   const softsusy::ComplexMatrix trans(phase_rotation(eigenvalues));
   u = trans * k.transpose();
   eigenvalues = eigenvalues.apply(fabs);
}

void Diagonalize2by2(const softsusy::DoubleMatrix& m, softsusy::DoubleMatrix& u,
                     softsusy::DoubleVector& eigenvalues)
{
   double theta;
   eigenvalues = m.sym2by2(theta);
   if (std::fabs(eigenvalues(1)) > std::fabs(eigenvalues(2))) {
      theta += 0.5 * M_PI;
      std::swap(eigenvalues(1), eigenvalues(2));
   }
   u = softsusy::rot2d(theta);
}

void Diagonalize2by2(const softsusy::DoubleMatrix& m, softsusy::ComplexMatrix& u,
                     softsusy::DoubleVector& eigenvalues)
{
   double theta;
   eigenvalues = m.sym2by2(theta);
   const softsusy::ComplexMatrix a(phase_rotation(eigenvalues).complexConjugate());
   u = a * softsusy::rot2d(theta);
   eigenvalues = eigenvalues.apply(fabs);
}

void AssociateOrderAbs(softsusy::DoubleMatrix& u, softsusy::DoubleMatrix& v, softsusy::DoubleVector& w)
{
#ifdef ENABLE_DEBUG
   if ((u.displayRows() != w.displayEnd()) ||
       (u.displayCols() != w.displayEnd()) ||
       (w.displayStart() != 1)) {
      std::ostringstream ii;
      ii << "Associated ordering incompatibility\n";
      throw DiagonalizationError(ii.str());
   }
#endif
   for (int i = w.displayStart(); i <= w.displayEnd(); ++i) {
      for (int j = i + 1; j <= w.displayEnd(); ++j) {
         if (std::fabs(w.display(i)) > std::fabs(w.display(j))) {
            w.swap(i, j);
            v.swapcols(i, j);
            u.swapcols(i, j);
         }
      }
   }
}

/**
 * diag = u m v^T
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param v mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const softsusy::DoubleMatrix& m, softsusy::DoubleMatrix& u,
                 softsusy::DoubleMatrix& v, softsusy::DoubleVector& eigenvalues)
{
#ifdef ENABLE_DEBUG
   const int c = m.displayCols();
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       v.displayCols() != c || v.displayRows() !=c ||
       u.displayCols() !=c || u.displayRows() !=c) {
      std::ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m
         << "with u" << u << "v " << v << "eigenvalues " << eigenvalues;
      throw DiagonalizationError(ii.str());
   }
#endif
   // Numerical routine replaces argument, so make a copy of elements
   u = m;
   diagonaliseSvd(u, eigenvalues, v);

   // Eigen::MatrixXd me(ToEigen(m)), ue(ToEigen(u)), ve(ToEigen(v));
   // Eigen::JacobiSVD<Eigen::MatrixXd> eigensolver(me, Eigen::ComputeThinU | 
   // u = ToSoftsusy(eigensolver.matrixU());
   // v = ToSoftsusy(eigensolver.matrixV());
   // eigenvalues = ToSoftsusy(eigensolver.singularValues());

   // Phase freedom - tends to give more familiar matrices
   u = -1.0 * u; v = -1.0 * v;

   AssociateOrderAbs(u, v, eigenvalues);

   v = v.transpose();
   u = u.transpose();
}

/**
 * diag = u^* m v^+
 *
 * @param m mass matrix
 * @param u mixing matrix
 * @param v mixing matrix
 * @param eigenvalues eigenvalues (sorted)
 */
void Diagonalize(const softsusy::DoubleMatrix& m, softsusy::ComplexMatrix& u,
                 softsusy::ComplexMatrix& v, softsusy::DoubleVector& eigenvalues)
{
#ifdef ENABLE_DEBUG
   const int c = m.displayCols();
   if (m.displayRows() != c || eigenvalues.displayEnd() != c ||
       v.displayCols() != c || v.displayRows() !=c || c > 10) {
      std::ostringstream ii;
      ii << "Error: Trying to diagonalise matrix \n" << m;
      throw DiagonalizationError(ii.str());
   }
#endif
   softsusy::DoubleMatrix realU(u.real()), realV(v.real());
   Diagonalize(m, realU, realV, eigenvalues);

   const softsusy::ComplexMatrix trans(phase_rotation(eigenvalues));
   u = trans * realU;
   v = trans * realV;
   eigenvalues = eigenvalues.apply(fabs);
}

void Diagonalize2by2(const softsusy::DoubleMatrix& m, softsusy::ComplexMatrix& u,
                     softsusy::ComplexMatrix& v, softsusy::DoubleVector& eigenvalues)
{
   double thetaL, thetaR;
   eigenvalues = m.asy2by2(thetaL, thetaR);
   positivise(thetaL, thetaR, eigenvalues, u, v);
   eigenvalues = eigenvalues.apply(fabs);
}

double MaxRelDiff(const softsusy::DoubleVector& a, const softsusy::DoubleVector& b)
{
   return a.compare(b);
}

void Symmetrize(softsusy::DoubleMatrix& m)
{
   const int r = m.displayRows();
   const int c = m.displayCols();
   for (int i = 1; i <= r; i++)
      for (int k = 1; k < i && k <= c; k++)
         m(i,k) = m(k,i);
}

namespace {
  inline double ZeroSqrt_d(double x) noexcept
  {
    return ZeroSqrt(x);
  }
} // anonymous namespace

softsusy::DoubleVector ZeroSqrt(const softsusy::DoubleVector& x)
{
   return x.apply(ZeroSqrt_d);
}

} // namespace flexiblesusy
