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

#include "wrappers.hpp"
#include "dilog.hpp"
#include "numerics2.hpp"

namespace flexiblesusy {

double AbsSqr(double z) noexcept
{
   return z * z;
}

double AbsSqr(const std::complex<double>& z) noexcept
{
   return std::norm(z);
}

double AbsSqrt(double x) noexcept
{
   return std::sqrt(std::fabs(x));
}

double ArcTan(double a) noexcept
{
   return std::atan(a);
}

double ArcSin(double a) noexcept
{
   return std::asin(a);
}

double ArcCos(double a) noexcept
{
   return std::acos(a);
}

double Arg(const std::complex<double>& z) noexcept
{
   return std::arg(z);
}

double Conj(double a) noexcept
{
   return a;
}

std::complex<double> Conj(const std::complex<double>& a) noexcept
{
   return std::conj(a);
}

double Tan(double a) noexcept
{
   return std::tan(a);
}

double Cot(double a) noexcept
{
   return 1./Tan(a);
}

double Cos(double x) noexcept
{
   return std::cos(x);
}

double Sin(double x) noexcept
{
   return std::sin(x);
}

double Sec(double x) noexcept
{
   return 1./Cos(x);
}

double Csc(double x) noexcept
{
   return 1./Sin(x);
}

int Delta(int i, int j) noexcept
{
   return i == j;
}

bool IsClose(double a, double b, double eps) noexcept
{
   return std::abs(a - b) < eps;
}

bool IsCloseRel(double a, double b, double eps) noexcept
{
   if (IsClose(a, b, std::numeric_limits<double>::epsilon()))
      return true;

   if (std::abs(a) < std::numeric_limits<double>::epsilon())
      return IsClose(a, b, eps);

   return std::abs((a - b)/a) < eps;
}

bool IsFinite(double x) noexcept
{
   return std::isfinite(x);
}

bool IsFinite(const std::complex<double>& x) noexcept
{
   return std::isfinite(x.real()) && std::isfinite(x.imag());
}

int KroneckerDelta(int i, int j) noexcept
{
   return i == j;
}

double Log(double a) noexcept
{
   return std::log(a);
}

std::complex<double> ComplexLog(double a) noexcept
{
   return fast_log(std::complex<double>(a,0.));
}

std::complex<double> ComplexLog(const std::complex<double>& z) noexcept
{
   return fast_log(z);
}

double FiniteLog(double a) noexcept
{
   const double l = std::log(a);
   return std::isfinite(l) ? l : 0.;
}

double MaxAbsValue(double x) noexcept
{
   return Abs(x);
}

double MaxAbsValue(const std::complex<double>& x) noexcept
{
   return Abs(x);
}

double MaxRelDiff(double a, double b)
{
   const double sTin = fabs(a);
   const double sTout = fabs(b);
   const double maxx = std::max(sTin, sTout);
   const double underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return std::abs((a - b) / maxx);
}

double MaxRelDiff(const std::complex<double>& a, const std::complex<double>& b)
{
   const double sTin = std::abs(a);
   const double sTout = std::abs(b);
   const double maxx = std::max(sTin, sTout);
   const double underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return std::abs(a - b) / maxx;
}

double PolyLog(int n, double z)
{
   if (n == 2)
      return dilog(z);
   throw SetupError("PolyLog(n!=2) not implemented");
}

std::complex<double> PolyLog(int n, const std::complex<double>& z)
{
   if (n == 2)
      return dilog(z);
   throw SetupError("PolyLog(n!=2) not implemented");
}

double Re(double x) noexcept
{
   return x;
}

double Re(const std::complex<double>& x) noexcept
{
   return std::real(x);
}

int Round(double a) noexcept
{
   return static_cast<int>(a >= 0. ? a + 0.5 : a - 0.5);
}

double Im(double) noexcept
{
   return 0.;
}

double Im(const std::complex<double>& x) noexcept
{
   return std::imag(x);
}

int Sign(double x) noexcept
{
   return (x >= 0.0 ? 1 : -1);
}

int Sign(int x) noexcept
{
   return (x >= 0 ? 1 : -1);
}

double SignedAbsSqrt(double a) noexcept
{
   return Sign(a) * AbsSqrt(a);
}

double Total(double a) noexcept
{
   return a;
}

std::complex<double> Total(const std::complex<double>& a) noexcept
{
   return a;
}

/// unit vector of length N into direction i
Eigen::VectorXd UnitVector(int N, int i)
{
   Eigen::VectorXd v = Eigen::VectorXd::Zero(N);
   v(i) = 1;

   return v;
}

/// unit matrix projector of size MxN into direction i, j
Eigen::MatrixXd MatrixProjector(int M, int N, int i, int j)
{
   Eigen::MatrixXd m = Eigen::MatrixXd::Zero(M,N);
   m(i,j) = 1;

   return m;
}

double ZeroSqrt(double x) noexcept
{
   return (x > 0.0 ? std::sqrt(x) : 0.0);
}

} // namespace flexiblesusy
