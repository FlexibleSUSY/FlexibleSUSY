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
#include "string_format.hpp"
#include "trilog.hpp"

#include <complex>
#include <cmath>

namespace flexiblesusy {

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
   return is_equal_rel(a, b, eps);
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

double MaxRelDiff(double a, double b) noexcept
{
   const double max = std::max(std::abs(a), std::abs(b));

   if (max < 1.0e-20) {
      return 0.0;
   }

   return std::abs((a - b) / max);
}

double MaxRelDiff(const std::complex<double>& a, const std::complex<double>& b) noexcept
{
   const double max = std::max(std::abs(a), std::abs(b));

   if (max < 1.0e-20) {
      return 0.0;
   }

   return std::abs((a - b) / max);
}

double PolyLog(int n, double z) noexcept
{
   return std::real(PolyLog(n, std::complex<double>(z, 0.0)));
}

std::complex<double> PolyLog(int n, const std::complex<double>& z) noexcept
{
   switch (n) {
   case 1: return -std::log(1.0 - z);
   case 2: return dilog(z);
   case 3: return trilog(z);
   default: break;
   }

   ERROR("PolyLog(n != 1|2|3) not implemented");

   return { 0.0, 0.0 };
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

#define DEFINE_ToString(type)                   \
   std::string ToString(type a)                 \
   {                                            \
      return flexiblesusy::to_string(a);        \
   }

DEFINE_ToString(char)
DEFINE_ToString(unsigned char)
DEFINE_ToString(unsigned short)
DEFINE_ToString(unsigned int)
DEFINE_ToString(unsigned long)
DEFINE_ToString(unsigned long long)
DEFINE_ToString(signed char)
DEFINE_ToString(signed short)
DEFINE_ToString(signed int)
DEFINE_ToString(signed long)
DEFINE_ToString(signed long long)
DEFINE_ToString(double)
DEFINE_ToString(const std::complex<double>&)

double Total(double a) noexcept
{
   return a;
}

std::complex<double> Total(const std::complex<double>& a) noexcept
{
   return a;
}

/// unit vector of length N into direction i
Eigen::VectorXd UnitVector(int N, int i) noexcept
{
   Eigen::VectorXd v = Eigen::VectorXd::Zero(N);
   v(i) = 1;

   return v;
}

/// unit matrix projector of size MxN into direction i, j
Eigen::MatrixXd MatrixProjector(int M, int N, int i, int j) noexcept
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
