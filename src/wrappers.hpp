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

#ifndef WRAPPERS_H
#define WRAPPERS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <Eigen/Core>

#include "eigen_tensor.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "if.hpp"
#include "sum.hpp"
#include "which.hpp"

namespace flexiblesusy {

static constexpr double Pi = M_PI;
static constexpr double oneOver16Pi = 1. / (16. * Pi);
static constexpr double oneOver16PiSqr = 1./(16. * Pi * Pi);
static constexpr double oneLoop = oneOver16PiSqr;
static constexpr double twoLoop = oneOver16PiSqr * oneOver16PiSqr;
static constexpr double threeLoop = oneOver16PiSqr * oneOver16PiSqr * oneOver16PiSqr;
static constexpr double fourLoop = twoLoop * twoLoop;
static constexpr double fiveLoop = oneLoop * fourLoop;
static constexpr bool True = true;

template <typename T>
T Abs(T a) noexcept
{
   return std::abs(a);
}

template <typename T>
T Abs(const std::complex<T>& z) noexcept
{
   return std::abs(z);
}

template <typename Scalar, int M, int N>
Eigen::Array<Scalar, M, N> Abs(const Eigen::Array<Scalar, M, N>& a)
{
   return a.cwiseAbs();
}

template <typename Scalar, int M, int N>
Eigen::Matrix<Scalar, M, N> Abs(const Eigen::Matrix<Scalar, M, N>& a)
{
   return a.cwiseAbs();
}

template <class T>
std::vector<T> Abs(std::vector<T> v) noexcept
{
   for (auto& e: v)
      e = Abs(e);
   return v;
}

double AbsSqr(double) noexcept;
double AbsSqr(const std::complex<double>&) noexcept;
double AbsSqrt(double) noexcept;

template <typename Derived>
Derived AbsSqrt(const Eigen::MatrixBase<Derived>& m)
{
   return m.cwiseAbs().cwiseSqrt();
}

template <typename Derived>
Derived AbsSqrt(const Eigen::ArrayBase<Derived>& m)
{
   return m.cwiseAbs().cwiseSqrt();
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

double ArcTan(double) noexcept;
double ArcSin(double) noexcept;
double ArcCos(double) noexcept;
double Arg(const std::complex<double>&) noexcept;

template <typename T>
constexpr T Cbrt(T a) noexcept
{
   return std::cbrt(a);
}

double Conj(double a) noexcept;
std::complex<double> Conj(const std::complex<double>& a) noexcept;

template<typename Scalar, int M, int N>
Eigen::Matrix<Scalar,M,N> Conj(const Eigen::Matrix<Scalar,M,N>& a) noexcept
{
   return a.conjugate();
}

template <class T>
T Conjugate(T a) noexcept
{
   return Conj(a);
}

template <typename T>
constexpr T Cube(T a) noexcept
{
   return a * a * a;
}

template <typename T>
T Exp(T z) noexcept
{
   return std::exp(z);
}

double Tan(double a) noexcept;
double Cot(double a) noexcept;
double Cos(double x) noexcept;
double Sin(double x) noexcept;
double Sec(double x) noexcept;
double Csc(double x) noexcept;
int Delta(int i, int j) noexcept;

#define FSFlagProblem(p) [&](){ (p); return 0.; }()
#define FSFlagWarning(p) [&](){ (p); return 0.; }()

bool IsClose(double, double, double eps = std::numeric_limits<double>::epsilon()) noexcept;
bool IsCloseRel(double, double, double eps = std::numeric_limits<double>::epsilon()) noexcept;
bool IsFinite(double) noexcept;
bool IsFinite(const std::complex<double>&) noexcept;

template <class Derived>
bool IsFinite(const Eigen::DenseBase<Derived>& m)
{
   // workaround for intel compiler / Eigen bug that causes unexpected
   // behavior from allFinite()
   const auto nr = m.rows();
   const auto nc = m.cols();

   for (int r = 0; r < nr; r++) {
      for (int c = 0; c < nc; c++) {
         if (!std::isfinite(m(r,c)))
            return false;
      }
   }

   return true;
}

int KroneckerDelta(int, int) noexcept;

template <class Derived>
typename Eigen::MatrixBase<Derived>::PlainObject Diag(const Eigen::MatrixBase<Derived>& m) noexcept
{
   static_assert(Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                 "Diag is only defined for squared matrices");

   typename Eigen::MatrixBase<Derived>::PlainObject diag(m);

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; ++i)
      for (int k = i + 1; k < Eigen::MatrixBase<Derived>::ColsAtCompileTime; ++k)
         diag(i,k) = 0.0;

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; ++i)
      for (int k = 0; k < i; ++k)
         diag(i,k) = 0.0;

   return diag;
}

std::complex<double> ComplexLog(double a) noexcept;
std::complex<double> ComplexLog(const std::complex<double>& z) noexcept;
double FiniteLog(double a) noexcept;

/**
 * Fills lower triangle of hermitian matrix from values
 * in upper triangle.
 *
 * @param m matrix
 */
template <typename Derived>
void Hermitianize(Eigen::PlainObjectBase<Derived>& m) noexcept
{
   static_assert(Eigen::PlainObjectBase<Derived>::RowsAtCompileTime ==
                 Eigen::PlainObjectBase<Derived>::ColsAtCompileTime,
                 "Hermitianize is only defined for squared matrices");

   for (int i = 0; i < Eigen::PlainObjectBase<Derived>::RowsAtCompileTime; i++) {
      for (int k = 0; k < i; k++) {
         m(i,k) = Conj(m(k,i));
      }
   }
}

///////////////////////// logger commands /////////////////////////

namespace {

template <typename Printer>
void PrintTo(std::stringstream& ostr, Printer&& printer)
{
   printer(ostr.str());
}

template <typename Printer, typename T0, typename... Ts>
void PrintTo(std::stringstream& ostr, Printer&& printer, T0&& v, Ts&&... vs)
{
   ostr << v;
   PrintTo(ostr, printer, std::forward<Ts>(vs)...);
}

} // anonymous namespace

/// print debug information to cerr
template<typename... Ts>
double PrintDEBUG(Ts&&... vs)
{
   std::stringstream ss;
   PrintTo(ss, [] (const std::string& s) { DEBUG_MSG(s); }, std::forward<Ts>(vs)...);
   return 0.;
}

/// print error to cerr
template<typename... Ts>
double PrintERROR(Ts&&... vs)
{
   std::stringstream ss;
   PrintTo(ss, [] (const std::string& s) { ERROR(s); }, std::forward<Ts>(vs)...);
   return 0.;
}

/// print error to cerr and stop program
template<typename... Ts>
double PrintFATAL(Ts&&... vs)
{
   std::stringstream ss;
   PrintTo(ss, [] (const std::string& s) { FATAL(s); }, std::forward<Ts>(vs)...);
   return 0.;
}

/// print an information message
template<typename... Ts>
double PrintINFO(Ts&&... vs)
{
   std::stringstream ss;
   PrintTo(ss, [] (const std::string& s) { INFO(s); }, std::forward<Ts>(vs)...);
   return 0.;
}

/// print verbose information to cerr
template<typename... Ts>
double PrintVERBOSE(Ts&&... vs)
{
   std::stringstream ss;
   PrintTo(ss, [] (const std::string& s) { VERBOSE_MSG(s); }, std::forward<Ts>(vs)...);
   return 0.;
}

/// print warning to cerr
template<typename... Ts>
double PrintWARNING(Ts&&... vs)
{
   std::stringstream ss;
   PrintTo(ss, [] (const std::string& s) { WARNING(s); }, std::forward<Ts>(vs)...);
   return 0.;
}

///////////////////////// end of logger commands /////////////////////////

double Log(double a) noexcept;

// MaxRelDiff //////////////////////////////////////////////////////////

double MaxRelDiff(double, double);

double MaxRelDiff(const std::complex<double>&, const std::complex<double>&);

template <class Derived>
auto MaxRelDiff(const Eigen::PlainObjectBase<Derived>& a,
                const Eigen::PlainObjectBase<Derived>& b)
   -> decltype(MaxRelDiff(a.data()[0], b.data()[0]))
{
   if (a.rows() != b.rows() || a.cols() != b.cols()) {
      throw SetupError("MaxRelDiff: Matrices/Vectors have different size!");
   }

   using Scalar_t = decltype(MaxRelDiff(a.data()[0], b.data()[0]));

   Scalar_t max = 0;

   for (Eigen::Index i = 0; i < a.size(); i++) {
      max = std::max(max, MaxRelDiff(a.data()[i], b.data()[i]));
   }

   return max;
}

// MaxAbsValue /////////////////////////////////////////////////////////

double MaxAbsValue(double x) noexcept;

double MaxAbsValue(const std::complex<double>& x) noexcept;

template <class Derived>
auto MaxAbsValue(const Eigen::MatrixBase<Derived>& x) -> decltype(x.cwiseAbs().maxCoeff())
{
   return x.cwiseAbs().maxCoeff();
}

template <class Derived>
auto MaxAbsValue(const Eigen::ArrayBase<Derived>& x) -> decltype(x.cwiseAbs().maxCoeff())
{
   return x.cwiseAbs().maxCoeff();
}

// Max /////////////////////////////////////////////////////////////////

template<typename T>
T Max(T&&t)
{
   return std::forward<T>(t);
}

template<typename T0, typename T1, typename... Ts>
typename std::common_type<T0, T1, Ts...>::type Max(T0&& val1, T1&& val2, Ts&&... vs)
{
   if (val2 < val1)
      return Max(val1, std::forward<Ts>(vs)...);
   else
      return Max(val2, std::forward<Ts>(vs)...);
}

template<typename T>
T Min(T&&t)
{
   return std::forward<T>(t);
}

template<typename T0, typename T1, typename... Ts>
typename std::common_type<T0, T1, Ts...>::type Min(T0&& val1, T1&& val2, Ts&&... vs)
{
   if (val2 < val1)
      return Min(val2, std::forward<Ts>(vs)...);
   else
      return Min(val1, std::forward<Ts>(vs)...);
}

int Sign(double x) noexcept;
int Sign(int x) noexcept;

template <typename T>
constexpr T Quad(T a) noexcept
{
   return a * a * a * a;
}

/// real polylogarithm
double PolyLog(int, double);

/// complex polylogarithm
std::complex<double> PolyLog(int, const std::complex<double>&);

template <typename Base, typename Exponent>
Base Power(Base base, Exponent exp) noexcept
{
   return std::pow(base, exp);
}

template <typename Base>
constexpr Base Power2(Base b) noexcept
{
   return b * b;
}

template <typename Base>
constexpr Base Power3(Base b) noexcept
{
   return b * b * b;
}

template <typename Base>
constexpr Base Power4(Base b) noexcept
{
   return Power2(Power2(b));
}

template <typename Base>
constexpr Base Power5(Base b) noexcept
{
   return Power4(b) * b;
}

template <typename Base>
constexpr Base Power6(Base b) noexcept
{
   return Power4(b) * Power2(b);
}

template <typename Base>
constexpr Base Power7(Base b) noexcept
{
   return Power6(b) * b;
}

template <typename Base>
constexpr Base Power8(Base b) noexcept
{
   return Power2(Power4(b));
}

template <typename Base>
constexpr Base Power9(Base b) noexcept
{
   return Power8(b) * b;
}

template <typename Base>
constexpr Base Power10(Base b) noexcept
{
   return Power8(b) * Power2(b);
}

template <typename Base>
constexpr Base Power11(Base b) noexcept
{
   return Power10(b) * b;
}

template <typename Base>
constexpr Base Power12(Base b) noexcept
{
   return Power2(Power6(b));
}

double Re(double) noexcept;
double Re(const std::complex<double>&) noexcept;

template<int M, int N>
Eigen::Matrix<double,M,N> Re(const Eigen::Matrix<double,M,N>& x)
{
   return x;
}

template<class Derived>
typename Eigen::Matrix<
   double,
   Eigen::MatrixBase<Derived>::RowsAtCompileTime,
   Eigen::MatrixBase<Derived>::ColsAtCompileTime>
Re(const Eigen::MatrixBase<Derived>& x)
{
   return x.real();
}

double Im(double) noexcept;
double Im(const std::complex<double>&) noexcept;

template<int M, int N>
Eigen::Matrix<double,M,N> Im(const Eigen::Matrix<double,M,N>&)
{
   return Eigen::Matrix<double,M,N>::Zero();
}

template<class Derived>
typename Eigen::Matrix<
   double,
   Eigen::MatrixBase<Derived>::RowsAtCompileTime,
   Eigen::MatrixBase<Derived>::ColsAtCompileTime>
Im(const Eigen::MatrixBase<Derived>& x)
{
   return x.imag();
}

template <typename T>
T RelDiff(T a, T b, T eps = std::numeric_limits<T>::epsilon()) noexcept
{
   const T max = std::max(a, b);

   if (std::abs(max) < eps)
      return T();

   return (a - b) / max;
}

int Round(double a) noexcept;

double SignedAbsSqrt(double a) noexcept;

template <typename Derived>
Derived SignedAbsSqrt(const Eigen::ArrayBase<Derived>& m)
{
   return m.unaryExpr([](double a) { return SignedAbsSqrt(a); });
}

template <class T, typename = typename std::enable_if<std::is_floating_point<T>::value,T>::type>
T Sqrt(T a) noexcept
{
   return std::sqrt(a);
}

template <class T, typename = typename std::enable_if<std::is_integral<T>::value,T>::type>
double Sqrt(T a) noexcept
{
   return std::sqrt(static_cast<double>(a));
}

template <typename Scalar, int M, int N>
Eigen::Array<Scalar, M, N> Sqrt(const Eigen::Array<Scalar, M, N>& m)
{
   return m.unaryExpr([](Scalar a){ return Sqrt(a); });
}

template <class T>
std::vector<T> Sqrt(std::vector<T> v)
{
   for (auto& e: v)
      e = Sqrt(e);
   return v;
}

template <typename T>
constexpr T Sqr(T a) noexcept
{
   return a * a;
}

template <typename Scalar, int M, int N>
Eigen::Array<Scalar, M, N> Sqr(const Eigen::Array<Scalar, M, N>& a)
{
   return a.unaryExpr([](Scalar a){ return Sqr(a); });
}

template <class T>
std::vector<T> Sqr(std::vector<T> v)
{
   for (auto& e: v)
      e = Sqr(e);
   return v;
}

#define DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(op)                     \
   template <typename T>                                                \
   std::complex<T> operator op(const std::complex<T>& lhs, int rhs)     \
   {                                                                    \
      return lhs op static_cast<T>(rhs);                                \
   }                                                                    \
                                                                        \
   template <typename T>                                                \
   std::complex<T> operator op(int lhs, const std::complex<T>& rhs)     \
   {                                                                    \
      return static_cast<T>(lhs) op rhs;                                \
   }

DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(*)
DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(/)
DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(+)
DEFINE_COMMUTATIVE_OPERATOR_COMPLEX_INT(-)

/**
 * Fills lower triangle of symmetric matrix from values in upper
 * triangle.
 *
 * @param m matrix
 */
template <typename Derived>
void Symmetrize(Eigen::PlainObjectBase<Derived>& m)
{
   static_assert(Eigen::PlainObjectBase<Derived>::RowsAtCompileTime ==
                 Eigen::PlainObjectBase<Derived>::ColsAtCompileTime,
                 "Symmetrize is only defined for squared matrices");

   for (int i = 0; i < Eigen::PlainObjectBase<Derived>::RowsAtCompileTime; i++) {
      for (int k = 0; k < i; k++) {
         m(i,k) = m(k,i);
      }
   }
}

#define UNITMATRIX(rows)             Eigen::Matrix<double,rows,rows>::Identity()
#define ZEROMATRIX(rows,cols)        Eigen::Matrix<double,rows,cols>::Zero()
#define ZEROTENSOR3(d1,d2,d3)        ZeroTensor3<double,d1,d2,d3>()
#define ZEROTENSOR4(d1,d2,d3,d4)     ZeroTensor4<double,d1,d2,d3,d4>()
#define ZEROVECTOR(rows)             Eigen::Matrix<double,rows,1>::Zero()
#define ZEROARRAY(rows)              Eigen::Array<double,rows,1>::Zero()
#define UNITMATRIXCOMPLEX(rows)      Eigen::Matrix<std::complex<double>,rows,rows>::Identity()
#define ZEROMATRIXCOMPLEX(rows,cols) Eigen::Matrix<std::complex<double>,rows,cols>::Zero()
#define ZEROVECTORCOMPLEX(rows)      Eigen::Matrix<std::complex<double>,rows,1>::Zero()
#define ZEROTENSOR3COMPLEX(d1,d2,d3) ZeroTensor3<std::complex<double>,d1,d2,d3>()
#define ZEROTENSOR4COMPLEX(d1,d2,d3,d4) ZeroTensor4<std::complex<double>,d1,d2,d3,d4>()
#define ZEROARRAYCOMPLEX(rows)       Eigen::Array<std::complex<double>,rows,1>::Zero()

// MxN matrix projection operator, which projects on the (X,Y)
// component
#define PROJECTOR Proj
#define DEFINE_PROJECTOR(M,N,X,Y)                                       \
   Eigen::Matrix<double,M,N> Proj(Eigen::Matrix<double,M,N>::Zero());   \
   Proj((X)-1,(Y)-1) = 1;

inline double FSThrow(const std::string& s)
{
   throw PhysicalError(s);
   return 0.;
}

template<class Scalar, int M>
Eigen::Matrix<Scalar,M,M> ToMatrix(const Eigen::Array<Scalar,M,1>& a)
{
   return Eigen::Matrix<Scalar,M,M>(a.matrix().asDiagonal());
}

template<class Scalar, int M, int N>
Eigen::Matrix<Scalar,M,N> ToMatrix(const Eigen::Matrix<Scalar,M,N>& a) noexcept
{
   return a;
}

// ToString ////////////////////////////////////////////////////////////

std::string ToString(char);

std::string ToString(unsigned char);
std::string ToString(unsigned short);
std::string ToString(unsigned int);
std::string ToString(unsigned long);
std::string ToString(unsigned long long);

std::string ToString(signed char);
std::string ToString(signed short);
std::string ToString(signed int);
std::string ToString(signed long);
std::string ToString(signed long long);

std::string ToString(double);
std::string ToString(const std::complex<double>&);

// Total ///////////////////////////////////////////////////////////////

double Total(double) noexcept;
std::complex<double> Total(const std::complex<double>&) noexcept;

template <class T>
T Total(const std::vector<T>& v)
{
   return std::accumulate(v.begin(), v.end(), T(0));
}

template <typename Scalar, int M, int N>
Scalar Total(const Eigen::Array<Scalar, M, N>& a)
{
   return a.sum();
}

template <typename Scalar, int M, int N>
Scalar Total(const Eigen::Matrix<Scalar, M, N>& a)
{
   return a.sum();
}

template <class Scalar, int M, int N>
Eigen::Array<Scalar,M,N> Total(const std::vector<Eigen::Array<Scalar,M,N> >& v)
{
   if (v.empty()) {
      Eigen::Array<Scalar,M,N> result(0,0);
      result.setZero();
      return result;
   }

   Eigen::Array<Scalar,M,N> result(v[0].rows(), v[0].cols());
   result.setZero();

   for (std::size_t i = 0; i < v.size(); i++)
      result += v[i];

   return result;
}

/// unit vector of length N into direction i
template <int N, int i, typename Scalar = double>
constexpr auto UnitVector() -> Eigen::Matrix<Scalar,N,1>
{
   return Eigen::Matrix<Scalar,N,1>::Unit(i);
}

/// unit vector of length N into direction i
template <int N, typename Scalar = double>
constexpr auto UnitVector(int i) -> Eigen::Matrix<Scalar,N,1>
{
   return Eigen::Matrix<Scalar,N,1>::Unit(i);
}

/// unit vector of length N into direction i
Eigen::VectorXd UnitVector(int N, int i);

/// matrix projector of size MxN into direction i, j
template <int M, int N, int i, int j, typename Scalar = double>
auto MatrixProjector() -> Eigen::Matrix<Scalar,M,N>
{
   Eigen::Matrix<Scalar,M,N> proj(Eigen::Matrix<Scalar,M,N>::Zero());
   proj(i,j) = 1;

   return proj;
}

/// matrix projector of size MxN into direction i, j
template <int M, int N, typename Scalar = double>
auto MatrixProjector(int i, int j) -> Eigen::Matrix<Scalar,M,N>
{
   Eigen::Matrix<Scalar,M,N> proj(Eigen::Matrix<Scalar,M,N>::Zero());
   proj(i,j) = 1;

   return proj;
}

/// unit matrix projector of size MxN into direction i, j
Eigen::MatrixXd MatrixProjector(int M, int N, int i, int j);

/// step function (0 for x < 0, 1 otherwise)
template <typename T>
constexpr int UnitStep(T x) noexcept
{
   return x < T() ? 0 : 1;
}

double ZeroSqrt(double x) noexcept;

template <typename Derived>
Derived ZeroSqrt(const Eigen::ArrayBase<Derived>& m)
{
   return m.unaryExpr([](double a){ return ZeroSqrt(a); });
}

} // namespace flexiblesusy

#endif
