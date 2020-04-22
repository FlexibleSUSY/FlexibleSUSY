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
#define BOOST_TEST_MODULE test_wrappers

#include <complex>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include "wrappers.hpp"
#include "stopwatch.hpp"
#include <boost/lexical_cast.hpp>

#if defined(__CYGWIN__) || defined(__FreeBSD__)
// std::to_string is missing on Cygwin and FreeBSD
// see http://stackoverflow.com/questions/22571838/gcc-4-8-1-stdto-string-error
#   undef  ENABLE_STD_TO_STRING
#else
#   define ENABLE_STD_TO_STRING 1
#endif

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_Constants )
{
   BOOST_CHECK_EQUAL(Pi            , 3.141592653589793238462643383279502884);
   BOOST_CHECK_EQUAL(oneOver16PiSqr, 0.006332573977646110715242466450607977432);
   BOOST_CHECK_EQUAL(oneLoop       , 0.006332573977646110715242466450607977432);
   BOOST_CHECK_EQUAL(twoLoop       , 0.00004010149318236068433262805963718239899);
   BOOST_CHECK_EQUAL(threeLoop     , 2.539456721913701894750580163364759355e-7);
   BOOST_CHECK_EQUAL(fourLoop      , 1.608129755454920445735606800552522082e-9);
   BOOST_CHECK_EQUAL(fiveLoop      , 1.018360064207223287777021556086771147e-11);
   BOOST_CHECK(True);
}

BOOST_AUTO_TEST_CASE( test_Abs )
{
   const std::complex<double> i(0.0, 1.0);

   BOOST_CHECK_EQUAL(Abs(-1), 1);
   BOOST_CHECK_EQUAL(Abs( 0), 0);
   BOOST_CHECK_EQUAL(Abs( 1), 1);

   BOOST_CHECK_EQUAL(Abs(-1.0), 1.0);
   BOOST_CHECK_EQUAL(Abs( 0.0), 0.0);
   BOOST_CHECK_EQUAL(Abs( 1.0), 1.0);

   BOOST_CHECK_EQUAL(Abs(i    ), 1.0);
   BOOST_CHECK_EQUAL(Abs(i*i  ), 1.0);
   BOOST_CHECK_EQUAL(Abs(1 + i), std::sqrt(2.0));

   Eigen::Matrix<double,2,2> m;
   m << -1.0, -2.0, -3.0, -4.0;

   BOOST_CHECK_EQUAL(Abs(m)(0,0), Abs(m(0,0)));
   BOOST_CHECK_EQUAL(Abs(m)(0,1), Abs(m(0,1)));
   BOOST_CHECK_EQUAL(Abs(m)(1,0), Abs(m(1,0)));
   BOOST_CHECK_EQUAL(Abs(m)(1,1), Abs(m(1,1)));

   Eigen::Array<double,2,2> a;
   a << -1.0, -2.0, -3.0, -4.0;

   BOOST_CHECK_EQUAL(Abs(a)(0,0), Abs(a(0,0)));
   BOOST_CHECK_EQUAL(Abs(a)(0,1), Abs(a(0,1)));
   BOOST_CHECK_EQUAL(Abs(a)(1,0), Abs(a(1,0)));
   BOOST_CHECK_EQUAL(Abs(a)(1,1), Abs(a(1,1)));
}

BOOST_AUTO_TEST_CASE( test_Delta )
{
   BOOST_CHECK_EQUAL(Delta(0,0), 1);
   BOOST_CHECK_EQUAL(Delta(1,1), 1);
   BOOST_CHECK_EQUAL(Delta(2,2), 1);
   BOOST_CHECK_EQUAL(Delta(3,3), 1);

   BOOST_CHECK_EQUAL(Delta(-1,-1), 1);
   BOOST_CHECK_EQUAL(Delta(-2,-2), 1);
   BOOST_CHECK_EQUAL(Delta(-3,-3), 1);

   BOOST_CHECK_EQUAL(Delta(0,1), 0);
   BOOST_CHECK_EQUAL(Delta(1,0), 0);
   BOOST_CHECK_EQUAL(Delta(0,2), 0);
   BOOST_CHECK_EQUAL(Delta(2,0), 0);

   BOOST_CHECK_EQUAL(Delta(0,-1), 0);
   BOOST_CHECK_EQUAL(Delta(-1,0), 0);
   BOOST_CHECK_EQUAL(Delta(0,-2), 0);
   BOOST_CHECK_EQUAL(Delta(-2,0), 0);
}

BOOST_AUTO_TEST_CASE( test_UnitStep )
{
   BOOST_CHECK_EQUAL(UnitStep(-1  ), 0);
   BOOST_CHECK_EQUAL(UnitStep(-2  ), 0);
   BOOST_CHECK_EQUAL(UnitStep(-0.5), 0);
   BOOST_CHECK_EQUAL(UnitStep(0   ), 1);
   BOOST_CHECK_EQUAL(UnitStep(0.5 ), 1);
   BOOST_CHECK_EQUAL(UnitStep(1   ), 1);
   BOOST_CHECK_EQUAL(UnitStep(2   ), 1);
}

BOOST_AUTO_TEST_CASE(test_Diag)
{
   Eigen::Matrix<double,3,3> m;
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         m(i,k) = (i+1) * (k+1);

   Eigen::Matrix<double,3,3> diag(Diag(m));

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++) {
         if (i == k)
            BOOST_CHECK_PREDICATE(std::not_equal_to<double>(), (diag(i,k))(0.));
         else
            BOOST_CHECK_EQUAL(diag(i,k), 0.);
      }
}

template <typename T>
std::string ToString_sstream(T a)
{
   std::ostringstream ostr;
   ostr << a;
   return ostr.str();
}

#ifdef ENABLE_STD_TO_STRING
template <typename T>
std::string ToString_to_string(T a)
{
   return std::to_string(a);
}
#endif

template <typename T>
std::string ToString_sprintf(T a)
{
   static const int buf_length = 20;
   char buf[buf_length];
   snprintf(buf, buf_length, "%i", a);
   return std::string(buf);
}

template <typename T>
std::string ToString_lexical_cast(T a)
{
   return boost::lexical_cast<std::string>(a);
}

#define MEASURE(type,number,iterations)                            \
   do {                                                            \
      Stopwatch stopwatch;                                         \
      double time = 0.;                                            \
      for (int i = 0; i < iterations; i++) {                       \
         stopwatch.start();                                        \
         ToString_##type(number);                                  \
         stopwatch.stop();                                         \
         time += stopwatch.get_time_in_seconds();                  \
      }                                                            \
      BOOST_TEST_MESSAGE("ToString via " #type ": " << time << " s");   \
   } while (0)

BOOST_AUTO_TEST_CASE(test_ToString)
{
   const int number = 123456;
   const int number_of_iterations = 1000000;

   MEASURE(sstream     , number, number_of_iterations);
#ifdef ENABLE_STD_TO_STRING
   MEASURE(to_string   , number, number_of_iterations);
#endif
   MEASURE(sprintf     , number, number_of_iterations);
   MEASURE(lexical_cast, number, number_of_iterations);
}

BOOST_AUTO_TEST_CASE(test_calculate_majorana_singlet_mass)
{
   double mass;
   std::complex<double> phase;

   mass = calculate_majorana_singlet_mass(100., phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), 1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_majorana_singlet_mass(-100., phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_SMALL(phase.real(), 1e-15);
   BOOST_CHECK_EQUAL(phase.imag(), 1.);

   mass = calculate_majorana_singlet_mass(std::complex<double>(100.,0.), phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), 1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_majorana_singlet_mass(std::complex<double>(-100.,0.), phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_SMALL(phase.real(), 1e-15);
   BOOST_CHECK_EQUAL(phase.imag(), 1.);

   mass = calculate_majorana_singlet_mass(std::complex<double>(1.,1.), phase);
   BOOST_CHECK_EQUAL(mass, sqrt(2.));
   BOOST_CHECK_EQUAL(std::abs(phase), 1.);
   BOOST_CHECK_CLOSE(std::arg(phase), 0.5 * Pi/4., 1e-13);
}

BOOST_AUTO_TEST_CASE(test_calculate_dirac_singlet_mass)
{
   double mass;
   std::complex<double> phase;

   mass = calculate_dirac_singlet_mass(100., phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), 1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_dirac_singlet_mass(-100., phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), -1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_dirac_singlet_mass(std::complex<double>(100.,0.), phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), 1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_dirac_singlet_mass(std::complex<double>(-100.,0.), phase);
   BOOST_CHECK_EQUAL(mass, 100.);
   BOOST_CHECK_EQUAL(phase.real(), -1.);
   BOOST_CHECK_SMALL(phase.imag(), 1e-15);

   mass = calculate_dirac_singlet_mass(std::complex<double>(1.,1.), phase);
   BOOST_CHECK_EQUAL(mass, sqrt(2.));
   BOOST_CHECK_EQUAL(std::abs(phase), 1.);
   BOOST_CHECK_CLOSE(std::arg(phase), Pi/4., 1e-13);
}

BOOST_AUTO_TEST_CASE(test_IF)
{
   BOOST_CHECK_EQUAL(IF(true , 1., 2.), 1.);
   BOOST_CHECK_EQUAL(IF(false, 1., 2.), 2.);
   BOOST_CHECK_EQUAL(IF(true , 1 , 2.), 1 );
   BOOST_CHECK_EQUAL(IF(false, 1., 2 ), 2 );

   BOOST_CHECK_EQUAL(std::real(IF(true , std::complex<double>(1.,1.), std::complex<double>(2.,2.))), 1.);
   BOOST_CHECK_EQUAL(std::real(IF(false, std::complex<double>(1.,1.), std::complex<double>(2.,2.))), 2.);

   BOOST_CHECK_EQUAL(std::real(IF(true , 1, std::complex<double>(2.,2.))), 1);
   BOOST_CHECK_EQUAL(std::real(IF(false, std::complex<double>(1.,1.), 2)), 2);
}

BOOST_AUTO_TEST_CASE(test_WHICH)
{
   // BOOST_CHECK_EQUAL(WHICH(true), 1.); // must not compile
   BOOST_CHECK_EQUAL(WHICH(true , 2.), 2.);
   BOOST_CHECK_EQUAL(WHICH(false, 2., true, 3.), 3.);
   BOOST_CHECK_EQUAL(WHICH(false, 2., false, 3., true, 4.), 4.);

   BOOST_CHECK_EQUAL(std::real(WHICH(true , std::complex<double>(2.,2.))), 2.);
   BOOST_CHECK_EQUAL(WHICH(false, std::complex<double>(2.,2.), true, 3.), 3.);
   BOOST_CHECK_EQUAL(WHICH(false, std::complex<double>(2.,2.), false, 3., true, 4.), 4.);

   BOOST_CHECK_EQUAL(WHICH(false, 0, true, 0.5), 0.5);
   BOOST_CHECK_EQUAL(WHICH(false, 0, false, 0, true, 0.5), 0.5);
}

BOOST_AUTO_TEST_CASE(test_MaxRelDiff)
{
   BOOST_CHECK_CLOSE(MaxRelDiff(0., 0.)  , 0. , 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(1., 0.)  , 1. , 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(-1., 0.) , 1. , 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(1., -1.) , 2. , 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(-1., -1.), 0. , 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(-1., -2.), 0.5, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(1., 2.), 0.5, 1e-10);

   Eigen::Matrix<double,1,1> M1, M2;
   M1(0) = 1.;
   M2(0) = -1.;

   BOOST_CHECK_CLOSE(MaxRelDiff(M1,M2), 2. , 1e-10);
}


BOOST_AUTO_TEST_CASE(test_MaxRelDiff_complex)
{
   const std::complex<double> unit(1.0, 0.0);

   BOOST_CHECK_CLOSE(MaxRelDiff( 0.0*unit,  0.0*unit), 0.0, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff( 1.0*unit,  0.0*unit), 1.0, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(-1.0*unit,  0.0*unit), 1.0, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff( 1.0*unit, -1.0*unit), 2.0, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(-1.0*unit, -1.0*unit), 0.0, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff(-1.0*unit, -2.0*unit), 0.5, 1e-10);
   BOOST_CHECK_CLOSE(MaxRelDiff( 1.0*unit,  2.0*unit), 0.5, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_MaxRelDiff_array)
{
   Eigen::Array<double,2,2> M1, M2;
   M1 << +1.0, +1.0, +1.0, +1.0;
   M2 << -1.0, -2.0, -3.0, -4.0;

   BOOST_CHECK_CLOSE(MaxRelDiff(M1,M2), 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_MaxRelDiff_matrix)
{
   Eigen::Matrix<double,2,2> M1, M2;
   M1 << +1.0, +1.0, +1.0, +1.0;
   M2 << -1.0, -2.0, -3.0, -4.0;

   BOOST_CHECK_CLOSE(MaxRelDiff(M1,M2), 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_MaxRelDiff_matrix_complex)
{
   Eigen::Matrix<std::complex<double>,2,2> M1, M2;
   M1 << +1.0, +1.0, +1.0, +1.0;
   M2 << -1.0, -2.0, -3.0, -4.0;

   BOOST_CHECK_CLOSE(MaxRelDiff(M1,M2), 2.0, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_MaxAbsValue)
{
   BOOST_CHECK_CLOSE(MaxAbsValue(0.)  , 0. , 1e-10);
   BOOST_CHECK_CLOSE(MaxAbsValue(1.)  , 1. , 1e-10);
   BOOST_CHECK_CLOSE(MaxAbsValue(-1.) , 1. , 1e-10);
   BOOST_CHECK_CLOSE(MaxAbsValue(std::complex<double>(0.,1.)), 1. , 1e-10);

   {
      Eigen::Matrix<double,3,1> M;
      M << -2.0, 0.0, 1.0;
      auto M2 = M;

      BOOST_CHECK_CLOSE(MaxAbsValue(M), 2.0, 1e-10);
      BOOST_CHECK_CLOSE(MaxAbsValue(M + M2), 4.0, 1e-10);
   }

   {
      Eigen::Matrix<std::complex<double>,2,1> M;
      M << std::complex<double>(0,1), std::complex<double>(0,2);
      auto M2 = M;

      BOOST_CHECK_CLOSE(MaxAbsValue(M), 2.0, 1e-10);
      BOOST_CHECK_CLOSE(MaxAbsValue(M + M2), 4.0, 1e-10);
   }

   {
      Eigen::Array<double,3,1> A;
      A << -2.0, 0.0, 1.0;
      auto A2 = A;

      BOOST_CHECK_CLOSE(MaxAbsValue(A), 2.0, 1e-10);
      BOOST_CHECK_CLOSE(MaxAbsValue(A + A2), 4.0, 1e-10);
   }

   {
      Eigen::Array<std::complex<double>,2,1> A;
      A << std::complex<double>(0,1), std::complex<double>(0,2);
      auto A2 = A;

      BOOST_CHECK_CLOSE(MaxAbsValue(A), 2.0, 1e-10);
      BOOST_CHECK_CLOSE(MaxAbsValue(A + A2), 4.0, 1e-10);
   }
}

BOOST_AUTO_TEST_CASE(test_Abs_Eigen_Array)
{
   Eigen::ArrayXd v(3);
   v(0) = 1.;
   v(1) = -1.;
   v(2) = 0.;

   Eigen::ArrayXd v_abs(Abs(v));
   BOOST_CHECK_CLOSE(v_abs(0), 1., 1e-10);
   BOOST_CHECK_CLOSE(v_abs(1), 1., 1e-10);
   BOOST_CHECK_CLOSE(v_abs(2), 0., 1e-10);
}

BOOST_AUTO_TEST_CASE(test_Sqr_Eigen_Array)
{
   Eigen::ArrayXd v(3);
   v(0) = 1.;
   v(1) = 2.;
   v(2) = 3.;

   Eigen::ArrayXd v_sqr(Sqr(v));
   BOOST_CHECK_CLOSE(v_sqr(0), 1., 1e-10);
   BOOST_CHECK_CLOSE(v_sqr(1), 4., 1e-10);
   BOOST_CHECK_CLOSE(v_sqr(2), 9., 1e-10);
}

BOOST_AUTO_TEST_CASE(test_Sqrt_overloads)
{
   BOOST_CHECK_EQUAL(Sqrt(2.f), std::sqrt(2.f));
   BOOST_CHECK_EQUAL(Sqrt(2.) , std::sqrt(2.));
   BOOST_CHECK_EQUAL(Sqrt(2.L), std::sqrt(2.L));
   BOOST_CHECK_EQUAL(Sqrt(2)  , std::sqrt(2.));
}

BOOST_AUTO_TEST_CASE(test_Sqrt_Eigen_Array)
{
   Eigen::ArrayXd v(3);
   v(0) = 1.;
   v(1) = 2.;
   v(2) = 3.;

   Eigen::ArrayXd v_sqrt(Sqrt(v));
   BOOST_CHECK_CLOSE(v_sqrt(0), 1., 1e-10);
   BOOST_CHECK_CLOSE(v_sqrt(1), std::sqrt(2.), 1e-10);
   BOOST_CHECK_CLOSE(v_sqrt(2), std::sqrt(3.), 1e-10);
}

BOOST_AUTO_TEST_CASE(test_Total_Array)
{
   Eigen::ArrayXd v(3);
   v(0) = 1.;
   v(1) = 2.;
   v(2) = 3.;

   BOOST_CHECK_CLOSE(Total(v), 6., 1e-10);
}

BOOST_AUTO_TEST_CASE(test_Re_scalars)
{
   BOOST_CHECK_EQUAL(Re(2.), 2.);
   BOOST_CHECK_EQUAL(Re(std::complex<double>(1., 0.)), 1.);
   BOOST_CHECK_EQUAL(Re(std::complex<double>(0., 3.)), 0.);
}

BOOST_AUTO_TEST_CASE(test_Re_Eigen_Matrix)
{
   Eigen::Matrix<double,2,2> m1;
   Eigen::Matrix<std::complex<double>,2,2> m2;

   m1 << 1.0, 2.0, -3.0, 1.0;
   m2 << std::complex<double>(1., 2.), std::complex<double>(-2., 1.),
      std::complex<double>(0., 0.), std::complex<double>(2., 0.);

   Eigen::Matrix<double,2,2> re1(Re(m1));
   Eigen::Matrix<double,2,2> re2(Re(m2));

   BOOST_CHECK_EQUAL(re1(0,0), 1.0);
   BOOST_CHECK_EQUAL(re1(0,1), 2.0);
   BOOST_CHECK_EQUAL(re1(1,0), -3.0);
   BOOST_CHECK_EQUAL(re1(1,1), 1.0);
   BOOST_CHECK_EQUAL(re2(0,0), 1.0);
   BOOST_CHECK_EQUAL(re2(0,1), -2.0);
   BOOST_CHECK_EQUAL(re2(1,0), 0.);
   BOOST_CHECK_EQUAL(re2(1,1), 2.0);
}

BOOST_AUTO_TEST_CASE(test_Im_scalars)
{
   BOOST_CHECK_EQUAL(Im(2.), 0.);
   BOOST_CHECK_EQUAL(Im(std::complex<double>(-3, 5.)), 5.);
   BOOST_CHECK_EQUAL(Im(std::complex<double>(3., 0.)), 0.);
}

BOOST_AUTO_TEST_CASE(test_Im_Eigen_Matrix)
{
   Eigen::Matrix<double,2,2> m1;
   Eigen::Matrix<std::complex<double>,2,2> m2;

   m1 << 1.0, 2.0, 3.0, 4.0;
   m2 << std::complex<double>(1., -2.), std::complex<double>(-2., 1.),
      std::complex<double>(0., 0.), std::complex<double>(0., 10.);

   Eigen::Matrix<double,2,2> im1(Im(m1));
   Eigen::Matrix<double,2,2> im2(Im(m2));

   BOOST_CHECK_EQUAL(im1(0,0), 0.);
   BOOST_CHECK_EQUAL(im1(0,1), 0.);
   BOOST_CHECK_EQUAL(im1(1,0), 0.);
   BOOST_CHECK_EQUAL(im1(1,1), 0.);
   BOOST_CHECK_EQUAL(im2(0,0), -2.0);
   BOOST_CHECK_EQUAL(im2(0,1), 1.0);
   BOOST_CHECK_EQUAL(im2(1,0), 0.);
   BOOST_CHECK_EQUAL(im2(1,1), 10.0);
}

BOOST_AUTO_TEST_CASE(test_Cube)
{
   BOOST_CHECK_EQUAL(Power(2.,3), Cube(2.));
   BOOST_CHECK_EQUAL(Power(3.,3), Cube(3.));
   BOOST_CHECK_EQUAL(Power(4.,3), Cube(4.));
   BOOST_CHECK_EQUAL(Power(2.,3), Power3(2.));
   BOOST_CHECK_EQUAL(Power(3.,3), Power3(3.));
   BOOST_CHECK_EQUAL(Power(4.,3), Power3(4.));
}

BOOST_AUTO_TEST_CASE(test_Quad)
{
   BOOST_CHECK_EQUAL(Power(2.,4), Quad(2.));
   BOOST_CHECK_EQUAL(Power(3.,4), Quad(3.));
   BOOST_CHECK_EQUAL(Power(4.,4), Quad(4.));
   BOOST_CHECK_EQUAL(Power(2.,4), Power4(2.));
   BOOST_CHECK_EQUAL(Power(3.,4), Power4(3.));
   BOOST_CHECK_EQUAL(Power(4.,4), Power4(4.));
}

BOOST_AUTO_TEST_CASE(test_Power5)
{
   BOOST_CHECK_EQUAL(Power(2.,5), Power5(2.));
   BOOST_CHECK_EQUAL(Power(3.,5), Power5(3.));
   BOOST_CHECK_EQUAL(Power(4.,5), Power5(4.));
}

BOOST_AUTO_TEST_CASE(test_Power6)
{
   BOOST_CHECK_EQUAL(Power(2.,6), Power6(2.));
   BOOST_CHECK_EQUAL(Power(3.,6), Power6(3.));
   BOOST_CHECK_EQUAL(Power(4.,6), Power6(4.));
}

BOOST_AUTO_TEST_CASE(test_Power7)
{
   BOOST_CHECK_EQUAL(Power(2.,7), Power7(2.));
   BOOST_CHECK_EQUAL(Power(3.,7), Power7(3.));
   BOOST_CHECK_EQUAL(Power(4.,7), Power7(4.));
}

BOOST_AUTO_TEST_CASE(test_Power8)
{
   BOOST_CHECK_EQUAL(Power(2.,8), Power8(2.));
   BOOST_CHECK_EQUAL(Power(3.,8), Power8(3.));
   BOOST_CHECK_EQUAL(Power(4.,8), Power8(4.));
}

BOOST_AUTO_TEST_CASE(test_Power9)
{
   BOOST_CHECK_EQUAL(Power(2.,9), Power9(2.));
   BOOST_CHECK_EQUAL(Power(3.,9), Power9(3.));
   BOOST_CHECK_EQUAL(Power(4.,9), Power9(4.));
}

BOOST_AUTO_TEST_CASE(test_Power10)
{
   BOOST_CHECK_EQUAL(Power(2.,10), Power10(2.));
   BOOST_CHECK_EQUAL(Power(3.,10), Power10(3.));
   BOOST_CHECK_EQUAL(Power(4.,10), Power10(4.));
}

BOOST_AUTO_TEST_CASE(test_Power11)
{
   BOOST_CHECK_EQUAL(Power(2.,11), Power11(2.));
   BOOST_CHECK_EQUAL(Power(3.,11), Power11(3.));
   BOOST_CHECK_EQUAL(Power(4.,11), Power11(4.));
}

BOOST_AUTO_TEST_CASE(test_Power12)
{
   BOOST_CHECK_EQUAL(Power(2.,12), Power12(2.));
   BOOST_CHECK_EQUAL(Power(3.,12), Power12(3.));
   BOOST_CHECK_EQUAL(Power(4.,12), Power12(4.));
}

template <typename Base, typename Exponent>
double time_Power(std::size_t N, Base b, Exponent e)
{
   Stopwatch sw;
   sw.start();
   for (std::size_t i = 0; i < N; i++)
      volatile auto r = Power(b,e);
   sw.stop();
   return sw.get_time_in_seconds();
}

#define DEF_FUN(f)                              \
   double time_##f(std::size_t N, double b)     \
   {                                            \
      Stopwatch sw;                             \
      sw.start();                               \
      for (std::size_t i = 0; i < N; i++)       \
         volatile auto r = f(b);                \
      sw.stop();                                \
      return sw.get_time_in_seconds();          \
   }

DEF_FUN(Sqr)
DEF_FUN(Cube)
DEF_FUN(Quad)
DEF_FUN(Power5)
DEF_FUN(Power6)
DEF_FUN(Power7)
DEF_FUN(Power8)
DEF_FUN(Power9)
DEF_FUN(Power10)
DEF_FUN(Power11)
DEF_FUN(Power12)

std::string format_line(int exponent, double t1, double t2)
{
   const double diff = MaxRelDiff(t1, t2) * 100;

   auto fmt = boost::format("%|8| | %|23.10| | %|14.10| | %|+8.3|%%")
              % exponent % t1 % t2 % diff;

   return fmt.str();
}

BOOST_AUTO_TEST_CASE(test_Power_benchmark)
{
   const std::size_t N = 1000000000;

   const double timed_Power_2  = time_Power(N, 3., 2);
   const double timed_Power_3  = time_Power(N, 3., 3);
   const double timed_Power_4  = time_Power(N, 3., 4);
   const double timed_Power_5  = time_Power(N, 3., 5);
   const double timed_Power_6  = time_Power(N, 3., 6);
   const double timed_Power_7  = time_Power(N, 3., 7);
   const double timed_Power_8  = time_Power(N, 3., 8);
   const double timed_Power_9  = time_Power(N, 3., 9);
   const double timed_Power_10 = time_Power(N, 3., 10);
   const double timed_Power_11 = time_Power(N, 3., 11);
   const double timed_Power_12 = time_Power(N, 3., 12);
   const double timed_Sqr      = time_Sqr(N, 3.);
   const double timed_Cube     = time_Cube(N, 3.);
   const double timed_Quad     = time_Quad(N, 3.);
   const double timed_Power5   = time_Power5(N, 3.);
   const double timed_Power6   = time_Power6(N, 3.);
   const double timed_Power7   = time_Power7(N, 3.);
   const double timed_Power8   = time_Power8(N, 3.);
   const double timed_Power9   = time_Power9(N, 3.);
   const double timed_Power10  = time_Power10(N, 3.);
   const double timed_Power11  = time_Power11(N, 3.);
   const double timed_Power12  = time_Power12(N, 3.);

   BOOST_TEST_MESSAGE("exponent | time/s (multiplication) | time/s (pow()) | rel. diff");
   BOOST_TEST_MESSAGE("---------------------------------------------------------------");
   BOOST_TEST_MESSAGE(format_line( 2, timed_Sqr    , timed_Power_2 ));
   BOOST_TEST_MESSAGE(format_line( 3, timed_Cube   , timed_Power_3 ));
   BOOST_TEST_MESSAGE(format_line( 4, timed_Quad   , timed_Power_4 ));
   BOOST_TEST_MESSAGE(format_line( 5, timed_Power5 , timed_Power_5 ));
   BOOST_TEST_MESSAGE(format_line( 6, timed_Power6 , timed_Power_6 ));
   BOOST_TEST_MESSAGE(format_line( 7, timed_Power7 , timed_Power_7 ));
   BOOST_TEST_MESSAGE(format_line( 8, timed_Power8 , timed_Power_8 ));
   BOOST_TEST_MESSAGE(format_line( 9, timed_Power9 , timed_Power_9 ));
   BOOST_TEST_MESSAGE(format_line(10, timed_Power10, timed_Power_10));
   BOOST_TEST_MESSAGE(format_line(11, timed_Power11, timed_Power_11));
   BOOST_TEST_MESSAGE(format_line(12, timed_Power12, timed_Power_12));
}

BOOST_AUTO_TEST_CASE(test_Min)
{
   BOOST_CHECK_EQUAL(Min(0.), 0.);
   BOOST_CHECK_EQUAL(Min(1.), 1.);
   BOOST_CHECK_EQUAL(Min(-1), -1);

   BOOST_CHECK_EQUAL(Min(0.,1.), 0.);
   BOOST_CHECK_EQUAL(Min(1.,0.), 0.);
   BOOST_CHECK_EQUAL(Min(1,0.), 0.);
   BOOST_CHECK_EQUAL(Min(1.,0), 0.);

   BOOST_CHECK_EQUAL(Min(-1,0,1), -1);
   BOOST_CHECK_EQUAL(Min(-1.,0.,1.), -1.);
}

BOOST_AUTO_TEST_CASE(test_Max)
{
   BOOST_CHECK_EQUAL(Max(0.), 0.);
   BOOST_CHECK_EQUAL(Max(1.), 1.);
   BOOST_CHECK_EQUAL(Max(-1), -1);

   BOOST_CHECK_EQUAL(Max(0.,1.), 1.);
   BOOST_CHECK_EQUAL(Max(1.,0.), 1.);
   BOOST_CHECK_EQUAL(Max(1,0.), 1.);
   BOOST_CHECK_EQUAL(Max(1.,0), 1.);

   BOOST_CHECK_EQUAL(Max(-1,0,1), 1);
   BOOST_CHECK_EQUAL(Max(-1.,0.,1.), 1.);
}

BOOST_AUTO_TEST_CASE(test_UnitVector)
{
   // 2-vector
   BOOST_CHECK_EQUAL((UnitVector<2,0>())(0), 1.);
   BOOST_CHECK_EQUAL((UnitVector<2,0>())(1), 0.);
   BOOST_CHECK_EQUAL((UnitVector<2,1>())(0), 0.);
   BOOST_CHECK_EQUAL((UnitVector<2,1>())(1), 1.);

   BOOST_CHECK_EQUAL(UnitVector<2>(0)(0), 1.);
   BOOST_CHECK_EQUAL(UnitVector<2>(0)(1), 0.);
   BOOST_CHECK_EQUAL(UnitVector<2>(1)(0), 0.);
   BOOST_CHECK_EQUAL(UnitVector<2>(1)(1), 1.);

   BOOST_CHECK_EQUAL(UnitVector(2,0)(0), 1.);
   BOOST_CHECK_EQUAL(UnitVector(2,0)(1), 0.);
   BOOST_CHECK_EQUAL(UnitVector(2,1)(0), 0.);
   BOOST_CHECK_EQUAL(UnitVector(2,1)(1), 1.);

   // 3-vector
   BOOST_CHECK_EQUAL((UnitVector<3,1>())(0), 0.);
   BOOST_CHECK_EQUAL((UnitVector<3,1>())(1), 1.);
   BOOST_CHECK_EQUAL((UnitVector<3,1>())(2), 0.);

   BOOST_CHECK_EQUAL(UnitVector<3>(1)(0), 0.);
   BOOST_CHECK_EQUAL(UnitVector<3>(1)(1), 1.);
   BOOST_CHECK_EQUAL(UnitVector<3>(1)(2), 0.);

   BOOST_CHECK_EQUAL(UnitVector(3,1)(0), 0.);
   BOOST_CHECK_EQUAL(UnitVector(3,1)(1), 1.);
   BOOST_CHECK_EQUAL(UnitVector(3,1)(2), 0.);
}

BOOST_AUTO_TEST_CASE(test_MatrixProjector)
{
   BOOST_CHECK_EQUAL((MatrixProjector<2,2,0,0>())(0,0), 1.);
   BOOST_CHECK_EQUAL((MatrixProjector<2,2,0,0>())(1,0), 0.);
   BOOST_CHECK_EQUAL((MatrixProjector<2,2,0,0>())(0,1), 0.);
   BOOST_CHECK_EQUAL((MatrixProjector<2,2,0,0>())(1,1), 0.);

   BOOST_CHECK_EQUAL((MatrixProjector<2,2>(0,0))(0,0), 1.);
   BOOST_CHECK_EQUAL((MatrixProjector<2,2>(0,0))(1,0), 0.);
   BOOST_CHECK_EQUAL((MatrixProjector<2,2>(0,0))(0,1), 0.);
   BOOST_CHECK_EQUAL((MatrixProjector<2,2>(0,0))(1,1), 0.);

   BOOST_CHECK_EQUAL(MatrixProjector(2,2,0,0)(0,0), 1.);
   BOOST_CHECK_EQUAL(MatrixProjector(2,2,0,0)(1,0), 0.);
   BOOST_CHECK_EQUAL(MatrixProjector(2,2,0,0)(0,1), 0.);
   BOOST_CHECK_EQUAL(MatrixProjector(2,2,0,0)(1,1), 0.);
}

BOOST_AUTO_TEST_CASE(test_Print_functions)
{
   const int num = 5;

   BOOST_CHECK_EQUAL(PrintVERBOSE("A verbose message: num = ", num), 0.);
   BOOST_CHECK_EQUAL(PrintDEBUG("A debug message: num = ", num), 0.);
   BOOST_CHECK_EQUAL(PrintINFO("An info message: num = ", num), 0.);
   BOOST_CHECK_EQUAL(PrintWARNING("A warning message: num = ", num), 0.);
   BOOST_CHECK_EQUAL(PrintERROR("An error message: num = ", num), 0.);
   BOOST_CHECK_THROW(PrintFATAL("A fatal message: num = ", num), flexiblesusy::FatalError);
}

BOOST_AUTO_TEST_CASE(test_FSThrow)
{
   BOOST_CHECK_THROW(FSThrow("msg"), flexiblesusy::PhysicalError);
}
