#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_legacy_diagonalization

#include <boost/test/unit_test.hpp>

#include "diagonalization.hpp"

#include <random>

using namespace softsusy;

DoubleMatrix random_real_matrix(int n, int m)
{
   static std::default_random_engine generator;
   static std::uniform_real_distribution<> o1(-3, 3);

   DoubleMatrix r(n, m);
   for (int i = 1; i <= n; ++i) {
      for (int j = 1; j <= n; ++j) {
         r(i, j) = o1(generator);
      }
   }
   return r;
}

BOOST_AUTO_TEST_CASE( test_svd )
{
   for (int n = 2; n <= 6; ++n) {
      DoubleMatrix  m(n,n);
      ComplexMatrix u(n,n);
      ComplexMatrix v(n,n);
      DoubleVector  s(n);
      ComplexMatrix diag(n,n);

      for (int count = 100; count; --count) {
         m = random_real_matrix(n,n);
         if (n == 2) {
            flexiblesusy::Diagonalize2by2(m, u, v, s);
         } else {
            flexiblesusy::Diagonalize(m, u, v, s);
         }
         diag = u.complexConjugate() * m * v.hermitianConjugate();

         for (int i = 1; i <= s.displayEnd(); ++i) {
            BOOST_CHECK(s(i) >= 0);
         }
         for (int i = 1; i <= diag.displayCols(); ++i) {
            for (int j = 1; j <= diag.displayRows(); ++j) {
               BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j?s(i):0)), 1e-13);
            }
         }
      }
   }
}

BOOST_AUTO_TEST_CASE( test_symmetric )
{
   for (int n = 2; n <= 6; ++n) {
      DoubleMatrix  m(n,n);
      ComplexMatrix u(n,n);
      DoubleVector  s(n);
      ComplexMatrix diag(n,n);

      for (int count = 100; count; --count) {
         m = random_real_matrix(n,n);
         m.symmetrise();
         if (n == 2) {
            flexiblesusy::Diagonalize2by2(m, u, s);
         } else {
            flexiblesusy::Diagonalize(m, u, s);
         }
         diag = u.complexConjugate() * m * u.hermitianConjugate();

         for (int i = 1; i <= s.displayEnd(); ++i) {
            BOOST_CHECK(s(i) >= 0);
         }
         for (int i = 1; i <= diag.displayCols(); ++i) {
            for (int j = 1; j <= diag.displayRows(); ++j) {
               BOOST_CHECK_SMALL(abs(diag(i,j) - (i==j?s(i):0)), 1e-13);
            }
         }
      }
   }
}
