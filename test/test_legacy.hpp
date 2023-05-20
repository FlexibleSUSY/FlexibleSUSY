
#ifndef TEST_LEGACY_H
#define TEST_LEGACY_H

#include "error_count.hpp"
#include "linalg.h"
#include "test.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

bool is_equal(const softsusy::DoubleMatrix& a, const softsusy::DoubleMatrix& b, double max_dev)
{
   if (a.displayRows() != b.displayRows()) {
      std::cout << "matrices have different number of rows\n";
      return false;
   }
   if (a.displayCols() != b.displayCols()) {
      std::cout << "matrices have different number of columns\n";
      return false;
   }
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         if (!is_equal(a(i,l), b(i,l), max_dev))
            return false;
      }
   }
   return true;
}

void check_equality(softsusy::Complex a, softsusy::Complex b,
                    const std::string& testMsg, double max_dev)
{
   std::ostringstream msg;
   msg << testMsg << " (real part)";
   check_equality(a.real(), b.real(), msg.str(), max_dev);
   msg.str("");
   msg << testMsg << " (imag part)";
   check_equality(a.imag(), b.imag(), msg.str(), max_dev);
}

void check_equality(const softsusy::DoubleVector& a, const softsusy::DoubleVector& b,
                    const std::string& testMsg, double max_dev)
{
   if (a.displayStart() != b.displayStart()) {
      std::cout << "test failed: " << testMsg << ": "
                << "vectors have different start value: "
                << "a.displayStart() == " << a.displayStart()
                << ", b.displayStart() == " << b.displayStart()
                << std::endl;
      gErrors++;
      return;
   }

   if (a.displayEnd() != b.displayEnd()) {
      std::cout << "test failed: " << testMsg << ": "
                << "vectors have different end value: "
                << "a.displayEnd() == " << a.displayEnd()
                << ", b.displayEnd() == " << b.displayEnd()
                << std::endl;
      gErrors++;
      return;
   }

   for (int i = a.displayStart(); i <= a.displayEnd(); ++i) {
      std::ostringstream element;
      element << testMsg << " [element " << i << "]";
      check_equality(a(i), b(i), element.str(), max_dev);
   }
}

void check_equality(const softsusy::DoubleMatrix& a, const softsusy::DoubleMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const softsusy::DoubleMatrix& b,
                    const Eigen::Matrix<double,3,3>& a,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 0; i < 3; ++i) {
      for (int l = 0; l < 3; ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i+1,l+1), element.str(), max_dev);
      }
   }
}

void check_equality(const softsusy::ComplexMatrix& a, const softsusy::ComplexMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const softsusy::DoubleMatrix& a, const softsusy::ComplexMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 1; i <= a.displayRows(); ++i) {
      for (int l = 1; l <= a.displayCols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "] ";
         check_equality(softsusy::Complex(a(i,l), 0.0), b(i,l), element.str(), max_dev);
      }
   }
}

void check_equality(const softsusy::ComplexMatrix& a, const softsusy::DoubleMatrix& b,
                    const std::string& testMsg, double max_dev)
{
   check_equality(b, a, testMsg, max_dev);
}

} // namespace flexiblesusy

#endif
