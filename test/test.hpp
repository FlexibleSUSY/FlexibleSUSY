
#ifndef TEST_H
#define TEST_H

#include "error_count.hpp"
#include "numerics2.hpp"
#include <cmath>
#include <iostream>
#include <type_traits>
#include <Eigen/Core>

namespace flexiblesusy {

static const double max_dev = 1.0e-12;

template <class DerivedA, class DerivedB>
bool is_equal(const Eigen::MatrixBase<DerivedA>& a,
              const Eigen::MatrixBase<DerivedB>& b, double max_dev)
{
   return (a - b).cwiseAbs().maxCoeff() <= max_dev;
}

template <int M, int N>
double max_diff(const Eigen::Matrix<double,M,N>& a,
                const Eigen::Matrix<double,M,N>& b)
{
   return (a - b).maxCoeff();
}

template <int M, int N>
double max_diff(const Eigen::Matrix<std::complex<double>,M,N>& a,
                const Eigen::Matrix<std::complex<double>,M,N>& b)
{
   return (a.cwiseAbs() - b.cwiseAbs()).maxCoeff();
}

template <typename T>
double max_diff(T a, T b)
{
   return a - b;
}

template <int M, int N>
double max_rel_diff(const Eigen::Matrix<double,M,N>& a,
                    const Eigen::Matrix<double,M,N>& b)
{
   return (a - b).maxCoeff() / a.maxCoeff();
}

template <int M, int N>
double max_rel_diff(const Eigen::Matrix<std::complex<double>,M,N>& a,
                    const Eigen::Matrix<std::complex<double>,M,N>& b)
{
   return (a.cwiseAbs() - b.cwiseAbs()).maxCoeff() / a.cwiseAbs().maxCoeff();
}

template <typename T>
double max_rel_diff(T a, T b)
{
   return (a - b) / a;
}

template <class T, class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
void check_equality(T a, T b, const std::string& testMsg, double max_dev)
{
   if (!is_equal(a, b, T(max_dev))) {
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b
                << " (diff.: " << max_diff(a,b) << ", rel. diff.: "
                << 100. * max_rel_diff(a,b) << "%)\n";
      gErrors++;
   }
}

template <typename T>
void check_equality(std::complex<T> a, std::complex<T> b, const std::string& testMsg, double max_dev)
{
   std::ostringstream msg;
   msg << testMsg << " (real part)";
   check_equality(a.real(), b.real(), msg.str(), max_dev);
   msg.str("");
   msg << testMsg << " (imag part)";
   check_equality(a.imag(), b.imag(), msg.str(), max_dev);
}

void check_equality(int a, int b, const std::string& testMsg, double)
{
   if (a != b) {
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b << ")\n";
      gErrors++;
   }
}

template <class DerivedA, class DerivedB>
void check_equality(const Eigen::MatrixBase<DerivedA>& a,
                    const Eigen::MatrixBase<DerivedB>& b,
                    const std::string& testMsg, double max_dev)
{
   assert(a.rows() == b.rows());
   assert(a.cols() == b.cols());

   for (int i = 0; i < a.rows(); ++i) {
      for (int l = 0; l < a.cols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

template <int NRows, int NCols>
void check_equality(const Eigen::Array<double,NRows,NCols>& a,
                    const Eigen::Array<double,NRows,NCols>& b,
                    const std::string& testMsg, double max_dev)
{
   for (int i = 0; i < NRows; ++i) {
      for (int l = 0; l < NCols; ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_equality(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

void check_condition(bool condition, const std::string& testMsg)
{
   if (!condition) {
      std::cout << "test failed: " << testMsg << "\n";
      gErrors++;
   }
}

template <typename T>
void check_greater_than(T a, T b, const std::string& testMsg)
{
   if (!(a > b)) {
      std::cout << "test failed: " << testMsg << ": " << a << " > "
                << b << "\n";
      gErrors++;
   }
}

template <typename T>
void check_relative_dev(T a, T b, const std::string& testMsg, T max_dev)
{
   if (!is_equal_rel(a, b, max_dev)) {
      std::cout << "test failed: " << testMsg << ": "
                << a << " != " << b
                << " (diff.: " << (a-b) << ", rel. diff.: "
                << 100. * (a-b)/a << "%)\n";
      gErrors++;
   }
}

template <class Derived>
void check_relative_dev(const Eigen::MatrixBase<Derived>& a,
                        const Eigen::MatrixBase<Derived>& b,
                        const std::string& testMsg,
                        typename Derived::Scalar max_dev)
{
   assert(a.rows() == b.rows());
   assert(a.cols() == b.cols());

   for (int i = 0; i < a.rows(); ++i) {
      for (int l = 0; l < a.cols(); ++l) {
         std::ostringstream element;
         element << testMsg << " [element " << i << "," << l << "]";
         check_relative_dev(a(i,l), b(i,l), element.str(), max_dev);
      }
   }
}

} // namespace flexiblesusy

#define S(x) #x
#define S_(x) S(x)
#define S__LINE__ S_(__LINE__)
#define TEST_EQUALITY(a, b) \
   flexiblesusy::check_equality(a, b, "line " S__LINE__ ": " #a " == " #b, max_dev)
#define TEST_CLOSE(a, b, dev) \
   flexiblesusy::check_equality(a, b, "line " S__LINE__ ": " #a " == " #b, dev)
#define TEST_GREATER(a, b) \
   flexiblesusy::check_greater_than(a, b, "line " S__LINE__ ": " #a " > " #b)
#define TEST_CLOSE_REL(a, b, dev) \
   flexiblesusy::check_relative_dev(a, b, "line " S__LINE__ ": " #a " =r= " #b, dev)
#define TEST(condition) \
   flexiblesusy::check_condition(condition, #condition);

#endif
