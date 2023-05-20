#ifndef RANDOM_MATRIX_H
#define RANDOM_MATRIX_H

#include "config.h"
#include "error.hpp"

#include <Eigen/QR>

#include <cmath>

#ifdef ENABLE_RANDOM
#include <random>

namespace flexiblesusy {

// generate a random complex element of the Ginibre ensemble
template <int Rows, int Cols, class Generator>
void random_ginibre_matrix(
   Eigen::Matrix<std::complex<double>, Rows, Cols>& m, Generator& generator)
{
   std::normal_distribution<double> dist(0., 1.);

   const auto rows = m.rows();
   const auto cols = m.cols();

   for (int j = 0; j < cols; ++j) {
      for (int i = 0; i < rows; ++i) {
         m(i, j) = std::complex<double>(dist(generator), dist(generator)) / std::sqrt(2.);
      }
   }
}

// generate a random member of the circular unitary ensemble using
// the algorithm presented in math-ph/0609050
template <int Rows, int Cols, class Generator>
void random_cue_matrix(
   Eigen::Matrix<std::complex<double>, Rows, Cols>& u, Generator& generator)
{
   const auto rows = u.rows();
   const auto cols = u.cols();

   Eigen::Matrix<std::complex<double>, Rows, Cols> z(rows, cols);
   random_ginibre_matrix(z, generator);

   Eigen::ColPivHouseholderQR<Eigen::Matrix<std::complex<double>, Rows, Cols> > cqr(
      z.colPivHouseholderQr());

   Eigen::Matrix<std::complex<double>, Rows, Cols> q_mat(cqr.householderQ());
   Eigen::Matrix<std::complex<double>, Rows, Cols> r_mat(
      cqr.matrixR().template triangularView<Eigen::Upper>());

   const auto calc_phase = [](const std::complex<double>& z) {
      double r = std::abs(z);
      return r == 0. ? 1 : z / r;
   };

   u = q_mat * (r_mat.diagonal().unaryExpr(calc_phase).asDiagonal());
}

} // namespace flexiblesusy

#else

namespace flexiblesusy {

class DisabledRandomError : public Error {
public:
   explicit DisabledRandomError(const std::string& msg_) : msg(msg_) {}
   virtual ~DisabledRandomError() = default;
   std::string what() const override { return msg; }

private:
   std::string msg;
};

template <int Rows, int Cols, class Generator>
void random_ginibre_matrix(
   Eigen::Matrix<std::complex<double>, Rows, Cols>& /* m */, Generator& /* generator */)
{
   throw DisabledRandomError(
      "cannot call random_ginibre_matrix, because <random> is not available");
}

template <int Rows, int Cols, class Generator>
void random_cue_matrix(
   Eigen::Matrix<std::complex<double>, Rows, Cols>& /* u */, Generator& /* generator */)
{
   throw DisabledRandomError(
      "cannot call random_ginibre_matrix, because <random> is not available");
}

} // namespace flexiblesusy

#endif

#endif
