// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#ifndef LI5_H
#define LI5_H

#include <complex>

namespace flexiblesusy {

/// complex polylogarithm with n=5
std::complex<double> Li5(const std::complex<double>&) noexcept;

/// complex polylogarithm with n=5 with long double precision
std::complex<long double> Li5(const std::complex<long double>&) noexcept;

} // namespace flexiblesusy

#endif
