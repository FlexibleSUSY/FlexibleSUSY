// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#ifndef LI_H
#define LI_H

#include <complex>
#include <cstdint>

namespace flexiblesusy {

/// complex polylogarithm for arbitrary integer n
std::complex<double> Li(int64_t n, const std::complex<double>&) noexcept;

} // namespace flexiblesusy

#endif
