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

#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include <complex>
#include "src/lowe.h"

namespace flexiblesusy {

// Eq.(2.6) of hep-ph/0503173
double calc_DeltaAH(double b) noexcept;

/// Eq.(2.6) of hep-ph/0503173
double calc_DeltaH(double b) noexcept;

/// Eq.(2.11) of hep-ph/0503173, 2-loop and higher order
double calc_Deltaqq(double alpha_s_red, double Nf) noexcept;

/// Eq.(2.31) of hep-ph/0503172, including edge cases
double RT(double x) noexcept;

// eq. 2.5, 2.7 & 2.8 of https://arxiv.org/pdf/hep-ph/0509189.pdf
std::complex<double> delta_hAA_2loopQCD_for_quark_loop(double mH, double mq, double mu) noexcept;

// eq. 2.17, 2.19 & 2.20 of https://arxiv.org/pdf/hep-ph/0509189.pdf
std::complex<double> delta_AhAA_2loopQCD_for_quark_loop(double mAh, double mq, double mu) noexcept;

std::complex<double> delta_AhAA_2loopQCD_for_squark_loop(double mAH, double msq, double mu) noexcept;

/* 2-loop QCD corrections to H->gamma gamma amplitude through scalar color triplet loop
 * from hep-ph/0611266 */

/* this functions could be improved to cover all range in r by linking
 * to CHAPLIN library (https://arxiv.org/abs/1106.5739) */
std::complex<double> delta_hAA_2loopQCD_for_squark_loop(double mH, double msq, double mu) noexcept;

// eq. 6 of https://arxiv.org/pdf/1109.5304.pdf
std::complex<double> f(double x) noexcept;

unsigned int number_of_active_flavours(softsusy::QedQcd const&, double m) noexcept;
double sm_up_quark_masses(softsusy::QedQcd const&, int);
double sm_down_quark_masses(softsusy::QedQcd const&, int);

} // namespace flexiblesusy

#endif
