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

#ifndef EFFECTIVE_COUPLINGS_H
#define EFFECTIVE_COUPLINGS_H

#include <complex>

namespace flexiblesusy {

namespace effective_couplings {

// loop functions of hep-ph/9504378, Eq. (53)
std::complex<double> scaling_function(double tau);
std::complex<double> AS0(double tau);
std::complex<double> AS12(double tau);
std::complex<double> AS1(double tau);
std::complex<double> AP12(double tau);

/// QCD corrections to the diphoton decay width
/// for scalars and pseudoscalars, at the scale
/// m_decay / 2
std::complex<double> scalar_diphoton_fermion_loop(
   double m_decay, double m_loop);
std::complex<double> pseudoscalar_diphoton_fermion_loop(
   double m_decay, double m_loop);

} // namespace effective_couplings

} // namespace flexiblesusy

#endif
