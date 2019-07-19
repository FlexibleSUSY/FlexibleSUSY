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

#ifndef DIAGONALIZATION_H
#define DIAGONALIZATION_H

#include "linalg.h"

namespace flexiblesusy {

/**
 * @class DiagonalizationError
 * @brief Diagonalization failed
 */
class DiagonalizationError : public Error {
public:
   explicit DiagonalizationError(const std::string& msg) : Error(msg) {}
   virtual ~DiagonalizationError() {}
};

softsusy::DoubleVector AbsSqrt(const softsusy::DoubleVector&);

void Diagonalize(const softsusy::DoubleMatrix&, softsusy::DoubleMatrix& , softsusy::DoubleVector&);
void Diagonalize(const softsusy::DoubleMatrix&, softsusy::ComplexMatrix&, softsusy::DoubleVector&);
void Diagonalize2by2(const softsusy::DoubleMatrix&, softsusy::DoubleMatrix& , softsusy::DoubleVector&);
void Diagonalize2by2(const softsusy::DoubleMatrix&, softsusy::ComplexMatrix&, softsusy::DoubleVector&);

// SVD
void Diagonalize(const softsusy::DoubleMatrix&, softsusy::DoubleMatrix& , softsusy::DoubleMatrix& , softsusy::DoubleVector&);
void Diagonalize(const softsusy::DoubleMatrix&, softsusy::ComplexMatrix&, softsusy::ComplexMatrix&, softsusy::DoubleVector&);
void Diagonalize2by2(const softsusy::DoubleMatrix&, softsusy::ComplexMatrix&, softsusy::ComplexMatrix&, softsusy::DoubleVector&);

double MaxRelDiff(const softsusy::DoubleVector&, const softsusy::DoubleVector&);

void Symmetrize(softsusy::DoubleMatrix&);

inline softsusy::DoubleMatrix Re(const softsusy::DoubleMatrix& m)
{
   return m;
}

inline softsusy::DoubleMatrix Re(const softsusy::ComplexMatrix& m)
{
   return m.real();
}

inline softsusy::ComplexMatrix Transpose(const softsusy::ComplexMatrix& m)
{
   return m.transpose();
}

inline softsusy::DoubleMatrix Transpose(const softsusy::DoubleMatrix& m)
{
   return m.transpose();
}

softsusy::DoubleVector ZeroSqrt(const softsusy::DoubleVector&);

} // namespace flexiblesusy

#endif
