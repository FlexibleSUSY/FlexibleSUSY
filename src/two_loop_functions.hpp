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

#ifndef TWO_LOOP_FUNCTIONS
#define TWO_LOOP_FUNCTIONS

#include "cextensions.hpp"

#ifndef PVATTR
#define PVATTR noexcept ATTR(const)
#endif

namespace flexiblesusy {
namespace two_loop_functions {

// tadpole functions

double TfSS(double, double, double, double) PVATTR;
double TfSSS(double, double, double, double) PVATTR;
double TfSSSS(double, double, double, double, double) PVATTR;
double TfSSFF(double, double, double, double, double) PVATTR;
double TfSSFbFb(double, double, double, double, double) PVATTR;
double TfFFFbS(double, double, double, double, double) PVATTR;
double TfFFbFS(double, double, double, double, double) PVATTR;
double TfFbFbFbS(double, double, double, double, double) PVATTR;
double TfSV(double, double) PVATTR;
double TfFV(double, double) PVATTR;

// self-energy functions

double WfSSSS(double, double, double, double, double, double) PVATTR;
double XfSSS(double, double, double, double, double) PVATTR;
double YfSSSS(double, double, double, double, double, double) PVATTR;
double ZfSSSS(double, double, double, double, double, double) PVATTR;
double SfSSS(double, double, double, double, double) PVATTR;
double UfSSSS(double, double, double, double, double, double) PVATTR;
double VfSSSSS(double, double, double, double, double, double, double) PVATTR;
double MfSSSSS(double, double, double, double, double, double, double) PVATTR;
double WfSSSV(double, double, double) PVATTR;
double MfSSSSV(double, double, double, double) PVATTR;
double WfSSFF(double, double, double, double, double, double) PVATTR;
double WfSSFbFb(double, double, double, double, double, double) PVATTR;
double MfFbFbFbFbS(double, double, double, double, double, double, double) PVATTR;
double MfFFbFbFS(double, double, double, double, double, double, double) PVATTR;
double MfFFbFFbS(double, double, double, double, double, double, double) PVATTR;
double MfFFFbFbS(double, double, double, double, double, double, double) PVATTR;
double MfFFFFS(double, double, double, double, double, double, double) PVATTR;
double MfSFbSFbFb(double, double, double, double, double, double, double) PVATTR;
double MfSFSFbF(double, double, double, double, double, double, double) PVATTR;
double MfSFSFFb(double, double, double, double, double, double, double) PVATTR;
double VfSSSFbFb(double, double, double, double, double, double, double) PVATTR;
double VfSSSFF(double, double, double, double, double, double, double) PVATTR;
double VfFbFbFbFbS(double, double, double, double, double, double, double) PVATTR;
double VfFbFFbFS(double, double, double, double, double, double, double) PVATTR;
double VfFbFFFbS(double, double, double, double, double, double, double) PVATTR;
double VfFFbFbFS(double, double, double, double, double, double, double) PVATTR;
double VfFFFbFbS(double, double, double, double, double, double, double) PVATTR;
double VfFFFFS(double, double, double, double, double, double, double) PVATTR;
double GfFFV(double, double, double, double) PVATTR;
double GfFbFbV(double, double, double, double) PVATTR;

} // namespace two_loop_functions
} // namespace flexiblesusy

#endif // TWO_LOOP_FUNCTIONS
