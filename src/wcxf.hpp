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

#ifndef WCXF_H
#define WCXF_H

#include <complex>
#include <string>
#include <utility>
#include <vector>

namespace flexiblesusy
{
namespace wcxf
{

/// names and values of Wilson coefficients
using Values = std::vector<std::pair<std::string, std::complex<double>>>;

/**
 * Writes Wilson coefficients of the weak effective theory (WET) in the
 * flavio basis to a file in the Wilson coefficient exchange format
 * (WCxf, arXiv:1712.05298), which flavio and wilson read.
 *
 * Several observables (b -> s gamma, b -> s l l) write to the same
 * file. If the file already contains Wilson coefficients in the same
 * EFT and basis at exactly the same scale, i.e. from another observable
 * of the same parameter point, the given values are added to them
 * (values of the same name are replaced). Otherwise the file is
 * overwritten.
 *
 * @param file_name name of the file
 * @param scale renormalization scale of the Wilson coefficients
 * @param values names and values of the Wilson coefficients
 */
void update_json_flavio(const std::string& file_name, double scale,
                        const Values& values);

} // namespace wcxf
} // namespace flexiblesusy

#endif
