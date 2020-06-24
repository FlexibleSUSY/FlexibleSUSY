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

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

namespace flexiblesusy {

namespace observables {

/// all observables supported by FlexibleSUSY
enum EObservables {
   a_muon,
   NUMBER_OF_OBSERVABLES
};

/// observable names
extern const char* const observable_names[NUMBER_OF_OBSERVABLES];

} // namespace observables

} // namespace flexiblesusy

#endif
