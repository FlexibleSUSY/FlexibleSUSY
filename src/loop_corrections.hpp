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

// File generated at @DateAndTime@

/**
 * @file loop_corrections.hpp
 * @brief contains struct for selection of loop corrections from the
 * literature
 */

#ifndef LOOP_CORRECTIONS_H
#define LOOP_CORRECTIONS_H

namespace flexiblesusy {

struct Loop_corrections {
   bool higgs_at_as{true};
   bool higgs_ab_as{true};
   bool higgs_at_at{true};
   bool higgs_atau_atau{true};
   bool higgs_at_as_as{true};
   bool higgs_ab_as_as{true};
   bool higgs_at_at_as{true};
   bool higgs_at_at_at{true};
   bool higgs_at_as_as_as{true};
   int higgs_3L_scheme{1};
   int top_qcd{1}; ///< top pole mass QCD corrections
};

} // namespace flexiblesusy

#endif
