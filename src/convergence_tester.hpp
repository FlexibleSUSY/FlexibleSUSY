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

#ifndef CONVERGENCE_TESTER_H
#define CONVERGENCE_TESTER_H

namespace flexiblesusy {

class Convergence_tester {
public:
   virtual ~Convergence_tester() = default;
   virtual bool accuracy_goal_reached() = 0;
   virtual int max_iterations() const = 0;
   virtual void restart() = 0;
};

} // namespace flexiblesusy

#endif
