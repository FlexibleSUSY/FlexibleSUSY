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

#ifndef COMPOSITE_CONVERGENCE_TESTER_H
#define COMPOSITE_CONVERGENCE_TESTER_H

#include "convergence_tester.hpp"
#include <vector>

namespace flexiblesusy {

/**
 * @class Composite_convergence_tester
 * @brief A composite convergence tester
 *
 * This class collects convergence testers that will be checked when
 * Composite_convergence_tester::accuracy_goal_reached() is called.
 */

class Composite_convergence_tester : public Convergence_tester {
public:
   Composite_convergence_tester() = default;
   virtual ~Composite_convergence_tester() = default;

   virtual bool accuracy_goal_reached() override;
   virtual int max_iterations() const override;
   virtual void restart() override;
   void add_convergence_tester(Convergence_tester*);

private:
   std::vector<Convergence_tester*> testers{};
};

} // namespace flexiblesusy

#endif
