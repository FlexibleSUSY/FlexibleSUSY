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

#include "observables.hpp"
#include <string>
#include <vector>

namespace flexiblesusy {

namespace observable_problems {

class Problem {
   virtual void clear() = 0;
   virtual bool have_problem() const = 0;
   virtual std::vector<std::string> get_problem_strings() const = 0;
};

/// a_muon problems
class Problem_a_muon : public Problem {
public:
   void clear() override;
   bool have_problem() const override;
   std::vector<std::string> get_problem_strings() const override;

   void flag_non_perturbative_running(double);
private:
   bool non_perturbative_running{false};
   double non_perturbative_running_to_scale{0.0};
};

} // namespace observable_problems

class Observable_problems : observable_problems::Problem {
public:
   void clear() override;
   bool have_problem() const override;
   std::vector<std::string> get_problem_strings() const override;

   observable_problems::Problem_a_muon a_muon{};
};

} // namespace flexiblesusy
