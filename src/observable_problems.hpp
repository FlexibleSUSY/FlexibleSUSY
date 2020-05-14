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

#ifndef OBSERVABLE_PROBLEMS_H
#define OBSERVABLE_PROBLEMS_H

#include "observables.hpp"

namespace flexiblesusy {

namespace observable_problems {

/// a_muon problems
class Problem_a_muon {
public:
   /// clears all problems
   void clear();
   /// returns true if there is a problem, false otherwise
   bool have_problem() const;
   /// copies problem strings to output iterator
   template <typename OutputIterator>
   void copy_problem_strings(OutputIterator oi) const {
      if (non_perturbative_running) {
         oi = "non-perturbative running";
         oi++;
      }
   }

   void flag_non_perturbative_running(double);
private:
   bool non_perturbative_running{false};
   double non_perturbative_running_to_scale{0.0};
};

} // namespace observable_problems

class Observable_problems {
public:
   /// clears all problems
   void clear();
   /// returns true if there is a problem, false otherwise
   bool have_problem() const;

   observable_problems::Problem_a_muon a_muon{};
};

} // namespace flexiblesusy

#endif
