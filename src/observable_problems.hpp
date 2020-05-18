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

namespace flexiblesusy {

namespace observable_problems {

/// general problems for all observables
class Problem_general {
public:
   /// clears all problems
   void clear();
   /// returns true if there is a problem, false otherwise
   bool have_problem() const;
   /// returns number of problems
   unsigned size() const;

   void flag_non_perturbative_running(double);
   bool have_non_perturbative_running() const;
   double get_non_perturbative_running_scale() const;

   void flag_thrown(const char*);
   bool have_thrown() const;
   const char* get_thrown_message() const;
private:
   bool thrown{false};
   const char* thrown_msg{nullptr};
   bool non_perturbative_running{false};
   double non_perturbative_running_to_scale{0.0};
};

/// a_muon problems
class Problem_a_muon {
public:
   /// clears all problems
   void clear();
   /// returns true if there is a problem, false otherwise
   bool have_problem() const;
   /// returns number of problems
   unsigned size() const;

   void flag_non_perturbative_running(double);
   bool have_non_perturbative_running() const;
   double get_non_perturbative_running_scale() const;
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
   /// returns number of problems
   unsigned size() const;

   /// general problems
   observable_problems::Problem_general general{};
   /// problems for a_muon
   observable_problems::Problem_a_muon a_muon{};
};

} // namespace flexiblesusy

#endif
