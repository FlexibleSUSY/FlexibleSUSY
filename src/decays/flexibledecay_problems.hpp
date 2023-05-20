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

#ifndef FLEXIBLEDECAY_PROBLEMS_H
#define FLEXIBLEDECAY_PROBLEMS_H

#include <vector>
#include <string>

namespace flexiblesusy {

class FlexibleDecay_problems {
public:
   void clear() {}
   bool have_problem() const { return false; }
   bool have_warning() const { return false; }
   std::vector<std::string> get_problem_strings() const {
      return error_strings;
   }
   std::vector<std::string> get_warning_strings() const {
      return warning_strings;
   }
   void add_warning(std::string const& warning) {warning_strings.push_back(warning);}
   void add_error(std::string const& error) {error_strings.push_back(error);}

private:
   std::string model_name{};
   std::vector<std::string> warning_strings {};
   std::vector<std::string> error_strings {};
};

} // namespace flexiblesusy

#endif
