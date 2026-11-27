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

/**
 * Errors and warnings of the decay calculation, written to DCINFO.
 *
 * An error that interrupted the decay calculation (add_fatal_error)
 * leaves the decay table incomplete, so it is a problem that
 * suppresses the decay output unless force output is set.  Other
 * errors (add_error), e.g. a channel set to zero because its width
 * isn't finite, are reported without suppressing the output.
 */
class FlexibleDecay_problems {
public:
   void clear()
   {
      warning_strings.clear();
      error_strings.clear();
      incomplete = false;
   }
   /// true if an error interrupted the decay calculation
   bool have_problem() const { return incomplete; }
   bool have_error() const { return !error_strings.empty(); }
   bool have_warning() const { return !warning_strings.empty(); }
   std::vector<std::string> get_problem_strings() const {
      return error_strings;
   }
   std::vector<std::string> get_warning_strings() const {
      return warning_strings;
   }
   void add_warning(std::string const& warning) {warning_strings.push_back(warning);}
   void add_error(std::string const& error) {error_strings.push_back(error);}
   /// adds an error that interrupted the decay calculation
   void add_fatal_error(std::string const& error)
   {
      error_strings.push_back(error);
      incomplete = true;
   }

private:
   std::vector<std::string> warning_strings {};
   std::vector<std::string> error_strings {};
   bool incomplete{false};
};

} // namespace flexiblesusy

#endif
