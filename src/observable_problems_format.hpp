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

#ifndef OBSERVABLE_PROBLEMS_FORMAT_H
#define OBSERVABLE_PROBLEMS_FORMAT_H

#include "observable_problems.hpp"
#include <iosfwd>
#include <string>
#include <boost/format.hpp>

namespace flexiblesusy {

namespace observable_problems {

template <typename T>
class SLHA_observable_problems_output_iterator {
public:
   SLHA_observable_problems_output_iterator(std::ostream& ostr_, int obs_idx_, int flag_)
      : ostr(ostr_), obs_idx(obs_idx_), flag(flag_) {}

   void operator=(const T& elem) {
      ostr << boost::format(" %5d %5d   %s\n") % obs_idx % flag % elem;
   }
   void operator++(int) {}
private:
   std::ostream& ostr;
   int obs_idx; ///< 1st index, observable index
   int flag;    ///< 2nd index, problem type (problem or warning)
};


/// copies problem strings to output iterator
template <typename OutputIterator>
void copy_problem_strings(const Problem_a_muon& p, OutputIterator oi)
{
   if (p.have_non_perturbative_running()) {
      oi = "non-perturbative running to scale "
         + std::to_string(p.get_non_perturbative_running_scale())
         + " GeV";
      oi++;
   }
}

} // namespace observable_problems

/// copies problem strings to output iterator
template <typename OutputIterator>
void copy_problem_strings(const Observable_problems& op, OutputIterator oi)
{
   copy_problem_strings(op.a_muon, oi);
}

} // namespace flexiblesusy

#endif
