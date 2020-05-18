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

#ifndef OBSERVABLE_PROBLEMS_FORMAT_SLHA_H
#define OBSERVABLE_PROBLEMS_FORMAT_SLHA_H

#include "observables.hpp"
#include "observable_problems.hpp"
#include "observable_problems_format.hpp"
#include "slha_format.hpp"

#include <string>

namespace flexiblesusy {

namespace observable_problems {

template <class OutputIterator>
class SLHA_output_iterator_adaptor {
public:
   SLHA_output_iterator_adaptor(OutputIterator& oi_, const char* obs_name_,
                                int obs_idx_, int flag_)
      : oi(oi_)
      , obs_name(obs_name_)
      , obs_idx(obs_idx_)
      , flag(flag_)
   {
   }

   void operator=(const std::string& elem) {
      oi = FORMAT_OBSINFO(obs_idx, flag, std::string(obs_name) + ": " + elem);
   }

   void operator++(int) { oi++; }

private:
   OutputIterator& oi;
   const char* obs_name{nullptr}; ///< name of observable
   int obs_idx{-1}; ///< 1st index, observable index
   int flag{-1};    ///< 2nd index, problem type (problem or warning)
};

} // namespace observable_problems


/// copies problem strings to output iterator
template <typename OutputIterator>
void format_problems_and_warnings(const Observable_problems& op, OutputIterator oi)
{
   const int problem_flag = 3;
   using SLHA_oi = observable_problems::SLHA_output_iterator_adaptor<OutputIterator>;

   {
      SLHA_oi slha_oi(oi,
                      "general",
                      0,
                      problem_flag);
      copy_problem_strings(op.general, slha_oi);
   }

   {
      SLHA_oi slha_oi(oi,
                      observables::observable_names[observables::a_muon],
                      observables::a_muon + 1,
                      problem_flag);
      copy_problem_strings(op.a_muon, slha_oi);
   }
}

} // namespace flexiblesusy

#endif
