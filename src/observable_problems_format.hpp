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

#include "observables.hpp"
#include "observable_problems.hpp"
#include <string>
#include <boost/format.hpp>

namespace flexiblesusy {

namespace observable_problems {

template <class OutputIterator>
class SLHA_output_iterator_adaptor {
public:
   SLHA_output_iterator_adaptor(OutputIterator& oi_) : oi(oi_) {}

   void set_observable_name(const char* obs_name_) { obs_name = obs_name_; }
   void set_observable_index(int obs_idx_) { obs_idx = obs_idx_; }
   void set_flag(int flag_) { flag = flag_; }

   template <typename T>
   void operator=(const T& elem) {
      oi = (boost::format(" %5d %5d   %s: %s") % obs_idx % flag % obs_name % elem).str();
   }
   void operator++(int) { oi++; }
private:
   OutputIterator& oi;
   const char* obs_name{nullptr}; ///< name of observable
   int obs_idx{-1}; ///< 1st index, observable index
   int flag{-1};    ///< 2nd index, problem type (problem or warning)
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


/// copies problem strings to output iterator
template <typename OutputIterator>
void format_problems_and_warnings(const Observable_problems& op, OutputIterator oi)
{
   observable_problems::SLHA_output_iterator_adaptor<OutputIterator> slha_oi(oi);

   slha_oi.set_observable_name(observables::observable_names[observables::a_muon]);
   slha_oi.set_observable_index(observables::a_muon + 1);
   slha_oi.set_flag(3); // problems have index 3
   copy_problem_strings(op.a_muon, slha_oi);
}

} // namespace flexiblesusy

#endif
