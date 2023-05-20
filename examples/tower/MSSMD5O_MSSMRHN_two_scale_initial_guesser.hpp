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

#ifndef MSSMD5O_MSSMRHN_TWO_SCALE_INITIAL_GUESSER_H
#define MSSMD5O_MSSMRHN_TWO_SCALE_INITIAL_GUESSER_H

#include "MSSMD5O_MSSMRHN_initial_guesser.hpp"
#include "MSSMD5O_input_parameters.hpp"
#include "MSSMD5O_two_scale_low_scale_constraint.hpp"
#include "MSSMD5O_two_scale_susy_scale_constraint.hpp"
#include "MSSMRHN_two_scale_high_scale_constraint.hpp"
#include "MSSMD5O_MSSMRHN_two_scale_matching.hpp"
#include "initial_guesser.hpp"

#include <sstream>

namespace flexiblesusy {

template <class T>
class MSSMD5O;

template <class T>
class MSSMRHN;

class Two_scale;

template<>
class MSSMD5O_MSSMRHN_initial_guesser<Two_scale> : public Initial_guesser {
public:
   MSSMD5O_MSSMRHN_initial_guesser(MSSMD5O<Two_scale>*, MSSMRHN<Two_scale>*,
				const MSSMD5O_input_parameters&,
                                const softsusy::QedQcd&,
				const MSSMD5O_low_scale_constraint<Two_scale>&,
				const MSSMD5O_susy_scale_constraint<Two_scale>&,
				const MSSMRHN_high_scale_constraint<Two_scale>&,
				const MSSMD5O_MSSMRHN_matching_up<Two_scale>&,
                                const MSSMD5O_MSSMRHN_matching_down<Two_scale>&);
   virtual ~MSSMD5O_MSSMRHN_initial_guesser();
   virtual void guess();

private:
   MSSMD5O<Two_scale>* model_1;
   MSSMRHN<Two_scale>* model_2;
   MSSMD5O_input_parameters input_pars;
   const softsusy::QedQcd qedqcd;
   MSSMD5O_low_scale_constraint<Two_scale> low_constraint_1;
   MSSMD5O_susy_scale_constraint<Two_scale> susy_constraint_1;
   MSSMRHN_high_scale_constraint<Two_scale> high_constraint_2;
   MSSMD5O_MSSMRHN_matching_up<Two_scale> matching_up;
   MSSMD5O_MSSMRHN_matching_down<Two_scale> matching_down;

   void guess_susy_parameters();
   void guess_soft_parameters();
};

} // namespace flexiblesusy

#endif
