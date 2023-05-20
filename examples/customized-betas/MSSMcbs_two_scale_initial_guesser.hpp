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

#ifndef MSSMCBS_TWO_SCALE_INITIAL_GUESSER_H
#define MSSMCBS_TWO_SCALE_INITIAL_GUESSER_H

#include "MSSMcbs_initial_guesser.hpp"
#include "MSSMcbs_two_scale_low_scale_constraint.hpp"
#include "CMSSM_two_scale_susy_scale_constraint.hpp"
#include "CMSSM_two_scale_high_scale_constraint.hpp"
#include "initial_guesser.hpp"

#include <sstream>

namespace flexiblesusy {

template <class T>
class MSSMcbs;

class Two_scale;

template<>
class MSSMcbs_initial_guesser<Two_scale> : public Initial_guesser {
public:
   MSSMcbs_initial_guesser(MSSMcbs<Two_scale>*,
                               const softsusy::QedQcd&,
                               const MSSMcbs_low_scale_constraint<Two_scale>&,
                               const CMSSM_susy_scale_constraint<Two_scale>&,
                               const CMSSM_high_scale_constraint<Two_scale>&);
   virtual ~MSSMcbs_initial_guesser();
   virtual void guess();

   void set_running_precision(double p) { running_precision = p; }

private:
   MSSMcbs<Two_scale>* model;
   softsusy::QedQcd qedqcd;
   double mu_guess;
   double mc_guess;
   double mt_guess;
   double md_guess;
   double ms_guess;
   double mb_guess;
   double me_guess;
   double mm_guess;
   double mtau_guess;
   double running_precision;
   MSSMcbs_low_scale_constraint<Two_scale> low_constraint;
   CMSSM_susy_scale_constraint<Two_scale> susy_constraint;
   CMSSM_high_scale_constraint<Two_scale> high_constraint;

   void guess_susy_parameters();
   void guess_soft_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
};

} // namespace flexiblesusy

#endif
