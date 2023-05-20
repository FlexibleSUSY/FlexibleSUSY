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

#ifndef SoftsusyMSSM_TWO_SCALE_INITIAL_GUESSER_H
#define SoftsusyMSSM_TWO_SCALE_INITIAL_GUESSER_H

#include "initial_guesser.hpp"
#include "SoftsusyMSSM_parameter_point.hpp"
#include "lowe_legacy.h"

namespace flexiblesusy {

template<class T> class SoftsusyMSSM;
class Two_scale;

class SoftsusyMSSM_low_scale_constraint;
class SoftsusyMSSM_susy_scale_constraint;
class SoftsusyMSSM_sugra_constraint;

class SoftsusyMSSM_initial_guesser : public Initial_guesser {
public:
   SoftsusyMSSM_initial_guesser(SoftsusyMSSM<Two_scale>*, const SoftsusyMSSM_parameter_point&,
                        const SoftsusyMSSM_low_scale_constraint&,
                        const SoftsusyMSSM_susy_scale_constraint&,
                        const SoftsusyMSSM_sugra_constraint&);
   virtual ~SoftsusyMSSM_initial_guesser();
   virtual void guess();
   void set_QedQcd(const softsusy::QedQcd_legacy& qedqcd) { oneset = qedqcd; }

private:
   SoftsusyMSSM<Two_scale>* mssm;     ///< Mssm model
   softsusy::QedQcd_legacy oneset;    ///< low-energy parameters
   SoftsusyMSSM_parameter_point pp;   ///< Mssm parameter point
   bool ewsbBCscale;          ///< EWSB at susy scale
};

}

#endif
