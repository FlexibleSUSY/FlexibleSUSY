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

#ifndef SoftsusyNMSSM_MSUSY_CONSTRAINT_H
#define SoftsusyNMSSM_MSUSY_CONSTRAINT_H

#include "single_scale_constraint.hpp"
#include "SoftsusyNMSSM_parameter_point.hpp"
#include "linalg.h"

namespace flexiblesusy {

class Two_scale;
template<class T> class SoftsusyNMSSM;

/**
 * @class SoftsusyNMSSM_susy_scale_constraint
 * @brief MSSM EWSB constraint at the Susy mass scale
 *
 * This class represents the electroweak symmetry breaking (EWSB)
 * constraint of the MSSM at the Susy mass scale MSusy.  The apply()
 * function calculates the MSSM \f$\overline{DR}\f$ parameters and
 * does the EWSB.
 */

class SoftsusyNMSSM_susy_scale_constraint : public Single_scale_constraint {
public:
   SoftsusyNMSSM_susy_scale_constraint(const SoftsusyNMSSM_parameter_point&);
   virtual ~SoftsusyNMSSM_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Model*);

private:
   SoftsusyNMSSM<Two_scale>* snmssm;
   double scale;
   SoftsusyNMSSM_parameter_point pp;   ///< SoftsusyNMSSM parameter point

   void update_scale();
};

}

#endif
