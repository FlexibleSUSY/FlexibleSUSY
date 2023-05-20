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

#ifndef SoftsusyMSSM_MZ_CONSTRAINT_H
#define SoftsusyMSSM_MZ_CONSTRAINT_H

#include "single_scale_constraint.hpp"
#include "SoftsusyMSSM_parameter_point.hpp"

namespace flexiblesusy {

class Two_scale;
template<class T> class SoftsusyMSSM;

/**
 * @class SoftsusyMSSM_low_scale_constraint
 * @brief MSSM low-energy constraint at the Z mass MZ
 *
 * This class represents the low-energy constraint of the MSSM at the
 * Z mass MZ.  The apply() function calculates the threshold
 * corrections to the gauge and Yukawa couplings.  It is assumed that
 * the MSSM model class is filled with the low-energy data set (see
 * MssmSoftsusy::setData).
 */

class SoftsusyMSSM_low_scale_constraint : public Single_scale_constraint {
public:
   SoftsusyMSSM_low_scale_constraint(const SoftsusyMSSM_parameter_point&);
   virtual ~SoftsusyMSSM_low_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Model*);

private:
   SoftsusyMSSM<Two_scale>* mssm;
   double scale;
   SoftsusyMSSM_parameter_point pp;   ///< Mssm parameter point

   void update_scale();
};

}

#endif
