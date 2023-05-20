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

#ifndef SoftsusyNMSSM_MZ_CONSTRAINT_H
#define SoftsusyNMSSM_MZ_CONSTRAINT_H

#include "single_scale_constraint.hpp"
#include "SoftsusyNMSSM_parameter_point.hpp"

namespace flexiblesusy {

class Two_scale;
template<class T> class SoftsusyNMSSM;

/**
 * @class SoftsusyNMSSM_low_scale_constraint
 * @brief NMSSM low-energy constraint at the Z mass MZ
 *
 * This class represents the low-energy constraint of the NMSSM at the
 * Z mass MZ.  The apply() function calculates the threshold
 * corrections to the gauge and Yukawa couplings.  It is assumed that
 * the NMSSM model class is filled with the low-energy data set (see
 * SoftsusyNMSSMSoftsusy::setData).
 */

class SoftsusyNMSSM_low_scale_constraint : public Single_scale_constraint {
public:
   SoftsusyNMSSM_low_scale_constraint(const SoftsusyNMSSM_parameter_point&);
   virtual ~SoftsusyNMSSM_low_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Model*);

private:
   SoftsusyNMSSM<Two_scale>* mssm;
   double scale;
   SoftsusyNMSSM_parameter_point pp;   ///< SoftsusyNMSSM parameter point

   void update_scale();
};

} // namespace flexiblesusy

#endif
