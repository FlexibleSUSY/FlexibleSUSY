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

#ifndef SoftsusyNMSSM_TWO_SCALE_CONVERGENCE_TESTER_H
#define SoftsusyNMSSM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "convergence_tester_drbar.hpp"
#include "SoftsusyNMSSM_two_scale.hpp"

namespace flexiblesusy {

class SoftsusyNMSSM_convergence_tester : public Convergence_tester_DRbar<SoftsusyNMSSM<Two_scale> > {
public:
   SoftsusyNMSSM_convergence_tester(SoftsusyNMSSM<Two_scale>*, double);
   virtual ~SoftsusyNMSSM_convergence_tester();

protected:
   virtual double max_rel_diff() const;

private:
   double sumTol(const SoftsusyNMSSM<Two_scale>&, const SoftsusyNMSSM<Two_scale>&) const;
};

}

#endif
