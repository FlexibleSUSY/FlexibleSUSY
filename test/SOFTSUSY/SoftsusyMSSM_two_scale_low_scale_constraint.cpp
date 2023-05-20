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

#include "SoftsusyMSSM_two_scale_low_scale_constraint.hpp"
#include "SoftsusyMSSM_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

SoftsusyMSSM_low_scale_constraint::SoftsusyMSSM_low_scale_constraint(const SoftsusyMSSM_parameter_point& pp_)
   : Single_scale_constraint()
   , mssm(NULL)
   , scale(softsusy::MZ)
   , pp(pp_)
{
}

SoftsusyMSSM_low_scale_constraint::~SoftsusyMSSM_low_scale_constraint()
{
}

void SoftsusyMSSM_low_scale_constraint::apply()
{
   assert(mssm && "Error: pointer to SoftsusyMSSM<Two_scale> cannot be zero");

   update_scale();
   mssm->sparticleThresholdCorrections(pp.tanBeta);
}

double SoftsusyMSSM_low_scale_constraint::get_scale() const
{
   return scale;
}

void SoftsusyMSSM_low_scale_constraint::set_model(Model* model)
{
   mssm = cast_model<SoftsusyMSSM<Two_scale>*>(model);
}

void SoftsusyMSSM_low_scale_constraint::update_scale()
{
   scale = mssm->displayMz();
}

}
