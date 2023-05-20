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

#include "SoftsusyMSSM_two_scale_susy_scale_constraint.hpp"
#include "SoftsusyMSSM_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

/**
 * Constructor
 *
 * @param pp_ Mssm parameter point
 */
SoftsusyMSSM_susy_scale_constraint::SoftsusyMSSM_susy_scale_constraint(const SoftsusyMSSM_parameter_point& pp_)
   : Single_scale_constraint()
   , mssm(NULL)
   , scale(pp_.msGuess)
   , pp(pp_)
{
}

SoftsusyMSSM_susy_scale_constraint::~SoftsusyMSSM_susy_scale_constraint()
{
}

void SoftsusyMSSM_susy_scale_constraint::apply()
{
   assert(mssm && "Error: pointer to SoftsusyMSSM<Two_scale> cannot be zero");

   mssm->calcDrBarPars();
   update_scale();
   const double mtrun = mssm->displayDrBarPars().mt;
   mssm->rewsb(pp.signMu, mtrun, pp.get_soft_pars());
}

double SoftsusyMSSM_susy_scale_constraint::get_scale() const
{
   return scale;
}

void SoftsusyMSSM_susy_scale_constraint::set_model(Model* model)
{
   mssm = cast_model<SoftsusyMSSM<Two_scale>*>(model);
}

void SoftsusyMSSM_susy_scale_constraint::update_scale()
{
   mssm->setMsusy(mssm->calcMs());
   scale = mssm->displayMsusy();
}

}
