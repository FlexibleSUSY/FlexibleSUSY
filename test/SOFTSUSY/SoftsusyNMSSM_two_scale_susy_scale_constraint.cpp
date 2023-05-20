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

#include "SoftsusyNMSSM_two_scale_susy_scale_constraint.hpp"
#include "SoftsusyNMSSM_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

/**
 * Constructor
 *
 * @param pp_ SoftsusyNMSSM parameter point
 */
SoftsusyNMSSM_susy_scale_constraint::SoftsusyNMSSM_susy_scale_constraint(const SoftsusyNMSSM_parameter_point& pp_)
   : Single_scale_constraint()
   , snmssm(NULL)
   , scale(pp_.msGuess)
   , pp(pp_)
{
}

SoftsusyNMSSM_susy_scale_constraint::~SoftsusyNMSSM_susy_scale_constraint()
{
}

void SoftsusyNMSSM_susy_scale_constraint::apply()
{
   assert(snmssm && "Error: pointer to SoftsusyNMSSM<Two_scale> cannot be zero");

   snmssm->calcDrBarPars();
   update_scale();
   const double mtrun = snmssm->displayDrBarPars().mt;
   snmssm->rewsb(pp.signMu, mtrun);
}

double SoftsusyNMSSM_susy_scale_constraint::get_scale() const
{
   return scale;
}

void SoftsusyNMSSM_susy_scale_constraint::set_model(Model* model)
{
   snmssm = cast_model<SoftsusyNMSSM<Two_scale>*>(model);
}

void SoftsusyNMSSM_susy_scale_constraint::update_scale()
{
   snmssm->setMsusy(snmssm->calcMs());
   scale = snmssm->displayMsusy();
}

}
