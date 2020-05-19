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


#include "standard_model_two_scale_low_scale_constraint.hpp"
#include "standard_model_two_scale_model.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "weinberg_angle.hpp"

#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace standard_model {

Standard_model_low_scale_constraint<Two_scale>::Standard_model_low_scale_constraint(
   StandardModel<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : Single_scale_constraint()
   , model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void Standard_model_low_scale_constraint<Two_scale>::apply()
{
   if (!model)
      throw SetupError("Error: Standard_model_low_scale_constraint::apply():"
                       " model pointer must not be zero");

   qedqcd.run_to(scale, 1.0e-5);
   model->calculate_DRbar_masses();

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mz_pole  = qedqcd.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_em > 0)
      delta_alpha_em = model->calculate_delta_alpha_em(alpha_em);

   if (model->get_thresholds() && model->get_threshold_corrections().alpha_s > 0)
      delta_alpha_s  = model->calculate_delta_alpha_s(alpha_s);

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);
   const double g1 = model->get_g1();
   const double g2 = model->get_g2();
   const double mZ = model->get_thresholds() && model->get_threshold_corrections().mz > 0 ?
      model->calculate_MVZ_DRbar(mz_pole) : mz_pole;
   const double theta_w = model->calculate_theta_w(qedqcd, alpha_em_drbar);

   double new_g1 = 1.2909944487358056*e_drbar*Sec(theta_w);
   double new_g2 = e_drbar*Csc(theta_w);

   if (IsFinite(new_g1)) {
      model->get_problems().unflag_non_perturbative_parameter(standard_model_info::g1);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         standard_model_info::g1, new_g1, get_scale());
      new_g1 = Electroweak_constants::g1;
   }

   if (IsFinite(new_g2)) {
      model->get_problems().unflag_non_perturbative_parameter(standard_model_info::g2);
   } else {
      model->get_problems().flag_non_perturbative_parameter(
         standard_model_info::g2, new_g2, get_scale());
      new_g2 = Electroweak_constants::g2;
   }

   model->set_v(Re((2*mZ)/Sqrt(0.6*Sqr(g1) + Sqr(g2))));
   model->calculate_Yu_DRbar(qedqcd);
   model->calculate_Yd_DRbar(qedqcd);
   model->calculate_Ye_DRbar(qedqcd);
   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(3.5449077018110318*Sqrt(alpha_s_drbar));

   if (model->get_thresholds() && model->get_threshold_corrections().sin_theta_w > 0)
      qedqcd.setPoleMW(model->recalculate_mw_pole(qedqcd.displayPoleMW()));
}

double Standard_model_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

void Standard_model_low_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<StandardModel<Two_scale>*>(model_);
}

void Standard_model_low_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& Standard_model_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void Standard_model_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void Standard_model_low_scale_constraint<Two_scale>::initialize()
{
   scale = qedqcd.displayPoleMZ();
}

} // namespace standard_model
} // namespace flexiblesusy
