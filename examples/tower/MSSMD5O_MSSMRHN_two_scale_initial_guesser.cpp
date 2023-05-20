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

#include "MSSMD5O_MSSMRHN_two_scale_initial_guesser.hpp"
#include "MSSMD5O_two_scale_model.hpp"
#include "MSSMRHN_two_scale_model.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <cassert>

namespace flexiblesusy {

#define INPUTPARAMETER(p) input_pars.p
#define MODEL1PARAMETER(p) model_1->get_##p()
#define MODEL2PARAMETER(p) model_2->get_##p()
#define LowEnergyConstant(p) Electroweak_constants::p
// #define MODEL model

MSSMD5O_MSSMRHN_initial_guesser<Two_scale>::MSSMD5O_MSSMRHN_initial_guesser(
   MSSMD5O<Two_scale>* model_1_, MSSMRHN<Two_scale>* model_2_,
   const MSSMD5O_input_parameters& input_pars_,
   const softsusy::QedQcd& qedqcd_,
   const MSSMD5O_low_scale_constraint<Two_scale>& low_constraint_1_,
   const MSSMD5O_susy_scale_constraint<Two_scale>& susy_constraint_1_,
   const MSSMRHN_high_scale_constraint<Two_scale>& high_constraint_2_,
   const MSSMD5O_MSSMRHN_matching_up<Two_scale>& matching_up_,
   const MSSMD5O_MSSMRHN_matching_down<Two_scale>& matching_down_
)
   : Initial_guesser()
   , model_1(model_1_), model_2(model_2_)
   , input_pars(input_pars_)
   , qedqcd(qedqcd_)
   , low_constraint_1(low_constraint_1_)
   , susy_constraint_1(susy_constraint_1_)
   , high_constraint_2(high_constraint_2_)
   , matching_up(matching_up_)
   , matching_down(matching_down_)
{
   assert(model_1 && model_2 && "MSSMD5O_MSSMRHN_initial_guesser: Error: pointers to models must not be zero");
}

MSSMD5O_MSSMRHN_initial_guesser<Two_scale>::~MSSMD5O_MSSMRHN_initial_guesser()
{
}

void MSSMD5O_MSSMRHN_initial_guesser<Two_scale>::guess()
{
   guess_susy_parameters();
   guess_soft_parameters();
}

void MSSMD5O_MSSMRHN_initial_guesser<Two_scale>::guess_susy_parameters()
{
   using namespace softsusy;

   QedQcd leAtMt(qedqcd);
   const double mtpole = leAtMt.displayPoleMt();

   // guess gauge couplings at mt
   const auto alpha_sm(leAtMt.guess_alpha_SM5(mtpole));

   model_1->set_g1(sqrt(4.0 * M_PI * alpha_sm(0)));
   model_1->set_g2(sqrt(4.0 * M_PI * alpha_sm(1)));
   model_1->set_g3(sqrt(4.0 * M_PI * alpha_sm(2)));
   model_1->set_scale(mtpole);
   model_1->set_loops(2);

   // apply user-defined initial guess at the low scale
   const auto TanBeta = INPUTPARAMETER(TanBeta);

   model_1->set_vd(LowEnergyConstant(vev)/Sqrt(1 + Sqr(TanBeta)));
   model_1->set_vu((TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta)));

   Eigen::Matrix<double,3,3> topDRbar(Eigen::Matrix<double,3,3>::Zero()),
      bottomDRbar(Eigen::Matrix<double,3,3>::Zero()),
      electronDRbar(Eigen::Matrix<double,3,3>::Zero());
   topDRbar(0,0)      = leAtMt.displayMass(mUp);
   topDRbar(1,1)      = leAtMt.displayMass(mCharm);
   topDRbar(2,2)      = leAtMt.displayMass(mTop) - 30.0;
   bottomDRbar(0,0)   = leAtMt.displayMass(mDown);
   bottomDRbar(1,1)   = leAtMt.displayMass(mStrange);
   bottomDRbar(2,2)   = leAtMt.displayMass(mBottom);
   electronDRbar(0,0) = leAtMt.displayMass(mElectron);
   electronDRbar(1,1) = leAtMt.displayMass(mMuon);
   electronDRbar(2,2) = leAtMt.displayMass(mTau);

   Eigen::Matrix<double,3,3> new_Yu, new_Yd, new_Ye;

   const auto vd = MODEL1PARAMETER(vd);
   const auto vu = MODEL1PARAMETER(vu);
   new_Yu = ((1.4142135623730951*topDRbar)/vu).transpose();
   new_Yd = ((1.4142135623730951*bottomDRbar)/vd).transpose();
   new_Ye = ((1.4142135623730951*electronDRbar)/vd).transpose();

   model_1->set_Yu(new_Yu);
   model_1->set_Yd(new_Yd);
   model_1->set_Ye(new_Ye);

   const auto mv1 = INPUTPARAMETER(mv1);
   const auto mv2 = INPUTPARAMETER(mv2);
   const auto mv3 = INPUTPARAMETER(mv3);
   const auto ThetaV12 = INPUTPARAMETER(ThetaV12);
   const auto ThetaV13 = INPUTPARAMETER(ThetaV13);
   const auto ThetaV23 = INPUTPARAMETER(ThetaV23);

   model_1->set_WOp(0, 0, (2*(mv1*Sqr(Cos(ThetaV12))*Sqr(Cos(ThetaV13)) + mv2*Sqr
      (Cos(ThetaV13))*Sqr(Sin(ThetaV12)) + mv3*Sqr(Sin(ThetaV13))))/Sqr(vu));
   model_1->set_WOp(0, 1, (2*(mv3*Cos(ThetaV13)*Sin(ThetaV13)*Sin(ThetaV23) + mv1
      *Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12)
      *Sin(ThetaV13)*Sin(ThetaV23)) + mv2*Cos(ThetaV13)*Sin(ThetaV12)*(Cos(
      ThetaV12)*Cos(ThetaV23) - Sin(ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23))))/Sqr(
      vu));
   model_1->set_WOp(0, 2, (2*(mv3*Cos(ThetaV13)*Cos(ThetaV23)*Sin(ThetaV13) + mv2
      *Cos(ThetaV13)*Sin(ThetaV12)*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(ThetaV13)) -
      Cos(ThetaV12)*Sin(ThetaV23)) + mv1*Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(
      ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) + Sin(ThetaV12)*Sin(ThetaV23))))/Sqr(
      vu));
   model_1->set_WOp(1, 0, (2*(mv3*Cos(ThetaV13)*Sin(ThetaV13)*Sin(ThetaV23) + mv1
      *Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12)
      *Sin(ThetaV13)*Sin(ThetaV23)) + mv2*Cos(ThetaV13)*Sin(ThetaV12)*(Cos(
      ThetaV12)*Cos(ThetaV23) - Sin(ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23))))/Sqr(
      vu));
   model_1->set_WOp(1, 1, (2*(mv3*Sqr(Cos(ThetaV13))*Sqr(Sin(ThetaV23)) + mv1*Sqr
      (-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23))
      + mv2*Sqr(Cos(ThetaV12)*Cos(ThetaV23) - Sin(ThetaV12)*Sin(ThetaV13)*Sin(
      ThetaV23))))/Sqr(vu));
   model_1->set_WOp(1, 2, (2*(mv1*(-(Cos(ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) +
      Sin(ThetaV12)*Sin(ThetaV23))*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12
      )*Sin(ThetaV13)*Sin(ThetaV23)) + mv2*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(
      ThetaV13)) - Cos(ThetaV12)*Sin(ThetaV23))*(Cos(ThetaV12)*Cos(ThetaV23) - Sin
      (ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23)) + mv3*Cos(ThetaV23)*Sin(ThetaV23)*
      Sqr(Cos(ThetaV13))))/Sqr(vu));
   model_1->set_WOp(2, 0, (2*(mv3*Cos(ThetaV13)*Cos(ThetaV23)*Sin(ThetaV13) + mv2
      *Cos(ThetaV13)*Sin(ThetaV12)*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(ThetaV13)) -
      Cos(ThetaV12)*Sin(ThetaV23)) + mv1*Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(
      ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) + Sin(ThetaV12)*Sin(ThetaV23))))/Sqr(
      vu));
   model_1->set_WOp(2, 1, (2*(mv1*(-(Cos(ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) +
      Sin(ThetaV12)*Sin(ThetaV23))*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12
      )*Sin(ThetaV13)*Sin(ThetaV23)) + mv2*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(
      ThetaV13)) - Cos(ThetaV12)*Sin(ThetaV23))*(Cos(ThetaV12)*Cos(ThetaV23) - Sin
      (ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23)) + mv3*Cos(ThetaV23)*Sin(ThetaV23)*
      Sqr(Cos(ThetaV13))))/Sqr(vu));
   model_1->set_WOp(2, 2, (2*(mv3*Sqr(Cos(ThetaV13))*Sqr(Cos(ThetaV23)) + mv2*Sqr
      (-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(ThetaV13)) - Cos(ThetaV12)*Sin(ThetaV23))
      + mv1*Sqr(-(Cos(ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) + Sin(ThetaV12)*Sin(
      ThetaV23))))/Sqr(vu));
}

void MSSMD5O_MSSMRHN_initial_guesser<Two_scale>::guess_soft_parameters()
{
   const double low_scale_guess_1 = low_constraint_1.get_initial_scale_guess();
   const double high_scale_guess_2 =
      high_constraint_2.get_initial_scale_guess();
   const double matching_scale_guess = matching_up.get_initial_scale_guess();

   model_1->run_to(matching_scale_guess);
   matching_up.set_models(model_1, model_2);
   matching_up.match();

   model_2->run_to(high_scale_guess_2);

   // apply high-scale constraint
   high_constraint_2.set_model(model_2);
   high_constraint_2.apply();

   // apply user-defined initial guess at the high scale
   model_2->set_Mu(1.);
   model_2->set_BMu(0.);

   model_2->run_to(matching_scale_guess);
   matching_down.set_models(model_2, model_1);
   matching_down.match();

   model_1->run_to(low_scale_guess_1);

   // apply EWSB constraint
   model_1->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model_1->calculate_DRbar_masses();
   model_1->set_thresholds(3);
   model_1->set_loops(2);
}

} // namespace flexiblesusy
