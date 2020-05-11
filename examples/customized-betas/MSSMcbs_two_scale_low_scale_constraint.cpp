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

#include "MSSMcbs_two_scale_low_scale_constraint.hpp"
#include "MSSMcbs_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME CMSSM<Two_scale>

MSSMcbs_low_scale_constraint<Two_scale>::MSSMcbs_low_scale_constraint()
   : Single_scale_constraint()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , qedqcd()
   , MZDRbar(0.)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
   , threshold_corrections_loop_order(2)
{
}

MSSMcbs_low_scale_constraint<Two_scale>::MSSMcbs_low_scale_constraint(
   MSSMcbs<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : Single_scale_constraint()
   , model(model_)
   , qedqcd(qedqcd_)
   , new_g1(0.)
   , new_g2(0.)
   , new_g3(0.)
{
   initialize();
}

MSSMcbs_low_scale_constraint<Two_scale>::~MSSMcbs_low_scale_constraint()
{
}

void MSSMcbs_low_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: MSSMcbs_low_scale_constraint:"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();
   calculate_DRbar_gauge_couplings();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
   MODEL->set_vd((2*MZDRbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
      )));
   MODEL->set_vu((2*MZDRbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta))));


   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);
}

double MSSMcbs_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MSSMcbs_low_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void MSSMcbs_low_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<MSSMcbs<Two_scale>*>(model_);
}

void MSSMcbs_low_scale_constraint<Two_scale>::set_sm_parameters(const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

void MSSMcbs_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   qedqcd = softsusy::QedQcd();
   MZDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
}

void MSSMcbs_low_scale_constraint<Two_scale>::initialize()
{
   initial_scale_guess = LowEnergyConstant(MZ);

   scale = initial_scale_guess;

   MZDRbar = 0.;
   new_g1 = 0.;
   new_g2 = 0.;
   new_g3 = 0.;
}

void MSSMcbs_low_scale_constraint<Two_scale>::update_scale()
{
   scale = LowEnergyConstant(MZ);


}

void MSSMcbs_low_scale_constraint<Two_scale>::calculate_DRbar_gauge_couplings()
{
   assert(qedqcd.get_scale() == get_scale() && "Error: low-energy data"
          " set is not defined at the same scale as the low-energy"
          " constraint.  You need to run the low-energy data set to this"
          " scale!");

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);

   double delta_alpha_em = 0.;
   double alS5DRbar_over_alS5MSbar = 1;
   double zeta_g_QCD_2 = 1;
   double zeta_g_SUSY_2 = 1;

   if (model->get_thresholds()) {
      delta_alpha_em = calculate_delta_alpha_em(alpha_em);
      if (model->get_thresholds() == 1) {
	 // one-loop level test against MSSM_low_scale_constraint
	 alS5DRbar_over_alS5MSbar = calculate_alS5DRbar_over_alS5MSbar(alpha_s);
	 zeta_g_QCD_2 = calculate_zeta_g_QCD_2(alpha_s);
	 zeta_g_SUSY_2 = calculate_zeta_g_SUSY_2(alpha_s);
	 alpha_s /= (1 - (alS5DRbar_over_alS5MSbar-1)
		       - (1/zeta_g_QCD_2-1) - (1/zeta_g_SUSY_2-1));
      }
      else {
      	 alS5DRbar_over_alS5MSbar = calculate_alS5DRbar_over_alS5MSbar(alpha_s);
      	 alpha_s *= alS5DRbar_over_alS5MSbar; // alS5MSbar -> alS5DRbar
      	 zeta_g_QCD_2 = calculate_zeta_g_QCD_2(alpha_s);
      	 alpha_s /= zeta_g_QCD_2;	      // alS5DRbar -> alS6DRbar
      	 zeta_g_SUSY_2 = calculate_zeta_g_SUSY_2(alpha_s);
      	 alpha_s /= zeta_g_SUSY_2;	      // alS6DRbar -> alS6DRbarMSSM
      }
   }

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s;
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);

   // interface variables
   MZDRbar
      = model->calculate_MVZ_DRbar(Electroweak_constants::MZ);
   const double MWDRbar
      = model->calculate_MVWm_DRbar(Electroweak_constants::MW);
   const double AlphaS = alpha_s_drbar;
   const double EDRbar = e_drbar;

   const double ThetaW = ArcSin(Sqrt(1 - Sqr(MWDRbar)/Sqr(MZDRbar)));
   new_g1 = 1.2909944487358056*EDRbar*Sec(ThetaW);
   new_g2 = EDRbar*Csc(ThetaW);
   new_g3 = 3.5449077018110318*Sqrt(AlphaS);

}

double MSSMcbs_low_scale_constraint<Two_scale>::calculate_delta_alpha_em(double alphaEm) const
{
   const double currentScale = model->get_scale();
   const auto MCha = MODELPARAMETER(MCha);
   const auto MFu = MODELPARAMETER(MFu);
   const auto MHpm = MODELPARAMETER(MHpm);
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSe = MODELPARAMETER(MSe);
   const auto MSu = MODELPARAMETER(MSu);

   const double delta_alpha_em_SM = 0.15915494309189535*alphaEm*(
      0.3333333333333333 - 1.7777777777777777*FiniteLog(Abs(MFu(2)/currentScale)))
      ;

   const double delta_alpha_em = 0.15915494309189535*alphaEm*(
      -1.3333333333333333*FiniteLog(Abs(MCha(0)/currentScale)) -
      1.3333333333333333*FiniteLog(Abs(MCha(1)/currentScale)) - 0.3333333333333333
      *FiniteLog(Abs(MHpm(1)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd
      (0)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(1)/currentScale))
      - 0.1111111111111111*FiniteLog(Abs(MSd(2)/currentScale)) -
      0.1111111111111111*FiniteLog(Abs(MSd(3)/currentScale)) - 0.1111111111111111*
      FiniteLog(Abs(MSd(4)/currentScale)) - 0.1111111111111111*FiniteLog(Abs(MSd(5
      )/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(0)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(1)/currentScale)) - 0.3333333333333333*
      FiniteLog(Abs(MSe(2)/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(3
      )/currentScale)) - 0.3333333333333333*FiniteLog(Abs(MSe(4)/currentScale)) -
      0.3333333333333333*FiniteLog(Abs(MSe(5)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSu(0)/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(1
      )/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(2)/currentScale)) -
      0.4444444444444444*FiniteLog(Abs(MSu(3)/currentScale)) - 0.4444444444444444*
      FiniteLog(Abs(MSu(4)/currentScale)) - 0.4444444444444444*FiniteLog(Abs(MSu(5
      )/currentScale)));

   return delta_alpha_em + delta_alpha_em_SM;

}

double MSSMcbs_low_scale_constraint<Two_scale>::calculate_alS5DRbar_over_alS5MSbar(double alS5MSbar) const
{
   // Eq. (7) in Harlander, Mihaila, Steinhauser, PRD72(2005)095009

   const int Nc = 3;
   const double CF = (Sqr(Nc) - 1.0) / (2*Nc);
   const int CA = Nc;
   const double T = 1.0/2;
   const int nf = 5;

   return 1 + alS5MSbar/Pi * CA/12
	    + (model->get_thresholds() < 2 ? 0 :
	       Sqr(alS5MSbar/Pi) * (11.0/72 * Sqr(CA) - 1.0/8 * CF * T * nf))
      ;
}

double MSSMcbs_low_scale_constraint<Two_scale>::calculate_zeta_g_QCD_2(double alS5DRbar) const
{
   const double currentScale = model->get_scale();
   const auto MFu = MODELPARAMETER(MFu);

   // Eqs. (19)-(20) in Harlander, Mihaila, Steinhauser, PRD72(2005)095009
   // expressed in powers of alpha_s^(nf-1) up to O(alpha_s^2)

   double Lt = FiniteLog(Sqr(currentScale/MFu(2)));
   double a1 = - Lt / 6;
   double a2 = - 5.0/24 - Lt * 19/24 + Sqr(Lt) / 36;
   double b1 = - a1;
   double b2 = - a2 + 2 * Sqr(a1);

   return 1 / (1 + b1 * alS5DRbar/Pi
	         + (model->get_thresholds() < 2 ? 0 : b2 * Sqr(alS5DRbar/Pi))
	      );
}

double MSSMcbs_low_scale_constraint<Two_scale>::calculate_zeta_g_SUSY_2(double alS6DRbar) const
{
   const double currentScale = model->get_scale();
   const auto MSd = MODELPARAMETER(MSd);
   const auto MSu = MODELPARAMETER(MSu);
   const auto MGlu = MODELPARAMETER(MGlu);

   const double delta_alpha_s = 0.15915494309189535*alS6DRbar*( - 2*FiniteLog(
      Abs(MGlu/currentScale)) - 0.16666666666666666*FiniteLog(Abs(MSd(0)
      /currentScale)) - 0.16666666666666666*FiniteLog(Abs(MSd(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSd(5)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(0)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(1)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(2)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(3)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(4)/currentScale)) -
      0.16666666666666666*FiniteLog(Abs(MSu(5)/currentScale)));

   // Eq. (26) in Harlander, Mihaila, Steinhauser, PRD72(2005)095009

   const int ns = 6;
   double M = Sqrt(MSu(0)*MSu(5));
   double LM = FiniteLog(Sqr(currentScale/M));
   return 1 /
      (1 + delta_alpha_s
         + (model->get_thresholds() < 2 ? 0 :
	    Sqr(alS6DRbar/Pi) * (
	    85.0/32 - 5.0/48 * ns
	    + LM * (3 - ns/12.0) + Sqr(LM) * Sqr(1.0/2 + ns/12.0)))
      );
}

void MSSMcbs_low_scale_constraint<Two_scale>::calculate_DRbar_yukawa_couplings()
{
   calculate_Yu_DRbar();
   calculate_Yd_DRbar();
   calculate_Ye_DRbar();
}

void MSSMcbs_low_scale_constraint<Two_scale>::calculate_Yu_DRbar()
{
   Eigen::Matrix<double,3,3> topDRbar(Eigen::Matrix<double,3,3>::Zero());
   topDRbar(0,0)      = qedqcd.displayMass(softsusy::mUp);
   topDRbar(1,1)      = qedqcd.displayMass(softsusy::mCharm);
   topDRbar(2,2)      = qedqcd.displayMass(softsusy::mTop);

   if (model->get_thresholds())
      topDRbar(2,2) = model->calculate_MFu_DRbar(qedqcd.displayPoleMt(), 2);

   const auto vu = MODELPARAMETER(vu);
   MODEL->set_Yu(((1.4142135623730951*topDRbar)/vu).transpose());

}

void MSSMcbs_low_scale_constraint<Two_scale>::calculate_Yd_DRbar()
{
   Eigen::Matrix<double,3,3> bottomDRbar(Eigen::Matrix<double,3,3>::Zero());
   bottomDRbar(0,0)   = qedqcd.displayMass(softsusy::mDown);
   bottomDRbar(1,1)   = qedqcd.displayMass(softsusy::mStrange);
   bottomDRbar(2,2)   = qedqcd.displayMass(softsusy::mBottom);

   if (model->get_thresholds())
      bottomDRbar(2,2) = model->calculate_MFd_DRbar(qedqcd.displayMass(softsusy::mBottom), 2);

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Yd(((1.4142135623730951*bottomDRbar)/vd).transpose());

}

void MSSMcbs_low_scale_constraint<Two_scale>::calculate_Ye_DRbar()
{
   Eigen::Matrix<double,3,3> electronDRbar(Eigen::Matrix<double,3,3>::Zero());
   electronDRbar(0,0) = qedqcd.displayMass(softsusy::mElectron);
   electronDRbar(1,1) = qedqcd.displayMass(softsusy::mMuon);
   electronDRbar(2,2) = qedqcd.displayMass(softsusy::mTau);

   if (model->get_thresholds())
      electronDRbar(2,2) = model->calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mTau), 2);

   const auto vd = MODELPARAMETER(vd);
   MODEL->set_Ye(((1.4142135623730951*electronDRbar)/vd).transpose());

}

} // namespace flexiblesusy
