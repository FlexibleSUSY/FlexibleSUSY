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

#include <cmath>
#include <Eigen/LU>
#include "numerics2.hpp"
#include "linalg2.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include "model.hpp" // for cast_model
#include "MSSMD5O_MSSMRHN_two_scale_matching.hpp"
#include "MSSMD5O_two_scale_model.hpp"
#include "MSSMRHN_two_scale_model.hpp"

#define LowEnergyConstant(p) Electroweak_constants::p

using namespace std;
namespace flexiblesusy {

MSSMD5O_MSSMRHN_matching::MSSMD5O_MSSMRHN_matching() :
    fixed_scale(0),
    lower(0), upper(0),
    inputPars()
{
    // no reasonable guess without input
    scale = initial_scale_guess = 0;
}

MSSMD5O_MSSMRHN_matching::MSSMD5O_MSSMRHN_matching(const MSSMD5O_input_parameters& inputPars_) :
    lower(0), upper(0),
    inputPars(inputPars_)
{
    make_initial_scale_guess();
}

namespace {

struct CompareSpectrum {
    CompareSpectrum(const Eigen::Array3d& s_) : s(s_) {}
    bool operator() (int i, int j) { return s[i] < s[j]; }
    const Eigen::Array3d& s;
};

}

void MSSMD5O_MSSMRHN_matching::invert_seesaw_formula
(const Eigen::Matrix3d& WOp, const Eigen::Vector3d& YvDiag,
 Eigen::Matrix3d& Yv, Eigen::Matrix3d& Mv)
{
    Eigen::Matrix3cd uh;
    Eigen::Array3d s;
    fs_diagonalize_symmetric(WOp, s, uh);
    Eigen::Matrix3d U = uh.adjoint().real();
    Eigen::Vector3d YvDiagInv(1, 1, 1);
    YvDiagInv.array() /= YvDiag.array();
    Eigen::Matrix3d YvInv = U * YvDiagInv.asDiagonal();

    Mv = (YvInv.transpose() * WOp * YvInv).inverse();
    Yv = YvDiag.asDiagonal() * U.adjoint();
}

void MSSMD5O_MSSMRHN_matching::match_low_to_high_scale_model()
{
    upper->set_Yd(lower->get_Yd());
    upper->set_Ye(lower->get_Ye());
    upper->set_Yu(lower->get_Yu());
    upper->set_Mu(lower->get_Mu());
    upper->set_g1(lower->get_g1());
    upper->set_g2(lower->get_g2());
    upper->set_g3(lower->get_g3());
    upper->set_vd(lower->get_vd());
    upper->set_vu(lower->get_vu());

    const auto& WOp = lower->get_WOp();
    const auto YvDiag1 = inputPars.YvDiag1;
    const auto YvDiag2 = inputPars.YvDiag2;
    const auto YvDiag3 = inputPars.YvDiag3;
    Eigen::Vector3d YvDiag;
    YvDiag << YvDiag1, YvDiag2, YvDiag3;
    Eigen::Matrix3d Yv;
    Eigen::Matrix3d Mv;
    invert_seesaw_formula(WOp, YvDiag, Yv, Mv);

    upper->set_Yv(Yv);
    upper->set_Mv(Mv);

    upper->set_TYd(lower->get_TYd());
    upper->set_TYe(lower->get_TYe());
    upper->set_TYu(lower->get_TYu());
    upper->set_BMu(lower->get_BMu());
    upper->set_mq2(lower->get_mq2());
    upper->set_ml2(lower->get_ml2());
    upper->set_mHd2(lower->get_mHd2());
    upper->set_mHu2(lower->get_mHu2());
    upper->set_md2(lower->get_md2());
    upper->set_mu2(lower->get_mu2());
    upper->set_me2(lower->get_me2());
    upper->set_MassB(lower->get_MassB());
    upper->set_MassWB(lower->get_MassWB());
    upper->set_MassG(lower->get_MassG());

    upper->set_scale(lower->get_scale());
}

void MSSMD5O_MSSMRHN_matching::match_high_to_low_scale_model()
{
    update_scale();

    lower->set_Yd(upper->get_Yd());
    lower->set_Ye(upper->get_Ye());
    lower->set_Yu(upper->get_Yu());
    lower->set_Mu(upper->get_Mu());
    lower->set_g1(upper->get_g1());
    lower->set_g2(upper->get_g2());
    lower->set_g3(upper->get_g3());
    lower->set_vd(upper->get_vd());
    lower->set_vu(upper->get_vu());

    const auto& Yv = upper->get_Yv();
    const auto& Mv = upper->get_Mv();
    lower->set_WOp(Yv.transpose() * Mv.inverse() * Yv);

    lower->set_TYd(upper->get_TYd());
    lower->set_TYe(upper->get_TYe());
    lower->set_TYu(upper->get_TYu());
    lower->set_BMu(upper->get_BMu());
    lower->set_mq2(upper->get_mq2());
    lower->set_ml2(upper->get_ml2());
    lower->set_mHd2(upper->get_mHd2());
    lower->set_mHu2(upper->get_mHu2());
    lower->set_md2(upper->get_md2());
    lower->set_mu2(upper->get_mu2());
    lower->set_me2(upper->get_me2());
    lower->set_MassB(upper->get_MassB());
    lower->set_MassWB(upper->get_MassWB());
    lower->set_MassG(upper->get_MassG());

    lower->set_scale(upper->get_scale());
}

double MSSMD5O_MSSMRHN_matching::get_scale() const
{
    return scale;
}

double MSSMD5O_MSSMRHN_matching::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void MSSMD5O_MSSMRHN_matching::set_models(Model *lower_, Model *upper_)
{
    lower = cast_model<MSSMD5O<Two_scale>*>(lower_);
    upper = cast_model<MSSMRHN<Two_scale>*>(upper_);
}

void MSSMD5O_MSSMRHN_matching::set_lower_input_parameters(const MSSMD5O_input_parameters& inputPars_)
{
   inputPars = inputPars_;
   make_initial_scale_guess();
}

void MSSMD5O_MSSMRHN_matching::make_initial_scale_guess()
{
    const auto TanBeta = inputPars.TanBeta;
    double vu = (TanBeta*LowEnergyConstant(vev))/Sqrt(1 + Sqr(TanBeta));
    const auto mv1 = inputPars.mv1;
    const auto mv2 = inputPars.mv2;
    const auto mv3 = inputPars.mv3;
    const auto ThetaV12 = inputPars.ThetaV12;
    const auto ThetaV13 = inputPars.ThetaV13;
    const auto ThetaV23 = inputPars.ThetaV23;

    Eigen::Matrix3d WOp;
    WOp <<
     (2*(mv1*Sqr(Cos(ThetaV12))*Sqr(Cos(ThetaV13)) + mv2*Sqr
      (Cos(ThetaV13))*Sqr(Sin(ThetaV12)) + mv3*Sqr(Sin(ThetaV13))))/Sqr(vu),
     (2*(mv3*Cos(ThetaV13)*Sin(ThetaV13)*Sin(ThetaV23) + mv1
      *Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12)
      *Sin(ThetaV13)*Sin(ThetaV23)) + mv2*Cos(ThetaV13)*Sin(ThetaV12)*(Cos(
      ThetaV12)*Cos(ThetaV23) - Sin(ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23))))/Sqr(
      vu),
     (2*(mv3*Cos(ThetaV13)*Cos(ThetaV23)*Sin(ThetaV13) + mv2
      *Cos(ThetaV13)*Sin(ThetaV12)*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(ThetaV13)) -
      Cos(ThetaV12)*Sin(ThetaV23)) + mv1*Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(
      ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) + Sin(ThetaV12)*Sin(ThetaV23))))/Sqr(
      vu),
     (2*(mv3*Cos(ThetaV13)*Sin(ThetaV13)*Sin(ThetaV23) + mv1
      *Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12)
      *Sin(ThetaV13)*Sin(ThetaV23)) + mv2*Cos(ThetaV13)*Sin(ThetaV12)*(Cos(
      ThetaV12)*Cos(ThetaV23) - Sin(ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23))))/Sqr(
      vu),
     (2*(mv3*Sqr(Cos(ThetaV13))*Sqr(Sin(ThetaV23)) + mv1*Sqr
      (-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23))
      + mv2*Sqr(Cos(ThetaV12)*Cos(ThetaV23) - Sin(ThetaV12)*Sin(ThetaV13)*Sin(
      ThetaV23))))/Sqr(vu),
     (2*(mv1*(-(Cos(ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) +
      Sin(ThetaV12)*Sin(ThetaV23))*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12
      )*Sin(ThetaV13)*Sin(ThetaV23)) + mv2*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(
      ThetaV13)) - Cos(ThetaV12)*Sin(ThetaV23))*(Cos(ThetaV12)*Cos(ThetaV23) - Sin
      (ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23)) + mv3*Cos(ThetaV23)*Sin(ThetaV23)*
      Sqr(Cos(ThetaV13))))/Sqr(vu),
     (2*(mv3*Cos(ThetaV13)*Cos(ThetaV23)*Sin(ThetaV13) + mv2
      *Cos(ThetaV13)*Sin(ThetaV12)*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(ThetaV13)) -
      Cos(ThetaV12)*Sin(ThetaV23)) + mv1*Cos(ThetaV12)*Cos(ThetaV13)*(-(Cos(
      ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) + Sin(ThetaV12)*Sin(ThetaV23))))/Sqr(
      vu),
     (2*(mv1*(-(Cos(ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) +
      Sin(ThetaV12)*Sin(ThetaV23))*(-(Cos(ThetaV23)*Sin(ThetaV12)) - Cos(ThetaV12
      )*Sin(ThetaV13)*Sin(ThetaV23)) + mv2*(-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(
      ThetaV13)) - Cos(ThetaV12)*Sin(ThetaV23))*(Cos(ThetaV12)*Cos(ThetaV23) - Sin
      (ThetaV12)*Sin(ThetaV13)*Sin(ThetaV23)) + mv3*Cos(ThetaV23)*Sin(ThetaV23)*
      Sqr(Cos(ThetaV13))))/Sqr(vu),
     (2*(mv3*Sqr(Cos(ThetaV13))*Sqr(Cos(ThetaV23)) + mv2*Sqr
      (-(Cos(ThetaV23)*Sin(ThetaV12)*Sin(ThetaV13)) - Cos(ThetaV12)*Sin(ThetaV23))
      + mv1*Sqr(-(Cos(ThetaV12)*Cos(ThetaV23)*Sin(ThetaV13)) + Sin(ThetaV12)*Sin(
      ThetaV23))))/Sqr(vu);

    const auto YvDiag1 = inputPars.YvDiag1;
    const auto YvDiag2 = inputPars.YvDiag2;
    const auto YvDiag3 = inputPars.YvDiag3;
    Eigen::Vector3d YvDiag;
    YvDiag << YvDiag1, YvDiag2, YvDiag3;
    Eigen::Matrix3d Yv;
    Eigen::Matrix3d Mv;
    invert_seesaw_formula(WOp, YvDiag, Yv, Mv);

    double RHN_scale = pow(abs(Mv.determinant()), 1.0/3);
    scale = initial_scale_guess = RHN_scale;
}

void MSSMD5O_MSSMRHN_matching::reset()
{
   scale = initial_scale_guess;
   fixed_scale = 0.;
   lower = NULL;
   upper = NULL;
}

void MSSMD5O_MSSMRHN_matching::update_scale()
{
    if (!is_zero(fixed_scale)) {
	scale = fixed_scale;
	return;
    }

    double RHN_scale = pow(abs(upper->get_Mv().determinant()), 1.0/3);
    scale = RHN_scale;
}

MSSMD5O_MSSMRHN_matching_up<Two_scale>::MSSMD5O_MSSMRHN_matching_up()
   : matching()
{}

MSSMD5O_MSSMRHN_matching_up<Two_scale>::MSSMD5O_MSSMRHN_matching_up(const MSSMD5O_input_parameters& input)
   : matching(input)
{}

void MSSMD5O_MSSMRHN_matching_up<Two_scale>::match()
{
   matching.match_low_to_high_scale_model();
}

double MSSMD5O_MSSMRHN_matching_up<Two_scale>::get_scale() const
{
   return matching.get_scale();
}

void MSSMD5O_MSSMRHN_matching_up<Two_scale>::set_models(
   Model *lower, Model *upper)
{
   matching.set_models(lower, upper);
}

double MSSMD5O_MSSMRHN_matching_up<Two_scale>::get_initial_scale_guess() const
{
   return matching.get_initial_scale_guess();
}

void MSSMD5O_MSSMRHN_matching_up<Two_scale>::set_lower_input_parameters(
   const MSSMD5O_input_parameters& input)
{
   matching.set_lower_input_parameters(input);
}

void MSSMD5O_MSSMRHN_matching_up<Two_scale>::set_scale(double scale)
{
   matching.set_scale(scale);
}

void MSSMD5O_MSSMRHN_matching_up<Two_scale>::reset()
{
   matching.reset();
}

MSSMD5O_MSSMRHN_matching_down<Two_scale>::MSSMD5O_MSSMRHN_matching_down()
   : matching()
{}

MSSMD5O_MSSMRHN_matching_down<Two_scale>::MSSMD5O_MSSMRHN_matching_down(const MSSMD5O_input_parameters& input)
   : matching(input)
{}

void MSSMD5O_MSSMRHN_matching_down<Two_scale>::match()
{
   matching.match_high_to_low_scale_model();
}

double MSSMD5O_MSSMRHN_matching_down<Two_scale>::get_scale() const
{
   return matching.get_scale();
}

void MSSMD5O_MSSMRHN_matching_down<Two_scale>::set_models(
   Model *upper, Model *lower)
{
   matching.set_models(lower, upper);
}

double MSSMD5O_MSSMRHN_matching_down<Two_scale>::get_initial_scale_guess() const
{
   return matching.get_initial_scale_guess();
}

void MSSMD5O_MSSMRHN_matching_down<Two_scale>::set_lower_input_parameters(
   const MSSMD5O_input_parameters& input)
{
   matching.set_lower_input_parameters(input);
}

void MSSMD5O_MSSMRHN_matching_down<Two_scale>::set_scale(double scale)
{
   matching.set_scale(scale);
}

void MSSMD5O_MSSMRHN_matching_down<Two_scale>::reset()
{
   matching.reset();
}

} // namespace flexiblesusy
