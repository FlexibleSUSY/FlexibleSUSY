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

#include "weinberg_angle.hpp"
#include "config.h"
#include "ew_input.hpp"
#include "logger.hpp"
#include "numerics.h"
#include "numerics2.hpp"
#include "sm_mw.hpp"
#include "wrappers.hpp"
#include "sum.hpp"
#include "loop_libraries/loop_library.hpp"

#include <cmath>
#include <complex>
#include <limits>

#define WARN_IF_ZERO(p,fun)                     \
   if (is_zero(p))                              \
      WARNING(#fun ": " #p " is zero!");

#define QUIT_IF(condition,fun)                  \
   if (condition)                               \
      FATAL(#fun ": " #condition);

#define ROOT2 Electroweak_constants::root2

#define MODELPARAMETER(p) model->get_##p()

namespace flexiblesusy {

namespace {

double sqr(double x) noexcept { return x*x; }

double calculate_G_fermi(double alpha_em, double mz_pole, double sin_2_cos_2, double delta_rhat) noexcept
{
   return Pi*alpha_em/(ROOT2*Sqr(mz_pole)*sin_2_cos_2*(1.0 - delta_rhat));
}

} // anonymous namespace

namespace weinberg_angle {

using namespace softsusy;

Weinberg_angle::Weinberg_angle(
   const standard_model::Standard_model* model_,
   const Sm_parameters& sm_parameters_)
   : model(model_)
   , sm_parameters(sm_parameters_)
{
}

void Weinberg_angle::set_number_of_iterations(int n)
{
   number_of_iterations = n;
}

void Weinberg_angle::set_number_of_loops(int n)
{
   number_of_loops = n;
}

void Weinberg_angle::set_precision_goal(double p)
{
   precision_goal = p;
}

double Weinberg_angle::get_rho_hat() const
{
   return rho_hat;
}

/**
 * Calculates the DR-bar weak mixing angle \f$\sin\hat{\theta}_W\f$ as
 * defined in Eq. (C.3) from hep-ph/9606211 given the Fermi constant,
 * the Z-boson pole mass and the DR-bar electromagnetic coupling as input
 * and taking the tree-level value of the \f$\hat{\rho}\f$ parameter into account.
 * Furthermore the W boson pole mass is determined from the final result.
 *
 * The function throws an exception of type NoSinThetaWConvergenceError if the
 * iterative procedure to determine the weak mixing angle does not converge.
 *
 * @param sinThetaW_start initial guess for the sine of the weak mixing angle
 *
 * @return sine of the DR-bar weak mixing angle (#1) and W pole mass (#2)
 */
double Weinberg_angle::calculate(double sinThetaW_start)
{
   const auto gY = model->get_g1() * standard_model_info::normalization_g1;
   const auto g2 = model->get_g2() * standard_model_info::normalization_g2;
   const double eDRbar     = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(eDRbar) / (4.0 * Pi);
   const double mw         = sm_parameters.mw_pole;
   const double mz         = sm_parameters.mz_pole;
   const double gfermi     = sm_parameters.fermi_constant;

   pizzt_MZ = calculate_self_energy_VZ(mz);
   piwwt_MW = calculate_self_energy_VWp(mw);
   piwwt_0  = calculate_self_energy_VWp(0.);

   int iteration = 0;
   bool not_converged = true;
   bool fudged = false;
   double sinThetaW_old = sinThetaW_start;
   double sinThetaW_new = sinThetaW_start;

   while (not_converged && iteration < number_of_iterations) {
      fudged = false;

      double deltaRhoHat = calculate_delta_rho_hat(sinThetaW_old);

      if (!std::isfinite(deltaRhoHat) || Abs(deltaRhoHat) >= 1.0) {
         fudged = true;
         deltaRhoHat = 0.;
      }

      const double rhohat_ratio = 1.0 / (1.0 - deltaRhoHat);

      double deltaRHat = calculate_delta_r_hat(rhohat_ratio, sinThetaW_old);

      if (!std::isfinite(deltaRHat) || Abs(deltaRHat) >= 1.0) {
         fudged = true;
         deltaRHat = 0.;
      }

      double sin2thetasqO4 = Pi * alphaDRbar /
         (ROOT2 * Sqr(mz) * gfermi * (1.0 - deltaRHat));

      if (sin2thetasqO4 >= 0.25) {
         fudged = true;
         sin2thetasqO4 = 0.25;
      }

      if (sin2thetasqO4 < 0.0) {
         fudged = true;
         sin2thetasqO4 = 0.0;
      }

      const double sin2theta = Sqrt(4.0 * sin2thetasqO4);
      const double theta = 0.5 * ArcSin(sin2theta);

      sinThetaW_new = Sin(theta);

      const double precision = Abs(sinThetaW_old / sinThetaW_new - 1.0);

      VERBOSE_MSG("\t\tIteration step " << iteration
                  << ": prec=" << precision
                  << " dRhoHat=" << deltaRhoHat
                  << " rhohat_ratio=" << rhohat_ratio
                  << " dRHat=" << deltaRHat
                  << " sinThetaW_new=" << sinThetaW_new
                  << " fudged = " << fudged);

      not_converged = precision >= precision_goal;

      sinThetaW_old = sinThetaW_new;
      iteration++;
   }

   if (fudged)
      throw NonPerturbativeSinThetaW();

   if (not_converged)
      throw NoSinThetaWConvergenceError(number_of_iterations, sinThetaW_new);

   return sinThetaW_new;
}

double Weinberg_angle::calculate_G_fermi()
{
   const double gY = model->get_g1() * standard_model_info::normalization_g1;
   const double g2 = model->get_g2() * standard_model_info::normalization_g2;
   const double e = gY*g2/Sqrt(Sqr(gY) + Sqr(g2));
   const double alpha_em_drbar = Sqr(e)/(4*Pi); // @todo(alex): may use input sm_parameters.alpha_em_drbar from qedqcd here
   const double mw_drbar = model->get_MVWp();
   const double mz_drbar = model->get_MVZ();
   const double cos_theta = mw_drbar/mz_drbar;
   const double theta = ArcCos(cos_theta);
   const double mz_pole = sm_parameters.mz_pole;
   const double sin_theta = Sin(theta);
   const double sin_2_cos_2 = Sqr(sin_theta*Cos(theta));

   if (number_of_loops <= 0) {
      return flexiblesusy::calculate_G_fermi(alpha_em_drbar, mz_pole, sin_2_cos_2, 0.0);
   }

   if (number_of_loops <= 1){
      const double delta_rho_hat = calculate_delta_rho_hat(sin_theta);
      const double rhohat_ratio = 1.0/(1.0 - delta_rho_hat); // Eq.(C.4) [arXiv:hep-ph/9606211]
      const double delta_r_hat = calculate_delta_r_hat(rhohat_ratio, sin_theta);
      return flexiblesusy::calculate_G_fermi(alpha_em_drbar, mz_pole, sin_2_cos_2, delta_r_hat);
   }

   // int iteration = 0;
   // bool not_converged = true;

   // double rhohat_1l = calculate_rho_hat(sin_theta, data, susy_contributions, 1);
   // double delta_r_hat_1l = calculate_delta_r(rhohat_1l, sin_theta, data, susy_contributions,1);
   // double gfermi_1l = Pi*alphaDRbar/(ROOT2*Sqr(mz_pole)*sin_2_cos_2*(1.0 - delta_r_hat_1l));

   // data.fermi_contant = gfermi_1l;
   // double gfermi_old = gfermi_1l;
   // double gfermi_new = gfermi_1l;

   // while (not_converged && iteration < number_of_iterations) {
   //    rhohat = calculate_rho_hat(sin_theta, data, susy_contributions, number_of_loops);
   //    delta_rhat = calculate_delta_r(rhohat, sin_theta, data, susy_contributions, number_of_loops);
   //    gfermi_new = Pi*alphaDRbar/(ROOT2*Sqr(mz_pole)*sin_2_cos_2*(1.0 - delta_rhat));
   //    const double precision = Abs(gfermi_old / gfermi_new - 1.0);
   //    not_converged = precision >= precision_goal;
   //    data.fermi_contant = gfermi_new;
   //    gfermi_old = gfermi_new;
   //    iteration++;
   // }

   // rho_hat = rhohat;
   // g_fermi = gfermi_new;

   // const int no_convergence_error = iteration == number_of_iterations;

   return 0.0;
}

/**
 * Calculates the vertex, box and external wave-function renormalization
 * corrections \f$\delta_{\text{VB}}\f$ for the specific model as e.g.
 * given in Eqs. (C.11)-(C.16), (C.20) from hep-ph/9606211 for the MSSM.
 *
 * @param rhohat_ratio = rhohat / rhohat_tree
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\delta_{\text{VB}}\f$
 */
double Weinberg_angle::calculate_delta_vb(double rhohat_ratio, double sinThetaW) const
{
   const double deltaVbSM = calculate_delta_vb_sm(sinThetaW);

   const double deltaVb = rhohat_ratio * deltaVbSM;

   return deltaVb;
}

/**
 * Calculates \f$\rho^{(2)}(r)\f$ as given in Eqs. (C.7)-(C.8) from
 * hep-ph/9606211 .
 *
 * @param r ratio of Higgs mass over top quark mass
 *
 * @return \f$\rho^{(2)}(r)\f$
 */
double Weinberg_angle::rho_2(double r)
{
   const double Pi2 = Pi * Pi;

   if (r <= std::numeric_limits<double>::epsilon()) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
      WARNING("rho_2: value of r is invalid: r = " << r);
      WARNING("-> setting 2-loop corrections ~ xt^2 to 0");
#endif
      return 0.;
   }

   if (r <= 1.9) {
      const double r2 = sqr(r);
      return 19.0 - 16.5 * r + 43.0 * r2 / 12.0 + 7.0 / 120.0 * r2 * r -
         Pi * sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * r2 + r2 * r /
                         256.0) - Pi2 * (2.0 - 2.0 * r + 0.5 * r2) -
         std::log(r) * (3.0 * r - 0.5 * r2);
   } else {
      const double rm1 = 1.0 / r, rm2 = sqr(rm1), rm3 = rm2 * rm1,
         rm4 = rm3 * rm1, rm5 = rm4 * rm1;
      return sqr(std::log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 - 168.0
                            * rm4 - 612.0 * rm5) -
         std::log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
                   8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5)
         + Pi2 * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
                      56.0 * rm4 - 204.0 * rm5)
         + 49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
         341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
   }
}

/**
 * Calculates 1-loop transverse Z boson self-energy
 * including the correction from usage of pole instead of DRbar top-quark mass.
 *
 * @param p momentum
 *
 * @return 1-loop transverse Z boson self-energy
 */
double Weinberg_angle::calculate_self_energy_VZ(double p) const
{
   const double mt      = sm_parameters.mt_pole;
   const double mtDRbar = model->get_MFu(2);
   const auto pizzt   = Re(model->self_energy_VZ_1loop(p));

   double pizzt_corrected = pizzt;

   if (model->get_thresholds() > 1) {
      pizzt_corrected =
         pizzt - calculate_self_energy_VZ_top(p, mtDRbar)
               + calculate_self_energy_VZ_top(p, mt);
   }

   return pizzt_corrected;
}

/**
 * Calculates 1-loop top-quark contribution to Z boson self-energy.
 *
 * @param p momentum
 * @param mt top-quark mass
 *
 * @return 1-loop top-quark contribution to Z boson self-energy
 */
double Weinberg_angle::calculate_self_energy_VZ_top(double p, double mt) const
{
   const double q  = model->get_scale();
   static constexpr double Nc = 3.0;
   const auto gY = model->get_g1() * standard_model_info::normalization_g1;
   const auto g2 = model->get_g2() * standard_model_info::normalization_g2;
   const double gY2 = Sqr(gY);
   const double g22 = Sqr(g2);
   const double sw2 = gY2 / (gY2 + g22);
   const double cw2 = 1.0 - sw2;
   const double guL = 0.5 - 2.0 * sw2 / 3.0;
   const double guR = 2.0 * sw2 / 3.0;

   const double self_energy_z_top =
      Nc * Sqr(g2) / cw2 * oneOver16PiSqr *
      (softsusy::hfn(p, mt, mt, q) * (Sqr(guL) + Sqr(guR)) -
       4.0 * guL * guR * Sqr(mt) * softsusy::b0(p, mt, mt, q));

   return self_energy_z_top;
}

/**
 * Calculates 1-loop transverse W boson self-energy
 * including the correction from usage of pole instead of DRbar top-quark mass.
 *
 * @param p momentum
 *
 * @return 1-loop transverse W boson self-energy
 */
double Weinberg_angle::calculate_self_energy_VWp(double p) const
{
   const double mt      = sm_parameters.mt_pole;
   const double mtDRbar = model->get_MFu(2);
   const auto piwwt   = Re(model->self_energy_VWp_1loop(p));

   double piwwt_corrected = piwwt;

   if (model->get_thresholds() > 1) {
      piwwt_corrected =
         piwwt - calculate_self_energy_VWp_top(p, mtDRbar)
               + calculate_self_energy_VWp_top(p, mt);
   }

   return piwwt_corrected;
}

/**
 * Calculates 1-loop top-quark contribution to W boson self-energy.
 *
 * @param p momentum
 * @param mt top-quark mass
 *
 * @return 1-loop top-quark contribution to W boson self-energy
 */
double Weinberg_angle::calculate_self_energy_VWp_top(double p, double mt) const
{
   const double q  = model->get_scale();
   const double mb = model->get_MFd(2);
   static constexpr double Nc = 3.0;
   const auto g2 = model->get_g2() * standard_model_info::normalization_g2;

   const double self_energy_w_top =
      0.5 * Nc * softsusy::hfn(p, mt, mb, q) * Sqr(g2) * oneOver16PiSqr;

   return self_energy_w_top;
}

/**
 * Calculates the \f$\Delta\hat{\rho}\f$ corrections as defined in
 * Eqs. (C.4), (C.6) from hep-ph/9606211 but with the dependency on
 * rhohat eliminated.
 *
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\Delta\hat{\rho}\f$ as defined in (C.4) and (C.6) from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_rho_hat(double sinThetaW) const
{
   const double gfermi = sm_parameters.fermi_constant;
   const double mw = sm_parameters.mw_pole;
   const double mz = sm_parameters.mz_pole;
   const double mt = sm_parameters.mt_pole;
   const double alphaS = sm_parameters.alpha_s;

   const double deltaRhoHat1Loop = number_of_loops > 0 ?
      1 - (1 + piwwt_MW / Sqr(mw)) / (1 + pizzt_MZ / Sqr(mz)) : 0.;

   std::complex<double> deltaRhoHat2LoopSM(0., 0.);

   if (number_of_loops > 1) {
      if (Abs(1. - model->get_scale() / mz) > 0.1) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
         WARNING("scale deviates from value mz_pole assumed in deltaRhoHat2LoopSM");
#endif
      }

      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto v = MODELPARAMETER(v);
      const auto Yu = MODELPARAMETER(Yu);
      const auto Vu = MODELPARAMETER(Vu);
      const auto Uu = MODELPARAMETER(Uu);
      const auto Mhh = MODELPARAMETER(Mhh);

      deltaRhoHat2LoopSM = (6.015223977354103e-6*((64*alphaS*Pi*Sqr(g1)*Sqr(g2)*(-
         2.24 + 1.262*Log(mt/mz) - (2.145*Sqr(mt))/Sqr(mw) - (0.85*Sqr(mz))/Sqr(mt)))
         /((0.6*Sqr(g1) + Sqr(g2))*Sqr(sinThetaW)) - 5*rho_2(Mhh/mt)*Sqr(gfermi)*Sqr(
         mt)*Sqr(v)*(Sqr(Abs(SUM(j2,0,2,Conj(Vu(2,j2))*SUM(j1,0,2,Conj(Uu(2,j1))*Yu(
         j1,j2))) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(2,j1))*Vu(2,j2)))) - Sqr
         (Abs(SUM(j2,0,2,Conj(Vu(2,j2))*SUM(j1,0,2,Conj(Uu(2,j1))*Yu(j1,j2))) + SUM(
         j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(2,j1))*Vu(2,j2)))))))/(1 + pizzt_MZ/Sqr(
         mz));   }

   const double deltaRhoHat2LoopSMreal = std::real(deltaRhoHat2LoopSM);

   const double deltaRhoHat = deltaRhoHat1Loop + deltaRhoHat2LoopSMreal;

   return deltaRhoHat;
}

/**
 * Calculates the \f$\Delta\hat{r}\f$ corrections as defined in
 * Eqs. (C.3), (C.5) from hep-ph/9606211 taking the tree-level
 * value of the rhohat parameter into account.
 *
 * @param rhohat_ratio = rhohat / rhohat_tree
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\Delta\hat{r}\f$ as defined in (C.3) and (C.5) from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_r_hat(double rhohat_ratio, double sinThetaW) const
{
   const double gfermi = sm_parameters.fermi_constant;
   const double mw = sm_parameters.mw_pole;
   const double mz = sm_parameters.mz_pole;
   const double mt = sm_parameters.mt_pole;
   const double alphaS = sm_parameters.alpha_s;

   const double dvb = number_of_loops > 0 ?
      calculate_delta_vb(rhohat_ratio, sinThetaW) : 0.;

   const double deltaRHat1Loop = number_of_loops > 0 ?
      rhohat_ratio * piwwt_0 / Sqr(mw) - pizzt_MZ / Sqr(mz) + dvb : 0.;

   std::complex<double> deltaRHat2LoopSM(0., 0.);

   if (number_of_loops > 1) {
      if (Abs(1. - model->get_scale() / mz) > 0.1) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
         WARNING("scale deviates from value mz_pole assumed in deltaRHat2LoopSM");
#endif
      }

      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto v = MODELPARAMETER(v);
      const auto Yu = MODELPARAMETER(Yu);
      const auto Vu = MODELPARAMETER(Vu);
      const auto Uu = MODELPARAMETER(Uu);
      const auto Mhh = MODELPARAMETER(Mhh);

      deltaRHat2LoopSM = 6.015223977354103e-6*((64*alphaS*Pi*Sqr(g1)*Sqr(g2)*(-0.224
         + 0.575*Log(mt/mz) + (2.145*Sqr(mt))/Sqr(mz) - (0.144*Sqr(mz))/Sqr(mt)))/((
         0.6*Sqr(g1) + Sqr(g2))*(1 - Sqr(sinThetaW))*Sqr(sinThetaW)) - 5*(-1 +
         deltaRHat1Loop)*rhohat_ratio*rho_2(Mhh/mt)*Sqr(gfermi)*Sqr(mt)*Sqr(v)*(Sqr(Abs
         (SUM(j2,0,2,Conj(Vu(2,j2))*SUM(j1,0,2,Conj(Uu(2,j1))*Yu(j1,j2))) - SUM(j2,0,
         2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(2,j1))*Vu(2,j2)))) - Sqr(Abs(SUM(j2,0,2,Conj
         (Vu(2,j2))*SUM(j1,0,2,Conj(Uu(2,j1))*Yu(j1,j2))) + SUM(j2,0,2,SUM(j1,0,2,
         Conj(Yu(j1,j2))*Uu(2,j1))*Vu(2,j2))))));   }

   const double deltaRHat2LoopSMreal = std::real(deltaRHat2LoopSM);

   const double deltaRHat = deltaRHat1Loop + deltaRHat2LoopSMreal;

   return deltaRHat;
}

/**
 * Calculates the Standard Model vertex and box corrections
 * \f$\delta_{\text{VB}}^{\text{SM}}\f$ as given in Eq. (C.12) from
 * hep-ph/9606211 taking the tree-level value of the rhohat parameter
 * into account.
 *
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\delta_{\text{VB}}^{\text{SM}}\f$ as defined in (C.12)
 * from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_vb_sm(double sinThetaW) const
{
   const double mz  = sm_parameters.mz_pole;
   const double mw  = sm_parameters.mw_pole;
   const double cw2 = Sqr(mw / mz);
   const double sw2 = 1.0 - cw2;
   const double sinThetaW2 = Sqr(sinThetaW);
   const double outcos2    = 1.0 - sinThetaW2;
   const double q   = model->get_scale();

   const auto gY = model->get_g1() * standard_model_info::normalization_g1;
   const auto g2 = model->get_g2() * standard_model_info::normalization_g2;
   const double eDRbar     = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(eDRbar) / (4.0 * Pi);

   const double deltaVbSM = alphaDRbar / (4.0 * Pi * sinThetaW2) *
      (6.0 + log(cw2) / sw2 *
       (3.5 - 2.5 * sw2 - sinThetaW2 * (5.0 - 1.5 * cw2 / outcos2))
       - 4. * Log(Sqr(mz/q)));

   return deltaVbSM;
}

} // namespace weinberg_angle

} // namespace flexiblesusy
