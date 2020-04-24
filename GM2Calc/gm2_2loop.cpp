// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#include "gm2_2loop_helpers.hpp"
#include "gm2_1loop_helpers.hpp"
#include "gm2_ffunctions.hpp"
#include "gm2_numerics.hpp"

#include <cmath>
#include <complex>

/**
 * \file gm2_2loop.cpp
 *
 * Contains functions necessary to calculate the SUSY contributions
 * for g-2 at the 2-loop level.
 */

namespace gm2calc {

namespace {

const double oneOver16PiSqr = 6.332573977646111e-3; // 1/(4 Pi)^2
const double root2 = 1.414213562373095; // Sqrt[2]

} // anonymous namespace

/**
 * \fn calculate_amu_2loop_non_tan_beta_resummed
 *
 * Calculates best 2-loop SUSY contribution to a_mu without tan(beta)
 * resummation.
 *
 * This function re-defines the muon Yukawa coupling in terms of the
 * tree-level relation with the muon pole mass, i.e. \f$y_\mu =
 * \frac{\sqrt{2} m_\mu^\text{pole}}{v_d}\f$.  Therefore, this
 * function does not use tan(beta) resummation.
 */
double calculate_amu_2loop_non_tan_beta_resummed(const MSSMNoFV_onshell& model)
{
   MSSMNoFV_onshell model_ytree(model);
   model_ytree.convert_to_non_tan_beta_resummed();

   return amu2LFSfapprox_non_tan_beta_resummed(model_ytree)
      + amu2LChipmPhotonic(model_ytree)
      + amu2LChi0Photonic(model_ytree)
      + amu2LaSferm(model_ytree)
      + amu2LaCha(model_ytree);
}

/**
 * \fn calculate_amu_2loop
 *
 * Calculates best 2-loop SUSY contribution to a_mu with tan(beta)
 * resummation.
 */
double calculate_amu_2loop(const MSSMNoFV_onshell& model)
{
   return amu2LFSfapprox(model)
      + amu2LChipmPhotonic(model)
      + amu2LChi0Photonic(model)
      + amu2LaSferm(model)
      + amu2LaCha(model);
}

// fermion/sfermion corrections, log-approximations

/**
 * Calculates \f$m_{SUSY}\f$, p.37 arxiv:1311.1775.
 * Finds minimum of special masses to normalize logarithms.
 */
double log_scale(const MSSMNoFV_onshell& model)
{
   return std::fmin(std::abs(model.get_MassB()),
           std::fmin(std::abs(model.get_MassWB()),
            std::fmin(std::abs(model.get_Mu()),
             std::fmin(std::sqrt(model.get_me2(1, 1)),
              std::sqrt(model.get_ml2(1, 1))))));
}

/**
 * Calculates \f$\Delta_{g_1}\f$, Eq (6.6a) arxiv:1311.1775.
 *
 * Contributions from 1st and 2nd generation sleptons have been
 * included in addition.
 */
double delta_g1(const MSSMNoFV_onshell& model)
{
   const double gY = model.get_gY();
   const Eigen::Matrix<double,3,3>& mu2(model.get_mu2());
   const Eigen::Matrix<double,3,3>& md2(model.get_md2());
   const Eigen::Matrix<double,3,3>& mq2(model.get_mq2());
   const Eigen::Matrix<double,3,3>& me2(model.get_me2());
   const Eigen::Matrix<double,3,3>& ml2(model.get_ml2());
   const double logscale = log_scale(model);

   return sqr(gY) * oneOver16PiSqr * 4. / 3. *
          (4. / 3. * std::log(std::sqrt(mu2(0, 0)) / logscale) +
           4. / 3. * std::log(std::sqrt(mu2(1, 1)) / logscale) +
           4. / 3. * std::log(std::sqrt(mu2(2, 2)) / logscale) +
           1. / 3. * std::log(std::sqrt(md2(0, 0)) / logscale) +
           1. / 3. * std::log(std::sqrt(md2(1, 1)) / logscale) +
           1. / 3. * std::log(std::sqrt(md2(2, 2)) / logscale) +
           1. / 6. * std::log(std::sqrt(mq2(0, 0)) / logscale) +
           1. / 6. * std::log(std::sqrt(mq2(1, 1)) / logscale) +
           1. / 6. * std::log(std::sqrt(mq2(2, 2)) / logscale) +
           std::log(std::sqrt(me2(0, 0)) / logscale) +
           std::log(std::sqrt(me2(1, 1)) / logscale) +
           std::log(std::sqrt(me2(2, 2)) / logscale) +
           0.5 * std::log(std::sqrt(ml2(0, 0)) / logscale) +
           0.5 * std::log(std::sqrt(ml2(1, 1)) / logscale) +
           0.5 * std::log(std::sqrt(ml2(2, 2)) / logscale));
}

/**
 * Calculates \f$\Delta_{\tilde{H}}\f$, Eq (6.6c) arxiv:1311.1775.
 */
double delta_yuk_higgsino(const MSSMNoFV_onshell& model)
{
   const double ytau = model.get_Ye(2, 2);
   const double ytop = model.get_Yu(2, 2);
   const double ybot = model.get_Yd(2, 2);
   const Eigen::Matrix<double,3,3>& mu2(model.get_mu2());
   const Eigen::Matrix<double,3,3>& md2(model.get_md2());
   const Eigen::Matrix<double,3,3>& mq2(model.get_mq2());
   const Eigen::Matrix<double,3,3>& me2(model.get_me2());
   const Eigen::Matrix<double,3,3>& ml2(model.get_ml2());
   const double logscale = log_scale(model);

   return oneOver16PiSqr * 0.5 *
          (3. * sqr(ytop) * std::log(std::sqrt(mu2(2, 2)) / logscale) +
           3. * sqr(ybot) * std::log(std::sqrt(md2(2, 2)) / logscale) +
           3. * (sqr(ytop) + sqr(ybot)) *
              std::log(std::sqrt(mq2(2, 2)) / logscale) +
           sqr(ytau) * (std::log(std::sqrt(me2(2, 2)) / logscale) +
                        std::log(std::sqrt(ml2(2, 2)) / logscale)));
}

/**
 * Calculates \f$\Delta_{\tilde{B}\tilde{H}}\f$, Eq (6.6d)
 * arxiv:1311.1775.
 */
double delta_yuk_bino_higgsino(const MSSMNoFV_onshell& model)
{
   const double ytop = model.get_Yu(2, 2);
   const Eigen::Matrix<double,3,3>& mu2(model.get_mu2());
   const Eigen::Matrix<double,3,3>& mq2(model.get_mq2());
   const double logscale = log_scale(model);

   return oneOver16PiSqr * sqr(ytop) *
          (-8. * std::log(std::sqrt(mu2(2, 2)) / logscale) +
           2. * std::log(std::sqrt(mq2(2, 2)) / logscale));
}

/**
 * Calculates \f$\Delta_{g_2}\f$, Eq (6.6b) arxiv:1311.1775.
 *
 * Contributions from 1st and 2nd generation sleptons have been
 * included in addition.
 */
double delta_g2(const MSSMNoFV_onshell& model)
{
   const double g2 = model.get_g2();
   const Eigen::Matrix<double,3,3>& mq2(model.get_mq2());
   const Eigen::Matrix<double,3,3>& ml2(model.get_ml2());
   const double logscale = log_scale(model);

   return sqr(g2) * oneOver16PiSqr * 4. / 3. *
          (1.5 * std::log(std::sqrt(mq2(0, 0)) / logscale) +
           1.5 * std::log(std::sqrt(mq2(1, 1)) / logscale) +
           1.5 * std::log(std::sqrt(mq2(2, 2)) / logscale) +
           0.5 * std::log(std::sqrt(ml2(0, 0)) / logscale) +
           0.5 * std::log(std::sqrt(ml2(1, 1)) / logscale) +
           0.5 * std::log(std::sqrt(ml2(2, 2)) / logscale));
}

/**
 * Calculates \f$\Delta_{\tilde{W}\tilde{H}}\f$, Eq (6.6e)
 * arxiv:1311.1775.
 */
double delta_yuk_wino_higgsino(const MSSMNoFV_onshell& model)
{
   const double ytop = model.get_Yu(2, 2);
   const Eigen::Matrix<double,3,3>& mq2(model.get_mq2());
   const double logscale = log_scale(model);

   return oneOver16PiSqr * (-6.) * sqr(ytop) *
          std::log(std::sqrt(mq2(2, 2)) / logscale);
}

/**
 * Calculates \f$\Delta_{t_\beta}\f$, Eq (6.6f) arxiv:1311.1775.
 */
double delta_tan_beta(const MSSMNoFV_onshell& model)
{
   const double ytau = model.get_Ye(2, 2);
   const double ytop = model.get_Yu(2, 2);
   const double ybot = model.get_Yd(2, 2);
   const double logscale = log_scale(model);
   const double Q = model.get_scale();

   return oneOver16PiSqr * (sqr(ytau) - 3. * sqr(ytop) + 3. * sqr(ybot)) *
          std::log(Q / logscale);
}

/**
 * Calculates 1st line of Eq (6.5) arxiv:1311.1775.
 */
double amu2LWHnu(const MSSMNoFV_onshell& model)
{
   return amu1LWHnu(model) *
          (0.015 + delta_g2(model) + delta_yuk_higgsino(model) +
           delta_yuk_wino_higgsino(model) + delta_tan_beta(model));
}

/**
 * Calculates 2nd line of Eq (6.5) arxiv:1311.1775.
 */
double amu2LWHmuL(const MSSMNoFV_onshell& model)
{
   return amu1LWHmuL(model) *
          (0.015 + delta_g2(model) + delta_yuk_higgsino(model) +
           delta_yuk_wino_higgsino(model) + delta_tan_beta(model));
}

/**
 * Calculates 3rd line of Eq (6.5) arxiv:1311.1775.
 */
double amu2LBHmuL(const MSSMNoFV_onshell& model)
{
   return amu1LBHmuL(model) *
          (0.015 + delta_g1(model) + delta_yuk_higgsino(model) +
           delta_yuk_bino_higgsino(model) + delta_tan_beta(model));
}

/**
 * Calculates 4th line of Eq (6.5) arxiv:1311.1775.
 */
double amu2LBHmuR(const MSSMNoFV_onshell& model)
{
   return amu1LBHmuR(model) *
          (0.04 + delta_g1(model) + delta_yuk_higgsino(model) +
           delta_yuk_bino_higgsino(model) + delta_tan_beta(model));
}

/**
 * Calculates 5th line of Eq (6.5) arxiv:1311.1775.
 */
double amu2LBmuLmuR(const MSSMNoFV_onshell& model)
{
   return amu1LBmuLmuR(model) *
          (0.03 + delta_g1(model) + delta_tan_beta(model));
}

/**
 * Calculates 2-loop leading log approximation for fermion-sfermion
 * loop contributions, Eq (6.5) arxiv:1311.1775.
 *
 * No tan(beta) resummation
 */
double amu2LFSfapprox_non_tan_beta_resummed(const MSSMNoFV_onshell& model)
{
   return amu2LWHnu(model) + amu2LWHmuL(model) + amu2LBHmuL(model) +
          amu2LBHmuR(model) + amu2LBmuLmuR(model);
}

/**
 * Calculates 2-loop leading log approximation for fermion-sfermion
 * loop contributions, Eq (6.5) arxiv:1311.1775.
 *
 * Includes tan(beta) resummation
 */
double amu2LFSfapprox(const MSSMNoFV_onshell& model)
{

   return amu2LFSfapprox_non_tan_beta_resummed(model) * tan_beta_cor(model);
}

// === photonic 2-loop corrections ===

/**
 * Calculates the photonic 2-loop contribution to the 1-loop chargino
 * diagram, Eq (35) arXiv:1003.5820.
 */
double amu2LChipmPhotonic(const MSSMNoFV_onshell& model)
{
   const double MM = model.get_MM();
   const Eigen::Array<double,2,1> AAC_(AAC(model));
   const Eigen::Array<double,2,1> BBC_(BBC(model));
   const Eigen::Array<double,2,1>& MCha(model.get_MCha());
   const double MSvmL = model.get_MSvmL();
   const Eigen::Array<double,2,1> x__k(x_k(model));
   const double Q = model.get_scale();

   double result = 0.;

   for (int k = 0; k < 2; k++) {
      result += (AAC_(k) * F1C(x__k(k)) / 12.
                 + MCha(k) / 3. * BBC_(k) / model.get_MM() * F2C(x__k(k)))
                    * 16. * std::log(MM / MSvmL)
                  - 47. * AAC_(k) * F3C(x__k(k)) / 72.
                  - 61. * MCha(k) / 9. * BBC_(k) / model.get_MM() * F4C(x__k(k))
                  - (0.5 * AAC_(k) * F1C(x__k(k))
                     + MCha(k) * BBC_(k) / model.get_MM() * F2C(x__k(k)))
                    * std::log(sqr(MSvmL / Q));
   }

   return  sqr(model.get_EL0()) * sqr(oneOver16PiSqr) * sqr(MM / MSvmL) * result;
}

/**
 * Calculates the photonic 2-loop contribution to the 1-loop
 * neutralino diagram, Eq (36) arXiv:1003.5820.
 */
double amu2LChi0Photonic(const MSSMNoFV_onshell& model)
{
   const double MM = model.get_MM();
   const Eigen::Matrix<double,4,2> AAN_(AAN(model));
   const Eigen::Matrix<double,4,2> BBN_(BBN(model));
   const Eigen::Array<double,4,1>& MNeu(model.get_MChi());
   const Eigen::Array<double,2,1>& MSmu(model.get_MSm());
   const Eigen::Matrix<double,4,2> x__im(x_im(model));
   const double Q = model.get_scale();

   double result = 0.;

   for (int i = 0; i < 4; ++i) {
      for (int m = 0; m < 2; ++m) {
         result += 1. / sqr(MSmu(m))
                    * ((- 1. / 12. * AAN_(i, m) * F1N(x__im(i, m))
                       - MNeu(i) / 6. * BBN_(i, m) / model.get_MM() * F2N(x__im(i, m)))
                        * 16. * std::log(MM / MSmu(m))
                      + 35. / 72. * AAN_(i, m) * F3N(x__im(i, m))
                      + 8. * MNeu(i) / 9. * BBN_(i, m) / model.get_MM() * F4N(x__im(i, m))
                      + (0.25 * AAN_(i, m) * F1N(x__im(i, m)))
                       * std::log(sqr(MSmu(m) / Q)) );
      }
   }

   return sqr(model.get_EL0()) * sqr(oneOver16PiSqr) * sqr(MM) * result;
}

// ==== amu2Loop_a corrections ====

/**
 * The following functions include resummation of 1/(1 + Delta_mu) within
 * the muon, tau and bottom Yukawa couplings.
 */

/**
 * Calculates CP-even Higgs mixing angle \f$\tan\alpha\f$.
 *
 * @note the result is < 0
 */
double tan_alpha(const MSSMNoFV_onshell& model)
{
   const double TB = model.get_TB();
   const double MZ = model.get_MZ();
   const double MA0 = model.get_MA0();
   const double tan2beta = 2. * TB / (1. - sqr(TB));
   const double tan2alpha = tan2beta * (sqr(MA0) + sqr(MZ)) / (sqr(MA0) - sqr(MZ));

   return -1.0 / tan2alpha - std::sqrt(1.0 / sqr(tan2alpha) + 1.0);
}

/**
 * Calculates \f$\lambda_{\mu}\f$, Eq (65), and
 * \f$\lambda_{\chi_k^+}\f$, Eq (66) arXiv:hep-ph/0609168.
 *
 * Row 0 contains \f$\lambda_{\chi_1^+}\f$, Eq. (66).
 * Row 1 contains \f$\lambda_{\chi_2^+}\f$, Eq. (66).
 * Row 2 contains \f$\lambda_{\mu}\f$, Eq. (65).
 *
 * includes tan(beta) resummation
 */
Eigen::Matrix<std::complex<double>,3,3> lambda_mu_cha(const MSSMNoFV_onshell& model)
{
   const Eigen::Array<double,2,1>& MCha(model.get_MCha());
   const double MW(model.get_MW());
   const double TB(model.get_TB());
   const double CB = 1. / std::sqrt(1. + sqr(TB));
   const double SB = std::sqrt(1. - sqr(CB));
   const double TA(tan_alpha(model));
   const double CA = 1. / std::sqrt(1. + sqr(TA));
   const double SA = - std::sqrt(1. - sqr(CA));
   const Eigen::Matrix<std::complex<double>,2,2>& U(model.get_UM());
   const Eigen::Matrix<std::complex<double>,2,2>& V(model.get_UP());
   const double one_over_cb_eff = root2 * model.get_Ye(1,1)
      * model.get_MW() / model.get_MM() / model.get_g2();

   Eigen::Matrix<std::complex<double>,3,3> result;

   for (int k = 0; k < 2; ++k) {
      result(k, 0) = root2 * MW / MCha(k)
                      * (U(k, 0) * V(k, 1) * CA + U(k, 1) * V(k, 0) * (-SA));
      result(k, 1) = root2 * MW / MCha(k)
                      * (U(k, 0) * V(k, 1) * SA + U(k, 1) * V(k, 0) * CA);
      result(k, 2) = root2 * MW / MCha(k)
                      * (U(k, 0) * V(k, 1) * (-CB) + U(k, 1) * V(k, 0) * (-SB));
   }
   result(2, 0) = -SA * one_over_cb_eff;
   result(2, 1) = CA * one_over_cb_eff;
   result(2, 2) = SB * one_over_cb_eff;

   return result;
}

/**
 * Calculates \f$\lambda_{\tilde{b}_i}\f$, Eq (67)
 * arXiv:hep-ph/0609168
 */
Eigen::Matrix<std::complex<double>,2,2> lambda_stop(const MSSMNoFV_onshell& model) {
   const double TB(model.get_TB());
   const double CB = 1. / std::sqrt(1. + sqr(TB));
   const double SB = std::sqrt(1. - sqr(CB));
   const double TA(tan_alpha(model));
   const double CA = 1. / std::sqrt(1. + sqr(TA));
   const double SA = - std::sqrt(1. - sqr(CA));
   const double MT(model.get_MT());
   const Eigen::Array<double,2,1>& MStop(model.get_MSt());
   const Eigen::Matrix<double,2,2>& UStop(model.get_USt());
   const double At(model.get_Au(2, 2));
   const double Mu(model.get_Mu());

   Eigen::Matrix<std::complex<double>,2,2> result;

   for (int i = 0; i < 2; ++i) {
      result(i, 0) = 2. * MT / (sqr(MStop(i)) * SB) * (Mu * SA + At * CA)
                      * std::conj(UStop(i, 0)) * UStop(i, 1);
      result(i, 1) = 2. * MT / (sqr(MStop(i)) * SB) * (Mu * (-CA) + At * SA)
                      * std::conj(UStop(i, 0)) * UStop(i, 1);
   }

   return result;
}

/**
 * Calculates \f$\lambda_{\tilde{b}_i}\f$, Eq (68)
 * arXiv:hep-ph/0609168
 *
 * includes tan(beta) resummation
 */
Eigen::Matrix<std::complex<double>,2,2> lambda_sbot(const MSSMNoFV_onshell& model)
{
   const double TA(tan_alpha(model));
   const double CA = 1. / std::sqrt(1. + sqr(TA));
   const double SA = - std::sqrt(1. - sqr(CA));
   const Eigen::Array<double,2,1>& MSbot(model.get_MSb());
   const Eigen::Matrix<double,2,2>& USbot(model.get_USb());
   const double Ab(model.get_Ad(2, 2));
   const double Mu(model.get_Mu());
   const double mb_over_cb_eff = root2 * model.get_Yd(2,2)
      * model.get_MW() / model.get_g2();

   Eigen::Matrix<std::complex<double>,2,2> result;

   for (int i = 0; i < 2; ++i) {
      result(i, 0) = 2. * mb_over_cb_eff / (sqr(MSbot(i))) * (- Mu * CA + Ab * (-SA))
                      * std::conj(USbot(i, 0)) * USbot(i, 1);
      result(i, 1) = 2. * mb_over_cb_eff / (sqr(MSbot(i))) * (- Mu * SA + Ab * CA)
                      * std::conj(USbot(i, 0)) * USbot(i, 1);
   }

   return result;
}

/**
 * Calculates \f$\lambda_{\tilde{t}_i}\f$, Eq (69)
 * arXiv:hep-ph/0609168
 *
 * includes tan(beta) resummation
 */
Eigen::Matrix<std::complex<double>,2,2> lambda_stau(const MSSMNoFV_onshell& model)
{
   const double TA(tan_alpha(model));
   const double CA = 1. / std::sqrt(1. + sqr(TA));
   const double SA = - std::sqrt(1. - sqr(CA));
   const Eigen::Array<double,2,1>& MStau(model.get_MStau());
   const Eigen::Matrix<double,2,2>& UStau(model.get_UStau());
   const double Al(model.get_Ae(2, 2));
   const double Mu(model.get_Mu());
   const double mtau_over_cb_eff = root2 * model.get_Ye()(2,2)
      * model.get_MW() / model.get_g2();

   Eigen::Matrix<std::complex<double>,2,2> result;

   for (int i = 0; i < 2; ++i) {
      result(i, 0) = 2. * mtau_over_cb_eff / (sqr(MStau(i))) * (- Mu * CA + Al * (-SA))
                      * std::conj(UStau(i, 0)) * UStau(i, 1);
      result(i, 1) = 2. * mtau_over_cb_eff / (sqr(MStau(i))) * (- Mu * SA + Al * CA)
                      * std::conj(UStau(i, 0)) * UStau(i, 1);
   }

   return result;
}

/**
 * Calculates 2-loop contribution to amu, where a sfermion loop has
 * been inserted into a 1-loop Standard Model diagram (photonic
 * Barr-Zee diagram \f$(\tilde{f}\gamma H)\f$), Eq (64)
 * arXiv:hep-ph/0609168.
 */
double amu2LaSferm(const MSSMNoFV_onshell& model)
{
   const double MM(model.get_MM());
   const double MW(model.get_MW());
   const double SW = std::sqrt(1. - sqr(MW / model.get_MZ()));
   const double EL(model.get_EL());
   const Eigen::Array<double,2,1>& m_stop(model.get_MSt());
   const Eigen::Array<double,2,1>& m_sbot(model.get_MSb());
   const Eigen::Array<double,2,1>& m_stau(model.get_MStau());
   const Eigen::Array<double,2,1>& M_higgs(model.get_Mhh());
   const Eigen::Matrix<std::complex<double>,3,3> lambda(lambda_mu_cha(model));
   Eigen::Array<std::complex<double>,2,1> lambda_mu;
   lambda_mu(0) = lambda(2, 0);
   lambda_mu(1) = lambda(2, 1);
   const Eigen::Matrix<std::complex<double>,2,2> lambdastop(lambda_stop(model));
   const Eigen::Matrix<std::complex<double>,2,2> lambdasbot(lambda_sbot(model));
   const Eigen::Matrix<std::complex<double>,2,2> lambdastau(lambda_stau(model));

   double result = 0.;
   double N_c = 3.;
   double Q = 2./3.;

   for (int i = 0; i < 2; ++i) {
      for (int s = 0; s < 2; ++s) {
         result += N_c * sqr(Q) * real(lambda_mu(s)
                   * lambdastop(i, s))
                   * f_sferm(sqr(m_stop(i) / M_higgs(s)));
      }
   }

   Q = -1./3.;
   for (int i = 0; i < 2; ++i) {
      for (int s = 0; s < 2; ++s) {
         result += N_c * sqr(Q) * real(lambda_mu(s)
                   * lambdasbot(i, s))
                   * f_sferm(sqr(m_sbot(i) / M_higgs(s)));
      }
   }

   N_c = 1.;
   Q = -1.;
   for (int i = 0; i < 2; ++i) {
      for (int s = 0; s < 2; ++s) {
         result += N_c * sqr(Q) * real(lambda_mu(s)
                   * lambdastau(i, s))
                   * f_sferm(sqr(m_stau(i) / M_higgs(s)));
      }
   }

   return result * 2. * sqr(oneOver16PiSqr * sqr(EL) * MM / (MW * SW));
}

/**
 * Calculates 2-loop contribution to amu, where a chargino loop has
 * been inserted into a 1-loop Standard Model diagram (photonic
 * Barr-Zee diagram \f$(\chi\gamma H)\f$), Eq (63)
 * arXiv:hep-ph/0609168.
 */
double amu2LaCha(const MSSMNoFV_onshell& model)
{
   const double MM(model.get_MM());
   const double MW(model.get_MW());
   const double MA0(model.get_MA0());
   const Eigen::Array<double,2,1>& M_higgs(model.get_Mhh());
   const double SW = std::sqrt(1. - sqr(MW / model.get_MZ()));
   const double EL(model.get_EL());
   const Eigen::Array<double,2,1>& m_cha(model.get_MCha());
   const Eigen::Matrix<std::complex<double>,3,3> lambda(lambda_mu_cha(model));

   double result = 0.;

   for (int k = 0; k < 2; ++k) {
      result += real(lambda(2, 2) * lambda(k, 2)) * f_PS(sqr(m_cha(k) / MA0));
      for (int s = 0; s < 2; ++s) {
         result += real(lambda(2, s) * lambda(k, s))
                   * f_S(sqr(m_cha(k) / M_higgs(s)));
      }
   }

   return result * 2. * sqr(oneOver16PiSqr * sqr(EL) * MM / (MW * SW));
}

} // namespace gm2calc
