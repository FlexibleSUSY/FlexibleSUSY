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

#include "gm2calc_interface.hpp"
#include "config.h"

/**
 * @file gm2calc_interface.cpp
 * @brief contains definitions of GM2Calc interface functions
 */

#ifdef ENABLE_GM2CALC

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2calc/THDM.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

namespace {

double calculate_e(double alpha) {
   return std::sqrt(4. * Pi * alpha);
}

double calculate_alpha(double e) {
   return e * e / (4. * Pi);
}

gm2calc::MSSMNoFV_onshell setup(const GM2Calc_MSSMNoFV_data& data)
{
   gm2calc::MSSMNoFV_onshell model;

   // fill couplings
   model.set_alpha_MZ(data.alpha_em_MZ);
   model.set_alpha_thompson(data.alpha_em_0);
   model.set_g3(calculate_e(data.alpha_s_MZ));

   // fill pole masses
   model.get_physical().MVZ     = data.MZ;    // 1L
   model.get_physical().MVWm    = data.MW;    // 1L
   model.get_physical().MFb     = data.mb_mb; // 1L
   model.get_physical().MFt     = data.MT;    // 1L
   model.get_physical().MFtau   = data.MTau;  // 1L
   model.get_physical().MFm     = data.MM;    // 1L
   model.get_physical().MAh(1)  = data.MA0;   // 2L
   model.get_physical().MSvmL   = data.MSvm;  // 1L
   model.get_physical().MSm     = data.MSm;   // 1L
   model.get_physical().MCha    = data.MCha;  // 1L
   model.get_physical().MChi    = data.MChi;  // 1L

   // fill DR-bar parameters
   model.set_scale(data.scale); // 2L
   model.set_TB(data.TB);       // 1L
   model.set_Mu(data.Mu);       // initial guess
   model.set_MassB(data.M1);    // initial guess
   model.set_MassWB(data.M2);   // initial guess
   model.set_MassG(data.M3);    // 2L
   model.set_mq2(data.mq2);     // 2L
   model.set_ml2(data.ml2);     // 1L
   model.set_md2(data.md2);     // 2L
   model.set_mu2(data.mu2);     // 2L
   model.set_me2(data.me2);     // 2L
   model.set_Au(data.Au);       // 2L
   model.set_Ad(data.Ad);       // 2L
   model.set_Ae(data.Ae);       // 2L

   // convert DR-bar parameters to on-shell
   model.convert_to_onshell();

   return model;
}

gm2calc::THDM setup(const GM2Calc_THDM_data& data)
{
   gm2calc::SM sm;
   sm.set_alpha_em_0(data.alpha_em_0);
   sm.set_alpha_em_mz(data.alpha_em_mz);
   sm.set_alpha_s_mz(data.alpha_s_mz);
   sm.set_mh(data.mh);
   sm.set_mw(data.mw);
   sm.set_mz(data.mz);
   sm.set_mu(data.mu);
   sm.set_md(data.md);
   sm.set_mv(data.mv);
   sm.set_ml(data.ml);
   sm.set_ckm(data.ckm);

   gm2calc::thdm::Gauge_basis basis;
   basis.yukawa_type = gm2calc::thdm::int_to_cpp_yukawa_type(data.yukawa_type);
   basis.lambda = data.lambda;
   basis.tan_beta = data.tan_beta;
   basis.m122 = data.m122;
   basis.zeta_u = data.zeta_u;
   basis.zeta_d = data.zeta_d;
   basis.zeta_l = data.zeta_l;
   basis.Delta_u = data.delta_u;
   basis.Delta_d = data.delta_d;
   basis.Delta_l = data.delta_l;
   basis.Pi_u = data.pi_u;
   basis.Pi_d = data.pi_d;
   basis.Pi_l = data.pi_l;

   return gm2calc::THDM(basis, sm);
}

} // anonymous namespace

/* ****************************** MSSM ****************************** */

/**
 * This function calculates \f$a_\mu\f$ up to the 2-loop level
 * including tan(beta) resummation in the MSSM with GM2Calc.
 *
 * If the calculation fails, 0 is returned.
 *
 * @param data data necessary for GM2Calc
 *
 * @return \f$a_\mu\f$ up to the 2-loop level w/ tan(beta) resummation
 */
double gm2calc_calculate_amu(const GM2Calc_MSSMNoFV_data& data)
{
   double amu = 0.;

   try {
      const gm2calc::MSSMNoFV_onshell model(setup(data));

      amu =
         + gm2calc::calculate_amu_1loop(model)
         + gm2calc::calculate_amu_2loop(model);
   } catch (const gm2calc::Error& e) {
      ERROR("GM2Calc: " << e.what());
   }

   return amu;
}

/**
 * This function calculates the uncertainty \f$\delta a_\mu\f$ of
 * \f$a_\mu\f$ at the 2-loop level in the MSSM by calling GM2Calc's
 * uncertainty estimation function.
 *
 * If the calculation fails, 0 is returned.
 *
 * @param data data necessary for GM2Calc
 *
 * @return \f$\delta a_\mu\f$ at the 2-loop level
 */
double gm2calc_calculate_amu_uncertainty(const GM2Calc_MSSMNoFV_data& data)
{
   double uncertainty = 0.;

   try {
      const gm2calc::MSSMNoFV_onshell model(setup(data));
      uncertainty = gm2calc::calculate_uncertainty_amu_2loop(model);
   } catch (const gm2calc::Error& e) {
      ERROR("GM2Calc: uncertainty estimation: " << e.what());
   }

   return uncertainty;
}

/* ****************************** THDM ****************************** */

/**
 * This function calculates \f$a_\mu\f$ up to the 2-loop level in the
 * THDM with GM2Calc.
 *
 * If the calculation fails, 0 is returned.
 *
 * @param data data necessary for GM2Calc
 *
 * @return \f$a_\mu\f$ up to the 2-loop level w/ tan(beta) resummation
 */
double gm2calc_calculate_amu(const GM2Calc_THDM_data& data)
{
   double amu = 0.;

   try {
      const gm2calc::THDM model = setup(data);

      amu =
         + gm2calc::calculate_amu_1loop(model)
         + gm2calc::calculate_amu_2loop(model);
   } catch (const gm2calc::Error& e) {
      ERROR("GM2Calc: " << e.what());
   }

   return amu;
}

/**
 * This function calculates the uncertainty \f$\delta a_\mu\f$ of
 * \f$a_\mu\f$ at the 2-loop level with GM2Calc's uncertainty
 * estimation function.
 *
 * If the calculation fails, 0 is returned.
 *
 * @param data data necessary for GM2Calc
 *
 * @return \f$\delta a_\mu\f$ at the 2-loop level
 */
double gm2calc_calculate_amu_uncertainty(const GM2Calc_THDM_data& data)
{
   double uncertainty = 0.;

   try {
      const gm2calc::THDM model = setup(data);
      uncertainty = gm2calc::calculate_uncertainty_amu_2loop(model);
   } catch (const gm2calc::Error& e) {
      ERROR("GM2Calc: uncertainty estimation: " << e.what());
   }

   return uncertainty;
}

} // namespace flexiblesusy

#else

namespace flexiblesusy {

/* ****************************** MSSM ****************************** */

double gm2calc_calculate_amu(const GM2Calc_MSSMNoFV_data&)
{
   return 0.;
}

double gm2calc_calculate_amu_uncertainty(const GM2Calc_MSSMNoFV_data&)
{
   return 0.;
}

/* ****************************** THDM ****************************** */

double gm2calc_calculate_amu(const GM2Calc_THDM_data&)
{
   return 0.;
}

double gm2calc_calculate_amu_uncertainty(const GM2Calc_THDM_data&)
{
   return 0.;
}

} // namespace flexiblesusy

#endif
