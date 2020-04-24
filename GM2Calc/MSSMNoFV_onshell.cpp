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

#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2calc/gm2_error.hpp"

#include "gm2_1loop_helpers.hpp"
#include "gm2_constants.h"
#include "gm2_log.hpp"
#include "gm2_mb.hpp"
#include "gm2_numerics.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <boost/math/tools/roots.hpp>

namespace gm2calc {

namespace {
   const double Pi = 3.141592653589793;
   const double eps = std::numeric_limits<double>::epsilon();
   const double root2 = 1.414213562373095; // Sqrt[2]

   /// fine-structure constant in the Thompson limit (Q = 0) from PDG (2019)
   const double ALPHA_EM_THOMPSON = GM2CALC_ALPHA_EM_THOMPSON;

   /// fine-structure constant at Q = MZ
   const double ALPHA_EM_MZ = GM2CALC_ALPHA_EM_MZ;

   /// calculates gauge coupling from alpha
   double calculate_e(double alpha) {
      return std::sqrt(4. * Pi * alpha);
   }

   /// calculates alpha from gauge coupling
   double calculate_alpha(double e) {
      return e * e / (4. * Pi);
   }

   /// prints Eigen Matrix/Array in pretty format
   template <class Derived>
   std::string pretty_print(const Eigen::DenseBase<Derived>& v)
   {
      const bool is_matrix = v.rows() > 1;

      std::ostringstream sstr;

      for (auto i = 0; i < v.rows(); i++) {
         if (is_matrix && i == 0) {
            sstr << '{';
         }

         for (auto k = 0; k < v.cols(); k++) {
            if (k == 0) {
               sstr << '{';
            }
            sstr << v(i,k);
            if (k + 1 == v.cols()) {
               sstr << '}';
            } else {
               sstr << ", ";
            }
         }

         if (is_matrix) {
            if (i + 1 == v.rows()) {
               sstr << '}';
            } else {
               sstr << ", ";
            }
         }
      }

      return sstr.str();
   }

} // anonymous namespace

namespace detail {

template <class Derived>
bool is_equal(const Eigen::ArrayBase<Derived>& a,
              const Eigen::ArrayBase<Derived>& b,
              double precision_goal)
{
   return (a - b).cwiseAbs().maxCoeff() < precision_goal;
}

template <class Derived>
bool is_zero(const Eigen::ArrayBase<Derived>& a, double eps)
{
   return a.cwiseAbs().maxCoeff() < eps;
}

bool is_zero(double a, double eps)
{
   return gm2calc::is_zero(a, eps);
}

/**
 * Returns index of most bino-like neutralino.  The function extracts
 * this information from the given neutralino mixing matrix.
 *
 * @param ZN neutralino mixing matrix
 */
template <class Derived>
unsigned find_bino_like_neutralino(const Eigen::MatrixBase<Derived>& ZN)
{
   unsigned max_bino;
   ZN.col(0).cwiseAbs().maxCoeff(&max_bino);

   return max_bino;
}

/**
 * Returns index of most right-handed smuon.  The function extracts
 * this information from the given smuon mixing matrix.
 *
 * @param ZM smuon mixing matrix
 */
template <class Derived>
unsigned find_right_like_smuon(const Eigen::MatrixBase<Derived>& ZM)
{
   return (std::abs(ZM(0,0)) > std::abs(ZM(0,1))) ? 1 : 0;
}

} // namespace detail

MSSMNoFV_onshell::MSSMNoFV_onshell()
   : EL(calculate_e(ALPHA_EM_MZ))
   , EL0(calculate_e(ALPHA_EM_THOMPSON))
{
   set_g3(calculate_e(GM2CALC_ALPHA_S_MZ));
   get_physical().MFt = GM2CALC_MT;
   get_physical().MFb = GM2CALC_MBMB;
   get_physical().MFe = GM2CALC_ME;
   get_physical().MFm  = GM2CALC_MM;
   get_physical().MFtau = GM2CALC_ML;
   get_physical().MVWm = GM2CALC_MW;
   get_physical().MVZ  = GM2CALC_MZ;
   set_scale(get_physical().MVZ);
}

void MSSMNoFV_onshell::set_alpha_MZ(double alpha)
{
   EL = calculate_e(alpha);
}

void MSSMNoFV_onshell::set_alpha_thompson(double alpha)
{
   EL0 = calculate_e(alpha);
}

void MSSMNoFV_onshell::set_TB(double tanb)
{
   const double vev = get_vev();
   const double rt = std::sqrt(1 + tanb*tanb);
   const double sinb = tanb / rt;
   const double cosb = 1.   / rt;
   set_vd(vev * cosb);
   set_vu(vev * sinb);
}

double MSSMNoFV_onshell::get_vev() const
{
   const double cW = get_MW() / get_MZ();
   const double g2 = get_EL() / std::sqrt(1. - cW*cW);
   const double vev = 2. * get_MW() / g2;

   return vev;
}

double MSSMNoFV_onshell::get_TB() const
{
   if (gm2calc::detail::is_zero(get_vd(), eps)) {
      throw EInvalidInput("down-type VEV vd = 0");
   }
   return get_vu() / get_vd();
}

/**
 * Converts the model parameters from the DR-bar scheme to the mixed
 * on-shell/DR-bar scheme.
 *
 * The function assumes that the physical struct is filled with SM and
 * SUSY pole masses and corresponding mixing matrices.  From these
 * quantities, the model parameters are calculated in the mixed
 * on-shell/DR-bar scheme, used to calculate amu.
 *
 * This function is intended to be used with input parameters in the
 * (SLHA-compatible) DR-bar scheme.
 *
 * @param precision accuracy goal for the conversion
 * @param max_iterations maximum number of iterations
 */
void MSSMNoFV_onshell::convert_to_onshell(
   double precision, unsigned max_iterations)
{
   check_input();
   calculate_DRbar_masses();
   copy_susy_masses_to_pole();

   calculate_mb_DRbar_MZ();
   convert_gauge_couplings();
   convert_BMu();
   convert_vev();
   convert_yukawa_couplings_treelevel();
   convert_Mu_M1_M2(precision, max_iterations);
   convert_yukawa_couplings(); // first guess of resummed yukawas
   convert_ml2();
   convert_yukawa_couplings();
   calculate_MSm(); // calculate MSm with new value of ml2(1,1) and ym
   convert_me2(precision, max_iterations);
   convert_yukawa_couplings();

   // final mass spectrum
   get_problems().clear_problems();
   calculate_DRbar_masses();
   check_problems();
}

/**
 * Calculates the SUSY mass spectrum in the mixed on-shell/DR-bar scheme.
 *
 * The function assumes that the physical struct is filled with SM
 * pole masses and the SUSY parameters are given in GM2Calc-specific
 * mixed on-shell/DR-bar scheme.  From these input quantities, the
 * corresponding (mixed on-shell/DR-bar) SUSY mass spectrum is
 * calculated.
 *
 * This function is inteded to be used with the input parameters given
 * in the GM2Calc-specific on-shell/DR-bar scheme.
 */
void MSSMNoFV_onshell::calculate_masses()
{
   check_input();
   calculate_mb_DRbar_MZ();
   convert_gauge_couplings();
   convert_BMu();
   convert_vev();
   convert_yukawa_couplings_treelevel();
   convert_yukawa_couplings(); // tan(beta) resummation in Yukawas

   // final mass spectrum
   get_problems().clear();
   calculate_DRbar_masses();
   check_problems();

   // copy SUSY masses to physical struct
   copy_susy_masses_to_pole();
   get_physical().MAh = get_MAh();
   get_physical().Mhh = get_Mhh();
}

/**
 * Converts parameters and masses to the case of no tan(beta)
 * resummation.
 */
void MSSMNoFV_onshell::convert_to_non_tan_beta_resummed()
{
   convert_yukawa_couplings_treelevel();
   calculate_DRbar_masses();
   check_problems();
}

void MSSMNoFV_onshell::check_input() const
{
#define WARN_OR_THROW_IF(cond,msg)                      \
   if (cond) {                                          \
      if (do_force_output()) {                          \
         WARNING(msg);                                  \
      } else {                                          \
         throw EInvalidInput(msg);                      \
      }                                                 \
   }

#define WARN_OR_THROW_IF_ZERO(mass,msg)                         \
   WARN_OR_THROW_IF(gm2calc::detail::is_zero(get_##mass(), eps), msg)

   const double MW = get_MW();
   const double MZ = get_MZ();
   const double TB = get_TB();

   WARN_OR_THROW_IF(MW >= MZ   , "MW >= MZ cannot be treated with GM2Calc");
   WARN_OR_THROW_IF_ZERO(MW    , "W mass is zero");
   WARN_OR_THROW_IF_ZERO(MZ    , "Z mass is zero");
   WARN_OR_THROW_IF_ZERO(MM    , "Muon mass is zero");
   WARN_OR_THROW_IF_ZERO(Mu    , "mu parameter is zero");
   WARN_OR_THROW_IF_ZERO(MassB , "Bino mass M1 is zero");
   WARN_OR_THROW_IF_ZERO(MassWB, "Wino mass M2 is zero");
   WARN_OR_THROW_IF_ZERO(TB    , "tan(beta) is zero");
   WARN_OR_THROW_IF(!std::isfinite(TB), "tan(beta) is infinite");

#undef WARN_OR_THROW_IF
#undef WARN_OR_THROW_IF_ZERO
}

void MSSMNoFV_onshell::check_problems() const
{
   if (get_problems().have_problem()) {
      if (!do_force_output()) {
         throw EPhysicalProblem(get_problems().get_problems());
      }
   }
   if (get_mu2().diagonal().minCoeff() < 0. ||
       get_md2().diagonal().minCoeff() < 0. ||
       get_mq2().diagonal().minCoeff() < 0. ||
       get_me2().diagonal().minCoeff() < 0. ||
       get_ml2().diagonal().minCoeff() < 0.) {
      if (!do_force_output()) {
         throw EInvalidInput("soft mass squared < 0");
      }
   }
   if (gm2calc::detail::is_zero(get_MCha(0), eps)) {
      if (!do_force_output()) {
         throw EInvalidInput("lightest chargino mass = 0");
      }
   }
}

void MSSMNoFV_onshell::copy_susy_masses_to_pole()
{
   // if pole masses are zero, interpret the tree-level masses as pole
   // masses

#define COPY_IF_ZERO_0(m)                                               \
   if (gm2calc::detail::is_zero(get_physical().m, eps)) {               \
      get_physical().m = get_##m();                                     \
   }
#define COPY_IF_ZERO_1(m,z)                                             \
   if (gm2calc::detail::is_zero(get_physical().m, eps)) {               \
      get_physical().m = get_##m();                                     \
      get_physical().z = get_##z();                                     \
   }
#define COPY_IF_ZERO_2(m,u,v)                                           \
   if (gm2calc::detail::is_zero(get_physical().m, eps)) {               \
      get_physical().m = get_##m();                                     \
      get_physical().u = get_##u();                                     \
      get_physical().v = get_##v();                                     \
   }

   COPY_IF_ZERO_1(MChi, ZN);
   COPY_IF_ZERO_2(MCha, UM, UP);
   COPY_IF_ZERO_0(MSveL);
   COPY_IF_ZERO_0(MSvmL);
   COPY_IF_ZERO_0(MSvtL);
   COPY_IF_ZERO_1(MSd, ZD);
   COPY_IF_ZERO_1(MSu, ZU);
   COPY_IF_ZERO_1(MSe, ZE);
   COPY_IF_ZERO_1(MSm, ZM);
   COPY_IF_ZERO_1(MStau, ZTau);
   COPY_IF_ZERO_1(MSs, ZS);
   COPY_IF_ZERO_1(MSc, ZC);
   COPY_IF_ZERO_1(MSb, ZB);
   COPY_IF_ZERO_1(MSt, ZT);

#undef COPY_IF_ZERO_0
#undef COPY_IF_ZERO_1
#undef COPY_IF_ZERO_2
}

void MSSMNoFV_onshell::convert_gauge_couplings()
{
   const double MW = get_MW(); // pole mass
   const double MZ = get_MZ(); // pole mass
   const double cW = MW / MZ;  // on-shell weak mixing angle

   set_g1(std::sqrt(5. / 3.) * EL / cW);
   set_g2(EL / std::sqrt(1. - sqr(cW)));
}

void MSSMNoFV_onshell::convert_BMu()
{
   const double TB = get_TB(); // DR-bar
   const double MA = get_MA0(); // pole mass
   const double sin2b = 2. * TB / (1. + sqr(TB));

   set_BMu(0.5 * sqr(MA) * sin2b);
}

void MSSMNoFV_onshell::convert_vev()
{
   const double TB = get_TB(); // DR-bar
   const double MW = get_MW(); // pole mass
   const double vev = 2. * MW / get_g2();

   set_vu(vev / std::sqrt(1. + 1. / sqr(TB)));
   set_vd(get_vu() / TB);
}

/**
 * Calculates mb(MZ) in DR-bar scheme.
 *
 * mb(MZ) DR-bar is calculated from mb(mb) MS-bar using the function
 * \a calculate_mb_SM5_DRbar .
 */
void MSSMNoFV_onshell::calculate_mb_DRbar_MZ()
{
   const double mb_mb = get_MBMB();
   const double alpha_s = calculate_alpha(get_g3());

   mb_DRbar_MZ = calculate_mb_SM5_DRbar(mb_mb, alpha_s, get_MZ());
}

void MSSMNoFV_onshell::convert_yukawa_couplings_treelevel()
{
   Eigen::Matrix<double,3,3> Ye_neu(Eigen::Matrix<double,3,3>::Zero());
   Ye_neu(0, 0) = root2 * get_ME() / get_vd();
   Ye_neu(1, 1) = root2 * get_MM() / get_vd();
   Ye_neu(2, 2) = root2 * get_ML() / get_vd();
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(Eigen::Matrix<double,3,3>::Zero());
   Yu_neu(0, 0) = root2 * get_MU() / get_vu();
   Yu_neu(1, 1) = root2 * get_MC() / get_vu();
   Yu_neu(2, 2) = root2 * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(Eigen::Matrix<double,3,3>::Zero());
   Yd_neu(0, 0) = root2 * get_MD() / get_vd();
   Yd_neu(1, 1) = root2 * get_MS() / get_vd();
   Yd_neu(2, 2) = root2 * get_MB() / get_vd();
   set_Yd(Yd_neu);

   // recalculate trilinear couplings with new Yukawas
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
}

/**
 * Calculate Yukawa couplings with resummed tan(beta) corrections
 *
 * @note The resummations depend on the model parameters ml2, me2, Mu,
 * MassB, MassWB, MassG.  Therefore, this routine needs to be called
 * after all these model parameters have been determined.
 */
void MSSMNoFV_onshell::convert_yukawa_couplings()
{
   Eigen::Matrix<double,3,3> Ye_neu(Eigen::Matrix<double,3,3>::Zero());
   Ye_neu(0, 0) = root2 * get_ME() / get_vd();
   Ye_neu(1, 1) = root2 * get_MM() / get_vd() / (1 + delta_mu_correction(*this));
   Ye_neu(2, 2) = root2 * get_ML() / get_vd() / (1 + delta_tau_correction(*this));
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(Eigen::Matrix<double,3,3>::Zero());
   Yu_neu(0, 0) = root2 * get_MU() / get_vu();
   Yu_neu(1, 1) = root2 * get_MC() / get_vu();
   Yu_neu(2, 2) = root2 * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(Eigen::Matrix<double,3,3>::Zero());
   Yd_neu(0, 0) = root2 * get_MD() / get_vd();
   Yd_neu(1, 1) = root2 * get_MS() / get_vd();
   Yd_neu(2, 2) = root2 * get_MB() / get_vd() / (1 + delta_bottom_correction(*this));
   set_Yd(Yd_neu);

   // recalculate trilinear couplings with new Yukawas
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
}

/**
 * Finds index of neutralino, which is most bino like.
 * @return index in neutralino mass multiplet
 */
unsigned MSSMNoFV_onshell::find_bino_like_neutralino()
{
   // try using pole mass mixing matrix ZN
   if (!gm2calc::is_zero(get_physical().ZN.cwiseAbs().maxCoeff(), eps)) {
      return detail::find_bino_like_neutralino(get_physical().ZN);
   }

   // try using DR mixing matrix ZN
   if (!gm2calc::is_zero(get_ZN().cwiseAbs().maxCoeff(), eps)) {
      calculate_MChi();
   }

   return detail::find_bino_like_neutralino(get_ZN());
}

/**
 * Determines the Mu parameter and the soft-breaking Bino and Wino
 * mass parameters from the two chargino pole masses and the most
 * bino-like neutralino pole mass.  The function uses a fixed-point
 * iteration.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 */
void MSSMNoFV_onshell::convert_Mu_M1_M2(
   double precision_goal,
   unsigned max_iterations)
{
   const auto bino_idx_pole = find_bino_like_neutralino();
   auto bino_idx_DR = detail::find_bino_like_neutralino(get_ZN());
   const auto MCha_goal(get_physical().MCha);
   auto MChi_goal(get_MChi());
   MChi_goal(bino_idx_DR) = get_physical().MChi(bino_idx_pole);

   if (verbose_output) {
      VERBOSE("Converting Mu, M1, M2 to on-shell scheme ...\n"
              "   Goal: MCha = " << pretty_print(MCha_goal.transpose())
              << ", MChi(" << bino_idx_DR << ") = " << MChi_goal(bino_idx_DR)
              << ", accuracy goal = " << precision_goal);
   }

   auto calc_precision = [&MCha_goal, &MChi_goal, this] (unsigned idx) {
      return std::max((MCha_goal - get_MCha()).cwiseAbs().maxCoeff(),
                      std::abs(MChi_goal(idx) - get_MChi(idx)));
   };

   double precision = calc_precision(bino_idx_DR);
   unsigned it = 0;

   while (precision > precision_goal && it < max_iterations) {

      const auto U(get_UM()); // neg. chargino mixing matrix
      const auto V(get_UP()); // pos. chargino mixing matrix
      const auto N(get_ZN()); // neutralino mixing matrix
      const Eigen::Matrix<double,2,2> X = (U.transpose() * MCha_goal.matrix().asDiagonal() * V).real();
      const Eigen::Matrix<double,4,4> Y = (N.transpose() * MChi_goal.matrix().asDiagonal() * N).real();

      if (!X.allFinite()) {
         WARNING("chargino mixing matrix contains NaNs");
         break;
      }

      if (!Y.allFinite()) {
         WARNING("neutralino mixing matrix contains NaNs");
         break;
      }

      set_MassB(Y(0,0));
      set_MassWB(X(0,0));
      set_Mu(X(1,1));

      calculate_MChi();
      calculate_MCha();

      bino_idx_DR = detail::find_bino_like_neutralino(get_ZN());

      MChi_goal = get_MChi();
      MChi_goal(bino_idx_DR) = get_physical().MChi(bino_idx_pole);

      const double old_precision = precision;
      precision = calc_precision(bino_idx_DR);

      if (gm2calc::is_equal(precision, old_precision, eps)) {
         if (verbose_output) {
            VERBOSE("   No improvement in last iteration step, stopping iteration ...");
         }
         break;
      }

      if (verbose_output) {
         VERBOSE("   Iteration " << it << ": Mu = " << get_Mu()
                 << ", M1 = " << get_MassB()
                 << ", M2 = " << get_MassWB()
                 << ", MCha = " << pretty_print(get_MCha().transpose())
                 << ", MChi(" << bino_idx_DR << ") = " << get_MChi(bino_idx_DR)
                 << ", accuracy = " << precision << " GeV");
      }

      it++;
   }

   calculate_DRbar_masses();

   if (precision > precision_goal) {
      get_problems().flag_no_convergence_Mu_MassB_MassWB(precision, max_iterations);
   } else {
      get_problems().unflag_no_convergence_Mu_MassB_MassWB();
   }

   if (verbose_output) {
      if (precision > precision_goal) {
         VERBOSE(
            "   DR-bar to on-shell conversion for Mu, M1 and M2 did"
            " not converge (reached absolute accuracy: " << precision <<
            " GeV, accuracy goal: " << precision_goal <<
            ", max. iterations: " << max_iterations << ")");
      }
      VERBOSE("   Achieved absolute accuracy: " << precision << " GeV");
   }
}

/**
 * Determines soft-breaking left-handed smuon mass parameter from the
 * muon sneutrino pole mass.
 */
void MSSMNoFV_onshell::convert_ml2()
{
   const double MSvmL_pole = get_physical().MSvmL;
   const double vd2 = sqr(get_vd());
   const double vu2 = sqr(get_vu());
   const double g12 = sqr(get_g1());
   const double g22 = sqr(get_g2());

   if (verbose_output) {
      VERBOSE("Converting msl(2,2) to on-shell scheme ...\n"
              "   Using MSvm_pole = " << MSvmL_pole);
   }

   // calculate ml2(1,1) from muon sneutrino pole mass
   const double ml211
      = sqr(MSvmL_pole) + 0.125*(0.6*g12*(vu2 - vd2) + g22*(vu2 - vd2));

   if (std::isfinite(ml211)) {
      set_ml2(1,1,ml211);
      calculate_MSvmL();
   } else {
      WARNING("msl(2,2) is NaN");
   }

   if (verbose_output) {
      VERBOSE("   New msl(2,2) = " << ml211 << ", new MSvmL = " << get_MSvmL());
   }
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 */
void MSSMNoFV_onshell::convert_me2(
   double precision_goal,
   unsigned max_iterations)
{
   double precision = convert_me2_fpi(precision_goal, max_iterations);

   if (precision > precision_goal) {
      precision = convert_me2_root(precision_goal, max_iterations);
   }

   if (precision > precision_goal) {
      get_problems().flag_no_convergence_me2(precision, max_iterations);
   } else {
      get_problems().unflag_no_convergence_me2();
   }
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a one-dimensional root-finding
 * algorithm.
 *
 * @param precision_goal precision goal for the root finding algorithm
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_root_modify(
   double precision_goal,
   unsigned max_iterations)
{
   class Difference_MSm {
   public:
      explicit Difference_MSm(const MSSMNoFV_onshell& model_)
         : model(model_) {}

      double operator()(double me211) {
         model.set_me2(1,1,me211);
         model.calculate_MSm();

         Eigen::Array<double,2,1> MSm_pole(model.get_physical().MSm);
         std::sort(MSm_pole.data(), MSm_pole.data() + MSm_pole.size());
         const int right_index = detail::find_right_like_smuon(model.get_ZM());

         return model.get_MSm(right_index) - MSm_pole(right_index);
      }
   private:
      MSSMNoFV_onshell model;
   };

   boost::uintmax_t it = max_iterations;

   // stopping criterion, given two brackets a, b
   auto Stop_crit = [precision_goal](double a, double b) -> bool {
      return gm2calc::is_equal(a,b,precision_goal);
   };

   if (verbose_output) {
      VERBOSE("Converting mse(2,2) to on-shell scheme with root finder ...");
   }

   // initial guess for the brackets of me2(1,1)
   const double initial_bracket = sqr(1e3 * get_MSm().cwiseAbs().maxCoeff());

   // find the root
   try {
      const std::pair<double,double> root = boost::math::tools::toms748_solve(
         Difference_MSm(*this), 0., initial_bracket, Stop_crit, it);
      set_me2(1, 1, 0.5*(root.first + root.second));
   } catch (const std::exception& e) {
      if (verbose_output) {
         VERBOSE("   DR-bar to on-shell conversion for mse failed with"
                 " root finder: " << e.what());
      }
   }

   calculate_MSm();

   const double precision = std::abs(Difference_MSm(*this)(get_me2(1,1)));

   if (verbose_output) {
      if (it >= max_iterations) {
         VERBOSE(
            "   DR-bar to on-shell conversion for mse did not converge with"
            " root finder (reached absolute accuracy: " << precision <<
            " GeV, accuracy goal: " << precision_goal <<
            ", max. iterations: " << max_iterations << ")");
      }
      VERBOSE("   Achieved absolute accuracy: " << precision << " GeV");
   }

   return precision;
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a one-dimensional root-finding
 * algorithm.
 *
 * If a NaN appears during the FPI, the function resets the
 * soft-breaking right-handed smuon mass parameter to the initial
 * value and returns the maximum double.
 *
 * @param precision_goal precision goal for the root finding algorithm
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_root(
   double precision_goal,
   unsigned max_iterations)
{
   const double me2_save = get_me2(1,1);
   const double precision = convert_me2_root_modify(precision_goal, max_iterations);

   if (!std::isfinite(precision) ||
       !std::isfinite(get_me2(1,1)) ||
       !get_MSm().allFinite() ||
       !get_ZM().allFinite()) {
      // reset
      set_me2(1,1,me2_save);
      calculate_MSm();

      return std::numeric_limits<double>::max();
   }

   return precision;
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a fixed-point iteration.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_fpi_modify(
   double precision_goal,
   unsigned max_iterations)
{
   // sorted pole masses
   Eigen::Array<double,2,1> MSm_pole_sorted(get_physical().MSm);
   std::sort(MSm_pole_sorted.data(),
             MSm_pole_sorted.data() + MSm_pole_sorted.size());

   Eigen::Array<double,2,1> MSm_goal(MSm_pole_sorted);

   int right_index = detail::find_right_like_smuon(get_ZM());

   if (verbose_output) {
      VERBOSE("Converting mse(2,2) to on-shell scheme with FPI ...\n"
              "   Goal: MSm(" << right_index << ") = "
              << MSm_goal(right_index));
   }

   auto calc_precision = [&MSm_goal, this] (unsigned idx) {
      return std::abs(get_MSm(idx) - MSm_goal(idx));
   };

   double precision = calc_precision(right_index);
   unsigned it = 0;

   while (precision > precision_goal && it < max_iterations) {
      const Eigen::Matrix<double,2,2> ZM(get_ZM()); // smuon mixing matrix
      const Eigen::Matrix<double,2,2> M(
         ZM.adjoint() * MSm_goal.square().matrix().asDiagonal() * ZM);

      const double vd2 = sqr(get_vd());
      const double vu2 = sqr(get_vu());
      const double g12 = sqr(get_g1());
      const double ymu2 = std::norm(Ye(1,1));

      const double me211 = M(1,1)
         - (0.5*ymu2*vd2 - 0.15*g12*vd2 + 0.15*g12*vu2);

      if (std::isfinite(me211)) {
         set_me2(1,1,me211);
         calculate_MSm();
      } else {
         WARNING("mse2(2,2) is NaN");
         break;
      }

      if (!std::isfinite(me211) || !get_MSm().allFinite() || !get_ZM().allFinite()) {
         if (verbose_output) {
            VERBOSE("   NaN appearing in DR-bar to on-shell conversion"
                    " for mse with FPI:\n"
                    "      mse2(2,2) = " << me211 <<
                    ", MSm = " << pretty_print(get_MSm().transpose()) <<
                    ", ZM = " << pretty_print(get_ZM()));
         }
         return std::numeric_limits<double>::max();
      }

      right_index = detail::find_right_like_smuon(get_ZM());

      MSm_goal = get_MSm();
      MSm_goal(right_index) = MSm_pole_sorted(right_index);

      const double old_precision = precision;
      precision = calc_precision(right_index);

      if (gm2calc::is_equal(precision, old_precision, eps)) {
         if (verbose_output) {
            VERBOSE("   No improvement in last iteration step, stopping iteration ...");
         }
         break;
      }

      if (verbose_output) {
         VERBOSE("   Iteration " << it << ": mse(2,2) = "
                 << signed_abs_sqrt(me211) << ", MSm(" << right_index
                 << ") = " << get_MSm(right_index)
                 << ", accuracy = " << precision);
      }

      it++;
   }

   if (verbose_output) {
      if (precision > precision_goal) {
         VERBOSE(
            "   DR-bar to on-shell conversion for mse did not converge with"
            " FPI (reached absolute accuracy: " << precision <<
            " GeV, accuracy goal: " << precision_goal <<
            ", max. iterations: " << max_iterations << ")");
      }
      VERBOSE("   Achieved absolute accuracy: " << precision << " GeV");
   }

   return precision;
}

/**
 * Determines soft-breaking right-handed smuon mass parameter from one
 * smuon pole mass.  The function uses a fixed-point iteration.
 *
 * If a NaN appears during the FPI, the function resets the
 * soft-breaking right-handed smuon mass parameter to the initial
 * value and returns the maximum double.
 *
 * @param precision_goal precision goal of iteration
 * @param max_iterations maximum number of iterations
 * @return achieved precision
 */
double MSSMNoFV_onshell::convert_me2_fpi(
   double precision_goal,
   unsigned max_iterations)
{
   const double me2_save = get_me2(1,1);
   const double precision = convert_me2_fpi_modify(precision_goal, max_iterations);

   if (!std::isfinite(precision) ||
       !std::isfinite(get_me2(1,1)) ||
       !get_MSm().allFinite() ||
       !get_ZM().allFinite()) {
      // reset
      set_me2(1,1,me2_save);
      calculate_MSm();

      return std::numeric_limits<double>::max();
   }

   return precision;
}

std::ostream& operator<<(std::ostream& os, const MSSMNoFV_onshell& model)
{
   auto sas = [] (double x) { return gm2calc::signed_abs_sqrt(x); };

   Eigen::Array<double,3,1> yd_non_res;
   yd_non_res << (root2 * model.get_MD() / model.get_vd()),
                 (root2 * model.get_MS() / model.get_vd()),
                 (root2 * model.get_MB() / model.get_vd());

   Eigen::Array<double,3,1> ye_non_res;
   ye_non_res << (root2 * model.get_ME() / model.get_vd()),
                 (root2 * model.get_MM() / model.get_vd()),
                 (root2 * model.get_ML() / model.get_vd());

   os <<
      "===============================\n"
      " GM2Calc masses and parameters \n"
      "===============================\n"
      "1/alpha(MZ)      = " << 1./calculate_alpha(model.get_EL()) << '\n' <<
      "1/alpha(0)       = " << 1./calculate_alpha(model.get_EL0()) << '\n' <<
      "alpha_s(MZ)      = " <<    calculate_alpha(model.get_g3()) << '\n' <<
      "MM pole          = " << model.get_MM() << " GeV\n" <<
      "MT               = " << model.get_MT() << " GeV\n" <<
      "mb(mb) MS-bar    = " << model.get_MBMB() << " GeV\n" <<
      "mb(MZ) DR-bar    = " << model.get_MB() << " GeV\n" <<
      "MTau             = " << model.get_ML() << " GeV\n" <<
      "MW pole          = " << model.get_MW() << " GeV\n" <<
      "MZ pole          = " << model.get_MZ() << " GeV\n" <<
      "MSm pole         = " << pretty_print(model.get_MSm().transpose()) << " GeV\n" <<
      "USm pole         = " << pretty_print(model.get_USm()) << '\n' <<
      "MSvm pole        = " << model.get_MSvmL() << " GeV\n" <<
      "MSb              = " << pretty_print(model.get_MSb().transpose()) << " GeV\n" <<
      "USb              = " << pretty_print(model.get_USb()) << '\n' <<
      "MSt              = " << pretty_print(model.get_MSt().transpose()) << " GeV\n" <<
      "USt              = " << pretty_print(model.get_USt()) << '\n' <<
      "MStau            = " << pretty_print(model.get_MStau().transpose()) << " GeV\n" <<
      "UStau            = " << pretty_print(model.get_UStau()) << '\n' <<
      "MCha pole        = " << pretty_print(model.get_MCha().transpose()) << " GeV\n" <<
      "UM pole          = " << pretty_print(model.get_UM()) << '\n' <<
      "UP pole          = " << pretty_print(model.get_UP()) << '\n' <<
      "MChi pole        = " << pretty_print(model.get_MChi().transpose()) << " GeV\n" <<
      "ZN pole          = {" << pretty_print(model.get_ZN().row(0)) << ",\n" <<
      "                    " << pretty_print(model.get_ZN().row(1)) << ",\n" <<
      "                    " << pretty_print(model.get_ZN().row(2)) << ",\n" <<
      "                    " << pretty_print(model.get_ZN().row(3)) << "}\n" <<
      "MA0              = " << model.get_MA0() << " GeV\n" <<
      "MH               = " << pretty_print(model.get_Mhh().transpose()) << " GeV\n" <<
      "tan(beta) DR-bar = " << model.get_TB() << '\n' <<
      "yu               = " << pretty_print(model.get_Yu().diagonal().transpose()) << '\n' <<
      "yd resummed      = " << pretty_print(model.get_Yd().diagonal().transpose()) << '\n' <<
      "ye resummed      = " << pretty_print(model.get_Ye().diagonal().transpose()) << '\n' <<
      "yd non res.      = " << pretty_print(yd_non_res.transpose()) << '\n' <<
      "ye non res.      = " << pretty_print(ye_non_res.transpose()) << '\n' <<
      "Mu on-shell      = " << model.get_Mu() << " GeV\n" <<
      "M1 on-shell      = " << model.get_MassB() << " GeV\n" <<
      "M2 on-shell      = " << model.get_MassWB() << " GeV\n" <<
      "M3               = " << model.get_MassG() << " GeV\n" <<
      "msl on-shell     = " << pretty_print(model.get_ml2().diagonal().transpose().unaryExpr(sas)) << " GeV\n" <<
      "mse on-shell     = " << pretty_print(model.get_me2().diagonal().transpose().unaryExpr(sas)) << " GeV\n" <<
      "msq              = " << pretty_print(model.get_mq2().diagonal().transpose().unaryExpr(sas)) << " GeV\n" <<
      "msu              = " << pretty_print(model.get_mu2().diagonal().transpose().unaryExpr(sas)) << " GeV\n" <<
      "msd              = " << pretty_print(model.get_md2().diagonal().transpose().unaryExpr(sas)) << " GeV\n" <<
      "Au               = " << pretty_print(model.get_Au().diagonal().transpose()) << " GeV\n" <<
      "Ad               = " << pretty_print(model.get_Ad().diagonal().transpose()) << " GeV\n" <<
      "Ae DR-bar        = " << pretty_print(model.get_Ae().diagonal().transpose()) << " GeV\n" <<
      "ren. scale       = " << model.get_scale() << " GeV\n"
      ;

   return os;
}

} // namespace gm2calc
