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

#ifndef WEINBERG_ANGLE_H
#define WEINBERG_ANGLE_H

#include "standard_model.hpp"
#include "ew_input.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

namespace weinberg_angle {

/**
 * @class Weinberg_angle
 * @brief Class to calculate the DR-bar weak mixing angle
 */
class Weinberg_angle {
public:
   /**
    * @class Sm_parameters
    * @brief SM parameters necessary for calculating the weak mixing angle
    */
   struct Sm_parameters {
      double fermi_constant{0.}; ///< Fermi constant
      double mw_pole{0.};        ///< W pole mass
      double mz_pole{0.};        ///< Z pole mass
      double mt_pole{0.};        ///< top quark pole mass
      double mh_pole{Electroweak_constants::MH}; ///< Higgs pole mass
      double alpha_s{0.};        ///< strong coupling at Q = mt_pole
      double alpha_s_mz{0.};     ///< strong coupling at Q = mz_pole
      double dalpha_s_5_had{Electroweak_constants::delta_alpha_s_5_had}; ///< 5-flavour hadronic contributions
   };

   Weinberg_angle();
   Weinberg_angle(const standard_model::Standard_model*, const Sm_parameters&);

   void set_number_of_iterations(int); ///< maximum number of iterations
   void set_number_of_loops(int);    ///< set number of loops
   void set_precision_goal(double);  ///< set precision goal
   double get_rho_hat() const;       ///< returns the rho parameter
   double get_sin_theta() const;     ///< returns sin(theta_w)

   /// calculates and returns the sine of the Weinberg angle
   double calculate(double sinThetaW_start = 0.48);
private:
   int number_of_iterations; ///< maximum number of iterations
   int number_of_loops;      ///< number of loops
   double precision_goal;         ///< precision goal
   double rho_hat;                ///< output rho-hat parameter
   double sin_theta;              ///< output sin(theta)

   const standard_model::Standard_model* model{nullptr}; ///< pointer to investigated model
   Sm_parameters sm_parameters{};     ///< SM parameters
   double calculate_self_energy_VZ(double p) const;
   double calculate_self_energy_VZ_top(double p, double mt) const;
   double calculate_self_energy_VWp(double p) const;
   double calculate_self_energy_VWp_top(double p, double mt) const;
   double calculate_delta_rho_hat(double sinThetaW) const;
   double pizzt_MZ{0.};               ///< transverse Z self-energy at p^2 = MZ^2
   double piwwt_MW{0.};               ///< transverse W self-energy at p^2 = MW^2
   double piwwt_0{0.};                ///< transverse W self-energy at p^2 = 0
   double calculate_delta_alpha_hat_bsm(double alpha_em) const;
   double calculate_delta_r_hat(double rhohat_ratio, double sinThetaW) const;
   double calculate_delta_vb(double rhohat_ratio, double sinThetaW) const;
   double calculate_delta_vb_sm(double sinThetaW) const;

   static double rho_2(double);
};

} // namespace weinberg_angle

} // namespace flexiblesusy

#endif
