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
    * @class Data
    * @brief Model parameters necessary for calculating weak mixing angle
    *
    * @attention The W and Z self-energies are assumed to be
    * calculated using the top quark pole mass, instead of the top
    * quark DR-bar mass.
    */
   struct Data {
      Data();

      double scale;                  ///< renormalization scale
      double alpha_em_drbar;         ///< alpha_em(MZ, DR-bar, SUSY)
      double fermi_contant;          ///< Fermi constant
      double self_energy_z_at_mz;    ///< self-energy Z at p = MZ, mt = mt_pole
      double self_energy_w_at_0;     ///< self-energy W at p = 0, mt = mt_pole
      double self_energy_w_at_mw;    ///< self-energy W at p = MW, mt = mt_pole
      double mw_pole;                ///< W pole mass
      double mz_pole;                ///< Z pole mass
      double mt_pole;                ///< top quark pole mass
      double mh_drbar;               ///< lightest CP-even Higgs DR-bar mass
      double hmix_12;                ///< CP-even Higgs mixing Cos(alpha)
      double msel_drbar;             ///< left-handed selectron DR-bar mass
      double msmul_drbar;            ///< left-handed smuon DR-bar mass
      double msve_drbar;             ///< electron-sneutrino DR-bar mass
      double msvm_drbar;             ///< muon-sneutrino DR-bar mass
      Eigen::ArrayXd mn_drbar;       ///< Neutralino DR-bar mass
      Eigen::ArrayXd mc_drbar;       ///< Chargino DR-bar mass
      Eigen::MatrixXcd zn;           ///< Neutralino mixing matrix
      Eigen::MatrixXcd um;           ///< Chargino mixing matrix
      Eigen::MatrixXcd up;           ///< Chargino mixing matrix
      double gY;                     ///< U(1)_Y gauge coupling
      double g2;                     ///< SU(2)_L gauge coupling
      double g3;                     ///< SU(3)_c gauge coupling
      double tan_beta;               ///< tan(beta) = vu / vd
   };

   struct Self_energy_data {
      Self_energy_data();
      double scale;                  ///< renormalization scale
      double mt_pole;                ///< top quark pole mass
      double mt_drbar;               ///< top quark DR-bar mass
      double mb_drbar;               ///< bottom quark DR-bar mass
      double gY;                     ///< U(1)_Y gauge coupling
      double g2;                     ///< SU(2)_L gauge coupling
   };

   Weinberg_angle();

   void enable_susy_contributions(); ///< enable susy contributions
   void disable_susy_contributions(); ///< disable susy contributions

   void set_data(const Data&);       ///< set data necessary for the calculation
   void set_number_of_iterations(int); ///< maximum number of iterations
   void set_number_of_loops(int);    ///< set number of loops
   void set_precision_goal(double);  ///< set precision goal
   double get_rho_hat() const;       ///< returns the rho parameter
   double get_sin_theta() const;     ///< returns sin(theta_w)

   /// calculates the sinus of the Weinberg angle
   int calculate(double rho_start = 1.0, double sin_start = 0.48);

   static double replace_mtop_in_self_energy_z(double, double, const Self_energy_data&);
   static double replace_mtop_in_self_energy_w(double, double, const Self_energy_data&);

private:
   int number_of_iterations; ///< maximum number of iterations
   int number_of_loops;      ///< number of loops
   double precision_goal;         ///< precision goal
   double rho_hat;                ///< output rho-hat parameter
   double sin_theta;              ///< output sin(theta)
   Data data;
   bool susy_contributions;       ///< model type

   static double calculate_delta_r(double, double, const Data&, bool add_susy_contributions = true, int number_of_loops = 2);
   static double calculate_delta_rho(double, double, const Data&, bool add_susy_contributions = true, int number_of_loops = 2);
   static double calculate_delta_vb(double, double, const Data&, bool add_susy_contributions = true);
   static double calculate_delta_vb_sm(double, double, const Data&);
   static double calculate_delta_vb_susy(double, const Data&);
   static double rho_2(double);

   static double calculate_self_energy_z_top(double, double, const Self_energy_data&);
   static double calculate_self_energy_w_top(double, double, const Self_energy_data&);
};

} // namespace weinberg_angle

} // namespace flexiblesusy

#endif
