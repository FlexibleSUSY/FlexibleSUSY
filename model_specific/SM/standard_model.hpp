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

/**
 * @file standard_model.hpp
 *
 * @brief contains class for Standard model running and self-energies
 *
 */

#ifndef STANDARD_MODEL_H
#define STANDARD_MODEL_H

#include "betafunction.hpp"
#include "standard_model_physical.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "error.hpp"
#include "names.hpp"
#include "problems.hpp"
#include "physical_input.hpp"

#include <array>
#include <iosfwd>
#include <string>

#include <Eigen/Core>

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {

class EWSB_solver;

namespace standard_model_info {

   enum Particles : int {VG, Hp, Fv, Ah, hh, VP, VZ, Fd, Fu, Fe, VWp,
      NUMBER_OF_PARTICLES};

   enum Parameters : int {g1, g2, g3, Lambdax, Yu0_0, Yu0_1, Yu0_2, Yu1_0,
      Yu1_1, Yu1_2, Yu2_0, Yu2_1, Yu2_2, Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2
      , Yd2_0, Yd2_1, Yd2_2, Ye0_0, Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0,
      Ye2_1, Ye2_2, mu2, v, NUMBER_OF_PARAMETERS};

   enum Mixings : int {ReVd00, ImVd00, ReVd01, ImVd01, ReVd02, ImVd02,
      ReVd10, ImVd10, ReVd11, ImVd11, ReVd12, ImVd12, ReVd20, ImVd20, ReVd21,
      ImVd21, ReVd22, ImVd22, ReUd00, ImUd00, ReUd01, ImUd01, ReUd02, ImUd02,
      ReUd10, ImUd10, ReUd11, ImUd11, ReUd12, ImUd12, ReUd20, ImUd20, ReUd21,
      ImUd21, ReUd22, ImUd22, ReVu00, ImVu00, ReVu01, ImVu01, ReVu02, ImVu02,
      ReVu10, ImVu10, ReVu11, ImVu11, ReVu12, ImVu12, ReVu20, ImVu20, ReVu21,
      ImVu21, ReVu22, ImVu22, ReUu00, ImUu00, ReUu01, ImUu01, ReUu02, ImUu02,
      ReUu10, ImUu10, ReUu11, ImUu11, ReUu12, ImUu12, ReUu20, ImUu20, ReUu21,
      ImUu21, ReUu22, ImUu22, ReVe00, ImVe00, ReVe01, ImVe01, ReVe02, ImVe02,
      ReVe10, ImVe10, ReVe11, ImVe11, ReVe12, ImVe12, ReVe20, ImVe20, ReVe21,
      ImVe21, ReVe22, ImVe22, ReUe00, ImUe00, ReUe01, ImUe01, ReUe02, ImUe02,
      ReUe10, ImUe10, ReUe11, ImUe11, ReUe12, ImUe12, ReUe20, ImUe20, ReUe21,
      ImUe21, ReUe22, ImUe22, NUMBER_OF_MIXINGS};

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = false;

   class Standard_model_particle_names : public Names {
   public:
      virtual ~Standard_model_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   class Standard_model_parameter_names : public Names {
   public:
      virtual ~Standard_model_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   const Standard_model_particle_names  particle_names_getter{};
   const Standard_model_parameter_names parameter_names_getter{};

} // namespace standard_model_info

namespace standard_model {

template <class T>
class StandardModel;

/**
 * @class Standard_model
 * @brief model class with routines for SM running and self-energies
 */
class Standard_model : public Beta_function {
public:

   Standard_model();
   Standard_model(double scale_, double loops_, double thresholds_
   , double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<
   double,3,3>& Yu_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<
   double,3,3>& Ye_, double mu2_, double v_);
   Standard_model(const Standard_model&) = default;
   Standard_model(Standard_model&&) = default;
   virtual ~Standard_model() = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   Standard_model& operator=(const Standard_model&) = default;
   Standard_model& operator=(Standard_model&&) = default;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 1;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(int);
   void set_loop_corrections(const Loop_corrections&);
   const Loop_corrections& get_loop_corrections() const;
   void set_threshold_corrections(const Threshold_corrections&);
   const Threshold_corrections& get_threshold_corrections() const;
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   void set_physical(const Standard_model_physical&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   const Standard_model_physical& get_physical() const;
   Standard_model_physical& get_physical();
   const Problems& get_problems() const;
   Problems& get_problems();
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   void print() const;
   void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&) override;

   Standard_model calc_beta() const;
   Standard_model calc_beta(int) const;
   void clear();
   void clear_running_parameters();
   void clear_DRbar_parameters();
   void clear_problems();

   void calculate_spectrum();
   std::string name() const;
   virtual void run_to(double scale, double eps = -1.0) override;
   void set_precision(double);
   double get_precision() const;

   void set_physical_input(const Physical_input& input_) { input = input_; }
   const Physical_input& get_physical_input() const { return input; }
   Physical_input& get_physical_input() { return input; }

   void initialise_from_input(const softsusy::QedQcd&);

   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_Lambdax(double Lambdax_) { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, double value) { Yu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, double value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, double value) { Ye(i,k) = value; }
   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }

   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_Lambdax() const { return Lambdax; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }
   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   double get_mu2() const { return mu2; }
   double get_v() const { return v; }

   double get_MVG() const { return MVG; }
   double get_MHp() const { return MHp; }
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   double get_MAh() const { return MAh; }
   double get_Mhh() const { return Mhh; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   double get_MVWp() const { return MVWp; }
   const Eigen::Array<double,2,1>& get_MVPVZ() const { return MVPVZ; }
   double get_MVPVZ(int i) const { return MVPVZ(i); }

   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const { return Vd; }
   const std::complex<double>& get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const { return Ud; }
   const std::complex<double>& get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const { return Vu; }
   const std::complex<double>& get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const { return Uu; }
   const std::complex<double>& get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const { return Ve; }
   const std::complex<double>& get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const { return Ue; }
   const std::complex<double>& get_Ue(int i, int k) const { return Ue(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const { return ZZ; }
   double get_ZZ(int i, int k) const { return ZZ(i,k); }

   double get_mass_matrix_VG() const;
   void calculate_MVG();
   double get_mass_matrix_Hp() const;
   void calculate_MHp();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   double get_mass_matrix_Ah() const;
   void calculate_MAh();
   double get_mass_matrix_hh() const;
   void calculate_Mhh();
   double get_mass_matrix_VP() const;
   void calculate_MVP();
   double get_mass_matrix_VZ() const;
   void calculate_MVZ();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   double get_mass_matrix_VWp() const;
   void calculate_MVWp();
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const;
   void calculate_MVPVZ();

   double get_ewsb_eq_hh_1() const;

   double CpconjHpHphh() const;
   double CpconjHpVWpVP() const;
   double CpconjHpVZVWp() const;
   double CpHpgWpCbargZ() const;
   double CpconjHpbargWpCgZ() const;
   double CpHpgZbargWp() const;
   double CpconjHpbargZgWp() const;
   double CpHpconjHpAhAh() const;
   double CpHpconjHphhhh() const;
   double CpHpconjHpconjHpHp() const;
   std::complex<double> CpconjHpVWpAh() const;
   double CpconjHpVWphh() const;
   double CpconjHpVPHp() const;
   double CpconjHpVZHp() const;
   double CpHpconjHpconjVWpVWp() const;
   std::complex<double> CpHpconjHpVZVZ() const;
   std::complex<double> CpconjHpbarFdFuPR(int gI1, int gI2) const;
   std::complex<double> CpconjHpbarFdFuPL(int gI1, int gI2) const;
   double CpconjHpbarFeFvPR(int , int ) const;
   std::complex<double> CpconjHpbarFeFvPL(int gI1, int gI2) const;
   double CpAhhhAh() const;
   std::complex<double> CpAhbargWpgWp() const;
   std::complex<double> CpAhbargWpCgWpC() const;
   double CpAhAhAhAh() const;
   double CpAhAhhhhh() const;
   double CpAhAhconjHpHp() const;
   std::complex<double> CpAhVZhh() const;
   std::complex<double> CpAhconjVWpHp() const;
   double CpAhAhconjVWpVWp() const;
   std::complex<double> CpAhAhVZVZ() const;
   std::complex<double> CpAhbarFdFdPR(int gI1, int gI2) const;
   std::complex<double> CpAhbarFdFdPL(int gI1, int gI2) const;
   std::complex<double> CpAhbarFeFePR(int gI1, int gI2) const;
   std::complex<double> CpAhbarFeFePL(int gI1, int gI2) const;
   std::complex<double> CpAhbarFuFuPR(int gI1, int gI2) const;
   std::complex<double> CpAhbarFuFuPL(int gI1, int gI2) const;
   double CphhAhAh() const;
   double Cphhhhhh() const;
   double CphhVZVZ() const;
   double CphhbargWpgWp() const;
   double CphhbargWpCgWpC() const;
   double CphhbargZgZ() const;
   double CphhconjHpHp() const;
   double CphhconjVWpVWp() const;
   double CphhhhAhAh() const;
   double Cphhhhhhhh() const;
   double CphhhhconjHpHp() const;
   std::complex<double> CphhVZAh() const;
   double CphhconjVWpHp() const;
   double CphhhhconjVWpVWp() const;
   std::complex<double> CphhhhVZVZ() const;
   std::complex<double> CphhbarFdFdPR(int gI1, int gI2) const;
   std::complex<double> CphhbarFdFdPL(int gI1, int gI2) const;
   std::complex<double> CphhbarFeFePR(int gI1, int gI2) const;
   std::complex<double> CphhbarFeFePL(int gI1, int gI2) const;
   std::complex<double> CphhbarFuFuPR(int gI1, int gI2) const;
   std::complex<double> CphhbarFuFuPL(int gI1, int gI2) const;
   std::complex<double> CpVZhhAh() const;
   double CpVZVZhh() const;
   double CpVZbargWpgWp() const;
   double CpVZbargWpCgWpC() const;
   double CpVZconjHpHp() const;
   double CpVZconjVWpHp() const;
   std::complex<double> CpVZVZAhAh() const;
   std::complex<double> CpVZVZhhhh() const;
   std::complex<double> CpVZVZconjHpHp() const;
   double CpVZconjVWpVWp() const;
   double CpVZbarFdFdPL(int gI1, int gI2) const;
   double CpVZbarFdFdPR(int gI1, int gI2) const;
   double CpVZbarFeFePL(int gI1, int gI2) const;
   double CpVZbarFeFePR(int gI1, int gI2) const;
   double CpVZbarFuFuPL(int gI1, int gI2) const;
   double CpVZbarFuFuPR(int gI1, int gI2) const;
   double CpVZbarFvFvPL(int gI1, int gI2) const;
   double CpVZbarFvFvPR(int , int ) const;
   double CpVZVZconjVWpVWp1() const;
   double CpVZVZconjVWpVWp2() const;
   double CpVZVZconjVWpVWp3() const;
   std::complex<double> CpconjVWpHpAh() const;
   double CpconjVWpHphh() const;
   double CpconjVWpVPHp() const;
   double CpconjVWpVWphh() const;
   double CpconjVWpVZHp() const;
   double CpconjVWpbargPgWp() const;
   double CpconjVWpbargWpCgP() const;
   double CpconjVWpbargWpCgZ() const;
   double CpconjVWpbargZgWp() const;
   double CpVWpconjVWpAhAh() const;
   double CpVWpconjVWphhhh() const;
   double CpVWpconjVWpconjHpHp() const;
   double CpconjVWpVWpVP() const;
   double CpconjVWpVZVWp() const;
   std::complex<double> CpconjVWpbarFdFuPL(int gI1, int gI2) const;
   double CpconjVWpbarFdFuPR(int , int ) const;
   std::complex<double> CpconjVWpbarFeFvPL(int gI1, int gI2) const;
   double CpconjVWpbarFeFvPR(int , int ) const;
   double CpVWpconjVWpVPVP1() const;
   double CpVWpconjVWpVPVP2() const;
   double CpVWpconjVWpVPVP3() const;
   double CpVWpconjVWpVZVZ1() const;
   double CpVWpconjVWpVZVZ2() const;
   double CpVWpconjVWpVZVZ3() const;
   double CpVWpconjVWpconjVWpVWp1() const;
   double CpVWpconjVWpconjVWpVWp2() const;
   double CpVWpconjVWpconjVWpVWp3() const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFdhhFdPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdhhFdPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFdVGFdPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdVGFdPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdVPFdPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdVPFdPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdVZFdPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdVZFdPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdconjHpFuPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdconjHpFuPR(int gO1, int gI2) const;
   double CpbarUFdconjVWpFuPR(int , int ) const;
   std::complex<double> CpbarUFdconjVWpFuPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFuhhFuPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuhhFuPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFuHpFdPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuHpFdPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFuVGFuPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuVGFuPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuVPFuPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuVPFuPL(int gO1, int gI2) const;
   double CpbarUFuVWpFdPR(int , int ) const;
   std::complex<double> CpbarUFuVWpFdPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuVZFuPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuVZFuPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFehhFePL(int gO2, int gI2) const;
   std::complex<double> CpbarUFehhFePR(int gO1, int gI2) const;
   std::complex<double> CpbarUFeVPFePR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeVPFePL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeVZFePR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeVZFePL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeconjHpFvPL(int gO2, int gI2) const;
   double CpbarUFeconjHpFvPR(int , int ) const;
   double CpbarUFeconjVWpFvPR(int , int ) const;
   double CpbarUFeconjVWpFvPL(int gO1, int gI2) const;
   std::complex<double> CpbarFdFdAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarFdFdAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarFdhhFdPL(int gO2, int gI2) const;
   std::complex<double> CpbarFdhhFdPR(int gO1, int gI2) const;
   double CpbarFdVZFdPR(int gO2, int gI2) const;
   double CpbarFdVZFdPL(int gO1, int gI2) const;
   std::complex<double> CpbarFdconjHpFuPL(int gO2, int gI2) const;
   std::complex<double> CpbarFdconjHpFuPR(int gO1, int gI2) const;
   double CpbarFdconjVWpFuPR(int , int ) const;
   std::complex<double> CpbarFdconjVWpFuPL(int gO1, int gI2) const;
   std::complex<double> CpbarFeFeAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarFeFeAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarFehhFePL(int gO2, int gI2) const;
   std::complex<double> CpbarFehhFePR(int gO1, int gI2) const;
   double CpbarFeVZFePR(int gO2, int gI2) const;
   double CpbarFeVZFePL(int gO1, int gI2) const;
   std::complex<double> CpbarFeconjHpFvPL(int gO2, int gI2) const;
   double CpbarFeconjHpFvPR(int , int ) const;
   double CpbarFeconjVWpFvPR(int , int ) const;
   std::complex<double> CpbarFeconjVWpFvPL(int gO1, int gI2) const;
   std::complex<double> CpbarFuFuAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarFuFuAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarFuhhFuPL(int gO2, int gI2) const;
   std::complex<double> CpbarFuhhFuPR(int gO1, int gI2) const;
   std::complex<double> CpbarFuHpFdPL(int gO2, int gI2) const;
   std::complex<double> CpbarFuHpFdPR(int gO1, int gI2) const;
   double CpbarFuVPFuPR(int gO2, int gI2) const;
   double CpbarFuVPFuPL(int gO1, int gI2) const;
   double CpbarFuVWpFdPR(int , int ) const;
   std::complex<double> CpbarFuVWpFdPL(int gO1, int gI2) const;
   double CpbarFuVZFuPR(int gO2, int gI2) const;
   double CpbarFuVZFuPL(int gO1, int gI2) const;
   std::complex<double> self_energy_Hp_1loop(double p ) const;
   std::complex<double> self_energy_Ah_1loop(double p ) const;
   std::complex<double> self_energy_hh_1loop(double p ) const;
   std::complex<double> self_energy_VZ_1loop(double p ) const;
   std::complex<double> self_energy_VWp_1loop(double p ) const;
   std::complex<double> self_energy_Fd_1loop_1(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_1(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fe_1loop_1(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const;

   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR(double p) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL(double p) const;

   std::complex<double> tadpole_hh_1loop() const;

   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole_equations() const;

   /// calculates Higgs 2-loop self-energy
   double self_energy_hh_2loop(double p) const;
   /// calculates Higgs 3-loop self-energy
   double self_energy_hh_3loop() const;
   /// calculates Higgs 4-loop self-energy
   double self_energy_hh_4loop() const;

   void calculate_MVG_pole();
   void calculate_MFv_pole();
   void calculate_Mhh_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFe_pole();
   void calculate_MVWp_pole();
   double calculate_MVWp_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MVP_DRbar(double);
   double calculate_MVZ_DRbar(double);
   double calculate_MVWp_DRbar(double);
   double calculate_Mhh_DRbar(double);

   double ThetaW() const;

   double calculate_delta_alpha_em(double alphaEm) const;
   double calculate_delta_alpha_s(double alphaS) const;
   void calculate_Lambdax_DRbar();
   double calculate_theta_w(const softsusy::QedQcd&, double alpha_em_drbar);
   void calculate_Yu_DRbar(const softsusy::QedQcd&);
   void calculate_Yd_DRbar(const softsusy::QedQcd&);
   void calculate_Ye_DRbar(const softsusy::QedQcd&);
   double recalculate_mw_pole(double);
   double max_rel_diff(const Standard_model& old) const;

protected:

   // Running parameters
   double g1{};
   double g2{};
   double g3{};
   double Lambdax{};
   Eigen::Matrix<double,3,3> Yu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Yd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Ye{Eigen::Matrix<double,3,3>::Zero()};
   double mu2{};
   double v{};

private:

   static const int numberOfParameters = 33;

   struct Beta_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};
      double traceYdAdjYuYuAdjYd{};
      double traceYdAdjYdYdAdjYdYdAdjYd{};
      double traceYdAdjYdYdAdjYuYuAdjYd{};
      double traceYdAdjYuYuAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYuYuAdjYu{};
   };
   void calc_beta_traces(Beta_traces&) const;

   double calc_beta_g1_one_loop(const Beta_traces&) const;
   double calc_beta_g1_two_loop(const Beta_traces&) const;
   double calc_beta_g1_three_loop(const Beta_traces&) const;
   double calc_beta_g2_one_loop(const Beta_traces&) const;
   double calc_beta_g2_two_loop(const Beta_traces&) const;
   double calc_beta_g2_three_loop(const Beta_traces&) const;
   double calc_beta_g3_one_loop(const Beta_traces&) const;
   double calc_beta_g3_two_loop(const Beta_traces&) const;
   double calc_beta_g3_three_loop(const Beta_traces&) const;
   double calc_beta_g3_four_loop(const Beta_traces&) const;
   double calc_beta_g3_five_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_one_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_two_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_three_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_four_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_one_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_two_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_three_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_four_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_one_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_two_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_three_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_one_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_two_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_three_loop(const Beta_traces&) const;
   double calc_beta_mu2_one_loop(const Beta_traces&) const;
   double calc_beta_mu2_two_loop(const Beta_traces&) const;
   double calc_beta_mu2_three_loop(const Beta_traces&) const;
   double calc_beta_v_one_loop(const Beta_traces&) const;
   double calc_beta_v_two_loop(const Beta_traces&) const;
   double calc_beta_v_three_loop(const Beta_traces&) const;

   using EWSB_vector_t = Eigen::Matrix<double,number_of_ewsb_equations,1>;

   class EEWSBStepFailed : public Error {
   public:
      EEWSBStepFailed() : Error("Could not perform EWSB step") {}
      virtual ~EEWSBStepFailed() = default;
   };

   int ewsb_loop_order{4};
   int pole_mass_loop_order{4};
   bool force_output{false};      ///< switch to force output of pole masses
   double precision{1e-4};        ///< RG running precision
   double ewsb_iteration_precision{1e-5};
   Standard_model_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{standard_model_info::model_name,
                     &standard_model_info::particle_names_getter,
                     &standard_model_info::parameter_names_getter};
   Loop_corrections loop_corrections{}; ///< used loop pole mass corrections
   Threshold_corrections threshold_corrections{}; ///< used low-energy threshold corrections
   Physical_input input{};

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(int);
   int solve_ewsb_iteratively_with(EWSB_solver*, const Eigen::Matrix<double, number_of_ewsb_equations, 1>&);
   int solve_ewsb_tree_level_custom();
   EWSB_vector_t ewsb_initial_guess();
   EWSB_vector_t ewsb_step() const;
   void copy_DRbar_masses_to_pole_masses();

   void initial_guess_for_parameters(const softsusy::QedQcd&);
   bool check_convergence(const Standard_model& old) const;

   // Passarino-Veltman loop functions
   double A0(double) const;
   double B0(double, double, double) const;
   double B1(double, double, double) const;
   double B00(double, double, double) const;
   double B22(double, double, double) const;
   double H0(double, double, double) const;
   double F0(double, double, double) const;
   double G0(double, double, double) const;

   // DR-bar masses
   double MVG{};
   double MHp{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double MAh{};
   double Mhh{};
   double MVP{};
   double MVZ{};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   double MVWp{};
   Eigen::Array<double,2,1> MVPVZ{Eigen::Array<double,2,1>::Zero()};

   // DR-bar mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};


};

std::ostream& operator<<(std::ostream&, const Standard_model&);

} // namespace standard_model

} // namespace flexiblesusy

#endif
