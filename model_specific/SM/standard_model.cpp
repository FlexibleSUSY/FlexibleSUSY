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
 * @file standard_model.cpp
 * @brief implementation of the Standard_model class
 *
 */

#include "standard_model.hpp"
#include "eigen_utils.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "config.h"
#include "loop_libraries/loop_library.hpp"
#include "raii.hpp"
#include "thread_pool.hpp"
#include "functors.hpp"
#include "ew_input.hpp"
#include "weinberg_angle.hpp"

#include "sm_twoloophiggs.hpp"
#include "sm_threeloophiggs.hpp"
#include "sm_fourloophiggs.hpp"
#include "sm_fourloop_as.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>

namespace flexiblesusy {

namespace standard_model_info {

   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {
      1, 1, 3, 1, 1, 1, 1, 3, 3, 3, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {
      "VG", "Hp", "Fv", "Ah", "hh", "VP", "VZ", "Fd", "Fu", "Fe", "VWp"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "H^+", R"(\nu)", "A^0", "h", R"(\gamma)", "Z", "d", "u", "e", "W^+"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {
      "g1", "g2", "g3",
      "Lambdax", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)", "Yu(1,1)", "Yu(1,2)",
      "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "Yd(0,0)", "Yd(0,1)", "Yd(0,2)", "Yd(1,0)"
      , "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)", "Ye(0,0)",
      "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)", "Ye(2,0)", "Ye(2,1)",
      "Ye(2,2)", "mu2", "v"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "Re(Vd(0,0))",
      "Im(Vd(0,0))", "Re(Vd(0,1))", "Im(Vd(0,1))", "Re(Vd(0,2))", "Im(Vd(0,2))",
      "Re(Vd(1,0))", "Im(Vd(1,0))", "Re(Vd(1,1))", "Im(Vd(1,1))", "Re(Vd(1,2))",
      "Im(Vd(1,2))", "Re(Vd(2,0))", "Im(Vd(2,0))", "Re(Vd(2,1))", "Im(Vd(2,1))",
      "Re(Vd(2,2))", "Im(Vd(2,2))", "Re(Ud(0,0))", "Im(Ud(0,0))", "Re(Ud(0,1))",
      "Im(Ud(0,1))", "Re(Ud(0,2))", "Im(Ud(0,2))", "Re(Ud(1,0))", "Im(Ud(1,0))",
      "Re(Ud(1,1))", "Im(Ud(1,1))", "Re(Ud(1,2))", "Im(Ud(1,2))", "Re(Ud(2,0))",
      "Im(Ud(2,0))", "Re(Ud(2,1))", "Im(Ud(2,1))", "Re(Ud(2,2))", "Im(Ud(2,2))",
      "Re(Vu(0,0))", "Im(Vu(0,0))", "Re(Vu(0,1))", "Im(Vu(0,1))", "Re(Vu(0,2))",
      "Im(Vu(0,2))", "Re(Vu(1,0))", "Im(Vu(1,0))", "Re(Vu(1,1))", "Im(Vu(1,1))",
      "Re(Vu(1,2))", "Im(Vu(1,2))", "Re(Vu(2,0))", "Im(Vu(2,0))", "Re(Vu(2,1))",
      "Im(Vu(2,1))", "Re(Vu(2,2))", "Im(Vu(2,2))", "Re(Uu(0,0))", "Im(Uu(0,0))",
      "Re(Uu(0,1))", "Im(Uu(0,1))", "Re(Uu(0,2))", "Im(Uu(0,2))", "Re(Uu(1,0))",
      "Im(Uu(1,0))", "Re(Uu(1,1))", "Im(Uu(1,1))", "Re(Uu(1,2))", "Im(Uu(1,2))",
      "Re(Uu(2,0))", "Im(Uu(2,0))", "Re(Uu(2,1))", "Im(Uu(2,1))", "Re(Uu(2,2))",
      "Im(Uu(2,2))", "Re(Ve(0,0))", "Im(Ve(0,0))", "Re(Ve(0,1))", "Im(Ve(0,1))",
      "Re(Ve(0,2))", "Im(Ve(0,2))", "Re(Ve(1,0))", "Im(Ve(1,0))", "Re(Ve(1,1))",
      "Im(Ve(1,1))", "Re(Ve(1,2))", "Im(Ve(1,2))", "Re(Ve(2,0))", "Im(Ve(2,0))",
      "Re(Ve(2,1))", "Im(Ve(2,1))", "Re(Ve(2,2))", "Im(Ve(2,2))", "Re(Ue(0,0))",
      "Im(Ue(0,0))", "Re(Ue(0,1))", "Im(Ue(0,1))", "Re(Ue(0,2))", "Im(Ue(0,2))",
      "Re(Ue(1,0))", "Im(Ue(1,0))", "Re(Ue(1,1))", "Im(Ue(1,1))", "Re(Ue(1,2))",
      "Im(Ue(1,2))", "Re(Ue(2,0))", "Im(Ue(2,0))", "Re(Ue(2,1))", "Im(Ue(2,1))",
      "Re(Ue(2,2))", "Im(Ue(2,2))"};

   const std::string model_name = "SM";

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                " << model_name << '\n'
      << "Is a low-energy model:     "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model: "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Number of multiplets:      " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:      " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                ";
   for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                ";
   for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }
   ostr << '\n';
}

} // namespace standard_model_info

namespace standard_model {

const int Standard_model::numberOfParameters;

#define PHYSICAL(parameter) physical.parameter

#define HIGGS_2LOOP_CORRECTION_AT_AS       loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS       loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT       loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU   loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION            loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS    loop_corrections.higgs_at_as_as
#define HIGGS_3LOOP_CORRECTION_AB_AS_AS    loop_corrections.higgs_ab_as_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AS    loop_corrections.higgs_at_at_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AT    loop_corrections.higgs_at_at_at
#define HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS loop_corrections.higgs_at_as_as_as

Standard_model::Standard_model()
{
   set_number_of_parameters(numberOfParameters);
   set_thresholds(3);
}

Standard_model::Standard_model(double scale_, double loops_, double thresholds_
   , double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<
   double,3,3>& Yu_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<
   double,3,3>& Ye_, double mu2_, double v_)
   : Beta_function()
   , g1(g1_), g2(g2_), g3(g3_), Lambdax(Lambdax_), Yu(Yu_), Yd(Yd_), Ye(Ye_)
   , mu2(mu2_), v(v_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

void Standard_model::do_force_output(bool flag)
{
   force_output = flag;
}

bool Standard_model::do_force_output() const
{
   return force_output;
}

void Standard_model::set_ewsb_loop_order(int loop_order)
{
   ewsb_loop_order = loop_order;
}

void Standard_model::set_loop_corrections(const Loop_corrections& loop_corrections_)
{
   loop_corrections = loop_corrections_;
}

const Loop_corrections& Standard_model::get_loop_corrections() const
{
   return loop_corrections;
}

void Standard_model::set_threshold_corrections(const Threshold_corrections& tc)
{
   threshold_corrections = tc;
}

const Threshold_corrections& Standard_model::get_threshold_corrections() const
{
   return threshold_corrections;
}

int Standard_model::get_number_of_ewsb_iterations() const
{
   return static_cast<int>(std::abs(-log10(ewsb_iteration_precision) * 10));
}

int Standard_model::get_number_of_mass_iterations() const
{
   return static_cast<int>(std::abs(-log10(precision) * 10));
}

void Standard_model::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
}

void Standard_model::set_pole_mass_loop_order(int loop_order)
{
   pole_mass_loop_order = loop_order;
}

int Standard_model::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void Standard_model::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double Standard_model::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double Standard_model::get_precision() const
{
   return precision;
}

double Standard_model::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const Standard_model_physical& Standard_model::get_physical() const
{
   return physical;
}

Standard_model_physical& Standard_model::get_physical()
{
   return physical;
}

void Standard_model::set_physical(const Standard_model_physical& physical_)
{
   physical = physical_;
}

const Problems& Standard_model::get_problems() const
{
   return problems;
}

Problems& Standard_model::get_problems()
{
   return problems;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double, Standard_model::number_of_ewsb_equations, 1> Standard_model::tadpole_equations() const
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole(
      Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop());
   }

   return tadpole;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int Standard_model::solve_ewsb_iteratively()
{
   auto ewsb_stepper = [this](const EWSB_vector_t& ewsb_pars) -> EWSB_vector_t {
      set_mu2(ewsb_pars(0));
      if (get_ewsb_loop_order() > 0)
         calculate_DRbar_masses();
      return ewsb_step();
   };
   auto tadpole_stepper = [this](const EWSB_vector_t& ewsb_pars) -> EWSB_vector_t {
      set_mu2(ewsb_pars(0));
      if (get_ewsb_loop_order() > 0)
         calculate_DRbar_masses();
      return tadpole_equations();
   };

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative<number_of_ewsb_equations> >(ewsb_stepper, get_number_of_ewsb_iterations(), fixed_point_iterator::Convergence_tester_relative<number_of_ewsb_equations>(ewsb_iteration_precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, get_number_of_ewsb_iterations(), ewsb_iteration_precision, Root_finder<number_of_ewsb_equations>::GSLHybridS)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, get_number_of_ewsb_iterations(), ewsb_iteration_precision, Root_finder<number_of_ewsb_equations>::GSLBroyden))
   };

   const auto x_init(ewsb_initial_guess());

   VERBOSE_MSG("\t\tSolving EWSB equations ...\n"
               "\t\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (auto& solver: solvers) {
      VERBOSE_MSG("\t\t\tStarting EWSB iteration using " << solver->name());
      status = solve_ewsb_iteratively_with(solver.get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\t\t\t" << solver->name() << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\t\t\t" << solver->name() << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\t\t\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int Standard_model::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const Eigen::Matrix<double, number_of_ewsb_equations, 1>& x_init
)
{
   const int status = solver->solve(x_init);
   const auto solution = solver->get_solution();

   mu2 = solution(0);


   return status;
}

int Standard_model::solve_ewsb_iteratively(int loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const auto save_loop_order_raii = make_raii_save(ewsb_loop_order);
   ewsb_loop_order = loop_order;
   return solve_ewsb_iteratively();
}


int Standard_model::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_mu2 = mu2;

   mu2 = Re(0.5*Lambdax*Sqr(v));

   const bool is_finite = IsFinite(mu2);

   if (!is_finite) {
      mu2 = old_mu2;
      error = 1;
   }


   return error;
}

int Standard_model::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int Standard_model::solve_ewsb()
{
   VERBOSE_MSG("\t\tSolving Standard model EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

Eigen::Matrix<double, Standard_model::number_of_ewsb_equations, 1> Standard_model::ewsb_initial_guess()
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> x_init;

   x_init[0] = mu2;


   return x_init;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @return new set of EWSB output parameters
 */
Eigen::Matrix<double, Standard_model::number_of_ewsb_equations, 1> Standard_model::ewsb_step() const
{
   std::array<double, number_of_ewsb_equations> tadpole{};
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh_1loop());
   }

   double mu2;

   mu2 = Re((0.5*(Cube(v)*Lambdax - 2*tadpole[0]))/v);

   const bool is_finite = IsFinite(mu2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = mu2;


   return ewsb_parameters;
}

void Standard_model::print() const
{
   this->print(std::cerr);
}

void Standard_model::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "Standard model Q = " << get_scale() << "GeV \n"
           "========================================\n";
   ostr << "running parameters:\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "v = " << v << '\n';
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 */

double Standard_model::A0(double m) const
{
   return Loop_library::get().A0(m*m, Sqr(get_scale())).real();
}

double Standard_model::B0(double p, double m1, double m2) const
{
   return Loop_library::get().B0(p*p, m1*m1, m2*m2, Sqr(get_scale())).real();
}

double Standard_model::B1(double p, double m1, double m2) const
{
   return Loop_library::get().B1(p*p, m1*m1, m2*m2, Sqr(get_scale())).real();
}

double Standard_model::B00(double p, double m1, double m2) const
{
   return Loop_library::get().B00(p*p, m1*m1, m2*m2, Sqr(get_scale())).real();
}

double Standard_model::B22(double p, double m1, double m2) const
{
   const double scl2 = Sqr(get_scale());
   return (
      Loop_library::get().B00(p*p, m1*m1, m2*m2, scl2) - Loop_library::get().A0(m1*m1, scl2)/4.0 -
      Loop_library::get().A0(m2*m2, scl2)/4.0
   ).real();
}

double Standard_model::H0(double p, double m1, double m2) const
{
   const double scl2 = Sqr(get_scale());
   return 4.0*Loop_library::get().B00(p*p, m1*m1, m2*m2, scl2).real() + Standard_model::G0(p, m1, m2);
}

double Standard_model::F0(double p, double m1, double m2) const
{
   const double scl2 = Sqr(get_scale());
   return (
      Loop_library::get().A0(m1*m1, scl2) - 2.0*Loop_library::get().A0(m2*m2, scl2)
      - (2.0*p*p + 2.0*m1*m1 - m2*m2) * Loop_library::get().B0(p*p, m1*m1, m2*m2, scl2)
   ).real();
}

double Standard_model::G0(double p, double m1, double m2) const
{
   const double scl2 = Sqr(get_scale());
   return (
      (p*p - m1*m1 - m2*m2) * Loop_library::get().B0(p*p, m1*m1, m2*m2, scl2)
      - Loop_library::get().A0(m1*m1, scl2) - Loop_library::get().A0(m2*m2, scl2)
   ).real();
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void Standard_model::calculate_DRbar_masses()
{
   const auto save_mu2_raii = make_raii_save(mu2);

   solve_ewsb_tree_level();

   calculate_MVPVZ();
   calculate_MVWp();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_MAh();
   calculate_MVP();
   calculate_Mhh();
   calculate_MVZ();
   calculate_MFv();
   calculate_MHp();
   calculate_MVG();
}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void Standard_model::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 9u));

   tp.run_task([this] () { calculate_MVG_pole(); });
   tp.run_task([this] () { calculate_MFv_pole(); });
   tp.run_task([this] () { calculate_Mhh_pole(); });
   tp.run_task([this] () { calculate_MVP_pole(); });
   tp.run_task([this] () { calculate_MVZ_pole(); });
   tp.run_task([this] () { calculate_MFd_pole(); });
   tp.run_task([this] () { calculate_MFu_pole(); });
   tp.run_task([this] () { calculate_MFe_pole(); });
   tp.run_task([this] () { calculate_MVWp_pole(); });

#else

   calculate_MVG_pole();
   calculate_MFv_pole();
   calculate_Mhh_pole();
   calculate_MVP_pole();
   calculate_MVZ_pole();
   calculate_MFd_pole();
   calculate_MFu_pole();
   calculate_MFe_pole();
   calculate_MVWp_pole();

#endif
}

void Standard_model::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MHp) = MHp;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MVWp) = MVWp;
   PHYSICAL(MVPVZ) = MVPVZ;
   PHYSICAL(ZZ) = ZZ;

}

/**
 * Checks the pole masses for tachyons
 */
void Standard_model::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mhh) < 0.) problems.flag_pole_tachyon(standard_model_info::hh);

}

/**
 * calculates spectrum for model once the MSbar parameters at
 * at low energies are known
 */
void Standard_model::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   if (pole_mass_loop_order == 0) {
      copy_DRbar_masses_to_pole_masses();
   }

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void Standard_model::clear_running_parameters()
{
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   mu2 = 0.;
   v = 0.;
}

void Standard_model::clear_DRbar_parameters()
{
   MVG = 0.;
   MHp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MAh = 0.;
   Mhh = 0.;
   MVP = 0.;
   MVZ = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWp = 0.;
   MVPVZ = Eigen::Matrix<double,2,1>::Zero();
   ZZ = Eigen::Matrix<double,2,2>::Zero();
}

void Standard_model::clear_problems()
{
   problems.unflag_all_tachyons();
}

void Standard_model::clear()
{
   reset();
   clear_running_parameters();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

std::string Standard_model::name() const
{
   return "Standard model";
}

void Standard_model::initialise_from_input(const softsusy::QedQcd& qedqcd_)
{
   const double scale = get_scale();
   auto qedqcd = qedqcd_;

   // initial guess
   qedqcd.to(qedqcd.displayPoleMZ());
   initial_guess_for_parameters(qedqcd);
   run_to(qedqcd.displayPoleMZ());

   // determine Standard model parameters iteratively
   Standard_model old(*this);

   bool converged = false;

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mz_pole  = qedqcd.displayPoleMZ();

   while (!converged) {
      calculate_DRbar_masses();

      double mz_run = mz_pole;

      if (get_thresholds() && threshold_corrections.mz > 0)
         mz_run = calculate_MVZ_DRbar(mz_pole);

      double delta_alpha_em = 0.;
      double delta_alpha_s  = 0.;

      if (get_thresholds() && threshold_corrections.alpha_em > 0)
         delta_alpha_em = calculate_delta_alpha_em(alpha_em);

      if (get_thresholds() && threshold_corrections.alpha_s > 0)
         delta_alpha_s  = calculate_delta_alpha_s(alpha_s);

      const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
      const double alpha_s_drbar  = alpha_s / (1.0 - delta_alpha_s);
      const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);
      const double theta_w_drbar  = calculate_theta_w(qedqcd, alpha_em_drbar);

      v = Re((2 * mz_run) / Sqrt(0.6 * Sqr(g1) + Sqr(g2)));

      calculate_Yu_DRbar(qedqcd);
      calculate_Yd_DRbar(qedqcd);
      calculate_Ye_DRbar(qedqcd);
      calculate_Lambdax_DRbar();

      solve_ewsb();

      g1 = 1.2909944487358056 * e_drbar * Sec(theta_w_drbar);
      g2 = e_drbar * Csc(theta_w_drbar);
      g3 = 3.5449077018110318 * Sqrt(alpha_s_drbar);

      if (IsFinite(g1)) {
         problems.unflag_non_perturbative_parameter(standard_model_info::g1);
      } else {
         problems.flag_non_perturbative_parameter(
            standard_model_info::g1, g1, get_scale());
         g1 = Electroweak_constants::g1;
      }

      if (IsFinite(g2)) {
         problems.unflag_non_perturbative_parameter(standard_model_info::g2);
      } else {
         problems.flag_non_perturbative_parameter(
            standard_model_info::g2, g2, get_scale());
         g2 = Electroweak_constants::g2;
      }

      if (get_thresholds() && threshold_corrections.sin_theta_w > 0)
         qedqcd.setPoleMW(recalculate_mw_pole(qedqcd.displayPoleMW()));

      converged = check_convergence(old);
      old = *this;
   }

   // run all SM parameters to the desired scale
   if (scale > 0.)
      run_to(scale);
}

void Standard_model::initial_guess_for_parameters(const softsusy::QedQcd& qedqcd)
{
   const double MH = input.get(Physical_input::mh_pole);
   const double mtpole = qedqcd.displayPoleMt();

   const double mu_guess = qedqcd.displayMass(softsusy::mUp);
   const double mc_guess = qedqcd.displayMass(softsusy::mCharm);
   const double mt_guess = get_thresholds() > 0 && threshold_corrections.mt > 0 ?
      qedqcd.displayMass(softsusy::mTop) - 30.0 :
      qedqcd.displayPoleMt();
   const double md_guess = qedqcd.displayMass(softsusy::mDown);
   const double ms_guess = qedqcd.displayMass(softsusy::mStrange);
   const double mb_guess = qedqcd.displayMass(softsusy::mBottom);
   const double me_guess = get_thresholds() > 0 ?
      qedqcd.displayMass(softsusy::mElectron) :
      qedqcd.displayPoleMel();
   const double mm_guess = get_thresholds() > 0 ?
      qedqcd.displayMass(softsusy::mMuon) :
      qedqcd.displayPoleMmuon();
   const double mtau_guess = qedqcd.displayMass(softsusy::mTau);

   // guess gauge couplings at mt
   const auto alpha_sm(qedqcd.guess_alpha_SM5(mtpole));

   g1 = Sqrt(4. * Pi * alpha_sm(0));
   g2 = Sqrt(4. * Pi * alpha_sm(1));
   g3 = Sqrt(4. * Pi * alpha_sm(2));

   set_scale(mtpole);

   v = Re(Electroweak_constants::vev);

   Eigen::Matrix<double,3,3> upQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;
   Yu = -((1.4142135623730951*upQuarksDRbar) / v).transpose();

   Eigen::Matrix<double,3,3> downQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;
   Yd = ((1.4142135623730951*downQuarksDRbar)/v).transpose();

   Eigen::Matrix<double,3,3> downLeptonsDRbar(Eigen::Matrix<double,3,3>::Zero());
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;
   Ye = ((1.4142135623730951*downLeptonsDRbar)/v).transpose();

   Lambdax = Sqr(MH) / Sqr(v);

   solve_ewsb_tree_level();
}

double Standard_model::calculate_delta_alpha_em(double alphaEm) const
{
   const double delta_alpha_em_SM = -0.28294212105225836*alphaEm*
      FiniteLog(Abs(MFu(2) / get_scale()));

   return delta_alpha_em_SM;
}

double Standard_model::calculate_delta_alpha_s(double alphaS) const
{
   const double delta_alpha_s_1loop = -0.1061032953945969*alphaS*
      FiniteLog(Abs(MFu(2) / get_scale()));

   double delta_alpha_s_2loop = 0.;
   double delta_alpha_s_3loop = 0.;
   double delta_alpha_s_4loop = 0.;

   if (get_thresholds() > 1 && threshold_corrections.alpha_s > 1) {
      sm_fourloop_as::Parameters pars;
      pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
      pars.mt   = MFu(2);
      pars.Q    = get_scale();

      const auto das_1L = sm_fourloop_as::delta_alpha_s_1loop_as(pars);
      const auto das_2L = sm_fourloop_as::delta_alpha_s_2loop_as_as(pars);

      delta_alpha_s_2loop = das_2L - Sqr(das_1L);
   }

   if (get_thresholds() > 2 && get_threshold_corrections().alpha_s > 2) {
      sm_fourloop_as::Parameters pars;
      pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
      pars.mt   = MFu(2);
      pars.Q    = get_scale();

      const auto das_1L = sm_fourloop_as::delta_alpha_s_1loop_as(pars);
      const auto das_2L = sm_fourloop_as::delta_alpha_s_2loop_as_as(pars);
      const auto das_3L = sm_fourloop_as::delta_alpha_s_3loop_as_as_as(pars);

      delta_alpha_s_3loop = das_3L + Power3(das_1L) - 2. * das_1L * das_2L;
   }

   if (get_thresholds() > 3 && get_threshold_corrections().alpha_s > 3) {
      sm_fourloop_as::Parameters pars;
      pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
      pars.mt   = MFu(2);
      pars.Q    = get_scale();

      const auto das_1L = sm_fourloop_as::delta_alpha_s_1loop_as(pars);
      const auto das_2L = sm_fourloop_as::delta_alpha_s_2loop_as_as(pars);
      const auto das_3L = sm_fourloop_as::delta_alpha_s_3loop_as_as_as(pars);
      const auto das_4L = sm_fourloop_as::delta_alpha_s_4loop_as_as_as_as(pars);

      delta_alpha_s_4loop = das_4L - 2. * das_1L * das_3L - Power2(das_2L) + 3.
         * Power2(das_1L) * das_2L - Power4(das_1L);
   }

   return delta_alpha_s_1loop + delta_alpha_s_2loop + delta_alpha_s_3loop
      + delta_alpha_s_4loop;

}

double Standard_model::calculate_theta_w(const softsusy::QedQcd& qedqcd, double alpha_em_drbar)
{
   double theta_w = 0.;

   using namespace weinberg_angle;

   const double scale               = get_scale();
   const double mw_pole             = qedqcd.displayPoleMW();
   const double mz_pole             = qedqcd.displayPoleMZ();
   const double mt_pole             = qedqcd.displayPoleMt();
   const double mt_drbar            = MFu(2);
   const double mb_drbar            = MFd(2);
   const double mh_drbar            = Mhh;
   const double gY                  = g1 * 0.7745966692414834;
   const double pizztMZ             = Re(self_energy_VZ_1loop(mz_pole));
   const double piwwt0              = Re(self_energy_VWp_1loop(0.));
   const double self_energy_w_at_mw = Re(self_energy_VWp_1loop(mw_pole));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY       = gY;
   se_data.g2       = g2;

   double pizztMZ_corrected = pizztMZ;
   double piwwtMW_corrected = self_energy_w_at_mw;
   double piwwt0_corrected  = piwwt0;

   if (get_thresholds() > 1 && threshold_corrections.sin_theta_w > 1) {
      pizztMZ_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole,
            se_data);
      piwwtMW_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(
            self_energy_w_at_mw, mw_pole, se_data);
      piwwt0_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0.,
            se_data);
   }

   Weinberg_angle::Data data;
   data.scale               = scale;
   data.alpha_em_drbar      = alpha_em_drbar;
   data.fermi_contant       = qedqcd.displayFermiConstant();
   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
   data.mw_pole             = mw_pole;
   data.mz_pole             = mz_pole;
   data.mt_pole             = mt_pole;
   data.mh_drbar            = mh_drbar;
   data.gY                  = gY;
   data.g2                  = g2;
   data.g3                  = g3;

   Weinberg_angle weinberg;
   weinberg.disable_susy_contributions();
   weinberg.set_number_of_loops(threshold_corrections.sin_theta_w);
   weinberg.set_data(data);

   const int error = weinberg.calculate();

   theta_w = ArcSin(weinberg.get_sin_theta());

   if (error)
      problems.flag_no_sinThetaW_convergence();
   else
      problems.unflag_no_sinThetaW_convergence();

   return theta_w;
}

void Standard_model::calculate_Yu_DRbar(const softsusy::QedQcd& qedqcd)
{
   Eigen::Matrix<double,3,3> upQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   upQuarksDRbar(0,0)      = qedqcd.displayMass(softsusy::mUp);
   upQuarksDRbar(1,1)      = qedqcd.displayMass(softsusy::mCharm);
   upQuarksDRbar(2,2)      = qedqcd.displayPoleMt();

   if (get_thresholds() && threshold_corrections.mt > 0)
      upQuarksDRbar(2,2) = calculate_MFu_DRbar(qedqcd.displayPoleMt(), 2);

   Yu = -((1.4142135623730951*upQuarksDRbar)/v).transpose();
}

void Standard_model::calculate_Yd_DRbar(const softsusy::QedQcd& qedqcd)
{
   Eigen::Matrix<double,3,3> downQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   downQuarksDRbar(0,0)   = qedqcd.displayMass(softsusy::mDown);
   downQuarksDRbar(1,1)   = qedqcd.displayMass(softsusy::mStrange);
   downQuarksDRbar(2,2)   = qedqcd.displayMass(softsusy::mBottom);

   if (get_thresholds() && threshold_corrections.mb > 0)
      downQuarksDRbar(2,2) = calculate_MFd_DRbar(qedqcd.displayMass(softsusy::mBottom), 2);

   Yd = ((1.4142135623730951*downQuarksDRbar)/v).transpose();
}

void Standard_model::calculate_Ye_DRbar(const softsusy::QedQcd& qedqcd)
{
   Eigen::Matrix<double,3,3> downLeptonsDRbar(Eigen::Matrix<double,3,3>::Zero());
   downLeptonsDRbar(0,0) = qedqcd.displayPoleMel();
   downLeptonsDRbar(1,1) = qedqcd.displayPoleMmuon();
   downLeptonsDRbar(2,2) = qedqcd.displayPoleMtau();

   if (get_thresholds()) {
      downLeptonsDRbar(0,0) = calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mElectron), 0);
      downLeptonsDRbar(1,1) = calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mMuon), 1);
   }

   if (get_thresholds() && threshold_corrections.mtau > 0)
      downLeptonsDRbar(2,2) = calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mTau), 2);

   Ye = ((1.4142135623730951*downLeptonsDRbar)/v).transpose();
}

void Standard_model::calculate_Lambdax_DRbar()
{
   double higgsDRbar = input.get(Physical_input::mh_pole);

   if (get_thresholds() && threshold_corrections.mh > 0)
      higgsDRbar = calculate_Mhh_DRbar(higgsDRbar);

   Lambdax = Sqr(higgsDRbar) / Sqr(v);
}

double Standard_model::recalculate_mw_pole(double mw_pole)
{
   calculate_MVWp();

   const double mw_drbar    = MVWp;
   const double mw_pole_sqr = Sqr(mw_drbar) - Re(self_energy_VWp_1loop(mw_pole));

   if (mw_pole_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VWp);

   return AbsSqrt(mw_pole_sqr);
}

double Standard_model::max_rel_diff(const Standard_model& old) const
{
   std::array<double, 12> diff{};

   diff[0] = MaxRelDiff(old.Mhh, Mhh);
   diff[1] = MaxRelDiff(old.MVZ, MVZ);
   for (int i = 0; i < 3; ++i)
      diff[i + 2] = MaxRelDiff(old.MFd(i), MFd(i));
   for (int i = 0; i < 3; ++i)
      diff[i + 5] = MaxRelDiff(old.MFu(i), MFu(i));
   for (int i = 0; i < 3; ++i)
      diff[i + 8] = MaxRelDiff(old.MFe(i), MFe(i));
   diff[11] = MaxRelDiff(old.MVWp, MVWp);

   return *std::max_element(diff.cbegin(), diff.cend());
}

bool Standard_model::check_convergence(const Standard_model& old) const
{
   return max_rel_diff(old) < precision;
}

Eigen::ArrayXd Standard_model::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g1;
   pars(1) = g2;
   pars(2) = g3;
   pars(3) = Lambdax;
   pars(4) = Yu(0,0);
   pars(5) = Yu(0,1);
   pars(6) = Yu(0,2);
   pars(7) = Yu(1,0);
   pars(8) = Yu(1,1);
   pars(9) = Yu(1,2);
   pars(10) = Yu(2,0);
   pars(11) = Yu(2,1);
   pars(12) = Yu(2,2);
   pars(13) = Yd(0,0);
   pars(14) = Yd(0,1);
   pars(15) = Yd(0,2);
   pars(16) = Yd(1,0);
   pars(17) = Yd(1,1);
   pars(18) = Yd(1,2);
   pars(19) = Yd(2,0);
   pars(20) = Yd(2,1);
   pars(21) = Yd(2,2);
   pars(22) = Ye(0,0);
   pars(23) = Ye(0,1);
   pars(24) = Ye(0,2);
   pars(25) = Ye(1,0);
   pars(26) = Ye(1,1);
   pars(27) = Ye(1,2);
   pars(28) = Ye(2,0);
   pars(29) = Ye(2,1);
   pars(30) = Ye(2,2);
   pars(31) = mu2;
   pars(32) = v;

   return pars;
}

void Standard_model::set(const Eigen::ArrayXd& pars)
{
   g1 = pars(0);
   g2 = pars(1);
   g3 = pars(2);
   Lambdax = pars(3);
   Yu(0,0) = pars(4);
   Yu(0,1) = pars(5);
   Yu(0,2) = pars(6);
   Yu(1,0) = pars(7);
   Yu(1,1) = pars(8);
   Yu(1,2) = pars(9);
   Yu(2,0) = pars(10);
   Yu(2,1) = pars(11);
   Yu(2,2) = pars(12);
   Yd(0,0) = pars(13);
   Yd(0,1) = pars(14);
   Yd(0,2) = pars(15);
   Yd(1,0) = pars(16);
   Yd(1,1) = pars(17);
   Yd(1,2) = pars(18);
   Yd(2,0) = pars(19);
   Yd(2,1) = pars(20);
   Yd(2,2) = pars(21);
   Ye(0,0) = pars(22);
   Ye(0,1) = pars(23);
   Ye(0,2) = pars(24);
   Ye(1,0) = pars(25);
   Ye(1,1) = pars(26);
   Ye(1,2) = pars(27);
   Ye(2,0) = pars(28);
   Ye(2,1) = pars(29);
   Ye(2,2) = pars(30);
   mu2 = pars(31);
   v = pars(32);
}

void Standard_model::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   Beta_function::run_to(scale, eps);
}

Eigen::ArrayXd Standard_model::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

Standard_model Standard_model::calc_beta(int loops) const
{
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_Lambdax = 0.;
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   double beta_mu2 = 0.;
   double beta_v = 0.;

   if (loops > 0) {
      Beta_traces traces;
      calc_beta_traces(traces);

      beta_g1 += calc_beta_g1_one_loop(traces);
      beta_g2 += calc_beta_g2_one_loop(traces);
      beta_g3 += calc_beta_g3_one_loop(traces);
      beta_Lambdax += calc_beta_Lambdax_one_loop(traces);
      beta_Yu += calc_beta_Yu_one_loop(traces);
      beta_Yd += calc_beta_Yd_one_loop(traces);
      beta_Ye += calc_beta_Ye_one_loop(traces);
      beta_mu2 += calc_beta_mu2_one_loop(traces);
      beta_v += calc_beta_v_one_loop(traces);

      if (loops > 1) {
         beta_g1 += calc_beta_g1_two_loop(traces);
         beta_g2 += calc_beta_g2_two_loop(traces);
         beta_g3 += calc_beta_g3_two_loop(traces);
         beta_Lambdax += calc_beta_Lambdax_two_loop(traces);
         beta_Yu += calc_beta_Yu_two_loop(traces);
         beta_Yd += calc_beta_Yd_two_loop(traces);
         beta_Ye += calc_beta_Ye_two_loop(traces);
         beta_mu2 += calc_beta_mu2_two_loop(traces);
         beta_v += calc_beta_v_two_loop(traces);

         if (loops > 2) {
            beta_g1 += calc_beta_g1_three_loop(traces);
            beta_g2 += calc_beta_g2_three_loop(traces);
            beta_g3 += calc_beta_g3_three_loop(traces);
            beta_Lambdax += calc_beta_Lambdax_three_loop(traces);
            beta_Yu += calc_beta_Yu_three_loop(traces);
            beta_Yd += calc_beta_Yd_three_loop(traces);
            beta_Ye += calc_beta_Ye_three_loop(traces);
            beta_mu2 += calc_beta_mu2_three_loop(traces);

            if (loops > 3) {
               beta_g3 += calc_beta_g3_four_loop(traces);
               beta_Lambdax += calc_beta_Lambdax_four_loop(traces);
               beta_Yu += calc_beta_Yu_four_loop(traces);

               if (loops > 4) {
                  beta_g3 += calc_beta_g3_five_loop(traces);
               }
            }
         }
      }
   }

   return Standard_model(get_scale(), loops, get_thresholds(),
                         beta_g1, beta_g2, beta_g3, beta_Lambdax, beta_Yu, beta_Yd, beta_Ye,
                         beta_mu2, beta_v);
}

Standard_model Standard_model::calc_beta() const
{
   return calc_beta(get_loops());
}

void Standard_model::calc_beta_traces(Beta_traces& traces) const
{
   const int loops = get_loops();

   if (loops > 0) {
      traces.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      traces.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      traces.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      traces.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      traces.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      traces.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());
   }

   if (loops > 1) {
      traces.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      traces.traceYdAdjYdYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint(
         )*Yd*Yd.adjoint()).trace());
      traces.traceYdAdjYdYdAdjYuYuAdjYd = Re((Yd*Yd.adjoint()*Yd*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace());
      traces.traceYdAdjYuYuAdjYdYdAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint(
         )*Yd*Yd.adjoint()).trace());
      traces.traceYdAdjYuYuAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace());
      traces.traceYeAdjYeYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint(
         )*Ye*Ye.adjoint()).trace());
      traces.traceYuAdjYuYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint(
         )*Yu*Yu.adjoint()).trace());
   }
}


double Standard_model::calc_beta_g1_one_loop(const Beta_traces&) const
{
   double beta_g1;

   beta_g1 = Re(4.1*oneOver16PiSqr*Cube(g1));

   return beta_g1;
}

double Standard_model::calc_beta_g1_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_g1;

   beta_g1 = Re(0.02*twoLoop*Cube(g1)*(199*Sqr(g1) + 5*(-5*traceYdAdjYd -
      15*traceYeAdjYe - 17*traceYuAdjYu + 27*Sqr(g2) + 88*Sqr(g3))));

   return beta_g1;
}

double Standard_model::calc_beta_g1_three_loop(const Beta_traces&) const
{
   double beta_g1;

   beta_g1 = Re(-0.000041666666666666665*threeLoop*Cube(g1)*(388613*Quad(
      g1) - 10*Sqr(g1)*(648*Lambdax + 1845*Sqr(g2) - 4384*Sqr(g3) - 8481*Sqr(Yu
      (2,2))) - 75*(3945*Quad(g2) - 6*Sqr(g2)*(-24*Lambdax + 32*Sqr(g3) + 785*
      Sqr(Yu(2,2))) + 4*(4752*Quad(g3) + 945*Quad(Yu(2,2)) - 36*Sqr(Lambdax) -
      464*Sqr(g3)*Sqr(Yu(2,2))))));

   return beta_g1;
}

double Standard_model::calc_beta_g2_one_loop(const Beta_traces&) const
{
   double beta_g2;

   beta_g2 = Re(-3.1666666666666665*oneOver16PiSqr*Cube(g2));

   return beta_g2;
}

double Standard_model::calc_beta_g2_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_g2;

   beta_g2 = Re(0.03333333333333333*twoLoop*Cube(g2)*(27*Sqr(g1) + 5*(-3*
      (3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu) + 35*Sqr(g2) + 72*Sqr(g3
      ))));

   return beta_g2;
}

double Standard_model::calc_beta_g2_three_loop(const Beta_traces&) const
{
   double beta_g2;

   beta_g2 = Re(0.000023148148148148147*threeLoop*Cube(g2)*(-151119*Quad(
      g1) + 270*Sqr(g1)*(24*Lambdax + 873*Sqr(g2) - 32*Sqr(g3) - 593*Sqr(Yu(2,2
      ))) + 25*(324953*Quad(g2) + 162*Sqr(g2)*(8*Lambdax + 416*Sqr(g3) - 243*
      Sqr(Yu(2,2))) + 108*(1296*Quad(g3) + 147*Quad(Yu(2,2)) - 12*Sqr(Lambdax)
      - 112*Sqr(g3)*Sqr(Yu(2,2))))));

   return beta_g2;
}

double Standard_model::calc_beta_g3_one_loop(const Beta_traces&) const
{
   double beta_g3;

   beta_g3 = Re(-7*oneOver16PiSqr*Cube(g3));

   return beta_g3;
}

double Standard_model::calc_beta_g3_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_g3;

   beta_g3 = Re(-0.1*twoLoop*Cube(g3)*(-11*Sqr(g1) + 5*(-9*Sqr(g2) + 4*(
      traceYdAdjYd + traceYuAdjYu + 13*Sqr(g3)))));

   return beta_g3;
}

double Standard_model::calc_beta_g3_three_loop(const Beta_traces&) const
{
   double beta_g3;

   beta_g3 = Re(0.008333333333333333*threeLoop*Cube(g3)*(-523*Quad(g1) +
      Sqr(g1)*(-9*Sqr(g2) + 616*Sqr(g3) - 303*Sqr(Yu(2,2))) + 15*(109*Quad(g2)
      + 3*Sqr(g2)*(56*Sqr(g3) - 31*Sqr(Yu(2,2))) + 20*(13*Quad(g3) + 6*Quad(Yu(
      2,2)) - 16*Sqr(g3)*Sqr(Yu(2,2))))));

   return beta_g3;
}

double Standard_model::calc_beta_g3_four_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   double beta_g3;

   beta_g3 = Re(-2472.2837425797156*Cube(g3)*Quad(oneOver16PiSqr)*(1.*
      Power6(g3) + 0.04569149546770327*Power6(Yu(2,2)) + 0.0060672647486441755*
      Lambdax*Quad(Yu(2,2)) - 0.12603833934188147*Quad(Yu(2,2))*Sqr(g3) +
      0.16927060578749137*Quad(g3)*Sqr(Yu(2,2)) - 0.003640358849186505*Sqr(
      Lambdax)*Sqr(Yu(2,2))));

   return beta_g3;
}

double Standard_model::calc_beta_g3_five_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   double beta_g3;

   beta_g3 = Re(271.4283824198132*Power11(g3)*Power5(oneOver16PiSqr));

   return beta_g3;
}

double Standard_model::calc_beta_Lambdax_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(0.27*Quad(g1) + 2.25*Quad(g2) - 9*
      Lambdax*Sqr(g2) + 0.9*Sqr(g1)*(-2*Lambdax + Sqr(g2)) + 4*(-3*
      traceYdAdjYdYdAdjYd - traceYeAdjYeYeAdjYe - 3*traceYuAdjYuYuAdjYu + 3*
      traceYdAdjYd*Lambdax + traceYeAdjYe*Lambdax + 3*traceYuAdjYu*Lambdax + 3*
      Sqr(Lambdax))));

   return beta_Lambdax;
}

double Standard_model::calc_beta_Lambdax_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      traces.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      traces.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      traces.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      traces.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      traces.traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      traces.traceYuAdjYuYuAdjYuYuAdjYu;

   double beta_Lambdax;

   beta_Lambdax = Re(twoLoop*(60*traceYdAdjYdYdAdjYdYdAdjYd - 24*
      traceYdAdjYdYdAdjYuYuAdjYd + 12*traceYdAdjYuYuAdjYdYdAdjYd - 12*
      traceYdAdjYuYuAdjYuYuAdjYd + 20*traceYeAdjYeYeAdjYeYeAdjYe + 60*
      traceYuAdjYuYuAdjYuYuAdjYu - 78*Cube(Lambdax) - 3*traceYdAdjYdYdAdjYd*
      Lambdax - 42*traceYdAdjYuYuAdjYd*Lambdax - traceYeAdjYeYeAdjYe*Lambdax -
      3*traceYuAdjYuYuAdjYu*Lambdax - 3.411*Power6(g1) + 38.125*Power6(g2) -
      0.125*(36*traceYdAdjYd + 12*traceYeAdjYe + 36*traceYuAdjYu + 73*Lambdax)*
      Quad(g2) + 1.5*Lambdax*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*
      traceYuAdjYu + 36*Lambdax)*Sqr(g2) - 0.015*Quad(g1)*(-60*traceYdAdjYd +
      300*traceYeAdjYe + 228*traceYuAdjYu - 629*Lambdax + 559*Sqr(g2)) - 64*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 64*traceYuAdjYuYuAdjYu*Sqr(g3) + 80*
      traceYdAdjYd*Lambdax*Sqr(g3) + 80*traceYuAdjYu*Lambdax*Sqr(g3) - 72*
      traceYdAdjYd*Sqr(Lambdax) - 24*traceYeAdjYe*Sqr(Lambdax) - 72*
      traceYuAdjYu*Sqr(Lambdax) - 0.025*Sqr(g1)*(289*Quad(g2) - 6*(36*
      traceYdAdjYd + 44*traceYeAdjYe + 84*traceYuAdjYu + 39*Lambdax)*Sqr(g2) -
      4*(16*traceYdAdjYdYdAdjYd - 48*traceYeAdjYeYeAdjYe - 32*
      traceYuAdjYuYuAdjYu + 25*traceYdAdjYd*Lambdax + 75*traceYeAdjYe*Lambdax +
      85*traceYuAdjYu*Lambdax + 108*Sqr(Lambdax)))));

   return beta_Lambdax;
}

double Standard_model::calc_beta_Lambdax_three_loop(const Beta_traces&) const
{
   double beta_Lambdax;

   const double beta_Lambdax_1 = Re(-6.031600412410262*threeLoop*(-
      9.339962685289576*Lambdax*Power6(g1) - 286.9827742168705*Lambdax*Power6(
      g2) + 74.07064511965402*Lambdax*Power6(Yd(2,2)) + 77.74412745450469*
      Lambdax*Power6(Ye(2,2)) + 1.*Power8(g1) + 75.66229767543491*Power8(g2) +
      161.25010404433053*Power8(Yd(2,2)) + 21.420306803784438*Power8(Ye(2,2)) -
      4.310550825785817*Quad(g1)*Quad(g2) - 249.16051734798086*Quad(Lambdax) +
      11.266458319444492*Quad(g1)*Quad(Yd(2,2)) - 10.53370494294467*Quad(g2)*
      Quad(Yd(2,2)) + 33.291890929823126*Quad(g3)*Quad(Yd(2,2)) -
      21.23162373590703*Quad(g1)*Quad(Ye(2,2)) - 3.262544766538025*Quad(g2)*
      Quad(Ye(2,2)) + 47.74852117315802*Quad(Yd(2,2))*Quad(Ye(2,2)) +
      6.423699789119175*Cube(Lambdax)*Sqr(g1) + 25.12726423496584*Power6(g2)*
      Sqr(g1) - 9.226669872878011*Power6(Yd(2,2))*Sqr(g1) - 17.462679807754416*
      Power6(Ye(2,2))*Sqr(g1) - 26.40696244250939*Lambdax*Quad(g2)*Sqr(g1) -
      11.990947297383288*Lambdax*Quad(Yd(2,2))*Sqr(g1) + 15.37217037999308*
      Lambdax*Quad(Ye(2,2))*Sqr(g1) + 32.118498945595874*Cube(Lambdax)*Sqr(g2)
      + 1.0231901443099605*Power6(g1)*Sqr(g2) - 49.16652864612754*Power6(Yd(2,2
      ))*Sqr(g2) - 16.388842882042518*Power6(Ye(2,2))*Sqr(g2) -
      20.47661959122464*Lambdax*Quad(g1)*Sqr(g2) + 1.8137171911699768*Lambdax*
      Quad(Yd(2,2))*Sqr(g2) - 2.3797101762657173*Lambdax*Quad(Ye(2,2))*Sqr(g2)
      + 25.306618248571045*Quad(Yd(2,2))*Sqr(g1)*Sqr(g2) + 22.872470753746928*
      Quad(Ye(2,2))*Sqr(g1)*Sqr(g2) - 0.4397999591948324*Power6(g1)*Sqr(g3) -
      4.997726809032187*Power6(g2)*Sqr(g3) - 166.12085657591095*Power6(Yd(2,2))
      *Sqr(g3) + 2.7790841962446198*Lambdax*Quad(g1)*Sqr(g3) +
      18.948301338031502*Lambdax*Quad(g2)*Sqr(g3) + 219.79763285742905*Lambdax*
      Quad(Yd(2,2))*Sqr(g3) - 0.9995453618064374*Quad(g2)*Sqr(g1)*Sqr(g3) -
      7.513283589518858*Quad(Yd(2,2))*Sqr(g1)*Sqr(g3) - 0.7329999319913874*Quad
      (g1)*Sqr(g2)*Sqr(g3) - 8.852950966952918*Quad(Yd(2,2))*Sqr(g2)*Sqr(g3) +
      30.759991811525552*Quad(g1)*Sqr(Lambdax) + 131.0233001667084*Quad(g2)*Sqr
      (Lambdax) - 293.1651465422671*Quad(Yd(2,2))*Sqr(Lambdax) -
      109.65884580737854*Quad(Ye(2,2))*Sqr(Lambdax) + 52.49683489237998*Sqr(g1)
      *Sqr(g2)*Sqr(Lambdax) - 72.36885240306762*Cube(Lambdax)*Sqr(Yd(2,2)) -
      3.23203546695533*Power6(g1)*Sqr(Yd(2,2)) - 41.447585739271474*Power6(g2)*
      Sqr(Yd(2,2)) + 24.620331229909603*Power6(Ye(2,2))*Sqr(Yd(2,2)) +
      11.460062483294552*Lambdax*Quad(g1)*Sqr(Yd(2,2)) + 105.99644858660962*
      Lambdax*Quad(g2)*Sqr(Yd(2,2)) - 118.365909843938*Lambdax*Quad(g3)*Sqr(Yd(
      2,2)) - 79.58086862193002*Lambdax*Quad(Ye(2,2))*Sqr(Yd(2,2)) -
      13.784493073708608*Quad(g2)*Sqr(g1)*Sqr(Yd(2,2)) - 8.381888765119909*Quad
      (g1)*Sqr(g2)*Sqr(Yd(2,2)) + 15.837899406921814*Lambdax*Sqr(g1)*Sqr(g2)*
      Sqr(Yd(2,2)) - 1.6285371043301784*Quad(g1)*Sqr(g3)*Sqr(Yd(2,2)) -
      10.91844691535382*Quad(g2)*Sqr(g3)*Sqr(Yd(2,2)) + 1.3873490852923183*
      Lambdax*Sqr(g1)*Sqr(g3)*Sqr(Yd(2,2)) - 5.021649329193744*Lambdax*Sqr(g2)*
      Sqr(g3)*Sqr(Yd(2,2)) - 5.958996138580129*Sqr(g1)*Sqr(g2)*Sqr(g3)*Sqr(Yd(2
      ,2)) + 12.5881938544474*Sqr(g1)*Sqr(Lambdax)*Sqr(Yd(2,2)) +
      59.60915139954888*Sqr(g2)*Sqr(Lambdax)*Sqr(Yd(2,2)) - 26.654542981504992*
      Sqr(g3)*Sqr(Lambdax)*Sqr(Yd(2,2)) - 24.12295080102254*Cube(Lambdax)*Sqr(
      Ye(2,2)) - 5.492026720736577*Power6(g1)*Sqr(Ye(2,2)) - 13.81586191309049*
      Power6(g2)*Sqr(Ye(2,2)) + 24.620331229909603*Power6(Yd(2,2))*Sqr(Ye(2,2))
      + 17.969663482497715*Lambdax*Quad(g1)*Sqr(Ye(2,2)) + 35.33214952886986*
      Lambdax*Quad(g2)*Sqr(Ye(2,2)) - 79.58086862193002*Lambdax*Quad(Yd(2,2))*
      Sqr(Ye(2,2)) - 1.5347452609060228*Quad(g2)*Sqr(g1)*Sqr(Ye(2,2)) -
      3.459017345205606*Quad(g1)*Sqr(g2)*Sqr(Ye(2,2)) - 6.687860766560821*
      Lambdax*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2)) + 1.9748527030237164*Sqr(g1)*Sqr(
      Lambdax)*Sqr(Ye(2,2)) + 19.86971713318296*Sqr(g2)*Sqr(Lambdax)*Sqr(Ye(2,2
      )) - 0.4078519516873914*Quad(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      0.746070643330594*Quad(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 1.7905695439934257
      *Lambdax*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 8.952847719967128*Lambdax*
      Sqr(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 0.4973804288870627*Sqr(g1)*Sqr(g2)*
      Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 35.81139087986851*Sqr(Lambdax)*Sqr(Yd(2,2))*
      Sqr(Ye(2,2)) - 7.372368699637704*Power6(g1)*Sqr(Yu(2,2)) -
      41.447585739271474*Power6(g2)*Sqr(Yu(2,2)) + 24.822556597226036*Lambdax*
      Quad(g1)*Sqr(Yu(2,2)) + 105.99644858660962*Lambdax*Quad(g2)*Sqr(Yu(2,2))
      - 118.365909843938*Lambdax*Quad(g3)*Sqr(Yu(2,2)) - 8.648442082445648*Quad
      (g2)*Sqr(g1)*Sqr(Yu(2,2)) - 7.047713089694501*Quad(g1)*Sqr(g2)*Sqr(Yu(2,2
      )) - 1.8617589799206342*Lambdax*Sqr(g1)*Sqr(g2)*Sqr(Yu(2,2)) -
      0.6735666808670179*Quad(g1)*Sqr(g3)*Sqr(Yu(2,2)) - 10.91844691535382*Quad
      (g2)*Sqr(g3)*Sqr(Yu(2,2)) - 5.787687768100881*Lambdax*Sqr(g1)*Sqr(g3)*Sqr
      (Yu(2,2)) - 5.021649329193744*Lambdax*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) -
      7.55061351101873*Sqr(g1)*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) +
      10.589103130834523*Sqr(g1)*Sqr(Lambdax)*Sqr(Yu(2,2)) + 59.60915139954888*
      Sqr(g2)*Sqr(Lambdax)*Sqr(Yu(2,2)) - 26.654542981504992*Sqr(g3)*Sqr(
      Lambdax)*Sqr(Yu(2,2))));
   const double beta_Lambdax_2 = Re(436.5*threeLoop*(-1.023515541010766*Lambdax
      *Power6(Yu(2,2)) - 2.228169974925502*Power8(Yu(2,2)) +
      0.14614411445535522*Quad(g1)*Quad(Yu(2,2)) + 0.14555578253854087*Quad(g2)
      *Quad(Yu(2,2)) - 0.4600306599364001*Quad(g3)*Quad(Yu(2,2)) +
      0.7931097093011755*Quad(Yd(2,2))*Quad(Yu(2,2)) - 0.6597938144329897*Quad(
      Ye(2,2))*Quad(Yu(2,2)) + 0.31092536460231723*Power6(Yu(2,2))*Sqr(g1) -
      0.09628707405242558*Lambdax*Quad(Yu(2,2))*Sqr(g1) + 0.6793879827233994*
      Power6(Yu(2,2))*Sqr(g2) - 0.025062124532088004*Lambdax*Quad(Yu(2,2))*Sqr(
      g2) - 0.6447330669284618*Quad(Yu(2,2))*Sqr(g1)*Sqr(g2) +
      2.2954745178309515*Power6(Yu(2,2))*Sqr(g3) - 3.0371855509499848*Lambdax*
      Quad(Yu(2,2))*Sqr(g3) + 0.16101073097774002*Quad(Yu(2,2))*Sqr(g1)*Sqr(g3)
      + 0.12233095693773434*Quad(Yu(2,2))*Sqr(g2)*Sqr(g3) + 4.050985151864039*
      Quad(Yu(2,2))*Sqr(Lambdax) - 1.2178606965749863*Power6(Yu(2,2))*Sqr(Yd(2,
      2)) + 4.458058162909423*Lambdax*Quad(Yu(2,2))*Sqr(Yd(2,2)) -
      0.10848390353339003*Quad(Yu(2,2))*Sqr(g1)*Sqr(Yd(2,2)) +
      0.13659793814432988*Quad(Yu(2,2))*Sqr(g2)*Sqr(Yd(2,2)) -
      0.5470674121572557*Quad(Yu(2,2))*Sqr(g3)*Sqr(Yd(2,2)) -
      0.3402061855670103*Power6(Yu(2,2))*Sqr(Ye(2,2)) + 1.099656357388316*
      Lambdax*Quad(Yu(2,2))*Sqr(Ye(2,2)) + 0.05154639175257732*Quad(Yu(2,2))*
      Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 1.*Cube(Lambdax)*Sqr(Yu(2,2)) -
      1.2178606965749863*Power6(Yd(2,2))*Sqr(Yu(2,2)) - 0.3402061855670103*
      Power6(Ye(2,2))*Sqr(Yu(2,2)) + 4.458058162909423*Lambdax*Quad(Yd(2,2))*
      Sqr(Yu(2,2)) + 1.099656357388316*Lambdax*Quad(Ye(2,2))*Sqr(Yu(2,2)) +
      0.040168042971726645*Quad(Yd(2,2))*Sqr(g1)*Sqr(Yu(2,2)) +
      0.13659793814432988*Quad(Yd(2,2))*Sqr(g2)*Sqr(Yu(2,2)) -
      0.5470674121572557*Quad(Yd(2,2))*Sqr(g3)*Sqr(Yu(2,2)) -
      0.0405119402990832*Quad(g1)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.2642469996319474*Quad(g2)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      1.7594501718213058*Quad(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.10996563573883161*Quad(Ye(2,2))*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.2194385717470402*Lambdax*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.31083128169247154*Lambdax*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.15977437983399062*Sqr(g1)*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.15302388409310907*Lambdax*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.9841691885756799*Sqr(g2)*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      2.1112878907901247*Sqr(Lambdax)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) +
      0.0963573883161512*Quad(g1)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.010309278350515464*Quad(g2)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.05154639175257732*Quad(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2)) -
      0.024742268041237116*Lambdax*Sqr(g1)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) -
      0.12371134020618557*Lambdax*Sqr(g2)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.03986254295532646*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) -
      0.4948453608247423*Sqr(Lambdax)*Sqr(Ye(2,2))*Sqr(Yu(2,2)) +
      0.09621993127147767*Lambdax*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2))));

   beta_Lambdax = beta_Lambdax_1 + beta_Lambdax_2;

   return beta_Lambdax;
}

double Standard_model::calc_beta_Lambdax_four_loop(const Beta_traces&) const
{
   double beta_Lambdax;

   beta_Lambdax = Re(16616.342695047893*Power6(g3)*Quad(oneOver16PiSqr)*Quad(Yu
      (2,2)));

   return beta_Lambdax;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 0.85*Sqr(g1) - 2.25*Sqr(g2) - 8*Sqr(g3)) - 1.5*(Yu*
      Yd.adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0016666666666666668*Yu*(1187*Quad(g1) + 5*Sqr(g1
      )*(75*traceYdAdjYd + 225*traceYeAdjYe + 255*traceYuAdjYu - 54*Sqr(g2) +
      152*Sqr(g3)) - 75*(46*Quad(g2) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Quad(g3) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax))))) + 0.0125*(-43*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yu*Yd.adjoint
      ()*Yd) + 0.0125*(223*Sqr(g1) + 675*Sqr(g2) + 20*(-3*(9*traceYdAdjYd + 3*
      traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 64*Sqr(g3)))*(Yu*Yu.adjoint(
      )*Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()
      *Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd + 1.5*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_three_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (16.09896843832796*threeLoop*(6.1494623322761255*PROJECTOR*Lambdax
      *Power5(Yu(2,2)) + 3.6401567229074585*PROJECTOR*Power7(Yu(2,2)) -
      1.5170193131852032*PROJECTOR*Cube(Yu(2,2))*Quad(g1) + 1.0553414347426693*
      PROJECTOR*Cube(Yu(2,2))*Quad(g2) + 22.59550733368146*PROJECTOR*Cube(Yu(2,
      2))*Quad(g3) + 2.821688161099187*PROJECTOR*Cube(Yu(2,2))*Quad(Yd(2,2)) +
      1.6072458368448965*PROJECTOR*Cube(Yu(2,2))*Quad(Ye(2,2)) -
      0.3944352102015494*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(g1) -
      1.8922019827975904*PROJECTOR*Power5(Yu(2,2))*Sqr(g1) - 2.096407613275952*
      PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(g2) - 6.184402459164058*PROJECTOR*
      Power5(Yu(2,2))*Sqr(g2) + 2.163462366965912*PROJECTOR*Cube(Yu(2,2))*Sqr(
      g1)*Sqr(g2) + 0.49692624907281824*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(g3)
      - 9.752177638054057*PROJECTOR*Power5(Yu(2,2))*Sqr(g3) +
      1.1226836416869554*PROJECTOR*Cube(Yu(2,2))*Sqr(g1)*Sqr(g3) +
      3.0045554008024813*PROJECTOR*Cube(Yu(2,2))*Sqr(g2)*Sqr(g3) +
      0.058233544813220885*PROJECTOR*Cube(Yu(2,2))*Sqr(Lambdax) +
      2.8883838227357557*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(Yd(2,2)) +
      2.8689726411313488*PROJECTOR*Power5(Yu(2,2))*Sqr(Yd(2,2)) -
      0.49957993142419754*PROJECTOR*Cube(Yu(2,2))*Sqr(g1)*Sqr(Yd(2,2)) -
      4.814159761919977*PROJECTOR*Cube(Yu(2,2))*Sqr(g2)*Sqr(Yd(2,2)) -
      0.7122084216159787*PROJECTOR*Cube(Yu(2,2))*Sqr(g3)*Sqr(Yd(2,2)) +
      0.9317367170115342*PROJECTOR*Cube(Yu(2,2))*Lambdax*Sqr(Ye(2,2)) +
      0.6522157019080739*PROJECTOR*Power5(Yu(2,2))*Sqr(Ye(2,2)) -
      0.4242586654541776*PROJECTOR*Cube(Yu(2,2))*Sqr(g1)*Sqr(Ye(2,2)) -
      1.9298448995321542*PROJECTOR*Cube(Yu(2,2))*Sqr(g2)*Sqr(Ye(2,2)) +
      0.1552894528352557*PROJECTOR*Cube(Yu(2,2))*Sqr(g3)*Sqr(Ye(2,2)) +
      0.21740523396935796*PROJECTOR*Cube(Yu(2,2))*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      0.2795210151034602*PROJECTOR*Cube(Lambdax)*Yu(2,2) + 1.*PROJECTOR*Power6(
      g1)*Yu(2,2) + 10.549072334227892*PROJECTOR*Power6(g2)*Yu(2,2) -
      38.47142449015523*PROJECTOR*Power6(g3)*Yu(2,2) + 2.187826890843716*
      PROJECTOR*Power6(Yd(2,2))*Yu(2,2) + 0.4996388893047732*PROJECTOR*Power6(
      Ye(2,2))*Yu(2,2) - 0.08455510706879672*PROJECTOR*Lambdax*Quad(g1)*Yu(2,2)
      - 0.33193120543535903*PROJECTOR*Lambdax*Quad(g2)*Yu(2,2) +
      0.4658683585057671*PROJECTOR*Lambdax*Quad(Yd(2,2))*Yu(2,2) +
      0.4658683585057671*PROJECTOR*Lambdax*Quad(Ye(2,2))*Yu(2,2) -
      0.29462264646735065*PROJECTOR*Quad(g2)*Sqr(g1)*Yu(2,2) -
      0.9376834535939306*PROJECTOR*Quad(g3)*Sqr(g1)*Yu(2,2) -
      0.23043972650846906*PROJECTOR*Quad(Yd(2,2))*Sqr(g1)*Yu(2,2) -
      0.5080205796037476*PROJECTOR*Quad(Ye(2,2))*Sqr(g1)*Yu(2,2) -
      0.2759050553623628*PROJECTOR*Quad(g1)*Sqr(g2)*Yu(2,2) + 4.575063689774406
      *PROJECTOR*Quad(g3)*Sqr(g2)*Yu(2,2) - 2.0795715998030686*PROJECTOR*Quad(
      Yd(2,2))*Sqr(g2)*Yu(2,2) - 0.5509041095110554*PROJECTOR*Quad(Ye(2,2))*Sqr
      (g2)*Yu(2,2) + 0.09084432990862458*PROJECTOR*Lambdax*Sqr(g1)*Sqr(g2)*Yu(2
      ,2) - 1.3863544090646365*PROJECTOR*Quad(g1)*Sqr(g3)*Yu(2,2) -
      1.3089127804653757*PROJECTOR*Quad(g2)*Sqr(g3)*Yu(2,2) +
      0.3148250285229063*PROJECTOR*Quad(Yd(2,2))*Sqr(g3)*Yu(2,2) -
      0.9969582872023416*PROJECTOR*Sqr(g1)*Sqr(g2)*Sqr(g3)*Yu(2,2) +
      0.1397605075517301*PROJECTOR*Sqr(g1)*Sqr(Lambdax)*Yu(2,2) +
      0.6988025377586506*PROJECTOR*Sqr(g2)*Sqr(Lambdax)*Yu(2,2) -
      0.20587853488920188*PROJECTOR*Quad(g1)*Sqr(Yd(2,2))*Yu(2,2) +
      0.40913682599408235*PROJECTOR*Quad(g2)*Sqr(Yd(2,2))*Yu(2,2) -
      12.757991577276114*PROJECTOR*Quad(g3)*Sqr(Yd(2,2))*Yu(2,2) +
      0.8230341000268552*PROJECTOR*Quad(Ye(2,2))*Sqr(Yd(2,2))*Yu(2,2) +
      0.5641039159322749*PROJECTOR*Sqr(g1)*Sqr(g2)*Sqr(Yd(2,2))*Yu(2,2) -
      1.3643639389175877*PROJECTOR*Sqr(g1)*Sqr(g3)*Sqr(Yd(2,2))*Yu(2,2) -
      8.902567024109379*PROJECTOR*Sqr(g2)*Sqr(g3)*Sqr(Yd(2,2))*Yu(2,2) -
      1.129730769376485*PROJECTOR*Sqr(Lambdax)*Sqr(Yd(2,2))*Yu(2,2) -
      1.5400233128895373*PROJECTOR*Quad(g1)*Sqr(Ye(2,2))*Yu(2,2) -
      0.9718846178821886*PROJECTOR*Quad(g2)*Sqr(Ye(2,2))*Yu(2,2) +
      1.3665471849502502*PROJECTOR*Quad(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2) -
      0.3364704668151931*PROJECTOR*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2))*Yu(2,2) -
      0.3494012688793253*PROJECTOR*Sqr(Lambdax)*Sqr(Ye(2,2))*Yu(2,2) -
      0.14904312779958145*PROJECTOR*Sqr(g1)*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2) -
      0.5159639826231229*PROJECTOR*Sqr(g2)*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2) -
      0.445163098127733*PROJECTOR*Sqr(g3)*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Yu(2,2))).
      real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_four_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (2308.1827208150457*PROJECTOR*Power8(g3)*Quad(oneOver16PiSqr)*Yu(2
      ,2)).real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yd_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(Sqr(g1) + 9*Sqr(g2) + 4*(-3*
      traceYdAdjYd - traceYeAdjYe - 3*traceYuAdjYu + 8*Sqr(g3))) + 1.5*(Yd*
      Yd.adjoint()*Yd) - 1.5*(Yd*Yu.adjoint()*Yu))).real();

   return beta_Yd;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yd_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0016666666666666668*Yd*(127*Quad(g1) + 5*Sqr(g1
      )*(-15*(5*traceYdAdjYd + 15*traceYeAdjYe + 17*traceYuAdjYu) + 162*Sqr(g2)
      - 248*Sqr(g3)) + 75*(46*Quad(g2) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Quad(g3) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax))))) + 0.0125*(187*Sqr(g1) + 675*Sqr(g2) + 20*(-3*(9*
      traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 64*Sqr(g3))
      )*(Yd*Yd.adjoint()*Yd) + 0.0125*(-79*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yd*Yu.adjoint
      ()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - Yd*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu - 0.25*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 2.75*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();

   return beta_Yd;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yd_three_loop(const Beta_traces&) const
{
   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (ZEROMATRIX(3,3)).real();

   return beta_Yd;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Ye_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 2.25*Sqr(g1) - 2.25*Sqr(g2)) + 1.5*(Ye*Ye.adjoint()*Ye)))
      .real();

   return beta_Ye;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Ye_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.005*Ye*(1371*Quad(g1) + 5*Sqr(g1)*(25*
      traceYdAdjYd + 75*traceYeAdjYe + 85*traceYuAdjYu + 54*Sqr(g2)) - 25*(46*
      Quad(g2) - 15*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*Sqr(g2) -
      2*(80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) - 3*(9*traceYdAdjYdYdAdjYd -
      2*traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2
      *Sqr(Lambdax))))) + 0.0375*(129*Sqr(g1) + 5*(-4*(9*traceYdAdjYd + 3*
      traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 45*Sqr(g2)))*(Ye*Ye.adjoint(
      )*Ye) + 1.5*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();

   return beta_Ye;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Ye_three_loop(const Beta_traces&) const
{
   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (ZEROMATRIX(3,3)).real();

   return beta_Ye;
}

double Standard_model::calc_beta_mu2_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_mu2;

   beta_mu2 = Re(oneOver16PiSqr*(2*mu2*(3*traceYdAdjYd + traceYeAdjYe + 3
      *(traceYuAdjYu + Lambdax)) - 0.9*mu2*Sqr(g1) - 4.5*mu2*Sqr(g2)));

   return beta_mu2;
}

double Standard_model::calc_beta_mu2_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   double beta_mu2;

   beta_mu2 = Re(0.0025*mu2*twoLoop*(1671*Quad(g1) + 10*Sqr(g1)*(50*
      traceYdAdjYd + 150*traceYeAdjYe + 170*traceYuAdjYu + 288*Lambdax + 45*Sqr
      (g2)) - 25*(145*Quad(g2) - 12*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*
      traceYuAdjYu + 48*Lambdax)*Sqr(g2) - 8*(80*(traceYdAdjYd + traceYuAdjYu)*
      Sqr(g3) - 3*(9*traceYdAdjYdYdAdjYd + 14*traceYdAdjYuYuAdjYd + 3*
      traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu + 24*traceYdAdjYd*Lambdax + 8
      *traceYeAdjYe*Lambdax + 24*traceYuAdjYu*Lambdax + 10*Sqr(Lambdax))))));

   return beta_mu2;
}

double Standard_model::calc_beta_mu2_three_loop(const Beta_traces&) const
{
   double beta_mu2;

   beta_mu2 = Re(2*mu2*threeLoop*(128.25*Cube(Lambdax) +
      8.378314604562993*Power6(g1) + 301.7235994495886*Power6(g2) +
      154.40506064218175*Power6(Yd(2,2)) + 3.4683535473939138*Power6(Ye(2,2)) +
      154.40506064218175*Power6(Yu(2,2)) + 173.69714554123618*Lambdax*Quad(Yd(
      2,2)) + 75.89904851374538*Lambdax*Quad(Ye(2,2)) + 173.69714554123618*
      Lambdax*Quad(Yu(2,2)) - 209.24048513745396*Quad(Yd(2,2))*Sqr(g3) -
      209.24048513745396*Quad(Yu(2,2))*Sqr(g3) + 178.48396765750311*Quad(g3)*
      Sqr(Yd(2,2)) + 72.*Quad(Ye(2,2))*Sqr(Yd(2,2)) + 296.2115485137454*Quad(Yu
      (2,2))*Sqr(Yd(2,2)) + 40.19238810996319*Lambdax*Sqr(g3)*Sqr(Yd(2,2)) +
      37.125*Sqr(Lambdax)*Sqr(Yd(2,2)) + 72.*Quad(Yd(2,2))*Sqr(Ye(2,2)) + 72.*
      Quad(Yu(2,2))*Sqr(Ye(2,2)) + 12.375*Sqr(Lambdax)*Sqr(Ye(2,2)) - 54.*
      Lambdax*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + Quad(g2)*(-32.25723415592714*Lambdax
      - 28.572145541236182*Sqr(g3) - 102.62651936694535*Sqr(Yd(2,2)) -
      34.208839788981784*Sqr(Ye(2,2)) - 102.62651936694535*Sqr(Yu(2,2))) + Quad
      (g1)*(-32.90278936623708*Lambdax + 9.778005607297105*Sqr(g2) -
      4.190581346047974*Sqr(g3) - 11.837888914984243*Sqr(Yd(2,2)) -
      20.542464319265452*Sqr(Ye(2,2)) - 27.721423523790055*Sqr(Yu(2,2))) +
      178.48396765750311*Quad(g3)*Sqr(Yu(2,2)) + 296.2115485137454*Quad(Yd(2,2)
      )*Sqr(Yu(2,2)) + 72.*Quad(Ye(2,2))*Sqr(Yu(2,2)) + 40.19238810996319*
      Lambdax*Sqr(g3)*Sqr(Yu(2,2)) + 37.125*Sqr(Lambdax)*Sqr(Yu(2,2)) -
      208.57214554123618*Lambdax*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 16.698731351660527
      *Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 54.*Lambdax*Sqr(Ye(2,2))*Sqr(Yu(2,2)
      ) + 10.5*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g2)*(
      -3.829281688145728*Quad(Yd(2,2)) + 3.223572770618091*Quad(Ye(2,2)) -
      3.829281688145728*Quad(Yu(2,2)) - 48.205536385309045*Sqr(Lambdax) -
      53.09857277061809*Lambdax*Sqr(Ye(2,2)) + Sqr(Yd(2,2))*(
      -159.29571831185424*Lambdax - 13.500000000000004*Sqr(Ye(2,2)) -
      62.83053638530905*Sqr(Yu(2,2))) - 159.29571831185424*Lambdax*Sqr(Yu(2,2))
      - 13.500000000000004*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g3)*(
      7.572145541236182*Sqr(Yd(2,2)) + 7.572145541236182*Sqr(Yu(2,2)))) + Sqr(
      g1)*(9.930778506201856*Quad(g2) + 11.916358216494471*Quad(Yd(2,2)) -
      15.051929108247235*Quad(Ye(2,2)) - 7.507690297250921*Quad(Yu(2,2)) -
      9.641107277061808*Sqr(Lambdax) - 32.86395336511993*Lambdax*Sqr(Yd(2,2)) -
      7.605285445876382*Lambdax*Sqr(Ye(2,2)) - 2.7*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      29.849524256872694*Lambdax*Sqr(Yu(2,2)) - 23.946234141895758*Sqr(Yd(2,2)
      )*Sqr(Yu(2,2)) - 2.7*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g3)*(
      -2.091983828751535*Sqr(Yd(2,2)) + 8.727254982244787*Sqr(Yu(2,2))) + Sqr(
      g2)*(-18.911535445876382*Lambdax - 8.109375*Sqr(Yd(2,2)) +
      10.602411385309043*Sqr(Ye(2,2)) + 11.470322300901756*Sqr(Yu(2,2))))));

   return beta_mu2;
}

double Standard_model::calc_beta_v_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_v;

   beta_v = Re(oneOver16PiSqr*(0.6*v*Sqr(g1) + v*(-3*traceYdAdjYd -
      traceYeAdjYe - 3*traceYuAdjYu + 3*Sqr(g2))));

   return beta_v;
}

double Standard_model::calc_beta_v_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   double beta_v;

   beta_v = Re(-0.00125*twoLoop*v*(1221*Quad(g1) + 10*Sqr(g1)*(122*
      traceYdAdjYd + 174*traceYeAdjYe + 242*traceYuAdjYu - 45*Sqr(g2)) - 25*(
      379*Quad(g2) - 108*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*Sqr(
      g2) - 8*(80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) - 3*(9*
      traceYdAdjYdYdAdjYd - 2*traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*
      traceYuAdjYuYuAdjYu - 2*Sqr(Lambdax))))));

   return beta_v;
}

double Standard_model::calc_beta_v_three_loop(const Beta_traces&) const
{
   double beta_v;

   beta_v = 0;

   return beta_v;
}

double Standard_model::get_mass_matrix_VG() const
{
   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void Standard_model::calculate_MVG()
{
   const auto mass_matrix_VG = get_mass_matrix_VG();
   MVG = mass_matrix_VG;
}

double Standard_model::get_mass_matrix_Hp() const
{
   const double mass_matrix_Hp = Re(-mu2 + 0.5*Lambdax*Sqr(v) + 0.25*Sqr(
      g2)*Sqr(v));

   return mass_matrix_Hp;
}

void Standard_model::calculate_MHp()
{
   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   MHp = mass_matrix_Hp;

   if (MHp < 0.) {
      problems.flag_running_tachyon(standard_model_info::Hp);
   }

   MHp = AbsSqrt(MHp);
}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,2) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void Standard_model::calculate_MFv()
{
   MFv.setConstant(0);
}

double Standard_model::get_mass_matrix_Ah() const
{
   const double mass_matrix_Ah = Re(0.25*(-4*mu2 + 2*Lambdax*Sqr(v) + Sqr
      (v)*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))));

   return mass_matrix_Ah;
}

void Standard_model::calculate_MAh()
{
   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   MAh = mass_matrix_Ah;

   if (MAh < 0.) {
      problems.flag_running_tachyon(standard_model_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

double Standard_model::get_mass_matrix_hh() const
{
   const double mass_matrix_hh = Re(-mu2 + 1.5*Lambdax*Sqr(v));

   return mass_matrix_hh;
}

void Standard_model::calculate_Mhh()
{
   const auto mass_matrix_hh = get_mass_matrix_hh();
   Mhh = mass_matrix_hh;

   if (Mhh < 0.) {
      problems.flag_running_tachyon(standard_model_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

double Standard_model::get_mass_matrix_VP() const
{
   const double mass_matrix_VP = Re(0);

   return mass_matrix_VP;
}

void Standard_model::calculate_MVP()
{
   const auto mass_matrix_VP = get_mass_matrix_VP();
   MVP = mass_matrix_VP;
}

double Standard_model::get_mass_matrix_VZ() const
{
   const double mass_matrix_VZ = Re(0.25*Sqr(v)*Sqr(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW())));

   return mass_matrix_VZ;
}

void Standard_model::calculate_MVZ()
{
   const auto mass_matrix_VZ = get_mass_matrix_VZ();
   MVZ = mass_matrix_VZ;

   if (MVZ < 0.) {
      problems.flag_running_tachyon(standard_model_info::VZ);
   }

   MVZ = AbsSqrt(MVZ);
}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void Standard_model::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(standard_model_info::Fd, eigenvalue_error > precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = -0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = -0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = -0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = -0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = -0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = -0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = -0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = -0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void Standard_model::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(standard_model_info::Fu, eigenvalue_error > precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void Standard_model::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(standard_model_info::Fe, eigenvalue_error > precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double Standard_model::get_mass_matrix_VWp() const
{
   const double mass_matrix_VWp = Re(0.25*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void Standard_model::calculate_MVWp()
{
   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   MVWp = mass_matrix_VWp;

   if (MVWp < 0.) {
      problems.flag_running_tachyon(standard_model_info::VWp);
   }

   MVWp = AbsSqrt(MVWp);
}

Eigen::Matrix<double,2,2> Standard_model::get_mass_matrix_VPVZ() const
{
   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void Standard_model::calculate_MVPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error
      );
   ZZ.transposeInPlace();
#else
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
   ZZ.transposeInPlace();
#endif


   MVPVZ = AbsSqrt(MVPVZ);
}


double Standard_model::get_ewsb_eq_hh_1() const
{
   double result = Re(-(mu2*v) + 0.5*Cube(v)*Lambdax);

   return result;
}



double Standard_model::CpconjHpHphh() const
{
   const double result = -(v*Lambdax);

   return result;
}

double Standard_model::CpconjHpVWpVP() const
{
   const double result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double Standard_model::CpconjHpVZVWp() const
{
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double Standard_model::CpHpgWpCbargZ() const
{
   const double result = 0.05*g2*v*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpconjHpbargWpCgZ() const
{
   const double result = -0.05*g2*v*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpHpgZbargWp() const
{
   const double result = -0.05*g2*v*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpconjHpbargZgWp() const
{
   const double result = 0.05*g2*v*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpHpconjHpAhAh() const
{
   const double result = -Lambdax;

   return result;
}

double Standard_model::CpHpconjHphhhh() const
{
   const double result = -Lambdax;

   return result;
}

double Standard_model::CpHpconjHpconjHpHp() const
{
   const double result = -2*Lambdax;

   return result;
}

std::complex<double> Standard_model::CpconjHpVWpAh() const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double Standard_model::CpconjHpVWphh() const
{
   const double result = -0.5*g2;

   return result;
}

double Standard_model::CpconjHpVPHp() const
{
   const double result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(
      ThetaW()));

   return result;
}

double Standard_model::CpconjHpVZHp() const
{
   const double result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(
      ThetaW()));

   return result;
}

double Standard_model::CpHpconjHpconjVWpVWp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpHpconjHpVZVZ() const
{
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW
      ())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW(
      ))));

   return result;
}

std::complex<double> Standard_model::CpconjHpbarFdFuPR(int gI1, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*
      Uu(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpconjHpbarFdFuPL(int gI1, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,
      Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

double Standard_model::CpconjHpbarFeFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpconjHpbarFeFvPL(int gI1, int gI2) const
{
   const std::complex<double> result = -SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,gI2))
      ;

   return result;
}

double Standard_model::CpAhhhAh() const
{
   const double result = -(v*Lambdax);

   return result;
}

std::complex<double> Standard_model::CpAhbargWpgWp() const
{
   const std::complex<double> result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpAhbargWpCgWpC() const
{
   const std::complex<double> result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return result;
}

double Standard_model::CpAhAhAhAh() const
{
   const double result = -3*Lambdax;

   return result;
}

double Standard_model::CpAhAhhhhh() const
{
   const double result = -Lambdax;

   return result;
}

double Standard_model::CpAhAhconjHpHp() const
{
   const double result = -Lambdax;

   return result;
}

std::complex<double> Standard_model::CpAhVZhh() const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpAhconjVWpHp() const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double Standard_model::CpAhAhconjVWpVWp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpAhAhVZVZ() const
{
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834
      *g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> Standard_model::CpAhbarFdFdPR(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(
      gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpAhbarFdFdPL(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*
      Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpAhbarFeFePR(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(
      gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpAhbarFeFePL(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*
      Ye(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpAhbarFuFuPR(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(
      gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpAhbarFuFuPL(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*
      Yu(j1,j2)));

   return result;
}

double Standard_model::CphhAhAh() const
{
   const double result = -(v*Lambdax);

   return result;
}

double Standard_model::Cphhhhhh() const
{
   const double result = -3*v*Lambdax;

   return result;
}

double Standard_model::CphhVZVZ() const
{
   const double result = 0.1*v*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double Standard_model::CphhbargWpgWp() const
{
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double Standard_model::CphhbargWpCgWpC() const
{
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double Standard_model::CphhbargZgZ() const
{
   const double result = -0.05*v*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double Standard_model::CphhconjHpHp() const
{
   const double result = -(v*Lambdax);

   return result;
}

double Standard_model::CphhconjVWpVWp() const
{
   const double result = 0.5*v*Sqr(g2);

   return result;
}

double Standard_model::CphhhhAhAh() const
{
   const double result = -Lambdax;

   return result;
}

double Standard_model::Cphhhhhhhh() const
{
   const double result = -3*Lambdax;

   return result;
}

double Standard_model::CphhhhconjHpHp() const
{
   const double result = -Lambdax;

   return result;
}

std::complex<double> Standard_model::CphhVZAh() const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CphhconjVWpHp() const
{
   const double result = 0.5*g2;

   return result;
}

double Standard_model::CphhhhconjVWpVWp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CphhhhVZVZ() const
{
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834
      *g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> Standard_model::CphhbarFdFdPR(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CphhbarFdFdPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(
      gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CphhbarFeFePR(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CphhbarFeFePL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(
      gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CphhbarFuFuPR(int gI1, int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CphhbarFuFuPL(int gI1, int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(
      gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpVZhhAh() const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpVZVZhh() const
{
   const double result = 0.1*v*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double Standard_model::CpVZbargWpgWp() const
{
   const double result = g2*Cos(ThetaW());

   return result;
}

double Standard_model::CpVZbargWpCgWpC() const
{
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double Standard_model::CpVZconjHpHp() const
{
   const double result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(
      ThetaW()));

   return result;
}

double Standard_model::CpVZconjVWpHp() const
{
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

std::complex<double> Standard_model::CpVZVZAhAh() const
{
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834
      *g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> Standard_model::CpVZVZhhhh() const
{
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834
      *g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> Standard_model::CpVZVZconjHpHp() const
{
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW
      ())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW(
      ))));

   return result;
}

double Standard_model::CpVZconjVWpVWp() const
{
   const double result = g2*Cos(ThetaW());

   return result;
}

double Standard_model::CpVZbarFdFdPL(int gI1, int gI2) const
{
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(15*g2*Cos
      (ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpVZbarFdFdPR(int gI1, int gI2) const
{
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpVZbarFeFePL(int gI1, int gI2) const
{
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpVZbarFeFePR(int gI1, int gI2) const
{
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpVZbarFuFuPL(int gI1, int gI2) const
{
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*
      Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpVZbarFuFuPR(int gI1, int gI2) const
{
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpVZbarFvFvPL(int gI1, int gI2) const
{
   const double result = -0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpVZbarFvFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpVZVZconjVWpVWp1() const
{
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpVZVZconjVWpVWp2() const
{
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpVZVZconjVWpVWp3() const
{
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpconjVWpHpAh() const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double Standard_model::CpconjVWpHphh() const
{
   const double result = 0.5*g2;

   return result;
}

double Standard_model::CpconjVWpVPHp() const
{
   const double result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double Standard_model::CpconjVWpVWphh() const
{
   const double result = 0.5*v*Sqr(g2);

   return result;
}

double Standard_model::CpconjVWpVZHp() const
{
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double Standard_model::CpconjVWpbargPgWp() const
{
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpconjVWpbargWpCgP() const
{
   const double result = g2*Sin(ThetaW());

   return result;
}

double Standard_model::CpconjVWpbargWpCgZ() const
{
   const double result = g2*Cos(ThetaW());

   return result;
}

double Standard_model::CpconjVWpbargZgWp() const
{
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpAhAh() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

double Standard_model::CpVWpconjVWphhhh() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

double Standard_model::CpVWpconjVWpconjHpHp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

double Standard_model::CpconjVWpVWpVP() const
{
   const double result = g2*Sin(ThetaW());

   return result;
}

double Standard_model::CpconjVWpVZVWp() const
{
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpconjVWpbarFdFuPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      Vu(gI2,j1))*Vd(gI1,j1));

   return result;
}

double Standard_model::CpconjVWpbarFdFuPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpconjVWpbarFeFvPL(int gI1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*Ve(gI1
      ,gI2),0);

   return result;
}

double Standard_model::CpconjVWpbarFeFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpVWpconjVWpVPVP1() const
{
   const double result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpVPVP2() const
{
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpVPVP3() const
{
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpVZVZ1() const
{
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpVZVZ2() const
{
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpVZVZ3() const
{
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpVWpconjVWpconjVWpVWp1() const
{
   const double result = -Sqr(g2);

   return result;
}

double Standard_model::CpVWpconjVWpconjVWpVWp2() const
{
   const double result = -Sqr(g2);

   return result;
}

double Standard_model::CpVWpconjVWpconjVWpVWp3() const
{
   const double result = 2*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI1,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(gI1,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdhhFdPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(Vd(gI2,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdhhFdPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Yd(j1,gO1))*Ud(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdVGFdPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Ud(gI2,gO2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdVGFdPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vd(gI2,gO1))),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdVPFdPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(
      ThetaW())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdVPFdPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(
      Vd(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdVZFdPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Sin(
      ThetaW())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdVZFdPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdconjHpFuPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-SUM(j2,0,2,Conj(Vu(gI2,j2))*
      Yd(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdconjHpFuPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yu(j1,gO1))*
      Uu(gI2,j1)),0);

   return result;
}

double Standard_model::CpbarUFdconjVWpFuPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarUFdconjVWpFuPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      Vu(gI2,gO1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI1,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(gI1,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuhhFuPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*SUM(j2,0,2
      ,Conj(Vu(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuhhFuPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.7071067811865475*SUM(j1,0,2
      ,Conj(Yu(j1,gO1))*Uu(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuHpFdPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-SUM(j2,0,2,Conj(Vd(gI2,j2))*
      Yu(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuHpFdPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yd(j1,gO1))*
      Ud(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuVGFuPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Uu(gI2,gO2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuVGFuPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vu(gI2,gO1))),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuVPFuPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.5163977794943222*g1*Cos(
      ThetaW())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuVPFuPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(
      Vu(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

double Standard_model::CpbarUFuVWpFdPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarUFuVWpFdPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      Vd(gI2,gO1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuVZFuPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Sin(
      ThetaW())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuVZFuPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Cos
      (ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI1,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI1,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFehhFePL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(Ve(gI2,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFehhFePR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Ye(j1,gO1))*Ue(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeVPFePR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.7745966692414834*g1*Cos(
      ThetaW())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeVPFePL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.3872983346207417*g1*Conj(Ve
      (gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeVZFePR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Sin(
      ThetaW())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeVZFePL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeconjHpFvPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-Ye(gO2,gI2),0);

   return result;
}

double Standard_model::CpbarUFeconjHpFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpbarUFeconjVWpFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpbarUFeconjVWpFvPL(int gO1, int gI2) const
{
   const double result = IF(gI2 < 3,-0.7071067811865475*g2*KroneckerDelta(gI2,
      gO1),0);

   return result;
}

std::complex<double> Standard_model::CpbarFdFdAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI1,j2))*SUM(j1,0,2,Conj(Ud(gO2,j1))*
      Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFdFdAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI1,j1))*Vd(
      gO1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFdhhFdPL(int gO2, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(
      gI2,j2))*SUM(j1,0,2,Conj(Ud(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFdhhFdPR(int gO1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gO1,j2));

   return result;
}

double Standard_model::CpbarFdVZFdPR(int gO2, int gI2) const
{
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI2,gO2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpbarFdVZFdPL(int gO1, int gI2) const
{
   const double result = 0.03333333333333333*KroneckerDelta(gI2,gO1)*(15*g2*Cos
      (ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpbarFdconjHpFuPL(int gO2, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,
      Conj(Ud(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFdconjHpFuPR(int gO1, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*
      Uu(gI2,j1))*Vd(gO1,j2));

   return result;
}

double Standard_model::CpbarFdconjVWpFuPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFdconjVWpFuPL(int gO1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      Vu(gI2,j1))*Vd(gO1,j1));

   return result;
}

std::complex<double> Standard_model::CpbarFeFeAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI1,j2))*SUM(j1,0,2,Conj(Ue(gO2,j1))*
      Ye(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFeFeAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI1,j1))*Ve(
      gO1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFehhFePL(int gO2, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(
      gI2,j2))*SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFehhFePR(int gO1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gO1,j2));

   return result;
}

double Standard_model::CpbarFeVZFePR(int gO2, int gI2) const
{
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI2,gO2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpbarFeVZFePL(int gO1, int gI2) const
{
   const double result = 0.1*KroneckerDelta(gI2,gO1)*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpbarFeconjHpFvPL(int gO2, int gI2) const
{
   const std::complex<double> result = -SUM(j1,0,2,Conj(Ue(gO2,j1))*Ye(j1,gI2))
      ;

   return result;
}

double Standard_model::CpbarFeconjHpFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpbarFeconjVWpFvPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFeconjVWpFvPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*Ve(gO1
      ,gI2),0);

   return result;
}

std::complex<double> Standard_model::CpbarFuFuAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI1,j2))*SUM(j1,0,2,Conj(Uu(gO2,j1))*
      Yu(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFuFuAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI1,j1))*Vu(
      gO1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFuhhFuPL(int gO2, int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,Conj(Vu(
      gI2,j2))*SUM(j1,0,2,Conj(Uu(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFuhhFuPR(int gO1, int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gO1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFuHpFdPL(int gO2, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,
      Conj(Uu(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFuHpFdPR(int gO1, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*
      Ud(gI2,j1))*Vu(gO1,j2));

   return result;
}

double Standard_model::CpbarFuVPFuPR(int gO2, int gI2) const
{
   const double result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(
      gI2,gO2);

   return result;
}

double Standard_model::CpbarFuVPFuPL(int gO1, int gI2) const
{
   const double result = -0.03333333333333333*KroneckerDelta(gI2,gO1)*(
      3.872983346207417*g1*Cos(ThetaW()) + 15*g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFuVWpFdPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFuVWpFdPL(int gO1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      Vd(gI2,j1))*Vu(gO1,j1));

   return result;
}

double Standard_model::CpbarFuVZFuPR(int gO2, int gI2) const
{
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI2,gO2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpbarFuVZFuPL(int gO1, int gI2) const
{
   const double result = 0.03333333333333333*KroneckerDelta(gI2,gO1)*(-15*g2*
      Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}


std::complex<double> Standard_model::self_energy_Hp_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjHpHphh())*B0(p,MHp,Mhh);
   result += 2*AbsSqr(CpconjHpVWpVP())*(-1 + 2*B0(p,0,MVWp));
   result += 2*AbsSqr(CpconjHpVZVWp())*(-1 + 2*B0(p,MVWp,MVZ));
   result += -0.5*A0(MAh)*CpHpconjHpAhAh();
   result += -(A0(MHp)*CpHpconjHpconjHpHp());
   result += -0.5*A0(Mhh)*CpHpconjHphhhh();
   result += -(B0(p,MVZ,MVWp)*CpconjHpbargWpCgZ()*CpHpgWpCbargZ());
   result += -(B0(p,MVWp,MVZ)*CpconjHpbargZgWp()*CpHpgZbargWp());
   result += AbsSqr(CpconjHpVWpAh())*F0(p,MAh,MVWp);
   result += AbsSqr(CpconjHpVWphh())*F0(p,Mhh,MVWp);
   result += AbsSqr(CpconjHpVPHp())*F0(p,MHp,0);
   result += AbsSqr(CpconjHpVZHp())*F0(p,MHp,MVZ);
   result += 4*A0(MVWp)*CpHpconjHpconjVWpVWp() - 2*CpHpconjHpconjVWpVWp()*Sqr(
      MVWp);
   result += CpHpconjHpVZVZ()*(2*A0(MVZ) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpconjHpbarFdFuPL(gI1,gI2)) +
      AbsSqr(CpconjHpbarFdFuPR(gI1,gI2)))*G0(p,MFd(gI1),MFu(gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpconjHpbarFeFvPL(gI1,gI2)) +
      AbsSqr(CpconjHpbarFeFvPR(gI1,gI2)))*G0(p,MFe(gI1),MFv(gI2))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(p,MFd(gI1),MFu(gI2))*(Conj(
      CpconjHpbarFdFuPR(gI1,gI2))*CpconjHpbarFdFuPL(gI1,gI2) + Conj(
      CpconjHpbarFdFuPL(gI1,gI2))*CpconjHpbarFdFuPR(gI1,gI2))*MFu(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(p,MFe(gI1),MFv(gI2))*(Conj(
      CpconjHpbarFeFvPR(gI1,gI2))*CpconjHpbarFeFvPL(gI1,gI2) + Conj(
      CpconjHpbarFeFvPL(gI1,gI2))*CpconjHpbarFeFvPR(gI1,gI2))*MFv(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Ah_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CpAhAhAhAh();
   result += -(A0(MHp)*CpAhAhconjHpHp());
   result += -0.5*A0(Mhh)*CpAhAhhhhh();
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpgWp()));
   result += AbsSqr(CpAhhhAh())*B0(p,Mhh,MAh);
   result += AbsSqr(CpAhVZhh())*F0(p,Mhh,MVZ);
   result += 2*AbsSqr(CpAhconjVWpHp())*F0(p,MHp,MVWp);
   result += 4*A0(MVWp)*CpAhAhconjVWpVWp() - 2*CpAhAhconjVWpVWp()*Sqr(MVWp);
   result += CpAhAhVZVZ()*(2*A0(MVZ) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpAhbarFdFdPL(gI1,gI2)) + AbsSqr
      (CpAhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpAhbarFeFePL(gI1,gI2)) + AbsSqr(
      CpAhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpAhbarFuFuPL(gI1,gI2)) + AbsSqr
      (CpAhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(p,MFd(gI1),MFd(gI2))*(Conj(
      CpAhbarFdFdPR(gI1,gI2))*CpAhbarFdFdPL(gI1,gI2) + Conj(CpAhbarFdFdPL(gI1,gI2)
      )*CpAhbarFdFdPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(p,MFe(gI1),MFe(gI2))*(Conj(
      CpAhbarFeFePR(gI1,gI2))*CpAhbarFeFePL(gI1,gI2) + Conj(CpAhbarFeFePL(gI1,gI2)
      )*CpAhbarFeFePR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(p,MFu(gI1),MFu(gI2))*(Conj(
      CpAhbarFuFuPR(gI1,gI2))*CpAhbarFuFuPL(gI1,gI2) + Conj(CpAhbarFuFuPL(gI1,gI2)
      )*CpAhbarFuFuPR(gI1,gI2))*MFu(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_hh_1loop(double p ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CphhAhAh())*B0(p,MAh,MAh);
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpgWp()));
   result += -(B0(p,MVZ,MVZ)*Sqr(CphhbargZgZ()));
   result += AbsSqr(CphhconjHpHp())*B0(p,MHp,MHp);
   result += 2*AbsSqr(CphhconjVWpVWp())*(-1 + 2*B0(p,MVWp,MVWp));
   result += -0.5*A0(MAh)*CphhhhAhAh();
   result += -(A0(MHp)*CphhhhconjHpHp());
   result += 0.5*AbsSqr(Cphhhhhh())*B0(p,Mhh,Mhh);
   result += -0.5*A0(Mhh)*Cphhhhhhhh();
   result += AbsSqr(CphhVZVZ())*(-1 + 2*B0(p,MVZ,MVZ));
   result += AbsSqr(CphhVZAh())*F0(p,MAh,MVZ);
   result += 2*AbsSqr(CphhconjVWpHp())*F0(p,MHp,MVWp);
   result += 4*A0(MVWp)*CphhhhconjVWpVWp() - 2*CphhhhconjVWpVWp()*Sqr(MVWp);
   result += CphhhhVZVZ()*(2*A0(MVZ) - Sqr(MVZ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CphhbarFdFdPL(gI1,gI2)) + AbsSqr
      (CphhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CphhbarFeFePL(gI1,gI2)) + AbsSqr(
      CphhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CphhbarFuFuPL(gI1,gI2)) + AbsSqr
      (CphhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(p,MFd(gI1),MFd(gI2))*(Conj(
      CphhbarFdFdPR(gI1,gI2))*CphhbarFdFdPL(gI1,gI2) + Conj(CphhbarFdFdPL(gI1,gI2)
      )*CphhbarFdFdPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(p,MFe(gI1),MFe(gI2))*(Conj(
      CphhbarFeFePR(gI1,gI2))*CphhbarFeFePL(gI1,gI2) + Conj(CphhbarFeFePL(gI1,gI2)
      )*CphhbarFeFePR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(p,MFu(gI1),MFu(gI2))*(Conj(
      CphhbarFuFuPR(gI1,gI2))*CphhbarFuFuPL(gI1,gI2) + Conj(CphhbarFuFuPL(gI1,gI2)
      )*CphhbarFuFuPR(gI1,gI2))*MFu(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_VZ_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWpCgWpC())*B00(p,MVWp,MVWp);
   result += AbsSqr(CpVZbargWpgWp())*B00(p,MVWp,MVWp);
   result += -4*AbsSqr(CpVZconjHpHp())*B00(p,MHp,MHp);
   result += 2*AbsSqr(CpVZconjVWpHp())*B0(p,MVWp,MHp);
   result += -4*AbsSqr(CpVZhhAh())*B00(p,MAh,Mhh);
   result += 0.5*A0(MAh)*CpVZVZAhAh();
   result += A0(MHp)*CpVZVZconjHpHp();
   result += -(A0(MVWp)*(4*CpVZVZconjVWpVWp1() + CpVZVZconjVWpVWp2() +
      CpVZVZconjVWpVWp3()));
   result += AbsSqr(CpVZVZhh())*B0(p,MVZ,Mhh);
   result += 0.5*A0(Mhh)*CpVZVZhhhh();
   result += 2*CpVZVZconjVWpVWp1()*Sqr(MVWp);
   result += -0.6666666666666666*AbsSqr(CpVZconjVWpVWp())*(3*A0(MVWp) + 15*B00(
      p,MVWp,MVWp) - 6*Sqr(MVWp) + Sqr(p) + 3*B0(p,MVWp,MVWp)*(Sqr(MVWp) + 2*Sqr(p
      )));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr
      (CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2)) + 4*B0(p,MFd(gI1),MFd(gI2)
      )*MFd(gI1)*MFd(gI2)*Re(Conj(CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2)))
      );
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
      CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2)) + 4*B0(p,MFe(gI1),MFe(gI2))
      *MFe(gI1)*MFe(gI2)*Re(Conj(CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2))))
      ;
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr
      (CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2)) + 4*B0(p,MFu(gI1),MFu(gI2)
      )*MFu(gI1)*MFu(gI2)*Re(Conj(CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2)))
      );
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
      CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2)) + 4*B0(p,MFv(gI1),MFv(gI2))
      *MFv(gI1)*MFv(gI2)*Re(Conj(CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2))))
      ;

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_VWp_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWpbargPgWp())*B00(p,MVWp,MVP);
   result += AbsSqr(CpconjVWpbargWpCgP())*B00(p,MVP,MVWp);
   result += AbsSqr(CpconjVWpbargWpCgZ())*B00(p,MVZ,MVWp);
   result += AbsSqr(CpconjVWpbargZgWp())*B00(p,MVWp,MVZ);
   result += -4*AbsSqr(CpconjVWpHpAh())*B00(p,MAh,MHp);
   result += -4*AbsSqr(CpconjVWpHphh())*B00(p,Mhh,MHp);
   result += AbsSqr(CpconjVWpVPHp())*B0(p,0,MHp);
   result += AbsSqr(CpconjVWpVWphh())*B0(p,MVWp,Mhh);
   result += AbsSqr(CpconjVWpVZHp())*B0(p,MVZ,MHp);
   result += 0.5*A0(MAh)*CpVWpconjVWpAhAh();
   result += A0(MHp)*CpVWpconjVWpconjHpHp();
   result += -(A0(MVWp)*(4*CpVWpconjVWpconjVWpVWp1() + CpVWpconjVWpconjVWpVWp2(
      ) + CpVWpconjVWpconjVWpVWp3()));
   result += 0.5*A0(Mhh)*CpVWpconjVWphhhh();
   result += 0;
   result += 2*CpVWpconjVWpconjVWpVWp1()*Sqr(MVWp);
   result += -0.3333333333333333*AbsSqr(CpconjVWpVWpVP())*(3*A0(MVWp) + 30*B00(
      p,MVWp,0) - 6*Sqr(MVWp) + 2*Sqr(p) + 3*B0(p,MVWp,0)*(Sqr(MVWp) + 4*Sqr(p)));
   result += -0.5*A0(MVZ)*(4*CpVWpconjVWpVZVZ1() + CpVWpconjVWpVZVZ2() +
      CpVWpconjVWpVZVZ3()) + CpVWpconjVWpVZVZ1()*Sqr(MVZ);
   result += -0.3333333333333333*AbsSqr(CpconjVWpVZVWp())*(3*A0(MVWp) + 3*A0(
      MVZ) + 30*B00(p,MVZ,MVWp) - 6*Sqr(MVWp) + 3*B0(p,MVZ,MVWp)*Sqr(MVWp) - 6*Sqr
      (MVZ) + 3*B0(p,MVZ,MVWp)*Sqr(MVZ) + 2*Sqr(p) + 12*B0(p,MVZ,MVWp)*Sqr(p));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpconjVWpbarFdFuPL(gI1,gI2)) +
      AbsSqr(CpconjVWpbarFdFuPR(gI1,gI2)))*H0(p,MFd(gI1),MFu(gI2)) + 4*B0(p,MFd(
      gI1),MFu(gI2))*MFd(gI1)*MFu(gI2)*Re(Conj(CpconjVWpbarFdFuPL(gI1,gI2))*
      CpconjVWpbarFdFuPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpconjVWpbarFeFvPL(gI1,gI2)) +
      AbsSqr(CpconjVWpbarFeFvPR(gI1,gI2)))*H0(p,MFe(gI1),MFv(gI2)) + 4*B0(p,MFe(
      gI1),MFv(gI2))*MFe(gI1)*MFv(gI2)*Re(Conj(CpconjVWpbarFeFvPL(gI1,gI2))*
      CpconjVWpbarFeFvPR(gI1,gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fd_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
      CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
      CpbarUFdhhFdPR(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),0))*Conj(
      CpbarUFdVGFdPR(gO2,gI2))*CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,
      gI2))*CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2
      ,gI2))*CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
      CpbarUFdconjHpFuPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),MVWp))*Conj(
      CpbarUFdconjVWpFuPR(gO2,gI2))*CpbarUFdconjVWpFuPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPR(gO2,gI1))*
      CpbarUFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPR(gO2,
      gI2))*CpbarUFdconjHpFuPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPL(
      gO2,gI2))*CpbarUFdconjVWpFuPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPR(gO2,gI2))*
      CpbarUFdhhFdPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),0))*Conj(
      CpbarUFdVGFdPL(gO2,gI2))*CpbarUFdVGFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPL(gO2,gI2)
      )*CpbarUFdVPFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPL(gO2,
      gI2))*CpbarUFdVZFdPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
      CpbarUFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,
      gI2))*CpbarUFdconjHpFuPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarUFdconjVWpFuPR(
      gO2,gI2))*CpbarUFdconjVWpFuPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
      CpbarUFdhhFdPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),0))*Conj(
      CpbarUFdVGFdPR(gO2,gI2))*CpbarUFdVGFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2)
      )*CpbarUFdVPFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,
      gI2))*CpbarUFdVZFdPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
      CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
      CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(
      gO2,gI2))*CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
      CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),0))*Conj(
      CpbarUFuVGFuPR(gO2,gI2))*CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,
      gI2))*CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2
      ,gI2))*CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
      CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
      CpbarUFuhhFuPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
      CpbarUFuHpFdPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(
      CpbarUFuVGFuPL(gO2,gI2))*CpbarUFuVGFuPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2)
      )*CpbarUFuVPFuPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPL(gO2,
      gI2))*CpbarUFuVWpFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,
      gI2))*CpbarUFuVZFuPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
      CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
      CpbarUFuhhFuPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
      CpbarUFuHpFdPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(
      CpbarUFuVGFuPR(gO2,gI2))*CpbarUFuVGFuPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2)
      )*CpbarUFuVPFuPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,
      gI2))*CpbarUFuVWpFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,
      gI2))*CpbarUFuVZFuPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
      CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
      CpbarUFehhFePR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,
      gI2))*CpbarUFeVPFePL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2
      ,gI2))*CpbarUFeVZFePL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,2,B0(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
      CpbarUFeconjHpFvPR(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFv(gI2),MVWp))*Conj(
      CpbarUFeconjVWpFvPR(gO2,gI2))*CpbarUFeconjVWpFvPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPR(gO2,gI1))*
      CpbarUFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPR(gO2,
      gI2))*CpbarUFeconjHpFvPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPL(
      gO2,gI2))*CpbarUFeconjVWpFvPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePR(gO2,gI2))*
      CpbarUFehhFePR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePL(gO2,gI2)
      )*CpbarUFeVPFePL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePL(gO2,
      gI2))*CpbarUFeVZFePL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
      CpbarUFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,
      gI2))*CpbarUFeconjHpFvPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarUFeconjVWpFvPR(
      gO2,gI2))*CpbarUFeconjVWpFvPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
      CpbarUFehhFePL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2)
      )*CpbarUFeVPFePR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,
      gI2))*CpbarUFeVZFePR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
      CpbarFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
      CpbarFdhhFdPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,
      gI2))*CpbarFdVZFdPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
      CpbarFdconjHpFuPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),MVWp))*Conj(
      CpbarFdconjVWpFuPR(gO2,gI2))*CpbarFdconjVWpFuPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPR(gO2,gI1))*
      CpbarFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPR(gO2,gI2
      ))*CpbarFdconjHpFuPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPL(
      gO2,gI2))*CpbarFdconjVWpFuPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPR(gO2,gI2))*
      CpbarFdhhFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPL(gO2,gI2
      ))*CpbarFdVZFdPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
      CpbarFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2
      ))*CpbarFdconjHpFuPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVWp))*Conj(CpbarFdconjVWpFuPR(
      gO2,gI2))*CpbarFdconjVWpFuPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
      CpbarFdhhFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2
      ))*CpbarFdVZFdPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
      CpbarFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
      CpbarFehhFePR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,
      gI2))*CpbarFeVZFePL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,2,B0(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
      CpbarFeconjHpFvPR(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFv(gI2),MVWp))*Conj(
      CpbarFeconjVWpFvPR(gO2,gI2))*CpbarFeconjVWpFvPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPR(gO2,gI1))*
      CpbarFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPR(gO2,gI2
      ))*CpbarFeconjHpFvPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPL(
      gO2,gI2))*CpbarFeconjVWpFvPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePR(gO2,gI2))*
      CpbarFehhFePR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePL(gO2,gI2
      ))*CpbarFeVZFePL(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
      CpbarFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2
      ))*CpbarFeconjHpFvPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFv(gI2),MVWp))*Conj(CpbarFeconjVWpFvPR(
      gO2,gI2))*CpbarFeconjVWpFvPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
      CpbarFehhFePL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2
      ))*CpbarFeVZFePR(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
      CpbarFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
      CpbarFuHpFdPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPR(
      gO2,gI2))*CpbarFuVWpFdPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
      CpbarFuhhFuPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,
      gI2))*CpbarFuVPFuPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,
      gI2))*CpbarFuVZFuPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPR(gO2,gI1))*
      CpbarFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPR(gO2,gI2))*
      CpbarFuhhFuPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPR(gO2,gI2))*
      CpbarFuHpFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPL(gO2,gI2))
      *CpbarFuVPFuPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPL(gO2,
      gI2))*CpbarFuVWpFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPL(gO2,gI2
      ))*CpbarFuVZFuPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
      CpbarFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
      CpbarFuhhFuPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
      CpbarFuHpFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))
      *CpbarFuVPFuPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarFuVWpFdPR(gO2,
      gI2))*CpbarFuVWpFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2
      ))*CpbarFuVZFuPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
      CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
      CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(
      gO2,gI2))*CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
      CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,
      gI2))*CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2
      ,gI2))*CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
      CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
      CpbarUFuhhFuPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
      CpbarUFuHpFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2)
      )*CpbarUFuVPFuPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPL(gO2,
      gI2))*CpbarUFuVWpFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,
      gI2))*CpbarUFuVZFuPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
      CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
      CpbarUFuhhFuPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
      CpbarUFuHpFdPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2)
      )*CpbarUFuVPFuPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFd(gI2),MVWp))*Conj(CpbarUFuVWpFdPR(gO2,
      gI2))*CpbarUFuVWpFdPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,
      gI2))*CpbarUFuVZFuPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::tadpole_hh_1loop() const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CphhAhAh();
   result += A0(MVWp)*CphhbargWpCgWpC();
   result += A0(MVWp)*CphhbargWpgWp();
   result += A0(MVZ)*CphhbargZgZ();
   result += -(A0(MHp)*CphhconjHpHp());
   result += -0.5*A0(Mhh)*Cphhhhhh();
   result += 4*A0(MVWp)*CphhconjVWpVWp() - 2*CphhconjVWpVWp()*Sqr(MVWp);
   result += CphhVZVZ()*(2*A0(MVZ) - Sqr(MVZ));
   result += 6*SUM(gI1,0,2,A0(MFd(gI1))*(CphhbarFdFdPL(gI1,gI1) + CphhbarFdFdPR
      (gI1,gI1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(MFe(gI1))*(CphhbarFeFePL(gI1,gI1) + CphhbarFeFePR
      (gI1,gI1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(MFu(gI1))*(CphhbarFuFuPL(gI1,gI1) + CphhbarFuFuPR
      (gI1,gI1))*MFu(gI1));

   return result * oneOver16PiSqr;

}




double Standard_model::self_energy_hh_2loop(double p) const
{
   using namespace flexiblesusy::sm_twoloophiggs;

   const double p2 = Sqr(p);
   const double mb = MFd(2);
   const double mt = MFu(2);
   const double mtau = MFe(2);
   const double yb = Yd(2,2);
   const double yt = Yu(2,2);
   const double ytau = Ye(2,2);
   const double gs = g3;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy -= delta_mh_2loop_at_at_sm(p2, scale, mt, yt, mb);
   }

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy -= delta_mh_2loop_at_as_sm(p2, scale, mt, yt, gs);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy -= delta_mh_2loop_ab_as_sm(p2, scale, mb, yb, gs);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy -= delta_mh_2loop_atau_atau_sm(p2, scale, mtau, ytau);
   }

   return self_energy;
}

double Standard_model::self_energy_hh_3loop() const
{
   using namespace flexiblesusy::sm_threeloophiggs;

   const double mt = MFu(2);
   const double yt = Yu(2,2);
   const double gs = g3;
   const double mh = Mhh;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_3LOOP_CORRECTION_AT_AT_AT) {
      self_energy -= delta_mh_3loop_at_at_at_sm(scale, mt, yt, mh);
   }

   if (HIGGS_3LOOP_CORRECTION_AT_AT_AS) {
      self_energy -= delta_mh_3loop_at_at_as_sm(scale, mt, yt, gs);
   }

   if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
      self_energy -= delta_mh_3loop_at_as_as_sm(scale, mt, yt, gs);
   }

   return self_energy;
}

double Standard_model::self_energy_hh_4loop() const
{
   using namespace flexiblesusy::sm_fourloophiggs;

   const double mt = MFu(2);
   const double yt = Yu(2,2);
   const double gs = g3;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS) {
      self_energy -= delta_mh_4loop_at_as_as_as_sm(scale, mt, yt, gs);
   }

   return self_energy;
}




void Standard_model::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void Standard_model::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void Standard_model::calculate_Mhh_pole()
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const double M_tree(get_mass_matrix_hh());
      const double p = old_Mhh;
      double self_energy = Re(self_energy_hh_1loop(p));
      if (pole_mass_loop_order > 1)
         self_energy += self_energy_hh_2loop(p);
      if (pole_mass_loop_order > 2)
         self_energy += self_energy_hh_3loop();
      if (pole_mass_loop_order > 3)
         self_energy += self_energy_hh_4loop();
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(Mhh) = SignedAbsSqrt(mass_sqr);

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(standard_model_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(standard_model_info::hh);
}

void Standard_model::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void Standard_model::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_VZ());
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void Standard_model::calculate_MFd_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFd(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(self_energy_Fd_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(self_energy_Fd_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(self_energy_Fd_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vd) mix_Vd;
      decltype(Ud) mix_Ud;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_Vd, mix_Ud, eigenvalue_error);
      problems.flag_bad_mass(standard_model_info::Fd, eigenvalue_error > precision
         * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_Vd, mix_Ud);
   #endif
      if (es == 0) {
         PHYSICAL(Vd) = mix_Vd;
         PHYSICAL(Ud) = mix_Ud;
      }
      PHYSICAL(MFd(es)) = Abs(eigen_values(es));
   }
}

void Standard_model::calculate_MFu_pole()
{
   // diagonalization with medium precision
   double qcd_1l = 0.;

   {
      const double currentScale = get_scale();
      qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(2))/Sqr(
         currentScale)))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = -0.005284774766427138*Quad(g3) - 0.0032348537833770956*
         Log(Sqr(currentScale)/Sqr(MFu(2)))*Quad(g3) - 0.0008822328500119351*
         Quad(g3)*Sqr(Log(Power(currentScale,2)/Sqr(MFu(2))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = -0.00003352082872926087*Power6(g3)*(35.702577217116016
         + 1.*Cube(Log(Sqr(currentScale)/Sqr(MFu(2)))) + 15.387410814884797*Log
         (Sqr(currentScale)/Sqr(MFu(2))) + 5.378787878787879*Sqr(Log(Power(
         currentScale,2)/Sqr(MFu(2)))));
   }

   double qcd_4l = 0.;

   if (pole_mass_loop_order > 3 && TOP_POLE_QCD_CORRECTION > 2) {
      const double currentScale = get_scale();
      qcd_4l = -1.6081297554549208e-9*Power8(g3)*(211681.74421123447 - 5638.*
         Cube(Log(Sqr(MFu(2))/Sqr(currentScale))) - 104673.38261571848*Log(Sqr(
         MFu(2))/Sqr(currentScale)) + 825.*Quad(Log(Sqr(MFu(2))/Sqr(
         currentScale))) + 22162.91142653778*Sqr(Log(Sqr(MFu(2))/Sqr(
         currentScale))));
   }

   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      for (int i1 = 0; i1 < 3; ++i1) {
         for (int i2 = 0; i2 < 3; ++i2) {
            if (i1 == 2 && i2 == 2) {
               self_energy_1(i1,i2)  = Re(
                  self_energy_Fu_1loop_1_heavy(p,i1,i2));
               self_energy_PL(i1,i2) = Re(
                  self_energy_Fu_1loop_PL_heavy(p,i1,i2));
               self_energy_PR(i1,i2) = Re(
                  self_energy_Fu_1loop_PR_heavy(p,i1,i2));
            } else {
               self_energy_1(i1,i2)  = Re(self_energy_Fu_1loop_1(p,
                  i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fu_1loop_PL(p
                  ,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fu_1loop_PR(p
                  ,i1,i2));
            }
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l + qcd_4l);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vu) mix_Vu;
      decltype(Uu) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(standard_model_info::Fu, eigenvalue_error > precision
         * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Vu) = mix_Vu;
         PHYSICAL(Uu) = mix_Uu;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void Standard_model::calculate_MFe_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFe(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(self_energy_Fe_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(self_energy_Fe_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(self_energy_Fe_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Ve;
      decltype(Ue) mix_Ue;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_Ve, mix_Ue, eigenvalue_error);
      problems.flag_bad_mass(standard_model_info::Fe, eigenvalue_error > precision
         * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_Ve, mix_Ue);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Ve;
         PHYSICAL(Ue) = mix_Ue;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}

void Standard_model::calculate_MVWp_pole()
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VWp))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_VWp());
   const double p = MVWp;
   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VWp);

   PHYSICAL(MVWp) = AbsSqrt(mass_sqr);
}

double Standard_model::calculate_MVWp_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VWp))
      return 0.;

   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = get_mass_matrix_VWp() - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VWp);

   return AbsSqrt(mass_sqr);
}

double Standard_model::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = get_mass_matrix_VZ() - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VZ);

   return AbsSqrt(mass_sqr);
}


double Standard_model::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double Standard_model::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar + self_energy_1 + m_sm_drbar *
      (self_energy_PL + self_energy_PR);

   return m_susy_drbar;
}

double Standard_model::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_1loop_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   const double qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(idx))
      /Sqr(currentScale)))*Sqr(g3);
   double qcd_2l = 0., qcd_3l = 0., qcd_4l = 0.;

   if (get_thresholds() > 1 && threshold_corrections.mt > 1) {
      const double q_2l = 2.2278607323533713e-6*Quad(g3)*(2372.129769909197 -
         1452.*Log(Sqr(MFu(idx))/Sqr(currentScale)) + 396.*Sqr(Log(Sqr(MFu(idx)
         )/Sqr(currentScale))));

      qcd_2l = -q_2l + qcd_1l * qcd_1l;
   }

   if (get_thresholds() > 2 && threshold_corrections.mt > 2) {
      qcd_3l = 1.0450439184830051e-10*Power6(g3)*(-8.404731799492894e6 + 48600.
         *Cube(Log(Sqr(MFu(idx))/Sqr(currentScale))) + 1.112325741480515e6*Log(
         Sqr(MFu(idx))/Sqr(currentScale)) - 208980.*Sqr(Log(Sqr(MFu(idx))/Sqr(
         currentScale))));
   }

   if (get_thresholds() > 3 && threshold_corrections.mt > 3) {
      qcd_4l = 1.6544544809207004e-13*Power8(g3)*(-1.5015630809970136e9 +
         3.1428e6*Cube(Log(Sqr(MFu(idx))/Sqr(currentScale))) +
         4.409910543513119e8*Log(Sqr(MFu(idx))/Sqr(currentScale)) - 826200.*
         Quad(Log(Sqr(MFu(idx))/Sqr(currentScale))) - 1.7811910793512717e7*Sqr(
         Log(Sqr(MFu(idx))/Sqr(currentScale))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l + qcd_4l);

   return m_susy_drbar;
}

double Standard_model::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree -
      self_energy_PL - self_energy_PR);

   return m_susy_drbar;
}

double Standard_model::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double Standard_model::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_running_tachyon(standard_model_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double Standard_model::calculate_MVWp_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWp_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_running_tachyon(standard_model_info::VWp);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double Standard_model::calculate_Mhh_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_hh_1loop(p));
   const double tadpole = Re(tadpole_hh_1loop());
   const double mass_sqr = Sqr(m_pole) + self_energy - tadpole/v;

   if (mass_sqr < 0.) {
      problems.flag_running_tachyon(standard_model_info::hh);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double Standard_model::ThetaW() const
{
   return ArcCos(Abs(ZZ(0,0)));
}

std::ostream& operator<<(std::ostream& ostr, const Standard_model& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace standard_model

} // namespace flexiblesusy
