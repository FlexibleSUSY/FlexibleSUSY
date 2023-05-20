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

#ifndef MSSMCBS_SPECTRUM_GENERATOR_H
#define MSSMCBS_SPECTRUM_GENERATOR_H

#include "CMSSM_info.hpp"
#include "MSSMcbs_two_scale_model.hpp"
#include "CMSSM_two_scale_ewsb_solver.hpp"
#include "CMSSM_two_scale_high_scale_constraint.hpp"
#include "CMSSM_two_scale_susy_scale_constraint.hpp"
#include "MSSMcbs_two_scale_low_scale_constraint.hpp"
#include "CMSSM_two_scale_convergence_tester.hpp"
#include "MSSMcbs_two_scale_initial_guesser.hpp"
#include "CMSSM_utilities.hpp"

#include "coupling_monitor.hpp"
#include "error.hpp"
#include "numerics2.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

namespace flexiblesusy {

template <class T>
class MSSMcbs_spectrum_generator {
public:
   MSSMcbs_spectrum_generator()
      : solver(), model()
      , high_scale_constraint()
      , susy_scale_constraint()
      , low_scale_constraint()
      , high_scale(0.)
      , susy_scale(0.)
      , low_scale(0.)
      , parameter_output_scale(0.)
      , precision_goal(1.0e-4)
      , max_iterations(0)
      , beta_loop_order(2)
      , threshold_corrections_loop_order(1)
      , calculate_sm_masses(false) {}
   ~MSSMcbs_spectrum_generator() {}

   double get_high_scale() const { return high_scale; }
   double get_susy_scale() const { return susy_scale; }
   double get_low_scale()  const { return low_scale;  }
   const MSSMcbs<T>& get_model() const { return model; }
   const Problems& get_problems() const { return model.get_problems(); }
   int get_exit_code() const { return get_problems().have_problem(); };
   void set_parameter_output_scale(double s) { parameter_output_scale = s; }
   void set_precision_goal(double precision_goal_) { precision_goal = precision_goal_; }
   void set_pole_mass_loop_order(int l) { model.set_pole_mass_loop_order(l); }
   void set_ewsb_loop_order(int l) { model.set_ewsb_loop_order(l); }
   void set_beta_loop_order(int l) { beta_loop_order = l; }
   void set_max_iterations(int n) { max_iterations = n; }
   void set_calculate_sm_masses(bool flag) { calculate_sm_masses = flag; }
   void set_threshold_corrections_loop_order(int t) { threshold_corrections_loop_order = t; }

   void run(const softsusy::QedQcd& qedqcd, const CMSSM_input_parameters& input);
   void write_running_couplings(const std::string& filename = "CMSSM_rge_running.dat") const;
   void write_spectrum(const std::string& filename = "CMSSM_spectrum.dat") const;

private:
   RGFlow<T> solver;
   MSSMcbs<T> model;
   CMSSM_high_scale_constraint<T> high_scale_constraint;
   CMSSM_susy_scale_constraint<T> susy_scale_constraint;
   MSSMcbs_low_scale_constraint<T> low_scale_constraint;
   double high_scale, susy_scale, low_scale;
   double parameter_output_scale; ///< output scale for running parameters
   double precision_goal; ///< precision goal
   int max_iterations; ///< maximum number of iterations
   int beta_loop_order; ///< beta-function loop order
   int threshold_corrections_loop_order; ///< threshold corrections loop order
   bool calculate_sm_masses; ///< calculate SM pole masses
};

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param qedqcd Standard Model input parameters
 * @param input model input parameters
 */
template <class T>
void MSSMcbs_spectrum_generator<T>::run(const softsusy::QedQcd& qedqcd,
					const CMSSM_input_parameters& input)
{
   model.clear();
   model.set_input_parameters(input);
   model.do_calculate_sm_pole_masses(calculate_sm_masses);
   model.set_loops(beta_loop_order);
   model.set_thresholds(threshold_corrections_loop_order);

   CMSSM_ewsb_solver<T> ewsb_solver;
   model.set_ewsb_solver(
      std::make_shared<CMSSM_ewsb_solver<T> >(ewsb_solver));

   high_scale_constraint.clear();
   susy_scale_constraint.clear();
   low_scale_constraint .clear();
   high_scale_constraint.set_model(&model);
   susy_scale_constraint.set_model(&model);
   low_scale_constraint .set_model(&model);
   low_scale_constraint .set_sm_parameters(qedqcd);
   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint .initialize();

   CMSSM_convergence_tester<T> convergence_tester(&model, precision_goal);
   if (max_iterations > 0)
      convergence_tester.set_max_iterations(max_iterations);

   MSSMcbs_initial_guesser<T> initial_guesser(&model, qedqcd,
                                                  low_scale_constraint,
                                                  susy_scale_constraint,
                                                  high_scale_constraint);

   Two_scale_increasing_precision precision(10.0, precision_goal);

   solver.reset();
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add(&low_scale_constraint, &model);
   solver.add(&high_scale_constraint, &model);
   solver.add(&susy_scale_constraint, &model);

   high_scale = susy_scale = low_scale = 0.;

   try {
      solver.solve();
      high_scale = high_scale_constraint.get_scale();
      susy_scale = susy_scale_constraint.get_scale();
      low_scale  = low_scale_constraint.get_scale();

      model.run_to(susy_scale);
      model.solve_ewsb();
      model.calculate_spectrum();

      // run to output scale (if scale > 0)
      if (!is_zero(parameter_output_scale)) {
         model.run_to(parameter_output_scale);
      }
   } catch (const NoConvergenceError&) {
      model.get_problems().flag_thrown("no convergence");
   } catch (const NonPerturbativeRunningError&) {
      model.get_problems().flag_no_perturbative();
   } catch (const Error& error) {
      model.get_problems().flag_thrown();
   } catch (const std::string& str) {
      model.get_problems().flag_thrown();
   } catch (const char* str) {
      model.get_problems().flag_thrown();
   } catch (const std::exception& error) {
      model.get_problems().flag_thrown();
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
template <class T>
void MSSMcbs_spectrum_generator<T>::write_running_couplings(const std::string& filename) const
{
   MSSMcbs<T> tmp_model(model);
   tmp_model.run_to(low_scale);

   // returns parameters at given scale
   auto data_getter = [&tmp_model](double scale) {
      tmp_model.run_to(scale);
      return CMSSM_parameter_getter::get_parameters(tmp_model);
   };

   std::vector<std::string> parameter_names(
      std::cbegin(CMSSM_info::parameter_names),
      std::cend(CMSSM_info::parameter_names));

   Coupling_monitor coupling_monitor(data_getter, parameter_names);
   coupling_monitor.run(low_scale, high_scale, 100, true);
   coupling_monitor.write_to_file(filename);
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template <class T>
void MSSMcbs_spectrum_generator<T>::write_spectrum(const std::string& filename) const
{
   CMSSM_spectrum_plotter plotter;
   plotter.extract_spectrum(model);
   plotter.write_to_file(filename);
}

} // namespace flexiblesusy

#endif
