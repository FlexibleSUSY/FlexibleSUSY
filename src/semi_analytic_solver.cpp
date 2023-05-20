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

#include "semi_analytic_solver.hpp"

#include "convergence_tester.hpp"
#include "error.hpp"
#include "initial_guesser.hpp"
#include "logger.hpp"
#include "model.hpp"
#include "single_scale_constraint.hpp"
#include "single_scale_matching.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

#include <algorithm>

/**
 * @file semi_analytic_solver.cpp
 * @brief contains the implementation of the RGFlow<Semi_analytic> class members
 */

namespace flexiblesusy {

/**
 * Add a model constraint for the inner iteration
 *
 * @param c constraint
 * @param m model
 */
void RGFlow<Semi_analytic>::add_inner(Single_scale_constraint* c, Model* m)
{
   if (!c) throw SetupError("constraint pointer is NULL");
   if (!m) throw SetupError("model pointer is NULL");
   inner_sliders.push_back(std::make_shared<Constraint_slider>(m, c));
}

/**
 * Adds a matching condition for the inner iteration.  This matching
 * condition matches the two models by calling match().  The two
 * models are passed to the set_models() function of the matching
 * condition.
 *
 * @param mc matching condition
 * @param m1 model 1
 * @param m2 model 2
 */
void RGFlow<Semi_analytic>::add_inner(Single_scale_matching* mc, Model* m1, Model* m2)
{
   if (!mc) throw SetupError("matching condition pointer is NULL");
   if (!m1) throw SetupError("model pointer 1 is NULL");
   if (!m2) throw SetupError("model pointer 2 is NULL");
   mc->set_models(m1, m2);
   inner_sliders.push_back(std::make_shared<Matching_slider>(m1, m2, mc));
}

/**
 * Add a model constraint for the outer iteration
 *
 * @param c constraint
 * @param m model
 */
void RGFlow<Semi_analytic>::add_outer(Single_scale_constraint* c, Model* m)
{
   if (!c) throw SetupError("constraint pointer is NULL");
   if (!m) throw SetupError("model pointer is NULL");
   outer_sliders.push_back(std::make_shared<Constraint_slider>(m, c));
}

/**
 * Adds a matching condition for the outer iteration.  This matching
 * condition matches the two models by calling match().  The two
 * models are passed to the set_models() function of the matching
 * condition.
 *
 * @param mc matching condition
 * @param m1 model 1
 * @param m2 model 2
 */
void RGFlow<Semi_analytic>::add_outer(Single_scale_matching* mc, Model* m1, Model* m2)
{
   if (!mc) throw SetupError("matching condition pointer is NULL");
   if (!m1) throw SetupError("model pointer 1 is NULL");
   if (!m2) throw SetupError("model pointer 2 is NULL");
   mc->set_models(m1, m2);
   outer_sliders.push_back(std::make_shared<Matching_slider>(m1, m2, mc));
}

/**
 * @brief Solves the boundary value problem.
 *
 * At first the initial_guess() is called.  Afterwards, the semianalytic
 * solution is applied: first, the function
 * iteratively runs the tower up and down and imposes the boundary
 * conditions.  The iteration stops if either the maximum number of
 * iterations is reached or the precision goal is achieved (defined by
 * the convergence_tester for the inner two-scale iteration). Then
 * the next, semianalytic step is applied. This is iterated until
 * either the maximum number of iterations is reached or the precision goal
 * is achieved (defined by the convergence_tester for the outer iteration).
 */
void RGFlow<Semi_analytic>::solve()
{
   check_setup();

   const int max_iterations = get_max_iterations();
   if ((inner_sliders.empty() && outer_sliders.empty())
       || max_iterations == 0)
      return;

   initial_guess();

   iteration = 0;
   bool accuracy_reached = false;
   // flag for non-convergence in iteration
   bool no_conv_inner_iteration = false;
   RGFlow<Two_scale> inner_solver;
   while (iteration < max_iterations && !accuracy_reached) {
      update_running_precision();
      // reset problem flags on each iteration step
      clear_problems();
      no_conv_inner_iteration = false;
      prepare_inner_iteration(inner_solver);
      try {
         inner_solver.solve();
      } catch (const NoConvergenceError&) {
         no_conv_inner_iteration = true;
      }
      run_outer_sliders();
      accuracy_reached = accuracy_goal_reached();
      ++iteration;
   }

   if (!accuracy_reached || no_conv_inner_iteration)
      throw NoConvergenceError(max_iterations);

   VERBOSE_MSG("convergence reached after " << iteration << " iterations");
}

/**
 * Sanity checks the models and boundary conditions
 */
void RGFlow<Semi_analytic>::check_setup() const
{
   if (!inner_convergence_tester) {
      throw SetupError("RGFlow<Semi_analytic>::Error: inner convergence tester must "
                       "not be NULL");
   }
   if (!outer_convergence_tester) {
      throw SetupError("RGFlow<Semi_analytic>::Error: outer convergence tester must "
                       "not be NULL");
   }
}

void RGFlow<Semi_analytic>::clear_problems()
{
   VERBOSE_MSG("> clearing problems ...");

   for (auto& s: outer_sliders)
      s->clear_problems();
}

/**
 * Initializes the two-scale solver used in the inner iteration.
 */
void RGFlow<Semi_analytic>::prepare_inner_iteration(RGFlow<Two_scale>& solver) const
{
   inner_convergence_tester->restart();

   solver.reset();
   solver.set_convergence_tester(inner_convergence_tester);
   solver.set_running_precision(running_precision_calculator);

   for (auto& s: inner_sliders) {
      s->clear_problems();
      s->add_constraint_to_solver(solver);
   }
}

/**
 * Does the initial guess by calling the guess() method of the initial
 * guesser (if given).
 */
void RGFlow<Semi_analytic>::initial_guess()
{
   if (initial_guesser)
      initial_guesser->guess();
}

void RGFlow<Semi_analytic>::run_outer_sliders()
{
   VERBOSE_MSG("> running all models (outer iteration " << iteration << ") ...");

   for (auto& s: outer_sliders) {
      s->set_precision(get_precision());
      s->slide();
   }

   VERBOSE_MSG("> running outer sliders finished");
}

/**
 * Returns the precision of the RG running.
 *
 * @return RG running precision
 */
double RGFlow<Semi_analytic>::get_precision()
{
   return running_precision;
}

/**
 * Recalculates the precision of the RG running using the user defined
 * Two_scale_running_precision_calculator class.
 */
void RGFlow<Semi_analytic>::update_running_precision()
{
   if (running_precision_calculator)
      running_precision = running_precision_calculator->get_precision(iteration);
}

/**
 * Returns the value returned by the accuracy_goal_reached() method of
 * the convergence tester.
 */
bool RGFlow<Semi_analytic>::accuracy_goal_reached() const
{
   return outer_convergence_tester->accuracy_goal_reached();
}

/**
 * Set the convergence tester to be used during the inner iteration.
 *
 * @param convergence_tester_ the convergence tester to be used
 */
void RGFlow<Semi_analytic>::set_inner_convergence_tester(Convergence_tester* convergence_tester_)
{
   inner_convergence_tester = convergence_tester_;
}

/**
 * Set the convergence tester to be used during the outer iteration.
 *
 * @param convergence_tester_ the convergence tester to be used
 */
void RGFlow<Semi_analytic>::set_outer_convergence_tester(Convergence_tester* convergence_tester_)
{
   outer_convergence_tester = convergence_tester_;
}

void RGFlow<Semi_analytic>::set_initial_guesser(Initial_guesser* ig)
{
   initial_guesser = ig;
}

/**
 * Set RG running precision calculator.
 *
 * @param rp running precision calculator
 */
void RGFlow<Semi_analytic>::set_running_precision(Two_scale_running_precision* rp)
{
   running_precision_calculator = rp;
}

/**
 * Returns the number of performed iterations
 * @return number of performed iterations
 */
int RGFlow<Semi_analytic>::number_of_iterations_done() const
{
   return iteration;
}

/**
 * Returns the maximum number of iterations set in the outer
 * convergence tester.
 */
int RGFlow<Semi_analytic>::get_max_iterations() const
{
   return outer_convergence_tester->max_iterations();
}

/**
 * Returns the pointer to the model at the given scale.
 *
 * @param scale scale for which corresponding model to return
 * @return model at scale
 */
Model* RGFlow<Semi_analytic>::get_model(double scale) const
{
   const std::vector<std::shared_ptr<Slider> > sorted_inner_sliders(sort_sliders(inner_sliders));
   const std::vector<std::shared_ptr<Slider> > sorted_outer_sliders(sort_sliders(outer_sliders));
   std::vector<std::shared_ptr<Slider> > sorted_sliders;
   std::merge(sorted_inner_sliders.begin(), sorted_inner_sliders.end(),
              sorted_outer_sliders.begin(), sorted_outer_sliders.end(),
              std::back_inserter(sorted_sliders),
              [](const std::shared_ptr<Slider>& s1,
                 const std::shared_ptr<Slider>& s2)
              { return s1->get_scale() < s2->get_scale(); });

   auto it = std::lower_bound(sorted_sliders.begin(), sorted_sliders.end(),
                              scale,
                              [](const std::shared_ptr<Slider>& s, double scale)
                              { return s->get_scale() < scale; });

   if (it == sorted_sliders.end())
      return nullptr;

   return (*it)->get_model();
}

/**
 * Returns the pointer to the model at the current scale.
 * @return model at current scale
 */
Model* RGFlow<Semi_analytic>::get_model() const
{
   return get_model(scale);
}

/**
 * @brief resets the solver to the initial condition
 *
 * The pointers to the models, matching conditions, convergence
 * testers, initial guesser, and running precision calculator are set
 * to zero.  The running precision is set to the default value 0.001.
 */
void RGFlow<Semi_analytic>::reset()
{
   inner_sliders.clear();
   outer_sliders.clear();

   iteration = 0;
   inner_convergence_tester = nullptr;
   outer_convergence_tester = nullptr;
   initial_guesser = nullptr;
   running_precision_calculator = nullptr;
   running_precision = 1.0e-3;
   scale = 0;
}

/**
 * Returns vector of sliders, sorted w.r.t. their scale.
 *
 * @return vector of sorted sliders
 */
std::vector<std::shared_ptr<RGFlow<Semi_analytic>::Slider> > RGFlow<Semi_analytic>::sort_sliders(const std::vector<std::shared_ptr<RGFlow<Semi_analytic>::Slider> >& sliders) const
{
   std::vector<std::shared_ptr<Slider> > sorted_sliders(sliders);

   std::sort(sorted_sliders.begin(), sorted_sliders.end(),
             [](const std::shared_ptr<Slider>& s1, const std::shared_ptr<Slider>& s2)
             { return s1->get_scale() < s2->get_scale(); });

   return sorted_sliders;
}

/**
 * Run the model tower to the given scale.
 *
 * @param scale_ scale to run to
 */
void RGFlow<Semi_analytic>::run_to(double scale_)
{
   scale = scale_;

   Model* model = get_model();

   if (model)
      model->run_to(scale);
}

/* Implementation of sliders */

void RGFlow<Semi_analytic>::Constraint_slider::clear_problems() {
   model->clear_problems();
}

Model* RGFlow<Semi_analytic>::Constraint_slider::get_model() {
   return model;
}

void RGFlow<Semi_analytic>::Constraint_slider::add_constraint_to_solver(
   RGFlow<Two_scale>& solver) {
   solver.add(constraint, model);
}

double RGFlow<Semi_analytic>::Constraint_slider::get_scale() {
   return constraint->get_scale();
}

void RGFlow<Semi_analytic>::Constraint_slider::slide() {
   VERBOSE_MSG("> \trunning " << model->name() << " to scale " << constraint->get_scale() << " GeV");
   model->run_to(constraint->get_scale());
   VERBOSE_MSG("> \tapplying " << constraint->name());
   constraint->apply();
}

void RGFlow<Semi_analytic>::Constraint_slider::set_precision(double p) {
   model->set_precision(p);
}

void RGFlow<Semi_analytic>::Matching_slider::clear_problems() {
   m1->clear_problems();
   m2->clear_problems();
}

Model* RGFlow<Semi_analytic>::Matching_slider::get_model() {
   return m1;
}

void RGFlow<Semi_analytic>::Matching_slider::add_constraint_to_solver(
   RGFlow<Two_scale>& solver) {
   solver.add(matching, m1, m2);
}

double RGFlow<Semi_analytic>::Matching_slider::get_scale() {
   return matching->get_scale();
}

void RGFlow<Semi_analytic>::Matching_slider::slide() {
   VERBOSE_MSG("> \trunning " << m1->name() << " to scale " << matching->get_scale() << " GeV");
   m1->run_to(matching->get_scale());
   VERBOSE_MSG("> \trunning " << m2->name() << " to scale " << matching->get_scale() << " GeV");
   m2->run_to(matching->get_scale());
   VERBOSE_MSG("> \tmatching " << m1->name() << " -> " << m2->name());
   matching->match();
}

void RGFlow<Semi_analytic>::Matching_slider::set_precision(double p) {
   m1->set_precision(p);
   m2->set_precision(p);
}

} // namespace flexiblesusy
