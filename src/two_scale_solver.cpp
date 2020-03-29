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

#include "two_scale_solver.hpp"

#include "convergence_tester.hpp"
#include "error.hpp"
#include "initial_guesser.hpp"
#include "logger.hpp"
#include "model.hpp"
#include "single_scale_constraint.hpp"
#include "single_scale_matching.hpp"
#include "two_scale_running_precision.hpp"

#include <algorithm>

/**
 * @file two_scale_solver.cpp
 * @brief contains the implementation of the RGFlow<Two_scale> class members
 */

namespace flexiblesusy {

/**
 * Adding a model constraint
 *
 * @param c constraint
 * @param m model
 */
void RGFlow<Two_scale>::add(Single_scale_constraint* c, Model* m)
{
   if (!c) throw SetupError("constraint pointer is NULL");
   if (!m) throw SetupError("model pointer is NULL");
   sliders.push_back(std::make_shared<Constraint_slider>(m, c));
}

/**
 * Adds a matching condition.  This matching condition matches the two
 * models by calling match().  The two models are passed to the
 * set_models() function of the matching condition.
 *
 * @param mc matching condition
 * @param m1 model 1
 * @param m2 model 2
 */
void RGFlow<Two_scale>::add(Single_scale_matching* mc, Model* m1, Model* m2)
{
   if (!mc) throw SetupError("matching condition pointer is NULL");
   if (!m1) throw SetupError("model pointer 1 is NULL");
   if (!m2) throw SetupError("model pointer 2 is NULL");
   mc->set_models(m1, m2);
   sliders.push_back(std::make_shared<Matching_slider>(m1, m2, mc));
}

/**
 * @brief Solves the boundary value problem.
 *
 * At first the initial_guess() is called.  Afterwards, the function
 * iteratively runs the tower up and down and imposes the boundary
 * conditions.  The iteration stops if either the maximum number of
 * iterations is reached or the precision goal is achieved (defined by
 * the convergence_tester).
 */
void RGFlow<Two_scale>::solve()
{
   check_setup();

   const int max_iterations = get_max_iterations();
   if (sliders.empty() || max_iterations == 0)
      return;

   initial_guess();

   iteration = 0;
   bool accuracy_reached = false;

   while (iteration < max_iterations && !accuracy_reached) {
      update_running_precision();
      clear_problems();
      run_sliders();
      accuracy_reached = accuracy_goal_reached();
      ++iteration;
   }

   if (!accuracy_reached)
      throw NoConvergenceError(max_iterations);

   VERBOSE_MSG("convergence reached after " << iteration << " iterations");
}

/**
 * Sanity checks the models and boundary conditions.
 */
void RGFlow<Two_scale>::check_setup() const
{
   if (!convergence_tester) {
      throw SetupError("RGFlow<Two_scale>::Error: convergence tester must "
                       "not be NULL");
   }
}

void RGFlow<Two_scale>::clear_problems()
{
   VERBOSE_MSG("> clearing problems ...");

   for (auto& s: sliders)
      s->clear_problems();
}

/**
 * Does the initial guess by calling the guess() method of the initial
 * guesser (if given).
 */
void RGFlow<Two_scale>::initial_guess()
{
   if (initial_guesser)
      initial_guesser->guess();
}

void RGFlow<Two_scale>::run_sliders()
{
   VERBOSE_MSG("> running all models (iteration " << iteration << ") ...");

   for (auto& s: sliders) {
      s->set_precision(get_precision());
      s->slide();
   }

   VERBOSE_MSG("> running sliders finished");
}

/**
 * Returns the precision of the RG running.
 *
 * @return RG running precision
 */
double RGFlow<Two_scale>::get_precision()
{
   return running_precision;
}

/**
 * Recalculates the precision of the RG running using the user defined
 * Two_scale_running_precision_calculator class.
 */
void RGFlow<Two_scale>::update_running_precision()
{
   if (running_precision_calculator)
      running_precision = running_precision_calculator->get_precision(iteration);
}

/**
 * Returns the value returned by the accuracy_goal_reached() method of
 * the convergence tester.
 */
bool RGFlow<Two_scale>::accuracy_goal_reached() const
{
   return convergence_tester->accuracy_goal_reached();
}

/**
 * Set the convergence tester to be used during the iteration.
 *
 * @param convergence_tester_ the convergence tester to be used
 */
void RGFlow<Two_scale>::set_convergence_tester(Convergence_tester* convergence_tester_)
{
   convergence_tester = convergence_tester_;
}

void RGFlow<Two_scale>::set_initial_guesser(Initial_guesser* ig)
{
   initial_guesser = ig;
}

/**
 * Set RG running precision calculator.
 *
 * @param rp running precision calculator
 */
void RGFlow<Two_scale>::set_running_precision(Two_scale_running_precision* rp)
{
   running_precision_calculator = rp;
}

/**
 * Returns the number of performed iterations
 * @return number of performed iterations
 */
int RGFlow<Two_scale>::number_of_iterations_done() const
{
   return iteration;
}

/**
 * Returns the maximum number of iterations set in the convergence
 * tester.
 */
int RGFlow<Two_scale>::get_max_iterations() const
{
   return convergence_tester->max_iterations();
}

/**
 * Returns the pointer to the model at the given scale.
 *
 * @param scale scale for which corresponding model to return
 * @return model at scale
 */
Model* RGFlow<Two_scale>::get_model(double scale) const
{
   const auto sorted_sliders = sort_sliders();

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
Model* RGFlow<Two_scale>::get_model() const
{
   return get_model(scale);
}

/**
 * @brief resets the solver to the initial condition
 *
 * The pointers to the models, matching conditions, convergence
 * tester, initial guesser, and running precision calculator are set
 * to zero.  The running precision is set to the default value 0.001.
 */
void RGFlow<Two_scale>::reset()
{
   sliders.clear();

   iteration = 0;
   convergence_tester = nullptr;
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
std::vector<std::shared_ptr<RGFlow<Two_scale>::Slider> > RGFlow<Two_scale>::sort_sliders() const
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
void RGFlow<Two_scale>::run_to(double scale_)
{
   scale = scale_;

   Model* model = get_model(scale_);

   if (model)
      model->run_to(scale_);
}

/* Implementation of sliders */

void RGFlow<Two_scale>::Constraint_slider::clear_problems() {
   model->clear_problems();
}

Model* RGFlow<Two_scale>::Constraint_slider::get_model() {
   return model;
}

double RGFlow<Two_scale>::Constraint_slider::get_scale() {
   return constraint->get_scale();
}

void RGFlow<Two_scale>::Constraint_slider::slide() {
   VERBOSE_MSG("> \trunning " << model->name() << " to scale " << constraint->get_scale() << " GeV");
   model->run_to(constraint->get_scale());
   VERBOSE_MSG("> \tapplying " << constraint->name());
   constraint->apply();
}

void RGFlow<Two_scale>::Constraint_slider::set_precision(double p) {
   model->set_precision(p);
}

void RGFlow<Two_scale>::Matching_slider::clear_problems() {
   m1->clear_problems();
   m2->clear_problems();
}

Model* RGFlow<Two_scale>::Matching_slider::get_model() {
   return m1;
}

double RGFlow<Two_scale>::Matching_slider::get_scale() {
   return matching->get_scale();
}

void RGFlow<Two_scale>::Matching_slider::slide() {
   VERBOSE_MSG("> \trunning " << m1->name() << " to scale " << matching->get_scale() << " GeV");
   m1->run_to(matching->get_scale());
   VERBOSE_MSG("> \trunning " << m2->name() << " to scale " << matching->get_scale() << " GeV");
   m2->run_to(matching->get_scale());
   VERBOSE_MSG("> \tmatching " << m1->name() << " -> " << m2->name());
   matching->match();
}

void RGFlow<Two_scale>::Matching_slider::set_precision(double p) {
   m1->set_precision(p);
   m2->set_precision(p);
}

} // namespace flexiblesusy
