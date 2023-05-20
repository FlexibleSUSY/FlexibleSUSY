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

#ifndef CONVERGENCE_TESTER_DRBAR_H
#define CONVERGENCE_TESTER_DRBAR_H

#include "convergence_tester.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"

#include <functional>
#include <limits>

namespace flexiblesusy {

template <class Model>
class Convergence_tester_DRbar : public Convergence_tester {
public:
   using Scale_getter = std::function<double()>;

   Convergence_tester_DRbar(const Model*, double, const Scale_getter& sg = Scale_getter());
   virtual ~Convergence_tester_DRbar() = default;

   virtual bool accuracy_goal_reached() override;
   virtual double get_accuracy_goal() const { return accuracy_goal; }
   virtual int max_iterations() const override { return max_it; }
   virtual void restart() override;

   double get_current_accuracy() const { return current_accuracy; }
   /// set maximum number of iterations
   void set_max_iterations(int it) { max_it = it; }
   /// set scale getter
   template <typename F>
   void set_scale_getter(F&& sg) { scale_getter = std::forward<F>(sg); }

protected:
   /// get current iteration number
   int get_iteration() const { return it_count; }
   /// get model at current iteration
   const Model& get_current_iteration_model() const { return current_model; }
   /// get model state during last iteration
   const Model& get_last_iteration_model() const { return last_iteration_model; }
   /// maximum relative difference to last iteration
   virtual double max_rel_diff() const = 0;

private:
   const Model* model{nullptr};             ///< pointer to model
   Model current_model{};                   ///< model state at current iteration
   Model last_iteration_model{};            ///< model state at last iteration
   Scale_getter scale_getter{};             ///< function to retrieve scale
   int it_count{0};                ///< iteration
   int max_it{40};                 ///< maximum number of iterations
   double accuracy_goal{1e-4};              ///< accuracy goal
   double current_accuracy{std::numeric_limits<double>::infinity()}; ///< current accuracy

   double scale_difference() const;         ///< absolute scale difference
   double rel_scale_difference() const;     ///< relative scale difference
   void run_to_scale();                     ///< runs models to comparison scale
};

template <class Model>
Convergence_tester_DRbar<Model>::Convergence_tester_DRbar
(const Model* model_, double accuracy_goal_, const Scale_getter& sg)
   : Convergence_tester()
   , model(model_)
   , scale_getter(sg)
   , max_it(static_cast<int>(-log10(accuracy_goal_) * 10))
   , accuracy_goal(accuracy_goal_)
{
   if (!model)
      throw SetupError("Convergence_tester_DRbar<Model>: "
                       "model pointer must not be zero!");
}

template <class Model>
bool Convergence_tester_DRbar<Model>::accuracy_goal_reached()
{
   current_model = *model;
   bool precision_reached = false;

   if (it_count > 0) {
      run_to_scale();
      const double scale_accuracy_goal = accuracy_goal * 16*Sqr(Pi);
      if (rel_scale_difference() < scale_accuracy_goal) {
         current_accuracy = max_rel_diff();
         precision_reached = current_accuracy < accuracy_goal;
         VERBOSE_MSG("Convergence_tester_DRbar: current accuracy = "
                     << current_accuracy
                     << ", accuracy goal = " << accuracy_goal);
      } else {
         VERBOSE_MSG("scale has changed by " << scale_difference()
                     << " GeV (" << rel_scale_difference()
                     << "), skipping parameter comparison");
      }
   }

   // save old model parameters
   last_iteration_model = current_model;
   ++it_count;

   return precision_reached;
}

template <class Model>
double Convergence_tester_DRbar<Model>::scale_difference() const
{
   return current_model.get_scale() - last_iteration_model.get_scale();
}

template <class Model>
double Convergence_tester_DRbar<Model>::rel_scale_difference()
   const
{
   const double diff = scale_difference();
   const double last_scale = last_iteration_model.get_scale();
   if (!is_zero(last_scale))
      return diff / last_scale;
   return std::numeric_limits<double>::infinity();
}

template <class Model>
void Convergence_tester_DRbar<Model>::run_to_scale()
{
   if (scale_getter) {
      const double scale = scale_getter();
      current_model.run_to(scale);
      current_model.calculate_DRbar_masses();
      last_iteration_model.run_to(scale);
      last_iteration_model.calculate_DRbar_masses();
   }
}

template <class Model>
void Convergence_tester_DRbar<Model>::restart()
{
   it_count = 0;
   current_model = Model();
   last_iteration_model = Model();
}

} // namespace flexiblesusy

#endif
