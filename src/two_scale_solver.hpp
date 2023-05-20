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

#ifndef TWO_SCALE_SOLVER_H
#define TWO_SCALE_SOLVER_H

#include "rg_flow.hpp"

#include <memory>
#include <vector>

/**
 * @file two_scale_solver.hpp
 * @brief contains the definition of the RGFlow<Two_scale> class
 */

namespace flexiblesusy {

class Convergence_tester;
class Initial_guesser;
class Model;
class Single_scale_constraint;
class Single_scale_matching;

class Two_scale;
class Two_scale_running_precision;

/**
 * @class RGFlow<Two_scale>
 * @brief Boundary condition solver (two-scale algorithm)
 *
 * This boundary condition solver uses the two-scale algorithm to
 * solve the boundary value problem: It uses RG running to iteratively
 * run the models to the boundary condition (constraint) scales and
 * imposes the constraints.
 *
 * To add constraints use the add() function.  Matching conditions are
 * added using the add_upwards() or add_downwards() functions,
 * depending on whether the low-scale model should be matched to the
 * high-scale one (add_upwards()) or vice versa (add_downwards()).
 * The added constraints and matching conditions are applied in their
 * given order.
 */

template<>
class RGFlow<Two_scale> {
public:
   RGFlow() = default;
   RGFlow(const RGFlow&) = delete;
   RGFlow(RGFlow&&) = delete;
   ~RGFlow() = default;
   RGFlow& operator=(const RGFlow&) = delete;
   RGFlow& operator=(RGFlow&&) = delete;

   /// add constraint
   void add(Single_scale_constraint*, Model*);
   /// add matching condition
   void add(Single_scale_matching*, Model*, Model*);
   /// get model at current scale
   Model* get_model() const;
   /// get number of used iterations
   int number_of_iterations_done() const;
   /// clear all internal data
   void reset();
   /// run model at given scale to given scale
   void run_to(double);
   /// set convergence tester
   void set_convergence_tester(Convergence_tester*);
   /// set running precision calculator
   void set_running_precision(Two_scale_running_precision*);
   /// set initial guesser
   void set_initial_guesser(Initial_guesser*);
   /// solves the boundary value problem
   void solve();

private:
   struct Slider {
   public:
      virtual ~Slider() {}
      virtual void clear_problems() {}
      virtual Model* get_model() = 0;
      virtual double get_scale() = 0;
      virtual void slide() {}
      virtual void set_precision(double) {}
   };

   struct Constraint_slider : public Slider {
   public:
      Constraint_slider(Model* m, Single_scale_constraint* c)
         : model(m), constraint(c) {}
      virtual ~Constraint_slider() {}
      virtual void clear_problems() override;
      virtual Model* get_model() override;
      virtual double get_scale() override;
      virtual void slide() override;
      virtual void set_precision(double) override;
   private:
      Model* model;
      Single_scale_constraint* constraint;
   };

   struct Matching_slider : public Slider {
   public:
      Matching_slider(Model* m1_, Model* m2_, Single_scale_matching* mc)
         : m1(m1_), m2(m2_), matching(mc) {}
      virtual ~Matching_slider() {}
      virtual void clear_problems() override;
      virtual Model* get_model() override;
      virtual double get_scale() override;
      virtual void slide() override;
      virtual void set_precision(double) override;
   private:
      Model *m1, *m2;
      Single_scale_matching* matching;
   };

   std::vector<std::shared_ptr<Slider> > sliders{}; ///< sliders to be run up and down
   int iteration{0};             ///< iteration number (starting at 0)
   Convergence_tester* convergence_tester{nullptr}; ///< the convergence tester
   Initial_guesser* initial_guesser{nullptr};       ///< does initial guess
   Two_scale_running_precision* running_precision_calculator{nullptr}; ///< RG running precision calculator
   double running_precision{1.0e-3};           ///< RG running precision
   double scale{0.};                           ///< current scale

   bool accuracy_goal_reached() const; ///< check if accuracy goal is reached
   void check_setup() const;           ///< check the setup
   void clear_problems();              ///< clear model problems
   int get_max_iterations() const; ///< returns max. number of iterations
   Model* get_model(double) const;     ///< returns model at given scale
   double get_precision();             ///< returns running precision
   void initial_guess();               ///< initial guess
   void run_sliders();                 ///< run all sliders
   std::vector<std::shared_ptr<Slider> > sort_sliders() const; ///< sort the sliders w.r.t. to scale
   void update_running_precision();    ///< update the RG running precision
};

} // namespace flexiblesusy

#endif
