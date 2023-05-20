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

#ifndef SEMI_ANALYTIC_SOLVER_H
#define SEMI_ANALYTIC_SOLVER_H

#include "rg_flow.hpp"

#include <memory>
#include <vector>

/**
 * @file semi_analytic_solver.hpp
 * @brief contains the definition of the RGFlow<Semi_analytic> class
 */

namespace flexiblesusy {

class Convergence_tester;
class Initial_guesser;
class Model;
class Single_scale_constraint;
class Single_scale_matching;

class Semi_analytic;
class Two_scale;
class Two_scale_running_precision;

/**
 * @class RGFlow<Semi_analytic>
 * @brief Boundary condition solver (semi-analytic algorithm)
 *
 * This boundary condition solver uses the semi-analytic algorithm to
 * solve the boundary value problem in two stages.  Firstly, a
 * two-scale iteration is performed to solve for a subset of the
 * parameters.  Once these are determined, semi-analytic solutions
 * for the remaining parameters are calculated at the relevant
 * scales.
 */

template<>
class RGFlow<Semi_analytic> {
public:
   /// Create empty semi-analytic solver
   /// The RG running precision is set to the default value of 0.001
   RGFlow() = default;
   RGFlow(const RGFlow&) = delete;
   RGFlow(RGFlow&&) = delete;
   ~RGFlow() = default;
   RGFlow& operator=(const RGFlow&) = delete;
   RGFlow& operator=(RGFlow&&) = delete;

   /// add inner constraint
   void add_inner(Single_scale_constraint*, Model*);
   /// add inner matching condition
   void add_inner(Single_scale_matching*, Model*, Model*);
   /// add outer constraint
   void add_outer(Single_scale_constraint*, Model*);
   /// add outer matching condition
   void add_outer(Single_scale_matching*, Model*, Model*);
   /// get model at current scale
   Model* get_model() const;
   /// get number of used iterations
   int number_of_iterations_done() const;
   /// clear all internal data
   void reset();
   /// run model at given scale to given scale
   void run_to(double);
   /// set convergence tester for inner two-scale iteration
   void set_inner_convergence_tester(Convergence_tester*);
   /// set convergence tester for overall iteration
   void set_outer_convergence_tester(Convergence_tester*);
   /// set running precision calculator
   void set_running_precision(Two_scale_running_precision*);
   /// set initial guesser
   void set_initial_guesser(Initial_guesser*);
   /// solves the boundary value problem
   void solve();

private:
   struct Slider {
   public:
      virtual ~Slider() = default;
      virtual void clear_problems() {}
      virtual Model* get_model() = 0;
      virtual void add_constraint_to_solver(RGFlow<Two_scale>&) = 0;
      virtual double get_scale() = 0;
      virtual void slide() {}
      virtual void set_precision(double) {}
   };

   struct Constraint_slider : public Slider {
   public:
      Constraint_slider(Model* m, Single_scale_constraint* c)
         : model(m), constraint(c) {}
      virtual ~Constraint_slider() = default;
      virtual void clear_problems() override;
      virtual Model* get_model() override;
      virtual void add_constraint_to_solver(RGFlow<Two_scale>&) override;
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
      virtual ~Matching_slider() = default;
      virtual void clear_problems() override;
      virtual Model* get_model() override;
      virtual void add_constraint_to_solver(RGFlow<Two_scale>&) override;
      virtual double get_scale() override;
      virtual void slide() override;
      virtual void set_precision(double) override;
   private:
      Model *m1, *m2;
      Single_scale_matching* matching;
   };

   std::vector<std::shared_ptr<Slider> > inner_sliders{}; ///< sliders to be used in the inner iteration
   std::vector<std::shared_ptr<Slider> > outer_sliders{}; ///< sliders to be used in the outer iteration
   int iteration{0};             ///< iteration number (starting at 0)
   Convergence_tester* inner_convergence_tester{nullptr}; ///< the convergence tester for the two-scale iteration
   Convergence_tester* outer_convergence_tester{nullptr}; ///< the convergence tester for the main iteration
   Initial_guesser* initial_guesser{nullptr};       ///< does initial guess
   Two_scale_running_precision* running_precision_calculator{nullptr}; ///< RG running precision calculator
   double running_precision{1.0e-3};           ///< RG running precision
   double scale{0.};                           ///< current scale

   bool accuracy_goal_reached() const; ///< check if accuracy goal is reached
   void check_setup() const;           ///< check the setup
   void clear_problems();              ///< clear model problems
   void prepare_inner_iteration(RGFlow<Two_scale>& solver) const; ///< set-up inner two-scale iteration
   int get_max_iterations() const; ///< returns max. number of iterations
   Model* get_model(double) const;     ///< returns model at given scale
   void initial_guess();               ///< initial guess
   double get_precision();             ///< returns running precision
   void update_running_precision();    ///< update the RG running precision
   std::vector<std::shared_ptr<Slider> > sort_sliders(
      const std::vector<std::shared_ptr<Slider> >&) const; ///< sort the outer sliders w.r.t. scale

   void run_outer_sliders();           ///< run all outer iteration sliders
};

} // namespace flexiblesusy

#endif
