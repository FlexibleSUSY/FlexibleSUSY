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

#ifndef SHOOTING_SOLVER_H
#define SHOOTING_SOLVER_H

#include "composite_root_finder.hpp"
#include "logger.hpp"
#include "model.hpp"
#include "single_scale_constraint.hpp"
#include "single_scale_matching.hpp"
#include <cstddef>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class Shooting_solver
 * @brief Solves the IVP using the shooting method
 */
template<std::size_t N>
class Shooting_solver {
public:
   using Vec_t = Eigen::Matrix<double, N, 1>;
   using Fun_t = std::function<Vec_t(const Vec_t&)>;
   using Pre_t = std::function<Vec_t()>;
   using Set_t = std::function<void(const Vec_t&)>;

   /// add constraint
   void add(Single_scale_constraint*, Model*);
   /// add matching condition
   void add(Single_scale_matching*, Model*, Model*);
   /// add predictor
   void add(const Pre_t& p) { predictor = p; }
   /// add parameter setter
   void add(const Set_t& s) { setter = s; }
   /// returns solution
   Vec_t get_solution() const { return solution; }
   /// set precision goal
   void set_precision(double p) { precision = p; }
   /// set maximum number of iterations
   void set_max_iterations(std::size_t n) { max_it = n; }
   /// clear all internal data
   void reset();
   /// solves the boundary value problem with initial guess
   void solve(const Vec_t&);

private:
   struct Slider {
   public:
      virtual ~Slider() {}
      virtual void clear_problems() {}
      virtual double get_scale() = 0;
      virtual void slide() {}
   };

   struct Constraint_slider : public Slider {
   public:
      Constraint_slider(Model* m, Single_scale_constraint* c)
         : model(m), constraint(c) {}
      virtual ~Constraint_slider() = default;
      virtual void clear_problems() override { model->clear_problems(); }
      virtual double get_scale() override { return constraint->get_scale(); }
      virtual void slide() override {
         const auto Q = get_scale();
         VERBOSE_MSG("> \trunning " << model->name() << " to scale " << Q << " GeV");
         model->run_to(Q);
         VERBOSE_MSG("> \tapplying " << constraint->name());
         constraint->apply();
      }
   private:
      Model* model;
      Single_scale_constraint* constraint;
   };

   struct Matching_slider : public Slider {
   public:
      Matching_slider(Model* m1_, Model* m2_, Single_scale_matching* mc)
         : m1(m1_), m2(m2_), matching(mc) {}
      virtual ~Matching_slider() = default;
      virtual void clear_problems() override {
         m1->clear_problems();
         m2->clear_problems();
      }
      virtual double get_scale() override { return matching->get_scale(); }
      virtual void slide() override {
         const auto Q = get_scale();
         VERBOSE_MSG("> \trunning " << m1->name() << " to scale " << Q << " GeV");
         m1->run_to(Q);
         VERBOSE_MSG("> \trunning " << m2->name() << " to scale " << Q << " GeV");
         m2->run_to(Q);
         VERBOSE_MSG("> \tmatching " << m1->name() << " -> " << m2->name());
         matching->match();
      }
   private:
      Model *m1, *m2;
      Single_scale_matching* matching;
   };

   std::vector<std::shared_ptr<Slider> > sliders{}; ///< sliders to be run
   Vec_t solution{Vec_t::Zero()}; ///< stored solution
   double precision{1e-3};        ///< running precision
   std::size_t max_it{100};       ///< maximum number of iterations
   Pre_t predictor{nullptr};      ///< calculates prediction
   Set_t setter{nullptr};         ///< sets parameters
};

/**
 * Adding a model constraint
 *
 * @param c constraint
 * @param m model
 */
template<std::size_t N>
void Shooting_solver<N>::add(Single_scale_constraint* c, Model* m)
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
template<std::size_t N>
void Shooting_solver<N>::add(Single_scale_matching* mc, Model* m1, Model* m2)
{
   if (!mc) throw SetupError("matching condition pointer is NULL");
   if (!m1) throw SetupError("model pointer 1 is NULL");
   if (!m2) throw SetupError("model pointer 2 is NULL");
   mc->set_models(m1, m2);
   sliders.push_back(std::make_shared<Matching_slider>(m1, m2, mc));
}

/**
 * Reset solver
 */
template<std::size_t N>
void Shooting_solver<N>::reset()
{
   sliders.clear();
   solution = Vec_t::Zero();
   precision = 1e-3;
   max_it = 100;
   predictor = nullptr;
   setter = nullptr;
}

/**
 * Solves the IVP by finding the root of the function defined by the
 * following procdure:
 *
 * 1. setting the parameters (using the setter)
 * 2. imposing all boundary and matching conditions
 * 3. calculating the output (calling the predictor)
 *
 * @param v0 initial parameter guess
 */
template<std::size_t N>
void Shooting_solver<N>::solve(const Vec_t& v0)
{
   if (!predictor)
      throw SetupError("No predictor set!");
   if (!setter)
      throw SetupError("No setter set!");

   Fun_t fun = [&] (const Vec_t& v) {
      setter(v);
      for (auto& s: sliders) {
         s->slide();
      }
      return predictor();
   };

   Composite_root_finder<N> crs;
   crs.set_precision(precision);
   crs.set_max_iterations(max_it);
   crs.solve(fun, v0);

   solution = crs.get_solution();
}

} // namespace flexiblesusy

#endif
