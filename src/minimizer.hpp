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

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <utility>
#include <Eigen/Core>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "error.hpp"
#include "ewsb_solver.hpp"
#include "logger.hpp"
#include "gsl_utils.hpp"
#include "gsl_vector.hpp"

namespace flexiblesusy {

/**
 * @class Minimizer
 * @brief Function minimizer
 *
 * The user has to provide the function to be minimized of the type
 * Function_t.  This function gets as arguments an Eigen vector of
 * lenght `dimension' and returns a double.
 *
 * Example:
 * @code
 * auto parabola = [](const Eigen::Matrix<double,2,1>& x) -> double {
 *    const double y = x(0);
 *    const double z = x(1);
 *    return Sqr(y - 5.0) + Sqr(z - 1.0);
 * };
 *
 * Minimizer<2> minimizer(parabola, 100, 1.0e-5);
 * const double start[2] = { 10, 10 };
 * const int status = minimizer.minimize(start);
 * @endcode
 */
template <std::size_t dimension>
class Minimizer : public EWSB_solver {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   using Function_t = std::function<double(const Vector_t&)>;
   enum Solver_type { GSLSimplex, GSLSimplex2, GSLSimplex2Rand };

   Minimizer() = default;
   template <typename F>
   Minimizer(F&&, std::size_t, double, Solver_type solver_type_ = GSLSimplex2);
   virtual ~Minimizer() = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   double get_minimum_value() const { return minimum_value; }
   template <typename F>
   void set_function(F&& f) { function = std::forward<F>(f); }
   void set_precision(double p) { precision = p; }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   void set_solver_type(Solver_type t) { solver_type = t; }
   int minimize(const Vector_t&);

   // EWSB_solver interface methods
   virtual std::string name() const override { return "Minimizer"; }
   virtual int solve(const Eigen::VectorXd&) override;
   virtual Eigen::VectorXd get_solution() const override { return minimum_point; }

private:
   std::size_t max_iterations{100};     ///< maximum number of iterations
   double precision{1.e-2};             ///< precision goal
   double minimum_value{0.};            ///< minimum function value found
   Vector_t minimum_point{Vector_t::Zero()}; ///< vector of minimum point
   Function_t function{nullptr};        ///< function to minimize
   Solver_type solver_type{GSLSimplex2};///< solver type

   void print_state(gsl_multimin_fminimizer*, std::size_t) const;
   static double gsl_function(const gsl_vector*, void*);
   const gsl_multimin_fminimizer_type* solver_type_to_gsl_pointer() const;
};

/**
 * Constructor
 *
 * @param function_ pointer to the function to minimize
 * @param max_iterations_ maximum number of iterations
 * @param precision_ precision goal
 * @param solver_type_ GSL multimin minimizer type
 */
template <std::size_t dimension>
template <typename F>
Minimizer<dimension>::Minimizer(
   F&& function_,
   std::size_t max_iterations_,
   double precision_,
   Solver_type solver_type_
)
   : max_iterations(max_iterations_)
   , precision(precision_)
   , function(std::forward<F>(function_))
   , solver_type(solver_type_)
{
}

/**
 * Start the minimization
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if minimum found)
 */
template <std::size_t dimension>
int Minimizer<dimension>::minimize(const Vector_t& start)
{
   if (!function)
      throw SetupError("Minimizer: function not callable");

   gsl_multimin_fminimizer *minimizer;
   gsl_multimin_function minex_func;

   GSL_vector min_point = to_GSL_vector(start);

   // Set initial step sizes
   GSL_vector step_size(dimension);
   step_size.set_all(1.0);

   // Initialize method and iterate
   minex_func.n = dimension;
   minex_func.f = gsl_function;
   minex_func.params = &function;

   minimizer = gsl_multimin_fminimizer_alloc(solver_type_to_gsl_pointer(), dimension);
   gsl_multimin_fminimizer_set(minimizer, &minex_func, min_point.raw(), step_size.raw());

   size_t iter = 0;
   int status;

   do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(minimizer);

      if (status)
         break;

      const double size = gsl_multimin_fminimizer_size(minimizer);
      status = gsl_multimin_test_size(size, precision);

      print_state(minimizer, iter);
   } while (status == GSL_CONTINUE && iter < max_iterations);

   VERBOSE_MSG("\t\t\tMinimization status = " << gsl_strerror(status));

   // save minimum point and function value
   minimum_point = to_eigen_vector<dimension>(minimizer->x);
   minimum_value = minimizer->fval;

   gsl_multimin_fminimizer_free(minimizer);

   return status;
}

/**
 * Print state of the minimizer
 *
 * @param minimizer minimizer
 * @param iteration iteration number
 */
template <std::size_t dimension>
void Minimizer<dimension>::print_state(gsl_multimin_fminimizer* minimizer,
                                               std::size_t iteration) const
{
   VERBOSE_MSG("\t\t\tIteration " << iteration
               << ": x = " << GSL_vector(minimizer->x)
               << ", f(x) = " << minimizer->fval);
}

template <std::size_t dimension>
int Minimizer<dimension>::solve(const Eigen::VectorXd& start)
{
   return (minimize(start) == GSL_SUCCESS ?
           EWSB_solver::SUCCESS : EWSB_solver::FAIL);
}

template <std::size_t dimension>
double Minimizer<dimension>::gsl_function(const gsl_vector* x, void* params)
{
   if (!is_finite(x))
      return std::numeric_limits<double>::max();

   Function_t* fun = static_cast<Function_t*>(params);
   const Vector_t arg(to_eigen_vector<dimension>(x));
   double result = std::numeric_limits<double>::max();

   try {
      result = (*fun)(arg);
   } catch (const flexiblesusy::Error&) {
   }

   return result;
}

template <std::size_t dimension>
const gsl_multimin_fminimizer_type* Minimizer<dimension>::solver_type_to_gsl_pointer() const
{
   switch (solver_type) {
   case GSLSimplex     : return gsl_multimin_fminimizer_nmsimplex;
   case GSLSimplex2    : return gsl_multimin_fminimizer_nmsimplex2;
   case GSLSimplex2Rand: return gsl_multimin_fminimizer_nmsimplex2rand;
   default:
      throw SetupError("Unknown minimizer solver type: "
                       + std::to_string(solver_type));
   }

   return nullptr;
}

} // namespace flexiblesusy

#endif
