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

#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <Eigen/Core>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "logger.hpp"
#include "error.hpp"
#include "ewsb_solver.hpp"
#include "gsl_utils.hpp"
#include "gsl_vector.hpp"

namespace flexiblesusy {

/**
 * @class Root_finder
 * @brief Function root finder
 *
 * The user has to provide the function (of which the root should be
 * found) of the type Function_t.  This function gets as arguments an
 * Eigen vector of lenght `dimension' and returns an Eigen vector of
 * the same length.
 *
 * Example:
 * @code
 * auto parabola = [](const Eigen::Matrix<double,2,1>& x) {
 *    const double y = x(0);
 *    const double z = x(1);
 *    Eigen::Matrix<double,2,1> f;
 *    f << y*(y - 5.0), z*(z - 1.0);
 *    return f;
 * };
 *
 * Root_finder<2> root_finder(parabola, 100, 1.0e-5);
 * const double start[2] = { 10, 10 };
 * const int status = root_finder.find_root(start);
 * @endcode
 */
template <std::size_t dimension>
class Root_finder : public EWSB_solver {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   using Function_t = std::function<Vector_t(const Vector_t&)>;
   enum Solver_type { GSLHybrid, GSLHybridS, GSLBroyden, GSLNewton };

   Root_finder() = default;
   Root_finder(const Function_t&, std::size_t, double, Solver_type solver_type_ = GSLHybrid);
   virtual ~Root_finder() = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   void set_function(const Function_t& f) { function = f; }
   void set_precision(double p) { precision = p; }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   void set_solver_type(Solver_type t) { solver_type = t; }
   int find_root(const Vector_t&);

   // EWSB_solver interface methods
   virtual std::string name() const override { return std::string("Root_finder<") + solver_type_name() + ">"; }
   virtual int solve(const Eigen::VectorXd&) override;
   virtual Eigen::VectorXd get_solution() const override { return root; }

private:
   std::size_t max_iterations{100};    ///< maximum number of iterations
   double precision{1.e-2};            ///< precision goal
   Vector_t root{Vector_t::Zero()};    ///< the root
   Function_t function{nullptr};       ///< function to minimize
   Solver_type solver_type{GSLHybrid}; ///< solver type

   void print_state(const gsl_multiroot_fsolver*, std::size_t) const;
   const char* solver_type_name() const;
   const gsl_multiroot_fsolver_type* solver_type_to_gsl_pointer() const;
   static int gsl_function(const gsl_vector*, void*, gsl_vector*);

   static bool is_finite(const Vector_t& v) {
      return std::any_of(v.data(), v.data() + v.size(),
                         [](double x) { return std::isfinite(x); });
   }
};

/**
 * Constructor
 *
 * @param function_ pointer to the function to minimize
 * @param max_iterations_ maximum number of iterations
 * @param precision_ precision goal
 * @param solver_type_ GSL multiroot solver type
 */
template <std::size_t dimension>
Root_finder<dimension>::Root_finder(
   const Function_t& function_,
   std::size_t max_iterations_,
   double precision_,
   Solver_type solver_type_
)
   : max_iterations(max_iterations_)
   , precision(precision_)
   , function(function_)
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
int Root_finder<dimension>::find_root(const Vector_t& start)
{
   if (!function)
      throw SetupError("Root_finder: function not callable");

   int status;
   std::size_t iter = 0;
   void* parameters = &function;
   gsl_multiroot_function f = {gsl_function, dimension, parameters};

   gsl_multiroot_fsolver* solver
      = gsl_multiroot_fsolver_alloc(solver_type_to_gsl_pointer(), dimension);

   if (!solver) {
      throw OutOfMemoryError(std::string("Cannot allocate gsl_multiroot_fsolver ") +
                             gsl_multiroot_fsolver_name(solver));
   }

#ifndef ENABLE_DEBUG
   gsl_set_error_handler_off();
#endif

   GSL_vector tmp_root = to_GSL_vector(start);

   gsl_multiroot_fsolver_set(solver, &f, tmp_root.raw());

   print_state(solver, iter);

   do {
      iter++;
      status = gsl_multiroot_fsolver_iterate(solver);

      print_state(solver, iter);

      if (status)   // check if solver is stuck
         break;

      status = gsl_multiroot_test_residual(solver->f, precision);
   } while (status == GSL_CONTINUE && iter < max_iterations);

   VERBOSE_MSG("\t\t\tRoot_finder status = " << gsl_strerror(status));

   root = to_eigen_vector<dimension>(solver->x);

   gsl_multiroot_fsolver_free(solver);

   return status;
}

/**
 * Print state of the root finder
 *
 * @param solver solver
 * @param iteration iteration number
 */
template <std::size_t dimension>
void Root_finder<dimension>::print_state(const gsl_multiroot_fsolver* solver,
                                         std::size_t iteration) const
{
   VERBOSE_MSG("\t\t\tIteration " << iteration
               << ": x = " << GSL_vector(solver->x)
               << ", f(x) = " << GSL_vector(solver->f));
}

template <std::size_t dimension>
int Root_finder<dimension>::gsl_function(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!flexiblesusy::is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   Function_t* fun = static_cast<Function_t*>(params);
   int status = GSL_SUCCESS;
   const Vector_t arg(to_eigen_vector<dimension>(x));
   Vector_t result;
   result.setConstant(std::numeric_limits<double>::max());

   try {
      result = (*fun)(arg);
      status = is_finite(result) ? GSL_SUCCESS : GSL_EDOM;
   } catch (const flexiblesusy::Error&) {
      status = GSL_EDOM;
   }

   // copy result -> f
   for (std::size_t i = 0; i < dimension; i++) {
      gsl_vector_set(f, i, result(i));
   }

   return status;
}

template <std::size_t dimension>
const char* Root_finder<dimension>::solver_type_name() const
{
   switch (solver_type) {
   case GSLHybrid : return "GSLHybrid";
   case GSLHybridS: return "GSLHybridS";
   case GSLBroyden: return "GSLBroyden";
   case GSLNewton : return "GSLNewton";
   default:
      throw SetupError("Unknown root solver type: "
                       + std::to_string(solver_type));
   }

   return "unknown";
}

template <std::size_t dimension>
const gsl_multiroot_fsolver_type* Root_finder<dimension>::solver_type_to_gsl_pointer() const
{
   switch (solver_type) {
   case GSLHybrid : return gsl_multiroot_fsolver_hybrid;
   case GSLHybridS: return gsl_multiroot_fsolver_hybrids;
   case GSLBroyden: return gsl_multiroot_fsolver_broyden;
   case GSLNewton : return gsl_multiroot_fsolver_dnewton;
   default:
      throw SetupError("Unknown root solver type: "
                       + std::to_string(solver_type));
   }

   return nullptr;
}

template <std::size_t dimension>
int Root_finder<dimension>::solve(const Eigen::VectorXd& start)
{
   return (find_root(start) == GSL_SUCCESS ?
           EWSB_solver::SUCCESS : EWSB_solver::FAIL);
}

} // namespace flexiblesusy

#endif
