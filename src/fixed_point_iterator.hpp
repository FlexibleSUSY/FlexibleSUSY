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

#ifndef FIXED_POINT_ITERATOR_H
#define FIXED_POINT_ITERATOR_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <Eigen/Core>

#include "logger.hpp"
#include "eigen_utils.hpp"
#include "error.hpp"
#include "ewsb_solver.hpp"

namespace flexiblesusy {

namespace fixed_point_iterator {

template <std::size_t dimension>
class Convergence_tester_absolute {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   enum Status : int { SUCCESS, CONTINUE };

   explicit Convergence_tester_absolute(double precision_ = 1.0e-2)
      : precision(precision_)
   {}

   const char* name() const { return "Convergence_tester_absolute"; }

   /**
    * Test whether the absolute value of the residual, defined by
    * \f$|a-b| = \sqrt{\sum_i (a_i - b_i)^2}\f$,
    * is less than the set precision.
    *
    * @param a first vector
    * @param b second vector
    * @return status code (SUCCESS or CONTINUE)
    */
   int operator()(const Vector_t& a, const Vector_t& b) const {
      const double res = std::sqrt((a - b).array().square().sum());
      return res < precision ? SUCCESS : CONTINUE;
   }

private:
   double precision;                 ///< precision goal
};

template <std::size_t dimension>
class Convergence_tester_relative {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   enum Status : int { SUCCESS, CONTINUE };

   explicit Convergence_tester_relative(double precision_ = 1.0e-2)
      : precision(precision_)
   {}

   const char* name() const { return "Convergence_tester_relative"; }

   /**
    * Test whether the relative difference is less than the set
    * precision. The relative difference test used here is carried out
    * by applying \a is_equal_rel() to each element of the vector.
    *
    * @param a first vector
    * @param b second vector
    * @return status code (SUCCESS or CONTINUE)
    */
   int operator()(const Vector_t& a, const Vector_t& b) const {
      return is_equal_rel(a, b, precision) ? SUCCESS : CONTINUE;
   }

private:
   double precision;                 ///< precision goal
};

template <std::size_t dimension>
class Convergence_tester_tadpole {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   using Function_t = std::function<Vector_t(const Vector_t&)>;
   enum Status : int { SUCCESS, CONTINUE };

   Convergence_tester_tadpole(double precision_,
                              const Function_t& tadpole_function_)
      : precision(precision_)
      , tadpole_function(tadpole_function_)
   {}

   const char* name() const { return "Convergence_tester_tadpole"; }

   /**
    * Test whether the relative difference is less than the set
    * precision. The relative difference test used here is carried out
    * by applying \a is_equal_rel() to each element of the vector. If the
    * relative difference is below the precision, it is tested whether
    * the tadpoles are below the precision. If the tadpoles are larger
    * than the precision, CONTINUE is returned.
    *
    * @param a first vector
    * @param b second vector
    * @return status code (SUCCESS or CONTINUE)
    */
   int operator()(const Vector_t& a, const Vector_t& b) const
   {
      if (!is_equal_rel(a, b, precision)) {
         return CONTINUE;
      }

      static const double eps = 10*std::pow(10., -std::numeric_limits<double>::digits10);

      if (is_equal_rel(a, b, eps)) {
         return SUCCESS;
      }

      return check_tadpoles(a);
   }

private:
   double precision;                 ///< precision goal
   const Function_t tadpole_function; ///< function to calculate tadpole

   int check_tadpoles(const Vector_t& x) const {
      const double res = x.cwiseAbs().sum();
      return res < precision ? SUCCESS : CONTINUE;
   }
};

} // namespace fixed_point_iterator

/**
 * @class Fixed_point_iterator
 * @brief Does fixed point iteration
 * @author Dylan Harries, Alexander Voigt
 * @tparam dimension dimension of function
 * @tparam Convergence_tester function for relative comparison
 *    of subsequent iteration steps
 *
 * The user has to provide the function (of which a fixed point should
 * be found) of the type \a Function_t. This function gets as
 * arguments a Eigen vector of length \a dimension and returns a
 * vector with the next point.
 *
 * @note The standard relative convergence criterion may not be
 * suitable in all situations, for example not when the iteration
 * converges slowly.  In this case subsequent steps are very close to
 * each other, but \f$x_n\f$ might not be close to the true fixed
 * point.
 *
 * @todo implement check for no progress towards solution
 */
template <std::size_t dimension, class Convergence_tester = fixed_point_iterator::Convergence_tester_relative<dimension>>
class Fixed_point_iterator : public EWSB_solver {
public:
   using Vector_t = Eigen::Matrix<double,dimension,1>;
   using Function_t = std::function<Vector_t(const Vector_t&)>;

   Fixed_point_iterator() = default;
   Fixed_point_iterator(const Function_t&, std::size_t, const Convergence_tester&);
   virtual ~Fixed_point_iterator() = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   void set_function(const Function_t& f) { function = f; }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   int find_fixed_point(const Eigen::VectorXd&);

   // EWSB_solver interface methods
   virtual std::string name() const override { return std::string("Fixed_point_iterator<") + convergence_tester.name() + ">"; }
   virtual int solve(const Eigen::VectorXd& start) override { return find_fixed_point(start); }
   virtual Eigen::VectorXd get_solution() const override { return fixed_point; }

private:
   std::size_t max_iterations{100};         ///< maximum number of iterations
   Vector_t xn{Vector_t::Zero()};           ///< current iteration point
   Vector_t fixed_point{Vector_t::Zero()};  ///< vector of fixed point estimate
   Function_t function{nullptr};            ///< function defining fixed point
   Convergence_tester convergence_tester{}; ///< convergence tester

   int fixed_point_iterator_iterate();
   void print_state(std::size_t) const;

   static bool is_finite(const Vector_t& v) {
      return std::any_of(v.data(), v.data() + v.size(),
                         [](double x) { return std::isfinite(x); });
   }
};

/**
 * Constructor
 *
 * @param function_ pointer to the function to find fixed point for
 * @param max_iterations_ maximum number of iterations
 * @param convergence_tester_ convergence tester
 */
template <std::size_t dimension, class Convergence_tester>
Fixed_point_iterator<dimension,Convergence_tester>::Fixed_point_iterator(
   const Function_t& function_,
   std::size_t max_iterations_,
   const Convergence_tester& convergence_tester_
)
   : max_iterations(max_iterations_)
   , function(function_)
   , convergence_tester(convergence_tester_)
{
}

/**
 * Start the iteration
 *
 * @param start starting point
 *
 * @return status code (SUCCESS if fixed point found)
 */
template <std::size_t dimension, class Convergence_tester>
int Fixed_point_iterator<dimension,Convergence_tester>::find_fixed_point(
   const Eigen::VectorXd& start
)
{
   if (!function) {
      throw SetupError("Fixed_point_iterator: function not callable");
   }

   int status;
   std::size_t iter = 0;

   fixed_point = xn = start;

   print_state(iter);

   do {
      iter++;
      status = fixed_point_iterator_iterate();

      print_state(iter);

      if (status == EWSB_solver::FAIL) {
         break;
      }

      status = convergence_tester(fixed_point, xn);

   } while (status == Convergence_tester::CONTINUE && iter < max_iterations);

   VERBOSE_MSG("\t\t\tFixed_point_iterator status = " << status);

   return status;
}

/**
 * Perform a single step of the fixed point iteration
 *
 * @return status
 */
template <std::size_t dimension, class Convergence_tester>
int Fixed_point_iterator<dimension,Convergence_tester>::fixed_point_iterator_iterate()
{
   int status = EWSB_solver::SUCCESS;
   xn = fixed_point;

   try {
      fixed_point = function(xn);
   } catch (const flexiblesusy::Error&) {
      status = EWSB_solver::FAIL;
   }

   // For safety, include a check for nans or infs
   if (!is_finite(fixed_point)) {
      status = EWSB_solver::FAIL;
   }

   return status;
}

/**
 * Print state of the fixed point iterator
 *
 * @param iteration iteration number
 */
template <std::size_t dimension, class Convergence_tester>
void Fixed_point_iterator<dimension,Convergence_tester>::print_state(std::size_t iteration) const
{
   VERBOSE_MSG("\t\t\tIteration n = " << iteration
               << ": x_{n} = [" << xn.transpose()
               << "], x_{n+1} = [" << fixed_point.transpose() << "]");
}

} // namespace flexiblesusy

#endif
