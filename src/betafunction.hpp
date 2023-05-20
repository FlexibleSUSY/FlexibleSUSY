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

/**
 * @file betafunction.hpp
 * @brief contains class Beta_function
 */

#ifndef BETAFUNCTION_H
#define BETAFUNCTION_H

#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class Beta_function
 * @brief beta function interface
 *
 * Beta_function is the abstract base class for the beta functions of
 * the parameter classes of the generated models.  It defines the
 * basic RG running interface.  The run() and run_to() functions use
 * the Runge-Kutta algorithm to integrate the RGEs up to a given
 * scale.
 */
class Beta_function {
public:
   Beta_function() = default;
   Beta_function(const Beta_function&) = default;
   Beta_function(Beta_function&&) = default;
   virtual ~Beta_function() = default;
   Beta_function& operator=(const Beta_function&) = default;
   Beta_function& operator=(Beta_function&&) = default;

   void set_scale(double s) { scale = s; }
   void set_number_of_parameters(int pars) { num_pars = pars; }
   void set_loops(int l) { loops = l; }
   void set_thresholds(int t) { thresholds = t; }
   void set_zero_threshold(double t) { zero_threshold = t; }

   double get_scale() const { return scale; }
   int get_number_of_parameters() const { return num_pars; }
   int get_loops() const { return loops; }
   int get_thresholds() const { return thresholds; }
   double get_zero_threshold() const { return zero_threshold; }

   void reset();

   virtual Eigen::ArrayXd get() const = 0;
   virtual void set(const Eigen::ArrayXd&) = 0;
   virtual Eigen::ArrayXd beta() const = 0;

   virtual void run(double, double, double eps = -1.0);
   virtual void run_to(double, double eps = -1.0);

private:
   int num_pars{0};              ///< number of parameters
   int loops{0};                 ///< to what loop order does the RG evolution run
   int thresholds{0};            ///< threshold correction loop order
   double scale{0.};             ///< current renormalization scale
   double tolerance{1.e-4};      ///< running tolerance
   double min_tolerance{1.e-11}; ///< minimum tolerance allowed
   double zero_threshold{1.e-11};///< threshold for treating values as zero

   Eigen::ArrayXd derivatives(double, const Eigen::ArrayXd&);
   double get_tolerance(double eps) const;
};

} // namespace flexiblesusy

#endif
