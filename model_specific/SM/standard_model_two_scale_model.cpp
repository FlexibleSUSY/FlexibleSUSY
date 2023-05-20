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
 * @file standard_model_two_scale_model.cpp
 * @brief implementation of the SM model class
 *
 * Contains the definition of the SM model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 */

#include "standard_model_two_scale_model.hpp"
#include "standard_model.hpp"

namespace flexiblesusy {
namespace standard_model {

using namespace standard_model_info;


StandardModel<Two_scale>::StandardModel()
   : Model()
   , Standard_model()
{
}

void StandardModel<Two_scale>::calculate_spectrum()
{
   Standard_model::calculate_spectrum();
}

void StandardModel<Two_scale>::clear_problems()
{
   Standard_model::clear_problems();
}

std::string StandardModel<Two_scale>::name() const
{
   return Standard_model::name();
}

void StandardModel<Two_scale>::run_to(double scale, double eps)
{
   Standard_model::run_to(scale, eps);
}

void StandardModel<Two_scale>::print(std::ostream& out) const
{
   Standard_model::print(out);
}

void StandardModel<Two_scale>::set_precision(double p)
{
   Standard_model::set_precision(p);
}

} // namespace standard_model
} // namespace flexiblesusy
