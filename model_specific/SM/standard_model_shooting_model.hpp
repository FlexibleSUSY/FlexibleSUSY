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
 * @file standard_model_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 */

#ifndef STANDARD_MODEL_SHOOTING_MODEL_H
#define STANDARD_MODEL_SHOOTING_MODEL_H

#include "model.hpp"
#include "standard_model.hpp"
#include "standard_model_two_scale_model.hpp"

namespace flexiblesusy {

class Shooting;
class Two_scale;

/**
 * @class StandardModel<Shooting>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */

namespace standard_model {

template<>
class StandardModel<Shooting> : public StandardModel<Two_scale> {};

} // namespace standard_model

} // namespace flexiblesusy

#endif
