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
 * @file HiggsTools_interface.hpp
 *
 * @brief contains interface to HiggsTools
 */

#ifndef HIGGSTOOLS_INTERFACE_H
#define HIGGSTOOLS_INTERFACE_H

#include "decays/decay.hpp"
#include "physical_input.hpp"
#include "lowe.h"
#include "spectrum_generator_settings.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

void call_HiggsTools(
   EffectiveCoupling_list const&,
   Physical_input const&,
   softsusy::QedQcd const&,
   Spectrum_generator_settings const&,
   FlexibleDecay_settings const&);

} // flexiblesusy

#endif
