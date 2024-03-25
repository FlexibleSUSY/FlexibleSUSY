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
 * @file experimental_constraints.hpp
 *
 * @brief contains Higgs effective coupling calculation and interfaces to
 * HiggsTools and Lilith
 */

#ifndef HIGGSTOOLS_INTERFACE_H
#define HIGGSTOOLS_INTERFACE_H

#include "config.h"

#include "decays/decay.hpp"
#include "physical_input.hpp"
#include "lowe.h"
#include "spectrum_generator_settings.hpp"
#include "decays/flexibledecay_settings.hpp"

#include <string>
#include <vector>
#include <optional>
#include <utility>

namespace flexiblesusy {

EffectiveCoupling_list get_normalized_effective_couplings(
  EffectiveCoupling_list const& bsm_input,
  Physical_input const& physical_input,
  softsusy::QedQcd const& qedqcd,
  Spectrum_generator_settings const& spectrum_generator_settings,
  FlexibleDecay_settings const& flexibledecay_settings);

struct SignalResult {
   std::size_t ndof {};
   double mhRef {};
   double chi2BSM {};
   // SM chi2 for Physical_input::mh_pole
   double chi2SM {};
};

#ifdef ENABLE_HIGGSTOOLS
std::tuple<SignalResult, std::vector<std::tuple<int, double, double, std::string>>> call_higgstools(
   EffectiveCoupling_list const&,
   Physical_input const&,
   std::string const& higgsbounds_dataset, std::string const& higgssignals_dataset);
#endif

/**
 *  @return pair containit -2 log likelihood and number of degrees of freedom
 */
#ifdef ENABLE_LILITH
std::optional<SignalResult>
call_lilith(EffectiveCoupling_list const&, Physical_input const&, std::string const&);
#endif

} // flexiblesusy

#endif
