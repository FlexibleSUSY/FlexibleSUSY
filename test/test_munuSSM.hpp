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

#ifndef TEST_munuSSM_H
#define TEST_munuSSM_H

#include <tuple>

#include "lowe.h"
#include "munuSSM_two_scale_spectrum_generator.hpp"

namespace flexiblesusy {

munuSSM_slha<munuSSM<Two_scale>>
setup_munuSSM(const munuSSM_input_parameters& input, const softsusy::QedQcd& qedqcd)
{
    Spectrum_generator_settings settings;
    settings.set(Spectrum_generator_settings::calculate_sm_masses, 0);
    settings.set(Spectrum_generator_settings::calculate_bsm_masses, 0);

    munuSSM_spectrum_generator<Two_scale> spectrum_generator;
    spectrum_generator.set_settings(settings);
    spectrum_generator.run(qedqcd, input);

    return std::get<0>(spectrum_generator.get_models_slha());
}

} // namespace flexiblesusy

#endif
