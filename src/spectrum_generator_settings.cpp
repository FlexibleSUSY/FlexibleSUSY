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

#include "spectrum_generator_settings.hpp"
#include "error.hpp"
#include "string_format.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace flexiblesusy {

namespace {
const std::array<std::string, Spectrum_generator_settings::NUMBER_OF_OPTIONS> descriptions = {
   "precision goal",
   "max. iterations (0 = automatic)",
   "solver (0 = all)",
   "calculate SM pole masses",
   "pole mass loop order",
   "EWSB loop order",
   "beta-functions loop order",
   "threshold corrections loop order",
   "Higgs 2-loop corrections O(alpha_t alpha_s)",
   "Higgs 2-loop corrections O(alpha_b alpha_s)",
   "Higgs 2-loop corrections O((alpha_t + alpha_b)^2)",
   "Higgs 2-loop corrections O(alpha_tau^2)",
   "force output",
   "Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)",
   "beta-function zero threshold",
   "calculate observables (a_muon, ...)",
   "force positive majorana masses",
   "pole mass renormalization scale (0 = SUSY scale)",
   "pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))",
   "EFT matching scale (0 = SUSY scale)",
   "EFT loop order for upwards matching",
   "EFT loop order for downwards matching",
   "EFT index of SM-like Higgs in the BSM model",
   "calculate BSM pole masses",
   "individual threshold correction loop orders",
   "Renormalization scheme for Higgs 3-loop corrections O(alpha_t alpha_s^2 + alpha_b alpha_s^2)",
   "Higgs 3-loop corrections O(alpha_t alpha_s^2)",
   "Higgs 3-loop corrections O(alpha_b alpha_s^2)",
   "Higgs 3-loop corrections O(alpha_t^2 alpha_s)",
   "Higgs 3-loop corrections O(alpha_t^3)",
   "Higgs 4-loop corrections O(alpha_t alpha_s^3)",
   "loop library type (0 = Softsusy)"
};

bool is_integer(double value)
{
   double intpart;
   return std::modf(value, &intpart) == 0.0;
}

void assert_bool(double value, const char* quantity)
{
   if (value != 0.0 && value != 1.0) {
      throw SetupError(std::string(quantity) + " must either 0 or 1");
   }
}

void assert_integer(double value, const char* quantity)
{
   if (!is_integer(value)) {
      throw SetupError(std::string(quantity) + " must be an integer");
   }
}

void assert_ge(double value, double lower_bound, const char* quantity)
{
   if (value < lower_bound) {
      throw SetupError(std::string(quantity) +
                       " must be greater than or equal to " +
                       flexiblesusy::to_string(lower_bound));
   }
}

void assert_gt(double value, double lower_bound, const char* quantity)
{
   if (value <= lower_bound) {
      throw SetupError(std::string(quantity) + " must be greater than " +
                       flexiblesusy::to_string(lower_bound));
   }
}

void assert_le(double value, double upper_bound, const char* quantity)
{
   if (value > upper_bound) {
      throw SetupError(std::string(quantity) +
                       " must be lower than or equal to " +
                       flexiblesusy::to_string(upper_bound));
   }
}

} // anonymous namespace

/**
 * Default constructor
 *
 * Calls reset() to initialize all spectrum generator settings to
 * their default values.
 */
Spectrum_generator_settings::Spectrum_generator_settings()
{
   reset();
}

double Spectrum_generator_settings::get(Settings o) const
{
   return values.at(o);
}

Spectrum_generator_settings::Settings_t Spectrum_generator_settings::get() const
{
   Settings_t s(&values[0]);
   return s;
}

std::string Spectrum_generator_settings::get_description(Settings o) const
{
   return descriptions.at(o);
}

void Spectrum_generator_settings::set(Settings o, double value)
{
   switch (o) {
   case precision: // 0 [double > 0]
      assert_gt(value, 0.0, descriptions.at(o).c_str());
      break;
   case max_iterations: // 1 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case solver: // 2 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case calculate_sm_masses: // 3 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case pole_mass_loop_order: // 4 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case ewsb_loop_order: // 5 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case beta_loop_order: // 6 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case threshold_corrections_loop_order: // 7 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case higgs_2loop_correction_at_as: // 8 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_2loop_correction_ab_as: // 9 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_2loop_correction_at_at: // 10 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_2loop_correction_atau_atau: // 11 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case force_output: // 12 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case top_pole_qcd_corrections: // 13 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case beta_zero_threshold: // 14 [double > 0]
      assert_ge(value, 0.0, descriptions.at(o).c_str());
      break;
   case calculate_observables: // 15 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case force_positive_masses: // 16 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case pole_mass_scale: // 17 [double >= 0]
      assert_ge(value, 0.0, descriptions.at(o).c_str());
      break;
   case eft_pole_mass_scale: // 18 [double >= 0]
      assert_ge(value, 0.0, descriptions.at(o).c_str());
      break;
   case eft_matching_scale: // 19 [double >= 0]
      assert_ge(value, 0.0, descriptions.at(o).c_str());
      break;
   case eft_matching_loop_order_up: // 20 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case eft_matching_loop_order_down: // 21 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case eft_higgs_index: // 22 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case calculate_bsm_masses: // 23 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case threshold_corrections: // 24 [int >= 0]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      break;
   case higgs_3loop_ren_scheme_atb_as2: // 25 [int >= 0 and <= 2]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      assert_le(value, 2, descriptions.at(o).c_str());
      break;
   case higgs_3loop_correction_at_as2: // 26 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_3loop_correction_ab_as2: // 27 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_3loop_correction_at2_as: // 28 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_3loop_correction_at3: // 29 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case higgs_4loop_correction_at_as3: // 30 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   default:
      break;
   }

   values.at(o) = value;
}

void Spectrum_generator_settings::set(const Spectrum_generator_settings::Settings_t& s)
{
   std::copy(s.data(), s.data() + s.size(), values.begin());
}

/**
 * Resets all spectrum generator settings to their defaults.
 *
 * | enum                             | possible values                                 | default value   |
 * |----------------------------------|-------------------------------------------------|-----------------|
 * | precision                        | any positive double                             | 1.0e-4          |
 * | max_iterations                   | any positive double                             | 0 (= automatic) |
 * | solver                           | 0 (all), 1 (two-scale) or 2 (semi-analytic)     | 0 (= all)       |
 * | calculate_sm_masses              | 0 (no) or 1 (yes)                               | 0 (= no)        |
 * | pole_mass_loop_order             | 0, 1, 2, 3, 4                                   | 4 (= 4-loop)    |
 * | ewsb_loop_order                  | 0, 1, 2, 3, 4                                   | 4 (= 4-loop)    |
 * | beta_loop_order                  | 0, 1, 2, 3, 4                                   | 4 (= 4-loop)    |
 * | threshold_corrections_loop_order | 0, 1, 2, 3, 4                                   | 4 (= 4-loop)    |
 * | higgs_2loop_correction_at_as     | 0, 1                                            | 1 (= enabled)   |
 * | higgs_2loop_correction_ab_as     | 0, 1                                            | 1 (= enabled)   |
 * | higgs_2loop_correction_at_at     | 0, 1                                            | 1 (= enabled)   |
 * | higgs_2loop_correction_atau_atau | 0, 1                                            | 1 (= enabled)   |
 * | force_output                     | 0 (no) or 1 (yes)                               | 0 (= no)        |
 * | top_pole_qcd_corrections         | 0 (1L), 1 (2L), 2 (3L)                          | 1 (= 2L QCD)    |
 * | beta_zero_threshold              | any positive double                             | 1.0e-11         |
 * | calculate_observables            | 0 (no) or 1 (yes)                               | 0 (= no)        |
 * | force_positive_masses            | 0 (no) or 1 (yes)                               | 0 (= no)        |
 * | pole_mass_scale                  | any positive double                             | 0 (= SUSY scale)|
 * | eft_pole_mass_scale              | any positive double                             | 0 (= minimum of {Mt, SUSY scale})|
 * | eft_matching_scale               | any positive double                             | 0 (= SUSY scale)|
 * | eft_matching_loop_order_up       | 0, 1, 2                                         | 2 (= 2-loop)    |
 * | eft_matching_loop_order_down     | 0, 1                                            | 1 (= 1-loop)    |
 * | eft_higgs_index                  | any integer >= 0                                | 0 (= lightest)  |
 * | calculate_bsm_masses             | 0 (no) or 1 (yes)                               | 1 (= yes)       |
 * | threshold_corrections            | positive integer                                | 124111421       |
 * | higgs_3loop_ren_scheme_atb_as2   | 0 (DR'), 1 (MDR'), 2 (H3m)                      | 0 (= DR')       |
 * | higgs_3loop_correction_at_as2    | 0, 1                                            | 1 (= enabled)   |
 * | higgs_3loop_correction_ab_as2    | 0, 1                                            | 1 (= enabled)   |
 * | higgs_3loop_correction_at2_as    | 0, 1                                            | 1 (= enabled)   |
 * | higgs_3loop_correction_at3       | 0, 1                                            | 1 (= enabled)   |
 * | higgs_4loop_correction_at_as3    | 0, 1                                            | 1 (= enabled)   |
 * | loop_library                     | 0(Softsusy),1(Collier),2(Looptools),3(fflite)   | 0 (= Softsusy)  |
 */
void Spectrum_generator_settings::reset()
{
   values[precision]             = 1.0e-4;
   values[max_iterations]        = 0.; // 0 = automatic
   values[solver]                = 0.; // 0 = all
   values[calculate_sm_masses]   = 0.; // 0 = false
   values[pole_mass_loop_order]  = 4.;
   values[ewsb_loop_order]       = 4.;
   values[beta_loop_order]       = 4.;
   values[threshold_corrections_loop_order] = 3.;
   values[higgs_2loop_correction_at_as]     = 1.;
   values[higgs_2loop_correction_ab_as]     = 1.;
   values[higgs_2loop_correction_at_at]     = 1.;
   values[higgs_2loop_correction_atau_atau] = 1.;
   values[force_output]                     = 0;
   values[calculate_sm_masses]   = 0.; // 0 = false
   values[top_pole_qcd_corrections]         = 1.;
   values[beta_zero_threshold]              = 1.0e-11;
   values[calculate_observables]            = 0;
   values[force_positive_masses]            = 0;
   values[pole_mass_scale]                  = 0;
   values[eft_pole_mass_scale]              = 0;
   values[eft_matching_scale]               = 0;
   values[eft_matching_loop_order_up]       = 2.;
   values[eft_matching_loop_order_down]     = 1.;
   values[eft_higgs_index]                  = 0;
   values[calculate_bsm_masses]             = 1.;
   values[threshold_corrections]            = Threshold_corrections().get();
   values[higgs_3loop_ren_scheme_atb_as2]   = 0.;
   values[higgs_3loop_correction_at_as2]    = 1.;
   values[higgs_3loop_correction_ab_as2]    = 1.;
   values[higgs_3loop_correction_at2_as]    = 1.;
   values[higgs_3loop_correction_at3]       = 1.;
   values[higgs_4loop_correction_at_as3]    = 1.;
   values[loop_library]                     = -1.; // -1 = (set via environment FLEXIBLESUSY_LOOP_LIBRARY)
}

Loop_corrections Spectrum_generator_settings::get_loop_corrections() const
{
   Loop_corrections loop_corrections;
   loop_corrections.higgs_at_as     = get(higgs_2loop_correction_at_as);
   loop_corrections.higgs_ab_as     = get(higgs_2loop_correction_ab_as);
   loop_corrections.higgs_at_at     = get(higgs_2loop_correction_at_at);
   loop_corrections.higgs_atau_atau = get(higgs_2loop_correction_atau_atau);
   loop_corrections.higgs_at_as_as  = get(higgs_3loop_correction_at_as2);
   loop_corrections.higgs_ab_as_as  = get(higgs_3loop_correction_ab_as2);
   loop_corrections.higgs_at_at_as  = get(higgs_3loop_correction_at2_as);
   loop_corrections.higgs_at_at_at  = get(higgs_3loop_correction_at3);
   loop_corrections.higgs_at_as_as_as   = get(higgs_4loop_correction_at_as3);
   loop_corrections.higgs_3L_scheme = get(higgs_3loop_ren_scheme_atb_as2);
   loop_corrections.top_qcd         = get(top_pole_qcd_corrections);

   return loop_corrections;
}

void Spectrum_generator_settings::set_loop_corrections(
   const Loop_corrections& loop_corrections)
{
   set(higgs_2loop_correction_at_as, loop_corrections.higgs_at_as);
   set(higgs_2loop_correction_ab_as, loop_corrections.higgs_ab_as);
   set(higgs_2loop_correction_at_at, loop_corrections.higgs_at_at);
   set(higgs_2loop_correction_atau_atau, loop_corrections.higgs_atau_atau);
   set(higgs_3loop_correction_at_as2, loop_corrections.higgs_at_as_as);
   set(higgs_3loop_correction_ab_as2, loop_corrections.higgs_ab_as_as);
   set(higgs_3loop_correction_at2_as, loop_corrections.higgs_at_at_as);
   set(higgs_3loop_correction_at3   , loop_corrections.higgs_at_at_at);
   set(higgs_4loop_correction_at_as3, loop_corrections.higgs_at_as_as_as);
   set(higgs_3loop_ren_scheme_atb_as2, loop_corrections.higgs_3L_scheme);
   set(top_pole_qcd_corrections, loop_corrections.top_qcd);
}

Threshold_corrections Spectrum_generator_settings::get_threshold_corrections() const
{
   return Threshold_corrections(get(threshold_corrections));
}

void Spectrum_generator_settings::set_threshold_corrections(const Threshold_corrections& tc)
{
   set(threshold_corrections, tc.get());
}

std::ostream& operator<<(std::ostream& ostr, const Spectrum_generator_settings& sgs)
{
   ostr << "(";

   for (int i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      ostr << sgs.get(static_cast<Spectrum_generator_settings::Settings>(i));
      if (i < Spectrum_generator_settings::NUMBER_OF_OPTIONS - 1)
         ostr << ", ";
   }

   ostr << ")";

   return ostr;
}

} // namespace flexiblesusy
