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

#ifndef FLEXIBLEDECAY_SETTINGS_H
#define FLEXIBLEDECAY_SETTINGS_H

#include <Eigen/Core>

namespace flexiblesusy {

class FlexibleDecay_settings {
public:
   /// FlexibleDecay settings
   enum Settings : int {
      calculate_decays,      ///< [0] calculate particle decays
      include_higher_order_corrections, ///< [1] include higher order corrections in decays
      offshell_VV_decays,    ///< [2]
      NUMBER_OF_OPTIONS      ///< number of possible options
   };

   using Settings_t = Eigen::Array<double,NUMBER_OF_OPTIONS,1>;

   FlexibleDecay_settings();

   double get(Settings) const; ///< get value of spectrum generator setting
   Settings_t get() const;     ///< get all spectrum generator settings
   std::string get_description(Settings) const; ///< get description of spectrum generator setting
   void set(Settings, double); ///< set value of spectrum generator setting
   void set(const Settings_t&);///< set all spectrum generator settings
   void reset();               ///< resets all settings to their defaults

private:
   std::array<double, NUMBER_OF_OPTIONS> values; ///< spectrum generator settings
};

std::ostream& operator<<(std::ostream&, const FlexibleDecay_settings&);

} // namespace flexiblesusy

#endif
