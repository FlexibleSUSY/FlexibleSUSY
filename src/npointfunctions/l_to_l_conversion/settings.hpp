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

#ifndef LTOLCONVERSION_SETTINGS_H
#define LTOLCONVERSION_SETTINGS_H

#include <array>
#include <Eigen/Core>

namespace flexiblesusy {

class LToLConversion_settings {
public:
   /// LToLConversion settings
   enum Settings : int {
      include_tensor_contribution,   ///< [0]
      include_gluonic_contribution,   ///< [1]
      NUMBER_OF_OPTIONS   ///< number of possible options
   };

   using Settings_t = Eigen::Array<double,NUMBER_OF_OPTIONS,1>;

   LToLConversion_settings();

   double get(Settings) const; ///< get value of setting
   Settings_t get() const;     ///< get all settings
   std::string get_description(Settings) const; ///< get description setting
   void set(Settings, double); ///< set value of setting
   void set(const Settings_t&);///< set all settings
   void reset();               ///< resets all settings to their defaults
   double tensor() const { return get(include_tensor_contribution); }

private:
   std::array<double, NUMBER_OF_OPTIONS> values; ///< settings
};

std::ostream& operator<<(std::ostream&, const LToLConversion_settings&);

} // namespace flexiblesusy

#endif
