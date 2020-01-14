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

#include "threshold_corrections.hpp"
#include "error.hpp"
#include "logger.hpp"
#include <string>
#include <iostream>

namespace flexiblesusy {

namespace {

/// returns digit [0-9] in flags at position pos
int get_digit(Threshold_corrections::Flags_t flags, int pos)
{
   if (pos < 0) {
      throw OutOfBoundsError(
         "get_digit: position ( " + std::to_string(pos) + ") must be positive");
   }

   int num = std::abs(static_cast<int>(flags));

   while (--pos >= 0) {
      num /= 10;
   }

   return static_cast<Threshold_corrections::Flags_t>(num % 10);
}

/// sets digit [0-9] in flags at position pos
void set_digit(Threshold_corrections::Flags_t& flags, int pos, int digit)
{
   if (pos < 0) {
      throw OutOfBoundsError(
         "set_digit: position ( " + std::to_string(pos) + ") must be positive");
   }

   if (digit > 9) {
      throw OutOfBoundsError(
         "set_digit: digit ( " + std::to_string(digit) + ") must be less or equal than 9.");
   }

   if (digit < 0) {
      WARNING("digit at position " << pos << " is negative (" << digit << ")."
              " I'm setting it to zero.");
      digit = 0;
   }

   const auto old_digit = get_digit(flags, pos);

   int dig = digit - old_digit;

   while (--pos >= 0) {
      dig *= 10;
   }

   flags += dig;
}

} // anonymous namespace

Threshold_corrections::Threshold_corrections(Flags_t flags)
{
   set(flags);
}

void Threshold_corrections::set(Flags_t flags)
{
   alpha_em    = get_digit(flags, static_cast<int>(Positions::alpha_em   ));
   sin_theta_w = get_digit(flags, static_cast<int>(Positions::sin_theta_w));
   alpha_s     = get_digit(flags, static_cast<int>(Positions::alpha_s    ));
   mz          = get_digit(flags, static_cast<int>(Positions::mz         ));
   mw          = get_digit(flags, static_cast<int>(Positions::mw         ));
   mh          = get_digit(flags, static_cast<int>(Positions::mh         ));
   mt          = get_digit(flags, static_cast<int>(Positions::mt         ));
   mb          = get_digit(flags, static_cast<int>(Positions::mb         ));
   mtau        = get_digit(flags, static_cast<int>(Positions::mtau       ));
}

Threshold_corrections::Flags_t Threshold_corrections::get() const
{
   Flags_t flags = 0;

   set_digit(flags, static_cast<int>(Positions::alpha_em   ), alpha_em);
   set_digit(flags, static_cast<int>(Positions::sin_theta_w), sin_theta_w);
   set_digit(flags, static_cast<int>(Positions::alpha_s    ), alpha_s);
   set_digit(flags, static_cast<int>(Positions::mz         ), mz);
   set_digit(flags, static_cast<int>(Positions::mw         ), mw);
   set_digit(flags, static_cast<int>(Positions::mh         ), mh);
   set_digit(flags, static_cast<int>(Positions::mt         ), mt);
   set_digit(flags, static_cast<int>(Positions::mb         ), mb);
   set_digit(flags, static_cast<int>(Positions::mtau       ), mtau);

   return flags;
}

std::ostream& operator<<(std::ostream& ostr, const Threshold_corrections& tc)
{
   ostr << '['
        << tc.mtau << ','
        << tc.mb << ','
        << tc.mt << ','
        << tc.mh << ','
        << tc.mw << ','
        << tc.mz << ','
        << tc.alpha_s << ','
        << tc.sin_theta_w << ','
        << tc.alpha_em << ']';

   return ostr;
}

} // namespace flexiblesusy
