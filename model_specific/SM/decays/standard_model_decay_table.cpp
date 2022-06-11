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


#include "standard_model_decay_table.hpp"

#include "error.hpp"

#include <iostream>
#include <sstream>

namespace flexiblesusy {

SM_decay_table::SM_decay_table()
   : decay_table({Decays_list(25)})
{
}

void SM_decay_table::clear()
{
   for (auto& particle: decay_table) {
      particle.clear();
   }
}

/*
void SM_decay_table::print(std::ostream& ostr) const
{
   for (const auto& particle: decay_table) {
      const auto pdg = particle.get_particle_id();
      for (const auto& decay: particle) {
         const auto& final_state = decay.second.get_final_state_particle_ids();
         const std::size_t final_state_size = final_state.size();
         ostr << "Partial width["
              << SM_info::get_particle_name_from_pdg(pdg)
              << " -> ";
         for (std::size_t i = 0; i < final_state_size; ++i) {
            ostr << SM_info::get_particle_name_from_pdg(final_state[i]);
            if (i < final_state_size - 1) {
               ostr << ", ";
            }
         }
         ostr << "] = " << decay.second.get_width() << '\n';
      }
   }
}

std::ostream& operator<<(std::ostream& ostr,
                         const SM_decay_table& table)
{
   table.print(ostr);
   return ostr;
}
*/

Decays_list& SM_decay_table::get_hh_decays()
{
   return decay_table[0];
}

const Decays_list& SM_decay_table::get_hh_decays() const
{
   return decay_table[0];
}

} // namespace flexiblesusy
