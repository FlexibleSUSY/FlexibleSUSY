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

#include <memory>

#include "config.h"
#include "generic_loop_library.hpp"
#include "loop_library_interface.hpp"
#include "collier.hpp"
#include "loop_tools.hpp"

namespace flexiblesusy {

int Generic_loop::type_ = -1;
std::unique_ptr<Loop_library_interface> Generic_loop::lib_;

void Generic_loop::setLibrary(int new_type) {
   if( Generic_loop::type_ == -1) {
      switch(new_type) {
         case 0 : Generic_loop::lib_ = std::make_unique<Collier>(); // @ToDo add SoftSUSY
                  Generic_loop::type_ = 0;
                  break;
#ifdef ENABLE_COLLIER
         case 1 : Generic_loop::lib_ = std::make_unique<Collier>();
                  Generic_loop::type_ = 1;
                  break;
#endif // ENABLE_COLLIER
#ifdef ENABLE_LOOPTOOLS
         case 2 : Generic_loop::lib_ = std::make_unique<Looptools>();
                  Generic_loop::type_ = 2;
                  break;
#endif // ENABLE_LOOPTOOLS
         default: throw std::invalid_argument("Unrecognized Generic Loop Library (check table inside FlexibleSUSY/src/spectrum_generator_settings.cpp)");
                  break;
      }
      Generic_loop::type_ = 1;
   }
}

Loop_library_interface& Generic_loop::library() {
   if( Generic_loop::type_ == -1) {
      Generic_loop::lib_ = std::make_unique<Collier>(); // @ToDo add SoftSUSY
      Generic_loop::type_ = 0;
   }
   return *Generic_loop::lib_;
}

} // namespace flexiblesusy
