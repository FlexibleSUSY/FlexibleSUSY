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
#include "loop_library.hpp"
#include "loop_library_interface.hpp"

#include "softsusy.hpp"

#ifdef ENABLE_COLLIER
#include "collier.hpp"
#define COLLIER_INFO ", 1 (=Collier)"
#else
#define COLLIER_INFO
#endif // ENABLE_COLLIER

#ifdef ENABLE_LOOPTOOLS
#include "looptools.hpp"
#define LOOPTOOLS_INFO ", 2 (=LoopTools)"
#else
#define LOOPTOOLS_INFO
#endif // ENABLE_LOOPTOOLS

#define STRINGIFY(X) #X
#define TOSTR(MACROS) STRINGIFY(MACROS)

namespace flexiblesusy {

int Loop::type_ = -1;
std::unique_ptr<looplibrary::Loop_library_interface> Loop::lib_;

void Loop::setLibrary(int new_type) {
   if( Loop::type_ == -1) {
      switch(new_type) {
         case 0 : Loop::lib_ = std::make_unique<looplibrary::Softsusy>();
                  Loop::type_ = 0;
                  break;
#ifdef ENABLE_COLLIER
         case 1 : Loop::lib_ = std::make_unique<looplibrary::Collier>();
                  Loop::type_ = 1;
                  break;
#endif // ENABLE_COLLIER
#ifdef ENABLE_LOOPTOOLS
         case 2 : Loop::lib_ = std::make_unique<looplibrary::Looptools>();
                  Loop::type_ = 2;
                  break;
#endif // ENABLE_LOOPTOOLS
         default: throw std::invalid_argument("Currently configured values are 0 (=Softsusy)" COLLIER_INFO LOOPTOOLS_INFO ".");
                  break;
      }
   }
}

looplibrary::Loop_library_interface& Loop::library() {
   if(Loop::type_ == -1) {
      Loop::lib_ = std::make_unique<looplibrary::Softsusy>();
      Loop::type_ = 0;
   }
   return *Loop::lib_;
}

} // namespace flexiblesusy
