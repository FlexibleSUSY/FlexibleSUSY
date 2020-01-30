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

#include "library_softsusy.hpp"

#ifdef ENABLE_COLLIER
#include "library_collier.hpp"
#define COLLIER_INFO ", 1 (=COLLIER)"
#else
#define COLLIER_INFO
#endif // ENABLE_COLLIER

#ifdef ENABLE_LOOPTOOLS
#include "library_looptools.hpp"
#define LOOPTOOLS_INFO ", 2 (=LoopTools)"
#else
#define LOOPTOOLS_INFO
#endif // ENABLE_LOOPTOOLS

#ifdef ENABLE_FFLITE
#include "library_fflite.hpp"
#define FFLITE_INFO ", 3 (=fflite)"
#else
#define FFLITE_INFO
#endif // ENABLE_FFLITE

#define STRINGIFY(X) #X
#define TOSTR(MACROS) STRINGIFY(MACROS)

namespace flexiblesusy {

int Looplibrary::type_ = -1;
std::unique_ptr<looplibrary::Loop_library_interface> Looplibrary::lib_;

void Looplibrary::set(int new_type) {
   if( Looplibrary::type_ == -1) {
      switch(new_type) {
         case 0 : Looplibrary::lib_ = std::make_unique<looplibrary::Softsusy>();
                  Looplibrary::type_ = 0;
                  break;
#ifdef ENABLE_COLLIER
         case 1 : Looplibrary::lib_ = std::make_unique<looplibrary::Collier>();
                  Looplibrary::type_ = 1;
                  break;
#endif // ENABLE_COLLIER
#ifdef ENABLE_LOOPTOOLS
         case 2 : Looplibrary::lib_ = std::make_unique<looplibrary::Looptools>();
                  Looplibrary::type_ = 2;
                  break;
#endif // ENABLE_LOOPTOOLS
#ifdef ENABLE_FFLITE
         case 3 : Looplibrary::lib_ = std::make_unique<looplibrary::Fflite>();
                  Looplibrary::type_ = 3;
                  break;
#endif // ENABLE_FFLITE
         default: throw std::invalid_argument("Currently configured values are 0 (=Softsusy)" COLLIER_INFO LOOPTOOLS_INFO FFLITE_INFO".");
                  break;
      }
   }
}

looplibrary::Loop_library_interface& Looplibrary::get() {
   if(Looplibrary::type_ == -1) throw std::logic_error("Loop library should be initialized before first usage.");
   return *Looplibrary::lib_;
}

} // namespace flexiblesusy
