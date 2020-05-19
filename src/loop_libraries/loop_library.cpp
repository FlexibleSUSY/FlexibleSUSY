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

#include <cstdlib>
#include <memory>

#include "config.h"
#include "logger.hpp"
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

namespace flexiblesusy
{

Loop_library::Library Loop_library::type_ = Loop_library::Library::Undefined;
std::unique_ptr<looplibrary::Loop_library_interface> Loop_library::lib_;

void Loop_library::set_default()
{
   Loop_library::lib_ = std::make_unique<looplibrary::Softsusy>();
   Loop_library::type_ = Loop_library::Library::Softsusy;
   VERBOSE_MSG("Enabling softsusy.");
}

void Loop_library::set(int new_type)
{
   if (Loop_library::type_ == Loop_library::Library::Undefined) {
      if (new_type == -1) {
         if (const char* flag = std::getenv("FLEXIBLESUSY_LOOP_LIBRARY")) {
            new_type = std::atoi(flag);
            VERBOSE_MSG("Setting looplibrary using environment variable "
                        << "FLEXIBLESUSY_LOOP_LIBRARY=" << new_type);
         }
      }
      switch (new_type) {
      case 0:
         Loop_library::set_default();
         break;
#ifdef ENABLE_COLLIER
      case 1:
         Loop_library::lib_ = std::make_unique<looplibrary::Collier>();
         Loop_library::type_ = Loop_library::Library::Collier;
         VERBOSE_MSG("Enabling COLLIER.");
         break;
#endif // ENABLE_COLLIER
#ifdef ENABLE_LOOPTOOLS
      case 2:
         Loop_library::lib_ = std::make_unique<looplibrary::Looptools>();
         Loop_library::type_ = Loop_library::Library::Looptools;
         VERBOSE_MSG("Enabling LoopTools.");
         break;
#endif // ENABLE_LOOPTOOLS
#ifdef ENABLE_FFLITE
      case 3:
         Loop_library::lib_ = std::make_unique<looplibrary::Fflite>();
         Loop_library::type_ = Loop_library::Library::Fflite;
         VERBOSE_MSG("Enabling FFlite.");
         break;
#endif // ENABLE_FFLITE
      default:
         ERROR("Warning: Check FlexibleSUSY[31]:\n"
               "Currently configured values are 0 (=softsusy)" COLLIER_INFO
                  LOOPTOOLS_INFO FFLITE_INFO ".\n"
               "Setting default library.");
         Loop_library::set_default();
         break;
      }
   }
}

Loop_library::Library Loop_library::get_type()
{
   return type_;
}

looplibrary::Loop_library_interface& Loop_library::get()
{
   if (Loop_library::type_ == Loop_library::Library::Undefined) {
      ERROR("Loop library should be initialized before first usage.\n"
            "Setting default library.");
      Loop_library::set(-1);
   }
   return *Loop_library::lib_;
}

} // namespace flexiblesusy
