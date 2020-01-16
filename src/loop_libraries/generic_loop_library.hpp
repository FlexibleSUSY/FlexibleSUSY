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

#include "config.h"
#ifdef ENABLE_GENERIC_LOOP_LIBRARY

#ifndef GENERIC_LOOP
#define GENERIC_LOOP

#include "loop_library_interface.hpp"

namespace flexiblesusy {

class Generic_loop {
public:
   static void setLibrary(int);
   static Loop_library_interface& library();

private:
   static int type_;
   static std::unique_ptr<Loop_library_interface> lib_;

   Generic_loop() {}
   Generic_loop(Generic_loop const&);
   void operator=(Generic_loop const&);
};

} // namespace flexiblesusy

#endif // GENERIC_LOOP
#endif // ENABLE_GENERIC_LOOP_LIBRARY
