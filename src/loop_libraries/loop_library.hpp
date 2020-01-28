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
#ifdef ENABLE_LOOP_LIBRARY

#ifndef LOOP
#define LOOP

#include <memory>
#include "loop_library_interface.hpp"

namespace flexiblesusy {

class Loop {
public:
   static void setLibrary(int);
   static looplibrary::Loop_library_interface& library();

private:
   static int type_;
   static std::unique_ptr<looplibrary::Loop_library_interface> lib_;

   Loop() {}
   Loop(Loop const&);
   void operator=(Loop const&);
};

} // namespace flexiblesusy

#endif // LOOP
#endif // ENABLE_LOOP_LIBRARY
