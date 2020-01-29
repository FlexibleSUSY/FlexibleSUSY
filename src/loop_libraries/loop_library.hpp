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

#ifndef LOOPLIBRARY_HPP
#define LOOPLIBRARY_HPP

#include <memory>
#include "loop_library_interface.hpp"

namespace flexiblesusy {

class Looplibrary {
public:
   static void set(int);
   static looplibrary::Loop_library_interface& get();

private:
   static int type_;
   static std::unique_ptr<looplibrary::Loop_library_interface> lib_;

   Looplibrary() {}
   Looplibrary(Looplibrary const&);
   void operator=(Looplibrary const&);
};

} // namespace flexiblesusy

#endif // LOOPLIBRARY_HPP
