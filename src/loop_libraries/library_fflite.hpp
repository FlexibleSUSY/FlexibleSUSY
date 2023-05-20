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

#ifndef LOOP_LIBRARY_FFLITE_H
#define LOOP_LIBRARY_FFLITE_H

#include "loop_library_interface.hpp"

#define REDEFINE(R, ARGS, NAME)                                                \
   std::complex<double> NAME ARGS noexcept override;

namespace flexiblesusy
{
namespace looplibrary
{
class Fflite : public Loop_library_interface
{
public:
   Fflite();
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (A_ARGS), A_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (B_ARGS), B_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (C_ARGS), C_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (D_ARGS), D_SEQ)
   void A(Acoeff_t&, A_ARGS) noexcept override;
   void B(Bcoeff_t&, B_ARGS) noexcept override;
   void C(Ccoeff_t&, C_ARGS) noexcept override;
   void D(Dcoeff_t&, D_ARGS) noexcept override;
   ~Fflite() noexcept override {};
};
} // namespace looplibrary
} // namespace flexiblesusy

#endif // LOOP_LIBRARY_FFLITE_H
