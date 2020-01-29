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

#ifndef LOOP_LIBRARY_INTERFACE
#define LOOP_LIBRARY_INTERFACE

#include <complex>

#define two_point_virtual(NAME)\
   virtual std::complex<double> NAME(\
      std::complex<double>,\
      std::complex<double>, std::complex<double>,\
      double) = 0;
#define three_point_virtual(NAME)\
   virtual std::complex<double> NAME(\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      double) = 0;
#define four_point_virtual(NAME)\
   virtual std::complex<double> NAME(\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>,\
      double) = 0;

namespace looplibrary
{

class Loop_library_interface
{
   public:
      virtual std::complex<double> A0(std::complex<double>, double) = 0;

      two_point_virtual(B0)
      two_point_virtual(B1)
      two_point_virtual(B00)

      three_point_virtual(C0)
      three_point_virtual(C1)
      three_point_virtual(C2)
      three_point_virtual(C00)
      three_point_virtual(C11)
      three_point_virtual(C12)
      three_point_virtual(C22)

      four_point_virtual(D0)
      four_point_virtual(D00)
      four_point_virtual(D1)
      four_point_virtual(D11)
      four_point_virtual(D12)
      four_point_virtual(D13)
      four_point_virtual(D2)
      four_point_virtual(D22)
      four_point_virtual(D23)
      four_point_virtual(D3)
      four_point_virtual(D33)

      virtual void A(
         std::complex<double> (&)[1],
         std::complex<double>,
         double) = 0;
      virtual void B(
         std::complex<double> (&)[2],
         std::complex<double>, std::complex<double>, std::complex<double>,
         double) = 0;
      virtual void C(
         std::complex<double> (&)[7],
         std::complex<double>, std::complex<double>, std::complex<double>,
         std::complex<double>, std::complex<double>, std::complex<double>,
         double) = 0;
      virtual void D(
         std::complex<double> (&)[11],
         std::complex<double>, std::complex<double>, std::complex<double>,
         std::complex<double>, std::complex<double>, std::complex<double>,
         std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>,
         double) = 0;
};

} // namespace looplibrary

#endif
