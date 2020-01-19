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

#include "loop_library_interface.hpp"

#define two_point(NAME)\
   std::complex<double> NAME(\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      double) noexcept;
#define three_point(NAME)\
   std::complex<double> NAME(\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      double) noexcept;
#define four_point(NAME)\
   std::complex<double> NAME(\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      std::complex<double>, std::complex<double>, std::complex<double>,\
      std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>,\
      double) noexcept;

namespace looplibrary {

class Collier : public Loop_library_interface
{
   private:
      double current_mu2_uv;
      void initialize() noexcept;
      void set_mu2_uv(double) noexcept;

   public:
      Collier() : current_mu2_uv(1.0) {
         initialize();
      }

      std::complex<double> A0(std::complex<double>, double) noexcept;

      two_point(B0)
      two_point(B1)

      three_point(C0)
      three_point(C1)
      three_point(C2)
      three_point(C00)
      three_point(C11)
      three_point(C12)
      three_point(C22)

      four_point(D0)
      four_point(D00)
      four_point(D1)
      four_point(D11)
      four_point(D12)
      four_point(D13)
      four_point(D2)
      four_point(D22)
      four_point(D23)
      four_point(D3)
      four_point(D33)

      void get_A(
         std::complex<double> (&)[1],
         std::complex<double>,
         double) noexcept;
      void get_B(
         std::complex<double> (&)[2],
         std::complex<double>, std::complex<double>, std::complex<double>,
         double) noexcept;
      void get_C(
         std::complex<double> (&)[7],
         std::complex<double>, std::complex<double>, std::complex<double>,
         std::complex<double>, std::complex<double>, std::complex<double>,
         double) noexcept;
      void get_D(
         std::complex<double> (&)[11],
         std::complex<double>, std::complex<double>, std::complex<double>,
         std::complex<double>, std::complex<double>, std::complex<double>,
         std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>,
         double) noexcept;
};

} // namespace looplibrary
