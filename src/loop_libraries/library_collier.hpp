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

#ifndef LOOP_LIBRARY_COLLIER_H
#define LOOP_LIBRARY_COLLIER_H

#include "loop_library_interface.hpp"

#include <array>
#include <cstddef>
#include <unordered_map>

#define REDEFINE(R, ARGS, NAME)                                                \
   std::complex<double> NAME ARGS noexcept override;

namespace flexiblesusy
{
namespace looplibrary
{

/**
 * Remembers tensor coefficients by the arguments passed to COLLIER
 * (momenta, masses and the renormalization scale).  The generated 1-loop
 * decay amplitudes ask for the same B and C tensors many times (MSSM:
 * 80% of the B and 30% of the C calls repeat).  COLLIER's own global
 * cache doesn't fit: it identifies integrals by their position in a fixed
 * sequence of calls.  The cache is cleared when it is full.  Like COLLIER
 * itself, it isn't thread-safe.
 */
template <std::size_t N, typename Coeff>
class Tensor_cache
{
public:
   using Key = std::array<double, N>;

   /// returns the cached coefficients, or nullptr
   const Coeff* find(const Key&) const noexcept;
   /// stores the coefficients; does nothing if memory can't be allocated
   void insert(const Key&, const Coeff&) noexcept;

private:
   static constexpr std::size_t max_size = 10000;

   struct Hash {
      std::size_t operator()(const Key&) const noexcept;
   };
   std::unordered_map<Key, Coeff, Hash> map{};
};

class Collier : public Loop_library_interface
{
private:
   double current_mu2_uv;
   Tensor_cache<6, Bcoeff_t> b_cache{};
   Tensor_cache<10, Ccoeff_t> c_cache{};
   Tensor_cache<15, Dcoeff_t> d_cache{};
   void initialize() noexcept;
   void set_mu2_uv(double) noexcept;

public:
   Collier() : current_mu2_uv(1.0) { initialize(); }
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (A_ARGS), A_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (B_ARGS), B_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (DB_ARGS), DB_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (C_ARGS), C_SEQ)
   BOOST_PP_SEQ_FOR_EACH(REDEFINE, (D_ARGS), D_SEQ)
   void A(Acoeff_t&, A_ARGS) noexcept override;
   void B(Bcoeff_t&, B_ARGS) noexcept override;
   void DB(DBcoeff_t&, DB_ARGS) noexcept override;
   void C(Ccoeff_t&, C_ARGS) noexcept override;
   void D(Dcoeff_t&, D_ARGS) noexcept override;
   ~Collier() noexcept override{};
};
} // namespace looplibrary
} // namespace flexiblesusy

#endif // LOOP_LIBRARY_COLLIER_H
