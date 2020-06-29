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

#ifndef FOR_EACH_H
#define FOR_EACH_H

#include <tuple>

namespace flexiblesusy {

namespace detail
{
   template<int... Is>
   struct seq {};

   template<int N, int... Is>
   struct gen_seq : gen_seq<N - 1, N - 1, Is...> {};

   template<int... Is>
   struct gen_seq<0, Is...> : seq<Is...> {};

   template<typename T, typename F, int... Is>
   void for_each(T&& t, F f, seq<Is...>)
   {
      auto l = { (f(std::get<Is>(t)), 0)... };
   }
} // namespace detail

/// applies f on each element of the tuple
template<typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...> const& t, F f)
{
   detail::for_each(t, f, detail::gen_seq<sizeof...(Ts)>());
}

} // namespace flexiblesusy

#endif
