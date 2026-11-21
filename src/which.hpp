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

#ifndef WHICH_H
#define WHICH_H

#include "pp_map.hpp"
#include <utility>
#include <type_traits>

namespace flexiblesusy {

#define LAMBDA(expr) [&](){ return expr; }

#define WHICH(...) lazy_which(MAP_LIST(LAMBDA, __VA_ARGS__))

template<typename If, typename Then>
auto lazy_which(If&& cif, Then&& cthen) -> decltype(cthen())
{
    return cif() ? cthen() : decltype(cthen()){};
}

template <typename If, typename Then, typename... Elses>
auto lazy_which(If&& cif, Then&& cthen, Elses&&... celses)
   -> std::common_type_t<decltype(cthen()), decltype(celses())...>
{
    return cif() ? cthen() : lazy_which(std::forward<Elses>(celses)...);
}

} // namespace flexiblesusy

#endif // sum_hpp
