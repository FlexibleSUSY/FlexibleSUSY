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

#ifndef RAII_H
#define RAII_H

#include <utility>

namespace flexiblesusy {

/**
 * @class RAII_save
 * @brief Saves value of variable and restores it at destruction
 */
template <typename T>
class RAII_save {
public:
   RAII_save(T& var_) noexcept : var(var_), value(var_) {}
   RAII_save(const RAII_save&) = delete;
   RAII_save(RAII_save&&) noexcept = default;
   ~RAII_save() { var = value; }
   RAII_save& operator=(const RAII_save&) = delete;
   RAII_save& operator=(RAII_save&& other) noexcept = default;

private:
   T& var;
   T value{};
};

template <typename T>
constexpr RAII_save<T> make_raii_save(T& var)
{
   return RAII_save<T>(var);
}

/**
 * @class RAII_guard
 * @brief Carries out provided clean-up actions at destruction
 */
template <typename F>
class RAII_guard {
public:
   RAII_guard(F f_) : clean_up(std::move(f_)) {}
   RAII_guard(const RAII_guard&) = delete;
   RAII_guard(RAII_guard&&) noexcept = default;
   ~RAII_guard() { clean_up(); }
   RAII_guard& operator=(const RAII_guard&) = delete;
   RAII_guard& operator=(RAII_guard&&) noexcept = default;
private:
   F clean_up;
};

template <typename F>
constexpr RAII_guard<F> make_raii_guard(F f)
{
   return RAII_guard<F>(std::move(f));
}

} // namespace flexiblesusy

#endif
