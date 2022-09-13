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

#ifndef SUM_H
#define SUM_H

#include <type_traits>
#include <cstddef>
#include <Eigen/Core>

namespace flexiblesusy {

#define SUM(...) (get_sum(__VA_ARGS__)(__VA_ARGS__))

#define get_sum(...) get_sum_macro(__VA_ARGS__, sum_user_t, sum_ptrdiff_t,)

#define get_sum_macro(_1, _2, _3, _4, _5, name, ...) name

#define sum_ptrdiff_t(idx, ini, fin, expr)      \
    sum<std::ptrdiff_t>((ini), (fin), [&](std::ptrdiff_t (idx)) { return (expr); })

#define sum_user_t(type, idx, ini, fin, expr)	\
    sum<type>((ini), (fin), [&](type (idx)) { return (expr); })

template<typename T>
struct is_eigen_type
{
    static constexpr auto value =
	std::is_base_of<Eigen::EigenBase<T>, T>::value;
};

template<typename Idx, typename Function, bool isEigenType>
struct EvalEigenXprImpl {
    static auto eval(Idx i, Function f) -> decltype(f(i)) {
	return f(i);
    }
};

template<typename Idx, typename Function>
struct EvalEigenXprImpl<Idx, Function, true> {
    static auto eval(Idx i, Function f) ->
	typename std::remove_reference<decltype(f(i).eval())>::type
    {
	return f(i).eval();
    }
};

template<typename Idx, typename Function>
auto EvalEigenXpr(Idx i, Function f) ->
    decltype(
	EvalEigenXprImpl<Idx, Function, is_eigen_type<decltype(f(i))>::value>::
	eval(i, f))
{
    return
	EvalEigenXprImpl<Idx, Function, is_eigen_type<decltype(f(i))>::value>::
	eval(i, f);
}

template<typename T, bool isEigenType>
struct create_zero {
    static const T zero() {
	return T();
    }
};

template<typename T>
struct create_zero<T, true> {
    static const T zero() {
	T z;
	z.setZero();
	return z;
    }
};

template<class Idx, class Function>
auto sum(Idx ini, Idx fin, Function f) -> decltype(EvalEigenXpr<Idx>(ini, f))
{
    using Evaled = decltype(EvalEigenXpr<Idx>(ini, f));
    using Acc = typename std::remove_cv<Evaled>::type;
    Acc s = create_zero<Acc, is_eigen_type<Acc>::value>::zero();
    for (Idx i = ini; i <= fin; i++) s += f(i);
    return s;
}

} // namespace flexiblesusy

#endif // sum_hpp
