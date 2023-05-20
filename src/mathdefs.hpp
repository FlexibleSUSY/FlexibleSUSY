#ifndef mathdefs_hpp
#define mathdefs_hpp

#include "wrappers.hpp"

#include <cmath>
#include <complex>
#include <limits>

namespace flexiblesusy {

using Real = double;
using Comp = std::complex<Real>;


constexpr double pi = Pi;
constexpr double r2 = 1.4142135623730950;
const Real epsilon = std::numeric_limits<Real>::epsilon();


template<unsigned p, bool even_p, class T>
struct _power;

template<unsigned p, class T>
struct _power<p, true, T> {
    static T f(T x) {
	return _power<p/2, (p/2) % 2 == 0, T>::f(x) *
	       _power<p/2, (p/2) % 2 == 0, T>::f(x);
    }
};

template<unsigned p, class T>
struct _power<p, false, T> {
    static T f(T x) {
	return _power<p-1, (p-1) % 2 == 0, T>::f(x) * x;
    }
};

template<class T>
struct _power<1, false, T> {
    static T f(T x) {
	return x;
    }
};

template<class T>
struct _power<0, true, T> {
    static T f(T /* x */) {
	return 1;
    }
};

template<unsigned p, class T> inline T Pow(T x)
{
    return _power<p, p % 2 == 0, T>::f(x);
}

template<class T> inline T Pow2(T x) { return Pow<2>(x); }
template<class T> inline T Pow3(T x) { return Pow<3>(x); }
template<class T> inline T Pow4(T x) { return Pow<4>(x); }
template<class T> inline T Pow5(T x) { return Pow<5>(x); }
template<class T> inline T sqr (T x) { return Pow<2>(x); }
template<class T> inline T cube(T x) { return Pow<3>(x); }

template<class T> inline T sign(T x)
{
    if (x == 0) return 0;
    else if (x > 0) return 1;
    else return -1;
}

} // namespace flexiblesusy

#endif // mathdefs_hpp
