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

#include <limits>
#include <cmath>
#include <cstdlib>
#include "derivative.hpp"
#include "logger.hpp"
#include "pv.hpp"

#ifdef ENABLE_LOOPTOOLS
#  include <clooptools.h>
#elif defined(ENABLE_FFLITE)
#  include "fflite.hpp"
#else
#  include "numerics.h"
#endif

namespace flexiblesusy {

namespace passarino_veltman {

using namespace std;

#ifdef ENABLE_LOOPTOOLS

namespace {

struct Initialize_looptools {
    Initialize_looptools() {
	ltini();
    }
    ~Initialize_looptools() {
	ltexi();
    }
} initialize_looptools;

const double deriv_eps = 1e-5; ///< epsilon for derivatives

} // anonymous namespace

complex<double> A0(double m2, double scl2) noexcept
{
    setmudim(scl2);
    return ::A0(m2);
}

complex<double> B0(double p2, double m2a, double m2b, double scl2) noexcept
{
    setmudim(scl2);
    return ::B0(p2, m2a, m2b);
}

complex<double> B1(double p2, double m2a, double m2b, double scl2) noexcept
{
    setmudim(scl2);
    return ::B1(p2, m2a, m2b);
}

complex<double> B00(double p2, double m2a, double m2b, double scl2) noexcept
{
    setmudim(scl2);
    return ::B00(p2, m2a, m2b);
}

std::complex<double> C0(double p2a, double p2b, double p2ab,
                        double m2a, double m2b, double m2c,
                        double scl2) noexcept
{
   setmudim(scl2);
   return ::C0(p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C1(double p2a, double p2b, double p2ab,
                        double m2a, double m2b, double m2c,
                        double scl2) noexcept
{
   setmudim(scl2);
   return ::C0i(cc1, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C2(double p2a, double p2b, double p2ab,
                        double m2a, double m2b, double m2c,
                        double scl2) noexcept
{
   setmudim(scl2);
   return ::C0i(cc2, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C00(double p2a, double p2b, double p2ab,
                         double m2a, double m2b, double m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0i(cc00, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C11(double p2a, double p2b, double p2ab,
                         double m2a, double m2b, double m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0i(cc11, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C12(double p2a, double p2b, double p2ab,
                         double m2a, double m2b, double m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0i(cc12, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C22(double p2a, double p2b, double p2ab,
                         double m2a, double m2b, double m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0i(cc22, p2a, p2b, p2ab, m2a, m2b, m2c);
}

complex<double> A0(complex<double> m2, double scl2) noexcept
{
    setmudim(scl2);
    return ::A0C(m2);
}

complex<double> B0
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    setmudim(scl2);
    return ::B0C(p2, m2a, m2b);
}

complex<double> B1
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    setmudim(scl2);
    return ::B1C(p2, m2a, m2b);
}

complex<double> B00
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    setmudim(scl2);
    return ::B00C(p2, m2a, m2b);
}

std::complex<double> C0(std::complex<double> p2a, std::complex<double> p2b,
                        std::complex<double> p2ab, std::complex<double> m2a,
                        std::complex<double> m2b, std::complex<double> m2c,
                        double scl2) noexcept
{
   setmudim(scl2);
   return ::C0C(p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C1(std::complex<double> p2a, std::complex<double> p2b,
                        std::complex<double> p2ab, std::complex<double> m2a,
                        std::complex<double> m2b, std::complex<double> m2c,
                        double scl2) noexcept
{
   setmudim(scl2);
   return ::C0iC(cc1, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C2(std::complex<double> p2a, std::complex<double> p2b,
                        std::complex<double> p2ab, std::complex<double> m2a,
                        std::complex<double> m2b, std::complex<double> m2c,
                        double scl2) noexcept
{
   setmudim(scl2);
   return ::C0iC(cc2, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C00(std::complex<double> p2a, std::complex<double> p2b,
                         std::complex<double> p2ab, std::complex<double> m2a,
                         std::complex<double> m2b, std::complex<double> m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0iC(cc00, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C11(std::complex<double> p2a, std::complex<double> p2b,
                         std::complex<double> p2ab, std::complex<double> m2a,
                         std::complex<double> m2b, std::complex<double> m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0iC(cc11, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C12(std::complex<double> p2a, std::complex<double> p2b,
                         std::complex<double> p2ab, std::complex<double> m2a,
                         std::complex<double> m2b, std::complex<double> m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0iC(cc12, p2a, p2b, p2ab, m2a, m2b, m2c);
}

std::complex<double> C22(std::complex<double> p2a, std::complex<double> p2b,
                         std::complex<double> p2ab, std::complex<double> m2a,
                         std::complex<double> m2b, std::complex<double> m2c,
                         double scl2) noexcept
{
   setmudim(scl2);
   return ::C0iC(cc22, p2a, p2b, p2ab, m2a, m2b, m2c);
}

complex<double> D1B0(double p2, double m2a, double m2b) noexcept
{
    const auto f = [m2a,m2b](double p2) { return B0(p2, m2a, m2b, 1.0); };
    return derivative_central<0>(f, p2, deriv_eps);
}

complex<double> D1B0
(complex<double> p2, complex<double> m2a, complex<double> m2b) noexcept
{
    const auto re_f = [p2,m2a,m2b](double re_p2) { return B0({re_p2, p2.imag()}, m2a, m2b, 1.0).real(); };
    const auto im_f = [p2,m2a,m2b](double re_p2) { return B0({re_p2, p2.imag()}, m2a, m2b, 1.0).imag(); };

    const auto dudx = derivative_central<0>(re_f, p2.real(), deriv_eps);
    const auto dvdx = derivative_central<0>(im_f, p2.real(), deriv_eps);

    return {dudx, dvdx};
}

#elif defined(ENABLE_FFLITE)

namespace {

const double nan = numeric_limits<double>::quiet_NaN();
const double deriv_eps = 1e-5; ///< epsilon for derivatives

// see src/include/ff.h in LoopTools
const double acc = 1e-13;
const double eps = 1e-22;
const complex<double> cIeps(0.0, 1e-50);

struct Initialize_looptools {
    Initialize_looptools() {
	ltini_();
    }
    ~Initialize_looptools() {
	ltexi_();
    }
} initialize_looptools;

template<class T> T sqr(T x) noexcept { return x*x; }
template<class T> T sign(T a, T b) noexcept { return b >= 0 ? std::abs(a) : - std::abs(a); }

namespace FF {

complex<double> A0(double m2, double scl2) noexcept
{
    complex<double> ca0;
    int ier;
    ljffxa0_(ca0, 0, scl2, m2, ier);
    return ca0;
}

complex<double> B0(double p2, double m2a, double m2b, double scl2) noexcept
{
    complex<double> cb0;
    int ier;
    ljffxb0_(cb0, 0, scl2, p2, m2a, m2b, ier);
    return cb0;
}

complex<double> B1(double p2, double m2a, double m2b, double scl2) noexcept
{
    complex<double> cb1;
    int ier;
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    double piDpj[9];
    ljffdot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffxb1_(cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb1;
}

complex<double> B00(double p2, double m2a, double m2b, double scl2) noexcept
{
    complex<double> cb2i[2];
    int ier;
    complex<double> cb1 = B1(p2, m2a, m2b, scl2);
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    double piDpj[9];
    ljffdot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffxb2p_(cb2i, cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb2i[1];
}

complex<double> A0
(complex<double> m2, double scl2) noexcept
{
    complex<double> ca0;
    int ier;
    ljffca0_(ca0, 0, scl2, m2, ier);
    return ca0;
}

complex<double> B0
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    complex<double> cb0;
    int ier;
    ljffcb0_(cb0, 0, scl2, p2, m2a, m2b, ier);
    return cb0;
}

complex<double> B1
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    complex<double> cb1;
    int ier;
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    complex<double> piDpj[9];
    ljffcot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffcb1_(cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb1;
}

complex<double> B00
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    complex<double> cb2i[2];
    int ier;
    complex<double> cb1 = B1(p2, m2a, m2b, scl2);
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    complex<double> piDpj[9];
    ljffcot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffcb2p_(cb2i, cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb2i[1];
}

complex<double> D1B0(double p2, double m2a, double m2b) noexcept
{
    const auto f = [m2a,m2b](double p2) { return B0(p2, m2a, m2b, 1.0); };
    return derivative_central<0>(f, p2, deriv_eps);
}

complex<double> D1B0
(complex<double> p2, complex<double> m2a, complex<double> m2b) noexcept
{
    const auto re_f = [p2,m2a,m2b](double re_p2) { return B0({re_p2, p2.imag()}, m2a, m2b, 1.0).real(); };
    const auto im_f = [p2,m2a,m2b](double re_p2) { return B0({re_p2, p2.imag()}, m2a, m2b, 1.0).imag(); };

    const auto dudx = derivative_central<0>(re_f, p2.real(), deriv_eps);
    const auto dvdx = derivative_central<0>(im_f, p2.real(), deriv_eps);

    return {dudx, dvdx};
}

} // namespace FF

namespace AD {

complex<double> fpv
(const int& n, const std::complex<double>& x, const std::complex<double>& y) noexcept
{
    complex<double> res; sub_fpv_(res, n, x, y); return res;
}

complex<double> yfpv
(const int& n, const std::complex<double>& x, const std::complex<double>& y) noexcept
{
    complex<double> res; sub_yfpv_(res, n, x, y); return res;
}

complex<double> fth
(const int& n, const std::complex<double>& x, const std::complex<double>& y) noexcept
{
    complex<double> res; sub_fth_(res, n, x, y); return res;
}

complex<double> xlogx(complex<double> x) noexcept
{
    return std::abs(x) == 0 ? 0.0 : x*log(x);
}

enum BFuncType { b0, b1, b00 };

// the two-point tensor coefficients from Ansgar Denner's bcanew.f
// adapted to the conventions of LoopTools
// translated to C++ for inclusion in FlexibleSUSY

template<BFuncType ft>
complex<double> B(double p2, double m2a, double m2b, double scl2) noexcept
{
    const double mudim = scl2;
    const double delta = 0;

    complex<double> x1, x2, y1, y2, r;
    complex<double> mu, f1, f2, g1, g2, a0;

    complex<double> B_bb0_, B_bb1_, B_bb00_;

    const double m1 = m2a;
    const double m2 = m2b;
    const double p = p2;
    const double dm = m1 - m2;

    // general case
    if ( std::abs(p) > eps*(m1 + m2) ) {
	r = sqrt(complex<double>(p*(p - m1 - m2) -
				 m1*(p - dm) - m2*(p + dm)));
	x1 = .5*(p + dm + r)/p;
	x2 = .5*(p + dm - r)/p;
	if ( std::abs(x2) > std::abs(x1) )
	    x1 = m1/(p*x2);
	else if ( std::abs(x1) > std::abs(x2) )
	    x2 = m1/(p*x1);
	x1 = x1 + sign(std::abs(x1), p)*cIeps;
	x2 = x2 - sign(std::abs(x2), p)*cIeps;

	y2 = .5*(p - dm + r)/p;
	y1 = .5*(p - dm - r)/p;
	if ( std::abs(y2) > std::abs(y1) )
	    y1 = m2/(p*y2);
	else if ( std::abs(y1) > std::abs(y2) )
	    y2 = m2/(p*y1);
	y1 = y1 - sign(std::abs(y1), p)*cIeps;
	y2 = y2 + sign(std::abs(y2), p)*cIeps;

	if ( std::abs(y1) > .5 && std::abs(y2) > .5 ) {
	    mu = log(m2/mudim) - delta;
	    B_bb0_ = -(mu + fpv(1, x1, y1) + fpv(1, x2, y2));
	    if (ft == b0) return B_bb0_;
	    B_bb1_ = 1/2.0*(mu + fpv(2, x1, y1) + fpv(2, x2, y2));
	    if (ft == b1) return B_bb1_;
	}
	else if ( std::abs(x1) < 10 && std::abs(x2) < 10 ) {
	    mu = log(p/mudim*(1.0 - cIeps)) - delta;
	    g1 = xlogx(y1);
	    f1 = xlogx(-x1) - g1 + 1.0;
	    g2 = xlogx(y2);
	    f2 = xlogx(-x2) - g2 + 1.0;
	    B_bb0_ = -(mu - f1 - f2);
	    if (ft == b0) return B_bb0_;
	    f1 = x1*f1 - g1 + 1/2.0;
	    f2 = x2*f2 - g2 + 1/2.0;
	    B_bb1_ = 1/2.0*(mu - f1 - f2);
	    if (ft == b1) return B_bb1_;
	}
	else if ( std::abs(x1) > .5 && std::abs(x2) > .5 ) {
	    mu = log(m1/mudim) - delta +
		fth(1, x1, y1) + fth(1, x2, y2);
	    B_bb0_ = -mu;
	    if (ft == b0) return B_bb0_;
	    mu = mu + fth(2, x1, y1) + fth(2, x2, y2);
	    B_bb1_ = 1/2.0*mu;
	    if (ft == b1) return B_bb1_;
	}
	else {
	    ERROR("Bcoeffb not defined for"
		  "  p = "  << p <<
		  "  m1 = " << m1 <<
		  "  m2 = " << m2);
	    B_bb0_ = nan;
	    if (ft == b0) return B_bb0_;
	    B_bb1_ = nan;
	    if (ft == b1) return B_bb1_;
	}

	a0 = 0;
	if ( m2 != 0 ) a0 = m2*(1 - log(m2/mudim) + delta);

	B_bb00_ = ((p + dm)*B_bb1_ +
		   2*m1*B_bb0_ + a0 + m1 + m2 - p/3.0)/6.0;
	if (ft == b00) return B_bb00_;
    }
    // zero momentum
    else if ( std::abs(dm) > acc*(m1 + m2) ) {
	x2 = m1/dm*(1.0 - cIeps);
	y2 = -m2/dm*(1.0 - cIeps);
	if ( std::abs(y2) > .5 ) {
	    mu = log(m2/mudim) - delta;
	    B_bb0_ = -(mu + fpv(1, x2, y2));
	    if (ft == b0) return B_bb0_;
	    B_bb1_ = 1/2.0*(mu + fpv(2, x2, y2));
	    if (ft == b1) return B_bb1_;
	    a0 = 0;
	    if ( m2 != 0 ) a0 = m2*(1 - log(m2/mudim) + delta);
	    B_bb00_ = (2.0*(m1*B_bb0_ + a0) + m1 + m2)/8.0;
	    if (ft == b00) return B_bb00_;
	}
	else {
	    mu = log(m1/mudim) - delta;
	    f1 = fpv(1, y2, x2);
	    B_bb0_ = -(mu + f1);
	    if (ft == b0) return B_bb0_;
	    B_bb1_ = 1/2.0*(mu + (1.0 + x2)*f1 + 1/2.0);
	    if (ft == b1) return B_bb1_;
	    a0 = 0;
	    if ( m1 != 0 ) a0 = m1*(1 - log(m1/mudim) + delta);
	    B_bb00_ = (2.0*(m2*B_bb0_ + a0) + m1 + m2)/8.0;
	    if (ft == b00) return B_bb00_;
	}
    }
    else {
	mu = log(m2/mudim) - delta;
	B_bb0_ = -mu;
	if (ft == b0) return B_bb0_;
	B_bb1_ = 1/2.0*mu;
	if (ft == b1) return B_bb1_;
	B_bb00_ = .5*m1*(1.0 - mu);
	if (ft == b00) return B_bb00_;
    }

    abort();			// this must not happen
    return 0.0;
}

} // namespace AD

} // namespace

complex<double> A0(double m2, double scl2) noexcept
{
    return FF::A0(m2, scl2);
}

// LoopTools evaluates B functions (by default) using Ansgar Denner's
// implementation since version 2.8
complex<double> B0(double p2, double m2a, double m2b, double scl2) noexcept
{
    // see src/B/Bcoeff.F in LoopTools
    if (fabs(p2) + fabs(m2a) + fabs(m2b) < eps) return 0.0;

    return AD::B<AD::b0>(p2, m2a, m2b, scl2);
    // return FF::B0(p2, m2a, m2b, scl2);
}

complex<double> B1(double p2, double m2a, double m2b, double scl2) noexcept
{
    // see src/B/Bcoeff.F in LoopTools
    if (fabs(p2) + fabs(m2a) + fabs(m2b) < eps) return 0.0;

    return AD::B<AD::b1>(p2, m2a, m2b, scl2);
    // return FF::B1(p2, m2a, m2b, scl2);
}

complex<double> B00(double p2, double m2a, double m2b, double scl2) noexcept
{
    // see src/B/Bcoeff.F in LoopTools
    if (fabs(p2) + fabs(m2a) + fabs(m2b) < eps) return 0.0;

    return AD::B<AD::b00>(p2, m2a, m2b, scl2);
    // return FF::B00(p2, m2a, m2b, scl2);
}

complex<double> A0
(complex<double> m2, double scl2) noexcept
{
    return FF::A0(m2, scl2);
}

complex<double> B0
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    return FF::B0(p2, m2a, m2b, scl2);
}

complex<double> B1
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    return FF::B1(p2, m2a, m2b, scl2);
}

complex<double> B00
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2) noexcept
{
    return FF::B00(p2, m2a, m2b, scl2);
}

complex<double> D1B0(double p2, double m2a, double m2b) noexcept
{
    return FF::D1B0(p2, m2a, m2b);
}

complex<double> D1B0
(complex<double> p2, complex<double> m2a, complex<double> m2b) noexcept
{
    return FF::D1B0(p2, m2a, m2b);
}

#endif // defined(ENABLE_FFLITE)

double ReA0(double m2, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return A0(m2, scl2).real();
#else
    return softsusy::a0(sqrt(m2), sqrt(scl2));
#endif
}

double ReB0(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B0(p2, m2a, m2b, scl2).real();
#else
    return softsusy::b0(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
#endif
}

double ReB1(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B1(p2, m2a, m2b, scl2).real();
#else
    return -softsusy::b1(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
#endif
}

double ReB00(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B00(p2, m2a, m2b, scl2).real();
#else
    return softsusy::b22(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
#endif
}

double ReB22(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B22(p2, m2a, m2b, scl2).real();
#else
    return ReB00(p2, m2a, m2b, scl2) - ReA0(m2a, scl2)/4 - ReA0(m2b, scl2)/4;
#endif
}

double ReH0(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return H0(p2, m2a, m2b, scl2).real();
#else
    return 4*ReB00(p2, m2a, m2b, scl2) + ReG0(p2, m2a, m2b, scl2);
#endif
}

double ReF0(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return F0(p2, m2a, m2b, scl2).real();
#else
    return ReA0(m2a, scl2) - 2*ReA0(m2b, scl2)
	   - (2*p2 + 2*m2a - m2b) * ReB0(p2, m2a, m2b, scl2);
#endif
}

double ReG0(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return G0(p2, m2a, m2b, scl2).real();
#else
    return (p2 - m2a - m2b) * ReB0(p2, m2a, m2b, scl2)
	   - ReA0(m2a, scl2) - ReA0(m2b, scl2);
#endif
}

double ReD1B0(double p2, double m2a, double m2b) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
   return D1B0(p2, m2a, m2b).real();
#else
   return softsusy::d1_b0(p2, m2a, m2b);
#endif
}

double ReD1F0(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return D1F0(p2, m2a, m2b, scl2).real();
#else
    return - (2.0*p2 + 2.0*m2a - m2b) * ReD1B0(p2, m2a, m2b)
       - 2.0 * ReB0(p2, m2a, m2b, scl2);
#endif
}

double ReD1G0(double p2, double m2a, double m2b, double scl2) noexcept
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return D1G0(p2, m2a, m2b, scl2).real();
#else
    return (p2 - m2a - m2b) * ReD1B0(p2, m2a, m2b)
       + ReB0(p2, m2a, m2b, scl2);
#endif
}

} // namespace passarino_veltman

} // namespace flexiblesusy
