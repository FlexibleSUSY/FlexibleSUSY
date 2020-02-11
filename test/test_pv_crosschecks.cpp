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

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <limits>
#include <vector>
#include <map>
#include <cstdlib>
#include <cmath>
#include "pv.hpp"
#include "wrappers.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pv_crosschecks

#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;

using namespace std;
using namespace flexiblesusy;

const double eps = numeric_limits<double>::min();

template<class T, class U>
struct LoopFunc {
    LoopFunc(size_t nargs_) : nargs(nargs_) {}
    virtual ~LoopFunc() {}
    const size_t nargs;
    virtual U eval(const vector<T>&, double scl2) = 0;
};

#define def_func(f, nargs, ...)				    \
template<class T, class U>				    \
struct f : LoopFunc<T, U> {				    \
    f() : LoopFunc<T, U>(nargs) {}			    \
    U eval(const vector<T>& a, double scl2)		    \
    { return passarino_veltman::f(__VA_ARGS__, scl2); }	    \
};

#define def_func3(f) def_func(f, 3, a[0], a[1], a[2])

#if defined(ENABLE_FFLITE) || defined(ENABLE_LOOPTOOLS)

def_func (  A0, 1, a[0])
def_func3(  B0 )
def_func3(  B1 )
def_func3(  B00)
def_func3(  B22)
def_func3(  H0 )
def_func3(  F0 )
def_func3(  G0 )

map<string, LoopFunc<complex<double>, complex<double> >*> ccfuncs;
map<string, LoopFunc<	     double , complex<double> >*> rcfuncs;

#endif

def_func (ReA0, 1, a[0])
def_func3(ReB0 )
def_func3(ReB1 )
def_func3(ReB00)
def_func3(ReB22)
def_func3(ReH0 )
def_func3(ReF0 )
def_func3(ReG0 )

map<string, LoopFunc<	     double ,         double  >*> rrfuncs;

struct Initialize_funcs {
    Initialize_funcs() {

#if defined(ENABLE_FFLITE) || defined(ENABLE_LOOPTOOLS)
	ccfuncs[  "A0" ] = new   A0 <complex<double>, complex<double> >;
	ccfuncs[  "B0" ] = new   B0 <complex<double>, complex<double> >;
	ccfuncs[  "B1" ] = new   B1 <complex<double>, complex<double> >;
	ccfuncs[  "B00"] = new   B00<complex<double>, complex<double> >;
	ccfuncs[  "B22"] = new   B22<complex<double>, complex<double> >;
	ccfuncs[  "F0" ] = new   F0 <complex<double>, complex<double> >;
	ccfuncs[  "G0" ] = new   G0 <complex<double>, complex<double> >;
	ccfuncs[  "H0" ] = new   H0 <complex<double>, complex<double> >;

	rcfuncs[  "A0" ] = new   A0 <        double , complex<double> >;
	rcfuncs[  "B0" ] = new   B0 <        double , complex<double> >;
	rcfuncs[  "B1" ] = new   B1 <        double , complex<double> >;
	rcfuncs[  "B00"] = new   B00<        double , complex<double> >;
	rcfuncs[  "B22"] = new   B22<        double , complex<double> >;
	rcfuncs[  "F0" ] = new   F0 <        double , complex<double> >;
	rcfuncs[  "G0" ] = new   G0 <        double , complex<double> >;
	rcfuncs[  "H0" ] = new   H0 <        double , complex<double> >;
#endif
	rrfuncs["ReA0" ] = new ReA0 <        double ,         double  >;
	rrfuncs["ReB0" ] = new ReB0 <        double ,         double  >;
	rrfuncs["ReB1" ] = new ReB1 <        double ,         double  >;
	rrfuncs["ReB00"] = new ReB00<        double ,         double  >;
	rrfuncs["ReB22"] = new ReB22<        double ,         double  >;
	rrfuncs["ReF0" ] = new ReF0 <        double ,         double  >;
	rrfuncs["ReG0" ] = new ReG0 <        double ,         double  >;
	rrfuncs["ReH0" ] = new ReH0 <        double ,         double  >;

    }
} initialize_funcs;

template<class T, class U> T random_m2_phase();

template<> double random_m2_phase<double, double>()
{
    return 1;
}

template<> double random_m2_phase<double, complex<double> >()
{
    return 1;
}

template<>
complex<double> random_m2_phase<complex<double>, complex<double> >()
{
    return polar(1.0, 2 * M_PI * rand() / RAND_MAX);
}

template<class T, class U> T random_m2(double order)
{
    double o = 4.0 * rand() / RAND_MAX - 2.0;
    return pow(10, order + o) * random_m2_phase<T, U>();
}

template<>
complex<double> random_m2<complex<double>, complex<double> >(double order)
{
    double o = 4.0 * rand() / RAND_MAX - 2.0;
    complex<double> m2 = pow(10, order + o) *
	random_m2_phase<complex<double>, complex<double> >();
    if (m2.real() < 0 && m2.imag() == 0)
	return complex<double>(m2.real(), -eps);
    else
	return m2;
}

template<class T, class U> T random_p2_phase();

template<> double random_p2_phase<double, double>()
{
    return 1;
}

template<> double random_p2_phase<double, complex<double> >()
{
    return 1 - rand() % 2 * 2;
}

template<>
complex<double> random_p2_phase<complex<double>, complex<double> >()
{
    return 1 - rand() % 2 * 2;
}

template<class T, class U> T random_p2(double order)
{
    double o = 4.0 * rand() / RAND_MAX - 2.0;
    return pow(10, order + o) * random_p2_phase<T, U>();
}

template<class T, class U>
void generate_points(map<string, LoopFunc<T, U>*>& funcs, size_t n)
{
    for (typename map<string, LoopFunc<T, U>*>::iterator f = funcs.begin();
	 f != funcs.end(); f++) {
	vector<T> args(f->second->nargs);
	for (size_t count = n; count; count--) {
	    double order = 40.0 * rand() / RAND_MAX;
	    size_t i = 0;
	    if (args.size() > 1) args[i++] = random_p2<T, U>(order);
	    for (; i < args.size(); i++) args[i] = random_m2<T, U>(order);
	    double s2 = random_m2<double, double>(order);
	    cout << "CHECK " << f->first;
	    for (size_t i = 0; i < args.size(); i++)
		cout << ' ' << args[i];
	    cout << ' ' << s2;
	    cout << ' ' << f->second->eval(args, s2) << '\n';
	}
    }
}

template<class T>
void check_close_fraction(T a, T b, double tol, double scl=1);

template<>
void check_close_fraction<double>(double a, double b, double tol, double scl)
{
    if      (a == 0)
	BOOST_CHECK_SMALL(b / scl, tol);
    else if (b == 0)
	BOOST_CHECK_SMALL(a / scl, tol);
    else
	BOOST_CHECK_CLOSE_FRACTION(a, b, tol);
}

template<>
void check_close_fraction<complex<double> >
(complex<double> a, complex<double> b, double tol, double)
{
    double scale = sqrt(norm(a) + norm(b));
    check_close_fraction(a.real(), b.real(), tol, scale);
    check_close_fraction(a.imag(), b.imag(), tol, scale);
}

template<class T, class U>
bool check_point(map<string, LoopFunc<T, U>*>& funcs, string& line, double tol)
{
    istringstream ls(line);
    string word;
    ls >> word; if (word != "CHECK") return false;

    ls >> word;
    typename map<string, LoopFunc<T, U>*>::iterator f = funcs.find(word);
    if (f == funcs.end()) return false;

    vector<T> args(f->second->nargs);
    for (size_t i = 0; i < args.size(); i++)
	if (!(ls >> args[i])) return false;
    double s2;
    U value_from_line;
    if (!(ls >> s2 >> value_from_line)) return false;

    U value_evaled = f->second->eval(args, s2);
    // cout << scientific << setprecision(17)
    // 	 << line << "  " << value_evaled << '\n';
    check_close_fraction(value_from_line, value_evaled, tol);

    return true;
}

#ifdef ENABLE_FFLITE

BOOST_AUTO_TEST_CASE(fflite_generate_points)
{
    cout << scientific << setprecision(17);
    const size_t n = 10000;
    generate_points<complex<double>, complex<double> >(ccfuncs, n);
    generate_points<        double , complex<double> >(rcfuncs, n);
    generate_points<        double ,         double  >(rrfuncs, n);
}

#elif defined(ENABLE_LOOPTOOLS)

BOOST_AUTO_TEST_CASE(looptools_check_points)
{
    string line;
    while (getline(cin, line))
	check_point(rrfuncs, line, 1e-15) ||
	check_point(rcfuncs, line, 1e-15) ||
	check_point(ccfuncs, line, 1e-15);
}

#else

BOOST_AUTO_TEST_CASE(softsusy_check_points)
{
    string line;
    while (getline(cin, line))
	check_point(rrfuncs, line, 3e-7);
}

#endif
