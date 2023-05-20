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

#ifndef LATTICE_MODEL_H
#define LATTICE_MODEL_H

#include <string>
#include <ostream>
#include <functional>
#include "mathdefs.hpp"

namespace flexiblesusy {

class Lattice;
template<class Method> class RGFlow;

struct Wilson {
    Wilson(size_t nvars) : width(nvars) {}
    virtual ~Wilson() = default;
    size_t width;		// 1 + number of Wilson coefficients
    // derivative wrt t == x[0]
    virtual Real  dx(const Real *x, size_t i) const = 0;
    // d dx[i] / d x[j] == d dx[i] / d y[j] / unit[j]
    virtual void ddx(const Real *x, size_t i, Real *ddx) const = 0;
};

struct ParWilson {
    ParWilson(size_t width_) : width(width_) {}
    virtual ~ParWilson() = default;
    size_t width;		// 1 + number of Wilson coefficients
    // derivative wrt t == x[0]
    virtual Real  dx(const Real a, const Real *x, size_t i) const = 0;
    // d dx[i] / d x[j] == d dx[i] / d y[j] / unit[j]
    virtual void ddx(const Real a, const Real *x, size_t i, Real *ddx) const=0;
};

class Lattice_model: public ParWilson {
public:
    Lattice_model(size_t width) : ParWilson(width) {}
    virtual ~Lattice_model() = default;
    virtual void init(RGFlow<Lattice> *flow, size_t theory)
    { f = flow; T = theory; }
    virtual void calculate_spectrum() = 0;
    virtual std::string name() const { return "unnamed"; }
    virtual int run_to(double, double eps = -1.0);
    virtual void print(std::ostream& out) const { out << "Model: " << name(); }
    friend std::ostream& operator<<
    (std::ostream& out, const Lattice_model& model) {
	model.print(out);
	return out;
    }

protected:
    RGFlow<Lattice> *f;
    size_t T;
};

class Lattice_translator {
public:
    template<class T>
    struct Var {
	Var(std::function<T()> get, std::function<void(T)> set) :
	    get_(get), set_(set) {}
	operator T() { return get_(); }
	T operator=(T value) { set_(value); return get_(); }
	T operator=(Var<T>& var) { return operator=(T(var)); }
	std::function<T()> get_;
	std::function<void(T)> set_;
    };

    Lattice_translator(RGFlow<Lattice> *flow, size_t theory, size_t site) :
	f(flow), T(theory), m(site)
	{}
    Real  u(size_t i);
    Real& y(size_t i);

private:
    RGFlow<Lattice> *f;
    size_t T;
    size_t m;
};

} // namespace flexiblesusy

#endif
