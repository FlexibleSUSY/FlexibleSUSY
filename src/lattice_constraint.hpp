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

#ifndef lattice_constraint_hpp
#define lattice_constraint_hpp


#include <Eigen/Dense>
#include "constraint.hpp"
#include "lattice_solver.hpp"

namespace flexiblesusy {

template <class T>
class Matching;

class Lattice;

class Lattice_constraint {
public:
    // Lattice_constraint() {}
    virtual ~Lattice_constraint() = default;
    virtual void init(RGFlow<Lattice> *flow) { f = flow; activate(); }
    virtual void deactivate();
    virtual void alloc_rows() = 0;
    virtual void free_rows() { rfree(); }
    virtual void operator()() = 0;
    virtual void relocate(const std::vector<std::vector<size_t>>& site_maps)=0;
    // size_t relocated_site(size_t m, size_t old_height, size_t new_height)
    // 	{ return size_t(0.5 + m*Real(new_height-1)/Real(old_height-1)); }
    RGFlow<Lattice> *f;
protected:
    virtual void activate();
    Real& A(size_t r, size_t T, size_t m, size_t j)
    { return f->A(rows[r]->rowSpec.realRow, T, m, j); }
    Real y(size_t T, size_t m, size_t i) const { return f->y(T, m, i); }
    Real& z(size_t r) { return f->z[rows[r]->rowSpec.realRow]; }
    Real u(size_t T, size_t i) const { return f->u(T, i); }
    Real x(size_t T, size_t m, size_t i) const { return f->x(T, m, i); }
    void ralloc(size_t nrows, size_t T, size_t m, size_t span);
    void rfree();
private:
    std::vector<RGFlow<Lattice>::EqRow *> rows;
};

using InterTheoryConstraint = Matching<Lattice>;

template<>
class Matching<Lattice> : public Lattice_constraint {
public:
    // InterTheoryConstraint() {}
    virtual void init
    (RGFlow<Lattice> *flow, size_t lower_theory)
    { Lattice_constraint::init(flow); TL = lower_theory; }
    virtual void relocate(const std::vector<std::vector<size_t>>&) {}
    using Lattice_constraint::init;
protected:
    Real& A(size_t r, size_t To, size_t j)
    { return Lattice_constraint::A(r, TL+To, m(To), j); }
    Real y(size_t To, size_t i) const
    { return Lattice_constraint::y(TL+To, m(To), i); }
    Real& z(size_t r) { return Lattice_constraint::z(r); }
    Real u(size_t To, size_t i) const
    { return Lattice_constraint::u(TL+To, i); }
    Real x(size_t To, size_t i) const
    { return Lattice_constraint::x(TL+To, m(To), i); }
    size_t m(size_t To) const { return To ? 0 : f->efts[TL].height - 1; }
    void ralloc(size_t nrows)
    { Lattice_constraint::ralloc(nrows, TL, m(0), 2); }
    size_t TL;
private:
};

class IntraTheoryConstraint : public Lattice_constraint {
public:
    // IntraTheoryConstraint() {}
    virtual void init
    (RGFlow<Lattice> *flow, size_t theory, size_t site, size_t span_) {
	Lattice_constraint::init(flow);
	T = theory; mbegin = site; span = span_;
    }
    virtual void relocate(const std::vector<size_t>& site_map) {
	size_t new_mbegin = site_map[mbegin];
	size_t new_span  = site_map[mbegin + span - 1] - new_mbegin + 1;
	mbegin = new_mbegin;
	span  = new_span;
    }
    virtual void relocate(const std::vector<std::vector<size_t>>& site_maps)
    { relocate(site_maps[T]); }
    size_t mbegin;
    using Lattice_constraint::init;
protected:
    Real& A(size_t r, size_t m, size_t j)
    { return Lattice_constraint::A(r, T, m, j); }
    Real y(size_t m, size_t i) const { return Lattice_constraint::y(T, m, i); }
    Real& z(size_t r) { return Lattice_constraint::z(r); }
    Real u(size_t i) const { return Lattice_constraint::u(T, i); }
    Real x(size_t m, size_t i) const { return Lattice_constraint::x(T, m, i); }
    void ralloc(size_t nrows, size_t m, size_t span)
    { Lattice_constraint::ralloc(nrows, T, m, span); }
    size_t T;
    size_t span;
private:
};

using SingleSiteConstraint = Constraint<Lattice>;

template<>
class Constraint<Lattice> : public IntraTheoryConstraint {
public:
    // SingleSiteConstraint() {}
    virtual void init(RGFlow<Lattice> *flow, size_t theory, size_t site)
    { IntraTheoryConstraint::init(flow, theory, site, 1); }
    virtual double get_scale() const { return std::exp(x(0))*f->scl0; }
    using IntraTheoryConstraint::init;
protected:
    Real& A(size_t r, size_t j)
    { return IntraTheoryConstraint::A(r, mbegin, j); }
    Real y(size_t i) const { return IntraTheoryConstraint::y(mbegin, i); }
    Real& z(size_t r) { return IntraTheoryConstraint::z(r); }
    Real x(size_t i) const { return IntraTheoryConstraint::x(mbegin, i); }
    void ralloc(size_t nrows)
    { IntraTheoryConstraint::ralloc(nrows, mbegin, 1); }
private:
};

class SingleSiteInterTheoryConstraint : public InterTheoryConstraint {
public:
    SingleSiteInterTheoryConstraint
    (SingleSiteConstraint *c, size_t T_offset) :
	ssc(c), To(T_offset) {}
    virtual void init(RGFlow<Lattice> *flow, size_t lower_theory) {
	InterTheoryConstraint::init(flow, lower_theory);
	ssc->init(f, TL+To, m(To));
    }
    void alloc_rows() { ssc->alloc_rows(); }
    void free_rows() { ssc->free_rows(); }
    void operator()() { (*ssc)(); }
    using InterTheoryConstraint::init;
private:
    SingleSiteConstraint *ssc;
    size_t To;
};

class AnySingleSiteConstraint : public SingleSiteConstraint {
public:
    AnySingleSiteConstraint
    (size_t nrows, std::function<void(AnySingleSiteConstraint *)> fxn) :
	nr(nrows), fxn_(fxn)
	{}
    void alloc_rows() { ralloc(nr); }
    void operator()() { fxn_(this); }
private:
    size_t nr;
    std::function<void(AnySingleSiteConstraint *)> fxn_;
};

class Lattice_RGE : public IntraTheoryConstraint {
public:
    // Lattice_RGE() {}
    void init(RGFlow<Lattice> *flow, size_t theory, size_t site, size_t span) {
	IntraTheoryConstraint::init(flow, theory, site, span);
	x.resize(f->efts[T].w->width);
	ddxm0i.resize(x.size());
	ddxm1i.resize(x.size());
    }
    void alloc_rows() {
	for (size_t n = 0; n < span - 1; n++)
	    ralloc(f->efts[T].w->width - 1, mbegin + n, 2);
    }
    void operator()();
    using IntraTheoryConstraint::init;

private:
    RVec x;
    Real dxm0i, dxm1i;
    RVec ddxm0i, ddxm1i;

    void calc_dxmi_ddxmi(size_t m, size_t i, Real& dxmi, RVec& ddxmi);
    void set_diff(size_t r, size_t m, size_t i);
};

class Lattice_RKRGE : public IntraTheoryConstraint {
public:
    template<class A, class V, class M>
    struct Adapter_ {
	Adapter_() : x(nullptr,0), D(nullptr,0,0) {}
	void set(A& xD, size_t width) {
	    v = &xD;
	    n = width;
	    new (&x) Eigen::Map<V>(v->data()  , n);
	    new (&D) Eigen::Map<M>(v->data()+n, n, n);
	}
	size_t n;
	Eigen::Map<V> x;
	Eigen::Map<M> D;
	A *v;
    };

    using Adapter = Adapter_<Eigen::ArrayXd, Eigen::VectorXd, Eigen::MatrixXd>;
    using const_Adapter = Adapter_<const Eigen::ArrayXd, const Eigen::VectorXd,
                                   const Eigen::MatrixXd>;

    // Lattice_RKRGE() {}
    void init(RGFlow<Lattice> *flow, size_t theory, size_t site, size_t span) {
	assert(span == 2);
	IntraTheoryConstraint::init(flow, theory, site, span);
	xD0.resize(f->efts[T].w->width * (1 + f->efts[T].w->width));
	xD1.resize(xD0.size());
	a0.set(xD0, f->efts[T].w->width);
	a1.set(xD1, f->efts[T].w->width);
	dx0.resize(f->efts[T].w->width);
	dx1.resize(f->efts[T].w->width);
    }
    void alloc_rows() {
	ralloc(f->efts[T].w->width - 1, mbegin, 2);
    }
    void operator()();
    using IntraTheoryConstraint::init;

private:
    int evolve_to(Real t_new, Adapter& a, Real eps = -1);

    Eigen::ArrayXd xD0, xD1;
    Adapter         a0,  a1;
    RVec           dx0, dx1;
};

class Uniform_dt : public IntraTheoryConstraint {
public:
    // Uniform_dt() {}
    void alloc_rows() {
	for (size_t n = 1; n < span - 1; n++)
	    ralloc(1, mbegin + n - 1, 3);
    }
    void operator()() {
	for (size_t n = 1; n < span - 1; n++) {
	    size_t m = mbegin + n;
	    size_t r = n - 1;
	    A(r,m-1,0) =  1;
	    A(r,m  ,0) = -2;
	    A(r,m+1,0) =  1;
	    z(r) = 0;
	}
    }
private:
};

class Match_t : public InterTheoryConstraint {
public:
    // Match_t() {}
    void alloc_rows() {
	ralloc(1);
    }
    void operator()() {
	A(0,0,0) =  u(0,0);
	A(0,1,0) = -u(1,0);
	z(0) = 0;
    }
};

class Unified_xi_xj : public AnySingleSiteConstraint {
public:
    Unified_xi_xj(size_t i, size_t j) :
	AnySingleSiteConstraint(1,
	    [=](AnySingleSiteConstraint *) {
		A(0,i) =  u(i);
		A(0,j) = -u(j);
		z(0) = 0;
	    })
	{}
};

class Fixed_x : public AnySingleSiteConstraint {
public:
    Fixed_x(size_t i, Real x) :
	AnySingleSiteConstraint(1,
	    [=](AnySingleSiteConstraint *) {
		A(0,i) = u(i);
		z(0) = x;
	    })
	{}
};

class Fixed_t : public AnySingleSiteConstraint {
public:
    Fixed_t(Real mu) :
	AnySingleSiteConstraint(1,
	    [=](AnySingleSiteConstraint *) {
		A(0,0) = u(0);
		z(0) = std::log(mu/f->scl0);
	    })
	{}
};

} // namespace flexiblesusy

#endif // lattice_constraint_hpp
