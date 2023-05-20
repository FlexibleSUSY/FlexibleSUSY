#ifndef LATTICE_FOREIGN_CONSTRAINT_H
#define LATTICE_FOREIGN_CONSTRAINT_H


#include "lattice_constraint.hpp"

namespace flexiblesusy {

class ForeignConstraint : public SingleSiteConstraint {
public:
    ForeignConstraint(size_t nrows) :
	SingleSiteConstraint(), nr(nrows) {}
    void init(RGFlow<Lattice> *flow, size_t theory, size_t site) {
	SingleSiteConstraint::init(flow, theory, site);
	x.resize(f->efts[T].w->width);
	row.resize(x.size());
    }
    void alloc_rows() { ralloc(nr); }
    using SingleSiteConstraint::init;
protected:
    void set_x() {
	for (size_t j = 0; j < x.size(); j++) x[j] = u(j)*y(j);
    }
    void copy_row(size_t r) {
	for (size_t j = 0; j < row.size(); j++) A(r,j) = row[j]*u(j);
	z(r) = rhs;
    }
    RVec x;
    RVec row;
    Real rhs;
private:
    size_t nr;
};

class ForeignMatching : public InterTheoryConstraint {
public:
    ForeignMatching(size_t nrows) :
	InterTheoryConstraint(), nr(nrows) {}
    void init(RGFlow<Lattice> *flow, size_t lower_theory) {
	InterTheoryConstraint::init(flow, lower_theory);
	w.resize(f->efts[TL  ].w->width);
	x.resize(f->efts[TL+1].w->width);
	row.resize(w.size() + x.size());
    }
    void alloc_rows() { ralloc(nr); }
    using InterTheoryConstraint::init;
protected:
    void set_w_x() {
	for (size_t j = 0; j < w.size(); j++) w[j] = u(0,j)*y(0,j);
	for (size_t j = 0; j < x.size(); j++) x[j] = u(1,j)*y(1,j);
    }
    void copy_row(size_t r) {
	RVec::const_iterator p = row.begin();
	for (size_t j = 0; j < w.size(); j++) A(r,0,j) = u(0,j) * *p++;
	for (size_t j = 0; j < x.size(); j++) A(r,1,j) = u(1,j) * *p++;
	z(r) = rhs;
    }
    RVec w;
    RVec x;
    RVec row;
    Real rhs;
private:
    size_t nr;
};

} // namespace flexiblesusy

#endif // LATTICE_FOREIGN_CONSTRAINT_H
