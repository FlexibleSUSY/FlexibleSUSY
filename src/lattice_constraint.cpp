#include <algorithm>

#include "lattice_constraint.hpp"
#include "rk.hpp"
#include "logger.hpp"

namespace flexiblesusy {

using namespace std;
using namespace Eigen;
using namespace runge_kutta;

void Lattice_constraint::activate()
{
    VERBOSE_MSG("registering elementary constraint " << this);
    f->elementary_constraints.insert(this);
}

void Lattice_constraint::deactivate()
{
    VERBOSE_MSG("deregistering elementary constraint " << this);
    f->elementary_constraints.erase(this);
}

void Lattice_constraint::ralloc
(size_t nrows, size_t T, size_t m, size_t span)
{
    VERBOSE_MSG("allocating " << nrows << " rows to " << this <<
		" at T=" << T << " m=" << m << " span=" << span);
    for (; nrows; nrows--)
	rows.push_back(f->ralloc(T, m, span));
}

void Lattice_constraint::rfree()
{
    VERBOSE_MSG("freeing " << rows.size() << " rows from " << this);
    for (auto r: rows) f->rfree(r);
    rows.clear();
}

void Lattice_RGE::operator()()
{
    for (size_t i = 1; i < x.size(); i++) { // x[0] untouched
	calc_dxmi_ddxmi(mbegin, i, dxm0i, ddxm0i);
	for (size_t n = 0; n < span - 1; n++) {
	    size_t m = mbegin + n;
	    calc_dxmi_ddxmi(m+1, i, dxm1i, ddxm1i);
	    // row order determined by alloc_rows()
	    set_diff(n * (f->efts[T].w->width-1) + i-1, m, i);
	    dxm0i = dxm1i;
	    swap(ddxm0i, ddxm1i); // cheaper than ddxm0i = ddxm1i
	}
    }
}

void Lattice_RGE::calc_dxmi_ddxmi(size_t m, size_t i, Real& dxmi, RVec& ddxmi)
{
    for (size_t j = 0; j < x.size(); j++) x[j] = y(m,j)*u(j);
    dxmi = f->efts[T].w->dx(f->a, &x[0], i);
    f->efts[T].w->ddx(f->a, &x[0], i, &ddxmi[0]);
}

void Lattice_RGE::set_diff(size_t r, size_t m, size_t i)
{
    A(r,m  ,0) = (ddxm0i[0]*(y(m+1,0)-y(m,0))*u(0)-(dxm1i+dxm0i))*u(0)/2;
    A(r,m+1,0) = (ddxm1i[0]*(y(m+1,0)-y(m,0))*u(0)+(dxm1i+dxm0i))*u(0)/2;
    z(r) = (y(m,0)*ddxm0i[0] + y(m+1,0)*ddxm1i[0])*u(0);
    for (size_t j = 1; j < x.size(); j++) {
	A(r,m  ,j) = (ddxm0i[j]*(y(m+1,0)-y(m,0))*u(0)/2+(i==j?1:0))*u(j);
	A(r,m+1,j) = (ddxm1i[j]*(y(m+1,0)-y(m,0))*u(0)/2-(i==j?1:0))*u(j);
	z(r) += (y(m,j)*ddxm0i[j] + y(m+1,j)*ddxm1i[j])*u(j);
    }
    z(r) *= (y(m+1,0)-y(m,0))*u(0)/2;
}

void Lattice_RKRGE::operator()()
{
    size_t m = mbegin;

    for (size_t j = 0; j < a0.n; j++) {
	a0.x(j) = y(m+0,j)*u(j);
	a1.x(j) = y(m+1,j)*u(j);
    }
    for (size_t j = 1; j < a0.n; j++) {
	dx0[j] = f->efts[T].w->dx(f->a, &a0.x(0), j);
	dx1[j] = f->efts[T].w->dx(f->a, &a1.x(0), j);
    }

    Real a_M = 0;
    Real t_M = (1-a_M)*a0.x(0) + a_M*a1.x(0);
    evolve_to(t_M, a0);
    evolve_to(t_M, a1);

    for (size_t i = 1; i < a0.n; i++) {
	size_t r = i - 1;
	z(r) = a1.x(i) - a0.x(i);
	for (size_t j = 1; j < a0.n; j++) {
	    A(r,m+0,0) -=  a0.D(i,j)*u(0)*dx0[j];
	    A(r,m+1,0) +=  a1.D(i,j)*u(0)*dx1[j];
	    A(r,m+0,j) =   a0.D(i,j)*u(j);
	    A(r,m+1,j) = - a1.D(i,j)*u(j);
	    z(r) += a0.D(i,j) * (u(j)*y(m+0,j) - u(0)*y(m+0,0)*dx0[j])
		  - a1.D(i,j) * (u(j)*y(m+1,j) - u(0)*y(m+1,0)*dx1[j]);
	}
    }
}

Real TOLERANCE = 1.0e-3;
const Real EPSTOL = 1.0e-11; ///< underflow accuracy

static int close(Real m1, Real m2, Real tol)
{
    return fabs(max(m1,m2)-min(m1,m2)) <= tol * max(fabs(m1),fabs(m2));
}

int Lattice_RKRGE::evolve_to(Real to, Adapter& a, Real eps)
{
    Real tol;
    if (eps < 0.0) tol = TOLERANCE;
    else if (eps < EPSTOL) tol = EPSTOL;
    else tol = eps;
    Real from = a.x(0);
    for (size_t i = 0; i < a.n; i++)
	for (size_t j = 0; j < a.n; j++)
	    a.D(i,j) = i==j ? 1 : 0;
    // from == to with high precision
    if (close(from, to, EPSTOL)) return 0;

    Real guess = (from - to) * 0.1; //first step size
    Real hmin = (from - to) * tol * 1.0e-5;

    RowVectorXd ddx(a.n);
    const_Adapter b;
    Adapter db;

    try {
       integrateOdes(*a.v, from, to, tol, guess, hmin,
	    [=,&ddx,&b,&db](Real, const ArrayXd& xD) -> ArrayXd {
		b.set(xD, a.n);

		ArrayXd dxD(xD.size());
		db.set(dxD, b.n);

		for (size_t i = 0; i < db.n; i++) {
		    db.x(i) = f->efts[T].w->dx(f->a, b.x.data(), i);

		    f->efts[T].w->ddx(f->a, b.x.data(), i, &ddx[0]);
		    db.D.block(i,0,1,db.n) = ddx * b.D;
		}

		return dxD;
	    }, odeStepper);

       a.x(0) = to;
    } catch(Error&) {
       return 1;
    }
    return 0;
}

} // namespace flexiblesusy
