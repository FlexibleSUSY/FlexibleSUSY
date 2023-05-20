#include "lattice_numerical_constraint.hpp"

namespace flexiblesusy {

using namespace std;


void NumericalConstraint::init
(RGFlow<Lattice> *flow, size_t theory, size_t site)
{
    ForeignConstraint::init(flow, theory, site);
    depends_on.resize(row.size());
    for (auto i: nonzeros) depends_on[i] = true;
    vector<size_t>().swap(nonzeros); // forces deallocation
}

void NumericalConstraint::operator()()
{
    set_x();
    Real BC0 = c(&x[0]);
    for (j = 0; j < row.size(); j++) {
	Real dc;
	Real xj = x[j];
	if (depends_on[j]) {
	    Real abserr;
	    gsl_deriv_central(&F_gsl, xj, deriv_epsilon*u(j), &dc, &abserr);
	    x[j] = xj;
	}
	else dc = 0;
	BC0 -= xj * (row[j] = dc);
    }
    rhs = -BC0;
    copy_row(0);
}

double NumericalConstraint::c_wrap(double xj, void *params)
{
    auto self = static_cast<NumericalConstraint *>(params);
    self->x[self->j] = xj;
    return self->c(&self->x[0]);
}

void NumericalMatching::init(RGFlow<Lattice> *flow, size_t lower_theory)
{
    ForeignMatching::init(flow, lower_theory);
    depends_on.resize(row.size());
    for (auto i: nonzerosL) depends_on[i         ] = true;
    for (auto i: nonzerosH) depends_on[i+w.size()] = true;
    vector<size_t>().swap(nonzerosL); // forces deallocation
    vector<size_t>().swap(nonzerosH); // forces deallocation
}

void NumericalMatching::operator()()
{
    set_w_x();
    Real MC0 = c(&w[0], &x[0]);
    for (j = 0; j < depends_on.size(); j++) {
	Real dc;
	Real wxj = wx(j);
	if (depends_on[j]) {
	    Real uj = j < w.size() ? u(0,j) : u(1,j-w.size());
	    Real abserr;
	    gsl_deriv_central(&F_gsl, wxj, deriv_epsilon*uj, &dc, &abserr);
	    wx(j) = wxj;
	}
	else dc = 0;
	MC0 -= wxj * (row[j] = dc);
    }
    rhs = -MC0;
    copy_row(0);
}

double NumericalMatching::c_wrap(double wxj, void *params)
{
    auto self = static_cast<NumericalMatching *>(params);
    self->wx(self->j) = wxj;
    return self->c(&self->w[0], &self->x[0]);
}

} // namespace flexiblesusy
