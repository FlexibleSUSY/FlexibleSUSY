#ifndef lattice_numerical_constraint_hpp
#define lattice_numerical_constraint_hpp


#include <gsl/gsl_deriv.h>

#include "lattice_foreign_constraint.hpp"

namespace flexiblesusy {

struct NumericalConstraintCommon {
    static constexpr Real default_epsilon = 1e-8;
};

class NumericalConstraint :
    public NumericalConstraintCommon,
    public ForeignConstraint {
public:
    NumericalConstraint
    (std::vector<size_t> dependence, Real epsilon = default_epsilon) :
	ForeignConstraint(1), nonzeros(dependence), deriv_epsilon(epsilon) {
	F_gsl.function = &c_wrap;
	F_gsl.params = this;
    }
    void init(RGFlow<Lattice> *flow, size_t theory, size_t site);
    void operator()();
    using ForeignConstraint::init;
protected:
    virtual Real c(const Real *x) const = 0;
private:
    std::vector<size_t> nonzeros;
    std::vector<bool> depends_on;
    gsl_function F_gsl;
    static double c_wrap(double xj, void *params);
    size_t j;
    Real deriv_epsilon;
};

class AnyNumericalConstraint : public NumericalConstraint {
public:
    AnyNumericalConstraint
    (std::vector<size_t> dependence,
     std::function<Real(const AnyNumericalConstraint *, const Real *x)> fxn,
     Real epsilon = default_epsilon) :
	NumericalConstraint(dependence, epsilon), fxn_(fxn)
	{}
protected:
    Real c(const Real *x) const { return fxn_(this, x); }
private:
    std::function<Real(const AnyNumericalConstraint *, const Real *x)> fxn_;
};

class NumericalMatching :
    public NumericalConstraintCommon,
    public ForeignMatching {
public:
    NumericalMatching(std::vector<size_t> depL, std::vector<size_t> depH,
		      Real epsilon = default_epsilon) :
	ForeignMatching(1), nonzerosL(depL), nonzerosH(depH),
	deriv_epsilon(epsilon) {
	F_gsl.function = &c_wrap;
	F_gsl.params = this;
    }
    void init(RGFlow<Lattice> *flow, size_t lower_theory);
    void operator()();
    using ForeignMatching::init;
protected:
    virtual Real c(const Real *w, const Real *x) const = 0;
private:
    std::vector<size_t> nonzerosL, nonzerosH;
    std::vector<bool> depends_on;
    gsl_function F_gsl;
    static double c_wrap(double wxj, void *params);
    size_t j;
    Real &wx(size_t j) { return j < w.size() ? w[j] : x[j - w.size()]; }
    Real deriv_epsilon;
};

class AnyNumericalMatching : public NumericalMatching {
public:
    AnyNumericalMatching
    (std::vector<size_t> depL, std::vector<size_t> depH,
     std::function<Real(const AnyNumericalMatching *,
			const Real *w, const Real *x)> fxn,
     Real epsilon = default_epsilon) :
	NumericalMatching(depL, depH, epsilon), fxn_(fxn)
	{}
protected:
    Real c(const Real *w, const Real *x) const { return fxn_(this, w, x); }
private:
    std::function<Real(const AnyNumericalMatching *,
		       const Real *w, const Real *x)> fxn_;
};

} // namespace flexiblesusy

#endif // lattice_numerical_constraint_hpp
