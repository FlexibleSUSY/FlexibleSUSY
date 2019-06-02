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

#ifndef lattice_solver_hpp
#define lattice_solver_hpp


#include <functional>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include "mathdefs.hpp"

#include "rg_flow.hpp"
#include "error.hpp"
#include "lattice_model.hpp"

namespace flexiblesusy {

#if 1
template <class T> class Constraint;
template <class T> class Matching;
template <class T> class Convergence_tester;
template <class T> class Initial_guesser;
#endif

class Lattice;
class Lattice_constraint;
// class SingleSiteConstraint;
// class InterTheoryConstraint;
using SingleSiteConstraint = Constraint<Lattice>;
using InterTheoryConstraint = Matching<Lattice>;

class Two_scale_running_precision;


using RVec = std::vector<Real>;


template<class T>
class band_matrix {
public:
    band_matrix(size_t N, size_t KL, size_t KU, size_t LD) :
	n(N), kl(KL), ku(KU), ld(LD) {
	assert(ld >= kl+ku+1);
	AB = new T[ld*n];
    }
    ~band_matrix() { delete[] AB; }
    T *pointer() { return AB; }
    void clear() { for (size_t i = 0; i < ld*n; i++) AB[i] = 0; }
    T& operator()(size_t i, size_t j) {
	assert(i < n && j < n && i+ku >= j && i <= j+kl);
	size_t idx = kl+ku+i-j + ld*j;
	assert(idx < ld*n);
	return AB[idx];
    }

private:
    size_t n, kl, ku, ld;
    T *AB;
};

			      
template<>
class RGFlow<Lattice> {
public:
    enum Inner_status { JUMPED=-2, ENDLESS=-1, CONVERGED=0 };

    class MemoryError : public Error {
    public:
	MemoryError(const std::string& message_) : Error(message_) {}
	virtual ~MemoryError() = default;
    };

    class SetupError : public Error {
    public:
	SetupError(const std::string& message_) : Error(message_) {}
	virtual ~SetupError() = default;
    };

    class NonInvertibleMatrixError : public Error {
    public:
	NonInvertibleMatrixError(const std::string& message_) :
	    Error(message_) {}
	virtual ~NonInvertibleMatrixError() = default;
    };

    class DivergenceError : public Error {
    public:
	DivergenceError(const std::string& message_) : Error(message_) {}
	virtual ~DivergenceError() = default;
    };

    class NoConvergenceError : public Error {
    public:
	NoConvergenceError(size_t number_of_iterations_)
	    : Error("RGFlow<Lattice>: no convergence")
            , number_of_iterations(number_of_iterations_) {}
	virtual ~NoConvergenceError() = default;
	std::string what_detailed() const override {
	    std::stringstream message;
	    message << "RGFlow<Lattice>::NoConvergenceError: no convergence"
		    << " after " << number_of_iterations << " iterations";
	    return message.str();
	}
    private:
	size_t number_of_iterations;
    };

    class NonPerturbativeRunningError : public Error {
    public:
	NonPerturbativeRunningError(Lattice_model* model_, double scale_)
	    : Error("non-perturbative RG running")
            , model(model_)
	    , scale(scale_)
	    {}
	virtual ~NonPerturbativeRunningError() = default;
	std::string what_detailed() const override;
    private:
	Lattice_model* model;
	double scale;
    };

    struct EFT {
	EFT(Lattice_model *model, RGFlow *flow) :
	    w(model), units(w->width, 1), f(flow)
	    {}
	Lattice_model *w;
	RVec units;		// x normalizations
	// std::vector<ptrdiff_t> r_offset;
	RGFlow *f;
	size_t height;
	size_t T;
	size_t offset;		// location within y_
    };

    union EqRow {
	EqRow *next;		// next free element
	struct {
	    size_t T, m;	// first site
	    size_t n;		// number of occupied sites
	    size_t realRow;	// row within A_
	} rowSpec;
    };

    // RGFlow(std::vector<EFT *> efts, std::vector<Constraint *> constraints,
    // 	   Real scale0, int verbose, std::ostream *lout);
    RGFlow();
    ~RGFlow();

    /// add a model and constraints
    /// order of constraints: ascending t
    void add_model(Lattice_model* model,
		   const std::vector<SingleSiteConstraint*>& constraints);
    /// add a model, constraints and matching condition
    /// order of constraints: ascending t
    void add_model(Lattice_model*,
		   InterTheoryConstraint *m = nullptr,
		   const std::vector<SingleSiteConstraint*>& constraints = std::vector<SingleSiteConstraint*>());
    /// add a model and up- and downwards constraints
    /// order of upward_constraints: ascending t
    /// order of downward_constraints: descending t
    void add_model(Lattice_model*,
		   const std::vector<SingleSiteConstraint*>& upward_constraints,
		   const std::vector<SingleSiteConstraint*>& downward_constraints);
    /// add a model, up- and downward constraints and matching condition
    /// order of upward_constraints: ascending t
    /// order of downward_constraints: descending t
    void add_model(Lattice_model*,
		   InterTheoryConstraint *m,
		   const std::vector<SingleSiteConstraint*>& upward_constraints,
		   const std::vector<SingleSiteConstraint*>& downward_constraints);
    /// clear all internal data
    void reset();
    /// set convergence tester
    void set_convergence_tester(Convergence_tester<Lattice>*);
    /// set running precision calculator
    /// TODO: replace Two_scale_running_precision by something lattice
    void set_running_precision(Two_scale_running_precision*);
    void set_initial_guesser(Initial_guesser<Lattice>*);

    void enable_hybrid() { hybrid = true; }
    void disable_multithreading() { multithreading = false; };
    void solve();
    friend std::ostream& operator<<(std::ostream &out, const RGFlow& flow);

    struct EFTspec : public EFT {
	EFTspec(Lattice_model *model,
		const std::vector<SingleSiteConstraint*>& cs,
		InterTheoryConstraint* m,
		RGFlow *flow);

	// Lattice_model *model;
	std::vector<SingleSiteConstraint*> constraints;
	InterTheoryConstraint *matching;
    };
    std::vector<EFTspec> efts;
    Initial_guesser<Lattice> *init_profile;
    Real tiny_dy;
    Real huge_dy;
    Real a;			// continuation parameter btw 0 & 1
    size_t max_a_steps;
    size_t max_iter;
    bool units_set;
    bool hybrid;
    Real scl0;			// t == log(mu/scl0)
    RVec y_;			// entire RG flow
    band_matrix<Real> *A_;
    RVec z;			// RHS of linear equations
    std::vector<EqRow> row_pool;
    EqRow *free_row_list_head;
    int N, KL, KU, LDA;		// for lapack
    std::vector<int> IPIV;	// ditto
    int verb;
    // std::ostream *log;
    size_t site_offset(size_t T, size_t m) const {
	assert(m < efts[T].height);
	return efts[T].offset + m*efts[T].w->width;
    }
    Real& A(size_t r, size_t T, size_t m, size_t i)
    { return (*A_)(r, site_offset(T,m) + i); }
    Real& y(size_t T, size_t m, size_t i) { return y_[site_offset(T,m) + i]; }
    Real  y(size_t T, size_t m, size_t i) const
    { return y_[site_offset(T,m) + i]; }
    Real  u(size_t T, size_t i) const { return efts[T].units[i]; }
    Real  x(size_t T, size_t m, size_t i) const { return y(T,m,i) * u(T,i); }
    std::vector<Lattice_constraint*> constraints;
    std::vector<size_t> teqidx;
    std::vector<size_t> rgeidx;
    std::unordered_set<Lattice_constraint*> elementary_constraints;
    bool multithreading;
    // as of v1.52 boost::thread gives valgrind the impression that
    // 8 bytes are not freed at exit
    // see http://www.cplusplus.com/forum/unices/83480/
    boost::thread_group *threads;
    boost::barrier *threads_begin, *threads_end;
    bool keep_threads;
    void set_units();
    void apply_constraints();
    void apply_constraints_thread(Lattice_constraint *c);
    void create_threads();
    void join_threads();
    Real maxdiff(const RVec& y0, const RVec& y1);
    void init_lattice();
    void increase_a();
    void increase_density();
    void rk_stage();
    Inner_status iterate();
    std::vector<std::vector<size_t>> refine_lattice();
    void enable_Runge_Kutta();
    void disable_Runge_Kutta();
    void resample(const std::vector<std::vector<size_t>>& site_maps);
    EqRow *ralloc(size_t T, size_t m, size_t span);
    void rfree(EqRow *r);
    void init_free_row_list();
    void sort_rows();
};

} // namespace flexiblesusy

#endif // lattice_solver_hpp
