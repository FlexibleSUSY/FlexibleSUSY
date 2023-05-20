#include <iomanip>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_spline.h>

#include "mathdefs.hpp"
#include "lattice_model.hpp"
#include "lattice_constraint.hpp"
#include "lattice_initial_guesser.hpp"
#include "lattice_solver.hpp"
#include "rk.hpp"
#include "logger.hpp"


namespace flexiblesusy {

using namespace std;


ostream& operator<<(ostream &out, const RGFlow<Lattice>& f)
{
    for (size_t T = 0; T < f.efts.size(); T++) {
	out << "EFT " << T << '\n';
	for (size_t m = 0; m < f.efts[T].height; m++) {
	    out << m;
	    for (size_t i = 0; i < f.efts[T].w->width; i++)
		out << ' ' << f.y(T,m,i)*f.efts[T].units[i];
	    out << '\n';
	}
    }
    return out;
}

int Lattice_model::run_to(double, double /* eps */)
{
    // TODO: slide scale pointer
    return 0;
}

RGFlow<Lattice>::EFTspec::EFTspec
(Lattice_model *model, const vector<SingleSiteConstraint*>& cs,
 InterTheoryConstraint* m, RGFlow *flow) :
    EFT(model, flow), constraints(cs), matching(m)
{
}

RGFlow<Lattice>::RGFlow() :
    tiny_dy(1e-2), huge_dy(100),
    max_a_steps(1024), max_iter(100),
    units_set(false), hybrid(false), scl0(1),
    multithreading(true)
{
}

RGFlow<Lattice>::~RGFlow()
{
    delete A_;
    for (auto i: teqidx) delete constraints[i];
    for (auto i: rgeidx) delete constraints[i];
}

void RGFlow<Lattice>::add_model
(Lattice_model *model, const vector<SingleSiteConstraint*>& constraints)
{
    add_model(model, nullptr, constraints);
}

void RGFlow<Lattice>::add_model
(Lattice_model *model, InterTheoryConstraint *matching,
 const vector<SingleSiteConstraint*>& constraints)
{
    efts.emplace_back(model, constraints, matching, this);
}

void RGFlow<Lattice>::add_model
(Lattice_model *model,
 const vector<SingleSiteConstraint*>& upward_constraints,
 const vector<SingleSiteConstraint*>& downward_constraints)
{
    add_model(model, nullptr, upward_constraints, downward_constraints);
}

void RGFlow<Lattice>::add_model
(Lattice_model *model,
 InterTheoryConstraint *matching,
 const vector<SingleSiteConstraint*>& us,
 const vector<SingleSiteConstraint*>& ds)
{
    using ci_t = vector<SingleSiteConstraint*>::const_iterator;

    vector<SingleSiteConstraint*> rs(ds.size());
    reverse_copy(ds.begin(), ds.end(), rs.begin());

    vector<SingleSiteConstraint*> combined;
    auto pu = us.begin();
    ci_t pr = rs.begin();

    while (pu != us.cend() || pr != rs.cend()) {
	auto pun = find(pu, us.cend(), pr != rs.cend() ? *pr : nullptr);
	auto prn = find(pr, rs.cend(), pu != us.cend() ? *pu : nullptr);

	if (pun == pu && prn == pr) {
	    assert(*pu == *pr);
	    combined.push_back(*pu);
	    pu++; pr++;
	}
	else if ((pr  == rs.cend() && pu  != us.cend()) ||
		 (prn == rs.cend() && pun != us.cend()))
	    while (pu != pun)
		combined.push_back(*pu++);
	else if ((pu  == us.cend() && pr  != rs.cend()) ||
		 (pun == us.cend() && prn != rs.cend()))
	    while (pr != prn)
		combined.push_back(*pr++);
	else {
	    stringstream msg;
	    msg << "RGFlow<Lattice>::Error: failed to combine upward "
		   "and downward constraints in model" << model->name();
	    throw SetupError(msg.str());
	}
    }
#ifdef ENABLE_VERBOSE
    {
	stringstream uss, rss, css;
	for (auto c: us	     ) uss << ' ' << c;
	for (auto c: rs	     ) rss << ' ' << c;
	for (auto c: combined) css << ' ' << c;
	VERBOSE_MSG("         upward   constraints:" << uss.str() << "\n"
		    "reversed downward constraints:" << rss.str() << "\n"
		    "         combined constraints:" << css.str());
    }
#endif
    add_model(model, matching, combined);
}

void RGFlow<Lattice>::reset()
{
    // TODO: do something reasonable
}

void RGFlow<Lattice>::set_convergence_tester(Convergence_tester<Lattice>*)
{
    // TODO: do something reasonable
}

void RGFlow<Lattice>::set_running_precision(Two_scale_running_precision*)
{
    // TODO: do something reasonable
}

void RGFlow<Lattice>::set_initial_guesser(Initial_guesser<Lattice>* guesser)
{
    init_profile = guesser;
    init_profile->init(this);
}

void RGFlow<Lattice>::solve()
{
    if (efts.empty())
	throw SetupError("RGFlow<Lattice>::Error: EFT tower empty");

    init_lattice();
    create_threads();
    increase_a();
    if (hybrid) rk_stage();
    else increase_density();
    join_threads();
}

void RGFlow<Lattice>::init_lattice()
{
    size_t max_width = 0;
    size_t offset = 0;
    for (size_t T = 0; T < efts.size(); T++) {
	efts[T].w->init(this, T);
	size_t height = !!T;
	for (auto c: efts[T].constraints)
	    c->init(this, T, height++);
	if (efts[T].matching) {
	    efts[T].matching->init(this, T);
	    height++;
	}
	efts[T].height = height;
	efts[T].T = T;
	efts[T].offset = offset;
	size_t width = efts[T].w->width;
	if (width > max_width) max_width = width;
	offset += width * height;
    }
    y_.resize(offset);
    z.resize(offset);
    row_pool.resize(offset);
    init_free_row_list();

    for (size_t T = 0; T < efts.size(); T++) {
	for (size_t m = 0; m < efts[T].height-1; m++) {
	    teqidx.push_back(constraints.size());
	    auto *teq = new Uniform_dt;
	    teq->init(this, T, m, 2);
	    constraints.push_back(teq);

	    rgeidx.push_back(constraints.size());
	    auto *rge = new Lattice_RGE;
	    rge->init(this, T, m, 2);
	    constraints.push_back(rge);
	}
	for (auto c: efts[T].constraints)
	    constraints.push_back(c);
	if (efts[T].matching)
	    constraints.push_back(efts[T].matching);
    }

    for (auto c: constraints) c->alloc_rows();
    sort_rows();

    N = y_.size();
    KL = KU = 2*max_width;
    LDA = 2*KL + KU + 1;
    IPIV.resize(N);
    A_ = new band_matrix<Real>(N, KL, KU, LDA);
    if (A_ == nullptr)
	throw MemoryError("RGFlow<Lattice>::Error: failed to allocate matrix");
}

void RGFlow<Lattice>::increase_a()
{
    enum { END, INIT, RUN } state = INIT;
    size_t a_steps = 1;

    VERBOSE_MSG("\n\nentering a-loop with a_steps=" << a_steps);
    while (state != END) {
	(*init_profile)();
	if (state == INIT) {
	    VERBOSE_MSG("initial RG flow:\n" << *this);
	    state = RUN;
	}

	for (size_t a_i = 0; a_i <= a_steps; a_i++) {
	    a = Real(a_i)/a_steps;
	    VERBOSE_MSG("\n\nentering Newton loop with a=" <<
			a_i << '/' << a_steps);
	    Inner_status s = iterate();
	    switch (s) {
	    case JUMPED:
		a_steps *= 2;
		if (a_steps <= max_a_steps) {
		    VERBOSE_MSG("RG flow jumped, restarting a-loop "
				"with doubled a_steps=" << a_steps);
		    goto end_outer;
		}
		else
		    throw DivergenceError
			("RGFlow<Lattice>::Error: Newton iterations diverged");
	    case ENDLESS:
		a_steps *= 2;
		if (a_steps <= max_a_steps) {
		    VERBOSE_MSG("Newton iterations do not converge, "
				"restarting a-loop "
				"with doubled a_steps=" << a_steps);
		    goto end_outer;
		}
		else
		    throw NoConvergenceError(max_iter);
	    default:
		VERBOSE_MSG(*this);
		break;
	    }
	}
	VERBOSE_MSG("exiting a-loop");
	state = END;
    end_outer:
	;
    }
}

void RGFlow<Lattice>::increase_density()
{
    VERBOSE_MSG("\n\nbeginning to increase lattice density");

    for (;;) {
	RVec old_y = y_;
	vector<size_t> old_heights(efts.size());
	vector<size_t> old_offsets(efts.size());
	for (size_t T = 0; T < efts.size(); T++) {
	    old_heights[T] = efts[T].height;
	    old_offsets[T] = efts[T].offset;
	}
	vector<vector<size_t>> site_maps = refine_lattice();
	stringstream heights;
#ifdef ENABLE_VERBOSE
	for (auto& e: efts) heights << ' ' << e.height;
	VERBOSE_MSG("\n\nentering Newton loop with heights:" <<
		    heights.str());
#endif
	Inner_status s = iterate();
	switch (s) {
	case JUMPED:
	    throw DivergenceError("RGFlow<Lattice>::Error: RG flow jumped, "
				  "this cannot happen");
	case ENDLESS:
	    throw NoConvergenceError(max_iter);
	default:
	    Real maxdy = 0;
	    for (size_t T = 0; T < efts.size(); T++)
		for (size_t m = 0; m < old_heights[T]; m++)
		    for (size_t i = 0; i < efts[T].w->width; i++) {
			Real diff =
			    fabs(y(T,site_maps[T][m],i) -
				 old_y[old_offsets[T]+m*efts[T].w->width+i]);
			if (diff > maxdy) maxdy = diff;
		    }
	    VERBOSE_MSG("maxdy=" << maxdy);
	    if (maxdy < tiny_dy) {
		VERBOSE_MSG("discretization fine enough");
		return;
	    }
	    break;
	}
    }
}

vector<vector<size_t>> RGFlow<Lattice>::refine_lattice()
{
    size_t T_top = efts.size() - 1;
    Real t_top = x(T_top,efts[T_top].height-1,0);
    Real t_bot = x(0,0,0);
    Real min_Dt = t_top - t_bot + 1;
    Real max_Dt = 0;
    for (size_t T = 0; T < efts.size(); T++) {
	for (size_t m = 0; m < efts[T].height-1; m++) {
	    Real Dt = x(T,m+1,0) - x(T,m,0);
	    if (Dt > max_Dt) max_Dt = Dt;
	    if (Dt < min_Dt) min_Dt = Dt;
	}
    }

    Real spacing = max_Dt > 2*min_Dt ? min_Dt : min_Dt/2;
    vector<vector<size_t>> site_maps(efts.size());
    for (size_t T = 0; T < efts.size(); T++) {
	size_t new_m = 0;
	site_maps[T].push_back(new_m);
	for (size_t m = 0; m < efts[T].height-1; m++) {
	    Real Dt = x(T,m+1,0) - x(T,m,0);
	    size_t q = Dt / spacing + 0.5;
	    site_maps[T].push_back(new_m += q);
	}
    }

    resample(site_maps);

    return site_maps;
}

void RGFlow<Lattice>::rk_stage()
{
    VERBOSE_MSG("\n\nentering Runge-Kutta stage");
    enable_Runge_Kutta();
    VERBOSE_MSG("\n\nentering Newton loop");
    Inner_status s = iterate();
    switch (s) {
    case JUMPED:
	throw DivergenceError("RGFlow<Lattice>::Error: RG flow jumped, "
			      "this cannot happen");
    case ENDLESS:
	throw NoConvergenceError(max_iter);
    default:
	break;
    }
    VERBOSE_MSG("exiting Runge-Kutta stage");
}

void RGFlow<Lattice>::enable_Runge_Kutta()
{
    join_threads();

    vector<vector<bool>> original(efts.size());
    for (size_t T = efts.size(); T--; ) {
	original[T].resize(efts[T].height);
	for (auto c: efts[T].constraints)
	    original[T][c->mbegin] = true;
	if (efts[T].matching)
	    original[T][efts[T].height-1] = original[T+1][0] = true;
    }

    vector<vector<size_t>> site_maps(efts.size());
    for (size_t T = 0; T < efts.size(); T++) {
	site_maps[T].resize(efts[T].height);
	for (size_t m = 0, new_m = 0; m < efts[T].height; m++) {
	    site_maps[T][m] = new_m;
	    if (original[T][m]) new_m++;
	}
    }

    resample(site_maps);

    VERBOSE_MSG("switching to Runge-Kutta RGEs");
    for (size_t i = 0, T = 0; T < efts.size(); T++)
	for (size_t m = 0; m < efts[T].height-1; m++, i++) {
	    constraints[rgeidx[i]]->free_rows();
	    constraints[rgeidx[i]]->deactivate();
	    delete constraints[rgeidx[i]];
	    auto *rkrge = new Lattice_RKRGE;
	    rkrge->init(this, T, m, 2);
	    constraints[rgeidx[i]] = rkrge;
	}
    for (auto i: rgeidx) constraints[i]->alloc_rows();
    sort_rows();

    create_threads();
}

void RGFlow<Lattice>::disable_Runge_Kutta()
{
    join_threads();

    VERBOSE_MSG("switching to difference RGEs");
    for (size_t i = 0, T = 0; T < efts.size(); T++)
	for (size_t m = 0; m < efts[T].height-1; m++, i++) {
	    constraints[rgeidx[i]]->free_rows();
	    constraints[rgeidx[i]]->deactivate();
	    delete constraints[rgeidx[i]];
	    auto *rge = new Lattice_RGE;
	    rge->init(this, T, m, 2);
	    constraints[rgeidx[i]] = rge;
	}
    for (auto i: rgeidx) constraints[i]->alloc_rows();
    sort_rows();

    create_threads();
}

void RGFlow<Lattice>::resample(const vector<vector<size_t>>& site_maps)
{
    vector<size_t> new_heights(efts.size());
    vector<size_t> new_offsets(efts.size());
    size_t offset = 0;
    for (size_t T = 0; T < efts.size(); T++) {
	new_heights[T] = site_maps[T][efts[T].height-1] + 1;
	new_offsets[T] = offset;
	offset += efts[T].w->width * new_heights[T];
    }
#ifdef ENABLE_VERBOSE
    for (size_t T = 0; T < efts.size(); T++) {
	stringstream maps;
	maps << "EFT " << T << ":";
	for (size_t m = 0; m < efts[T].height; m++)
	    maps << ' ' << m << "->" << site_maps[T][m];
	VERBOSE_MSG(maps.str());
    }
#endif

    RVec new_y(offset);
    for (size_t T = 0; T < efts.size(); T++) {
	RVec y0(efts[T].height);
	for (size_t m = 0; m < efts[T].height; m++) {
	    size_t new_m = site_maps[T][m];
	    y0[m] = new_y[new_offsets[T] + new_m*efts[T].w->width] = y(T,m,0);
	}
	for (size_t m = 0; m < efts[T].height-1; m++) {
	    size_t q = site_maps[T][m+1] - site_maps[T][m];
	    for (size_t n = 1; n < q; n++) {
		size_t new_m = site_maps[T][m] + n;
		new_y[new_offsets[T] + new_m*efts[T].w->width] =
		    ((q-n)*y(T,m,0) + n*y(T,m+1,0)) / q;
	    }
	}
	RVec yi(efts[T].height);
	for (size_t i = 1; i < efts[T].w->width; i++) {
	    gsl_interp_accel *acc = gsl_interp_accel_alloc();
	    const gsl_interp_type *itype =
		efts[T].height > 2 ? gsl_interp_cspline : gsl_interp_linear;
	    gsl_spline *spl = gsl_spline_alloc(itype, efts[T].height);
	    for (size_t m = 0; m < efts[T].height; m++) yi[m] = y(T,m,i);
	    gsl_spline_init(spl, &y0[0], &yi[0], efts[T].height);
	    for (size_t m = 0; m < new_heights[T]; m++)
		new_y[new_offsets[T] + m*efts[T].w->width + i] =
		    gsl_spline_eval
		    (spl, new_y[new_offsets[T] + m*efts[T].w->width], acc);
	    gsl_spline_free(spl);
	    gsl_interp_accel_free(acc);
	}
    }

    for (auto c: constraints) c->free_rows();
    y_ = new_y;
    z.resize(offset);
    row_pool.resize(offset);
    init_free_row_list();
    for (auto c: constraints) c->relocate(site_maps);
    for (size_t T = 0; T < efts.size(); T++) {
	efts[T].height = new_heights[T];
	efts[T].offset = new_offsets[T];
    }
    for (auto c: constraints) c->alloc_rows();
    sort_rows();

    N = y_.size();
    IPIV.resize(N);
    delete A_;
    A_ = new band_matrix<Real>(N, KL, KU, LDA);
    if (A_ == nullptr)
	throw MemoryError("RGFlow<Lattice>::Error: failed to allocate matrix");
}

void RGFlow<Lattice>::set_units()
{
    for (size_t T = 0; T < efts.size(); T++) {
	for (size_t i = 0; i < efts[T].w->width; i++) {
	    Real miny = y(T,0,i);
	    Real maxy = miny;
	    for (size_t m = 0; m < efts[T].height; m++) {
		if (y(T,m,i) < miny) miny = y(T,m,i);
		if (y(T,m,i) > maxy) maxy = y(T,m,i);
	    }
	    efts[T].units[i] = max(maxy-miny,max(fabs(maxy),fabs(miny)));
	    if (efts[T].units[i] == 0) {
		VERBOSE_MSG("no hint about unit of x[" << i << "] in EFT "
			    << T << ", defaulting to 1");
		efts[T].units[i] = 1;
	    }
	    for (size_t m = 0; m < efts[T].height; m++)
		y(T,m,i) /= efts[T].units[i];
	}
    }
}

void RGFlow<Lattice>::apply_constraints()
{
    A_->clear();

    if (multithreading) {
	threads_begin->wait();
	(**elementary_constraints.begin())();
	threads_end->wait();
    }
    else
	for (auto c: elementary_constraints) (*c)();

    // for (auto c: constraints) (*c)();
}

void RGFlow<Lattice>::apply_constraints_thread(Lattice_constraint *c)
{
    while (threads_begin->wait(), keep_threads) {
	(*c)();
	threads_end->wait();
    }
}

void RGFlow<Lattice>::create_threads()
{
    if (!multithreading) return;

    threads_begin = new boost::barrier(elementary_constraints.size());
    threads_end   = new boost::barrier(elementary_constraints.size());
    keep_threads = true;
    threads = new boost::thread_group;
    for (auto c: elementary_constraints)
	if (c != *elementary_constraints.begin())
	    threads->add_thread(
		new boost::thread
		(&RGFlow<Lattice>::apply_constraints_thread, this, c));
    VERBOSE_MSG("launched " << threads->size() << " subthreads");
}

void RGFlow<Lattice>::join_threads()
{
    if (!multithreading) return;

    keep_threads = false;
    threads_begin->wait();
    threads->join_all();
    VERBOSE_MSG("joined " << threads->size() << " subthreads");
    delete threads;
    delete threads_begin;
    delete threads_end;
}

Real RGFlow<Lattice>::maxdiff(const RVec& y0, const RVec& y1)
{
    assert(y0.size() == y1.size());
    Real max = 0;
    auto p = y0.begin();
    auto q = y1.begin();
    while (p < y0.end()) {
	Real diff = fabs(*p++ - *q++);
	if (diff > max) max = diff;
    }
    return max;
}

RGFlow<Lattice>::EqRow *RGFlow<Lattice>::ralloc
(size_t T, size_t m, size_t span)
{
    EqRow *r = free_row_list_head;
    assert(r != nullptr);
    if (r != nullptr) {
	free_row_list_head = r->next;
	r->rowSpec = { T, m, span, 0/* to be determined by sort_rows() */ };
    }
    return r;
}

void RGFlow<Lattice>::rfree(RGFlow<Lattice>::EqRow *r)
{
    r->next = free_row_list_head;
    free_row_list_head = r;
}

void RGFlow<Lattice>::sort_rows()
{
    if (free_row_list_head != nullptr)
	throw SetupError("RGFlow<Lattice>::Error: constraints not enough");
    vector<EqRow *> layout(row_pool.size());
    for (size_t r = 0; r < layout.size(); r++)
	layout[r] = &row_pool[r];
    sort(layout.begin(), layout.end(),
	 [](EqRow *a, EqRow *b) -> bool {
	     if (a->rowSpec.T < b->rowSpec.T) return true;
	     if (a->rowSpec.T > b->rowSpec.T) return false;
	     if (a->rowSpec.m < b->rowSpec.m) return true;
	     if (a->rowSpec.m > b->rowSpec.m) return false;
	     if (a->rowSpec.n < b->rowSpec.n) return true;
	     return false;
	 });
    for (size_t r = 0; r < layout.size(); r++)
	layout[r]->rowSpec.realRow = r;
}

void RGFlow<Lattice>::init_free_row_list()
{
    for (size_t r = 0; r < row_pool.size() - 1; r++)
	row_pool[r].next = &row_pool[r+1];
    row_pool.back().next = nullptr;
    free_row_list_head = &row_pool[0];
}

extern "C" void dgbsv_
(const int& N, const int& KL, const int& KU, const int& NRHS,
 double *AB, const int& LDAB, int *IPIV, double *B, const int& LDB, int *INFO);

RGFlow<Lattice>::Inner_status RGFlow<Lattice>::iterate()
{
    if (!units_set) { set_units(); units_set = true; }

    // const Real epsdy = 1e-8;
    // const Real jumpdy = .5;	// TODO: find a better criterion
    size_t iter = 0;
    enum { END, INIT } state = INIT;

    int NRHS = 1;
    int LDB = N;
    int INFO;

    while (iter < max_iter && state != END) {
	iter++;
	apply_constraints();
#ifdef DUMP_MATRIX
	cout << setprecision(15) << scientific;
	for (size_t i = 0; i < y_.size(); i++) {
	    for (size_t j = 0; j < y_.size(); j++) {
		if (j) cout << " ";
		if (abs(signed(i)-signed(j)) < KL) cout << (*A_)(i,j);
		else cout << 0;
	    }
	    cout << "\n";
	}
#endif
	dgbsv_(N,KL,KU,NRHS,A_->pointer(),LDA,&IPIV[0],&z[0],LDB,&INFO);
	assert(INFO >= 0);
	if (INFO > 0) {
	    stringstream msg;
	    msg << "RGFlow<Lattice>::Error: failed to solve equations, "
		   "DGBSV returned INFO = " << INFO;
	    throw NonInvertibleMatrixError(msg.str());
	}
	Real maxdy = maxdiff(y_, z);
	VERBOSE_MSG("iter=" << iter << " maxdy=" << maxdy);
	if (maxdy < tiny_dy)
	    state = END;
	else if (maxdy > huge_dy)
	    return JUMPED;
	y_ = z;
    }

    if (state == END) {
	VERBOSE_MSG("converged after " << iter << " Newton iterations");
	return CONVERGED;
    }
    else
	return ENDLESS;
}

} // namespace flexiblesusy
