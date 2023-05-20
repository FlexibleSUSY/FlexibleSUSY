#include "models/MSSM/MSSM_input_parameters.hpp"
#include "models/MSSM/MSSM_two_scale_model.hpp"
#include "models/MSSM/MSSM_two_scale_high_scale_constraint.hpp"
#include "models/MSSM/MSSM_two_scale_susy_scale_constraint.hpp"
#include "models/MSSM/MSSM_two_scale_low_scale_constraint.hpp"
#include "models/MSSM/MSSM_two_scale_convergence_tester.hpp"
#include "models/MSSM/MSSM_two_scale_initial_guesser.hpp"
#include "models/MSSM/MSSM_utilities.hpp"

#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"
#include "coupling_monitor.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <vector>

#include <cstdlib>

#include "consts.hpp"
#include "SM.hpp"
#include "fmssm_oneloop.hpp"
#include "fmssm_lattice.hpp"
#include "lattice_initial_guesser.hpp"
#include "fmssm_lattice_mx_constraint.hpp"
#include "fmssm_lattice_msusy_constraint.hpp"
#include "fmssm_lattice_mz_constraint.hpp"
#include "lattice_convergence_tester.hpp"
#include "lattice_solver.hpp"
#include "logger.hpp"


namespace flexiblesusy {

template<>
class MSSM<Lattice> : public Fmssm<Lattice> {
public:
    MSSM() : Fmssm<Lattice>()
       , problems(&MSSM_info::particle_names_getter,
                  &MSSM_info::parameter_names_getter)
	{}
    void set_input(const MSSM_input_parameters&) {}
    void do_calculate_sm_pole_masses(bool) {}
    const Problems& get_problems() const { return problems; }
    Problems& get_problems() { return problems; }
private:
    Problems problems;
};

// auxiliary class for initializing own members before the base class
// see http://www.boost.org/doc/libs/1_53_0/libs/utility/base_from_member.html
struct MSSM_high_scale_constraint_ {
    MSSM_high_scale_constraint_() : mxc(), mhc(), mgc(), mfc(), tfc() {}
    Fmssm_constraint_on_mx mxc;
    Fmssm_constraint_on_higgs_masses mhc;
    Fmssm_constraint_on_gaugino_masses mgc;
    Fmssm_constraint_on_sfermion_masses mfc;
    Fmssm_constraint_real_trilinear_factors tfc;
};

template<>
class MSSM_high_scale_constraint<Lattice> :
    public MSSM_high_scale_constraint_,
    public CompoundConstraint<Lattice> {
public:
    MSSM_high_scale_constraint() :
	CompoundConstraint<Lattice>::CompoundConstraint
	({&mxc, &mhc, &mgc, &mfc, &tfc})
	{}
    MSSM_high_scale_constraint(const MSSM_input_parameters& i) :
	CompoundConstraint<Lattice>::CompoundConstraint
	({&mxc, &mhc, &mgc, &mfc, &tfc})
    { set_input_parameters(i); }
    void set_input_parameters(const MSSM_input_parameters& i) {
	double m0  = i.m0   *GeV;
	double m12 = i.m12  *GeV;
	double a0  = i.Azero*GeV;
	mgc.M1 = mgc.M2 = mgc.M3 = m12;
	mfc.m2Q << sqr(m0), 0, 0,
	    	   0, sqr(m0), 0,
	    	   0, 0, sqr(m0);
	mfc.m2U = mfc.m2D = mfc.m2L = mfc.m2E = mfc.m2Q;
	mhc.m2Hu_fin = mhc.m2Hd_fin = sqr(m0);
	// to be set by initial guesser
	// mhc.m2Hu_ini = m2Hu_ini;
	// mhc.m2Hd_ini = m2Hd_ini;
	tfc.Au << a0, a0, a0,
	    	  a0, a0, a0,
	    	  a0, a0, a0;
	tfc.Ad = tfc.Ae = tfc.Au;
    }
    void reset() {}		///< reset to initial state
};

#if 0
// Au, Ad, Ae complex -> trilinear equations become trivial
// resulting in a singular matrix when a0 = 0
// see subroutine Fmssm_trilinear_factors in fmssm_lattice_constraints.f
template<>
class MSSM_high_scale_constraint<Lattice> : public Fmssm_mx_constraint {
public:
    MSSM_high_scale_constraint(const MSSM_input_parameters& i) {
	double m0  = i.m0   *GeV;
	double m12 = i.m12  *GeV;
	double a0  = i.Azero*GeV;
	mgc.M1 = mgc.M2 = mgc.M3 = m12;
	mfc.m2Q << sqr(m0), 0, 0,
	    	   0, sqr(m0), 0,
	    	   0, 0, sqr(m0);
	mfc.m2U = mfc.m2D = mfc.m2L = mfc.m2E = mfc.m2Q;
	mhc.m2Hu_fin = mhc.m2Hd_fin = sqr(m0);
	// to be set by initial guesser
	// mhc.m2Hu_ini = m2Hu_ini;
	// mhc.m2Hd_ini = m2Hd_ini;
	tfc.Au << a0, a0, a0,
	    	  a0, a0, a0,
	    	  a0, a0, a0;
	tfc.Ad = tfc.Ae = tfc.Au;
    }
};
#endif

template<>
class MSSM_susy_scale_constraint<Lattice> : public Fmssm_msusy_constraint {
public:
    MSSM_susy_scale_constraint() :
	Fmssm_msusy_constraint(10)
	{}
    MSSM_susy_scale_constraint(const MSSM_input_parameters& input) :
	Fmssm_msusy_constraint(input.TanBeta)
	{}
    void set_input_parameters(const MSSM_input_parameters& i) {
	Real beta = atan(i.TanBeta);
	Real vu = vv * sin(beta);
	Real vd = vv * cos(beta);
	msc.vu = ewsb.vu = vu;
	msc.vd = ewsb.vd = vd;
    }
    void reset() {}		///< reset to initial state
};

template<>
class MSSM_low_scale_constraint<Lattice> : public Fmssm_mz_constraint {
public:
    MSSM_low_scale_constraint() :
	Fmssm_mz_constraint(10)
	{}
    MSSM_low_scale_constraint(const MSSM_input_parameters& input) :
	Fmssm_mz_constraint(input.TanBeta)
	{}
    void set_input_parameters(const MSSM_input_parameters& i) {
	Real beta = atan(i.TanBeta);
	Real vu = vv * sin(beta);
	Real vd = vv * cos(beta);

	CM33 VCKM = standard_VCKM(60*deg);

	CM33 MUMW;
	MUMW << muMW, 0,    0,
	        0,    mcMW, 0,
	        0,    0,    mtMW;
	CM33 MDMW;
	MDMW << mdMW, 0,    0,
	        0,    msMW, 0,
	        0,    0,    mbMW;
	CM33 MEMW;
	MEMW << me,   0,    0,
	        0,    mmu,  0,
	        0,    0,    mtau;

	ycs.Yu = VCKM.transpose() * MUMW / vu;
	ycs.Yd = MDMW / vd;
	ycs.Ye = MEMW / vd;
    }
    void reset() {}		///< reset to initial state
};

template<>
class MSSM_initial_guesser<Lattice> : public Initial_guesser<Lattice> {
public:
    MSSM_initial_guesser<Lattice>
    (void *,
     const MSSM_input_parameters& i,
     const QedQcd&,
     MSSM_low_scale_constraint<Lattice>& mzc_,
     MSSM_susy_scale_constraint<Lattice>& msc_,
     MSSM_high_scale_constraint<Lattice>& mxc_) :
	mu(i.SignMu*i.m0*GeV),	// sign wish enters here
	b(sqr(i.m0*GeV)),	// something of the right order
	mzc(mzc_),
	msc(msc_),
	mxc(mxc_)
    {}

    void operator()() {
	size_t fmssm_height = f->efts[0].height;

	Real gY = mzc.gcs.g1 * sqrt(3/5.0);
	Real g2 = mzc.gcs.g2;
	Comp Yt = mzc.ycs.Yu(2,2); // neglecting CKM rotation
	Comp At = mxc.tfc.Au(2,2);
	Real vu = msc.ewsb.vu;
	Real vd = msc.ewsb.vd;
	// invent Higgs soft masses meeting EWSB conditions for a=0
	mxc.mhc.m2Hu_ini = b*vd/vu - sqr(mu) +
	    (sqr(g2)+sqr(gY))*(Pow4(vd)-Pow4(vu))/(4*(sqr(vd)+sqr(vu)));
	mxc.mhc.m2Hd_ini = b*vu/vd - sqr(mu) +
	    (sqr(g2)+sqr(gY))*(Pow4(vu)-Pow4(vd))/(4*(sqr(vd)+sqr(vu)));

	Real m2Q33 = real(mxc.mfc.m2Q(2,2));
	Real m2U33 = real(mxc.mfc.m2U(2,2));

	Real M2stLL =
	    m2Q33+sqr(vu)*norm(Yt)+(3*sqr(g2)-sqr(gY))*(sqr(vd)-sqr(vu))/12;
	Real M2stRR = m2U33+sqr(vu)*norm(Yt)+sqr(gY)*(sqr(vd)-sqr(vu))/3;
	Comp M2stLR = vu*conj(At*Yt) - vd*conj(Yt)*mu;
	Real detM2st = M2stLL*M2stRR - norm(M2stLR);
	Real MS_guess = pow(detM2st, 1/4.0);

	Real mxGuess = MX1L();
	Real t0 = log(mZ/f->scl0);
	Real t1 = log(MS_guess/f->scl0);
	Real t2 = log(mxGuess/f->scl0);

	const Fmssm<Lattice>& fmssm =
	    *static_cast<const Fmssm<Lattice>*>(f->efts[0].w);

	// when a = 0, no Wilson coefficient runs except for g1, g2, g3

	fmssm(0).m2Hu() = mxc.mhc.m2Hu_ini;
	fmssm(0).m2Hd() = mxc.mhc.m2Hd_ini;
	fmssm(0).mu()   = mu;
	fmssm(0).b()    = b;
	fmssm(0).M1()   = mxc.mgc.M1;
	fmssm(0).M2()   = mxc.mgc.M2;
	fmssm(0).M3()   = mxc.mgc.M3;
	for (size_t i = 0; i < 3; i++)
	    for (size_t j = 0; j < 3; j++) {
		fmssm(0).Yu(i,j)  = mzc.ycs.Yu(i,j);
		fmssm(0).Yd(i,j)  = mzc.ycs.Yd(i,j);
		fmssm(0).Ye(i,j)  = mzc.ycs.Ye(i,j);
		fmssm(0).m2Q(i,j) = mxc.mfc.m2Q(i,j);
		fmssm(0).m2U(i,j) = mxc.mfc.m2U(i,j);
		fmssm(0).m2D(i,j) = mxc.mfc.m2D(i,j);
		fmssm(0).m2L(i,j) = mxc.mfc.m2L(i,j);
		fmssm(0).m2E(i,j) = mxc.mfc.m2E(i,j);
		fmssm(0).TAu(i,j) = mxc.tfc.Au(i,j)*mzc.ycs.Yu(i,j);
		fmssm(0).TAd(i,j) = mxc.tfc.Ad(i,j)*mzc.ycs.Yd(i,j);
		fmssm(0).TAe(i,j) = mxc.tfc.Ae(i,j)*mzc.ycs.Ye(i,j);
	    }
	for (size_t m = 1; m < fmssm_height; m++)
	    for (size_t i = 0; i < fmssm.width; i++) f->y(0,m,i) = f->y(0,0,i);

	fmssm(0).t()  = t0;
	fmssm(0).g1() = g1L1(exp(t0)*f->scl0);
	fmssm(0).g2() = g1L2(exp(t0)*f->scl0);
	fmssm(0).g3() = g1L3(exp(t0)*f->scl0);
	fmssm(1).t()  = t1;
	fmssm(1).g1() = g1L1(exp(t1)*f->scl0);
	fmssm(1).g2() = g1L2(exp(t1)*f->scl0);
	fmssm(1).g3() = g1L3(exp(t1)*f->scl0);
	fmssm(2).t()  = t2;
	fmssm(2).g1() = g1L1(exp(t2)*f->scl0);
	fmssm(2).g2() = g1L2(exp(t2)*f->scl0);
	fmssm(2).g3() = g1L3(exp(t2)*f->scl0);

	// TODO: set units
    }

private:
    Real mu;
    Real b;
    MSSM_low_scale_constraint<Lattice>& mzc;
    MSSM_susy_scale_constraint<Lattice>& msc;
    MSSM_high_scale_constraint<Lattice>& mxc;
};

template<>
class MSSM_convergence_tester<Lattice> : public Convergence_tester<Lattice> {
public:
   MSSM_convergence_tester(MSSM<Lattice>*, double accuracy_goal) {}
   virtual ~MSSM_convergence_tester();
   void set_max_iterations(int) {}; ///< set maximum number of iterations
};

template<class T>
class MSSM_runner {
public:
   MSSM_runner()
      : solver(), model()
      , high_scale_constraint()
      , susy_scale_constraint()
      , low_scale_constraint()
      , high_scale(0.), susy_scale(0.), low_scale(0.)
      , precision_goal(1.0e-5)
      , max_iterations(0)
      , calculate_sm_masses(false) {}
   ~MSSM_runner() {}

   double get_high_scale() const { return high_scale; }
   double get_susy_scale() const { return susy_scale; }
   double get_low_scale()  const { return low_scale;  }
   const MSSM<T>& get_model() const { return model; }
   const Problems& get_problems() const { return model.get_problems(); }
   void set_precision_goal(double precision_goal_) { precision_goal = precision_goal_; }
   void set_max_iterations(int n) { max_iterations = n; }
   void set_calculate_sm_masses(bool flag) { calculate_sm_masses = flag; }

   void run(const QedQcd& qedqcd, const MSSM_input_parameters& input);
   void write_running_couplings(const std::string& filename = "MSSM_rge_running.dat") const;
   void write_spectrum(const std::string& filename = "MSSM_spectrum.dat") const;

private:
   RGFlow<T> solver;
   MSSM<T> model;
   MSSM_high_scale_constraint<T> high_scale_constraint;
   MSSM_susy_scale_constraint<T> susy_scale_constraint;
   MSSM_low_scale_constraint<T>  low_scale_constraint;
   double high_scale, susy_scale, low_scale;
   double precision_goal; ///< precision goal
   int max_iterations; ///< maximum number of iterations
   bool calculate_sm_masses; ///< calculate SM pole masses
};

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param qedqcd Standard Model input parameters
 * @param input model input parameters
 */
template<class T>
void MSSM_runner<T>::run(const QedQcd& qedqcd,
			 const MSSM_input_parameters& input)
{
   high_scale_constraint.reset();
   susy_scale_constraint.reset();
   low_scale_constraint.reset();
   high_scale_constraint.set_input_parameters(input);
   susy_scale_constraint.set_input_parameters(input);
   low_scale_constraint .set_input_parameters(input);

   std::vector<Constraint<T>*> upward_constraints;
   upward_constraints.push_back(&low_scale_constraint);
   upward_constraints.push_back(&high_scale_constraint);

   std::vector<Constraint<T>*> downward_constraints;
   downward_constraints.push_back(&high_scale_constraint);
   downward_constraints.push_back(&susy_scale_constraint);
   downward_constraints.push_back(&low_scale_constraint);

   model.set_input(input);
   model.do_calculate_sm_pole_masses(calculate_sm_masses);

   MSSM_convergence_tester<T> convergence_tester(&model, precision_goal);
   if (max_iterations > 0)
      convergence_tester.set_max_iterations(max_iterations);

   MSSM_initial_guesser<T> initial_guesser(&model, input, qedqcd,
					 low_scale_constraint,
					 susy_scale_constraint,
					 high_scale_constraint);
   Two_scale_increasing_precision precision(10.0, precision_goal);

   solver.reset();
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_model(&model, upward_constraints, downward_constraints);

   high_scale = susy_scale = low_scale = 0.;

   try {
      solver.solve();
      high_scale = high_scale_constraint.get_scale();
      susy_scale = susy_scale_constraint.get_scale();
      low_scale  = low_scale_constraint.get_scale();

      if (model.run_to(susy_scale))
         throw NonPerturbativeRunningError(susy_scale);

      model.calculate_spectrum();

      if (model.run_to(low_scale))
         throw NonPerturbativeRunningError(low_scale);
   } catch (const NoConvergenceError& error) {
      model.get_problems().flag_thrown("no convergence");
      ERROR(error.what());
   } catch (const NonPerturbativeRunningError& error) {
      model.get_problems().flag_no_perturbative();
      ERROR(error.what());
   } catch (const Error& error) {
      model.get_problems().flag_thrown();
      ERROR(error.what());
   } catch (const std::string& str) {
      model.get_problems().flag_thrown();
      ERROR(str);
   } catch (const char* str) {
      model.get_problems().flag_thrown();
      ERROR(str);
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
template<class T>
void MSSM_runner<T>::write_running_couplings(const std::string& filename) const
{
#if 0
   // FIXME: impossible for lattice method to simulate since tmp_model
   // by itself would be stateless unless linked to an RGFlow
   // POSSIBLE SOLUTION: duplicate entire RGFlow
   MSSM<T> tmp_model(model);
   const int error = tmp_model.run_to(low_scale);
   if (error) {
      ERROR("MSSM_runner::write_running_couplings: run to scale "
            << low_scale << " failed");
      return;
   }

   MSSM_parameter_getter parameter_getter;
   Coupling_monitor<MSSM<T>, MSSM_parameter_getter>
      coupling_monitor(tmp_model, parameter_getter);

   coupling_monitor.run(low_scale, high_scale, 100, true);
   coupling_monitor.write_to_file(filename);
#endif
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template<class T>
void MSSM_runner<T>::write_spectrum(const std::string& filename) const
{
#if 0
   MSSM_spectrum_plotter plotter;
   plotter.extract_spectrum(model);
   plotter.write_to_file(filename);
#endif
}

} // namespace flexiblesusy

template<class T>
int main_(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef T algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_exit())
      return options.status();

   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   MSSM_slha_io slha_io;
   Spectrum_generator_settings spectrum_generator_settings;
   QedQcd qedqcd;
   MSSM_input_parameters input;

   if (!slha_input_file.empty()) {
      slha_io.read_from_file(slha_input_file);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
      slha_io.fill(spectrum_generator_settings);
   }
   qedqcd.toMz(); // run SM fermion masses to MZ

   MSSM_runner<algorithm_type> runner;
   runner.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   runner.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   runner.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);

   runner.run(qedqcd, input);

   if (runner.get_problems().have_problem()) {
      runner.get_problems().print();
   } else {
      // runner.write_running_couplings();
      // runner.write_spectrum();
   }

#if 0
   slha_io.set_spinfo(runner.get_problems());
   slha_io.set_spectrum(runner.get_model());
   slha_io.write_to_file(slha_output_file);
#endif
   return 0;
}

using namespace flexiblesusy;

int main(int argc, char *argv[])
{
    if (argc > 1)
	switch (argv[1][0]) {
	case 't':
	    return main_<Two_scale>(argc-1, (const char **)(argv+1));
	case 'l':
	    return main_<Lattice>  (argc-1, (const char **)(argv+1));
	}
    return 1;
}
