
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "loop_libraries/loop_library.hpp"
#include "test_SM.hpp"
#include "wrappers.hpp"
#include "SM_two_scale_model.hpp"
#include "standard_model.hpp"
#include "config.h"

#define SARAH_VERSION_AT_LEAST(x,y,z) (                                 \
      (SARAH_MAJOR == 0 && SARAH_MINOR == 0 && SARAH_PATCH == 0) ||     \
      (SARAH_MAJOR > x || (SARAH_MAJOR >= x &&                          \
                           (SARAH_MINOR > y || (SARAH_MINOR >= y &&     \
                                                SARAH_PATCH >= z))))    \
      )

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_one_loop_Higgs_masses )
{
   Loop_library::set(0);
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   // set mu2 which is not in agreement with EWSB
   const double mu2 = 100.;
   m.set_mu2(mu2);

   m.calculate_DRbar_masses();

   // check that mu2 was not changed
   BOOST_CHECK_CLOSE(m.get_mu2(), mu2, 1.0e-10);

   // check Higgs tree-level mass
   const double Mhh(m.get_Mhh());
   const double lambda = m.get_Lambdax();
   const double v = m.get_v();
   const double hh_tree = Sqrt(lambda * Sqr(v));

   BOOST_CHECK_CLOSE(Mhh, hh_tree, 1.0e-12);

   m.set_pole_mass_loop_order(1);
   m.do_calculate_sm_pole_masses(true);
   m.solve_ewsb_one_loop();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   // check that mu2 was changed
   BOOST_CHECK(m.get_mu2() != mu2);

   // check Higgs pole mass
   const double Mhh_1l(m.get_physical().Mhh);

   // Note: in SARAH v4.13.0 the definition of mu2 was changed by an
   // overall sign
#if SARAH_VERSION_AT_LEAST(4,13,0)
   const double hh_1l = Sqrt(m.get_mu2() + 1.5*lambda*Sqr(v) - Re(m.self_energy_hh_1loop(Mhh_1l)));
#else
   const double hh_1l = Sqrt(-m.get_mu2() + 1.5*lambda*Sqr(v) - Re(m.self_energy_hh_1loop(Mhh_1l)));
#endif

   BOOST_CHECK_CLOSE(Mhh_1l, hh_1l, 2.0e-4);

   // check that tree-level Higgs mass has not changed
   BOOST_CHECK_CLOSE(m.get_Mhh(), hh_tree, 1.0e-12);
}

BOOST_AUTO_TEST_CASE( test_SM_heavy_top_self_energy )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.calculate_DRbar_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double p = 100.;

   const Eigen::Array<double,3,1> MFu(m.get_MFu());
   const double g3 = m.get_g3();
   const double scale = m.get_scale();

   Eigen::Matrix<std::complex<double>,3,3> se_t;
   Eigen::Matrix<std::complex<double>,3,3> se_t_no_gluon;

   // top self-energies with and without gluon
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         se_t(i,k) = m.self_energy_Fu_1loop_1(p,i,k);
         se_t_no_gluon(i,k) = m.self_energy_Fu_1loop_1_heavy_rotated(p,i,k);
      }
   }

   // adding gluon contrbution
   Eigen::Matrix<std::complex<double>,3,3> se_t_check(se_t_no_gluon);

   for (int i = 0; i < 3; i++) {
      const double gluon_contrib =
         -5.333333333333333 *
         (-0.5 + Loop_library::get().B0(Sqr(p),Sqr(MFu(i)),0,Sqr(scale)).real())
            * Sqr(g3) * MFu(i);

      se_t_check(i,i) += gluon_contrib * oneOver16PiSqr;
   }

   const Eigen::Matrix<std::complex<double>,3,3> Uu(m.get_Uu());
   const Eigen::Matrix<std::complex<double>,3,3> Vu(m.get_Vu());

   se_t_check = Vu.transpose() * se_t_check * Uu;

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE(Re(se_t(i,k)), Re(se_t_check(i,k)), 1.0e-10);
         BOOST_CHECK_CLOSE(Im(se_t(i,k)), Im(se_t_check(i,k)), 1.0e-10);
      }
   }

}

#define CHECK_CLOSE_1(mass,eps)                                         \
   do {                                                                 \
      BOOST_CHECK_CLOSE(m.get_physical().mass, sm.get_physical().mass, (eps)); \
      BOOST_TEST_MESSAGE(#mass " = " << m.get_physical().mass << " = " << sm.get_physical().mass); \
   } while (false);

#define CHECK_CLOSE_N(mass,eps,N)                                       \
   do {                                                                 \
      for (int i = 0; i < N; i++) {                                     \
         BOOST_CHECK_CLOSE(m.get_physical().mass(i), sm.get_physical().mass(i), (eps)); \
         BOOST_TEST_MESSAGE(#mass "(" << i << ") = " << m.get_physical().mass(i) << " = " << sm.get_physical().mass(i)); \
      }                                                                 \
   } while (false);

BOOST_AUTO_TEST_CASE( test_SM_one_loop_masses )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   standard_model::Standard_model sm;
   setup_SM_const(m, input);
   setup_SM_const(sm, input);

   m.set_thresholds(1);
   sm.set_thresholds(1);

   m.do_calculate_sm_pole_masses(true);

   m.calculate_DRbar_masses();
   sm.calculate_DRbar_masses();

   m.solve_ewsb_one_loop();
   sm.solve_ewsb_one_loop();

   m.calculate_pole_masses();
   sm.calculate_pole_masses();

   const double eps = 1e-13;

   CHECK_CLOSE_1(MVP, eps);
   CHECK_CLOSE_1(MVG, eps);
   CHECK_CLOSE_1(MVZ, eps);
   CHECK_CLOSE_1(MVWp, eps);
   CHECK_CLOSE_1(MAh, eps);
   CHECK_CLOSE_1(MHp, eps);
   CHECK_CLOSE_1(Mhh, eps);
   CHECK_CLOSE_N(MFv, eps, 3);
   CHECK_CLOSE_N(MFu, eps, 3);
   CHECK_CLOSE_N(MFd, eps, 3);
   CHECK_CLOSE_N(MFe, eps, 3);

   const double MW = m.get_physical().MVWp;
   const double MZ = m.get_physical().MVZ;

   // specific pole mass functions for gauge bosons
   BOOST_CHECK_CLOSE(m.calculate_MVWp_pole(MW), sm.calculate_MVWp_pole(MW), eps);
   BOOST_CHECK_CLOSE(m.calculate_MVWp_pole(MW), MW, 0.02);

   BOOST_CHECK_CLOSE(m.calculate_MVZ_pole(MZ) , sm.calculate_MVZ_pole(MZ) , eps);
   BOOST_CHECK_CLOSE(m.calculate_MVZ_pole(MZ) , MZ, 0.02);

   // specific DR-bar mass functions
   BOOST_CHECK_CLOSE(m.calculate_MVP_DRbar(0), sm.calculate_MVP_DRbar(0), eps);
   BOOST_CHECK_CLOSE(m.calculate_MVP_DRbar(0), 0., eps);

   BOOST_CHECK_CLOSE(m.calculate_MVZ_DRbar(MZ), sm.calculate_MVZ_DRbar(MZ), eps);
   BOOST_CHECK_CLOSE(m.calculate_MVZ_DRbar(MZ), m.get_MVZ(), 0.01);

   BOOST_CHECK_CLOSE(m.calculate_MVWp_DRbar(MW), sm.calculate_MVWp_DRbar(MW), eps);
   BOOST_CHECK_CLOSE(m.calculate_MVWp_DRbar(MW), m.get_MVWp(), 0.02);

   for (int i = 0; i < 3; i++) {
      const double eps = 1e-5, eps2 = 0.0002;

      {
         const double p = 0.;
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFv_DRbar(p,i), sm.calculate_MFv_DRbar(p,i), eps);
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFv_DRbar(p,i), m.get_MFv(i), eps2);
      }
      {
         const double p = m.get_MFe(i);
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFe_DRbar(p,i), sm.calculate_MFe_DRbar(p,i), eps);
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFe_DRbar(p,i), m.get_MFe(i), eps2);
      }
      {
         const double p = m.get_physical().MFu(i);
         const double eps2 = i < 2 ? 0.4 : 0.006;
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFu_DRbar(p,i), sm.calculate_MFu_DRbar(p,i), eps);
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFu_DRbar(p,i), m.get_MFu(i), eps2);
      }
      {
         const double p = m.get_MFd(i);
         const double eps2 = i < 2 ? 0.4 : 0.005;
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFd_DRbar(p,i), sm.calculate_MFd_DRbar(p,i), eps);
         BOOST_CHECK_CLOSE_FRACTION(m.calculate_MFd_DRbar(p,i), m.get_MFd(i), eps2);
      }
   }
}
