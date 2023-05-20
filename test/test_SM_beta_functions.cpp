
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_beta_functions

#include <boost/test/unit_test.hpp>

#include "SM_two_scale_model.hpp"
#include "standard_model.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

namespace {

std::pair<SM_mass_eigenstates, standard_model::Standard_model> make_point()
{
   SM_mass_eigenstates m1;
   standard_model::Standard_model m2;

   const double vev = 246.;
   const double g1 = 0.2;
   const double g2 = 0.3;
   const double g3 = 0.4;

   m1.set_scale(91.);
   m1.set_v(vev);
   m1.set_mu2(vev*vev);
   m1.set_g1(g1);
   m1.set_g2(g2);
   m1.set_g3(g3);
   m1.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m1.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m1.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);
   m1.set_Lambdax(0.2);

   m2.set(m1.get());

   m1.set_loops(5);
   m2.set_loops(5);

   return std::make_pair(m1, m2);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_SM_beta_functions_literature )
{
   SM_input_parameters input;
   SM<Two_scale> m;

   const double gut_norm = Sqrt(0.6); // GUT normalization

   input.LambdaIN = 0.1;

   const double vev = 246.;
   const double g1 = 0.2;
   const double g2 = 0.3;
   const double g3 = 0.4;

   m.set_scale(91.);
   m.set_v(vev);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.set_loops(1);

   const SM_soft_parameters beta(m.calc_beta());

   const double beta_g1 = oneOver16PiSqr * Power(m.get_g1(),3) * Sqr(gut_norm) * (41./6.);
   const double beta_g2 = oneOver16PiSqr * Power(m.get_g2(),3) * (-19./6.);
   const double beta_g3 = oneOver16PiSqr * Power(m.get_g3(),3) * (-7.);

   BOOST_CHECK_CLOSE_FRACTION(beta.get_g1(), beta_g1, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(beta.get_g2(), beta_g2, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(beta.get_g3(), beta_g3, 1.0e-10);
}

void compare_SM_pars(const SM_soft_parameters& m1, const standard_model::Standard_model& m2)
{
   const double eps = 1e-15;

   BOOST_CHECK_CLOSE_FRACTION(m1.get_g1()     , m2.get_g1()     , eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_g2()     , m2.get_g2()     , eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_g3()     , m2.get_g3()     , eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_v()      , m2.get_v()      , eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_mu2()    , m2.get_mu2()    , eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Lambdax(), m2.get_Lambdax(), eps);

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE_FRACTION(m1.get_Yu(i,k), m2.get_Yu(i,k), eps);
         BOOST_CHECK_CLOSE_FRACTION(m1.get_Yd(i,k), m2.get_Yd(i,k), eps);
         BOOST_CHECK_CLOSE_FRACTION(m1.get_Ye(i,k), m2.get_Ye(i,k), eps);
      }
   }

}

BOOST_AUTO_TEST_CASE( test_SM_beta_functions )
{
   const auto m = make_point();

   compare_SM_pars(m.first             , m.second             );
   compare_SM_pars(m.first.calc_beta(0), m.second.calc_beta(0));
   compare_SM_pars(m.first.calc_beta(1), m.second.calc_beta(1));
   compare_SM_pars(m.first.calc_beta(2), m.second.calc_beta(2));
   compare_SM_pars(m.first.calc_beta(3), m.second.calc_beta(3));
   compare_SM_pars(m.first.calc_beta(4), m.second.calc_beta(4));
   compare_SM_pars(m.first.calc_beta(5), m.second.calc_beta(5));
}
