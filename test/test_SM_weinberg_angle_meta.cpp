
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_weinberg_angle_meta

#include <boost/test/unit_test.hpp>

#include "SM_mass_eigenstates.hpp"
#include "test_SM.hpp"

#define private public

#include "weinberg_angle.hpp"
#include "SM_weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

BOOST_AUTO_TEST_CASE( test_rho_2 )
{
   double r;

   r = 0.1;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                                 SM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 1.8;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                                 SM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 1.9;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                                 SM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 2.0;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                                 SM_weinberg_angle::rho_2(r), 1.0e-10);

   r = 2.1;
   BOOST_CHECK_CLOSE_FRACTION(Weinberg_angle::rho_2(r),
                                 SM_weinberg_angle::rho_2(r), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   standard_model::Standard_model sm;
   SM_mass_eigenstates SM;

   setup_SM_const(sm, input);
   setup_SM_const(SM, input);

   sm.set_thresholds(2);
   sm.calculate_DRbar_masses();
   SM.set_thresholds(2);
   SM.calculate_DRbar_masses();

   const double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   SM_weinberg_angle::Sm_parameters SM_parameters;
   SM_parameters.fermi_constant = Electroweak_constants::gfermi;
   SM_parameters.mw_pole = Electroweak_constants::MW;
   SM_parameters.mz_pole = Electroweak_constants::MZ;
   SM_parameters.mt_pole = 165.0;
   SM_parameters.alpha_s = 0.1176;

   Weinberg_angle wein_sm(&sm, sm_parameters);
   SM_weinberg_angle wein_SM(&SM, SM_parameters);

   const double delta_vb_sm =
      wein_sm.calculate_delta_vb(outrho, outsin);
   const double delta_vb_SM = wein_SM.calculate_delta_vb(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_vb_sm, delta_vb_SM, 1.0e-16);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   standard_model::Standard_model sm;
   SM_mass_eigenstates SM;

   setup_SM_const(sm, input);
   setup_SM_const(SM, input);

   sm.set_thresholds(2);
   sm.calculate_DRbar_masses();
   SM.set_thresholds(2);
   SM.calculate_DRbar_masses();

   const double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   SM_weinberg_angle::Sm_parameters SM_parameters;
   SM_parameters.fermi_constant = Electroweak_constants::gfermi;
   SM_parameters.mw_pole = Electroweak_constants::MW;
   SM_parameters.mz_pole = Electroweak_constants::MZ;
   SM_parameters.mt_pole = 165.0;
   SM_parameters.alpha_s = 0.1176;

   Weinberg_angle wein_sm(&sm, sm_parameters);
   SM_weinberg_angle wein_SM(&SM, SM_parameters);

   // initialize self-energies
   wein_sm.pizzt_MZ = wein_sm.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein_sm.piwwt_MW = wein_sm.calculate_self_energy_VWp(Electroweak_constants::MW);
   wein_sm.piwwt_0  = wein_sm.calculate_self_energy_VWp(0.);
   wein_SM.pizzt_MZ = wein_SM.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein_SM.piwwt_MW = wein_SM.calculate_self_energy_VWp(Electroweak_constants::MW);
   wein_SM.piwwt_0  = wein_SM.calculate_self_energy_VWp(0.);

   const double delta_r_sm = wein_sm.calculate_delta_r_hat(outrho, outsin);
   const double delta_r_SM = wein_SM.calculate_delta_r_hat(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_r_sm, delta_r_SM, 8e-16);
}

BOOST_AUTO_TEST_CASE( test_sin_theta )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   standard_model::Standard_model sm;
   SM_mass_eigenstates SM;

   setup_SM_const(sm, input);
   setup_SM_const(SM, input);

   sm.set_thresholds(2);
   sm.calculate_DRbar_masses();
   SM.set_thresholds(2);
   SM.calculate_DRbar_masses();

   static constexpr double tol = 1.0e-10;
   static constexpr int maxTries = 20;
   static constexpr double rho_start = 1.0;
   static constexpr double sin_start = 0.48;

   Weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   Weinberg_angle wein_sm(&sm, sm_parameters);
   wein_sm.set_number_of_iterations(maxTries);
   wein_sm.set_precision_goal(tol);
   double sin_theta_sm;
   BOOST_REQUIRE_NO_THROW(sin_theta_sm = wein_sm.calculate(sin_start));

   SM_weinberg_angle::Sm_parameters SM_parameters;
   SM_parameters.fermi_constant = Electroweak_constants::gfermi;
   SM_parameters.mw_pole = Electroweak_constants::MW;
   SM_parameters.mz_pole = Electroweak_constants::MZ;
   SM_parameters.mt_pole = 165.0;
   SM_parameters.alpha_s = 0.1176;
   SM_weinberg_angle wein_SM(&SM, SM_parameters);
   wein_SM.set_number_of_iterations(maxTries);
   wein_SM.set_precision_goal(tol);
   double sin_theta_SM;
   BOOST_REQUIRE_NO_THROW(sin_theta_SM = wein_SM.calculate(sin_start).first);

   BOOST_CHECK_CLOSE_FRACTION(sin_theta_sm, sin_theta_SM, 1.0e-16);
}
