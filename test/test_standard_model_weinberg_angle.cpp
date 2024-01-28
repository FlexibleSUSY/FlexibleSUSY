#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle

#include <boost/test/unit_test.hpp>

#include "conversion.hpp"
#include "SM_two_scale_model.hpp"
#include "standard_model.hpp"
#include "test_SM.hpp"
#include "logger.hpp"

#define private public

#include "weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   static constexpr double outsin = 0.48;

   standard_model::Standard_model sm;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(sm, input);

   Weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   Weinberg_angle wein(&sm, sm_parameters);

   const double fs_delta_vb =
      wein.calculate_delta_vb_sm(outsin);

   BOOST_CHECK_CLOSE_FRACTION(fs_delta_vb, 0.0096639639254316526, 1.0e-16);
}

BOOST_AUTO_TEST_CASE( test_delta_r_hat )
{
   standard_model::Standard_model sm;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(sm, input);

   Weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   Weinberg_angle wein(&sm, sm_parameters);

   static constexpr double rho = 1.0;
   static constexpr double sinw = 0.48;

   wein.set_number_of_loops(0);
   const double fs_delta_r_0l =
      wein.calculate_delta_r_hat(rho, sinw);

   wein.set_number_of_loops(1);
   const double fs_delta_r_1l =
      wein.calculate_delta_r_hat(rho, sinw);

   wein.set_number_of_loops(2);
   const double fs_delta_r_2l =
      wein.calculate_delta_r_hat(rho, sinw);

   BOOST_CHECK_SMALL(fs_delta_r_0l, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_delta_r_1l, 0.0096639639254316526, 1.0e-16);
   BOOST_CHECK_CLOSE_FRACTION(fs_delta_r_2l, 0.010595894357437782, 1.0e-16);
}

BOOST_AUTO_TEST_CASE( test_delta_rho_hat )
{
   standard_model::Standard_model sm;
   SM_input_parameters input;
   input.LambdaIN = 0.1;

   setup_SM_const(sm, input);

   Weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   Weinberg_angle wein(&sm, sm_parameters);

   static constexpr double sinw = 0.48;

   wein.pizzt_MZ = wein.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein.piwwt_MW = wein.calculate_self_energy_VWp(Electroweak_constants::MW);
   wein.piwwt_0  = wein.calculate_self_energy_VWp(0.);

   wein.set_number_of_loops(0);
   const double fs_delta_r_0l =
      wein.calculate_delta_rho_hat(sinw);

   wein.set_number_of_loops(1);
   const double fs_delta_r_1l =
      wein.calculate_delta_rho_hat(sinw);

   wein.set_number_of_loops(2);
   const double fs_delta_r_2l =
      wein.calculate_delta_rho_hat(sinw);

   BOOST_CHECK_SMALL(fs_delta_r_0l, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_delta_r_1l, 0.035080116163659691, 1.0e-16);
   BOOST_CHECK_CLOSE_FRACTION(fs_delta_r_2l, 0.034020320922642386, 1.0e-16);
}

