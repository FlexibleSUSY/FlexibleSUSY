
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

void setup_data(const SM_mass_eigenstates& model,
                Weinberg_angle::Data& data)
{
   const double mw_pole = Electroweak_constants::MW;
   const double mz_pole = Electroweak_constants::MZ;
   // use mt_drbar instead of mt_pole such that tests work up to high precision
   const double mt_pole = 165.0;
   const double scale   = model.get_scale();
   const double gY      = model.get_g1() * Sqrt(0.6);
   const double g2      = model.get_g2();
   const double e_drbar = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));

   BOOST_CHECK_EQUAL(scale, mz_pole);

   data.fermi_contant  = Electroweak_constants::gfermi;
   data.mw_pole        = mw_pole;
   data.mz_pole        = mz_pole;
   data.mt_pole        = mt_pole;
   data.scale          = scale;
   data.gY             = gY;
   data.g2             = g2;
   data.g3             = model.get_g3();
   data.alpha_em_drbar = Sqr(e_drbar) / (4.0 * Pi);
   data.mh_drbar       = model.get_Mhh();

   const double pizztMZ = Re(model.self_energy_VZ_1loop(mz_pole));
   const double piwwtMW = Re(model.self_energy_VWp_1loop(mw_pole));
   const double piwwt0  = Re(model.self_energy_VWp_1loop(0.));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = model.get_MFu(2);
   se_data.mb_drbar = model.get_MFd(2);
   se_data.gY       = gY;
   se_data.g2       = g2;

   const double pizztMZ_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole, se_data);

   const double piwwtMW_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwtMW, mw_pole, se_data);

   const double piwwt0_corrected =
      Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0., se_data);

   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
}

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
   SM_mass_eigenstates model;
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   setup_SM_const(model, input);
   model.set_thresholds(2);
   model.calculate_DRbar_masses();

   const double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   const double delta_vb_1 =
      Weinberg_angle::calculate_delta_vb(outrho, outsin, data, false);

   SM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   SM_weinberg_angle wein(&model, sm_parameters);
   const double delta_vb_2 = wein.calculate_delta_vb(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_vb_1, delta_vb_2, 1.0e-6);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   SM_mass_eigenstates model;
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   setup_SM_const(model, input);
   model.set_thresholds(2);
   model.calculate_DRbar_masses();

   const double outrho = 1.0, outsin = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   const double delta_r_1 =
      Weinberg_angle::calculate_delta_r(outrho, outsin, data, false);

   SM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   SM_weinberg_angle wein(&model, sm_parameters);
   // initialize self-energies
   wein.pizzt_MZ = wein.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein.piwwt_MW = wein.calculate_self_energy_VWp(Electroweak_constants::MW);
   wein.piwwt_0  = wein.calculate_self_energy_VWp(0.);
   const double delta_r_2 = wein.calculate_delta_r_hat(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_r_1, delta_r_2, 2.0e-6);
}

BOOST_AUTO_TEST_CASE( test_sin_theta )
{
   SM_mass_eigenstates model;
   SM_input_parameters input;
   input.LambdaIN = 0.25;

   setup_SM_const(model, input);
   model.set_thresholds(2);
   model.calculate_DRbar_masses();

   const double tol = 1.0e-10;
   const int maxTries = 20;
   const double rho_start = 1.0, sin_start = 0.48;

   Weinberg_angle::Data data;
   setup_data(model, data);
   Weinberg_angle wein;
   wein.disable_susy_contributions();
   wein.set_data(data);
   wein.set_number_of_iterations(maxTries);
   wein.set_precision_goal(tol);
   const int error = wein.calculate(rho_start, sin_start);
   const double sin_theta_1 = wein.get_sin_theta();
   BOOST_REQUIRE(error == 0);

   SM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = 165.0;
   sm_parameters.alpha_s = 0.1176;
   SM_weinberg_angle pwein(&model, sm_parameters);
   pwein.set_number_of_iterations(maxTries);
   pwein.set_precision_goal(tol);
   double sin_theta_2;
   BOOST_REQUIRE_NO_THROW(sin_theta_2 = pwein.calculate(sin_start).first);

   BOOST_CHECK_CLOSE_FRACTION(sin_theta_1, sin_theta_2, 6.0e-9);
}
