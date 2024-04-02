
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMNoFV_weinberg_angle_meta

#include <boost/test/unit_test.hpp>

#include "CMSSMNoFV_mass_eigenstates.hpp"
#include "CMSSM_mass_eigenstates.hpp"
#include "test_CMSSMNoFV.hpp"

#define private public

#include "CMSSMNoFV_weinberg_angle.hpp"
#include "CMSSM_weinberg_angle.hpp"

using namespace flexiblesusy;

void setup_CMSSM_models_non_3rd_gen(CMSSMNoFV_mass_eigenstates& m1,
                                    CMSSM_mass_eigenstates& m2,
                                    const CMSSMNoFV_input_parameters& input)
{
   setup_CMSSM_models(m1, m2, input);
   m1.set_thresholds(2);
   m2.set_thresholds(2);

   const double ymu = 0.001;

   Eigen::Matrix<double,3,3> Ye1 = m1.get_Ye();
   Ye1(1,1) = ymu;
   m1.set_Ye(Ye1);
   Eigen::Matrix<double,3,3> Ye2 = m2.get_Ye();
   Ye2(1,1) = ymu;
   m2.set_Ye(Ye2);
}

BOOST_AUTO_TEST_CASE( test_delta_vb )
{
   CMSSMNoFV_mass_eigenstates m1;
   CMSSM_mass_eigenstates m2;
   CMSSMNoFV_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_models_non_3rd_gen(m1, m2, input);

   const double outrho = 1.0, outsin = 0.48;

   CMSSMNoFV_weinberg_angle::Parameters parameters_1;
   parameters_1.fermi_constant = Electroweak_constants::gfermi;
   parameters_1.mw_pole = Electroweak_constants::MW;
   parameters_1.mz_pole = Electroweak_constants::MZ;
   parameters_1.mt_pole = 165.0;
   parameters_1.alpha_s = 0.1176;
   CMSSMNoFV_weinberg_angle wein1(&m1, parameters_1);
   const double delta_vb_1 = wein1.calculate_delta_vb(outrho, outsin);

   CMSSM_weinberg_angle::Parameters parameters_2;
   parameters_2.fermi_constant = Electroweak_constants::gfermi;
   parameters_2.mw_pole = Electroweak_constants::MW;
   parameters_2.mz_pole = Electroweak_constants::MZ;
   parameters_2.mt_pole = 165.0;
   parameters_2.alpha_s = 0.1176;
   CMSSM_weinberg_angle wein2(&m2, parameters_2);
   const double delta_vb_2 = wein2.calculate_delta_vb(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_vb_1, delta_vb_2, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_r )
{
   CMSSMNoFV_mass_eigenstates m1;
   CMSSM_mass_eigenstates m2;
   CMSSMNoFV_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_models_non_3rd_gen(m1, m2, input);

   const double outrho = 1.0, outsin = 0.48;

   CMSSMNoFV_weinberg_angle::Parameters parameters_1;
   parameters_1.fermi_constant = Electroweak_constants::gfermi;
   parameters_1.mw_pole = Electroweak_constants::MW;
   parameters_1.mz_pole = Electroweak_constants::MZ;
   parameters_1.mt_pole = 165.0;
   parameters_1.alpha_s = 0.1176;
   CMSSMNoFV_weinberg_angle wein1(&m1, parameters_1);
   // initialize self-energies
   wein1.pizzt_MZ = wein1.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein1.piwwt_MW = wein1.calculate_self_energy_VWm(Electroweak_constants::MW);
   wein1.piwwt_0  = wein1.calculate_self_energy_VWm(0.);
   const double delta_r_1 = wein1.calculate_delta_r_hat(outrho, outsin);

   CMSSM_weinberg_angle::Parameters parameters_2;
   parameters_2.fermi_constant = Electroweak_constants::gfermi;
   parameters_2.mw_pole = Electroweak_constants::MW;
   parameters_2.mz_pole = Electroweak_constants::MZ;
   parameters_2.mt_pole = 165.0;
   parameters_2.alpha_s = 0.1176;
   CMSSM_weinberg_angle wein2(&m2, parameters_2);
   // initialize self-energies
   wein2.pizzt_MZ = wein2.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein2.piwwt_MW = wein2.calculate_self_energy_VWm(Electroweak_constants::MW);
   wein2.piwwt_0  = wein2.calculate_self_energy_VWm(0.);
   const double delta_r_2 = wein2.calculate_delta_r_hat(outrho, outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_r_1, delta_r_2, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_delta_rho )
{
   CMSSMNoFV_mass_eigenstates m1;
   CMSSM_mass_eigenstates m2;
   CMSSMNoFV_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_models_non_3rd_gen(m1, m2, input);

   const double outsin = 0.48;

   CMSSMNoFV_weinberg_angle::Parameters parameters_1;
   parameters_1.fermi_constant = Electroweak_constants::gfermi;
   parameters_1.mw_pole = Electroweak_constants::MW;
   parameters_1.mz_pole = Electroweak_constants::MZ;
   parameters_1.mt_pole = 165.0;
   parameters_1.alpha_s = 0.1176;
   CMSSMNoFV_weinberg_angle wein1(&m1, parameters_1);
   // initialize self-energies
   wein1.pizzt_MZ = wein1.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein1.piwwt_MW = wein1.calculate_self_energy_VWm(Electroweak_constants::MW);
   wein1.piwwt_0  = wein1.calculate_self_energy_VWm(0.);
   const double delta_rho_1 = wein1.calculate_delta_rho_hat(outsin);

   CMSSM_weinberg_angle::Parameters parameters_2;
   parameters_2.fermi_constant = Electroweak_constants::gfermi;
   parameters_2.mw_pole = Electroweak_constants::MW;
   parameters_2.mz_pole = Electroweak_constants::MZ;
   parameters_2.mt_pole = 165.0;
   parameters_2.alpha_s = 0.1176;
   CMSSM_weinberg_angle wein2(&m2, parameters_2);
   // initialize self-energies
   wein2.pizzt_MZ = wein2.calculate_self_energy_VZ(Electroweak_constants::MZ);
   wein2.piwwt_MW = wein2.calculate_self_energy_VWm(Electroweak_constants::MW);
   wein2.piwwt_0  = wein2.calculate_self_energy_VWm(0.);
   const double delta_rho_2 = wein2.calculate_delta_rho_hat(outsin);

   BOOST_CHECK_CLOSE_FRACTION(delta_rho_1, delta_rho_2, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_sin_theta )
{
   CMSSMNoFV_mass_eigenstates m1;
   CMSSM_mass_eigenstates m2;
   CMSSMNoFV_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 20.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_models_non_3rd_gen(m1, m2, input);

   const double tol = 1.0e-10;
   const int maxTries = 20;
   const double sin_start = 0.48;

   CMSSMNoFV_weinberg_angle::Parameters parameters_1;
   parameters_1.fermi_constant = Electroweak_constants::gfermi;
   parameters_1.mw_pole = Electroweak_constants::MW;
   parameters_1.mz_pole = Electroweak_constants::MZ;
   parameters_1.mt_pole = 165.0;
   parameters_1.alpha_s = 0.1176;
   CMSSMNoFV_weinberg_angle wein1(&m1, parameters_1);
   wein1.set_number_of_iterations(maxTries);
   wein1.set_precision_goal(tol);
   double sin_theta_1;
   BOOST_REQUIRE_NO_THROW(sin_theta_1 = wein1.calculate(sin_start).first);

   CMSSM_weinberg_angle::Parameters parameters_2;
   parameters_2.fermi_constant = Electroweak_constants::gfermi;
   parameters_2.mw_pole = Electroweak_constants::MW;
   parameters_2.mz_pole = Electroweak_constants::MZ;
   parameters_2.mt_pole = 165.0;
   parameters_2.alpha_s = 0.1176;
   CMSSM_weinberg_angle wein2(&m2, parameters_2);
   wein2.set_number_of_iterations(maxTries);
   wein2.set_precision_goal(tol);
   double sin_theta_2;
   BOOST_REQUIRE_NO_THROW(sin_theta_2 = wein2.calculate(sin_start).first);

   BOOST_CHECK_CLOSE_FRACTION(sin_theta_1, sin_theta_2, 1.0e-10);
}
