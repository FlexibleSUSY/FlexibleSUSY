#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCPVSemiAnalytic_semi_analytic_solutions

#include <boost/test/unit_test.hpp>

#include "test.hpp"
#include "test_CMSSMCPVSemiAnalytic.hpp"
#include "CMSSMCPVSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMCPVSemiAnalytic_coefficients )
{
   CMSSMCPVSemiAnalytic_input_parameters input;
   CMSSMCPVSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMCPVSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMCPVSemiAnalytic(model, values);

   CMSSMCPVSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CMSSMCPVSemiAnalytic_mass_eigenstates coeffs_model(model);
   solns.evaluate_solutions(coeffs_model);

   BOOST_CHECK_CLOSE_FRACTION(Re(model.get_MassB()), Re(coeffs_model.get_MassB()), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(Im(model.get_MassB()), Im(coeffs_model.get_MassB()), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(Re(model.get_MassWB()), Re(coeffs_model.get_MassWB()), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(Im(model.get_MassWB()), Im(coeffs_model.get_MassWB()), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(Re(model.get_MassG()), Re(coeffs_model.get_MassG()), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(Im(model.get_MassG()), Im(coeffs_model.get_MassG()), 1.0e-3);

   BOOST_CHECK_CLOSE_FRACTION(Re(model.get_BMu()), Re(coeffs_model.get_BMu()), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(Im(model.get_BMu()), Im(coeffs_model.get_BMu()), 1.0e-3);

   BOOST_CHECK_CLOSE_FRACTION(model.get_mHd2(), coeffs_model.get_mHd2(), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_mHu2(), coeffs_model.get_mHu2(), 1.0e-3);

   TEST_CLOSE_REL(model.get_TYu().real(), coeffs_model.get_TYu().real(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYu().imag(), coeffs_model.get_TYu().imag(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYd().real(), coeffs_model.get_TYd().real(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYd().imag(), coeffs_model.get_TYd().imag(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYe().real(), coeffs_model.get_TYe().real(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYe().imag(), coeffs_model.get_TYe().imag(), 1.0e-3);

   TEST_CLOSE_REL(model.get_mq2().real(), coeffs_model.get_mq2().real(), 1.0e-2);
   TEST_CLOSE_REL(model.get_mq2().imag(), coeffs_model.get_mq2().imag(), 1.0e-2);
   TEST_CLOSE_REL(model.get_mu2().real(), coeffs_model.get_mu2().real(), 1.0e-2);
   TEST_CLOSE_REL(model.get_mu2().imag(), coeffs_model.get_mu2().imag(), 1.0e-2);
   TEST_CLOSE_REL(model.get_md2().real(), coeffs_model.get_md2().real(), 1.0e-2);
   TEST_CLOSE_REL(model.get_md2().imag(), coeffs_model.get_md2().imag(), 1.0e-2);
   TEST_CLOSE_REL(model.get_ml2().real(), coeffs_model.get_ml2().real(), 1.0e-2);
   TEST_CLOSE_REL(model.get_ml2().imag(), coeffs_model.get_ml2().imag(), 1.0e-2);
   TEST_CLOSE_REL(model.get_me2().real(), coeffs_model.get_me2().real(), 1.0e-2);
   TEST_CLOSE_REL(model.get_me2().imag(), coeffs_model.get_me2().imag(), 1.0e-2);

   BOOST_CHECK_EQUAL(get_errors(), 0);
}
