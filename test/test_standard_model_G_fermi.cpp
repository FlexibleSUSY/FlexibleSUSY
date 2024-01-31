#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_standard_model_G_fermi

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "standard_model.hpp"
#include "lowe.h"

using namespace flexiblesusy;


double calculate_G_fermi(const softsusy::QedQcd& qedqcd, int loops)
{
   Threshold_corrections tc;
   tc.sin_theta_w = loops;

   BOOST_TEST_MESSAGE("initialize SM ...");
   standard_model::Standard_model sm;
   sm.set_pole_mass_loop_order(loops);
   sm.set_ewsb_loop_order(loops);
   sm.set_threshold_corrections(tc);

   sm.initialise_from_input(qedqcd);
   sm.calculate_DRbar_masses();
   sm.solve_ewsb();
   sm.calculate_pole_masses();

   BOOST_TEST_MESSAGE("calculate G_Fermi ...");
   return sm.calculate_G_fermi(qedqcd);
}


BOOST_AUTO_TEST_CASE( test_G_fermi )
{
   softsusy::QedQcd qedqcd;
   const double G_fermi_input = qedqcd.displayFermiConstant();
   const double G_fermi_0l = calculate_G_fermi(qedqcd, 0);
   const double G_fermi_1l = calculate_G_fermi(qedqcd, 1);
   const double G_fermi_2l = calculate_G_fermi(qedqcd, 2);

   BOOST_TEST_MESSAGE("G_fermi_input = " << G_fermi_input);
   BOOST_TEST_MESSAGE("G_fermi_0l = " << G_fermi_0l);
   BOOST_TEST_MESSAGE("G_fermi_1l = " << G_fermi_1l);
   BOOST_TEST_MESSAGE("G_fermi_2l = " << G_fermi_2l);

   BOOST_CHECK_CLOSE_FRACTION(G_fermi_input, G_fermi_0l, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(G_fermi_input, G_fermi_1l, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(G_fermi_input, G_fermi_2l, 1.0e-5);
}
