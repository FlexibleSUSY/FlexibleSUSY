#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_standard_model_G_fermi

#include <boost/test/unit_test.hpp>

#define private public

#include "config.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "standard_model.hpp"
#include "lowe.h"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_G_fermi )
{
   Threshold_corrections tc;
   tc.sin_theta_w = 0; // no loop corrections to sin(theta_W)

   standard_model::Standard_model sm;
   sm.set_pole_mass_loop_order(0);
   sm.set_ewsb_loop_order(0);
   sm.set_threshold_corrections(tc);

   softsusy::QedQcd qedqcd;
   sm.initialise_from_input(qedqcd);
   sm.calculate_DRbar_masses();

   const double G_fermi_input = qedqcd.displayFermiConstant();
   const double G_fermi_pred  = sm.calculate_G_fermi(qedqcd);

   const double gY                  = sm.get_g1() * 0.7745966692414834;
   const double g2                  = sm.get_g2();
   const double e                   = gY*g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alpha_em_drbar      = Sqr(e) / (4*Pi);
   const double cos_theta           = sm.get_MVWp()/sm.get_MVZ();
   const double theta_w             = ArcCos(cos_theta);
   const double sin_theta           = Sin(theta_w);
   const double MZ2pole             = Sqr(qedqcd.displayPoleMZ());
   const double sqrt2               = 1.4142135623730950;
   const double Pi                  = 3.1415926535897932;

   const double G_fermi_0l = Pi * alpha_em_drbar / (AbsSqr(sin_theta*cos_theta) * sqrt2 * MZ2pole);

   BOOST_TEST_MESSAGE("G_fermi_input = " << G_fermi_input);
   BOOST_TEST_MESSAGE("G_fermi_0l = " << G_fermi_0l);
   BOOST_TEST_MESSAGE("G_fermi_pred = " << G_fermi_pred);

   BOOST_CHECK_CLOSE_FRACTION(G_fermi_input, G_fermi_0l, 1.0e-2);

   BOOST_CHECK_CLOSE_FRACTION(G_fermi_input, G_fermi_pred, 1.0e-5);
}
