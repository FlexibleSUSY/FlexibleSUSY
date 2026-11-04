
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_self_energies

#include <boost/test/unit_test.hpp>

#include "test_NMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;


BOOST_AUTO_TEST_CASE( test_NMSSM_self_energy_neutral_higgs )
{
   NMSSM_input_parameters input;
   input.m0 = 300.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM_const(m, s, input);
   m.do_force_output(true);

   // initial guess
   m.set_Kappa(0.1);
   m.set_vS(5000.);
   m.set_ms2(-Sqr(input.m0));
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   s.setKappa(m.get_Kappa());
   s.setSvev(m.get_vS());
   s.setMsSquared(m.get_ms2());
   s.setMh1Squared(m.get_mHd2());
   s.setMh2Squared(m.get_mHu2());

   s.calcDrBarPars();
   m.calculate_DRbar_masses();

   // check tree-level
   DoubleVector hh_ss(s.displayDrBarPars().mh0);
   DoubleVector Ah_ss(s.displayDrBarPars().mA0);
   Eigen::Array<double,3,1> hh_fs(m.get_Mhh());
   Eigen::Array<double,3,1> Ah_fs(m.get_MAh());

   BOOST_CHECK_CLOSE_FRACTION(hh_ss(1), hh_fs(0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(2), hh_fs(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(3), hh_fs(2), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(Ah_ss(1), Ah_fs(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(Ah_ss(2), Ah_fs(2), 1.0e-10);

   // solve EWSB
   int loop_order = 2;
   softsusy::numRewsbLoops = loop_order;
   softsusy::numHiggsMassLoops = loop_order;
   m.set_ewsb_loop_order(loop_order);
   m.set_pole_mass_loop_order(loop_order);

   ensure_n_loop_ewsb(m, loop_order);
   ensure_n_loop_ewsb(s, loop_order);

   const double kappa_ss = s.displayKappa();
   const double vS_ss    = s.displaySvev();
   const double ms2_ss   = s.displayMsSquared();

   const double kappa_fs = m.get_Kappa();
   const double vS_fs    = m.get_vS();
   const double ms2_fs   = m.get_ms2();

   BOOST_CHECK_CLOSE_FRACTION(kappa_ss, kappa_fs, 5.5e-9);
   BOOST_CHECK_CLOSE_FRACTION(vS_ss   , vS_fs   , 3.6e-8);
   BOOST_CHECK_CLOSE_FRACTION(ms2_ss  , ms2_fs  , 7.5e-8);

   const int accuracy = loop_order;
   const double piWWT = 0., pizztMS = 0.;

   // setting initial conditions for the iteration
   softsusy::sPhysical physical(s.displayDrBarPars());
   physical.mh0(1) = s.displayDrBarPars().mh0(1);
   physical.mh0(2) = s.displayDrBarPars().mh0(2);
   physical.mh0(3) = s.displayDrBarPars().mh0(3);
   physical.mA0(1) = s.displayDrBarPars().mA0(1);
   physical.mA0(2) = s.displayDrBarPars().mA0(2);

   s.setPhys(physical); // initialization
   s.higgs(accuracy, piWWT, pizztMS, physical); // does one iteration
   s.setPhys(physical); // now modified

   m.set_precision(std::pow(10,-0.1)); // performs 1 iteration
   m.calculate_Mhh_pole();
   m.calculate_MAh_pole();

   hh_ss = s.displayPhys().mh0;
   hh_fs = m.get_physical().Mhh;
   Ah_ss = s.displayPhys().mA0;
   Ah_fs = m.get_physical().MAh;

   BOOST_CHECK_CLOSE_FRACTION(hh_ss(1), hh_fs(0), 2.0e-04);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(2), hh_fs(1), 4.0e-04);
   BOOST_CHECK_CLOSE_FRACTION(hh_ss(3), hh_fs(2), 7.0e-04);
   BOOST_CHECK_CLOSE_FRACTION(Ah_ss(1), Ah_fs(1), 9.0e-05);
   BOOST_CHECK_CLOSE_FRACTION(Ah_ss(2), Ah_fs(2), 4.0e-04);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_self_energy_photon_Z )
{
   NMSSM_input_parameters input;
   input.m0 = 300.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM_const(m, s, input);
   m.do_force_output(true);

   // initial guess
   m.set_Kappa(0.1);
   m.set_vS(5000.);
   m.set_ms2(-Sqr(input.m0));
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   m.calculate_DRbar_masses();

   const double p = 100;

   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VPVZ_1loop(p).real(), 109.34557661080734, 7.0e-15);
   BOOST_CHECK_SMALL(m.self_energy_VPVZ_1loop(p).imag(), 0.);
}
