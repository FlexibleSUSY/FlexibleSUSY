#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_standard_model_hh_deriv


#include <boost/test/unit_test.hpp>
#include <iomanip>

#include "derivative.hpp"
#include "wrappers.hpp"
#include "standard_model.hpp"
#include "sm_twoloophiggs.hpp"
#include "lowe.h"

using namespace flexiblesusy;

namespace {

/// flavour integers
constexpr int top    = 0;
constexpr int bottom = 1;
constexpr int tau    = 2;

/**
 * Function returns 1-loop correction to Higgs mass at:
 *
 * O(at)   if flavour == 1
 * O(ab)   if flavour == 2
 * O(atau) if flavour == 3
 */
double one_loop_correction(standard_model::Standard_model& sm, double y,
                           double v, double p, int flavour)
{
   switch (flavour) {
   case top:
      sm.set_Yu(2, 2, y);
      break;
   case bottom:
      sm.set_Yd(2, 2, y);
      break;
   case tau:
      sm.set_Ye(2, 2, y);
      break;
   }

   sm.set_v(v);
   sm.calculate_DRbar_masses();

   const double self_energy = Re(sm.self_energy_hh_1loop(p));
   const double tadpole = Re(sm.tadpole_hh_1loop() / v);
   return -self_energy + tadpole;
}

} // anonymous namespace

/**
 * Test for routines which return the derivatives w.r.t. x, x \in { y_{t,b,\tau}, v, p^2 }
 * of the 1-loop correction at O(at + ab + atau).
 *
 */
BOOST_AUTO_TEST_CASE( test_yuk_derivative )
{

   standard_model::Standard_model sm;
   softsusy::QedQcd qedqcd;
   sm.initialise_from_input(qedqcd);
   sm.set_g1(0.);
   sm.set_g2(0.);
   sm.set_g3(0.);
   sm.set_Lambdax(0.);

   Eigen::Matrix<double,3,3> zero_mat;
   zero_mat.setZero();

   auto sm_yt = sm;
   sm_yt.set_Yd(zero_mat);
   sm_yt.set_Ye(zero_mat);

   auto sm_yb = sm;
   sm_yb.set_Yu(zero_mat);
   sm_yb.set_Ye(zero_mat);

   auto sm_ytau = sm;
   sm_ytau.set_Yu(zero_mat);
   sm_ytau.set_Yd(zero_mat);

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         if (i == 2 && k == 2) {
            continue;
         }
         sm_yt.set_Yu(i, k, 0.);
         sm_yb.set_Yd(i, k, 0.);
         sm_ytau.set_Ye(i, k, 0.);
      }
   }

   sm_yt.calculate_DRbar_masses();
   sm_yb.calculate_DRbar_masses();
   sm_ytau.calculate_DRbar_masses();

   const double over_sqrt2 = 1./1.4142135623730950;
   const double v          = sm.get_v();
   const double yt         = sm.get_Yu(2,2);
   const double yb         = sm.get_Yd(2,2);
   const double ytau       = sm.get_Ye(2,2);
   const double p          =  0.;
   const double Q          = sm.get_scale();

   /*
    * Test for derivatives of the higgs 1-loop corrrection w.r.t. Yukawa couplings y_{t,b,tau}.
    *
    */

   ///  O( at )
   auto fyt1 = [&](double yt) { return one_loop_correction(sm_yt, yt, v, p, top);                                };
   auto fyt2 = [&](double yt) { return sm_twoloophiggs::delta_mh_1loop_at_sm(p, Q, yt*v*over_sqrt2, yt);         };
   auto dfyt = [&](double yt) { return sm_twoloophiggs::delta_mh_1loop_at_sm_deriv_yt(p, Q, yt*v*over_sqrt2, yt);};

   BOOST_CHECK_CLOSE_FRACTION(fyt2(yt), fyt1(yt), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fyt1,yt), dfyt(yt), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fyt2,yt), dfyt(yt), 1e-6);

   ///  O( ab )
   auto fyb1 = [&](double yb) { return one_loop_correction(sm_yb, yb, v, p, bottom);                             };
   auto fyb2 = [&](double yb) { return sm_twoloophiggs::delta_mh_1loop_ab_sm(p, Q, yb*v*over_sqrt2, yb);         };
   auto dfyb = [&](double yb) { return sm_twoloophiggs::delta_mh_1loop_ab_sm_deriv_yb(p, Q, yb*v*over_sqrt2, yb);};

   BOOST_CHECK_CLOSE_FRACTION(fyb2(yb), fyb1(yb), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fyb1,yb), dfyb(yb), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fyb2,yb), dfyb(yb), 1e-6);

   ///  O( atau )
   auto fytau1 = [&](double ytau) { return one_loop_correction(sm_ytau, ytau, v, p, tau);                          };
   auto fytau2 = [&](double ytau) { return sm_twoloophiggs::delta_mh_1loop_atau_sm(p, Q, ytau*v*over_sqrt2, ytau); };
   auto dfytau = [&](double ytau) { return sm_twoloophiggs::delta_mh_1loop_atau_sm_deriv_ytau(p, Q, ytau*v*over_sqrt2, ytau);};

   BOOST_CHECK_CLOSE_FRACTION(fytau2(ytau), fytau1(ytau), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fytau1,ytau), dfytau(ytau), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fytau2,ytau), dfytau(ytau), 1e-6);


   /*
    * Test for derivatives of the higgs 1-loop corrrection w.r.t. VEV.
    *
    */

   ///  O( at )
   auto fytv1 = [&](double v) { return one_loop_correction(sm_yt, yt, v, p, top);                               };
   auto fytv2 = [&](double v) { return sm_twoloophiggs::delta_mh_1loop_at_sm(p, Q, yt*v*over_sqrt2, yt);        };
   auto dfytv = [&](double v) { return sm_twoloophiggs::delta_mh_1loop_at_sm_deriv_v(p, Q, yt*v*over_sqrt2, yt);};

   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fytv1,v), dfytv(v), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fytv2,v), dfytv(v), 1e-6);

   ///  O( ab)
   auto fybv1 = [&](double v) { return one_loop_correction(sm_yb, yb, v, p, bottom);                            };
   auto fybv2 = [&](double v) { return sm_twoloophiggs::delta_mh_1loop_ab_sm(p, Q, yb*v*over_sqrt2, yb);        };
   auto dfybv = [&](double v) { return sm_twoloophiggs::delta_mh_1loop_ab_sm_deriv_v(p, Q, yb*v*over_sqrt2, yb);};

   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fybv1,v), dfybv(v), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fybv2,v), dfybv(v), 1e-6);

   ///  O( atau)
   auto fytauv1 = [&](double v) { return one_loop_correction(sm_ytau, ytau, v, p, tau);                                 };
   auto fytauv2 = [&](double v) { return sm_twoloophiggs::delta_mh_1loop_atau_sm(p, Q, ytau*v*over_sqrt2, ytau);        };
   auto dfytauv = [&](double v) { return sm_twoloophiggs::delta_mh_1loop_atau_sm_deriv_v(p, Q, ytau*v*over_sqrt2, ytau);};

   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fytauv1,v), dfytauv(v), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(fytauv2,v), dfytauv(v), 1e-6);


   /*
    * Test for derivativesof the higgs 1-loop corrrection w.r.t. p^2 .
    * To succeed the test, FS was configured with fflite library.
    *
    */

   ///  O( at )
   auto fpt1 = [&](double p2) { return one_loop_correction(sm_yt, yt, v, Sqrt(p2), top);};
   auto fpt2 = [&](double p2) { return sm_twoloophiggs::delta_mh_1loop_at_sm(Sqrt(p2), Q, yt*v*over_sqrt2, yt);};
   auto dfpt = [&](double p2) { return sm_twoloophiggs::delta_mh_1loop_at_sm_deriv_p2(AbsSqrt(p2), Q, yt*v*over_sqrt2, yt);};
   BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(fpt1,p, 1e-5), dfpt(p), 2e-4);
   BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(fpt2,p, 1e-5), dfpt(p), 1e-4);

   ///  O( ab )
   auto fpb1 = [&](double p2) { return one_loop_correction(sm_yb, yb, v, Sqrt(p2), bottom);};
   auto fpb2 = [&](double p2) { return sm_twoloophiggs::delta_mh_1loop_ab_sm(Sqrt(p2), Q, yb*v*over_sqrt2, yb);};
   auto dfpb = [&](double p2) { return sm_twoloophiggs::delta_mh_1loop_ab_sm_deriv_p2(AbsSqrt(p2), Q, yb*v*over_sqrt2, yb);};

   BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(fpb1,p, 1e-5), dfpb(p), 1e-4);
   BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(fpb2,p, 1e-5), dfpb(p), 1e-4);

   ///  O( atau )
   auto fptau1 = [&](double p2) { return one_loop_correction(sm_ytau, ytau, v, Sqrt(p2), tau);};
   auto fptau2 = [&](double p2) { return sm_twoloophiggs::delta_mh_1loop_atau_sm(Sqrt(p2), Q, ytau*v*over_sqrt2, ytau);};
   auto dfptau = [&](double p2) { return sm_twoloophiggs::delta_mh_1loop_atau_sm_deriv_p2(AbsSqrt(p2), Q, ytau*v*over_sqrt2, ytau);};

   BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(fptau1,p, 1e-5), dfptau(p), 1e-4);
   BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(fptau2,p, 1e-5), dfptau(p), 1e-4);

}
