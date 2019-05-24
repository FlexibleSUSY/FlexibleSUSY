#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_higgs_loop_corrections

#include <boost/test/unit_test.hpp>
#include "sm_twoloop_mt.hpp"
#include <complex>

#ifdef TSIL_SIZE_DOUBLE
typedef std::complex<double> TSIL_COMPLEXCPP;
#else
typedef std::complex<long double> TSIL_COMPLEXCPP;
#endif

BOOST_AUTO_TEST_CASE( test_fixed_values )
{
   using namespace flexiblesusy::sm_twoloop_mt;

   const TSIL_REAL PI = 3.1415926535897932384626433832795L;
   const TSIL_REAL k = 1.0L/4.0L/4.0L/PI/PI;
   const TSIL_REAL k2 = k*k;
   const TSIL_REAL g3 = 1.166;
   const TSIL_REAL g32 = g3*g3;
   const TSIL_REAL yt = 0.938L;
   const TSIL_REAL t = 26955.07L;
   const TSIL_REAL h = 15647.51L;
   const TSIL_REAL T = 30046.76L;
   const TSIL_REAL s = T;
   const TSIL_REAL qq = 30046.76L;
   const TSIL_REAL Mt = 173.34L;
   const TSIL_REAL mt = 164.18L;
   const TSIL_REAL eps = 1.0e-15L;

   // values obtained with Mathematica
   const TSIL_COMPLEXCPP delta1QCDmat = {11.53533409904775L, 0.0L};
   const TSIL_COMPLEXCPP SigmaSmat = {-108.6728541205668L, 1.508316325898297e-30L};
   const TSIL_COMPLEXCPP SigmaRmat = {0.4197494929184663L, 3.26021741300443e-33L};
   const TSIL_COMPLEXCPP SigmaLmat = {0.883555612555213L, 0.6910278616762643L};   
   const TSIL_COMPLEXCPP SigmaSMtmat = {-122.3083900408013L, -23.3476986375348L};
   const TSIL_COMPLEXCPP SigmaRMtmat = {0.458895852022693L, 0.06744581488147738L};
   const TSIL_COMPLEXCPP SigmaLMtmat = {0.898817852022693L, 0.7584736765577417L};   

   BOOST_CHECK_CLOSE_FRACTION(-0.5*k*g32*std::real(delta1QCDmat),
                              std::real(delta_mt_1loop_as(g3, t, qq)),
                              eps);

   BOOST_CHECK_CLOSE_FRACTION(k*std::real(SigmaSmat),
                              std::real(delta_mt_1loop_at_S(yt, t, h, t, qq)),
                              eps);

   BOOST_CHECK_CLOSE_FRACTION(k*std::real(SigmaRmat),
                              std::real(delta_mt_1loop_at_R(yt, t, h, t, qq)),
                              eps);

   BOOST_CHECK_CLOSE_FRACTION(k*std::real(SigmaLmat),
                              std::real(delta_mt_1loop_at_L(yt, t, h, t, qq)),
                              eps);

   BOOST_CHECK_CLOSE_FRACTION(k*std::real(SigmaSMtmat),
                              std::real(delta_mt_1loop_at_S(yt, t, h, T, qq)),
                              eps);

   BOOST_CHECK_CLOSE_FRACTION(k*std::real(SigmaRMtmat),
                              std::real(delta_mt_1loop_at_R(yt, t, h, T, qq)),
                              eps);

   BOOST_CHECK_CLOSE_FRACTION(k*std::real(SigmaLMtmat),
                              std::real(delta_mt_1loop_at_L(yt, t, h, T, qq)),
                              eps);

   {
      BOOST_TEST_MESSAGE("Testing loop corrections in FlexibleSUSY convention ...");

      const TSIL_COMPLEXCPP Sigma1S = {-0.7745269280201692, -0.1478510288299765};
      const TSIL_COMPLEXCPP Sigma1L_Mt = {0.9866219058664667, 0.832567736200295};
      const TSIL_COMPLEXCPP Sigma1R_Mt = {0.5037246413141052, 0.07403448682214211};
      const TSIL_COMPLEXCPP deltam1tQCD_Mt = {-8.607486434418441, 0.};
      const TSIL_COMPLEXCPP deltamt2QCD_Mt = {-1.381698521871175, 0.};
      const TSIL_COMPLEXCPP deltamt2mixed_Mt = {-0.183268360083289, -0.3363619645269211};
      const TSIL_COMPLEXCPP deltamt2Higgs_Mt = {0.03383958167495036, 0.02128586972099847};
      const TSIL_COMPLEXCPP Sigma2SHiggs = {-0.001811823606387967, -0.01758830494571239};
      const TSIL_COMPLEXCPP Sigma2Smixed = {0.09298641046896923, 0.1148468902046284};

      const TSIL_REAL mt_FS_fixed =
         std::real(Mt
                   + Sigma1S
                   + Sigma2Smixed
                   + Sigma2SHiggs
                   + deltam1tQCD_Mt
                   + deltamt2QCD_Mt
                   + deltamt2mixed_Mt
                   + deltamt2Higgs_Mt
            );

      const TSIL_REAL mt_FS_num =
         std::real(Mt
                   + delta_mt_1loop_at_S(yt, t, h, T, qq)
                   + delta_mt_2loop_as_at_S_flexiblesusy(g3, yt, t, h, T, qq)
                   + delta_mt_2loop_at_at_S_flexiblesusy(yt, t, h, T, qq)
                   + Mt*(
                      + delta_mt_1loop_as(g3, t, qq)
                      + delta_mt_2loop_as_as_flexiblesusy(g3, t, qq)
                      + delta_mt_2loop_as_at_LR_flexiblesusy(g3, yt, t, h, T, qq)
                      + delta_mt_2loop_at_at_LR_flexiblesusy(yt, t, h, T, qq)
                      )
            );

      BOOST_CHECK_CLOSE_FRACTION(mt_FS_fixed, mt_FS_num, eps);
   }

   {
      BOOST_TEST_MESSAGE("Testing loop corrections in SPheno convention ...");

      const TSIL_COMPLEXCPP Sigma1S = { -0.7745269280201692, -0.1478510288299765 };
      const TSIL_COMPLEXCPP Sigma1L_mt = { 0.9344847381167445, 0.788571425691499 };
      const TSIL_COMPLEXCPP Sigma1R_mt = { 0.4771057552264324, 0.07012219941421075 };
      const TSIL_COMPLEXCPP deltamt1QCD_mt = { -8.152631376501787, 0. };
      const TSIL_COMPLEXCPP deltamt2QCDspheno = { -1.713516368402346, 0. };
      const TSIL_COMPLEXCPP deltamt2mixedspheno_mt= { -0.01437079735664272, -0.1323887306086537 };
      const TSIL_COMPLEXCPP deltamt2Higgsspheno_mt = { 0.03211725198423005, -0.003263467271381239 };
      const TSIL_COMPLEXCPP Sigma2SHiggsspheno = { -0.004237440728220745, -0.003011460529185966 };
      const TSIL_COMPLEXCPP Sigma2Smixedspheno = { 0.03417266899507037, -4.742968698291361e-34 };

      const TSIL_REAL mt_SP_fixed =
         std::real(Mt
                   + Sigma1S
                   + Sigma2Smixedspheno
                   + Sigma2SHiggsspheno
                   + deltamt1QCD_mt
                   + deltamt2QCDspheno
                   + deltamt2mixedspheno_mt
                   + deltamt2Higgsspheno_mt
            );

      const TSIL_REAL mt_SP_num =
         std::real(Mt
                   + delta_mt_1loop_at_S(yt, t, h, T, qq)
                   + delta_mt_2loop_as_at_S_spheno(g3, yt, t, h, T, qq)
                   + delta_mt_2loop_at_at_S_spheno(yt, t, h, T, qq)
                   + mt*(
                      + delta_mt_1loop_as(g3, t, qq)
                      + delta_mt_2loop_as_as_spheno(g3, t, qq)
                      + delta_mt_2loop_as_at_LR_spheno(g3, yt, t, h, T, qq)
                      + delta_mt_2loop_at_at_LR_spheno(yt, t, h, T, qq)
                      )
            );

      BOOST_CHECK_CLOSE_FRACTION(mt_SP_fixed, mt_SP_num, eps);
   }
}
