#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_mass_eigenstates_decoupling_scheme

#include <boost/test/unit_test.hpp>

#include "test_CMSSM.hpp"
#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_mass_eigenstates_decoupling_scheme.hpp"
#include "standard_model.hpp"

#include <memory>
#include <utility>

#define COMPARE_0(a,b,par,eps)                                          \
   do {                                                                 \
      BOOST_CHECK_CLOSE(std::real(a.get_##par()), std::real(b.get_##par()), eps); \
      BOOST_CHECK_CLOSE(std::imag(a.get_##par()), std::imag(b.get_##par()), eps); \
   } while (false)

#define COMPARE_1(a,b,par,N,eps)                                        \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         BOOST_CHECK_CLOSE(std::real(a.get_##par()(i)), std::real(b.get_##par()(i)), eps); \
         BOOST_CHECK_CLOSE(std::imag(a.get_##par()(i)), std::imag(b.get_##par()(i)), eps); \
      }                                                                 \
   } while (false)

#define COMPARE_2(a,b,par,N,M,eps)                                      \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         for (int k = 0; k < (M); k++) {                                \
            BOOST_CHECK_CLOSE(std::real(a.get_##par()(i,k)), std::real(b.get_##par()(i,k)), eps); \
            BOOST_CHECK_CLOSE(std::imag(a.get_##par()(i,k)), std::imag(b.get_##par()(i,k)), eps); \
         }                                                              \
      }                                                                 \
   } while (false)

#define COMPARE_POLE_0(a,b,par,eps)                                     \
   do {                                                                 \
      BOOST_CHECK_CLOSE(std::real(a.par), std::real(b.par), eps); \
      BOOST_CHECK_CLOSE(std::imag(a.par), std::imag(b.par), eps); \
   } while (false)

#define COMPARE_POLE_1(a,b,par,N,eps)                                   \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         BOOST_CHECK_CLOSE(std::real(a.par(i)), std::real(b.par(i)), eps); \
         BOOST_CHECK_CLOSE(std::imag(a.par(i)), std::imag(b.par(i)), eps); \
      }                                                                 \
   } while (false)

#define COMPARE_POLE_2(a,b,par,N,M,eps)                                 \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         for (int k = 0; k < (M); k++) {                                \
            BOOST_CHECK_CLOSE(std::real(a.par(i,k)), std::real(b.par(i,k)), eps); \
            BOOST_CHECK_CLOSE(std::imag(a.par(i,k)), std::imag(b.par(i,k)), eps); \
         }                                                              \
      }                                                                 \
   } while (false)

using namespace flexiblesusy;

CMSSM_mass_eigenstates make_model(const CMSSM_input_parameters& input)
{
   return setup_CMSSM(input);
}

standard_model::Standard_model make_sm(const CMSSM_input_parameters& input)
{
   Eigen::Matrix<double,3,3> Yu, Yd, Ye;
   Yu << 0, 0, 0, 0, 0, 0, 0, 0, 0.8;
   Ye << 0, 0, 0, 0, 0, 0, 0, 0, 0.1;
   Yd << 0, 0, 0, 0, 0, 0, 0, 0, 0.3;

   standard_model::Standard_model sm(
      91. , // scale
      0   , // loops
      0   , // thresholds
      0.45, // g1
      0.61, // g2
      1.06, // g3
      0.24, // lambda
      Yu  ,
      Yd  ,
      Ye  ,
      0.  , // mu2
      245.  // v
   );

   sm.solve_ewsb_tree_level();
   sm.calculate_DRbar_masses();

   return sm;
}

/// checks that fill routine works correctly
BOOST_AUTO_TEST_CASE( test_CMSSM_mass_eigenstates_conversion )
{
   const auto eps = std::numeric_limits<double>::epsilon();

   CMSSM_input_parameters input;
   input.m0 = 250.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   auto model = make_model(input);

   model.solve_ewsb_equations_tree_level();
   model.calculate_tree_level_mass_spectrum();
   model.calculate_pole_mass_spectrum();

   CMSSM_mass_eigenstates_decoupling_scheme dec(input);
   dec.fill_from(model);

   // parameters
   COMPARE_0(model, dec, g1, eps);
   COMPARE_0(model, dec, g2, eps);
   COMPARE_0(model, dec, g3, eps);
   COMPARE_0(model, dec, vu, eps);
   COMPARE_0(model, dec, vd, eps);
   COMPARE_2(model, dec, Yu, 3, 3, eps);
   COMPARE_2(model, dec, Yd, 3, 3, eps);
   COMPARE_2(model, dec, Ye, 3, 3, eps);

   // tree-level masses
   COMPARE_0(model, dec, MVP, eps);
   COMPARE_0(model, dec, MVZ, eps);
   COMPARE_0(model, dec, MVWm, eps);
   COMPARE_0(model, dec, MVG, eps);
   COMPARE_1(model, dec, Mhh, 2, eps);
   COMPARE_1(model, dec, MHpm, 2, eps);
   COMPARE_1(model, dec, MAh, 2, eps);
   COMPARE_1(model, dec, MFu, 3, eps);
   COMPARE_1(model, dec, MFd, 3, eps);
   COMPARE_1(model, dec, MFe, 3, eps);
   COMPARE_1(model, dec, MFv, 3, eps);
   COMPARE_1(model, dec, MSu, 6, eps);
   COMPARE_1(model, dec, MSd, 6, eps);
   COMPARE_1(model, dec, MSe, 6, eps);
   COMPARE_1(model, dec, MSv, 3, eps);
   COMPARE_1(model, dec, MChi, 4, eps);
   COMPARE_1(model, dec, MCha, 2, eps);
   COMPARE_0(model, dec, MGlu, eps);

   COMPARE_2(model, dec, ZD, 6, 6, eps);
   COMPARE_2(model, dec, ZV, 3, 3, eps);
   COMPARE_2(model, dec, ZU, 6, 6, eps);
   COMPARE_2(model, dec, ZE, 6, 6, eps);
   COMPARE_2(model, dec, ZH, 2, 2, eps);
   COMPARE_2(model, dec, ZA, 2, 2, eps);
   COMPARE_2(model, dec, ZP, 2, 2, eps);
   COMPARE_2(model, dec, ZN, 4, 4, eps);
   COMPARE_2(model, dec, UM, 2, 2, eps);
   COMPARE_2(model, dec, UP, 2, 2, eps);
   COMPARE_2(model, dec, ZEL, 3, 3, eps);
   COMPARE_2(model, dec, ZER, 3, 3, eps);
   COMPARE_2(model, dec, ZDL, 3, 3, eps);
   COMPARE_2(model, dec, ZDR, 3, 3, eps);
   COMPARE_2(model, dec, ZUL, 3, 3, eps);
   COMPARE_2(model, dec, ZUR, 3, 3, eps);
   COMPARE_2(model, dec, ZZ, 2, 2, eps);

   // pole masses
   const auto pole_1 = model.get_physical();
   const auto pole_2 = dec.get_physical();

   COMPARE_POLE_0(pole_1, pole_2, MVP, eps);
   COMPARE_POLE_0(pole_1, pole_2, MVZ, eps);
   COMPARE_POLE_0(pole_1, pole_2, MVWm, eps);
   COMPARE_POLE_0(pole_1, pole_2, MVG, eps);
   COMPARE_POLE_1(pole_1, pole_2, Mhh, 2, eps);
   COMPARE_POLE_1(pole_1, pole_2, MHpm, 2, eps);
   COMPARE_POLE_1(pole_1, pole_2, MAh, 2, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFu, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFd, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFe, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFv, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MSu, 6, eps);
   COMPARE_POLE_1(pole_1, pole_2, MSd, 6, eps);
   COMPARE_POLE_1(pole_1, pole_2, MSe, 6, eps);
   COMPARE_POLE_1(pole_1, pole_2, MSv, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MChi, 4, eps);
   COMPARE_POLE_1(pole_1, pole_2, MCha, 2, eps);
   COMPARE_POLE_0(pole_1, pole_2, MGlu, eps);

   COMPARE_POLE_2(pole_1, pole_2, ZD, 6, 6, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZV, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZU, 6, 6, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZE, 6, 6, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZH, 2, 2, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZA, 2, 2, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZP, 2, 2, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZN, 4, 4, eps);
   COMPARE_POLE_2(pole_1, pole_2, UM, 2, 2, eps);
   COMPARE_POLE_2(pole_1, pole_2, UP, 2, 2, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZEL, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZER, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZDL, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZDR, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZUL, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZUR, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZZ, 2, 2, eps);
}

/// checks that SM-like parameters are equal after matching
BOOST_AUTO_TEST_CASE( test_CMSSM_mass_eigenstates_decoupling_scheme_matching )
{
   const auto eps = 1e-12;

   CMSSM_input_parameters input;
   input.m0 = 250.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   CMSSM_mass_eigenstates_decoupling_scheme dec(input);

   {
      auto model = make_model(input);

      model.solve_ewsb_equations_tree_level();
      model.calculate_tree_level_mass_spectrum();
      model.calculate_pole_mass_spectrum();

      dec.fill_from(model);
   }

   auto sm = make_sm(input);

   dec.fill_from(sm);

   const double vu = dec.get_vu();
   const double vd = dec.get_vd();
   const double v = std::sqrt(vu*vu + vd*vd);
   const double beta = std::atan(vu/vd);
   const double sb = std::sin(beta);
   const double cb = std::cos(beta);

   // check that SM-like parameters are equal
   COMPARE_0(sm, dec, g1, eps);
   COMPARE_0(sm, dec, g2, eps);
   COMPARE_0(sm, dec, g3, eps);

   BOOST_CHECK_CLOSE(sm.get_v(), v, eps);

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE(sm.get_Yu(i,k), dec.get_Yu(i,k)*sb, eps);
         BOOST_CHECK_CLOSE(sm.get_Yd(i,k), dec.get_Yd(i,k)*cb, eps);
         BOOST_CHECK_CLOSE(sm.get_Ye(i,k), dec.get_Ye(i,k)*cb, eps);
      }
   }
}
