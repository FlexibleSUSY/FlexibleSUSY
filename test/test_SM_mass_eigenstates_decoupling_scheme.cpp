#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_mass_eigenstates_decoupling_scheme

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_mass_eigenstates_decoupling_scheme.hpp"
#include "standard_model.hpp"

#include <memory>
#include <utility>

#define COMPARE_0(a,b,par,eps)                                          \
   do {                                                                 \
      BOOST_CHECK_CLOSE(std::real(a->get_##par()), std::real(b->get_##par()), eps); \
      BOOST_CHECK_CLOSE(std::imag(a->get_##par()), std::imag(b->get_##par()), eps); \
   } while (false)

#define COMPARE_1(a,b,par,N,eps)                                        \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         BOOST_CHECK_CLOSE(std::real(a->get_##par()(i)), std::real(b->get_##par()(i)), eps); \
         BOOST_CHECK_CLOSE(std::imag(a->get_##par()(i)), std::imag(b->get_##par()(i)), eps); \
      }                                                                 \
   } while (false)

#define COMPARE_2(a,b,par,N,M,eps)                                      \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         for (int k = 0; k < (M); k++) {                                \
            BOOST_CHECK_CLOSE(std::real(a->get_##par()(i,k)), std::real(b->get_##par()(i,k)), eps); \
            BOOST_CHECK_CLOSE(std::imag(a->get_##par()(i,k)), std::imag(b->get_##par()(i,k)), eps); \
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

using Model_ifc_ptr = std::unique_ptr<SM_mass_eigenstates_interface>;
using Model_ifc_ptrs = std::pair<Model_ifc_ptr, Model_ifc_ptr>;
using Model_ptrs = std::pair<std::unique_ptr<SM_mass_eigenstates>, std::unique_ptr<SM_mass_eigenstates_decoupling_scheme>>;

Model_ptrs make_model_ptrs(const SM_input_parameters& input)
{
   auto me = std::make_unique<SM_mass_eigenstates>();
   setup_SM_const(*me, input);

   auto dec = std::make_unique<SM_mass_eigenstates_decoupling_scheme>();
   dec->set_g1(me->get_g1());
   dec->set_g2(me->get_g2());
   dec->set_g3(me->get_g3());
   dec->set_Yu(me->get_Yu());
   dec->set_Yd(me->get_Yd());
   dec->set_Ye(me->get_Ye());
   dec->set_Lambdax(me->get_Lambdax());
   dec->set_v(me->get_v());
   // dec->set_mu2(me->get_mu2()); // determine mu2 from EWSB conditions
   // dec->set_scale(me->get_scale()); // scale is not defined for decoupling scheme

   return std::make_pair(std::move(me), std::move(dec));
}

Model_ifc_ptrs make_model_ifc_ptrs(const SM_input_parameters& input)
{
   auto ptrs = make_model_ptrs(input);

   return std::make_pair(std::move(ptrs.first), std::move(ptrs.second));
}

std::unique_ptr<standard_model::Standard_model> make_sm(
   const SM_input_parameters& input)
{
   Eigen::Matrix<double,3,3> Yu, Yd, Ye;
   Yu << 0, 0, 0, 0, 0, 0, 0, 0, 0.8;
   Ye << 0, 0, 0, 0, 0, 0, 0, 0, 0.1;
   Yd << 0, 0, 0, 0, 0, 0, 0, 0, 0.3;

   auto sm = std::make_unique<standard_model::Standard_model>(
      91.           , // scale
      0             , // loops
      0             , // thresholds
      0.45          , // g1
      0.61          , // g2
      1.06          , // g3
      input.LambdaIN, // lambda
      Yu            ,
      Yd            ,
      Ye            ,
      0.            , // mu2
      245.            // v
   );

   sm->solve_ewsb_tree_level();
   sm->calculate_DRbar_masses();

   return std::move(sm);
}

BOOST_AUTO_TEST_CASE( test_SM_mass_eigenstates_conversion )
{
   const auto eps = std::numeric_limits<double>::epsilon();

   SM_input_parameters input;
   input.LambdaIN = 0.24;
   input.Qin = 91.0;
   input.QEWSB = 173.34;

   auto models = make_model_ptrs(input);
   auto model = std::move(std::get<0>(models));

   model->solve_ewsb_equations_tree_level();
   model->calculate_tree_level_mass_spectrum();
   model->calculate_pole_mass_spectrum();

   auto dec = std::make_unique<SM_mass_eigenstates_decoupling_scheme>(input);
   dec->fill_from(*model);

   // parameters
   COMPARE_0(model, dec, g1, eps);
   COMPARE_0(model, dec, g2, eps);
   COMPARE_0(model, dec, g3, eps);
   COMPARE_0(model, dec, Lambdax, eps);
   COMPARE_0(model, dec, v, eps);
   COMPARE_0(model, dec, mu2, eps);
   COMPARE_2(model, dec, Yu, 3, 3, eps);
   COMPARE_2(model, dec, Yd, 3, 3, eps);
   COMPARE_2(model, dec, Ye, 3, 3, eps);

   // tree-level masses
   COMPARE_0(model, dec, MVP, eps);
   COMPARE_0(model, dec, MVZ, eps);
   COMPARE_0(model, dec, MVWp, eps);
   COMPARE_0(model, dec, MVG, eps);
   COMPARE_0(model, dec, Mhh, eps);
   COMPARE_0(model, dec, MHp, eps);
   COMPARE_0(model, dec, MAh, eps);
   COMPARE_1(model, dec, MFu, 3, eps);
   COMPARE_1(model, dec, MFd, 3, eps);
   COMPARE_1(model, dec, MFe, 3, eps);
   COMPARE_1(model, dec, MFv, 3, eps);

   COMPARE_2(model, dec, Vd, 3, 3, eps);
   COMPARE_2(model, dec, Ud, 3, 3, eps);
   COMPARE_2(model, dec, Vu, 3, 3, eps);
   COMPARE_2(model, dec, Uu, 3, 3, eps);
   COMPARE_2(model, dec, Ve, 3, 3, eps);
   COMPARE_2(model, dec, Ue, 3, 3, eps);
   COMPARE_2(model, dec, ZZ, 2, 2, eps);

   // pole masses
   const auto pole_1 = model->get_physical();
   const auto pole_2 = dec->get_physical();

   COMPARE_POLE_0(pole_1, pole_2, MVP, eps);
   COMPARE_POLE_0(pole_1, pole_2, MVZ, eps);
   COMPARE_POLE_0(pole_1, pole_2, MVWp, eps);
   COMPARE_POLE_0(pole_1, pole_2, MVG, eps);
   COMPARE_POLE_0(pole_1, pole_2, Mhh, eps);
   COMPARE_POLE_0(pole_1, pole_2, MHp, eps);
   COMPARE_POLE_0(pole_1, pole_2, MAh, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFu, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFd, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFe, 3, eps);
   COMPARE_POLE_1(pole_1, pole_2, MFv, 3, eps);

   COMPARE_POLE_2(pole_1, pole_2, Vd, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, Ud, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, Vu, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, Uu, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, Ve, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, Ue, 3, 3, eps);
   COMPARE_POLE_2(pole_1, pole_2, ZZ, 2, 2, eps);
}

BOOST_AUTO_TEST_CASE( test_SM_mass_eigenstates_decoupling_scheme_matching )
{
   const auto eps = 1e-13;

   SM_input_parameters input;
   input.LambdaIN = 0.24;
   input.Qin = 91.0;
   input.QEWSB = 173.34;

   auto dec = std::make_unique<SM_mass_eigenstates_decoupling_scheme>(input);

   {
      auto models = make_model_ptrs(input);
      auto model = std::move(std::get<0>(models));

      model->solve_ewsb_equations_tree_level();
      model->calculate_tree_level_mass_spectrum();
      model->calculate_pole_mass_spectrum();

      dec->fill_from(*model);
   }

   auto sm = make_sm(input);

   dec->fill_from(*sm);

   // parameters
   COMPARE_0(sm, dec, g1, eps);
   COMPARE_0(sm, dec, g2, eps);
   COMPARE_0(sm, dec, g3, eps);
   // COMPARE_0(sm, dec, Lambdax, eps);
   COMPARE_0(sm, dec, v, eps);
   COMPARE_2(sm, dec, Yu, 3, 3, eps);
   COMPARE_2(sm, dec, Yd, 3, 3, eps);
   COMPARE_2(sm, dec, Ye, 3, 3, eps);
}
