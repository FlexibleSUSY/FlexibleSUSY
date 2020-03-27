#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_mass_eigenstates_interface

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_mass_eigenstates_decoupling_scheme.hpp"

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

#define COMPARE_POLE_PTR_0(a,b,par,eps)                                 \
   do {                                                                 \
      BOOST_CHECK_CLOSE(std::real(a.par), std::real(b->get_##par()), eps); \
      BOOST_CHECK_CLOSE(std::imag(a.par), std::imag(b->get_##par()), eps); \
   } while (false)

#define COMPARE_POLE_PTR_1(a,b,par,N,eps)                               \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         BOOST_CHECK_CLOSE(std::real(a.par(i)), std::real(b->get_##par()(i)), eps); \
         BOOST_CHECK_CLOSE(std::imag(a.par(i)), std::imag(b->get_##par()(i)), eps); \
      }                                                                 \
   } while (false)

#define COMPARE_POLE_PTR_2(a,b,par,N,M,eps)                             \
   do {                                                                 \
      for (int i = 0; i < (N); i++) {                                   \
         for (int k = 0; k < (M); k++) {                                \
            BOOST_CHECK_CLOSE(std::real(a.par(i,k)), std::real(b->get_##par()(i,k)), eps); \
            BOOST_CHECK_CLOSE(std::imag(a.par(i,k)), std::imag(b->get_##par()(i,k)), eps); \
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

BOOST_AUTO_TEST_CASE( test_SM_mass_eigenstates )
{
   const auto eps = std::numeric_limits<double>::epsilon();

   SM_input_parameters input;
   input.LambdaIN = 0.24;
   input.Qin = 91.0;
   input.QEWSB = 173.34;

   auto models = make_model_ifc_ptrs(input);
   auto model_1 = std::move(std::get<0>(models));
   auto model_2 = std::move(std::get<1>(models));

   model_1->solve_ewsb_equations_tree_level();
   model_1->calculate_tree_level_mass_spectrum();
   model_1->calculate_pole_mass_spectrum();
   model_2->solve_ewsb_equations_tree_level();
   model_2->calculate_tree_level_mass_spectrum();
   model_2->calculate_pole_mass_spectrum();

   // parameters
   COMPARE_0(model_1, model_2, g1, eps);
   COMPARE_0(model_1, model_2, g2, eps);
   COMPARE_0(model_1, model_2, g3, eps);
   COMPARE_0(model_1, model_2, Lambdax, eps);
   COMPARE_0(model_1, model_2, v, eps);
   COMPARE_0(model_1, model_2, mu2, eps);
   COMPARE_2(model_1, model_2, Yu, 3, 3, eps);
   COMPARE_2(model_1, model_2, Yd, 3, 3, eps);
   COMPARE_2(model_1, model_2, Ye, 3, 3, eps);

   // tree-level masses
   COMPARE_0(model_1, model_2, MVP, eps);
   COMPARE_0(model_1, model_2, MVZ, eps);
   COMPARE_0(model_1, model_2, MVWp, eps);
   COMPARE_0(model_1, model_2, MVG, eps);
   COMPARE_0(model_1, model_2, Mhh, eps);
   COMPARE_0(model_1, model_2, MHp, eps);
   COMPARE_0(model_1, model_2, MAh, eps);
   COMPARE_1(model_1, model_2, MFu, 3, eps);
   COMPARE_1(model_1, model_2, MFd, 3, eps);
   COMPARE_1(model_1, model_2, MFe, 3, eps);
   COMPARE_1(model_1, model_2, MFv, 3, eps);

   COMPARE_2(model_1, model_2, Vd, 3, 3, eps);
   COMPARE_2(model_1, model_2, Ud, 3, 3, eps);
   COMPARE_2(model_1, model_2, Vu, 3, 3, eps);
   COMPARE_2(model_1, model_2, Uu, 3, 3, eps);
   COMPARE_2(model_1, model_2, Ve, 3, 3, eps);
   COMPARE_2(model_1, model_2, Ue, 3, 3, eps);
   COMPARE_2(model_1, model_2, ZZ, 2, 2, eps);

   // pole masses
   const auto pole_2 = model_2->get_physical();

   COMPARE_POLE_PTR_0(pole_2, model_2, MVP, eps);
   COMPARE_POLE_PTR_0(pole_2, model_2, MVZ, eps);
   COMPARE_POLE_PTR_0(pole_2, model_2, MVWp, eps);
   COMPARE_POLE_PTR_0(pole_2, model_2, MVG, eps);
   COMPARE_POLE_PTR_0(pole_2, model_2, Mhh, eps);
   COMPARE_POLE_PTR_0(pole_2, model_2, MHp, eps);
   COMPARE_POLE_PTR_0(pole_2, model_2, MAh, eps);
   COMPARE_POLE_PTR_1(pole_2, model_2, MFu, 3, eps);
   COMPARE_POLE_PTR_1(pole_2, model_2, MFd, 3, eps);
   COMPARE_POLE_PTR_1(pole_2, model_2, MFe, 3, eps);
   COMPARE_POLE_PTR_1(pole_2, model_2, MFv, 3, eps);

   COMPARE_POLE_PTR_2(pole_2, model_2, Vd, 3, 3, eps);
   COMPARE_POLE_PTR_2(pole_2, model_2, Ud, 3, 3, eps);
   COMPARE_POLE_PTR_2(pole_2, model_2, Vu, 3, 3, eps);
   COMPARE_POLE_PTR_2(pole_2, model_2, Uu, 3, 3, eps);
   COMPARE_POLE_PTR_2(pole_2, model_2, Ve, 3, 3, eps);
   COMPARE_POLE_PTR_2(pole_2, model_2, Ue, 3, 3, eps);
   COMPARE_POLE_PTR_2(pole_2, model_2, ZZ, 2, 2, eps);
}
