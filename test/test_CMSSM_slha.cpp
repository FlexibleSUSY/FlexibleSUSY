
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_slha

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_model_slha.hpp"
#include "wrappers.hpp"
#include "ckm.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_two_scale_slha_cctor )
{
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4 * Pi * alpha1);
   const double g2 = Sqrt(4 * Pi * alpha2);
   const double g3 = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   CMSSM<Two_scale> model(input);
   model.set_scale(91);
   model.set_loops(1);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   model.calculate_spectrum();

   BOOST_REQUIRE(!model.get_problems().have_problem());

   // check that model is in non-SLHA convention
   BOOST_CHECK_GT(model.get_physical().MChi.maxCoeff(), 0.);
   BOOST_CHECK_GT(model.get_physical().MChi.minCoeff(), 0.);
   BOOST_CHECK_GT(model.get_physical().ZN.imag().maxCoeff(), 0.);

   // fill SLHA wrapper class
   // automatic conversion to SLHA happens here
   CMSSM_slha slha_model(model);

   // check that model wrapper is in SLHA convention
   BOOST_CHECK_GT(slha_model.get_physical_slha().MChi.maxCoeff(), 0.);
   BOOST_CHECK_LT(slha_model.get_physical_slha().MChi.minCoeff(), 0.);
   BOOST_CHECK_EQUAL(slha_model.get_physical_slha().ZN.imag().maxCoeff(), 0.);

   // no automatic conversion
   CMSSM_slha slha_model_not_converted(model, false);

   // check that model is in non-SLHA convention
   BOOST_CHECK_GT(slha_model_not_converted.get_physical().MChi.maxCoeff(), 0.);
   BOOST_CHECK_GT(slha_model_not_converted.get_physical().MChi.minCoeff(), 0.);
   BOOST_CHECK_GT(slha_model_not_converted.get_physical().ZN.imag().maxCoeff(), 0.);
}

BOOST_AUTO_TEST_CASE( test_CMSSM_two_scale_slha_calculate_spectrum )
{
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4 * Pi * alpha1);
   const double g2 = Sqrt(4 * Pi * alpha2);
   const double g3 = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   CMSSM_slha model(input);
   model.set_scale(91);
   model.set_loops(1);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   model.calculate_spectrum(); // conversion happens here

   BOOST_REQUIRE(!model.get_problems().have_problem());

   // check that model class is in non-SLHA convention
   BOOST_CHECK_GT(model.get_physical().MChi.maxCoeff(), 0.);
   BOOST_CHECK_GT(model.get_physical().MChi.minCoeff(), 0.);
   BOOST_CHECK_GT(model.get_physical().ZN.imag().maxCoeff(), 0.);

   // check that model wrapper class is in SLHA convention
   BOOST_CHECK_GT(model.get_physical_slha().MChi.maxCoeff(), 0.);
   BOOST_CHECK_LT(model.get_physical_slha().MChi.minCoeff(), 0.);
   BOOST_CHECK_EQUAL(model.get_physical_slha().ZN.imag().maxCoeff(), 0.);
}

BOOST_AUTO_TEST_CASE( test_CMSSM_two_scale_slha_diagonal_yukawas )
{
   CMSSM_input_parameters input;
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4 * Pi * alpha1);
   const double g2 = Sqrt(4 * Pi * alpha2);
   const double g3 = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         Yu(i,k) = 0.1 * (i+1) * (k+1);
         Yd(i,k) = 0.1 * (i+1) * (k+1);
         Ye(i,k) = 0.1 * (i+1) * (k+1);
      }
   }
   Yu <<  0.1, 0.2, 0.3,
         -0.2, 0.4, 0.6,
          0.3, 0.6, 165.0   * root2 / (vev * sinBeta);
   Yd <<  0.1, 0.2, 0.3,
         -0.2, 0.4, 0.6,
          0.3, 0.6, 2.9     * root2 / (vev * sinBeta);
   Ye <<  0.1, 0.2, 0.3,
         -0.2, 0.4, 0.6,
          0.3, 0.6, 1.77699 * root2 / (vev * sinBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   CMSSM_slha model(input);
   model.set_scale(91);
   model.set_loops(1);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0 * Yu);
   model.set_ml2(mm0 * Ye);
   model.set_md2(mm0 * Yd);
   model.set_mu2(mm0 * Yu);
   model.set_me2(mm0 * Ye);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   model.solve_ewsb_tree_level();
   model.calculate_spectrum(); // conversion happens here

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK(model.get_Yu(i,k) != 0.);
         BOOST_CHECK(model.get_Yd(i,k) != 0.);
         BOOST_CHECK(model.get_Ye(i,k) != 0.);

         BOOST_CHECK(model.get_TYu(i,k) != 0.);
         BOOST_CHECK(model.get_TYd(i,k) != 0.);
         BOOST_CHECK(model.get_TYe(i,k) != 0.);
      }
   }

   // check that SLHA-compliant Yukawa couplings are diagonal
   for (int i = 0; i < 3; i++) {
      BOOST_CHECK_GT(model.get_Yu_slha(i), 0.);
      BOOST_CHECK_GT(model.get_Yd_slha(i), 0.);
      BOOST_CHECK_GT(model.get_Ye_slha(i), 0.);
   }

   // check SLHA-compliant trilinear couplings
   const Eigen::Matrix<std::complex<double>,3,3> ZDL_slha(model.get_ZDL_slha());
   const Eigen::Matrix<std::complex<double>,3,3> ZUL_slha(model.get_ZUL_slha());
   const Eigen::Matrix<std::complex<double>,3,3> ZDR_slha(model.get_ZDR_slha());
   const Eigen::Matrix<std::complex<double>,3,3> ZUR_slha(model.get_ZUR_slha());
   const Eigen::Matrix<std::complex<double>,3,3> ZEL_slha(model.get_ZEL_slha());
   const Eigen::Matrix<std::complex<double>,3,3> ZER_slha(model.get_ZER_slha());

   // slha convention
   const Eigen::Matrix<double,3,3> TYu_slha(model.get_TYu_slha());
   const Eigen::Matrix<double,3,3> TYd_slha(model.get_TYd_slha());
   const Eigen::Matrix<double,3,3> TYe_slha(model.get_TYe_slha());

   // non-slha
   const Eigen::Matrix<double,3,3> TYu((ZUR_slha.transpose() * TYu_slha * ZUL_slha).real());
   const Eigen::Matrix<double,3,3> TYd((ZDR_slha.transpose() * TYd_slha * ZDL_slha).real());
   const Eigen::Matrix<double,3,3> TYe((ZER_slha.transpose() * TYe_slha * ZEL_slha).real());

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE_FRACTION(model.get_TYu(i,k), TYu(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(model.get_TYd(i,k), TYd(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(model.get_TYe(i,k), TYe(i,k), 1.0e-10);
      }
   }

   // slha convention
   const Eigen::Matrix<double,3,3> mq2_slha(model.get_mq2_slha());
   const Eigen::Matrix<double,3,3> ml2_slha(model.get_ml2_slha());
   const Eigen::Matrix<double,3,3> mu2_slha(model.get_mu2_slha());
   const Eigen::Matrix<double,3,3> md2_slha(model.get_md2_slha());
   const Eigen::Matrix<double,3,3> me2_slha(model.get_me2_slha());

   // non-slha
   const Eigen::Matrix<double,3,3> mq2((ZDL_slha.adjoint() * mq2_slha * ZDL_slha).real());
   const Eigen::Matrix<double,3,3> ml2((ZEL_slha.adjoint() * ml2_slha * ZEL_slha).real());
   const Eigen::Matrix<double,3,3> mu2((ZUR_slha.transpose() * mu2_slha * ZUR_slha.conjugate()).real());
   const Eigen::Matrix<double,3,3> md2((ZDR_slha.transpose() * md2_slha * ZDR_slha.conjugate()).real());
   const Eigen::Matrix<double,3,3> me2((ZER_slha.transpose() * me2_slha * ZER_slha.conjugate()).real());

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE_FRACTION(model.get_mq2(i,k), mq2(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(model.get_ml2(i,k), ml2(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(model.get_mu2(i,k), mu2(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(model.get_md2(i,k), md2(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(model.get_me2(i,k), me2(i,k), 1.0e-10);
      }
   }

   // check CKM matrix
   const Eigen::Matrix<std::complex<double>,3,3> ckm_matrix(ZUL_slha * ZDL_slha.adjoint());
   const Eigen::Matrix<std::complex<double>,3,3> ckm_slha(model.get_ckm_matrix());

   BOOST_TEST_MESSAGE("ckm_matrix =\n" << ckm_matrix);
   BOOST_TEST_MESSAGE("ckm_slha =\n" << ckm_slha);

   // check that SLHA mixing matrix is in PDG convention
   BOOST_CHECK(Re(ckm_slha(0,0)) > 0.);
   BOOST_CHECK(Re(ckm_slha(1,1)) > 0.);
   BOOST_CHECK(Re(ckm_slha(2,2)) > 0.);
   BOOST_CHECK(Re(ckm_slha(0,1)) > 0.);
   BOOST_CHECK(Re(ckm_slha(1,2)) > 0.);

   // check that SLHA mixing matrix is in PDG convention
   BOOST_CHECK(Re(ckm_matrix(0,1)) > 0.);
   BOOST_CHECK(Re(ckm_matrix(1,2)) > 0.);
   BOOST_CHECK(Re(ckm_matrix(0,0)) > 0.);
   BOOST_CHECK(Re(ckm_matrix(1,1)) > 0.);
   BOOST_CHECK(Re(ckm_matrix(2,2)) > 0.);

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE_FRACTION(Re(ckm_matrix(i,k)), Re(ckm_slha(i,k)), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(Im(ckm_matrix(i,k)), Im(ckm_slha(i,k)), 1.0e-10);
      }
   }
}
