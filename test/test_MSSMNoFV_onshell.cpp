
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMNoFV_onshell

#include <boost/test/unit_test.hpp>

#include "test_CMSSMNoFV.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/MSSMNoFV_onshell.hpp"

using namespace flexiblesusy;
using namespace gm2calc;

void setup(MSSMNoFV_onshell& model)
{
   const double TanBeta = 10;
   const double SignMu = 1;
   const double m0 = 125;
   const double m12 = 500;
   const double Azero = 0;
   const double sinBeta = sin(atan(TanBeta));
   const double cosBeta = cos(atan(TanBeta));
   const double vev = 245.;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double root2 = Sqrt(2.0);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(1,1) = 0.105   * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   Eigen::Matrix<double,3,3> mm0_offset(mm0);
   mm0_offset(0,0) = Sqr(100.);
   mm0_offset(1,1) = Sqr(110.);
   mm0_offset(2,2) = Sqr(120.);

   model.set_scale(1000);
   model.set_g1(Electroweak_constants::g1);
   model.set_g2(Electroweak_constants::g2);
   model.set_g3(Electroweak_constants::g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(m12);
   model.set_MassG(m12);
   model.set_MassWB(m12 - 20);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0_offset);
   model.set_mu2(mm0_offset);
   model.set_me2(mm0_offset);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(-Sqr(m0));
   model.set_TYu(Azero * Yu);
   model.set_TYd(Azero * Yd);
   model.set_TYe(Azero * Ye);
   model.set_vu(vu);
   model.set_vd(vd);
   model.set_Mu(SignMu * 150);
   model.set_BMu(10000);

   model.solve_ewsb_tree_level();
   model.calculate_DRbar_masses();
   model.copy_DRbar_masses_to_pole_masses();

   model.get_physical().MChi(0) = model.get_MChi(0) - 2;
   model.get_physical().MChi(1) = model.get_MChi(1) - 2;
   model.get_physical().MChi(2) = model.get_MChi(2) - 3;
   model.get_physical().MChi(3) = model.get_MChi(3) - 5;
   model.get_physical().MCha(0) = model.get_MCha(0) - 2;
   model.get_physical().MCha(1) = model.get_MCha(1) - 3;
   model.get_physical().MSm(0) = model.get_MSm(0) - 10;
   model.get_physical().MSm(1) = model.get_MSm(1) - 9;
   model.get_physical().MVWm = Electroweak_constants::MW;
   model.get_physical().MVZ  = Electroweak_constants::MZ;
}

BOOST_AUTO_TEST_CASE( test_DRbar_os_conversion )
{
   MSSMNoFV_onshell model;
   setup(model);

   MSSMNoFV_onshell osmodel(model);

   BOOST_TEST_MESSAGE(osmodel);

   BOOST_CHECK_GT(std::abs(osmodel.get_MCha(0) - osmodel.get_physical().MCha(0)), 0.8);
   BOOST_CHECK_GT(std::abs(osmodel.get_MCha(1) - osmodel.get_physical().MCha(1)), 2.);
   BOOST_CHECK_GT(std::abs(osmodel.get_MChi(0) - osmodel.get_physical().MChi(0)), 1.);
   BOOST_CHECK_GT(std::abs(osmodel.get_MSm(0) - osmodel.get_physical().MSm(0)), 5.);
   BOOST_CHECK_GT(std::abs(osmodel.get_MSm(1) - osmodel.get_physical().MSm(1)), 5.);

   BOOST_TEST_MESSAGE("MChi_pole = " << osmodel.get_physical().MChi.transpose());
   BOOST_TEST_MESSAGE("MChi = " << osmodel.get_MChi().transpose());
   BOOST_TEST_MESSAGE("MSm_pole = " << osmodel.get_physical().MSm.transpose());
   BOOST_TEST_MESSAGE("MSm = " << osmodel.get_MSm().transpose());

   try {
      osmodel.set_verbose_output(true);
      osmodel.convert_to_onshell();
   } catch (const gm2calc::Error& e) {
      BOOST_TEST_MESSAGE(osmodel);
      BOOST_FAIL(e.what());
   }

   BOOST_TEST_MESSAGE("MChi = " << osmodel.get_MChi().transpose());
   BOOST_TEST_MESSAGE("MSm = " << osmodel.get_MSm().transpose());

   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MCha(0), osmodel.get_physical().MCha(0), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MCha(1), osmodel.get_physical().MCha(1), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MChi(0), osmodel.get_physical().MChi(0), 1e-4);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MSm(0), osmodel.get_physical().MSm(0), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MSvmL(), osmodel.get_physical().MSveL, 1e-5);
}
