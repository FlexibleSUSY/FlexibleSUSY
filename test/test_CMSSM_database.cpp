#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_database

#include <boost/test/unit_test.hpp>

#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_input_parameters.hpp"
#include "CMSSM_utilities.hpp"
#include "ew_input.hpp"
#include "CMSSM_observables.hpp"
#include "physical_input.hpp"
#include "test.hpp"
#include "lowe.h"
#include "wrappers.hpp"

using namespace flexiblesusy;

void setup_CMSSM(CMSSM_mass_eigenstates& m, const CMSSM_input_parameters& input)
{
   constexpr double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   constexpr double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   constexpr double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(2);
   m.set_thresholds(3);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   m.do_force_output(true);

   m.calculate_DRbar_masses();
   m.solve_ewsb();
   m.calculate_spectrum();
}

template <class T, class P>
void test_parameter_equality(const T& a, const P& b, double eps)
{
   TEST_CLOSE(a.get_loops(), b.get_loops(), eps);
   TEST_CLOSE(a.get_scale(), b.get_scale(), eps);
   TEST_CLOSE(a.get_thresholds(), b.get_thresholds(), eps);
   TEST_CLOSE(a.get_ewsb_loop_order(), b.get_ewsb_loop_order(), eps);
   TEST_CLOSE(a.get_pole_mass_loop_order(), b.get_pole_mass_loop_order(), eps);
   TEST_CLOSE(a.get_precision(), b.get_precision(), eps);
   TEST_CLOSE(a.get_ewsb_iteration_precision(), b.get_ewsb_iteration_precision(), eps);

   TEST_CLOSE(a.get_number_of_parameters(), b.get_number_of_parameters(), eps);

   TEST_CLOSE(a.get_g1(), b.get_g1(), eps);
   TEST_CLOSE(a.get_g2(), b.get_g2(), eps);
   TEST_CLOSE(a.get_g3(), b.get_g3(), eps);

   TEST_CLOSE(a.get_Yu(), b.get_Yu(), eps);
   TEST_CLOSE(a.get_Yd(), b.get_Yd(), eps);
   TEST_CLOSE(a.get_Ye(), b.get_Ye(), eps);

   TEST_CLOSE(a.get_MassB(), b.get_MassB(), eps);
   TEST_CLOSE(a.get_MassWB(), b.get_MassWB(), eps);
   TEST_CLOSE(a.get_MassG(), b.get_MassG(), eps);

   TEST_CLOSE(a.get_mHd2(), b.get_mHd2(), eps);
   TEST_CLOSE(a.get_mHu2(), b.get_mHu2(), eps);
   TEST_CLOSE(a.get_mq2(), b.get_mq2(), eps);
   TEST_CLOSE(a.get_mu2(), b.get_mu2(), eps);
   TEST_CLOSE(a.get_md2(), b.get_md2(), eps);
   TEST_CLOSE(a.get_ml2(), b.get_ml2(), eps);
   TEST_CLOSE(a.get_me2(), b.get_me2(), eps);

   TEST_CLOSE(a.get_TYu(), b.get_TYu(), eps);
   TEST_CLOSE(a.get_TYd(), b.get_TYd(), eps);
   TEST_CLOSE(a.get_TYe(), b.get_TYe(), eps);

   TEST_CLOSE(a.get_Mu(), b.get_Mu(), eps);
   TEST_CLOSE(a.get_BMu(), b.get_BMu(), 1e-10);

   TEST_CLOSE(a.get_vu(), b.get_vu(), eps);
   TEST_CLOSE(a.get_vd(), b.get_vd(), eps);
}

template <class T, class P>
void test_mass_equality(const T& a, const P& b, double eps)
{
   TEST_CLOSE(a.Mhh, b.Mhh, eps);
   TEST_CLOSE(a.MAh, b.MAh, eps);
   TEST_CLOSE(a.MHpm, b.MHpm, eps);
   TEST_CLOSE(a.MCha, b.MCha, eps);
   TEST_CLOSE(a.MChi, b.MChi, eps);
   TEST_CLOSE(a.MGlu, b.MGlu, eps);
   TEST_CLOSE(a.MSu, b.MSu, eps);
   TEST_CLOSE(a.MSd, b.MSd, eps);
   TEST_CLOSE(a.MSe, b.MSe, eps);
   TEST_CLOSE(a.MSv, b.MSv, eps);
   TEST_CLOSE(a.MFu, b.MFu, eps);
   TEST_CLOSE(a.MFd, b.MFd, eps);
   TEST_CLOSE(a.MFe, b.MFe, eps);
   TEST_CLOSE(a.MFv, b.MFv, eps);
   TEST_CLOSE(a.MVP, b.MVP, eps);
   TEST_CLOSE(a.MVG, b.MVG, eps);
   TEST_CLOSE(a.MVWm, b.MVWm, eps);
   TEST_CLOSE(a.MVZ, b.MVZ, eps);

   TEST_CLOSE(a.ZH, b.ZH, eps);
   TEST_CLOSE(a.ZA, b.ZA, eps);
   TEST_CLOSE(a.ZP, b.ZP, eps);
   TEST_CLOSE(a.ZN, b.ZN, eps);
   TEST_CLOSE(a.UM, b.UM, eps);
   TEST_CLOSE(a.UP, b.UP, eps);
   TEST_CLOSE(a.ZU, b.ZU, eps);
   TEST_CLOSE(a.ZD, b.ZD, eps);
   TEST_CLOSE(a.ZE, b.ZE, eps);
   TEST_CLOSE(a.ZV, b.ZV, eps);
   TEST_CLOSE(a.ZEL, b.ZEL, eps);
   TEST_CLOSE(a.ZER, b.ZER, eps);
   TEST_CLOSE(a.ZUL, b.ZUL, eps);
   TEST_CLOSE(a.ZUR, b.ZUR, eps);
   TEST_CLOSE(a.ZDL, b.ZDL, eps);
   TEST_CLOSE(a.ZDR, b.ZDR, eps);
}

template <class T, class P>
void test_input_parameter_equality(const T& a, const P& b, double eps)
{
   TEST_CLOSE(a.m0, b.m0, eps);
   TEST_CLOSE(a.m12, b.m12, eps);
   TEST_CLOSE(a.TanBeta, b.TanBeta, eps);
   TEST_CLOSE(a.SignMu, b.SignMu, eps);
   TEST_CLOSE(a.Azero, b.Azero, eps);
}

BOOST_AUTO_TEST_CASE( test_CMSSM_read_write )
{
   CMSSM_observables obs1, obs2;
   obs1.a_muon = 2e-9;
   Physical_input physical_input1, physical_input2;
   physical_input1.set(Physical_input::alpha_em_0, 0.1);
   softsusy::QedQcd qedqcd1, qedqcd2;
   qedqcd1.setAlpha(softsusy::ALPHA, 1./127.);
   qedqcd1.toMz();

   CMSSM_input_parameters pp;
   pp.m0 = 125.;
   pp.m12 = 500.;
   pp.TanBeta = 10.;
   pp.SignMu = 1;
   pp.Azero = 0.1;

   CMSSM_mass_eigenstates model(pp);
   setup_CMSSM(model, pp);

   const std::string db_file("test/test_CMSSM_database.db");

   BOOST_TEST_MESSAGE("writing to database " << db_file);
   CMSSM_database::to_database(db_file, model, &qedqcd1, &physical_input1, &obs1);

   BOOST_TEST_MESSAGE("reading from database " << db_file);
   const CMSSM_mass_eigenstates tmp(CMSSM_database::from_database(db_file, -1, &qedqcd2, &physical_input2, &obs2));

   constexpr double eps = 1e-10;

   gErrors = 0;
   test_parameter_equality(model, tmp, eps);
   test_mass_equality(model.get_physical(), tmp.get_physical(), eps);
   test_input_parameter_equality(model.get_input(), tmp.get_input(), eps);
   BOOST_REQUIRE(gErrors == 0);
   BOOST_REQUIRE((qedqcd1.display_input() - qedqcd2.display_input()).cwiseAbs().maxCoeff() < 1e-10);
   BOOST_REQUIRE((physical_input1.get() - physical_input2.get()).cwiseAbs().maxCoeff() < 1e-10);
   BOOST_REQUIRE((obs1.get() - obs2.get()).cwiseAbs().maxCoeff() < 1e-10);
}
