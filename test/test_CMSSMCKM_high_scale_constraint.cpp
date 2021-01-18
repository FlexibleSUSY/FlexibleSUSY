
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCKM_high_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_legacy.hpp"
#include <functional>
#include <Eigen/Dense>

#include "flavoursoft.h"
#include "CMSSMCKM_two_scale_model.hpp"
#include "CMSSMCKM_two_scale_high_scale_constraint.hpp"
#include "CMSSMCKM_input_parameters.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "conversion.hpp"
#include "def.h"

using namespace flexiblesusy;
using namespace softsusy;

void setup(CMSSMCKM_input_parameters& input, DoubleVector& input2)
{
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double vev = 246.0;
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> mm0;
   mm0 << Sqr(130), 200     , 100,
          200     , Sqr(170), 300,
          100     , 300     , Sqr(200);

   input.MassBInput = M12;
   input.MassGInput = M12;
   input.MassWBInput = M12;

   input.mq2Input = mm0;
   input.ml2Input = mm0;
   input.md2Input = mm0;
   input.mu2Input = mm0;
   input.me2Input = mm0;

   input.mHd2IN = Sqr(125);
   input.mHu2IN = Sqr(150);

   input.TYuInput << 100, 2  , 3,
                     4  , 500, 6,
                     7  , 8  , 900;
   input.TYdInput = input.TYuInput;
   input.TYeInput = input.TYuInput;

   // fill Softsusy input parameter vector
   int count = 1;

   // gauginos
   input2(count++) = input.MassBInput;
   input2(count++) = input.MassWBInput;
   input2(count++) = input.MassGInput;

   // soft squared masses
   for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++)
         input2(count++) = input.mq2Input(i-1,j-1);

   for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++)
         input2(count++) = input.mu2Input(i-1,j-1);

   for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++)
         input2(count++) = input.md2Input(i-1,j-1);

   for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++)
         input2(count++) = input.ml2Input(i-1,j-1);

   for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++)
         input2(count++) = input.me2Input(i-1,j-1);

  for (int i=1; i<=3; i++)
    for (int j=1; j<=3; j++)
       input2(count++) = input.TYuInput(i-1,j-1);

  for (int i=1; i<=3; i++)
    for (int j=1; j<=3; j++)
       input2(count++) = input.TYdInput(i-1,j-1);

  for (int i=1; i<=3; i++)
    for (int j=1; j<=3; j++)
       input2(count++) = input.TYeInput(i-1,j-1);

  input2(63) = input.mHd2IN;
  input2(64) = input.mHu2IN;
}

void setup(CMSSMCKM<Two_scale>& m, FlavourMssmSoftsusy& s)
{
   Eigen::Matrix<std::complex<double>,3,3>
      Yu(Eigen::Matrix<std::complex<double>,3,3>::Zero()),
      Yd(Eigen::Matrix<std::complex<double>,3,3>::Zero()),
      Ye(Eigen::Matrix<std::complex<double>,3,3>::Zero());

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double vev = 246.0;
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Yu << 1, 0.2, 0.1,
         0, 2  , 0.2,
         0, 0  , 3;

   Yd << 4, 0.3, 0.2,
         0, 5  , 0.3,
         0, 0  , 6;

   Ye << 7, 0.4, 0.1,
         0, 8  , 0.4,
         0, 0  , 9;

   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);

   s.setYukawaMatrix(YU, ToDoubleMatrix(Yu.transpose().real()));
   s.setYukawaMatrix(YD, ToDoubleMatrix(Yd.transpose().real()));
   s.setYukawaMatrix(YE, ToDoubleMatrix(Ye.transpose().real()));
   s.setTanb(tanBeta);
   s.setHvev(vev);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);
}

BOOST_AUTO_TEST_CASE( test_high_scale_constraint )
{
   CMSSMCKM<Two_scale> m;
   FlavourMssmSoftsusy s;
   CMSSMCKM_input_parameters input;
   DoubleVector inputParameters(64);
   setup(input, inputParameters);
   setup(m, s);

   m.set_input_parameters(input);
   CMSSMCKM_high_scale_constraint<Two_scale> constraint(&m);

   try {
      m.calculate_DRbar_masses();
      s.calcDrBarPars();
   } catch (const std::string& str) {
      BOOST_FAIL("exception thrown: " << str);
   } catch (...) {
      BOOST_FAIL("exception thrown");
   }

   // precondition: all parameters are equal
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHd2(), s.displayMh1Squared(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHu2(), s.displayMh2Squared(), 1.0e-10);

   TEST_CLOSE(ToDoubleMatrix(m.get_mq2().real()), s.displaySoftMassSquared(mQl), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_ml2().real()), s.displaySoftMassSquared(mLl), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_md2().real()), s.displaySoftMassSquared(mDr), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_mu2().real()), s.displaySoftMassSquared(mUr), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_me2().real()), s.displaySoftMassSquared(mEr), 1.0e-6);

   BOOST_CHECK_CLOSE_FRACTION(m.get_MassB() , s.displayGaugino(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassG() , s.displayGaugino(2), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassWB(), s.displayGaugino(3), 1.0e-10);

   TEST_CLOSE(ToDoubleMatrix(m.get_TYu().real()), s.displayTrilinear(UA).transpose(), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_TYd().real()), s.displayTrilinear(DA).transpose(), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_TYe().real()), s.displayTrilinear(EA).transpose(), 1.0e-6);

   for (int i = 0; i < 27; i++)
      slha2setTrilinear[i] = true;

   constraint.apply();
   flavourBcs(s, inputParameters);

   // postcondition: all parameters are equal
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHd2(), s.displayMh1Squared(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHu2(), s.displayMh2Squared(), 1.0e-10);

   TEST_CLOSE(ToDoubleMatrix(m.get_mq2().real()), s.displaySoftMassSquared(mQl), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_ml2().real()), s.displaySoftMassSquared(mLl), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_md2().real()), s.displaySoftMassSquared(mDr), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_mu2().real()), s.displaySoftMassSquared(mUr), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_me2().real()), s.displaySoftMassSquared(mEr), 1.0e-6);

   BOOST_CHECK_CLOSE_FRACTION(m.get_MassB() , s.displayGaugino(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassG() , s.displayGaugino(2), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassWB(), s.displayGaugino(3), 1.0e-10);

   TEST_CLOSE(ToDoubleMatrix(m.get_TYu().real()), s.displayTrilinear(UA).transpose(), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_TYd().real()), s.displayTrilinear(DA).transpose(), 1.0e-6);
   TEST_CLOSE(ToDoubleMatrix(m.get_TYe().real()), s.displayTrilinear(EA).transpose(), 1.0e-6);

   BOOST_CHECK_EQUAL(gErrors, 0);
}
