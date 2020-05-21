
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCKM_tree_level_spectrum

#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>

#include "wrappers.hpp"
#include "conversion.hpp"
#include "flavoursoft.h"
#include "CMSSMCKM_two_scale_model.hpp"
#include "CMSSMCKM_model_slha.hpp"
#include "CMSSMCKM_two_scale_high_scale_constraint.hpp"
#include "CMSSMCKM_input_parameters.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void setup(CMSSMCKM_input_parameters& input, DoubleVector& input2)
{
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> mm0;
   mm0 << Sqr(330), 200     , 100,
          200     , Sqr(370), 300,
          100     , 300     , Sqr(350);

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
                     7  , 8  , 700;
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

void setup(CMSSMCKM_mass_eigenstates& m, FlavourMssmSoftsusy& s)
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
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Yu << 0.1, 0.1, 0.1,
         0  , 0.3, 0.2,
         0  , 0  , 0.9;

   Yd << 0.4, 0.3, 0.2,
         0  , 0.5, 0.3,
         0  , 0  , 0.6;

   Ye << 0.7, 0.4, 0.1,
         0  , 0.8, 0.4,
         0  , 0  , 0.9;

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

   // CHECK: Quark Yukawa couplings are transposed
   s.setYukawaMatrix(YU, ToDoubleMatrix(Yu.transpose().real()));
   s.setYukawaMatrix(YD, ToDoubleMatrix(Yd.transpose().real()));
   s.setYukawaMatrix(YE, ToDoubleMatrix(Ye.real()));
   s.setTanb(tanBeta);
   s.setHvev(vev);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);

   // fill soft-breaking parameters
   CMSSMCKM_input_parameters input;
   DoubleVector inputParameters(64);
   setup(input, inputParameters);

   m.set_TYe(input.TYeInput.cast<std::complex<double> >());
   m.set_TYd(input.TYdInput.cast<std::complex<double> >());
   m.set_TYu(input.TYuInput.cast<std::complex<double> >());
   m.set_mHd2(input.mHd2IN);
   m.set_mHu2(input.mHu2IN);
   m.set_mq2(input.mq2Input.cast<std::complex<double> >());
   m.set_ml2(input.ml2Input.cast<std::complex<double> >());
   m.set_md2(input.md2Input.cast<std::complex<double> >());
   m.set_mu2(input.mu2Input.cast<std::complex<double> >());
   m.set_me2(input.me2Input.cast<std::complex<double> >());
   m.set_MassB(input.MassBInput);
   m.set_MassWB(input.MassWBInput);
   m.set_MassG(input.MassGInput);

   // CHECK: trilinear couplings are transposed
   s.setGauginoMass(1, input.MassBInput);
   s.setGauginoMass(2, input.MassWBInput);
   s.setGauginoMass(3, input.MassGInput);
   s.setSoftMassMatrix(mQl, ToDoubleMatrix(input.mq2Input));
   s.setSoftMassMatrix(mLl, ToDoubleMatrix(input.ml2Input));
   s.setSoftMassMatrix(mDr, ToDoubleMatrix(input.md2Input));
   s.setSoftMassMatrix(mUr, ToDoubleMatrix(input.mu2Input));
   s.setSoftMassMatrix(mEr, ToDoubleMatrix(input.me2Input));
   s.setMh1Squared(input.mHd2IN);
   s.setMh2Squared(input.mHu2IN);
   s.setTrilinearMatrix(EA, ToDoubleMatrix(input.TYeInput.transpose()));
   s.setTrilinearMatrix(DA, ToDoubleMatrix(input.TYdInput.transpose()));
   s.setTrilinearMatrix(UA, ToDoubleMatrix(input.TYuInput.transpose()));
   s.setMw(s.displayMwRun());
}

DoubleMatrix get_CKM(const FlavourMssmSoftsusy& s)
{
   DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3);
   DoubleVector yu(3), yd(3);

   s.displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
   s.displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
   ckmNormalise(Vu, Vd, Uu, Ud);

   DoubleMatrix CKM(Vu.transpose() * Vd);

   return CKM;
}

BOOST_AUTO_TEST_CASE( test_CMSSMCKM_tree_level_masses )
{
   // This test compares the DR-bar spectrum of Softsusy (in the
   // super-CKM basis) with the DR-bar masses of CMSSMCKM.

   softsusy::PRINTOUT = 2;
   FlavourMssmSoftsusy s;
   CMSSMCKM<Two_scale> m0;
   m0.do_force_output(true);
   setup(m0, s);

   m0.calculate_DRbar_masses();
   s.calcDrBarPars();

   CMSSMCKM_slha m(m0); // converts to SLHA-2

   // re-set model parameters to super-CKM basis
   m.set_Yu(m.get_Yu_slha().matrix().cast<std::complex<double> >().asDiagonal());
   m.set_Yd(m.get_Yd_slha().matrix().cast<std::complex<double> >().asDiagonal());
   m.set_Ye(m.get_Ye_slha().matrix().cast<std::complex<double> >().asDiagonal());

   m.set_mq2(m.get_mq2_slha().cast<std::complex<double> >());
   m.set_ml2(m.get_ml2_slha().cast<std::complex<double> >());
   m.set_md2(m.get_md2_slha().cast<std::complex<double> >());
   m.set_mu2(m.get_mu2_slha().cast<std::complex<double> >());
   m.set_me2(m.get_me2_slha().cast<std::complex<double> >());

   m.set_TYu(m.get_TYu_slha().cast<std::complex<double> >());
   m.set_TYd(m.get_TYd_slha().cast<std::complex<double> >());
   m.set_TYe(m.get_TYe_slha().cast<std::complex<double> >());

   BOOST_TEST_MESSAGE("super-CKM: Yu = " << m.get_Yu());
   BOOST_TEST_MESSAGE("super-CKM: Yd = " << m.get_Yd());
   BOOST_TEST_MESSAGE("super-CKM: Ye = " << m.get_Ye());
   BOOST_TEST_MESSAGE("super-CKM: mq2 = " << m.get_mq2());
   BOOST_TEST_MESSAGE("super-CKM: mu2 = " << m.get_mu2());
   BOOST_TEST_MESSAGE("super-CKM: md2 = " << m.get_md2());
   BOOST_TEST_MESSAGE("super-CKM: ml2 = " << m.get_ml2());
   BOOST_TEST_MESSAGE("super-CKM: me2 = " << m.get_me2());
   BOOST_TEST_MESSAGE("super-CKM: TYu = " << m.get_TYu());
   BOOST_TEST_MESSAGE("super-CKM: TYd = " << m.get_TYd());
   BOOST_TEST_MESSAGE("super-CKM: TYe = " << m.get_TYe());

   m.calculate_DRbar_masses();

   // Now CMSSMCKM down-type sfermion masses should be in super-CKM
   // basis.

   // neutral CP even Higgs
   const DoubleVector hh(ToDoubleVector(m.get_Mhh()));
   const DoubleVector mh0(s.displayDrBarPars().mh0);
   BOOST_CHECK_CLOSE(hh(1), mh0(1), 1.0e-12);
   BOOST_CHECK_CLOSE(hh(2), mh0(2), 1.0e-12);

   // neutral CP odd Higgs
   const DoubleVector Ah(ToDoubleVector(m.get_MAh()));
   const DoubleVector mA0(s.displayDrBarPars().mA0);
   const double mz = s.displayMzRun();
   BOOST_CHECK_CLOSE(Ah(1), mz    , 1.0e-12);
   BOOST_CHECK_CLOSE(Ah(2), mA0(1), 1.0e-12);

   // charged Higgs
   const DoubleVector Hpm(ToDoubleVector(m.get_MHpm()));
   const double mHpm = s.displayDrBarPars().mHpm;
   const double mw = s.displayMwRun();
   BOOST_CHECK_CLOSE(Hpm(1), mw  , 1.0e-12);
   BOOST_CHECK_CLOSE(Hpm(2), mHpm, 1.0e-12);

   // neutralinos
   const DoubleVector mneut(ToDoubleVector(m.get_MChi()));
   const DoubleVector mnBpmz = s.displayDrBarPars().mnBpmz;
   BOOST_CHECK_CLOSE(mneut(1), mnBpmz(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(2), mnBpmz(2), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(3), mnBpmz(3), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(4), mnBpmz(4), 1.0e-12);

   // charginos
   const DoubleVector mch(ToDoubleVector(m.get_MCha()));
   const DoubleVector mchBpmz = s.displayDrBarPars().mchBpmz;
   BOOST_CHECK_CLOSE(mch(1), mchBpmz(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mch(2), mchBpmz(2), 1.0e-12);

   // photon, W, Z, gluon mass
   const double vp = m.get_MVP();
   const double vz = m.get_MVZ();
   const double vw = m.get_MVWm();
   const double vg = m.get_MVG();
   BOOST_CHECK_EQUAL(vp, 0.0);
   BOOST_CHECK_CLOSE(vz, mz, 1.0e-12);
   BOOST_CHECK_CLOSE(vw, mw, 1.0e-12);
   BOOST_CHECK_EQUAL(vg, 0.0);

   // down-type squarks
   const DoubleVector Sd(ToDoubleVector(m.get_MSd()));
   const DoubleVector md(s.displayDrBarPars().md.flatten().sort());
   BOOST_CHECK_CLOSE(Sd(1), md(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(2), md(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(3), md(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(4), md(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(5), md(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(6), md(6), 1.0e-12);

   // down-type sleptons
   const DoubleVector Se(ToDoubleVector(m.get_MSe()));
   const DoubleVector me(s.displayDrBarPars().me.flatten().sort());
   BOOST_CHECK_CLOSE(Se(1), me(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(2), me(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(3), me(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(4), me(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(5), me(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(6), me(6), 1.0e-12);

   // gluinos
   const double MGlu = m.get_MGlu();
   const double mGluino = s.displayDrBarPars().mGluino;
   BOOST_CHECK_CLOSE(MGlu, mGluino, 1.0e-12);

   // neutrinos
   const DoubleVector MFv(ToDoubleVector(m.get_MFv()));
   BOOST_CHECK_EQUAL(MFv(1), 0.0);
   BOOST_CHECK_EQUAL(MFv(2), 0.0);
   BOOST_CHECK_EQUAL(MFv(3), 0.0);

   // leptons
   const DoubleVector MFe(ToDoubleVector(m.get_MFe()));
   BOOST_CHECK_CLOSE(MFe(3), s.displayDrBarPars().mtau, 1.0e-12);

   // ups
   const DoubleVector MFu(ToDoubleVector(m.get_MFu()));
   BOOST_CHECK_CLOSE(MFu(3), s.displayDrBarPars().mt, 1.0e-12);

   // downs
   const DoubleVector MFd(ToDoubleVector(m.get_MFd()));
   BOOST_CHECK_CLOSE(MFd(3), s.displayDrBarPars().mb, 1.0e-12);

   // compare CKM matrices

   const auto VCKM = m.get_ckm_matrix();
   const auto VCKM_SS = get_CKM(s);

   BOOST_TEST_MESSAGE("CKM(FS): " << VCKM);
   BOOST_TEST_MESSAGE("CKM(SS): " << VCKM_SS);

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE(std::real(VCKM(i,k)), std::real(VCKM_SS(i+1,k+1)), 1e-10);
         BOOST_CHECK_SMALL(std::imag(VCKM(i,k)), 1e-7);
         BOOST_CHECK_SMALL(std::imag(VCKM_SS(i+1,k+1)), 1e-10);
      }
   }

   // Now set CMSSMCKM up-type sfermion masses to super-CKM basis, as
   // in Eq. (11) of SLHA-2.

   const auto mq2_super_CKM = m.get_mq2();
   const auto mq2_up = VCKM * mq2_super_CKM * VCKM.adjoint();

   m.set_mq2(mq2_up);

   BOOST_TEST_MESSAGE("SLHA-2: mq2 = " << m.get_mq2());
   BOOST_TEST_MESSAGE("SLHA-2: mu2 = " << m.get_mu2());

   m.calculate_DRbar_masses();

   // up-type squarks
   const DoubleVector Su(ToDoubleVector(m.get_MSu()));
   const DoubleVector mu(s.displayDrBarPars().mu.flatten().sort());
   BOOST_CHECK_CLOSE(Su(1), mu(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(2), mu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(3), mu(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(4), mu(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(5), mu(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(6), mu(6), 1.0e-12);
}
