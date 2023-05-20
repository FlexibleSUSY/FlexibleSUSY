
#ifndef TEST_CMSSMCPVSEMIANALYTIC_H
#define TEST_CMSSMCPVSEMIANALYTIC_H

#include "CMSSMCPVSemiAnalytic_input_parameters.hpp"
#include "CMSSMCPVSemiAnalytic_mass_eigenstates.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

struct Boundary_values {
   double m12{};
   double Imm12{};
   double Azero{};
   double ImAzero{};
   double m0Sq{};
   double ReBMu0{};
   double ImBMu0{};
   std::complex<double> Mu;
};

void setup_high_scale_CMSSMCPVSemiAnalytic_const(
   flexiblesusy::CMSSMCPVSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   using namespace flexiblesusy;

   model.get_input().m12 = values.m12;
   model.get_input().Imm12 = values.Imm12;
   model.get_input().Azero = values.Azero;
   model.get_input().ImAzero = values.ImAzero;

   model.set_TYu(model.get_Yu() * values.Azero);
   model.set_TYd(model.get_Yd() * values.Azero);
   model.set_TYe(model.get_Ye() * values.Azero);

   model.set_mq2((values.m0Sq * UNITMATRIX(3)).cast<std::complex<double> >());
   model.set_mu2((values.m0Sq * UNITMATRIX(3)).cast<std::complex<double> >());
   model.set_md2((values.m0Sq * UNITMATRIX(3)).cast<std::complex<double> >());
   model.set_ml2((values.m0Sq * UNITMATRIX(3)).cast<std::complex<double> >());
   model.set_me2((values.m0Sq * UNITMATRIX(3)).cast<std::complex<double> >());
   model.set_mHd2(Re(values.m0Sq));
   model.set_mHu2(Re(values.m0Sq));

   model.set_MassB(values.m12 + std::complex<double>(0.,1.) * values.Imm12);
   model.set_MassWB(values.m12 + std::complex<double>(0.,1.) * values.Imm12);
   model.set_MassG(values.m12 + std::complex<double>(0.,1.) * values.Imm12);

   model.set_Mu(values.Mu);
   model.set_BMu(values.ReBMu0 + std::complex<double>(0.,1.) * values.ImBMu0);

   model.set_m0Sq(values.m0Sq);
   model.set_ReBMu0(values.ReBMu0);
   model.set_ImBMu0(values.ImBMu0);
   model.set_MuBV(values.Mu);
}

void setup_high_scale_CMSSMCPVSemiAnalytic(
   flexiblesusy::CMSSMCPVSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   using namespace flexiblesusy;

   values.m12 = 300.0;
   values.Imm12 = 0.;
   values.Azero = -100.;
   values.ImAzero = 0.;
   values.m0Sq = Sqr(1000.);
   values.ReBMu0 = 1.0e4;
   values.ImBMu0 = 0.;
   values.Mu = model.get_Mu();

   setup_high_scale_CMSSMCPVSemiAnalytic_const(model, values);
}

void setup_CMSSMCPVSemiAnalytic_const(
   flexiblesusy::CMSSMCPVSemiAnalytic_mass_eigenstates& m,
   const flexiblesusy::CMSSMCPVSemiAnalytic_input_parameters& input)
{
   using namespace flexiblesusy;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5.0 * ALPHAMZ / (3.0 * (1.0 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4.0 * Pi * alpha1);
   const double g2 = Sqrt(4.0 * Pi * alpha2);
   const double g3 = Sqrt(4.0 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = Sin(ArcTan(tanBeta));
   const double cosBeta = Cos(ArcTan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m12;
   const double a0 = input.Azero;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const std::complex<double> susyMu = input.AbsMu *
      (input.CosPhiMu + std::complex<double>(0.,1.) * input.SinPhiMu);
   const double BMu = Sqr(2.0 * Abs(susyMu));
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<std::complex<double>,3,3> Yu(
      Eigen::Matrix<std::complex<double>,3,3>::Zero());
   Eigen::Matrix<std::complex<double>,3,3> Yd(
      Eigen::Matrix<std::complex<double>,3,3>::Zero());
   Eigen::Matrix<std::complex<double>,3,3> Ye(
      Eigen::Matrix<std::complex<double>,3,3>::Zero());
   Eigen::Matrix<std::complex<double>,3,3> mm0(
      Eigen::Matrix<std::complex<double>,3,3>::Zero());

   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);

   mm0 = Sqr(m0) * Eigen::Matrix<std::complex<double>,3,3>::Identity();

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(1);
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
   m.set_eta(input.etaInput);
}

void setup_CMSSMCPVSemiAnalytic(
   flexiblesusy::CMSSMCPVSemiAnalytic_mass_eigenstates& m,
   flexiblesusy::CMSSMCPVSemiAnalytic_input_parameters& input)
{
   input.m12 = 500.;
   input.Azero = -100.;
   input.TanBeta = 10.0;
   input.AbsMu = 300.;
   input.CosPhiMu = 1.;
   input.etaInput = 0.1;

   setup_CMSSMCPVSemiAnalytic_const(m, input);
}

#endif
