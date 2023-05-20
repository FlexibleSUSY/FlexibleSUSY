
#ifndef TEST_lowNUHMSSMSEMIANALYTIC_H
#define TEST_lowNUHMSSMSEMIANALYTIC_H

#include "lowNUHMSSMSemiAnalytic_input_parameters.hpp"
#include "lowNUHMSSMSemiAnalytic_mass_eigenstates.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

struct Boundary_values {
   double m0{};
   double m12{};
   double Azero{};
   double mHd2{};
   double mHu2{};
   double Mu{};
   double BMu{};
};

void setup_susy_scale_lowNUHMSSMSemiAnalytic_const(
   flexiblesusy::lowNUHMSSMSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   using namespace flexiblesusy;

   model.get_input().m0 = values.m0;
   model.get_input().m12 = values.m12;
   model.get_input().Azero = values.Azero;
   model.get_input().MuInput = values.Mu;
   model.get_input().BMuInput = values.BMu;

   model.set_TYu(model.get_Yu() * values.Azero);
   model.set_TYd(model.get_Yd() * values.Azero);
   model.set_TYe(model.get_Ye() * values.Azero);

   model.set_mq2(Sqr(values.m0) * UNITMATRIX(3));
   model.set_mu2(Sqr(values.m0) * UNITMATRIX(3));
   model.set_md2(Sqr(values.m0) * UNITMATRIX(3));
   model.set_ml2(Sqr(values.m0) * UNITMATRIX(3));
   model.set_me2(Sqr(values.m0) * UNITMATRIX(3));
   model.set_mHd2(values.mHd2);
   model.set_mHu2(values.mHu2);

   model.set_MassB(values.m12);
   model.set_MassWB(values.m12);
   model.set_MassG(values.m12);

   model.set_Mu(values.Mu);
   model.set_BMu(values.BMu);

   model.set_mHd20(values.mHd2);
   model.set_mHu20(values.mHu2);
   model.set_MuBV(values.Mu);
}

void setup_susy_scale_lowNUHMSSMSemiAnalytic(
   flexiblesusy::lowNUHMSSMSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   using namespace flexiblesusy;

   values.m12 = model.get_input().m12;
   values.Azero = model.get_input().Azero;
   values.m0 = model.get_input().m0;
   values.BMu = model.get_BMu();
   values.Mu = model.get_Mu();
   values.mHd2 = model.get_mHd2();
   values.mHu2 = model.get_mHu2();

   setup_susy_scale_lowNUHMSSMSemiAnalytic_const(model, values);
}

void setup_lowNUHMSSMSemiAnalytic_const(
   flexiblesusy::lowNUHMSSMSemiAnalytic_mass_eigenstates& m,
   const flexiblesusy::lowNUHMSSMSemiAnalytic_input_parameters& input)
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
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double mHd2 = Sqr(input.m0);
   const double mHu2 = -Sqr(0.5 * input.m0);
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.MuInput;
   const double BMu = input.BMuInput;
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Ye(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> mm0(Eigen::Matrix<double,3,3>::Zero());

   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);

   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

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
   m.set_mHd2(mHd2);
   m.set_mHu2(mHu2);
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);
}

void setup_lowNUHMSSMSemiAnalytic(flexiblesusy::lowNUHMSSMSemiAnalytic_mass_eigenstates& m,
                             flexiblesusy::lowNUHMSSMSemiAnalytic_input_parameters& input)
{
   input.m0 = 800.;
   input.m12 = 500.;
   input.Azero = 0.;
   input.TanBeta = 10.0;
   input.MuInput = 130.;
   input.BMuInput = 3000.;

   setup_lowNUHMSSMSemiAnalytic_const(m, input);
}

#endif
