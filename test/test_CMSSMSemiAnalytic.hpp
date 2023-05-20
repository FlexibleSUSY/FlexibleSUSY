
#ifndef TEST_CMSSMSEMIANALYTIC_H
#define TEST_CMSSMSEMIANALYTIC_H

#include "CMSSMSemiAnalytic_input_parameters.hpp"
#include "CMSSMSemiAnalytic_mass_eigenstates.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

struct Boundary_values {
   double m12{};
   double Azero{};
   double m0Sq{};
   double BMu0{};
   double Mu{};
};

void setup_high_scale_CMSSMSemiAnalytic_const(
   flexiblesusy::CMSSMSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   model.get_input().m12 = values.m12;
   model.get_input().Azero = values.Azero;

   model.set_TYu(model.get_Yu() * values.Azero);
   model.set_TYd(model.get_Yd() * values.Azero);
   model.set_TYe(model.get_Ye() * values.Azero);

   model.set_mq2(values.m0Sq * UNITMATRIX(3));
   model.set_mu2(values.m0Sq * UNITMATRIX(3));
   model.set_md2(values.m0Sq * UNITMATRIX(3));
   model.set_ml2(values.m0Sq * UNITMATRIX(3));
   model.set_me2(values.m0Sq * UNITMATRIX(3));
   model.set_mHd2(values.m0Sq);
   model.set_mHu2(values.m0Sq);

   model.set_MassB(values.m12);
   model.set_MassWB(values.m12);
   model.set_MassG(values.m12);

   model.set_Mu(values.Mu);
   model.set_BMu(values.BMu0);

   model.set_m0Sq(values.m0Sq);
   model.set_BMu0(values.BMu0);
   model.set_MuBV(values.Mu);
}

void setup_high_scale_CMSSMSemiAnalytic(
   flexiblesusy::CMSSMSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   using namespace flexiblesusy;

   values.m12 = 300.0;
   values.Azero = -200.;
   values.m0Sq = Sqr(250);
   values.BMu0 = -1.e4;
   values.Mu = model.get_Mu();

   setup_high_scale_CMSSMSemiAnalytic_const(model, values);
}

void setup_CMSSMSemiAnalytic_const(flexiblesusy::CMSSMSemiAnalytic_mass_eigenstates& m,
                                   const flexiblesusy::CMSSMSemiAnalytic_input_parameters& input)
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
   const double susyMu = input.MuInput;
   const double BMu = Sqr(2.0 * susyMu);
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
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);
}

void setup_CMSSMSemiAnalytic(flexiblesusy::CMSSMSemiAnalytic_mass_eigenstates& m,
                             flexiblesusy::CMSSMSemiAnalytic_input_parameters& input)
{
   input.m12 = 300.;
   input.Azero = -200.;
   input.TanBeta = 7.0;
   input.MuInput = 500.;

   setup_CMSSMSemiAnalytic_const(m, input);
}

#endif
