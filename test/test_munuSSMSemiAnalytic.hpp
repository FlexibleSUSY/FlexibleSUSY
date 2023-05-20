
#ifndef TEST_munuSSMSemiAnalytic_H
#define TEST_munuSSMSemiAnalytic_H

#include "munuSSMSemiAnalytic_input_parameters.hpp"
#include "munuSSMSemiAnalytic_mass_eigenstates.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

struct Boundary_values {
   double m12{};
   double Azero{};
   double m0{};
   double mHd20{};
   double mHu20{};
   double mv20{};
   Eigen::Matrix<double,3,1> mlHd20{Eigen::Matrix<double,3,1>::Zero()};
};

void setup_high_scale_munuSSMSemiAnalytic_const(
   flexiblesusy::munuSSMSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   using namespace flexiblesusy;

   const auto m12 = values.m12;
   const auto Azero = values.Azero;
   const auto m0 = values.m0;
   const auto mHd20 = values.mHd20;
   const auto mHu20 = values.mHu20;
   const auto mv20 = values.mv20;
   const auto mlHd20 = values.mlHd20;
   const auto LambdaInput = model.get_input().LambdaInput;
   const auto Ye = model.get_Ye();
   const auto Yd = model.get_Yd();
   const auto Yu = model.get_Yu();
   const auto Yv = model.get_Yv();
   const auto Kappa = model.get_Kappa();

   model.get_input().m12 = m12;
   model.get_input().Azero = Azero;
   model.get_input().m0 = m0;

   model.set_TYe((Azero*Ye).real());
   model.set_TYd((Azero*Yd).real());
   model.set_TYu((Azero*Yu).real());
   model.set_TYv((Azero*Yv).real());
   model.set_mq2((Sqr(m0)*UNITMATRIX(3)).real());
   model.set_ml2((Sqr(m0)*UNITMATRIX(3)).real());
   model.set_md2((Sqr(m0)*UNITMATRIX(3)).real());
   model.set_mu2((Sqr(m0)*UNITMATRIX(3)).real());
   model.set_me2((Sqr(m0)*UNITMATRIX(3)).real());
   model.set_mHd2(Re(mHd20));
   model.set_mHu2(Re(mHu20));
   model.set_mv2(Re(mv20));
   model.set_mlHd2((mlHd20).real());
   model.set_TKappa(Re(Azero*Kappa));
   model.set_TLambdax(Re(Azero*LambdaInput));
   model.set_MassB(Re(m12));
   model.set_MassWB(Re(m12));
   model.set_MassG(Re(m12));

   model.set_mHd20(mHd20);
   model.set_mHu20(mHu20);
   model.set_mv20(mv20);
   model.set_mlHd20(mlHd20);
}

void setup_high_scale_munuSSMSemiAnalytic(
   flexiblesusy::munuSSMSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   values.m12 = 300.;
   values.Azero = -10.;
   values.m0 = 500.;
   values.mHd20 = 1.17385044e5;
   values.mHu20 = 5.28526767e5;
   values.mv20 = 1.66738568e4;
   values.mlHd20(0) = -1.18767390e5;
   values.mlHd20(1) = -1.18764617e5;
   values.mlHd20(2) = -1.17971011e5;

   setup_high_scale_munuSSMSemiAnalytic_const(model, values);
}

void setup_munuSSMSemiAnalytic_const(
   flexiblesusy::munuSSMSemiAnalytic_mass_eigenstates& m,
   const flexiblesusy::munuSSMSemiAnalytic_input_parameters& input)
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
   const double lambda = input.LambdaInput;
   const double kappa = input.KappaInput;
   const auto Yv = input.YvInput;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double vr = input.vRInput;
   const double vl1 = input.vL1Input;
   const double vl2 = input.vL2Input;
   const double vl3 = input.vL3Input;
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
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
   m.set_Yv(Yv);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vR(vr);
   m.set_vL(0, vl1);
   m.set_vL(1, vl2);
   m.set_vL(2, vl3);

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
   m.set_mv2(Sqr(m0));
   m.set_mlHd2(0, Sqr(m0));
   m.set_mlHd2(1, Sqr(m0));
   m.set_mlHd2(2, Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_TYv(a0 * Yv);
   m.set_TLambdax(a0 * lambda);
   m.set_TKappa(a0 * kappa);
}

void setup_munuSSMSemiAnalytic(flexiblesusy::munuSSMSemiAnalytic_mass_eigenstates& m,
                               flexiblesusy::munuSSMSemiAnalytic_input_parameters& input)
{
   input.TanBeta = 10.;
   input.vRInput = 10.;
   input.vL1Input = 10.;
   input.vL2Input = 10.;
   input.vL3Input = 10.;
   input.m0 = 500.;
   input.m12 = 300.;
   input.Azero = -10.;
   input.LambdaInput = 0.2;
   input.KappaInput = 0.1;
   input.YvInput(0) = 0.01;
   input.YvInput(1) = 0.01;
   input.YvInput(2) = 0.01;

   setup_munuSSMSemiAnalytic_const(m, input);
}

#endif
