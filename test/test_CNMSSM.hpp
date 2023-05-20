#ifndef TEST_CNMSSM_H
#define TEST_CNMSSM_H

#include "CNMSSM_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

struct Boundary_values {
   double m12{};
   double Azero{};
   double m0Sq{};
};

void setup_high_scale_CNMSSM_const(
   flexiblesusy::CNMSSM_mass_eigenstates& model,
   const Boundary_values& values)
{
   model.get_input().Azero = values.Azero;
   model.get_input().m12 = values.m12;
   model.get_input().LambdaInput = model.get_Lambdax();

   model.set_TYu(model.get_Yu() * values.Azero);
   model.set_TYd(model.get_Yd() * values.Azero);
   model.set_TYe(model.get_Ye() * values.Azero);
   model.set_TLambdax(model.get_Lambdax() * values.Azero);
   model.set_TKappa(model.get_Kappa() * values.Azero);

   model.set_mq2(values.m0Sq * UNITMATRIX(3));
   model.set_mu2(values.m0Sq * UNITMATRIX(3));
   model.set_md2(values.m0Sq * UNITMATRIX(3));
   model.set_ml2(values.m0Sq * UNITMATRIX(3));
   model.set_me2(values.m0Sq * UNITMATRIX(3));
   model.set_mHd2(values.m0Sq);
   model.set_mHu2(values.m0Sq);
   model.set_ms2(values.m0Sq);

   model.set_MassB(values.m12);
   model.set_MassWB(values.m12);
   model.set_MassG(values.m12);

   model.set_m0Sq(values.m0Sq);
}

void setup_high_scale_CNMSSM(
   flexiblesusy::CNMSSM_mass_eigenstates& model,
   Boundary_values& values)
{
   using namespace flexiblesusy;

   values.m12 = 133.33;
   values.Azero = -300.;
   values.m0Sq = Sqr(40.);

   setup_high_scale_CNMSSM_const(model, values);
}

void setup_CNMSSM_const(flexiblesusy::CNMSSM_mass_eigenstates& m,
                        const flexiblesusy::CNMSSM_input_parameters& input)
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
   const double lambda = input.LambdaInput;
   const double kappa = 0.0297;
   const double tanBeta = input.TanBeta;
   const double sinBeta = Sin(ArcTan(tanBeta));
   const double cosBeta = Cos(ArcTan(tanBeta));
   const double M12 = input.m12;
   const double m0 = 0.333 * input.m12;
   const double a0 = input.Azero;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double vS = 7000.;
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Ye(Eigen::Matrix<double,3,3>::Zero());

   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);

   Eigen::Matrix<double,3,3> mm0(Sqr(m0) * Eigen::Matrix<double,3,3>::Identity());

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
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
   m.set_ms2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_TLambdax(a0 * lambda);
   m.set_TKappa(a0 * kappa);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vS(vS);
}

void setup_CNMSSM(flexiblesusy::CNMSSM_mass_eigenstates& m,
                  flexiblesusy::CNMSSM_input_parameters& input)
{
   input.m12 = 133.33;
   input.TanBeta = 10.;
   input.Azero = -300.;
   input.LambdaInput = -0.05;
   input.SignvS = 1;

   setup_CNMSSM_const(m, input);
}

#endif
