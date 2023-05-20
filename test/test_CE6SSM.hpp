#ifndef TEST_CE6SSM_H
#define TEST_CE6SSM_H

#include "CE6SSM_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

struct Boundary_values {
   double m12{};
   double Azero{};
   double m0Sq{};
   double MuPr{};
   double BMuPr{};
   double Lambda{};
   double Lambda12{};
   double Kappa{};
};

void setup_high_scale_CE6SSM_const(
   flexiblesusy::CE6SSM_mass_eigenstates& model,
   const Boundary_values& values)
{
   model.get_input().BMuPrimeInput = values.BMuPr;

   model.set_Lambdax(values.Lambda);
   model.set_Lambda12(values.Lambda12 * UNITMATRIX(2));
   model.set_Kappa(values.Kappa * UNITMATRIX(3));

   model.set_MuPr(values.MuPr);
   model.set_BMuPr(values.BMuPr);

   model.set_TYu(model.get_Yu() * values.Azero);
   model.set_TYd(model.get_Yd() * values.Azero);
   model.set_TYe(model.get_Ye() * values.Azero);
   model.set_TLambdax(values.Lambda * values.Azero);
   model.set_TKappa(values.Kappa * values.Azero * UNITMATRIX(3));
   model.set_TLambda12(values.Lambda12 * values.Azero * UNITMATRIX(2));

   model.set_mq2(values.m0Sq * UNITMATRIX(3));
   model.set_mu2(values.m0Sq * UNITMATRIX(3));
   model.set_md2(values.m0Sq * UNITMATRIX(3));
   model.set_ml2(values.m0Sq * UNITMATRIX(3));
   model.set_me2(values.m0Sq * UNITMATRIX(3));
   model.set_mDx2(values.m0Sq * UNITMATRIX(3));
   model.set_mDxbar2(values.m0Sq * UNITMATRIX(3));
   model.set_mH1I2(values.m0Sq * UNITMATRIX(2));
   model.set_mH2I2(values.m0Sq * UNITMATRIX(2));
   model.set_msI2(values.m0Sq * UNITMATRIX(2));
   model.set_mHd2(values.m0Sq);
   model.set_mHu2(values.m0Sq);
   model.set_ms2(values.m0Sq);
   model.set_mHp2(values.m0Sq);
   model.set_mHpbar2(values.m0Sq);

   model.set_MassB(values.m12);
   model.set_MassBp(values.m12);
   model.set_MassWB(values.m12);
   model.set_MassG(values.m12);

   model.set_m0Sq(values.m0Sq);
   model.set_m12(values.m12);
   model.set_Azero(values.Azero);
   model.set_MuPrBV(values.MuPr);
}

void setup_high_scale_CE6SSM(
   flexiblesusy::CE6SSM_mass_eigenstates& model,
   Boundary_values& values)
{
   using namespace flexiblesusy;

   values.m12 = 1000;
   values.Azero = -1000.;
   values.m0Sq = 6.25e6;
   values.MuPr = 10000.;
   values.BMuPr = 10000.;
   values.Lambda = model.get_input().LambdaInput;
   values.Lambda12 = model.get_input().Lambda12Input;
   values.Kappa = model.get_input().KappaInput;

   setup_high_scale_CE6SSM_const(model, values);
}

void setup_CE6SSM_const(flexiblesusy::CE6SSM_mass_eigenstates& m,
                        const flexiblesusy::CE6SSM_input_parameters& input)
{
   using namespace flexiblesusy;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5.0 * ALPHAMZ / (3.0 * (1.0 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4.0 * Pi * alpha1);
   const double gN = 1.01 * g1;
   const double g2 = Sqrt(4.0 * Pi * alpha2);
   const double g3 = Sqrt(4.0 * Pi * ALPHASMZ);
   const double lambda = input.LambdaInput;
   const double tanBeta = input.TanBeta;
   const double sinBeta = Sin(ArcTan(tanBeta));
   const double cosBeta = Cos(ArcTan(tanBeta));
   const double M12 = 300.;
   const double m0 = input.vsInput;
   const double a0 = -1000.;
   const double MuPr = input.MuPrimeInput;
   const double BMuPr = input.BMuPrimeInput;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double vs = input.vsInput;
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Ye(Eigen::Matrix<double,3,3>::Zero());

   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);

   Eigen::Matrix<double,3,3> kappa(Eigen::Matrix<double,3,3>::Zero());
   kappa(0,0) = input.KappaInput;
   kappa(1,1) = input.KappaInput;
   kappa(2,2) = input.KappaInput;

   Eigen::Matrix<double,2,2> lambda12(Eigen::Matrix<double,2,2>::Zero());
   lambda12(0,0) = input.Lambda12Input;
   lambda12(1,1) = input.Lambda12Input;

   Eigen::Matrix<double,3,3> mm033(Sqr(m0) * Eigen::Matrix<double,3,3>::Identity());
   Eigen::Matrix<double,2,2> mm022(Sqr(m0) * Eigen::Matrix<double,2,2>::Identity());

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(1);

   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_gN(gN);

   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);

   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
   m.set_Lambda12(lambda12);

   m.set_MuPr(MuPr);
   m.set_BMuPr(BMuPr);

   m.set_MassB(M12);
   m.set_MassBp(M12);
   m.set_MassWB(M12);
   m.set_MassG(M12);

   m.set_mq2(mm033);
   m.set_ml2(mm033);
   m.set_md2(mm033);
   m.set_mu2(mm033);
   m.set_me2(mm033);
   m.set_mDx2(mm033);
   m.set_mDxbar2(mm033);

   m.set_mH1I2(mm022);
   m.set_mH2I2(mm022);
   m.set_msI2(mm022);

   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_ms2(Sqr(m0));
   m.set_mHp2(Sqr(m0));
   m.set_mHpbar2(Sqr(m0));

   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);

   m.set_TLambdax(a0 * lambda);
   m.set_TKappa(a0 * kappa);
   m.set_TLambda12(a0 * lambda12);

   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vs(vs);

}

void setup_CE6SSM(flexiblesusy::CE6SSM_mass_eigenstates& m,
                  flexiblesusy::CE6SSM_input_parameters& input)
{
   input.TanBeta = 10.;
   input.LambdaInput = 0.2;
   input.KappaInput = 0.15;
   input.MuPrimeInput = 10000.;
   input.BMuPrimeInput = 10000.;
   input.vsInput = 6000.;
   input.Lambda12Input = 0.2;
   input.m0SqGuess = flexiblesusy::Sqr(6000.);
   input.m12Guess = 600.;
   input.AzeroGuess = 600.;

   setup_CE6SSM_const(m, input);
}

#endif
