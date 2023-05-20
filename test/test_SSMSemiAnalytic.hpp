#ifndef TEST_SSMSemiAnalytic_H
#define TEST_SSMSemiAnalytic_H

#include "SSMSemiAnalytic_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

struct Boundary_values {
   double K1{};
   double K2{};
   double Kappa{};
   double MS0{};
   double mu20{};
};

void setup_high_scale_SSMSemiAnalytic_const(
   flexiblesusy::SSMSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   using namespace flexiblesusy;

   SSMSemiAnalytic_input_parameters& input(model.get_input());

   input.K1input = values.K1;
   input.K2input = values.K2;
   input.Kappainput = values.Kappa;

   model.set_K1(values.K1);
   model.set_K2(values.K2);
   model.set_Kappa(values.Kappa);
   model.set_MS(values.MS0);
   model.set_mu2(values.mu20);

   model.set_MS0(values.MS0);
   model.set_mu20(values.mu20);
}

void setup_high_scale_SSMSemiAnalytic(
   flexiblesusy::SSMSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   values.K1 = -100.;
   values.K2 = 0.1;
   values.Kappa = 100.;
   values.MS0 = 100.;
   values.mu20 = 8000.;

   setup_high_scale_SSMSemiAnalytic_const(model, values);
}

void setup_SSMSemiAnalytic_const(
   flexiblesusy::SSMSemiAnalytic_mass_eigenstates& m,
   const flexiblesusy::SSMSemiAnalytic_input_parameters& input)
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
   const double lambda = input.Lambdainput;
   const double lambdas = input.LambdaSinput;
   const double kappa = input.Kappainput;
   const double k1 = input.K1input;
   const double k2 = input.K2input;
   const double vs = input.vSInput;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Ye(Eigen::Matrix<double,3,3>::Zero());

   Yu(2,2) = -165.0   * root2 / vev;
   Yd(2,2) = 2.9     * root2 / vev;
   Ye(2,2) = 1.77699 * root2 / vev;

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
   m.set_LambdaS(lambdas);
   m.set_K2(k2);

   m.set_Kappa(kappa);
   m.set_K1(k1);
   m.set_MS(Sqr(vev));
   m.set_mu2(Sqr(vev));
   m.set_v(vev);
   m.set_vS(vs);
}

void setup_SSMSemiAnalytic(
   flexiblesusy::SSMSemiAnalytic_mass_eigenstates& m,
   flexiblesusy::SSMSemiAnalytic_input_parameters& input)
{
   input.Qin = 1000.;
   input.QEWSB = 173.3;
   input.Lambdainput = 0.19;
   input.LambdaSinput = 0.1;
   input.Kappainput = 100.;
   input.K1input = -100.;
   input.K2input = 0.1;
   input.vSInput = 3.;

   setup_SSMSemiAnalytic_const(m, input);
}

#endif
