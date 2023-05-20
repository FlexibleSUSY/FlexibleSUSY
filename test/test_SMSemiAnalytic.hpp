#ifndef TEST_SMSemiAnalytic_H
#define TEST_SMSemiAnalytic_H

#include "SMSemiAnalytic_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

struct Boundary_values {
   double mu20{};
};

void setup_high_scale_SMSemiAnalytic_const(
   flexiblesusy::SMSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   model.set_mu2(values.mu20);
   model.set_mu20(values.mu20);
}

void setup_high_scale_SMSemiAnalytic(
   flexiblesusy::SMSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   values.mu20 = 8000.;

   setup_high_scale_SMSemiAnalytic_const(model, values);
}

void setup_SMSemiAnalytic_const(
   flexiblesusy::SMSemiAnalytic_mass_eigenstates& m,
   const flexiblesusy::SMSemiAnalytic_input_parameters& input)
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
   const double lambda = input.LambdaIN;
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

   m.set_mu2(Sqr(vev));
   m.set_v(vev);
}

void setup_SMSemiAnalytic(
   flexiblesusy::SMSemiAnalytic_mass_eigenstates& m,
   flexiblesusy::SMSemiAnalytic_input_parameters& input)
{
   input.LambdaIN = 0.192;
   input.Qin = 1000.;
   input.QEWSB = 173.3;

   setup_SMSemiAnalytic_const(m, input);
}

#endif
