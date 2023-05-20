#ifndef TEST_THDMIIEWSBAtMZSemiAnalytic_H
#define TEST_THDMIIEWSBAtMZSemiAnalytic_H

#include "THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

struct Boundary_values {
   double M112{};
   double M122{};
   double M222{};
};

void setup_high_scale_THDMIIEWSBAtMZSemiAnalytic_const(
   flexiblesusy::THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   using namespace flexiblesusy;

   THDMIIEWSBAtMZSemiAnalytic_input_parameters& input(model.get_input());

   input.M122IN = values.M122;

   model.set_M112(values.M112);
   model.set_M122(values.M122);
   model.set_M222(values.M222);

   model.set_M112IN(values.M112);
   model.set_M222IN(values.M222);
}

void setup_high_scale_THDMIIEWSBAtMZSemiAnalytic(
   flexiblesusy::THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates& model,
   Boundary_values& values)
{
   values.M112 = 20000.;
   values.M122 = 1000.;
   values.M222 = -7000.;

   setup_high_scale_THDMIIEWSBAtMZSemiAnalytic_const(model, values);
}

void setup_THDMIIEWSBAtMZSemiAnalytic_const(
   flexiblesusy::THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates& m,
   const flexiblesusy::THDMIIEWSBAtMZSemiAnalytic_input_parameters& input)
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
   const double root2 = Sqrt(2.0);
   const double tb = input.TanBeta;
   const double cb = Cos(ArcTan(tb));
   const double sb = Sin(ArcTan(tb));
   const double vev = 246.0;
   const double v1 = vev * cb / root2;
   const double v2 = vev * sb / root2;
   const double scale = Electroweak_constants::MZ;
   const double lam1 = input.Lambda1IN;
   const double lam2 = input.Lambda2IN;
   const double lam3 = input.Lambda3IN;
   const double lam4 = input.Lambda4IN;
   const double lam5 = input.Lambda5IN;
   const double lam6 = input.Lambda6IN;
   const double lam7 = input.Lambda7IN;
   const double m122 = input.M122IN;

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

   m.set_Lambda1(lam1);
   m.set_Lambda2(lam2);
   m.set_Lambda3(lam3);
   m.set_Lambda4(lam4);
   m.set_Lambda5(lam5);
   m.set_Lambda6(lam6);
   m.set_Lambda7(lam7);

   m.set_M112(Sqr(vev));
   m.set_M122(m122);
   m.set_M222(-Sqr(vev));

   m.set_v1(v1);
   m.set_v2(v2);
}

void setup_THDMIIEWSBAtMZSemiAnalytic(
   flexiblesusy::THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates& m,
   flexiblesusy::THDMIIEWSBAtMZSemiAnalytic_input_parameters& input)
{
   input.Lambda1IN = 2.;
   input.Lambda2IN = 0.1;
   input.Lambda3IN = 0.5;
   input.Lambda4IN = 0.8;
   input.Lambda5IN = -2.;
   input.Lambda6IN = 0.;
   input.Lambda7IN = 0.;
   input.M122IN = 1000.;
   input.TanBeta = 10.;
   input.Qin = 126.;

   setup_THDMIIEWSBAtMZSemiAnalytic_const(m, input);
}

#endif
