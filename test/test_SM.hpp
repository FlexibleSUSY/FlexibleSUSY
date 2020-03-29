
#ifndef TEST_SM_H
#define TEST_SM_H

#include "wrappers.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;

template <class TSM, class TInput>
void setup_SM_const(TSM& m, const TInput& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double lambda = input.LambdaIN;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero());
   Yu(0,0) = 0.001;
   Yu(1,1) = 0.001;
   Yu(2,2) = 165.0   * root2 / vev;
   Yd(0,0) = 0.001;
   Yd(1,1) = 0.001;
   Yd(2,2) = 2.9     * root2 / vev;
   Ye(0,0) = 5.11e-4 * root2 / vev;
   Ye(1,1) = 0.105  * root2 / vev;
   Ye(2,2) = 1.77699 * root2 / vev;

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_Lambdax(lambda);
   m.set_v(vev);

   m.solve_ewsb_tree_level();
}

#endif
