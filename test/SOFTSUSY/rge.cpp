
/** \file rge.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: general abstract renormalisation group equation object

*/

#include "rge.h"

#include <cmath>

namespace softsusy {

static RGE * tempRge;

// runto/run functions return >0 if there's a problem with the running
int RGE::runto(double x2, double eps) {
  double tol;
  if (eps < 0.0) tol = TOLERANCE;
  else if (eps < EPSTOL) tol = EPSTOL;
  else tol = eps;
  double x1 = this -> displayMu();
  return run(x1, x2, tol);
}

int RGE::run(double x1, double x2, double eps) {
  double tol;

  if (eps < 0.0) tol = TOLERANCE;
  else if (eps < EPSTOL) tol = EPSTOL;
  else tol = eps;

  tempRge = this;
  DoubleVector y(this -> display());
  int err = callRK(x1, x2, y, allDerivs, tol);
  if (err == 0) {
    this -> set(y);
    tempRge -> setMu(x2);
  }
  return err;
}

// Runge-Kutte user defined routine: given log renorm scale x and the dependent
// variables y of an RGE, will calculate the derivitives dydx.
DoubleVector allDerivs(double x, const DoubleVector & y)
{
  tempRge->setMu(std::exp(x));
  tempRge->set(y);
  return tempRge->beta();
}

//Does the actual calling of Runge Kutta: default precision is TOLERANCE defined in
//def.h
//Returns >0 if there's a problem with the running
int RGE::callRK(double x1, double x2, DoubleVector & v,
		DoubleVector (*derivs)(double, const DoubleVector &), 
		double eps) {
  using std::fabs;
  using std::log;

  double tol;
  if (eps < 0.0) tol = TOLERANCE;
  else if (eps < EPSTOL) tol = EPSTOL;
  else tol = eps;
  // x1 == x2 with high precision
  if (close(fabs(x1), fabs(x2), EPSTOL)) return 0;

  // RGE in terms of natural log of renormalisation scale
  double from = log(fabs(x1));
  double to = log(fabs(x2));

  double guess = (from - to) * 0.1; //first step size
  double hmin = (from - to) * tol * 1.0e-5; 

  int err =
    integrateOdes(v, from, to, tol, guess, hmin, derivs, odeStepper);
  
  setMu(x2);
  return err;
}

} // namespace softsusy
