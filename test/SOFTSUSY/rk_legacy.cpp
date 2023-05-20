
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "utils.h"
#include "def.h"
#include "rk_legacy.hpp"

using namespace softsusy;

// returns >0 if there's a problem:
int integrateOdes(DoubleVector & ystart, double from, double to, double eps,
	      double h1, double hmin, 
	      DoubleVector (*derivs)(double, const DoubleVector &),
	      int (*rkqs)(DoubleVector & y, const DoubleVector & dydx, double
			   *x, double htry, double eps, DoubleVector & yscal,
			   double *hdid, double *hnext, 
			   DoubleVector (*derivs)(double, const DoubleVector
						  &)) ) {  
  int nvar =  ystart.displayEnd();
  int nstp, i;
  double x, hnext, hdid, h;
  DoubleVector yscal(nvar), y(ystart), dydx(nvar);
  
  x = from;
  h = sign(h1, to - from);
  
  const int MAXSTP = 400;
  const double TINY = 1.0e-16;

  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    dydx = (*derivs)(x, y);
    for (i = 1; i <= nvar; i++)
      yscal(i) = fabs(y(i)) + fabs(dydx(i) * h) + TINY;
    if ((x + h - to) * (x + h - from) > 0.0) h = to - x;
    int smallStep = (*rkqs)(y, dydx, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (smallStep) return 1;

    if ((x - to) * (to - from) >= 0.0) {
      for (i = 1; i<= nvar; i++) ystart(i) = y(i);
      return 0;
    }
      
    if (fabs(hnext) <= hmin) {
      nstp = MAXSTP; // bail out
      if (PRINTOUT > 1) {
	std::cout << "Step size too small in rk_legacy.cpp:integrateOdes\n";
	std::cout << "**********x = " << x << "*********\n";
	for (i = 1;i<= nvar;i++) 
	  std::cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	    ") = " << dydx(i) << std::endl;
	std::cout.flush();
      }
    }
    
    h = hnext;
  }
  
  if (PRINTOUT > 1) {
    std::cout << "Bailed out of rk_legacy.cpp:too many steps in integrateOdes\n";
    std::cout << "**********x = " << x << "*********\n";
    for (i = 1;i<= nvar;i++) 
      std::cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	") = " << dydx(i) << std::endl;
    std::cout.flush();
  }
  
  return 1;
}

int odeStepper(DoubleVector & y, const DoubleVector & dydx, double *x, double
		htry, double eps, DoubleVector & yscal, double *hdid, 
		double *hnext,		
		DoubleVector (*derivs)(double, const DoubleVector &))
{
  const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;

  int i, n = y.displayEnd();
  double errmax, h, htemp, xnew;
  
  DoubleVector yerr(n), ytemp(n);
  h = htry;
  for (;;) {
    rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, derivs);
    errmax = 0.0;
    for (i = 1; i<= n;i++) errmax = maximum(errmax, fabs(yerr(i) / yscal(i)));
    errmax  /= eps;
    if (errmax <= 1.0) break;
    htemp = SAFETY * h * pow(errmax, PSHRNK);
    h = (h >= 0.0 ? maximum(htemp ,0.1 * h) : minimum(htemp, 0.1 * h));
    xnew = (*x) + h;
    if (xnew == *x) 
      {
	if (PRINTOUT) {
	std::cout << "At x = " << *x;
	std::cout << ",stepsize underflow in odeStepper" << std::flush << std::endl;
	}
	return 1;
      }
  }
  if (errmax > ERRCON) *hnext = SAFETY * h * pow(errmax,PGROW);
  else *hnext = 5.0 * h;
  *x += (*hdid = h);
  y = ytemp;
  return 0;
}

void rungeKuttaStep(const DoubleVector & y, const DoubleVector & dydx, 
	     double x, double h, DoubleVector & yout, DoubleVector & yerr, 
	     DoubleVector (*derivs)(double, const DoubleVector &)) {
  int i;
  const static double a2 = 0.2,a3 = 0.3,a4 = 0.6,a5 = 1.0,a6 = 0.875,b21 =
    0.2,b31 = 3.0 / 40.0,b32 = 9.0 / 40.0,b41 = 0.3,b42 = -0.9,b43 = 1.2,
    b51 = -11.0 / 54.0, b52 = 2.5,b53 = -70.0 / 27.0,b54 = 35.0 / 27.0,
    b61 = 1631.0 / 55296.0,b62 = 175.0 / 512.0,b63 = 575.0 / 13824.0,
    b64 = 44275.0 / 110592.0,b65 = 253.0 / 4096.0,c1 = 37.0 / 378.0,
    c3 = 250.0 / 621.0,c4 = 125.0 / 594.0,c6 = 512.0 / 1771.0,
    dc5 = -277.00 / 14336.0;
  const double dc1 = c1-2825.0 / 27648.0,dc3 = c3-18575.0 / 48384.0,
    dc4 = c4-13525.0 / 55296.0,dc6 = c6-0.25;
  
  int n = y.displayEnd();

  DoubleVector ytemp(b21 * h * dydx.display() + y.display());
  DoubleVector ak2((*derivs)(x + a2 * h, ytemp));

  // Allowing piece-wise calculating of ytemp for speed reasons
  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b31 * dydx.display(i) + b32 * ak2(i));
  DoubleVector ak3((*derivs)(x + a3 * h, ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b41 * dydx.display(i) + b42 * ak2(i) + b43
				   * ak3(i));
  DoubleVector ak4((*derivs)(x+a4*h,ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b51 * dydx.display(i) + b52 * ak2(i) + b53
				   * ak3(i) + b54 * ak4(i));
  DoubleVector ak5((*derivs)(x + a5 * h, ytemp));

  for (i = 1; i<= n; i++)
    ytemp(i) = y.display(i) + h * (b61 * dydx.display(i) + b62 * ak2(i) + b63
				   * ak3(i) + b64 * ak4(i) + b65 * ak5(i));
  DoubleVector ak6((*derivs)(x + a6 * h, ytemp));

  for (i = 1; i<= n; i++)
    yout(i) = y.display(i) + h * (c1 * dydx.display(i) + c3 * ak3(i) + c4 *
				  ak4(i) + c6 * ak6(i));
  for (i = 1; i<= n; i++)
    yerr(i) = h * (dc1 * dydx.display(i) + dc3 * ak3(i) + 
		   dc4 * ak4(i) + dc5 * ak5(i) + dc6 * ak6(i));
}
