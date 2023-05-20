/** \file numerics_legacy.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief numerical routines - differential equation solver, differentiator
   and function minimiser for instance
*/

#ifndef NUMERICS_LEGACY_H
#define NUMERICS_LEGACY_H

#include "linalg.h"
#include "mycomplex.h"

namespace softsusy {

/// func is user-supplied, h is an estimate of what step-size to start with
/// and err returns error flags
double calcDerivative(double (*func)(double), 
		     double x, double h, double *err);

/// f is user-defined function, minimum value returned in xmin. Based on a
/// golden section search
double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin);

double dilog(double x);

/// These three functions are for the calculation of 2-loop log pieces of g-2
/// of the muon
double fps(double z);
double fs(double z);
double ffbar(double z);

Complex dilog(const Complex& x);

/// Returns a 3 by 3 real mixing matrix. Input angles are standard CKM
/// parameterisation. If the phase d is not zero, the result is only an
/// approximation to the full complex matrix: see SOFTSUSY manual for
/// details. 
DoubleMatrix display3x3RealMixing(double theta12, double theta13,
				  double theta23, double d);

/// useful for 2-loop mb/mt corrections
double fin(double mm1, double mm2);
double den(double a, int b); /// 1/a^b
/// Given r and qt, carry out a Jacobi rotation on rows i and i+1 of each
/// matrix. a and b are parameters of the rotation: $\cos
/// \theta=a/\sqrt{a^2+b^2}, \sin \theta = b / \sqrt{a^2 + b^2}$.
void rotate(DoubleMatrix& r, DoubleMatrix& qt, int n, int i, float a, float b);

} // namespace softsusy

#endif
