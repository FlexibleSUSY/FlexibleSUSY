
/** \file utils.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief A few handy bits and pieces - little mathematical functions 
   and the like
*/

#ifndef UTILS_H
#define UTILS_H

#include "mycomplex.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include "def.h"
#include <stdio.h>
#include <string>
#include <exception>
/// Error handling for softsusy
#include <csignal>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <signal.h>
#include <dirent.h>

/// returns either sqrt(f) for f>0 or 0 otherwise
inline double zeroSqrt(double f){ if (f > 0.) return sqrt(f); 
  else return softsusy::EPSTOL;
}
/// Gives exponent of the largest number: either imaginary or real part
double frexp(const softsusy::Complex & c, int * i);
/// Exception handler - will even trap errors in the fortran portion
void FPE_ExceptionHandler(int nErrType);
/// Standard theta function: 1 is a>0, 0 otherwise
int theta(double a);
/// Sets number of decimal places to display in scientific output
void outputCharacteristics(int);
/// square of a number
inline double sqr(double a) { return a * a; }
/// maximum of two numbers
inline double maximum(double a, double b) { return ((a > b) ? a : b); }
/// minimum of a and b
inline double minimum(double a, double b) { return ((a < b) ? a : b); }
/// minimum of a and b
inline int minimum(int a, int b) { return ((a < b) ? a : b); }
/// Finds fractional difference between |a| and |b|
double toleranceCheck(double sTin, double sTout);

inline double minimum(double a, double b, double c) {
   if(a <= b && a <= c) return a;
   else if(b <= a && b <= c) return b;
   else return c;
}

inline double maximum(double a, double b, double c) {
   if(a >= b && a >= c) return a;
   else if(b >= a && b >= c) return b;
   else return c;
}

inline double middle(double a, double b, double c) {
   if((a >= b && a <= c) ||(a >= c && a <= b)   ) return a;
   else if((b >= a && b <= c) ||(b >= c && b <= a) ) return b;
   else  return c;
}

inline int massorder(double & a, double & b, double & c) {
   double anew, bnew, cnew;
   anew = minimum(a,b,c);
   cnew= maximum(a,b,c);
   bnew = middle(a,b,c);

   a = anew;
   b = bnew;
   c = cnew;
   return 0;
}

/// checks if ABSOLUTE (or squared) values are closer than tol, or both
/// numbers are smaller than EPSTOL
bool close(double m1, double m2, double tol);

/// Returns |a| with sign of b in front
inline double sign(double a, double b) 
{ return ((b) >= 0.0 ? fabs(a) : -fabs(a)); }

/// gives the sign (+-1) of x
inline int sgn(double x)
{ return (x >= 0.0 ? 1 : -1); }

/// Outputs a space if greater than zero, a minus otherwise.
/// Useful for outputting negative numbers in rows
void printRow(std::ostream & o, double x);

/// Returns true if f's a nan. Unfortunately, the behaviour of this function
/// is not standardised across all platforms yet. 
bool testNan(double f);

#endif

