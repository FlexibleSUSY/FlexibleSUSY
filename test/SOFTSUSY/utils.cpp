
/** \file utils.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

 */

#include "utils.h"

double frexp(const softsusy::Complex & c, int * i) {
  int a, b;
  frexp(c.real(), &a);
  frexp(c.imag(), &b);
  if (std::abs(a) > std::abs(b)) return *i = a;
  *i = b;
  return 0.;
}

void FPE_ExceptionHandler(int /* nErrType */) {
  throw "SIGFPE"; ///< This reverts back to softsusy code!
}

int theta(double a) {
  int temp = 0;
  if (a > 0.0) temp = 1;
  return temp;
}

// Just sets precision and format of outputs
void outputCharacteristics(int n) {
  std::cin.setf(std::ios::scientific, std::ios::floatfield);
  std::cin.precision(n);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.precision(n);
  std::cerr.setf(std::ios::scientific, std::ios::floatfield);
  std::cerr.precision(n);
}

// Finds fractional difference between |a| and |b|
double toleranceCheck(double a, double b) {
  double sTin = fabs(a), sTout = fabs(b);
  double maxx = maximum(sTin, sTout);

  const double underflow = 1.0e-20;

  if (maxx < underflow) return 0.0;
  return fabs(1.0 - minimum(sTin, sTout) / maxx);
}

// Outputs a space if greater than zero, a minus otherwise.
// Useful for outputting negative numbers in rows
void printRow(std::ostream & out, double x) {

  // make it return a character when you've worked out the equivalent of printf

  double underflow = 1.0e-120;
  if (fabs(x) < underflow) x = 0.0; // Traps -0.0
  if (x >= 0.0) out << " " << x;
  else out << x;
}

bool testNan(double f) {
  return (f != f);
}

bool close(double m1, double m2, double tol) {
  const double mmax = fabs(maximum(fabs(m1), fabs(m2)));
  const double mmin = fabs(minimum(fabs(m1), fabs(m2)));
  const double max_tol = tol * mmax;
  if (max_tol == 0.0 && mmax != 0.0 && tol != 0.0)
    return (mmax - mmin <= tol);

  return (mmax - mmin <= max_tol);
}
