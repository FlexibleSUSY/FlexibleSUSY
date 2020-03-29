/** \file numerics_legacy.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "numerics_legacy.h"
#include "dilog.hpp"
#include "error.hpp"
#include "utils.h"
#include <cmath>

namespace flexiblesusy {

class SoftsusyNumericsError : public Error {
public:
   explicit SoftsusyNumericsError(std::string msg)
      : Error(msg) {}
   virtual ~SoftsusyNumericsError() {}
};

} // namespace flexiblesusy

namespace softsusy {

namespace {
inline void shft2(double & a, double & b, double c) { a = b; b = c; }
inline void shft3(double & a, double & b, double & c, double d) {
   a = b; b = c; c = d;
}
} // anonymous namespace

double calcDerivative(double (*func)(double), double x, double h, double
		      *err){
  const double CON = 1.4, CON2 = CON * CON, BIG = 1.0e30, 
    SAFE = 2.0; 
  const int NTAB = 10;
  
  int i, j;
  double errt, fac, hh, ans = 0.0;
  
  if (h == 0.0) throw flexiblesusy::SoftsusyNumericsError("h must be nonzero in numerics.cpp:calcDerivative");


  DoubleMatrix a(NTAB, NTAB);
  hh = h;
  a(1, 1) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
  *err = BIG;
  for (i=2; i<=NTAB; i++) {
    hh /= CON;
    a(1, i) = ((*func)(x + hh) - (*func)(x - hh)) / (2.0 * hh);
    fac = CON2;
    for (j=2; j<=i; j++) {
      a(j, i) = (a(j-1, i) * fac - a(j-1, i-1)) / (fac - 1.0);
      fac = CON2 * fac;
      errt = maximum(fabs(a(j, i) - a(j-1, i)), fabs(a(j, i) - a(j-1, i-1)));
      if (errt <= *err) {
	*err = errt;
	ans = a(j, i);
      }
    }
    if (fabs(a(i, i) - a(i-1, i-1)) >= SAFE * (*err)) break;
  }

  return ans;
}

double findMinimum(double ax, double bx, double cx, double (*f)(double),
		   double tol, double *xmin)
{
  const double R = 0.61803399, C = 1.0 - R;
  double f1, f2, x0, x1, x2, x3;
  
  x0 = ax; 
  x3 = cx; 
  if (fabs(cx - bx) > fabs(bx - ax)) {
    x1 = bx; 
    x2 = bx + C * (cx - bx); 
  } else {
    x2 = bx; 
    x1 = bx - C * (bx - ax); 
  }
  f1 = (*f)(x1); 
  f2 = (*f)(x2); 
  while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {
    if (f2 < f1) {
      shft3(x0, x1, x2, R * x1 + C * x3);
      shft2(f1, f2, (*f)(x2));
    } else {
      shft3(x3, x2, x1, R * x2 + C * x0);
      shft2(f2, f1, (*f)(x1));
	}
  }
  if (f1 < f2) {
    *xmin = x1; 
    return f1; 
  } else {
    *xmin = x2; 
    return f2; 
  }
}

double dilog(double x) {
  return flexiblesusy::dilog(x);
}

Complex dilog(const Complex& x) {
  return flexiblesusy::dilog(x);
}

double fps(double z) {
  if (z < 0.25) {
    double y = std::sqrt(1.0 - 4.0 * z);
    return 2.0 * z / y * 
      (dilog(1.0 - (1.0 - y) / (2.0 * z)) - 
       dilog(1.0 - (1.0 + y) / (2.0 * z)));
  }
  Complex zz(z);
  Complex y = sqrt(1.0 - 4.0 * zz);

  Complex ans = 2.0 * zz / y * 
    (dilog(1.0 - (1.0 - y) / (2.0 * zz)) - 
     dilog(1.0 - (1.0 + y) / (2.0 * zz)));

  /// answer should always be real
  if (ans.imag() > EPSTOL) throw flexiblesusy::SoftsusyNumericsError("Error in fps");
  return ans.real();
}

double fs(double z) {
  return (2.0 * z - 1) * fps(z) - 2.0 * z * (2.0 + std::log(z));
}

double ffbar(double z) {
  return z * 0.5 * (2.0 + std::log(z) - fps(z));
}

DoubleMatrix display3x3RealMixing(double theta12, double theta13, 
				  double theta23, double d) {
  Complex eID(cos(d), sin(d));
  double s12 = sin(theta12);
  double s13 = sin(theta13);
  double s23 = sin(theta23);
  
  double c12 = cos(theta12);
  double c13 = cos(theta13);
  double c23 = cos(theta23);

  DoubleMatrix ckmMatrix(3, 3);
  ckmMatrix(1, 1) = c12 * c13;      
  ckmMatrix(1, 2) = s12 * c13; 

  /// phase factor e^i delta: we'll set it to + or - 1 depending on the sign
  /// of s13
  int pf = 1;
  if (s13 < 0.) pf = -1;
  
  ckmMatrix(1, 3) = pf * s13;
  ckmMatrix(2, 1) = (-s12 * c23 - pf * c12 * s23 * s13);
  ckmMatrix(2, 2) = (c12 * c23 - pf * s12 * s23 * s13);
  ckmMatrix(2, 3) = s23 * c13; 
  ckmMatrix(3, 1) = (s12 * s23 - pf * c12 * c23 * s13); 
  ckmMatrix(3, 2) = (-c12 * s23 - pf * s12 * c23 * s13); 
  ckmMatrix(3, 3) = c23 * c13;

  return ckmMatrix;
  
}

double den(double a, int b) {
  double aa = a;
  for (int i=1; i<b; i++) aa = aa * a;
  return 1. / aa;
}

double fin(double mm1, double mm2)
{
  using std::log;

  if (mm1>mm2)
    return (-3.5 - (7.*mm2)/(2.*mm1) +
	    dilog(mm2/mm1) - (mm2*dilog(mm2/mm1))/mm1 -
	    (3.*mm2*log(mm1))/mm1 - log(mm1)*log(mm1 - mm2) +
	    (mm2*log(mm1)*log(mm1 - mm2))/mm1 + (3.*mm2*log(mm2))/mm1 -
	    log(mm1)*log(mm2) + (2.*mm2*log(mm1)*log(mm2))/mm1 +
	    log(mm1 - mm2)*log(mm2) - (mm2*log(mm1 - mm2)*log(mm2))/mm1 - 
	    sqr(PI) * 0.25 +
	    (mm2*sqr(PI))/(12.*mm1) + sqr(log(mm1)) - 
	    (3.*mm2*sqr(log(mm1)))/(2.*mm1) -
	    (mm2*sqr(log(mm2)))/(2.*mm1));
  else if (mm1<mm2)
    return (-3.5 - (7.*mm2)/(2.*mm1) -
	    dilog(mm1/mm2) + (mm2*dilog(mm1/mm2))/mm1 -
	    (3*mm2*log(mm1))/mm1 + (3*mm2*log(mm2))/mm1 +
	    (mm2*log(mm1)*log(mm2))/mm1 - log(mm1)*
	    log(-mm1 + mm2)+ (mm2*log(mm1)*log(-mm1 + mm2))/mm1 
	    + log(mm2)*log(-mm1 + mm2) -
	    (mm2*log(mm2)*log(-mm1 + mm2))/mm1 + sqr(PI) / 12. -
	    (mm2*sqr(PI))/(4.*mm1) + sqr(log(mm1)) * 0.5 - 
	    (mm2*sqr(log(mm1)))/mm1 -
	    sqr(log(mm2))/2);
  else return 7.-sqr(PI)/6.;  
}


void rotate(DoubleMatrix & r, DoubleMatrix & qt, int n, int i, float a, 
	    float b)
{
  using std::sqrt;

  int j;
  double c,fact,s,w,y;
  
  if (a == 0.0) {
    c = 0.0;
    s = (b >= 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) {
    fact = b / a;
    c = sign(1.0/sqrt(1.0+(fact*fact)),a);
    s = fact*c;
  } else {
    fact = a / b;
    s = sign(1.0 / sqrt(1.0 + (fact * fact)), b);
    c = fact * s;
  }
  for (j=i; j<=n; j++) {
    y = r(i, j);
    w = r(i+1, j);
    r(i, j) = c * y - s * w;
    r(i+1, j) = s * y + c * w;
  }
  for (j=1; j<=n; j++) {
    y = qt(i, j);
    w = qt(i+1, j);
    qt(i, j) = c * y - s * w;
    qt(i+1, j) = s * y + c * w;
  }
}

} // namespace softsusy
