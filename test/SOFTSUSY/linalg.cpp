
/** \file linalg.cpp 
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach 
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
*/

#include "linalg.h"
#include "utils.h"
#include "def.h"
#include "error.hpp"
#include "diagonalization.hpp"
#include <cmath>
using std::abs;
using std::cos;
using std::sin;
#include <sstream>
using std::ostringstream;
#include <istream>
using std::istream;
using std::ostream;

namespace softsusy {

/************************************
 * 
 *   *** DOUBLE VECTOR CLASS ***
 *
 ************************************/

void DoubleVector::append(const DoubleVector& v)
{
   const size_t old_size = x.size();
   std::valarray<double> y(x); // old vector
   x.resize(old_size + v.size());
   end = start + x.size() - 1;
   for (size_t i = 0; i < y.size(); ++i)
      x[i] = y[i];
   for (size_t i = 0; i < v.size(); ++i)
      x[i + old_size] = v(i + v.displayStart());
}

/*
 *  ELEMENT ACCESS
 */

// Changes the length of a vector - copies as many elements of old one as
// possible, and fills any extra up with zeroes
void DoubleVector::setEnd(int e) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (e < start) {
    ostringstream ii;
    ii << "DoubleVector.setEnd called with new end " << e << " before start " << start;
    ii << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  if (e != end) {
    std::valarray<double> t(0.0, e-start+1);

    if (e < end) {
      t = x[std::slice(0,t.size(),1)];
    } else {
      t[std::slice(0,x.size(),1)] = x;
    }

    x.resize(t.size(),0.0);
    x = t;
    end = e;
  }
}

// Most negative / smallest element
double DoubleVector::min(int & p) const {
  const double & m = x.min();
  for (int i=displayStart(); i<=displayEnd();  ++i)
    if (display(i) == m) {
      p = i;
      break;
    }
  return m;
}

void DoubleVector::swap(int i, int j) {
  double m = display(j);
  (*this)(j) = display(i);
  (*this)(i) = m;
}

/*
 *  STANDARD INPUT / OUTPUT
 */

namespace { // NUMDisplay to be used in this file only
  const int NUMDisplay = 5;
}

ostream & operator <<(ostream &left, const DoubleVector &v) {
  left << "(" << v.displayStart() << "," << v.displayEnd() << "):\n";
  for(int i=v.displayStart(); i<=v.displayEnd(); ++i) {
    if (v.display(i) >= 0.0) left << " ";
    left << v.display(i);
    if ((i % NUMDisplay) == 0 || (i == v.displayEnd())) 
      left << '\n';
    else
      left << " ";
  }
  return left;
}

istream & operator>>(istream & left, DoubleVector &v) {
  DoubleVector temp(v.displayStart(), v.displayEnd());
  for (int i=temp.displayStart(); i<=temp.displayEnd(); ++i) {
    left >> temp(i); 
  }
  v = temp;
  return left;
}

int DoubleVector::closest(double a) const {
  double dist = 6.66e66; 
  int position = -666;
  for (int i = displayStart(); i <= displayEnd(); ++i) {
    double propDist = fabs(display(i) - a);
    if (propDist < dist) {
      dist = propDist;
      position = i;
    }
  }
  return position;
}

void DoubleVector::fillArray(double* array, unsigned offset) const
{
   softsusy::fillArray(x, array, offset);
}







/************************************
 * 
 *   *** DOUBLE MATRIX CLASS ***
 *
 ************************************/

/*
 *  CONSTRUCTORS
 */

// Makes diagonal square matrix out of v
DoubleMatrix::DoubleMatrix(const DoubleVector &v)
  : rows(v.displayEnd()), cols(v.displayEnd()), x(0.0, v.displayEnd() * v.displayEnd()) { 
  diag() = v.x;
}

/*
 *  ELEMENT ACCESS
 */

const DoubleMatrix & DoubleMatrix::operator=(double v) {
  // square matrix?
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols)  {
    ostringstream ii;
    ii << "DoubleMatrix = double overload; not square\n";
    ii << *this << "=\n" << v;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  std::valarray<double> diagvec(v,rows);
  x = 0.0; // clear all matrix entries
  diag() = diagvec;
  return *this;
}

double DoubleMatrix::min(int & k, int & l) const {
  const double & m = x.min();
  for (int i=1; i<=displayRows();  ++i)
    for (int j=1; j<=displayCols(); ++j)
      if (display(i, j) == m) {
	k = i; l = j;
	break;
      }
  return m;
}

double DoubleMatrix::max(int & k, int & l) const {
  const double & m = x.max();
  for (int i=1; i<=displayRows();  ++i)
    for (int j=1; j<=displayCols(); ++j)
      if (display(i, j) == m) {
	k = i; l = j;
	break;
      }
  return m;
}

void DoubleMatrix::swaprows(int i, int j) {
  std::valarray<double> temp(row(i));
  row(i) = row(j);
  row(j) = temp;
}

void DoubleMatrix::swapcols(int i, int j) {
  std::valarray<double> temp(col(i));
  col(i) = col(j);
  col(j) = temp;
}

void DoubleMatrix::setCols(int numberOfCols)
{
   if (numberOfCols == cols)
      return;

   std::valarray<double> old(x);
   x.resize(rows * numberOfCols, 0.0);
   if (numberOfCols > cols) {
      for (std::size_t i = 0; i < old.size(); ++i)
         x[i + (i / cols) * (numberOfCols - cols)] = old[i];
   } else {
      for (std::size_t i = 0; i < x.size(); ++i)
         x[i] = old[i + (i / numberOfCols) * (cols - numberOfCols)];
   }
   cols = numberOfCols;
}

void DoubleMatrix::setRows(int numberOfRows)
{
   if (numberOfRows == rows)
      return;

   std::valarray<double> old(x);
   x.resize(cols * numberOfRows, 0.0);
   const std::size_t minSize = std::min(old.size(), x.size());
   for (std::size_t i = 0; i < minSize; ++i)
      x[i] = old[i];
   rows = numberOfRows;
}

void DoubleMatrix::resize(int numberOfRows, int numberOfCols)
{
   setCols(numberOfCols);
   setRows(numberOfRows);
}

double DoubleMatrix::trace() const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols)  {
    ostringstream ii;
    ii << "trace of non-square DoubleMatrix\n" << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  double sum = 0.0;
  for (int i = 1; i<=rows; ++i)
    sum += elmt(i,i);
  return sum;
}

DoubleMatrix DoubleMatrix::transpose() const {
  DoubleMatrix temp(cols,rows);
  for (int i=1; i<=rows; ++i)
     for (int k=1; k<=cols; ++k)
        temp(k,i) = elmt(i,k);
  return temp;
}

bool DoubleMatrix::testNan() const {
  for (int i=1; i<=rows; i++)
    for (int j=1; j<=cols; j++)
      if (display(i, j) != display(i, j)) return true;
  return false;
}

DoubleVector DoubleMatrix::flatten() const {
   DoubleVector v(x.size());
   for (size_t i = 0; i < x.size(); ++i)
      v(i + 1) = x[i];
   return v;
}

/*
 *  NUMERICAL DIAGONALIZATION ROUTINES ETC.
 */

// Perform on a mixing matrix O such that
// diag = O^double M O. It will put the diag evals in abs ascending order amd
// change O accordingly
void DoubleMatrix::associateOrderAbs(DoubleVector &v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if ((displayRows() != v.displayEnd()) ||
      (displayCols() != v.displayEnd()) ||
      (v.displayStart() != 1)) {
    ostringstream ii;
    ii << "Associated ordering incompatibility\n";
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  for (int i=v.displayStart(); i<=v.displayEnd(); ++i)
    for (int j=i+1; j<=v.displayEnd(); ++j)
      if (abs(v.display(i)) > abs(v.display(j))) {
	v.swap(i,j);
	swapcols(i,j);
      }
}

// Perform such that
// diag = U^T double M V. It will put the diag evals in abs ascending order amd
// change U,V accordingly
void DoubleMatrix::associateOrderAbs(DoubleMatrix & u, 
				     DoubleMatrix & v, 
				     DoubleVector &w) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if ((u.displayRows() != w.displayEnd()) ||
      (u.displayCols() != w.displayEnd()) ||
      (w.displayStart() != 1)) {
    ostringstream ii;
    ii << "Associated ordering incompatibility\n";
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  for (int i=w.displayStart(); i<=w.displayEnd(); ++i)
    for (int j=i+1; j<=w.displayEnd(); ++j)
      if (abs(w.display(i)) > abs(w.display(j))) {
	w.swap(i, j);
	v.swapcols(i, j);
	u.swapcols(i, j);
      }
}


void DoubleMatrix::symmetrise() {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "Error: symmetrising rectangular matrix " << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  for (int i=1; i<=rows; ++i)
     for (int k = 1; k < i && k <= cols; ++k)
        elmt(i,k) = elmt(k,i);
}

// Gives sum of difference between two matrices
double DoubleMatrix::compare(const DoubleMatrix & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != a.displayRows() || displayCols() !=
      a.displayCols()) {
    ostringstream ii;
    ii << "Error: comparing matrices of different sizes" << *this << 
      " and " << a;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  return abs(x - a.x).sum();
}

double DoubleMatrix::compare(const ComplexMatrix & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != a.displayRows() || displayCols() !=
      a.displayCols()) {
    ostringstream ii;
    ii << "Error: comparing matrices of different sizes" << *this << 
      " and " << a;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  ComplexMatrix temp(*this);
  return abs(abs(a.x) - abs(temp.x)).sum().real();
}

//Returns diagonal values ONLY of a matrix
DoubleVector DoubleMatrix::diagVals() const {
  DoubleVector temp(1, displayCols());//what's this for?   W(1, displayCols());
  DoubleMatrix U(displayRows(), displayCols());
  DoubleMatrix V(displayRows(), displayCols());
  diagonalise(U, V, temp);
  return temp;
}

// Returns inverse of a matrix
DoubleMatrix DoubleMatrix::inverse() const {
  DoubleMatrix ans(*this);
#ifdef ARRAY_BOUNDS_CHECKING
  if (ans.displayCols() != ans.displayRows()) { 
    ostringstream ii;
    ii << "Error in DoubleMatrix::invert, trying to invert non-square";
    ii << " matrix " << *this << '\n';
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  DoubleMatrix u(ans.displayRows(), ans.displayRows());
  DoubleMatrix v(ans.displayRows(), ans.displayRows());
  DoubleVector w(ans.displayRows());
  ans.diagonalise(u, v, w);
  
  const double underflow = 1.0e-16;
  for (int i=1; i<=ans.displayRows(); ++i) {
    if (w(i) > underflow) w(i) = 1 / w(i); // Leave zero eigenvalues alone
    else {
      ostringstream ii;
      ii << "Problem in DoubleMatrix::inverse, " << *this 
	 << "non-invertible\n";
      throw flexiblesusy::DiagonalizationError(ii.str());
    }
  }
  DoubleMatrix temp(w);
  return DoubleMatrix(v * temp * u.transpose());
}

// A = U.W.V^T where W is a matrix of the eigenvalues
// therefore W = U^T A V. Uses SVD but is not yet set up for really
// singular matrices
// You should really set up more/faster routines for symmetric case etc
double DoubleMatrix::diagonalise(DoubleMatrix & u, DoubleMatrix & v, 
				 DoubleVector & w) const {
#ifdef ARRAY_BOUNDS_CHECKING
  const int c = displayCols();
  if (displayRows() != c || w.displayEnd() != c || v.displayCols() != c
      || v.displayRows() !=c || u.displayCols() !=c || u.displayRows()
      !=c) {
    ostringstream ii;
    ii << "Error: Trying to diagonalise matrix \n" << *this 
       << "with u" << u << "v " << v << "w " << w;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  // Numerical routine replaces argument, so make a copy of elements
  u = *this;
  diagonaliseSvd(u, w, v);
  
  // Phase freedom - tends to give more familiar matrices
  u = -1.0 * u; v = -1.0 * v;
  
  associateOrderAbs(u, v, w);
  // Calculate residue: difference in two matrices, divided by absolute value
  // of largest eigenvalue
  DoubleMatrix temp(w);
  return compare(DoubleMatrix(u * temp * v.transpose())) / w.max();
}

// For SYMMETRIC MATRICES ONLY!
// A = V.W.V^T where W is a matrix of the eigenvalue therefore W =
// V^T A V. 
double DoubleMatrix::diagonaliseSym(DoubleMatrix & v, DoubleVector & w) const {
  const int c = displayCols();
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != c || w.displayEnd() != c || v.displayCols() != c
      || v.displayRows() !=c || c > 10) {
    ostringstream ii;
    ii << "Error: Trying to diagonalise matrix \n" << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  // Numerical recipes routine replaces argument, so make a copy of elements
  DoubleMatrix a(*this);
  int nrot;
  diagonaliseJac(a, c, w, v, &nrot);
  
  // So Far, V is such that the eigenvalues are in some random order.
  // We must now re-order the rows of V to get the eigenvalues in the
  // correct order (taken to be increasing abs(w)).
  v.associateOrderAbs(w);
  
  // Calculate residue: difference in two matrices, divided by absolute
  // value of largest eigenvalue
  DoubleMatrix temp(w); 
  return compare(DoubleMatrix(v * temp * v.transpose())) / w.max();
}

// For SYMMETRIC MATRICES ONLY!
// A = V.W.V^double where W is a matrix of the eigenvalue therefore W =
// V^double A V. all diagonal values are made positive!
double DoubleMatrix::diagonaliseSym(ComplexMatrix & v, DoubleVector & w) const {
  const int c = displayCols();
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != c || w.displayEnd() != c || v.displayCols() != c
      || v.displayRows() !=c || c > 10) {
    ostringstream ii;
    ii << "Error: Trying to diagonalise matrix \n" << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  // Numerical recipes routine replaces argument, so make a copy of elements
  DoubleMatrix a(*this), k(c, c);
  int nrot;
  diagonaliseJac(a, c, w, k, &nrot);
  
  // So Far, V is such that the eigenvalues are in some random order.
  // We must now re-order the rows of V to get the eigenvalues in the
  // correct order (taken to be increasing abs(w)).
  k.associateOrderAbs(w);
  
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  ComplexMatrix m(c, c);

  for (int i=1; i<=c; ++i) 
    if (w(i) < 0.0) {
      for (int j=1; j<=c; ++j) v(i, j) = Complex(0.0, k(i, j));
      w(i) = - w(i);
    }
    else {
      for (int j=1; j<=c; ++j) v(i, j) = Complex(k(i, j), 0.0);
    }
  
  // Calculate residue: difference in two matrices, divided by absolute
  // value of largest eigenvalue
  ComplexMatrix temp(c, c); 
  for (int i=1; i<=c; ++i) temp(i, i) = w(i);
  //this seems to be a relic
  //cout << v * temp * v.transpose();
  return compare(ComplexMatrix(v.transpose() * temp * v)) / w.max();
}

// Special case that's often used: eigenvalues of 2 by 2 symmetric matrix
// A(2,1) is assumed to be equal to A(1,2). theta is the mixing angle:
// [ cos theta   sin theta   ]    A    [cos theta   -sin theta ]  = [ m1 0 ]
// [ -sin theta  cos theta   ]         [sin theta    cos theta ]    [ 0 m2 ]
// where m1 is NOT NECESSARILY > m2
DoubleVector DoubleMatrix::sym2by2(double & theta) const  {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() !=2 || displayCols() !=2) {
    ostringstream ii;
    ii << "Called sym2by2 with matrix of dimension ("
       << displayRows() << ", " << displayCols() << ")\n";
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
  if (display(1, 2) != display(2, 1)) {
    ostringstream ii;
    ii << "Called non symmetric sym2by2:" << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  
  DoubleVector temp(1, 2);
  double sumTol = abs(abs(display(1, 1)) - abs(display(2, 2))) /
    maximum(abs(display(1, 1)), abs(display(2, 2)));
  
  double t2t = 2.0 * display(1, 2) / (display(1, 1) - display(2, 2));
  
  if (sumTol > EPSTOL || t2t < EPSTOL) theta = atan(t2t) * 0.5; 
  else {
    if (t2t > 0.0) 
      theta = PI * 0.25;
    else 
      theta = -PI * 0.25;
  }
  
  DoubleMatrix mm(rot2d(theta) * (*this) * rot2d(theta).transpose());
  
  temp(1) = mm(1, 1);
  temp(2) = mm(2, 2);
  
#ifdef ARRAY_BOUNDS_CHECKING
  double maxtol = maximum(abs(temp(1)), abs(temp(2))) * TOLERANCE;
  if (abs(mm(2, 1)) > maxtol || abs(mm(1, 2)) > maxtol) {
    ostringstream ii;
    ii << "Inaccurate diagonalisation in LINALG::sym2by2\n";
    ii << "Diagonalising " << * this;
    ii << "Found diagonalised matrix to be " << mm;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  
  return temp; 
}

// 2 by 2 asymmetric matrices
// [ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ]  = diag
// [ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
DoubleVector DoubleMatrix::asy2by2(double & thetaL, double & thetaR) const {
   DoubleVector temp1(2), temp2(2);
  DoubleMatrix e(*this);
  
  DoubleMatrix h(e.transpose() * e);
  temp1 = (h.sym2by2(thetaR)).apply(zeroSqrt);
  if (temp1(1) > temp1(2))
    thetaR = thetaR + PI * 0.5;
  
  h = e * e.transpose();
  temp2 = (h.sym2by2(thetaL)).apply(zeroSqrt);
  
  // Did the eigenvalues come out in the wrong order? 
  // If so, swap the order of the second case by adding pi/2 to the angle.
  if (temp2(1) > temp2(2)) {
    thetaL = thetaL + PI * 0.5;
  }
  
  DoubleMatrix mm(2, 2);
  mm = rot2d(thetaL) * e * rot2d(-thetaR);
  
#ifdef ARRAY_BOUNDS_CHECKING
  double maxtol = maximum(abs(temp1(1)), abs(temp1(2))) * TOLERANCE;
  if (abs(mm(2, 1)) > maxtol || abs(mm(1, 2)) > maxtol) {
    outputCharacteristics(15);
    ostringstream ii;
    ii << "Inaccurate diagonalisation in LINALG::asy2by2\n";
    ii << "temp1=" << temp1 << " temp2=" << temp2 << '\n';
    ii << "diagonalising " << e << " where thetaR=" << thetaR 
       << " thetaL=" << thetaL << '\n';
    ii << "m=" << mm;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  
  temp1(1) = mm(1, 1); temp1(2) = mm(2, 2);
  
  return temp1;
}

/*
 *  STANDARD INPUT / OUTPUT
 */

istream & operator>>(istream & left, DoubleMatrix &m) {
  DoubleMatrix x(m.displayRows(), m.displayCols());
  for (int i=1; i<=m.displayRows(); ++i)
    for (int j=1; j<=m.displayCols(); ++j)
      left >> x(i, j);
  m = x;  
  return left;
}

ostream & operator <<(ostream &left, const DoubleMatrix &v) {
  static const double underflow = 1.0e-120;
  double value;
  
  left << "(" << v.displayRows() << "," <<
    v.displayCols() << "):\n";
  for(int i=1; i<=v.displayRows(); ++i) {
    for(int j=1; j<=v.displayCols(); ++j) {
      value = v.display(i, j);
      if (abs(value) < underflow) value = 0.0; // Traps -0.0
      if (value >= 0.0) left << " "; 
      left << value << " ";
    }
    left << '\n';
  }
  return left;
}

/*
 *  NON-MEMBER DIAGONALIZATION ROUTINES ETC.
 */

// 2 dimensional rotation matrix M
// Returns U = [ cos theta    sin theta ]
//             [-sin theta    cos theta ]
DoubleMatrix rot2d(double theta) {
  DoubleMatrix u(2, 2);  
  u(1, 1) = cos(theta); u(2, 2) = u(1, 1);
  u(1, 2) = sin(theta); u(2, 1) = - u(1, 2);
  return u;
}

// Returns mixing matrix
// [ -sin(theta)  cos(theta) ]
// [  cos(theta)  sin(theta) ] --
// takes into account different conventions for the ordering of a.
DoubleMatrix rot2dTwist(double theta) {
  DoubleMatrix n(2, 2);
  n(1, 1) = -sin(theta);
  n(1, 2) = cos(theta);
  n(2, 1) = cos(theta); 
  n(2, 2) = sin(theta);
  return n;
}

// LCT: 3 dimensional rotation matrix U
// Returns U = [ cos theta  sin theta   0 ]
//             [ sin theta  -cos theta  0 ]
//             [ 0          0           1 ]
DoubleMatrix rot3d(double theta) {
  DoubleMatrix u(3, 3);
  u(1, 1) = -cos(theta); u(2, 2) = -u(1, 1);
  u(1, 2) = sin(theta);  u(2, 1) = u(1, 2);
  u(3, 3) = 1.0;
  return u;
}

// Redefines mixing matrices to be complex such that diagonal values are
// positive for a 2 by 2: if  
// [ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ]  = diag
// [ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
// as given by asy2by2!
// then u^* A v^+ = mdiagpositive
// IT'S BEEN CHECKED!
void positivise(double thetaL, double thetaR, const DoubleVector & diag,
		ComplexMatrix & u, ComplexMatrix & v) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (diag.displayStart()!=1 || diag.displayEnd()!=2 ||
      u.displayRows()!=2 || u.displayCols()!=2 ||
      v.displayRows()!=2 || v.displayCols()!=2) {
    ostringstream ii;
    ii << "DoubleMatrix::positivise currentl only available for 2 by 2"
       << "matrices, not " << diag << u << v;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  
  ComplexMatrix a(2, 2);
  for (int i=1; i<=2; ++i)
     if (diag.display(i) < 0.0)
        a(i, i) = Complex(0.0, 1.0);
     else
        a(i, i) = Complex(1.0, 0.0);
  
  u = a * rot2d(thetaL);
  v = a * rot2d(thetaR);
  return;
}

// ---------------- diagonalisation algorithms --------------------
void diagonaliseSvd(DoubleMatrix & a, DoubleVector & w, DoubleMatrix & v) {
  if (a.testNan()) {
    ostringstream ii;
    ii << "Nans present in linalg.cpp:diagonaliseSvd: diagonalising\n" 
       << a << std::endl;
    throw flexiblesusy::DiagonalizationError(ii.str());
  }
  
  int n = a.displayCols(), m = a.displayRows();
  int flag, i, its, j, jj, k, l = 0, nm = 0; 
  double anorm, c, f, g, h, s, scale, x, y, z; 
  
  DoubleVector rv1(a.displayCols());
  
  g = scale = anorm = 0.0; 
  for (i = 1; i<= n; ++i) {
    l = i + 1; 
    rv1(i) = scale * g; 
    g = s = scale = 0.0; 
    if (i <=  m) {
      for (k = i; k <= m; k++) scale += abs(a(k, i)); 
      if (scale) {
	for (k = i; k <= m; k++) {
	  a(k, i) /=  scale; 
	  s +=  a(k, i) * a(k, i); 
	}
	f = a(i, i); 
	g = -sign(std::sqrt(s), f); 
	h = f * g - s; 
	a(i, i) = f-g; 
	for (j = l; j <= n; ++j) {
	  for (s = 0.0, k = i; k <= m; k++) s +=  a(k, i) * a(k, j); 
	  f = s / h; 
	  for (k = i; k <= m; k++) a(k, j) +=  f * a(k, i); 
	}
	for (k = i; k <= m; k++) a(k, i)  *=  scale; 
      }
    }
    w(i) = scale * g; 
    g = s = scale = 0.0; 
    if (i  <=  m && i !=  n) {
      for (k = l; k <= n; k++) scale += abs(a(i, k)); 
      if (scale) {
	for (k = l; k <= n; k++) {
	  a(i, k) /=  scale; 
	  s +=  a(i, k) * a(i, k); 
	}
	f = a(i, l); 
	g  =  -sign(std::sqrt(s), f); 
	h = f * g-s; 
	a(i, l) = f - g; 
	for (k = l; k <= n; k++) rv1(k) = a(i, k) / h; 
	for (j = l; j <= m; ++j) {
	  for (s = 0.0, k = l; k <= n; k++) s +=  a(j, k) * a(i, k); 
	  for (k = l; k <= n; k++) a(j, k) +=  s * rv1(k); 
	}
	for (k = l; k <= n; k++) a(i, k)  *=  scale; 
      }
    }
    anorm = maximum(anorm, (abs(w(i)) + abs(rv1(i)))); 
  }
  for (i = n; i >= 1; i--) {
    if (i < n) {
      if (g) {
	for (j = l; j <= n; ++j)
	  v(j, i) = (a(i, j) / a(i, l)) / g; 
	for (j = l; j <= n; ++j) {
	  for (s = 0.0, k = l; k <= n; k++) s += a(i, k) * v(k, j); 
	  for (k = l; k <= n; k++) v(k, j) += s * v(k, i); 
	}
      }
      for (j = l; j <= n; ++j) v(i, j) = v(j, i) = 0.0; 
    }
    v(i, i) = 1.0; 
    g = rv1(i); 
    l = i; 
  }
  for (i = minimum(m, n); i >= 1; i--) {
    l = i + 1; 
    g = w(i); 
    for (j = l; j <= n; ++j) a(i, j) = 0.0; 
    if (g) {
      g = 1.0 / g; 
      for (j = l; j <= n; ++j) {
	for (s = 0.0, k = l; k <= m; k++) s +=  a(k, i) * a(k, j); 
	f = (s / a(i, i)) * g; 
	for (k = i; k <= m; k++) a(k, j) +=  f * a(k, i); 
      }
      for (j = i; j <= m; ++j) a(j, i) *=  g; 
    } else for (j = i; j <= m; ++j) a(j, i) = 0.0; 
    ++a(i, i); 
  }
  for (k = n; k >= 1; k--) {
    for (its = 1; its <= 30; its++) {
      flag = 1; 
      for (l = k; l >= 1; l--) {
	nm = l - 1; 
	if (abs(rv1(l))+anorm  ==  anorm) {
	  flag = 0; 
	  break; 
	}
	if (abs(w(nm))+anorm  ==  anorm) break; 
      }
      if (flag) {
	c = 0.0; 
	s = 1.0; 
	for (i = l; i <= k; ++i) {
	  f = s * rv1(i); 
	  rv1(i) = c * rv1(i); 
	  if (abs(f)+anorm  ==  anorm) break; 
	  g = w(i); 
	  h = pythagoras(f, g); 
	  w(i) = h; 
	  h = 1.0 / h; 
	  c = g * h; 
	  s  =  -f * h; 
	  for (j = 1; j <= m; ++j) {
	    y = a(j, nm); 
	    z = a(j, i); 
	    a(j, nm) = y * c + z * s; 
	    a(j, i) = z * c - y * s; 
	  }
	}
      }
      z = w(k); 
      if (l  ==  k) {
	if (z < 0.0) {
	  w(k) = -z; 
	  for (j = 1; j <= n; ++j) v(j, k) = -v(j, k); 
	}
	break; 
      }
      if (its  ==  30 || l < 1 || k < 1) 
	throw flexiblesusy::DiagonalizationError("no convergence in 30 diagonaliseSvd iterations");
      
      x = w(l); 
      nm = k - 1; 
      y = w(nm); 
      g = rv1(nm); 
      h = rv1(k); 
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y); 
      g = pythagoras(f, 1.0); 
      f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x; 
      c = s = 1.0; 
      for (j = l; j <= nm; ++j) {
	i = j+1; 
	g = rv1(i); 
	y = w(i); 
	h = s * g; 
	g = c * g; 
	z = pythagoras(f, h); 
	rv1(j) = z; 
	c = f / z; 
	s = h / z; 
	f = x * c + g * s; 
	g  =  g * c - x * s; 
	h = y * s; 
	y  *=  c; 
	for (jj = 1; jj <= n; ++jj) {
	  x = v(jj, j); 
	  z = v(jj, i); 
	  v(jj, j) = x * c + z * s; 
	  v(jj, i) = z * c - x * s; 
	}
	z = pythagoras(f, h); 
	w(j) = z; 
	if (z) {
	  z = 1.0 / z; 
	  c = f * z; 
	  s = h * z; 
	}
	f = c * g+s * y; 
	x = c * y-s * g; 
	for (jj = 1; jj <= m; ++jj) {
	  y = a(jj, j); 
	  z = a(jj, i); 
	  a(jj, j) = y * c + z * s; 
	  a(jj, i) = z * c - y * s; 
	}
      }
      rv1(l) = 0.0; 
      rv1(k) = f; 
      w(k) = x; 
    }
  }
}

double pythagoras(double a,  double b) {
  using std::sqrt;
  double absa, absb; 
  absa = abs(a); 
  absb = abs(b); 
  if (absa > absb) 
    return absa  *  sqrt(1.0 + sqr(absb / absa)); 
  else 
    return (absb  == 0.0 ? 0.0 : absb * sqrt(1.0 + sqr(absa / absb))); 
}

namespace { // helper function only
  inline void rotate(DoubleMatrix & a, int i, int j, int k, int l, 
		     double s, double tau) {
    double g = a(i, j);
    double h = a(k, l);
    a(i, j) = g - s * (h + g * tau);
    a(k, l) = h + s * (g - h * tau);
  }
}

void diagonaliseJac(DoubleMatrix & a,  int n,  DoubleVector & d,  
		    DoubleMatrix & v,  int *nrot) {
  using std::sqrt;

  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c; 
  
  DoubleVector b(n), z(n); 
  
  for (ip=1; ip<=n; ip++) {
    for (iq=1; iq<=n; iq++) v(ip,  iq) = 0.0; 
    v(ip,  ip) = 1.0; 
  }
  for (ip=1; ip<=n; ip++) {
    b(ip) = d(ip) = a(ip,  ip); 
    z(ip) = 0.0; 
  }
  *nrot = 0; 
  for (i=1; i<=50; ++i) {
    sm = 0.0; 
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++)
	sm += abs(a(ip,  iq)); 
    }
    if (sm == 0.0) return; 
    
    if (i < 4)
      tresh = 0.2  *sm/(n*n); 
    else
      tresh = 0.0; 
    for (ip=1; ip<=n-1; ip++) {
      for (iq=ip+1; iq<=n; iq++) {
	g=100.0 * abs(a(ip,  iq)); 
	if (i > 4 && abs(d(ip)) + g == abs(d(ip))
	    && abs(d(iq)) + g == abs(d(iq)))
	  a(ip,  iq) = 0.0; 
	else if (abs(a(ip,  iq)) > tresh) {
	  h = d(iq) - d(ip); 
	  if (abs(h) + g == abs(h))
	    t=(a(ip,  iq)) / h; 
	  else {
	    theta = 0.5 * h / (a(ip,  iq)); 
	    t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta)); 
	    if (theta < 0.0) t = -t; 
	  }
	  c = 1.0 / sqrt(1 + t * t); 
	  s = t * c; 
	  tau = s / (1.0 + c); 
	  h = t * a(ip,  iq); 
	  z(ip) -= h; 
	  z(iq) += h; 
	  d(ip) -= h; 
	  d(iq) += h; 
	  a(ip,  iq) = 0.0; 
	  for (j=1; j<=ip-1; ++j) rotate(a, j, ip, j, iq, s, tau); 
	  for (j=ip+1; j<=iq-1; ++j) rotate(a, ip, j, j, iq, s, tau);
	  for (j=iq+1; j<=n; ++j) rotate(a, ip, j, iq, j, s, tau);
	  for (j=1; j<=n; ++j) rotate(v, j, ip, j, iq, s, tau);
	  ++(*nrot); 
	}
      }
    }
    for (ip=1; ip<=n; ip++) {
      b(ip) += z(ip); 
      d(ip) = b(ip); 
      z(ip) = 0.0; 
    }
  }
  ostringstream ii;
  if (PRINTOUT) {
    ii << "Too many iterations in routine diagonaliseJac: diagonalising\n" 
       << a << "to" << d << " with " << v;
  }
  throw flexiblesusy::DiagonalizationError(ii.str());
}







/************************************
 * 
 *   *** COMPLEX VECTOR CLASS ***
 *
 ************************************/

/*
 *  ELEMENT ACCESS
 */

// Changes the length of a vector - copies as many elements of old one as
// possible, and fills any extra up with zeroes
void ComplexVector::setEnd(int e) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (e < start) {
    ostringstream ii;
    ii << "ComplexVector.setEnd called with incompatible length " << e << '\n'; 
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  if (e != end) {
    std::valarray<Complex> t(0.0, e-start+1);

    if (e < end) {
      t = x[std::slice(0,t.size(),1)];
    } else {
      t[std::slice(0,x.size(),1)] = x;
    }

    x.resize(t.size(),0.0);
    x = t;
    end = e;
  }
}

// Smallest absolute value of element
Complex ComplexVector::min(int & p) const {
  const Complex & m = x.min();
  for (int i=displayStart(); i<=displayEnd();  ++i)
    if (display(i) == m) {
      p = i;
      break;
    }
  return m;
}

void ComplexVector::swap(int i, int j) {
  Complex m = display(j);
  (*this)(j) = display(i);
  (*this)(i) = m;
}

/*
 *  STANDARD INPUT / OUTPUT
 */

ostream & operator <<(ostream &left, const ComplexVector &v) {
  const int NUMDisplay = 5;
  left << "(" << v.displayStart() << "," << v.displayEnd() << "):\n";
  if (v.displayEnd() > NUMDisplay) 
    left << '\n';
  else
    left << " ";
  for(int i=v.displayStart(); i<=v.displayEnd(); ++i) {
    if ((i % NUMDisplay) == 1) left << "(" << i << ")=";
    left << v.display(i);
    if ((v.displayEnd() % NUMDisplay) == 0 || (i == v.displayEnd())) 
      left << '\n';
    else
      left << " ";
  }
  return left;
}

istream & operator>>(istream & left, ComplexVector &v) {
  ComplexVector temp(v.displayStart(), v.displayEnd()); // mismatch with double version
  for (int i=v.displayStart(); i<=v.displayEnd(); ++i)
    left >> temp(i);
  v = temp;
  
  return left;
}









/************************************
 * 
 *   *** COMPLEX MATRIX CLASS ***
 *
 ************************************/

/*
 *  CONSTRUCTORS
 */

ComplexMatrix::ComplexMatrix(const DoubleMatrix & m)
  : rows(m.displayRows()), cols(m.displayCols()), x(0.0, m.x.size()) {
  for (int i=0; i < rows*cols; ++i)
    x[i] = m.x[i];
}

// Makes diagonal square matrix out of W
ComplexMatrix::ComplexMatrix(const ComplexVector &v) 
  : rows(v.displayEnd()), cols(v.displayEnd()), x(0.0, v.displayEnd() * v.displayEnd()) { 
  diag() = v.x;
}

/*
 *  ELEMENT ACCESS
 */

const ComplexMatrix & ComplexMatrix::operator=(const Complex &v) {
  // square matrix?
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "ComplexMatrix = Complex overload; must be square\n";
    ii << *this << "=\n" << v;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  std::valarray<Complex> diagvec(v,rows);
  x = 0.0; // reset all matrix entries
  diag() = diagvec;
  return *this;
}

Complex ComplexMatrix::min(int & k, int & l) const {
  const Complex & m = x.min();
  int i, j; 
  for (i=1; i<=displayRows();  ++i)
    for (j=1; j<=displayCols(); ++j)
      if (display(i, j) == m) {
	k = i; l = j;
	break;
      }
  return m;
}


void ComplexMatrix::swaprows(int i, int j) {
  std::valarray<Complex> temp(row(i));
  row(i) = row(j);
  row(j) = temp;
}

void ComplexMatrix::swapcols(int i,int j) {
  std::valarray<Complex> temp(col(i));
  col(i) = col(j);
  col(j) = temp;
}


void ComplexMatrix::setCols(int numberOfCols)
{
   if (numberOfCols == cols)
      return;

   std::valarray<Complex> old(x);
   x.resize(rows * numberOfCols, 0.0);
   if (numberOfCols > cols) {
      for (std::size_t i = 0; i < old.size(); ++i)
         x[i + (i / cols) * (numberOfCols - cols)] = old[i];
   } else {
      for (std::size_t i = 0; i < x.size(); ++i)
         x[i] = old[i + (i / numberOfCols) * (cols - numberOfCols)];
   }
   cols = numberOfCols;
}

void ComplexMatrix::setRows(int numberOfRows)
{
   if (numberOfRows == rows)
      return;

   std::valarray<Complex> old(x);
   x.resize(cols * numberOfRows, 0.0);
   const std::size_t minSize = std::min(old.size(), x.size());
   for (std::size_t i = 0; i < minSize; ++i)
      x[i] = old[i];
   rows = numberOfRows;
}


void ComplexMatrix::resize(int numberOfRows, int numberOfCols)
{
   setCols(numberOfCols);
   setRows(numberOfRows);
}


Complex ComplexMatrix::trace() const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols)  {
    ostringstream ii;
    ii << "Complexrace of non-square matrix\n" << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  Complex sum(0.0);
  for(int i=1; i<=rows; ++i) 
    sum += elmt(i, i);
  return sum;
}

ComplexMatrix ComplexMatrix::transpose() const {
  ComplexMatrix temp(cols,rows);
  for (int i=1; i<=rows; ++i)
     for (int k=1; k<=cols; ++k)
        temp(k,i) = elmt(i,k);
  return temp;
}

ComplexMatrix ComplexMatrix::hermitianConjugate() const {
  ComplexMatrix temp(cols,rows);
  for (int i=1; i<=rows; ++i)
     for (int k=1; k<=cols; ++k)
        temp(k,i) = conj(elmt(i,k));
  return temp;
}

namespace {
   Complex conjugator(const Complex& z) {
      return Complex(std::conj(static_cast<std::complex<double> >(z)));
   }
}

ComplexMatrix ComplexMatrix::complexConjugate() const { 
  return ComplexMatrix(x.apply(conjugator),rows,cols);
}

DoubleMatrix ComplexMatrix::real() const
{
   std::valarray<double> reals(x.size());
   for (std::size_t i = 0; i < x.size(); ++i)
      reals[i] = x[i].real();
   return DoubleMatrix(reals, rows, cols);
}

DoubleMatrix ComplexMatrix::imag() const
{
   std::valarray<double> imags(x.size());
   for (std::size_t i = 0; i < x.size(); ++i)
      imags[i] = x[i].imag();
   return DoubleMatrix(imags, rows, cols);
}

/*
 *  NUMERICAL DIAGONALIZATION ROUTINES ETC.
 */

void ComplexMatrix::symmetrise() {
#ifdef ARRAY_BOUNDS_CHECKING
  if (rows != cols) {
    ostringstream ii;
    ii << "Error: symmetrising rectangular matrix " << *this;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  for (int i=1; i<=rows; ++i)
     for (int k = 1; k < i && k <= cols; ++k)
        elmt(i,k) = elmt(k,i);
}

// Gives sum of difference between two matrices
double ComplexMatrix::compare(const ComplexMatrix & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayRows() != a.displayRows() || displayCols() !=
      a.displayCols()) {
    ostringstream ii;
    ii << "Error: comparing matrices of different sizes" << *this << 
      " and " << a;
    throw flexiblesusy::OutOfBoundsError(ii.str());
  }
#endif
  return abs(x - a.x).sum().real();
}

/*
 *  STANDARD INPUT / OUTPUT
 */

istream & operator>>(istream & left, ComplexMatrix &m) {
  int i, j;
  ComplexMatrix x(m.displayRows(), m.displayCols());
  for (i=1; i<=m.displayRows(); ++i)
    for (j=1; j<=m.displayCols(); ++j)
      left >> x(i, j);
  m = x;
  
  return left;
}

ostream & operator <<(ostream &left, const ComplexMatrix &v) {
  left << "(" << v.displayRows() << "," <<
    v.displayCols() << "):\n";
  for(int i=1; i<=v.displayRows(); ++i) {
    for(int j=1; j<=v.displayCols(); ++j) {
      if (v.display(i, j).real() >= 0.0)
	left << " ";  
      left << v.display(i, j) << " ";
    }
    left << '\n';
  }
  return left;
}


double DoubleVector::nmin(int & p) const {
  double m1 = display(displayStart()), m2;
  int p1 = displayStart(), p2;
  
  if ( display(displayStart()+1) < m1) {
    m2 = m1;
    p2 = p1;
    
    m1 = display(displayStart()+1);
    p1 = displayStart()+1;
  }
  else {
    m2 = display(displayStart()+1);
    p2 = displayStart()+1;
  }
  
  
  int i; 
  for (i=displayStart() + 2; i<=displayEnd();  i++){
    if (display(i) < m1) {
	m2 = m1;
	p2 = p1;
	
	m1 = display(i);
	p1 = i;
      }
    else if (display(i) < m2)  {
	m2 = display(i);
	p2 = i;
      }
  }
  
  p =  p2;
  return m2;
}


double DoubleMatrix::nmin(int & k, int & l) const {
  double m1 = display(1, 1), m2;
  int k1 = 1, l1 = 1, k2, l2;

  if (display(1, 2) < m1){
    m2 = m1;
    k2 = k1;
    l2 = l1;
    
    m1 = display(1, 2);
    k1 = 1;
    l1 = 2; 
  }
  else{
    m2 = display(1, 2);
    k2 = 1;
    l2 = 2; 
  }
  int i, j;
  for (i=1; i<=displayRows();  i++)
    for (j=1; j<=displayCols(); j++){
      if (display(i, j) < m1) {
	m2 = m1;
	k2 = k1;
	l2 = l1;

	m1 = display(i, j);
	k1 = i; l1 = j;
      }
      else if (display(i, j) < m2){
	m2 = display(i, j);
	k2 = i; l2 = j;
      }
    }
  k = k2; l = l2;  
  return m2;
}

double DoubleVector::compare(const DoubleVector & a) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (displayStart() != a.displayStart() || displayEnd() != a.displayEnd())
    throw flexiblesusy::OutOfBoundsError("In DoubleVector::compare, two vectors are of incompatible size.");
#endif
  
  DoubleVector sumTol(display());

  int i; for (i=displayStart(); i<=displayEnd(); i++) {
    sumTol(i) = fabs(1.0 - minimum(display(i), a.display(i)) 
		     / maximum(display(i), a.display(i)));
  }

  return sumTol.max();
}

void DoubleMatrix::fillArray(double* array, unsigned offset) const
{
   softsusy::fillArray(x, array, offset);
}

double DoubleVector::norm() const {
  using std::sqrt;
  double f = 0.0;
  int i; for (i=start; i<=end; i++) f += x[i] * x[i];
  return sqrt(f);				   
}

// Division of elements of vectors
DoubleVector DoubleVector::divide(DoubleVector const &v) const {
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != v.displayStart() || end != v.displayEnd())
    {
      ostringstream ii;
      ii << "DoubleVector * overload; incompatible lengths\n";
      ii << *this << "*\n" << v;
      throw flexiblesusy::OutOfBoundsError(ii.str());
    }
#endif
  DoubleVector temp(v);
  for(int i=start; i<=end; i++) temp(i) = x[i] / v.display(i);
  return temp;
}

DoubleVector DoubleMatrix::vectorfy () const {
  int length = displayCols() * displayRows();
  DoubleVector ans(length);
  int i,j; for (i=1; i<displayCols(); i++)
    for (j=1; j<displayRows(); j++)
      ans((i-1)*j + j) = display(j, i);
  return ans;
}

double DoubleVector::average() const {
  double f = 0.0;
  int i; for (i=start; i<=end; i++) f += x[i];
  return f / double(end);
}

/// fill array from valarray, starting at offset
void fillArray(const std::valarray<double>& x, double* array, unsigned offset)
{
   for (unsigned i = 0; i < x.size(); ++i)
      array[i + offset] = x[i];
}

const ComplexMatrix& ComplexMatrix::operator+=(const DoubleMatrix& other)
{
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != other.displayStart() || end != other.displayEnd())
    {
      ostringstream ii;
      ii << "ComplexMatrix::operator+(const DoubleMatrix&) incompatible lengths\n";
      ii << *this << "*\n" << other;
      throw flexiblesusy::OutOfBoundsError(ii.str());
    }
#endif
  for (std::size_t i = 0; i < x.size(); i++)
     x[i] += other.x[i];
  return *this;
}

ComplexMatrix ComplexMatrix::operator+(const DoubleMatrix& other)
{
#ifdef ARRAY_BOUNDS_CHECKING
  if (start != other.displayStart() || end != other.displayEnd())
    {
      ostringstream ii;
      ii << "ComplexMatrix::operator+=(const DoubleMatrix&) incompatible lengths\n";
      ii << *this << "*\n" << other;
      throw flexiblesusy::OutOfBoundsError(ii.str());
    }
#endif
  ComplexMatrix temp(*this);
  return temp += other;
}

} // namespace softsusy
