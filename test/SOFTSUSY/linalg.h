
/** \file linalg.h
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach 
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://projects.hepforge.org/softsusy/
    \brief DoubleVector and DoubleMatrix classes of doubles and 
    operations between them, complexified copies also. 

*/

#ifndef LINALG_H
#define LINALG_H

#include <vector>

#include "mycomplex.h"
#include "xpr-vector.h"
#include "xpr-matrix.h"
#include <exception>
#include <iosfwd>
#include <valarray>
#include <sstream>
namespace softsusy {
/************************************
 * 
 *   *** DOUBLE VECTOR CLASS ***
 *
 ************************************/

/// DoubleVector is of variable length, and contains double precision
class DoubleVector : public Indexable<double, DoubleVector> {
private:
  int start, end; ///< vector dimensions
  std::valarray<double> x; ///< valarray of values held in x
  ///Conversion from a valarray to a double vector with given start and end
  DoubleVector(const std::valarray<double> & valarr, int s, int e)
    : Indexable<double,DoubleVector>(), start(s), end(e), x(valarr) {}
  friend class DoubleMatrix;

public:

  /*
   *  CONSTRUCTORS
   */


  explicit DoubleVector(int e) ///< Default dimensions start at 1, end at e
    : Indexable<double,DoubleVector>(), start(1), end(e), x(0.0,e) {}
  DoubleVector(int s, int e) ///< Dimensions start at s, end at e
    : Indexable<double,DoubleVector>(), start(s), end(e), x(0.0,e-s+1) {}
  DoubleVector(double a, double b)
    : Indexable<double,DoubleVector>(), start(1), end(2), x(0.0,2) {
    x[0] = a; x[1] = b; }
  ~DoubleVector() {}

  /// this is the only required operator to make this class work with ETs
  template<class E>
  DoubleVector & operator=(const Xpr<double,E> & x) {
    *this = copy_from(x);
    return *this;
  }

  DoubleVector& operator=(const DoubleVector& other) {
    x.resize(other.x.size());
    x = other.x;
    start = other.start;
    end = other.end;
    return *this;
  }

  DoubleVector& operator=(double value) {
    for (std::size_t i = 0; i < x.size(); ++i)
      x[i] = value;
    return *this;
  }

  template <typename E>
  DoubleVector(const Xpr<double,E> & v)
    : Indexable<double,DoubleVector>(), 
      start(v.displayStart()), end(v.displayEnd()), 
      x(0.0,v.displayEnd()-v.displayStart()+1)
  {
    assign_from(v);
  }

  void append(const DoubleVector&); ///< append another vector

  /*
   *  ELEMENT ACCESS
   */

  double & operator() (int i); ///< Reference a single element
  double operator() (int i) const { return x[i-start]; }
  double display(int i) const { return x[i-start]; }///< display one element
  void set(int i, double f);///< sets ith element

  int displayStart() const { return start; }///< start of dimension
  int displayEnd() const { return end; }///< returns end of dimension
  const DoubleVector & display() const { return *this; }///< returns whole thing
  double compare(const DoubleVector & a) const; 

  /// Changes the length of a vector - copies as many elements of old one as
  /// possible, and fills any extra up with zeroes
  void setEnd(int e);


  double norm() const;
  /// applies fn to every element in a vector
  DoubleVector apply(double (*fn)(double)) const { 
    return DoubleVector(x.apply(fn),start,end); 
  }
  double nmin(int & p) const;
  double min() const { return x.min(); } ///< minimum element in vector
  double max() const { return x.max(); }///< maximum element in vector

  double max(int & p) const {
    double m = double(0);
    int i; 
    for (i=displayStart(); i<=displayEnd();  i++)
      if (display(i) > m) {
	m = display(i);
	p = i;
      }
    return m;
  }

  DoubleVector abs() const {
    return DoubleVector(std::abs(x), start, end);
  }

  double absmax(int & p) const {
    return abs().max(p);
  }

  double min(int & p) const;///< minimum element in vector
  double sum() const { return x.sum(); }
  /// Returns sum of absolute values of all elements
  double sumElements() const { return std::abs(x).sum(); }  
  std::size_t size() const { return x.size(); }
  
  void swap(int i, int j); ///< swaps ith and jth elements of a vector
  
  double dot(const DoubleVector & v) const { ///< dot product
    return (x * v.x).sum();
  }

  DoubleVector sort() const {
    DoubleVector temp(*this), saveIt(*this);
    
    int j, pos;
    for (j=start; j<=end; j++) {
      double d = temp.min(pos);
      if (d < 1.0e66) saveIt(j) = d;
      
      temp(pos) = 6.66e66; // take that element out of the min calculation
    }
    return saveIt;
  }
  
  /// Acts on each, element by element, dividing each element by v(i)
  DoubleVector divide(DoubleVector const &v) const;
  /// Average element of the vector
  double average() const; ///< average element
  /// Closest element of the vector to the input
  int closest(double a) const; 
  /// fill array, starting at offset
  void fillArray(double* array, unsigned offset = 0) const;
  bool operator==(const DoubleVector& o) const {
     if (displayStart() != o.displayStart())
        return false;
     if (displayEnd() != o.displayEnd())
        return false;
     for (int j=start; j<=end; ++j)
        if (display(j) != o.display(j))
           return false;
     return true;
  }
};


/*
 *  STANDARD INPUT / OUTPUT
 */

/// prints a vector out formatted, maximum of 5 elements per line
std::ostream & operator<<(std::ostream &left, const DoubleVector &V);
/// inputs a vector
std::istream & operator>>(std::istream & left, DoubleVector &V);


/*
 *  INLINE FUNCTION DEFINITIONS
 */

inline double & DoubleVector::operator() (int i) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i>end || i<start) {
      std::ostringstream ii;
      ii << "Trying to access " << i << "th element of DoubleVector\n";
      ii << "start " << start << " end " << end;
      ii << *this;
      throw std::runtime_error(ii.str());
  }
#endif
  return x[i-start];
}

inline void DoubleVector::set(int i, double f) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < start || i > end) 
    {
      std::ostringstream ii;
      ii << "Cannot access " << i << "th element of vector " << *this;
      throw std::runtime_error(ii.str());
    }
#endif
  x[i-start] = f;
}


/************************************
 * 
 *   *** DOUBLE MATRIX CLASS ***
 *
 ************************************/


class ComplexVector; class ComplexMatrix;
/// Matrix from 1..rows, 1..cols of double values
class DoubleMatrix : public MatIndexable<double, DoubleMatrix> {
private:
  int rows, cols;
  std::valarray<double> x; ///< Matrix is stored sequentially, row after row
  
  /// for easy access to single rows and columns
  /// not very fast; for speedy access use for-loops over the elements
  // for usage see e.g. http://www.roguewave.com/support/docs/sourcepro/stdlibref/2-10.html
  std::slice_array<double>  row(int i) { return x[std::slice((i-1)*cols,cols,1)]; }
  std::slice_array<double>  col(int i) { return x[std::slice(i-1,rows,cols)]; }
  const std::valarray<double>  row(int i) const {
    return std::valarray<double>(x[std::slice((i-1)*cols,cols,1)]); 
  }
  const std::valarray<double>  col(int i) const {
    return std::valarray<double>(x[std::slice(i-1,rows,cols)]);
  }

  std::slice_array<double>  diag() { return x[std::slice(0,rows,cols+1)]; }
  const std::valarray<double>  diag() const {
    return std::valarray<double>(x[std::slice(0,rows,cols+1)]);
  }

  ///fast access to single matrix elements
  double & elmt(int i, int j) { return x[(i-1) * cols + (j-1)]; }
  double elmt(int i, int j) const { return x[(i-1) * cols + (j-1)]; }

  ///Constructor straight from another valarray
  DoubleMatrix(const std::valarray<double> & valarr, int r, int c)
    : MatIndexable<double,DoubleMatrix>(), rows(r), cols(c), x(valarr) {}

  friend class ComplexMatrix;
  
public:

  /*
   *  CONSTRUCTORS
   */

  /// Constructor for matrix of zeroes 1..r,1..c
  DoubleMatrix(int r, int c) 
    : MatIndexable<double,DoubleMatrix>(), rows(r), cols(c), x(0.0, r*c) {}
  /// Makes diagonal square matrix: diagonal elements are from v, in order
  explicit DoubleMatrix(const DoubleVector &v);  
  ~DoubleMatrix() {}



  /// this is the only required operator to make this class work with ETs
  template<class E>
  DoubleMatrix & operator=(const MatXpr<double,E> & x) {
    *this = copy_from(x);
    return *this;
  }

  DoubleMatrix& operator=(const DoubleMatrix& other) {
    x.resize(other.x.size());
    x = other.x;
    rows = other.rows;
    cols = other.cols;
    return *this;
  }

  template <typename E>
  DoubleMatrix(const MatXpr<double, E> & m)
    : MatIndexable<double,DoubleMatrix>(), 
      rows(m.displayRows()), cols(m.displayCols()), x(0.0, rows*cols) 
  {
    assign_from(m);
  }

  /*
   *  ELEMENT ACCESS
   */
  double operator() (int i, int j) const; ///< to reference one element 
  double & operator() (int i, int j); ///< to reference one element 
  double display(int i, int j) const;///< ijth element
  double sum() const { return x.sum(); }
  double sumElements() const { return abs(x).sum(); }  
  /// Routines for outputting size of matrix
  int displayRows() const { return rows; };
  int displayCols() const { return cols; };
  const DoubleMatrix & display() const { return *this; }///< whole matrix returned

  DoubleVector displayRow(int i) const {
    DoubleVector temp(displayCols());
    int j; for (j=1; j<=displayCols(); j++) temp(j) = display(i, j);
    return temp;
  }

  DoubleVector displayCol(int i) const {
    DoubleVector temp(displayRows());
    int j; for (j=1; j<=displayRows(); j++) temp(j) = display(j, i);
    return temp;
  }


  void operator+=(DoubleMatrix & f) {
    x += f.x;
  }

  double nmin(int & k, int & l) const;

  /// Applies fn to every element of a matrix
  DoubleMatrix apply(double (*fn)(double)) const { 
    return DoubleMatrix(x.apply(fn),rows,cols); 
  }

  /// Sets diagonal entries equal to v, rest are 0
  const DoubleMatrix & operator=(double v);  

  double min(int & k, int & l) const; ///< minimum element
  double max(int & k, int & l) const; ///< maximum element
  /// Obvious elementary row/column operations
  void swaprows(int i, int j);
  void swapcols(int i,int j);
  std::size_t size() const { return x.size(); }
  /// change number of columns (Warning: can be slow because it internally copys a std::valarray<double>)
  void setCols(int);
  /// change number of rows (Warning: can be slow because it internally copys a std::valarray<double>)
  void setRows(int);
  /// resize matrix (Warning: can be slow because it internally copys a std::valarray<double>)
  void resize(int, int);

  double trace() const;///< trace must only be performed on a square matrix
  DoubleMatrix transpose() const; ///< can be any size
  bool testNan() const; ///< tests that all elements are not Nans 
  DoubleVector flatten() const; ///< matrix to vector, row by row

  /*
   *  NUMERICAL DIAGONALIZATION ROUTINES ETC.
   */

  /// Perform on a mixing matrix O such that
  /// diag v = O^double M O. It will put the diag evals in abs ascending order 
  /// and change O accordingly
  void associateOrderAbs(DoubleVector &v);
  // Perform on asymmetric matrices M such that
  /// diag = U^T double M V. It will put the diag evals in abs ascending order 
  /// and changes U,V accordingly
  void associateOrderAbs(DoubleMatrix & u, DoubleMatrix & v, DoubleVector &w) 
    const;
  /// fills in bottom-left hand corner of a matrix
  void symmetrise(); 
  /// Sums up the absolute value of the difference in each element of two
  /// matrices 
  double compare(const DoubleMatrix & a) const;
  double compare(const ComplexMatrix & a) const;
  /// Just returns a vector filled with the diagonal values 
  DoubleVector diagVals() const;
  DoubleMatrix inverse() const; ///< returns inverse of a matrix

  /// Uses singular SVD algorithm:
  /// \f$ A = U.W.V^T \f$ where W is a matrix of the eigenvalues,
  /// therefore \f$ W = U^T A V \f$.
  double diagonalise(DoubleMatrix & u, DoubleMatrix & v, DoubleVector & w) const;
  /// For SYMMETRIC MATRICES ONLY!
  /// \f$ A = V.W.V^T \f$ where W is a matrix of the eigenvalue therefore 
  /// \f$ W = V^T A V \f$. 
  double diagonaliseSym(DoubleMatrix & v, DoubleVector & w) const;
  double diagonaliseSym(ComplexMatrix & v, DoubleVector & w) const;

  /// Special case that's often used: eigenvalues of 2 by 2 symmetric matrix
  /// A(2,1) is assumed to be equal to A(1,2). theta is the mixing angle:
  /// m1 is NOT NECESSARILY > m2 
  // -[ cos theta    sin theta ]   A   [ cos theta -sin theta ] = diagonal
  // -[ -sin theta   cos theta ]       [ sin theta  cos theta ]
  DoubleVector sym2by2(double & theta) const; 
  /// 2 by 2 asymmetric matrices
  // -[ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ] = diagonal
  // -[ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
  DoubleVector asy2by2(double & thetaL, double & thetaR) const;
  /// converts a list of entries into a vector
  DoubleVector vectorfy() const;

  /// fill array, starting at offset
  void fillArray(double* array, unsigned offset = 0) const;
  bool operator==(const DoubleMatrix& o) const {
     if (displayRows() != o.displayRows())
        return false;
     if (displayCols() != o.displayCols())
        return false;
     for (int i = 1; i <= displayRows(); ++i)
        for (int k = 1; k <= displayCols(); ++k)
           if (display(i,k) != o.display(i,k))
              return false;
     return true;
  }
};

/*
 *  STANDARD INPUT / OUTPUT
 */

std::istream & operator>>(std::istream & left, DoubleMatrix &M);
std::ostream & operator <<(std::ostream &left, const DoubleMatrix &V);

/*
 *  NON-MEMBER DIAGONALIZATION ROUTINES ETC.
 */

/// Returns a 2x2 orthogonal matrix of rotation by angle theta
// [ cos theta    sin theta ]
// [-sin theta    cos theta ]
DoubleMatrix rot2d(double theta);
/// Returns 2 by 2 orthogonal mixing matrix
// [ -sin(theta)  cos(theta) ]
// [  cos(theta)  sin(theta) ] --
DoubleMatrix rot2dTwist(double theta);

/// LCT: Returns a 3x3 orthogonal matrix of rotation by angle theta.
/// Used in rotating CP-odd Higgs matrix
// [ -cos theta  sin theta  0 ]
// [ sin theta   cos theta  0 ]
// [ 0           0          1 ]
DoubleMatrix rot3d(double theta);

/// Redefines mixing matrices to be complex such that diagonal values are
/// positive for a 2 by 2: 
// [ cos thetaL    sin thetaL ]   A   [ cos thetaR -sin thetaR ]  = diag
// [ -sin thetaL   cos thetaL ]       [ sin thetaR  cos thetaR ]
// as given by asy2by2!
/// \f$ u^* A v^+ \f$ = mdiagpositive
void positivise(double thetaL, double thetaR, const DoubleVector & diag,
		  ComplexMatrix & u, ComplexMatrix & v);
void positivise(double theta, const DoubleVector& diag, ComplexMatrix& u);
void positivise(double theta, const DoubleVector& diag, DoubleMatrix& u);
/// Diagonalisation routines
void diagonaliseSvd(DoubleMatrix & a, DoubleVector & w, DoubleMatrix & v);
double pythagoras(double a, double b);
void diagonaliseJac(DoubleMatrix & a,  int n,  DoubleVector & d, DoubleMatrix
		    & v,  int *nrot);


/// fill array from valarray, starting at offset
void fillArray(const std::valarray<double>&, double*, unsigned offset = 0);

/*
 *  INLINE FUNCTION DEFINITIONS
 */

inline double DoubleMatrix::operator() (int i, int j) const { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i > rows || j > cols || i < 1 || j < 1) { 
    std::ostringstream ii;
    ii << "Trying to access " << i << "," << j << 
      "th element of DoubleMatrix\n" << *this;
    throw std::runtime_error(ii.str());
  }
#endif
  return elmt(i,j); 
}  

inline double &  DoubleMatrix::operator() (int i, int j) { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i > rows || j > cols || i < 1 || j < 1) { 
    std::ostringstream ii;
    ii << "Trying to access " << i << "," << j << 
      "th element of DoubleMatrix\n" << *this;
    throw std::runtime_error(ii.str());
  }
#endif
  return elmt(i,j); 
}  

inline double DoubleMatrix::display(int i, int j) const { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < 1 || i > displayRows() || j < 1 || j > displayCols()) { 
    std::ostringstream ii;
    ii << "Error: Requested display (" << i << "," << j << 
      ")th element of matrix " << *this; 
    throw std::runtime_error(ii.str());
  }
#endif
  return elmt(i,j); 
}  

/************************************
 * 
 *   *** COMPLEX VECTOR CLASS ***
 *
 ************************************/

/// Vector of double complex values
class ComplexVector : public Indexable<Complex, ComplexVector> {
private:
  int start, end; ///< Dimensions of vector
  std::valarray<Complex> x; ///< valarray of contents of vector
  ///Constructor directly from valarray, with given start and end
  ComplexVector(const std::valarray<Complex> & valarr, int s, int e)
    : Indexable<Complex, ComplexVector>(), start(s), end(e), x(valarr) {}  
  friend class ComplexMatrix;


public:

  /*
   *  CONSTRUCTORS
   */

  explicit ComplexVector(int e)///< Default dimensions start at 1, end at e
    : Indexable<Complex, ComplexVector>(), start(1), end(e), x(0.0,e) {}
  ComplexVector(int s, int e)///< Dimension starts at s and ends at e
    : Indexable<Complex, ComplexVector>(), start(s), end(e), x(0.0,e-s+1) {}
  ~ComplexVector() {}

  /// this is the only required operator to make this class work with ETs
  template<class E>
  ComplexVector & operator=(const Xpr<Complex,E> & x) {
    *this = copy_from(x);
    return *this;
  }

  ComplexVector& operator=(const ComplexVector& other) {
    x.resize(other.x.size());
    x = other.x;
    start = other.start;
    end = other.end;
    return *this;
  }

  template <typename E>
  ComplexVector(const Xpr<Complex,E> & v)
    : Indexable<Complex,ComplexVector>(), 
      start(v.displayStart()), end(v.displayEnd()), 
      x(0.0,v.displayEnd()-v.displayStart()+1)
  {
    assign_from(v);
  }


  /*
   *  ELEMENT ACCESS
   */

  Complex & operator() (int i);///< reference one element
  Complex operator() (int i) const { return x[i-start]; }
  Complex display(int i) const { return x[i-start]; }///<display ith element
  void set(int i, Complex f);///< set ith element to f

  int displayStart() const { return start; }///< displays start of dimension
  int displayEnd() const { return end; }///< displays end of dimension
  const ComplexVector & display() const { return *this; }///< displays whole vector

  /// Changes the length of a vector - copies as many elements of old one as
  /// possible, and fills any extra up with zeroes
  void setEnd(int e);

  /// Apply fn to every element
  ComplexVector apply(Complex (*fn)(Complex)) const { 
    return ComplexVector(x.apply(fn),start,end);
  }
  
  /// maximum absolute value, determined by definition of
  /// Complex::operator<(const Complex &) in mycomplex.h
  Complex max() const { return x.max(); }///< maximum absolute value
  Complex min(int & p) const;///< smallest absolute element
  
  void swap(int i, int j); ///< swap ith and jth elements

};

/*
 *  STANDARD INPUT / OUTPUT
 */

/// Formatted output
std::ostream & operator<<(std::ostream &left, const ComplexVector &V);
/// Formatted input
std::istream & operator>>(std::istream & left, ComplexVector &V);

/*
 *  INLINE FUNCTION DEFINITIONS
 */

inline Complex & ComplexVector::operator() (int i) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i>end || i<start) {
    std::ostringstream ii;
    ii << "Trying to access " << i << "th element of DoubleVector\n";
    ii << "start " << start << " end " << end;
    ii << *this;
    throw std::runtime_error(ii.str());
  }
#endif
  return x[i-start];
}

inline void ComplexVector::set(int i, Complex f) {
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < start || i > end) {
    std::ostringstream ii;
    ii << "Cannot access " << i << "th element of vector " << *this;
    throw std::runtime_error(ii.str());
  }
#endif
  x[i-start] = f;
}

/************************************
 * 
 *   *** COMPLEX MATRIX CLASS ***
 *
 ************************************/

/// matrix of complex double values dimensions (rows x cols)
class ComplexMatrix : public MatIndexable<Complex, ComplexMatrix> {
private:
  int rows, cols;///< gives size of matrix as 1..rows 1..cols
  std::valarray<Complex> x;///< data stored sequentially in this valarray, row after row
  
  /// for easy access to single rows and columns
  /// not very fast; for speedy access use for-loops over the elements
  // for usage see e.g. http://www.roguewave.com/support/docs/sourcepro/stdlibref/2-10.html
  std::slice_array<Complex>  row(int i) { return x[std::slice((i-1)*cols,cols,1)]; }
  std::slice_array<Complex>  col(int i) { return x[std::slice(i-1,rows,cols)]; }
  const std::valarray<Complex>  row(int i) const {
    return std::valarray<Complex>(x[std::slice((i-1)*cols,cols,1)]); 
  }
  const std::valarray<Complex>  col(int i) const {
    return std::valarray<Complex>(x[std::slice(i-1,rows,cols)]);
  }

  std::slice_array<Complex>  diag() { return x[std::slice(0,rows,cols+1)]; }
  const std::valarray<Complex>  diag() const {
    return std::valarray<Complex>(x[std::slice(0,rows,cols+1)]);
  }

  ///fast access to single matrix elements
  Complex & elmt(int i, int j) { return x[(i-1)*cols+(j-1)]; }
  const Complex elmt(int i, int j) const { return x[(i-1)*cols+(j-1)]; }

  ///Constructor straight from a valarray
  ComplexMatrix(const std::valarray<Complex> & valarr, int r, int c)
    : MatIndexable<Complex,ComplexMatrix>(), rows(r), cols(c), x(valarr) {}

  friend class DoubleMatrix;
  
public:

  /*
   *  CONSTRUCTORS
   */

  /// Default constructor: full of zeroes
  ComplexMatrix(int r, int c) 
    : MatIndexable<Complex,ComplexMatrix>(), rows(r), cols(c), x(0.0, r*c) {}
  /// Constructor sets matrix equal to m
  explicit ComplexMatrix(const DoubleMatrix & m);
  /// Constructor sets square matrix's diagonal values equal to elements of v
  explicit ComplexMatrix(const ComplexVector &v);  
  ~ComplexMatrix() {}

  /// this is the only required operator to make this class work with ETs
  template<class E>
  ComplexMatrix & operator=(const MatXpr<Complex,E> & x) {
    *this = copy_from(x);
    return *this;
  }

  ComplexMatrix& operator=(const ComplexMatrix& other) {
    x.resize(other.x.size());
    x = other.x;
    rows = other.rows;
    cols = other.cols;
    return *this;
  }

  template <typename E>
  ComplexMatrix(const MatXpr<Complex, E> & m)
    : MatIndexable<Complex,ComplexMatrix>(), 
      rows(m.displayRows()), cols(m.displayCols()), x(0.0, rows*cols) 
  {
    assign_from(m);
  }
  
  /*
   *  ELEMENT ACCESS
   */

  Complex & operator() (int i, int j); ///< Returns ijth element
  Complex operator() (int i, int j) const;
  Complex display(int i, int j) const;///< returns ijth element

  int displayRows() const { return rows; }; 
  int displayCols() const { return cols; };
  const ComplexMatrix & display() const { return *this; };///< returns whole matrix

  /// Sets diagonal entries equal to v, rest are 0
  const ComplexMatrix & operator=(const Complex &v);  
  const ComplexMatrix & operator+=(const DoubleMatrix&);
  ComplexMatrix operator+(const DoubleMatrix&);

  Complex min(int & k, int & l) const; ///< smallest absolute element
  /// Obvious elementary row/column operations
  void swaprows(int i, int j);///< Swaps row i with row j
  void swapcols(int i,int j);///< Swaps column i with column j

  /// change number of columns (Warning: can be slow because it internally copys a std::valarray<double>)
  void setCols(int);
  /// change number of rows (Warning: can be slow because it internally copys a std::valarray<double>)
  void setRows(int);
  /// resize matrix (Warning: can be slow because it internally copys a std::valarray<double>)
  void resize(int, int);

  Complex trace() const;  
  ComplexMatrix transpose() const;
  ComplexMatrix hermitianConjugate() const;
  ComplexMatrix complexConjugate() const;
  DoubleMatrix real() const;
  DoubleMatrix imag() const;
  
  /*
   *  NUMERICAL DIAGONALIZATION ROUTINES ETC.
   */

  /// Fills in lower bottom half of a square matrix copying the top right
  void symmetrise();
  /// Returns the sum of the modulus of the difference of each element
  double compare(const ComplexMatrix & a) const;
};

/*
 *  STANDARD INPUT / OUTPUT
 */

/// Formatted input
std::istream & operator>>(std::istream & left, ComplexMatrix &M);
/// Formatted output
std::ostream & operator<<(std::ostream &left, const ComplexMatrix &V);

/*
 *  INLINE FUNCTION DEFINITIONS
 */

inline Complex & ComplexMatrix::operator() (int i, int j) { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i > rows || j > cols || i < 0 || j < 0) { 
    std::ostringstream ii;
    ii << "Trying to access " << i << "," << j << 
      "th element of ComplexMatrix\n" << *this;
    throw std::runtime_error(ii.str());
  }
#endif
  return elmt(i,j); 
}  

inline Complex ComplexMatrix::operator() (int i, int j) const { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i > rows || j > cols || i < 0 || j < 0) { 
    std::ostringstream ii;
    ii << "Trying to access " << i << "," << j << 
      "th element of ComplexMatrix\n" << *this;
    throw std::runtime_error(ii.str());
  }
#endif
  return elmt(i,j); 
}  

inline Complex ComplexMatrix::display(int i, int j) const { 
#ifdef ARRAY_BOUNDS_CHECKING
  if (i < 1 || i > displayRows() || j < 1 || j > displayCols()) { 
    std::ostringstream ii;
    ii << "Error: Requested display (" << i << "," << j << 
	")th element of matrix " << *this; 
    throw std::runtime_error(ii.str());
    }
#endif
  return elmt(i,j); 
}  
}//namespace softsusy
#endif
