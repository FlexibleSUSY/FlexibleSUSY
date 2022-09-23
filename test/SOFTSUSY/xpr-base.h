/** \file xpr-base.h
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach, this part by D Grellscheid
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://projects.hepforge.org/softsusy/

    \brief Symbolic under-object for vectors/matrices for speed upgrade. The
    linear algebra operations are often set up at compile time rather than at
    run-time
*/



#ifndef XPR_BASE_H
#define XPR_BASE_H

#include <exception>
#include <iostream>
#include <math.h>

/****************** VECTOR **************************/

/// Under-object for speed upgrade of linear algebra
template<typename T, typename E>
class Xpr //  : public Indexable<T,Xpr> // inheriting doesn't work. Pity.
{
  E e;
public:
  Xpr( const E& e_ ) : e(e_) {}
  T operator() ( int n ) const { return e(n); }
  int displayStart() const { return e.displayStart(); }
  int displayEnd() const { return e.displayEnd(); }
};

/// const version of under-object for speed upgrade of linear algebra
template<typename T, typename Container>
class ConstRef {
  const Container& c;
public:
  ConstRef( const Container& c_) : c(c_) {}
  T operator() ( int n ) const { return c(n); }
  int displayStart() const { return c.displayStart(); }
  int displayEnd() const { return c.displayEnd(); }
};

/// binary operator on two xpr objects
template<typename T, typename A, typename B, typename Op>
class XprBinOp {
  A a;
  B b;
public:
  XprBinOp( const A& a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int n ) const
  {
    return Op::apply( a(n), b(n) );
  }
  int displayStart() const { return a.displayStart(); }
  int displayEnd() const { return a.displayEnd(); }
};

/// scalar operator between two xpr object
template<typename T, typename B, typename Op>
class XprScalOp {
  T a;
  B b;
public:
  XprScalOp( T a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int n ) const
  {
    return Op::apply( a, b(n) );
  }
  int displayStart() const { return b.displayStart(); }
  int displayEnd() const { return b.displayEnd(); }
};

/// provides indexing of xpr type objects
template<typename T, typename C>
class Indexable {
public:
  explicit Indexable() {}
  T operator() ( int n ) const
    {
      return static_cast<const C*>(this)->operator()(n);
    }

  int displayStart() const
  {
    return static_cast<const C*>(this)->displayStart();
  }
  
  int displayEnd() const
  {
    return static_cast<const C*>(this)->displayEnd();
  }

  template<class E>
  C& assign_from(const Xpr<T,E>& x )
  {
    C *me = static_cast<C*>(this);
    for(int i=me->displayStart(); i <= me->displayEnd(); ++i)
      me->operator()(i) = x(i);
    return *me;
  }

  template<class E>
  C copy_from(const Xpr<T,E>& x )
  {
    C *meptr = static_cast<C*>(this);
    C me(meptr->displayStart(), meptr->displayEnd());
    for(int i=me.displayStart(); i <= me.displayEnd(); ++i)
      me.operator()(i) = x(i);
    return me;
  }
 
  template<class E>
  Indexable& operator+=( const Xpr<T,E>& x )
  {
    C *me = static_cast<C*>(this);
    for ( int i=me->displayStart(); i <= me->displayEnd(); ++i)
      me->operator() (i) += x(i);
    return *this;
  }
  
  template<class E>
  Indexable& operator-=( const Xpr<T,E>& x )
  {
    C *me = static_cast<C*>(this);
    for ( int i=me->displayStart(); i <= me->displayEnd(); ++i)
      me->operator() (i) -= x(i);
    return *this;
  }
  
  Indexable& operator*=( T x )
  {
    C *me = static_cast<C*>(this);
    for ( int i=me->displayStart(); i <= me->displayEnd(); ++i)
      me->operator() (i) *= x;
    return *this;
  }

  // ...
};

/***************** operations ********************/

/// Add two xpr objects
template<typename T>
class OpAdd {
public:
  static inline T apply( const T& a, const T& b )
  {
    return a + b;
  }
};

/// subtract two xpr objects
template<typename T>
class OpSubtract {
public:
  static inline T apply( const T& a, const T& b )
  {
    return a - b;
  }
};

/// multiply two xpr objects
template<typename T>
class OpMultiply {
public:
  static inline T apply( const T& a, const T& b ) {
      return a * b;
  }
};

/// take the negative of an xpr object
template<typename T>
class OpNegate {
public:
  static inline T apply( const T& a )
  {
    return -a ;
  }
};


/************** MATRIX ********************/
/// matrix form of xpr object
template<typename T, typename E>
class MatXpr
{
  E e;
public:
  MatXpr( const E& e_ ) : e(e_) {}
  T operator() ( int m, int n ) const { return e(m,n); }
  int displayRows() const { return e.displayRows(); }
  int displayCols() const { return e.displayCols(); }
  T trace() const { 
    T retval = T();
    for (int i=1; i<= displayRows(); ++i)
      retval += this->operator()(i,i);
    return retval;
  }
};


/// const version of xpr matrix
template<typename T, typename Container>
class MatConstRef {
  const Container& c;
public:
  MatConstRef( const Container& c_) : c(c_) {}
  T operator() ( int m, int n ) const { return c(m,n); }
  int displayRows() const { return c.displayRows(); }
  int displayCols() const { return c.displayCols(); }
};

/// binary operator between two matrices
template<typename T, typename A, typename B, typename Op>
class MatXprBinOp {
  A a;
  B b;
public:
  MatXprBinOp( const A& a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int m, int n ) const
  {
    return Op::apply( a(m,n), b(m,n) );
  }
  int displayRows() const { return a.displayRows(); }
  int displayCols() const { return a.displayCols(); }
};

/// outer form of xpr matrix operator
template<typename T, typename A, typename B, typename Op>
class MatXprOuterOp {
  A a;
  B b;
public:
  MatXprOuterOp( const A& a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int m, int n ) const
  {
    return Op::apply( a(m), b(n) );
  }
  int displayRows() const { return a.displayEnd(); }
  int displayCols() const { return b.displayEnd(); }
};

/// xpr Matrix multiplication
template<typename T, typename A, typename B>
class MatXprMatMultOp {
  A a;
  B b;
public:
  MatXprMatMultOp( const A& a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int m, int n ) const
  {
    T value = T();
    for (int i=1; i<= a.displayCols(); ++i) {
#ifdef ARRAY_BOUNDS_CHECKING
      int n2, n3;
      frexp(a(m,i), &n2); frexp(b(i,n), &n3);
      if (n2+n3>=1024 || n2+n3 <= -1024) { 
	throw std::runtime_error("Multiply * overload; overflow\n");
	if (n2+n3<=-1024) return 1e-290;
	else return 1e290;
      } else 
#endif
	value += a(m,i) * b(i,n);
    }
    return value;
  }
  int displayRows() const { return a.displayRows(); }
  int displayCols() const { return b.displayCols(); }
};

/// another xpr matrix multiplication
template<typename T, typename A, typename B>
class XprMatMultOp {
  A a;
  B b;
public:
  XprMatMultOp( const A& a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int n ) const
  {
    T value = T();
    for (int i=1; i<= a.displayCols(); ++i) {
#ifdef ARRAY_BOUNDS_CHECKING
      int n2, n3;
      frexp(a(n,i), &n2); frexp(b(i), &n3);
      if (n2+n3>=1024 || n2+n3 <= -1024) { 
	throw std::runtime_error("Multiply * overload; overflow\n");
	if (n2+n3<=-1024) return 1e-290;
	else return 1e290;
      } else 
#endif
      value += a(n,i) * b(i);
    }
    return value;
  }
  int displayStart() const { return 1; }
  int displayEnd() const { return a.displayRows(); }
};

/// unary operator on an xpr matrix
template<typename T, typename A, typename Op>
class MatXprUnaryOp {
  A a;
public:
  MatXprUnaryOp( const A& a_)
    : a(a_) {}
  T operator() ( int m, int n ) const
  {
    return Op::apply( a(m,n) );
  }
  int displayRows() const { return a.displayRows(); }
  int displayCols() const { return a.displayCols(); }
};

/// scalar operator on xpr matrix
template<typename T, typename B, typename Op>
class MatXprScalOp {
  T a;
  B b;
public:
  MatXprScalOp( T a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int m, int n ) const
  {
    return Op::apply( a, b(m,n) );
  }
  int displayRows() const { return b.displayRows(); }
  int displayCols() const { return b.displayCols(); }
};

/// scalar operator on diagonal parts of xpr matrix 
template<typename T, typename B, typename Op>
class MatXprScalDiagOp1 {
  T a;
  B b;
public:
  MatXprScalDiagOp1( T a_, const B& b_ )
    : a(a_), b(b_) {}
  T operator() ( int m, int n ) const
  {
    return (m == n) ? Op::apply( a, b(m,n) ) : b(m,n);
  }
  int displayRows() const { return b.displayRows(); }
  int displayCols() const { return b.displayCols(); }
};

/// scalar operator on diagonal parts of xpr matrix - different place for const
template<typename T, typename B, typename Op>
class MatXprScalDiagOp2 {
  T a;
  B b;
public:
  MatXprScalDiagOp2(const B& b_, T a_ )
    : a(a_), b(b_) {}
  T operator() ( int m, int n ) const
  {
    return (m == n) ? Op::apply( b(m,n) , a ) : b(m,n);
  }
  int displayRows() const { return b.displayRows(); }
  int displayCols() const { return b.displayCols(); }
};

/// Indexable form of xpr matrix
template<typename T, typename C>
class MatIndexable {
public:
  explicit MatIndexable() {}
  T operator() ( int m, int n ) const
    {
      return static_cast<const C*>(this)->operator()(m,n);
    }

  int displayRows() const
  {
    return static_cast<const C*>(this)->displayRows();
  }
  
  
  int displayCols() const
  {
    return static_cast<const C*>(this)->displayCols();
  }

  template<class E>
  C& assign_from(const MatXpr<T,E>& x )
  {
    C *me = static_cast<C*>(this);
    for(int m=1; m <= me->displayRows(); ++m)
      for (int n=1; n <= me->displayCols(); ++n)
	me->operator()(m,n) = x(m,n);
    return *me;
  }

  template<class E>
  C copy_from(const MatXpr<T,E>& x )
  {
    C * meptr = static_cast<C*>(this);
    C me(meptr->displayRows(), meptr->displayCols());
    for(int m=1; m <= me.displayRows(); ++m)
      for (int n=1; n <= me.displayCols(); ++n)
	me.operator()(m,n) = x(m,n);
    return me;
  }
  
  template<class E>
  MatIndexable& operator+=( const MatXpr<T,E>& x )
  {
    C *me = static_cast<C*>(this);
    for(int m=1; m <= me->displayRows(); ++m)
      for (int n=1; n <= me->displayCols(); ++n)
	me->operator()(m,n) += x(m,n);
    return *this;
  }
  

  MatIndexable& operator*=(T x)
  {
    C *me = static_cast<C*>(this);
    for(int m=1; m <= me->displayRows(); ++m)
      for (int n=1; n <= me->displayCols(); ++n)
	me->operator()(m,n) *= x;
    return *this;
  }

};




#endif
