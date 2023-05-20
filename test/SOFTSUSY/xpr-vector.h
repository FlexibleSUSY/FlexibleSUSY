/** \file xpr-vector.h
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach, this part by D Grellscheid
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://projects.hepforge.org/softsusy/

    \brief Symbolic object for vectors for speed upgrade
*/

#ifndef XPR_VECTOR_H
#define XPR_VECTOR_H

#include "xpr-base.h"

/****** DOT PRODUCT *** temporary *****/
/// temporary xpr A.B
template<typename T, typename A, typename B >
T dot(const Indexable<T,A>& a, const Indexable<T,B>& b )
{
  T retval = T();
  for (int i = a.displayStart(); i <= a.displayEnd(); ++i)
    retval += a(i) * b(i);
  return retval;
}

/// temporary xpr A.B
template<typename T, typename A, typename B >
T dot(const Indexable<T,A>& a, const Xpr<T,B>& b )
{
  T retval = T();
  for (int i = a.displayStart(); i <= a.displayEnd(); ++i)
    retval += a(i) * b(i);
  return retval;
}

/// temporary xpr A.B
template<typename T, typename A, typename B >
T dot(const Xpr<T,A>& a, const Indexable<T,B>& b )
{
  T retval = T();
  for (int i = b.displayStart(); i <= b.displayEnd(); ++i)
    retval += a(i) * b(i);
  return retval;
}

/// temporary xpr A.B
template<typename T, typename A, typename B >
T dot(const Xpr<T,A>& a, const Xpr<T,B>& b )
{
  T retval = T();
  for (int i = b.displayStart(); i <= b.displayEnd(); ++i)
    retval += a(i) * b(i);
  return retval;
}

/******************* ADDITION ******************************/

/// temporary xpr A+B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpAdd<T> > >
operator+( const Indexable<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpAdd<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// temporary xpr A+B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpAdd<T> > >
operator+( const Indexable<T,A>& a, const Xpr<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpAdd<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/// temporary xpr A+B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpAdd<T> > >
operator+( const Xpr<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpAdd<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// temporary xpr A+B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpAdd<T> > >
operator+( const Xpr<T,A>& a, const Xpr<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpAdd<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/**************************** SUBTRACTION ***************************/

/// temporary xpr A-B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpSubtract<T> > >
operator-( const Indexable<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// temporary xpr A-B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpSubtract<T> > >
operator-( const Indexable<T,A>& a, const Xpr<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/// temporary xpr A-B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpSubtract<T> > >
operator-( const Xpr<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// temporary xpr A-B
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpSubtract<T> > >
operator-( const Xpr<T,A>& a, const Xpr<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/*********** ELEMENT WISE MULTIPLY ***************************/
/// xpr A*x, for each element
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpMultiply<T> > >
operator*( const Indexable<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// xpr A*x, for each element
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpMultiply<T> > >
operator*( const Indexable<T,A>& a, const Xpr<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/// xpr A*x, for each element
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpMultiply<T> > >
operator*( const Xpr<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// xpr A*x, for each element
template<typename T, typename A, typename B>
Xpr< T, 
     XprBinOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpMultiply<T> > >
operator*( const Xpr<T,A>& a, const Xpr<T,B>& b )
{
  typedef XprBinOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/********* SCALAR MULTIPLY ***************/
/// xpr A*x, x scalar
template<typename T, typename B>
Xpr< T, 
     XprScalOp<T, 
               ConstRef<T,Indexable<T,B> >, 
	      OpMultiply<T> > >
operator*( T a, const Indexable<T,B>& b )
{
  typedef XprScalOp< T, 
    ConstRef<T,Indexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( a,
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// xpr A*x, x scalar
template<typename T, typename B>
Xpr< T, 
     XprScalOp<T, 
               ConstRef<T,Indexable<T,B> >, 
	      OpMultiply<T> > >
operator*( const Indexable<T,B>& b, T a )
{
  typedef XprScalOp< T, 
    ConstRef<T,Indexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( a,
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// xpr A*x, x scalar
template<typename T, typename B>
Xpr< T, 
     XprScalOp<T, 
               ConstRef<T,Xpr<T,B> >, 
	      OpMultiply<T> > >
operator*( T a, const Xpr<T,B>& b )
{
  typedef XprScalOp< T, 
    ConstRef<T,Xpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( a,
				ConstRef<T,Xpr<T,B> >(b) ));
}

/// xpr A*x, x scalar
template<typename T, typename B>
Xpr< T, 
     XprScalOp<T, 
               ConstRef<T,Xpr<T,B> >, 
	      OpMultiply<T> > >
operator*( const Xpr<T,B>& b, T a )
{
  typedef XprScalOp< T, 
    ConstRef<T,Xpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( a,
				ConstRef<T,Xpr<T,B> >(b) ));
}

#endif
