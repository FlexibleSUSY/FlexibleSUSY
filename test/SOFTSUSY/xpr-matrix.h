/** \file xpr-matrix.h
    - Project:     SOFTSUSY 
    - Author:      Ben Allanach, this part by D Grellscheid
    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    - Webpage:     http://projects.hepforge.org/softsusy/

    \brief Symbolic object for matrices for speed upgrade
*/

#ifndef XPR_MATRIX_H
#define XPR_MATRIX_H

#include "xpr-base.h"
#include "mycomplex.h"

#define Complex softsusy::Complex

/******************* ADDITION ******************************/
/// add two xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      MatConstRef<T,MatIndexable<T,B> >, 
	      OpAdd<T> > >
operator+( const MatIndexable<T,A>& a, const MatIndexable<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatIndexable<T,A> >,
    MatConstRef<T,MatIndexable<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// add two xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      MatConstRef<T,MatIndexable<T,B> >, 
	      OpAdd<T> > >
operator+( const MatXpr<T,A>& a, const MatIndexable<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatXpr<T,A> >,
    MatConstRef<T,MatIndexable<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// add two xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      MatConstRef<T,MatXpr<T,B> >, 
	      OpAdd<T> > >
operator+( const MatIndexable<T,A>& a, const MatXpr<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatIndexable<T,A> >,
    MatConstRef<T,MatXpr<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// add two xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      MatConstRef<T,MatXpr<T,B> >, 
	      OpAdd<T> > >
operator+( const MatXpr<T,A>& a, const MatXpr<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatXpr<T,A> >,
    MatConstRef<T,MatXpr<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/****** ADDITION TO TRACE *******************/
/// scalar + matrix
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp1<T, 
               MatConstRef<T,MatIndexable<T,B> >, 
	      OpAdd<T> > >
operator+( T a, const MatIndexable<T,B>& b )
{
  typedef MatXprScalDiagOp1< T, 
    MatConstRef<T,MatIndexable<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// scalar + matrix
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp2<T, 
               MatConstRef<T,MatIndexable<T,B> >, 
	      OpAdd<T> > >
operator+( const MatIndexable<T,B>& b, T a )
{
  typedef MatXprScalDiagOp2< T, 
    MatConstRef<T,MatIndexable<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( 
				MatConstRef<T,MatIndexable<T,B> >(b) , a));
}

/// scalar + matrix
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp1<T, 
               MatConstRef<T,MatXpr<T,B> >, 
	      OpAdd<T> > >
operator+( T a, const MatXpr<T,B>& b )
{
  typedef MatXprScalDiagOp1< T, 
    MatConstRef<T,MatXpr<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// scalar + matrix
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp2<T, 
               MatConstRef<T,MatXpr<T,B> >, 
	      OpAdd<T> > >
operator+( const MatXpr<T,B>& b, T a )
{
  typedef MatXprScalDiagOp2< T, 
    MatConstRef<T,MatXpr<T,B> >, 
    OpAdd<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( 
				MatConstRef<T,MatXpr<T,B> >(b) , a));
}



/******************* SUBTRACTION ******************************/
/// xpr matrix subtraction
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      MatConstRef<T,MatIndexable<T,B> >, 
	      OpSubtract<T> > >
operator-( const MatIndexable<T,A>& a, const MatIndexable<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatIndexable<T,A> >,
    MatConstRef<T,MatIndexable<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// xpr matrix subtraction
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      MatConstRef<T,MatIndexable<T,B> >, 
	      OpSubtract<T> > >
operator-( const MatXpr<T,A>& a, const MatIndexable<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatXpr<T,A> >,
    MatConstRef<T,MatIndexable<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// xpr matrix subtraction
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      MatConstRef<T,MatXpr<T,B> >, 
	      OpSubtract<T> > >
operator-( const MatIndexable<T,A>& a, const MatXpr<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatIndexable<T,A> >,
    MatConstRef<T,MatXpr<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// xpr matrix subtraction
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprBinOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      MatConstRef<T,MatXpr<T,B> >, 
	      OpSubtract<T> > >
operator-( const MatXpr<T,A>& a, const MatXpr<T,B>& b )
{
  typedef MatXprBinOp< T, MatConstRef<T,MatXpr<T,A> >,
    MatConstRef<T,MatXpr<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/************* OUTER PRODUCT ************/
/// outer product of xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprOuterOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpMultiply<T> > >
outerProduct( const Xpr<T,A>& a, const Xpr<T,B>& b )
{
  typedef MatXprOuterOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/// outer product of xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprOuterOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Xpr<T,B> >, 
	      OpMultiply<T> > >
outerProduct( const Indexable<T,A>& a, const Xpr<T,B>& b )
{
  typedef MatXprOuterOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Xpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Xpr<T,B> >(b) ));
}

/// outer product of xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprOuterOp<T, 
	      ConstRef<T,Xpr<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpMultiply<T> > >
outerProduct( const Xpr<T,A>& a, const Indexable<T,B>& b )
{
  typedef MatXprOuterOp< T, ConstRef<T,Xpr<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( ConstRef<T,Xpr<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// outer product of xpr matrices
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprOuterOp<T, 
	      ConstRef<T,Indexable<T,A> >,
	      ConstRef<T,Indexable<T,B> >, 
	      OpMultiply<T> > >
outerProduct( const Indexable<T,A>& a, const Indexable<T,B>& b )
{
  typedef MatXprOuterOp< T, ConstRef<T,Indexable<T,A> >,
    ConstRef<T,Indexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( ConstRef<T,Indexable<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/************* UNARY NEGATIVE **************/
/// negative of xpr matrices
template<typename T, typename A>
MatXpr< T, 
     MatXprUnaryOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      OpNegate<T> > >
operator-( const MatXpr<T,A>& a )
{
  typedef MatXprUnaryOp< T, MatConstRef<T,MatXpr<T,A> >,
    OpNegate<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a) ));
}

/// negative of xpr matrices
template<typename T, typename A>
MatXpr< T, 
     MatXprUnaryOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      OpNegate<T> > >
operator-( const MatIndexable<T,A>& a )
{
  typedef MatXprUnaryOp< T, MatConstRef<T,MatIndexable<T,A> >,
    OpNegate<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a) ));
}

/****** SUBTRACTION FROM TRACE *******************/
/// xpr matrices - some number
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp1<T, 
               MatConstRef<T,MatIndexable<T,B> >, 
	      OpSubtract<T> > >
operator-( T a, const MatIndexable<T,B>& b )
{
  typedef MatXprScalDiagOp1< T, 
    MatConstRef<T,MatIndexable<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// xpr matrices - some number
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp2<T, 
               MatConstRef<T,MatIndexable<T,B> >, 
	      OpSubtract<T> > >
operator-( const MatIndexable<T,B>& b, T a )
{
  typedef MatXprScalDiagOp2< T, 
    MatConstRef<T,MatIndexable<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,B> >(b) , a));
}

/// xpr matrices - some number
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp1<T, 
               MatConstRef<T,MatXpr<T,B> >, 
	      OpSubtract<T> > >
operator-( T a, const MatXpr<T,B>& b )
{
  typedef MatXprScalDiagOp1< T, 
    MatConstRef<T,MatXpr<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// xpr matrices - some number
template<typename T, typename B>
MatXpr< T, 
     MatXprScalDiagOp2<T, 
               MatConstRef<T,MatXpr<T,B> >, 
	      OpSubtract<T> > >
operator-( const MatXpr<T,B>& b, T a )
{
  typedef MatXprScalDiagOp2< T, 
    MatConstRef<T,MatXpr<T,B> >, 
    OpSubtract<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT(	MatConstRef<T,MatXpr<T,B> >(b) , a ));
}

/// multiply xpr matrix by a scalar
template<typename T, typename B>
MatXpr< T, 
     MatXprScalOp<T, 
               MatConstRef<T,MatIndexable<T,B> >, 
	      OpMultiply<T> > >
operator*( T a, const MatIndexable<T,B>& b )
{
  typedef MatXprScalOp< T, 
    MatConstRef<T,MatIndexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >( ExprT( a,	MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// multiply xpr matrix by a scalar
template<typename T, typename B>
MatXpr< T,
     MatXprScalOp<T,
               MatConstRef<T,MatIndexable<T,B> >,
	      OpMultiply<T> > >
operator*( int a, const MatIndexable<T,B>& b )
{
  typedef MatXprScalOp< T,
    MatConstRef<T,MatIndexable<T,B> >,
    OpMultiply<T> > ExprT;

  return MatXpr< T, ExprT >( ExprT( static_cast<T>(a), MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// multiply xpr complex matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatIndexable<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( double a, const MatIndexable<Complex,B>& b )
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatIndexable<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatIndexable<Complex,B> >(b) ));
}

/// multiply xpr complex matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatIndexable<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( const std::complex<double>& a, const MatIndexable<Complex,B>& b )
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatIndexable<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatIndexable<Complex,B> >(b) ));
}

/// multiply xpr matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatXpr<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( double a, const MatXpr<Complex,B>& b )
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatXpr<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatXpr<Complex,B> >(b) ));
}

/// multiply xpr matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatXpr<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( const std::complex<double>& a, const MatXpr<Complex,B>& b )
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatXpr<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatXpr<Complex,B> >(b) ));
}

/// multiply xpr matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatIndexable<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( const MatIndexable<Complex,B>& b , double a)
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatIndexable<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatIndexable<Complex,B> >(b) ));
}

/// multiply xpr matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatIndexable<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( const MatIndexable<Complex,B>& b , const std::complex<double>& a)
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatIndexable<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatIndexable<Complex,B> >(b) ));
}

/// multiply xpr complex matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatXpr<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( const MatXpr<Complex,B>& b , double a)
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatXpr<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatXpr<Complex,B> >(b) ));
}

/// multiply xpr complex matrix by a scalar
template<typename B>
MatXpr< Complex, 
     MatXprScalOp<Complex, 
               MatConstRef<Complex,MatXpr<Complex,B> >, 
	      OpMultiply<Complex> > >
operator*( const MatXpr<Complex,B>& b , const std::complex<double>& a)
{
  typedef MatXprScalOp< Complex, 
    MatConstRef<Complex,MatXpr<Complex,B> >, 
    OpMultiply<Complex> > ExprT;
  
  return MatXpr< Complex, ExprT >( ExprT( a,	MatConstRef<Complex,MatXpr<Complex,B> >(b) ));
}


/// multiply xpr matrix by scalar
template<typename T, typename B>
MatXpr< T, 
     MatXprScalOp<T, 
               MatConstRef<T,MatIndexable<T,B> >, 
	      OpMultiply<T> > >
operator*( const MatIndexable<T,B>& b, T a )
{
  typedef MatXprScalOp< T, 
    MatConstRef<T,MatIndexable<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// multiply xpr matrix by scalar
template<typename T, typename B>
MatXpr< T,
     MatXprScalOp<T,
               MatConstRef<T,MatIndexable<T,B> >,
	      OpMultiply<T> > >
operator*( const MatIndexable<T,B>& b, int a )
{
  typedef MatXprScalOp< T,
    MatConstRef<T,MatIndexable<T,B> >,
    OpMultiply<T> > ExprT;

  return MatXpr< T, ExprT >(
                             ExprT( static_cast<T>(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// multiply xpr matrix by scalar
template<typename T, typename B>
MatXpr< T, 
     MatXprScalOp<T, 
               MatConstRef<T,MatXpr<T,B> >, 
	      OpMultiply<T> > >
operator*( T a, const MatXpr<T,B>& b )
{
  typedef MatXprScalOp< T, 
    MatConstRef<T,MatXpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// multiply xpr matrix by scalar
template<typename T, typename B>
MatXpr< T,
     MatXprScalOp<T,
               MatConstRef<T,MatXpr<T,B> >,
	      OpMultiply<T> > >
operator*( int a, const MatXpr<T,B>& b )
{
  typedef MatXprScalOp< T,
    MatConstRef<T,MatXpr<T,B> >,
    OpMultiply<T> > ExprT;

  return MatXpr< T, ExprT >(
                             ExprT( static_cast<T>(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// multiply xpr matrix by scalar
template<typename T, typename B>
MatXpr< T, 
     MatXprScalOp<T, 
               MatConstRef<T,MatXpr<T,B> >, 
	      OpMultiply<T> > >
operator*( const MatXpr<T,B>& b, T a )
{
  typedef MatXprScalOp< T, 
    MatConstRef<T,MatXpr<T,B> >, 
    OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( a,
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// multiply xpr matrix by scalar
template<typename T, typename B>
MatXpr< T,
     MatXprScalOp<T,
               MatConstRef<T,MatXpr<T,B> >,
	      OpMultiply<T> > >
operator*( const MatXpr<T,B>& b, int a )
{
  typedef MatXprScalOp< T,
    MatConstRef<T,MatXpr<T,B> >,
    OpMultiply<T> > ExprT;

  return MatXpr< T, ExprT >(
                             ExprT( static_cast<T>(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/****** DIVIDE BY SCALAR ***************/
/// divide xpr matrix by scalar
template<typename T, typename B >
MatXpr< T, MatXprScalOp<T, MatConstRef<T,MatIndexable<T,B> >, OpMultiply<T> > >
operator/( const MatIndexable<T,B>& b, T a )
{
  typedef MatXprScalOp< T, MatConstRef<T,MatIndexable<T,B> >, OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(ExprT (1.0/a, MatConstRef<T,MatIndexable<T,B> >(b))  );
}


/// divide xpr matrix by scalar
template<typename T, typename B >
MatXpr< T, MatXprScalOp<T, MatConstRef<T,MatXpr<T,B> >, OpMultiply<T> > >
operator/( const MatXpr<T,B>& b, T a )
{
  typedef MatXprScalOp< T, MatConstRef<T,MatXpr<T,B> >, OpMultiply<T> > ExprT;
  
  return MatXpr< T, ExprT >(ExprT (1.0/a, MatConstRef<T,MatXpr<T,B> >(b))  );
}




/************ MATRIX MULTIPLY ***************/
/// xpr matrix multiplication
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprMatMultOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      MatConstRef<T,MatIndexable<T,B> > > >
operator*( const MatIndexable<T,A>& a, const MatIndexable<T,B>& b )
{
  typedef MatXprMatMultOp< T, MatConstRef<T,MatIndexable<T,A> >,
    MatConstRef<T,MatIndexable<T,B> > > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// xpr matrix multiplication
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprMatMultOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      MatConstRef<T,MatIndexable<T,B> > > >
operator*( const MatXpr<T,A>& a, const MatIndexable<T,B>& b )
{
  typedef MatXprMatMultOp< T, MatConstRef<T,MatXpr<T,A> >,
    MatConstRef<T,MatIndexable<T,B> > > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a),
				MatConstRef<T,MatIndexable<T,B> >(b) ));
}

/// xpr matrix multiplication
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprMatMultOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      MatConstRef<T,MatXpr<T,B> > > >
operator*( const MatIndexable<T,A>& a, const MatXpr<T,B>& b )
{
  typedef MatXprMatMultOp< T, MatConstRef<T,MatIndexable<T,A> >,
    MatConstRef<T,MatXpr<T,B> > > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// xpr matrix multiplication
template<typename T, typename A, typename B>
MatXpr< T, 
     MatXprMatMultOp<T, 
	      MatConstRef<T,MatXpr<T,A> >,
	      MatConstRef<T,MatXpr<T,B> > > >
operator*( const MatXpr<T,A>& a, const MatXpr<T,B>& b )
{
  typedef MatXprMatMultOp< T, MatConstRef<T,MatXpr<T,A> >,
    MatConstRef<T,MatXpr<T,B> > > ExprT;
  
  return MatXpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatXpr<T,A> >(a),
				MatConstRef<T,MatXpr<T,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<Complex,MatIndexable<Complex,A> >,
	      MatConstRef<double,MatIndexable<double,B> > > >
operator*( const MatIndexable<Complex,A>& a, const MatIndexable<double,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<Complex,MatIndexable<Complex,A> >,
    MatConstRef<double,MatIndexable<double,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<Complex,MatIndexable<Complex,A> >(a),
				MatConstRef<double,MatIndexable<double,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<Complex,MatXpr<Complex,A> >,
	      MatConstRef<double,MatIndexable<double,B> > > >
operator*( const MatXpr<Complex,A>& a, const MatIndexable<double,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<Complex,MatXpr<Complex,A> >,
    MatConstRef<double,MatIndexable<double,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<Complex,MatXpr<Complex,A> >(a),
				MatConstRef<double,MatIndexable<double,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<Complex,MatIndexable<Complex,A> >,
	      MatConstRef<double,MatXpr<double,B> > > >
operator*( const MatIndexable<Complex,A>& a, const MatXpr<double,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<Complex,MatIndexable<Complex,A> >,
    MatConstRef<double,MatXpr<double,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<Complex,MatIndexable<Complex,A> >(a),
				MatConstRef<double,MatXpr<double,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<Complex,MatXpr<Complex,A> >,
	      MatConstRef<double,MatXpr<double,B> > > >
operator*( const MatXpr<Complex,A>& a, const MatXpr<double,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<Complex,MatXpr<Complex,A> >,
    MatConstRef<double,MatXpr<double,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<Complex,MatXpr<Complex,A> >(a),
				MatConstRef<double,MatXpr<double,B> >(b) ));
}

// special 2
/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<double,MatIndexable<double,A> >,
	      MatConstRef<Complex,MatIndexable<Complex,B> > > >
operator*( const MatIndexable<double,A>& a, const MatIndexable<Complex,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<double,MatIndexable<double,A> >,
    MatConstRef<Complex,MatIndexable<Complex,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<double,MatIndexable<double,A> >(a),
				MatConstRef<Complex,MatIndexable<Complex,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<double,MatXpr<double,A> >,
	      MatConstRef<Complex,MatIndexable<Complex,B> > > >
operator*( const MatXpr<double,A>& a, const MatIndexable<Complex,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<double,MatXpr<double,A> >,
    MatConstRef<Complex,MatIndexable<Complex,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<double,MatXpr<double,A> >(a),
				MatConstRef<Complex,MatIndexable<Complex,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<double,MatIndexable<double,A> >,
	      MatConstRef<Complex,MatXpr<Complex,B> > > >
operator*( const MatIndexable<double,A>& a, const MatXpr<Complex,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<double,MatIndexable<double,A> >,
    MatConstRef<Complex,MatXpr<Complex,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<double,MatIndexable<double,A> >(a),
				MatConstRef<Complex,MatXpr<Complex,B> >(b) ));
}

/// xpr matrix multiplication specialisation
template<typename A, typename B>
MatXpr< Complex, 
     MatXprMatMultOp<Complex, 
	      MatConstRef<double,MatXpr<double,A> >,
	      MatConstRef<Complex,MatXpr<Complex,B> > > >
operator*( const MatXpr<double,A>& a, const MatXpr<Complex,B>& b )
{
  typedef MatXprMatMultOp< Complex, MatConstRef<double,MatXpr<double,A> >,
    MatConstRef<Complex,MatXpr<Complex,B> > > ExprT;
  
  return MatXpr< Complex, ExprT >(
			 ExprT( MatConstRef<double,MatXpr<double,A> >(a),
				MatConstRef<Complex,MatXpr<Complex,B> >(b) ));
}

/************ MATRIX * VECTOR MULTIPLY ***************/
/// xpr matrix multiplication by vector
template<typename T, typename A, typename B>
Xpr< T, 
     XprMatMultOp<T, 
	      MatConstRef<T,MatIndexable<T,A> >,
	      ConstRef<T,Indexable<T,B> > > >
operator*( const MatIndexable<T,A>& a, const Indexable<T,B>& b )
{
  typedef XprMatMultOp< T, MatConstRef<T,MatIndexable<T,A> >,
    ConstRef<T,Indexable<T,B> > > ExprT;
  
  return Xpr< T, ExprT >(
			 ExprT( MatConstRef<T,MatIndexable<T,A> >(a),
				ConstRef<T,Indexable<T,B> >(b) ));
}

/// xpr matrix multiplication by vector specialization
template<typename A, typename B>
Xpr< Complex, 
     XprMatMultOp<Complex, 
	      MatConstRef<Complex,MatIndexable<Complex,A> >,
	      ConstRef<double,Indexable<double,B> > > >
operator*( const MatIndexable<Complex,A>& a, const Indexable<double,B>& b )
{
  typedef XprMatMultOp< Complex, MatConstRef<Complex,MatIndexable<Complex,A> >,
    ConstRef<double,Indexable<double,B> > > ExprT;
  
  return Xpr< Complex, ExprT >(
			 ExprT( MatConstRef<Complex,MatIndexable<Complex,A> >(a),
				ConstRef<double,Indexable<double,B> >(b) ));
}

/// xpr matrix multiplication by vector specialization
template<typename A, typename B>
Xpr< Complex, 
     XprMatMultOp<Complex, 
	      MatConstRef<double,MatIndexable<double,A> >,
	      ConstRef<Complex,Indexable<Complex,B> > > >
operator*( const MatIndexable<double,A>& a, const Indexable<Complex,B>& b )
{
  typedef XprMatMultOp< Complex, MatConstRef<double,MatIndexable<double,A> >,
    ConstRef<Complex,Indexable<Complex,B> > > ExprT;
  
  return Xpr< Complex, ExprT >(
			 ExprT( MatConstRef<double,MatIndexable<double,A> >(a),
				ConstRef<Complex,Indexable<Complex,B> >(b) ));
}

#undef Complex

#endif
