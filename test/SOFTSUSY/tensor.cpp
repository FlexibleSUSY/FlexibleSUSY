
/** \file tensor.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: tensor class - 3 index objects encoded as matrices

*/

#include "tensor.h"

namespace softsusy {

void Tensor::setTensor(const DoubleMatrix & B1, const DoubleMatrix & B2,
		      const DoubleMatrix & B3) {
  threecheck(B1, B2, B3);
  A1 = B1;
  A2 = B2;
  A3 = B3;
}

Tensor::Tensor()
  : A1(3,3), A2(3,3), A3(3,3) {
  int i, j;
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      A1(i, j) = 0.0;
      A2(i, j) = 0.0;
      A3(i, j) = 0.0;
    }
}

Tensor::Tensor(const Tensor &T)
  : A1(T.A1), A2(T.A2), A3(T.A3) {}

double & Tensor::operator() (int i, int j, int k) {
  switch(i) {
  case 1: return A1(j, k); break;
  case 2: return A2(j, k); break;
  case 3: return A3(j, k); break;
  default: 
    std::ostringstream ii;
    ii << "Trying to access " << i << j << k << "th element of tensor";
    ii << *this;
    throw ii.str();
    break;
  }
}

void Tensor::threecheck(const DoubleMatrix &B1, const DoubleMatrix &B2,
			const DoubleMatrix &B3) {
  if (B1.displayRows() != 3 || B1.displayCols() != 3 ||
      B2.displayRows() != 3 || B2.displayCols() != 3 ||
      B3.displayRows() != 3 || B3.displayCols() != 3) {
    std::ostringstream ii;
    ii << "Incorrect tensor initialisation with \n" << B1 << B2 << B3;
    throw ii.str();
  }
}

DoubleMatrix & Tensor::operator() (int i) {
  switch(i) {
  case 1: return A1; break;
  case 2: return A2; break;
  case 3: return A3; break;
  default:  
    std::ostringstream ii;
    ii << "Tensor index out of range: " << i << "\n";
    throw ii.str();
    break;
  }
}

const DoubleMatrix & Tensor::display(int i) const { 
  switch (i) { 
  case 1: return A1; break;
  case 2: return A2; break;
  case 3: return A3; break;
  default:  
    std::ostringstream ii;
    ii << "Tensor index out of range: " << i << "\n";
    throw ii.str();
    break;
  }
}

double Tensor::display(int i, int j, int k) const {
  double ans;
  switch(i) {
  case 1: ans = A1.display(j, k); break;
  case 2: ans = A2.display(j, k); break;
  case 3: ans = A3.display(j, k); break;
  default:  
    std::ostringstream ii;
    ii << "Tensor index out of range: " << i << "\n";
    throw ii.str();
    break;
  }
  return ans;
}

void Tensor::set(int i, int j, int k, double f) {
  switch(i) {
  case 1: A1(j, k) = f; break;
  case 2: A2(j, k) = f; break;
  case 3: A3(j, k) = f; break;
  default:  
    std::ostringstream ii;
    ii << "Tensor index out of range: " << i << "\n";
    throw ii.str();
    break;
  }
}

Tensor Tensor::transpose() const {
  Tensor temp;
  int i;
  for (i=1;i<=3;i++) temp(i) = display(i).transpose();
  return temp;
}

// Does V^i = T^jij where the input l is the position of i in T
DoubleVector Tensor::trace(int l) const {
  DoubleVector V(3);
  int j, k;
  for (k=1;k<=3;k++)
    for (j=1;j<=3;j++) 
      switch(l)	{
      case 1: V(k) += display(k, j, j); break;
      case 2: V(k) += display(j, k, j); break;
      case 3: V(k) += display(j, j, k); break;
      default:
	std::ostringstream ii;
	ii << "Attempting to trace over indices other than " << l << " of " 
	   << *this << "\n";
	throw ii.str();
	break;
      }
  return V;
}


// Outputs T^ijk V_i summing over ith (1st 2nd or 3rd) index. eg i = 
// 1: M_ij = T^kij V_k
// 2: M_ij = T^jki V_k
// 3: M_ij = T^ijk V_k
// so the matrix indices are just in order from L-R AFTER summed index.
DoubleMatrix Tensor::dotProd(const DoubleVector & V, int i) const {
  DoubleMatrix result(3, 3);

  int j, k, l;
  for (j=1; j<=3; j++)
    for (k=1; k<=3; k++)
      switch(i) {
      case 1: 
	for (l=1; l<=3; l++) 
	  result(j, k) += display(l, j, k) * V.display(l);
	break; 
      case 2:
	for (l=1; l<=3; l++) 
	  result(j, k) += display(k, l, j) * V.display(l); 
	break;
      case 3:
	for (l=1; l<=3; l++) 
	  result(j, k) += display(j, k, l) * V.display(l); 
	break;
      default:
	std::ostringstream ii;
	ii << "sum out of range in dot product " << *this << "*" 
	   << V << "(" << V.displayStart() << "," << V.displayEnd() << ")"
	   << " on " << i << "th index.\n"; 
	throw ii.str();
      }
  return result;
}

Tensor Tensor::operator*(double f) const {
  Tensor T(*this);
  int i; for(i=1; i<=3; i++)
    T(i) = T(i) * f;
  return T;
}

Tensor Tensor::operator/(double f) const {
  Tensor T(*this);
  return T * (1.0 / f);
}

// Swaps the other indices apart ith one eg i=1: T^ijk'=T^ikj
Tensor Tensor::swap(int i) {
  
  Tensor T;
  int j,k,l;
  for(j=1; j<=3; j++)
    for(k=1; k<=3; k++)
      for(l=1; l<=3; l++) {
	switch(i) {
	case 1: T(j, k, l) = display(j, l, k); break;
	case 2: T(j, k, l) = display(l, k, j); break;
	case 3: T(j, k, l) = display(k, j, l); break;
	default: 
	  std::ostringstream ii;
	  ii << "Trying to swap around " << i << "th element of tensor";
	  throw ii.str();
	  break;
	}
      }
  return T;
}

// T^ijk = T^ijl M_lk
Tensor Tensor::operator*(const DoubleMatrix & M) const {
  if (M.displayRows() != 3 || M.displayCols() !=3) {
    std::ostringstream ii;
    ii << "Tensor " << *this << " * matrix" << M << "of wrong size.\n";
    throw ii.str();
  }
  Tensor T;
  int i; 
  for (i=1; i<=3; i++) T(i) = display(i) * M;
  return T;
}

// Does T^ijk = T^ljk M_li
Tensor Tensor::product(const DoubleMatrix & M) const {
  if (M.displayRows() != 3 || M.displayCols() !=3) {
    std::ostringstream ii;
    ii << "Tensor " << *this << " * matrix" << M << "of wrong size.\n";
    throw ii.str();
  }
  Tensor T;
  int i,j,k,l; 
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++)
      for (k=1; k<=3; k++) 
	for (l=1; l<=3; l++) 
	  T(i, j, k) = T(i, j, k) + display(l, j, k) * M.display(l, i);
  return T;
}

// outputs tensors IF they're elements sum to more than tol
void Tensor::checkOut(double tol) const {
  // Tolerance: if matrices have elements of less than this number, it won't
  // display them!
  //  double tol = 1.0e-14;
  int i; 

  for (i=1; i<=3; i++)
    if (display(i).sumElements() > tol)
	std::cout << "Matrix " << i << display(i);
}

std::ostream & operator << (std::ostream &left, const Tensor &T) {
  int i; 
  for (i=1; i<=3; i++)
    left << "Matrix " << i << T.display(i);
  return left;
}

Tensor Tensor::operator+(const Tensor &T) const {
  Tensor temp;
  int i; for(i=1; i<=3; i++) temp(i) = T.display(i) + display(i);
  return temp;
}

Tensor Tensor::operator-(const Tensor &T) const {
  Tensor temp;
  int i; for(i=1; i<=3; i++) temp(i) = display(i) - T.display(i);
  return temp;
}

// l labels the position the vector index goes in.
// After that, indices are cyclic.
Tensor outerProduct(const DoubleVector &V, const DoubleMatrix & M, int l) {
  Tensor temp;
  int i, j, k;
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++)
      for (k=1; k<=3; k++) {
	switch(l) {
	case 1: temp(i, j, k) = V.display(i) * M.display(j, k); break;
	case 2: temp(i, j, k) = V.display(j) * M.display(k, i); break;
	case 3: temp(i, j, k) = V.display(k) * M.display(i, j); break;
	default: 
	  std::ostringstream ii;
	  ii << "Trying to outer product " << l << "th element of tensor";
	  throw ii.str();
	  break;
	}
      }
  return temp;
}

// M_ij = A^ikl B^jlk
DoubleMatrix matrixify(const Tensor & A, const Tensor &B) {
  DoubleMatrix temp(3, 3);
  int i, j; 
  for (i=1; i<=3; i++) 
    for (j=1; j<=3; j++)
      temp(i, j) +=  (A.display(i) * B.display(j)).trace();
  return temp;
}

// M_ij = A^lik B^lkj
DoubleMatrix sumProd(const Tensor & A, const Tensor & B) {
  DoubleMatrix temp(3, 3);
  int i; for(i=1; i<=3; i++) temp = temp + A.display(i) * B.display(i);
  return temp;
}

// T^ijk = f * T^ijk
Tensor operator*(double f, const Tensor &M) {
  Tensor T(M);
  int i; for(i=1; i<=3; i++)
    T(i) = T(i) * f;
  return T;
}

// T^kij=M^kl T^lij
Tensor Tensor::raise(const DoubleMatrix & M) const {
  Tensor temp;
  int i,k,j,l; 
    for(k=1; k<=3; k++)
    for(i=1; i<=3; i++) 
      for(j=1; j<=3; j++)
	for(l=1; l<=3; l++)
	  temp(k, i, j) += M.display(k, l) * this -> display(l, i, j);
    return temp;
}

// T^kij = M_il T^klj
Tensor operator*(const DoubleMatrix &M, const Tensor &T) {
  if (M.displayRows() !=3 || M.displayCols() != 3) {
    std::ostringstream ii;
    ii << "DoubleMatrix " << M << " * Tensor" << T << "of wrong size\n";
    throw ii.str();
  }
  Tensor temp;
  int i,j,k,l; 
  for(k=1; k<=3; k++)
    for(i=1; i<=3; i++) 
      for(j=1; j<=3; j++)
	for(l=1; l<=3; l++)
	  temp(k, i, j) += M.display(i, l) * T.display(k, l, j);
  return temp;  
}


} // namespace softsusy
