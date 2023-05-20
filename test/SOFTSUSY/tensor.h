
/** \file tensor.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief three-index tensor for containing information on R-parity
   violating couplings. Also contains linear algebra functions for interaction
   with matrices and vectors.
*/

#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <cmath>
#include <fstream>
#include "linalg.h"
//#include "rge.h"
//#include "lowe_legacy.h"
//#include "susy.h"

namespace softsusy {

class Tensor;
std::ostream & operator <<(std::ostream &left, const Tensor &);

/// Three-index tensor for containing information on RPV couplings. 
class Tensor {
/// Also contains linear algebra functions for interaction
/// with matrices and vectors. 
/// Specifically: a vector of matrices, specially designed for R-parity
/// violating interactions ie dimensions (3,3,3) 
private:
  DoubleMatrix A1, A2, A3; ///< Three matrices held in memory
  
  /// sets contents equal to the matrices given as arguments
  void setTensor(const DoubleMatrix &, const DoubleMatrix &, const
		 DoubleMatrix &); 
  
public:
  Tensor(); ///< Constructor fills object with zeroes by default
  /// Constructor sets object to be equal to another
  Tensor(const Tensor &);
  
  /// Checks that the three input matrices have dimension 3x3
  void threecheck(const DoubleMatrix &, const DoubleMatrix &, 
		  const DoubleMatrix &); 
  
  /// References a single element of the tensor
  double & operator () (int, int, int);
  /// References a single matrix of the tensor
  DoubleMatrix & operator () (int);
  /// Returns a single matrix of the tensor
  const DoubleMatrix & display(int) const;
  /// Returns a single element of the tensor
  double display(int, int, int) const;
  /// outputs tensors to standard input IF they're elements sum to more than
  /// tol 
  void checkOut(double) const;
  /// Sets a single element of tensor \f$ T_{ijk} \f$=f
  void set(int i, int j, int k, double f);

  /// Does \f$ V^i = T^{jij} \f$ where the input l is the position of i in T
  DoubleVector trace(int) const;
  /// Transposes each matrix held in the tensor
  Tensor transpose() const;
  /// Outputs \f$ T^{ijk} V_i \f$ summing over ith (1st 2nd or 3rd) index. 
  /// eg i=1: \f$ M_{ij} = T^{kij} V_k \f$
  /// 2: \f$ M_{ij} = T^{jki} V_k \f$
  /// 3: \f$ M_{ij} = T^{ijk} V_k \f$
  /// so the matrix indices are just in order from L-R AFTER summed index.
  DoubleMatrix dotProd(const DoubleVector & v, int i) const;
  /// Multiplies all matrices in tensor by a double
  Tensor operator*(double) const;
  /// Divides all matrices in tensor by a double
  Tensor operator/(double) const;
  /// \f$ T^{ijk} = T^{ijl} M_{lk} \f$
  Tensor operator*(const DoubleMatrix &) const;
  /// Adds all matrices between two tensors
  Tensor operator+(const Tensor &) const;
  /// Subtracts all matrices between two tensors
  Tensor operator-(const Tensor &) const;
  /// Does \f$ T^{ijk} = T^{ljk} M_{li} \f$
  Tensor product(const DoubleMatrix &) const;
  /// Swaps the other indices apart ith one eg i=1: \f$ {T^{ijk}}'=T^{ikj} \f$
  Tensor swap(int);
  /// \f$ T^kij=M^kl T^lij \f$
  Tensor raise(const DoubleMatrix & M) const;
};
/// Multiplies all matrices in tensor by a double
Tensor operator*(double, const Tensor &);
/// \f$ T^{kij} = M_{il} T^{klj} \f$
Tensor operator*(const DoubleMatrix &, const Tensor &);
/// Make a tensor out of a vector and a matrix. l labels the position the
/// vector index goes in. After that, indices are cyclic.
Tensor outerProduct(const DoubleVector &, const DoubleMatrix &, int);
/// \f$ M_{ij} = A^{ikl} B^{jlk} \f$
DoubleMatrix matrixify(const Tensor &, const Tensor &);
// \f$ M_{ij} = A^{lik} B^{lkj} \f$
DoubleMatrix sumProd(const Tensor & A, const Tensor & B);

} // namespace softsusy

#endif
