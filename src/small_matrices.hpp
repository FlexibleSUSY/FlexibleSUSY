#ifndef small_matrices_hpp
#define small_matrices_hpp


#include <Eigen/Dense>
#include "mathdefs.hpp"

namespace flexiblesusy {

// NOTE: the following matrices are under column-major storage scheme
// see http://eigen.tuxfamily.org/dox/TopicStorageOrders.html

using RM22 = Eigen::Matrix<Real, 2, 2>;
using RM33 = Eigen::Matrix<Real, 3, 3>;
using RM44 = Eigen::Matrix<Real, 4, 4>;
using RM66 = Eigen::Matrix<Real, 6, 6>;
using CM22 = Eigen::Matrix<Comp, 2, 2>;
using CM33 = Eigen::Matrix<Comp, 3, 3>;
using CM44 = Eigen::Matrix<Comp, 4, 4>;
using CM66 = Eigen::Matrix<Comp, 6, 6>;
using RVe2 = Eigen::Matrix<Real, 2, 1>;
using RVe3 = Eigen::Matrix<Real, 3, 1>;
using RVe4 = Eigen::Matrix<Real, 4, 1>;
using RVe6 = Eigen::Matrix<Real, 6, 1>;

} // namespace flexiblesusy

#endif // small_matrices_hpp
