// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "unitarity.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <vector>

namespace flexiblesusy {

namespace {

/// groups the indices that are connected by a non-zero matrix element
std::vector<std::vector<Eigen::Index>> find_blocks(Eigen::MatrixXcd const& m)
{
   const Eigen::Index n = m.rows();
   std::vector<Eigen::Index> parent(n);
   std::iota(parent.begin(), parent.end(), Eigen::Index{0});

   const auto find_root = [&parent](Eigen::Index i) {
      while (parent[i] != i) {
         parent[i] = parent[parent[i]];
         i = parent[i];
      }
      return i;
   };

   for (Eigen::Index col = 0; col < n; ++col) {
      for (Eigen::Index row = 0; row < n; ++row) {
         if (row != col && m(row, col) != std::complex<double>(0., 0.)) {
            parent[find_root(row)] = find_root(col);
         }
      }
   }

   std::vector<std::vector<Eigen::Index>> members(n);
   for (Eigen::Index i = 0; i < n; ++i) {
      members[find_root(i)].push_back(i);
   }

   members.erase(
      std::remove_if(members.begin(), members.end(),
                     [](auto const& block) { return block.empty(); }),
      members.end());

   return members;
}

} // anonymous namespace

/**
 * Returns max(|Re(eigenvalues(m))|), or NaN if m has non-finite elements.
 *
 * A reordering of the rows and columns brings m to a block-diagonal form, whose
 * eigenvalues are those of the blocks.  This holds for every square matrix, so
 * m need not be Hermitian.  Scattering matrices have many small blocks, and the
 * cost of the eigenvalue calculation grows with the third power of the size.
 */
double max_abs_re_eigenvalue(Eigen::MatrixXcd const& m)
{
   if (!m.allFinite()) {
      return std::numeric_limits<double>::quiet_NaN();
   }

   double result = 0.;

   for (auto const& block : find_blocks(m)) {
      const auto size = static_cast<Eigen::Index>(block.size());

      if (size == 1) {
         result = std::max(result, std::abs(m(block[0], block[0]).real()));
         continue;
      }

      Eigen::MatrixXcd sub(size, size);
      for (Eigen::Index col = 0; col < size; ++col) {
         for (Eigen::Index row = 0; row < size; ++row) {
            sub(row, col) = m(block[row], block[col]);
         }
      }

      result = std::max(result, sub.eigenvalues().real().cwiseAbs().maxCoeff());
   }

   return result;
}

} // namespace flexiblesusy
