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

#ifndef COUPLING_MONITOR_H
#define COUPLING_MONITOR_H

#include <vector>
#include <functional>
#include <string>

#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class Coupling_monitor
 * @brief records model parameters at different scales
 */
class Coupling_monitor {
public:
   /// function to return parameter values at given scale
   using Data_getter = std::function<Eigen::ArrayXd(double)>;
   /// tuple of scale and couplings
   using Tuple = std::pair<double, Eigen::ArrayXd>;

   Coupling_monitor(const Data_getter&, const std::vector<std::string>&);

   /// get couplings at all scales
   void run(double, double, int number_of_steps = 20, bool include_endpoint = false);
   /// get maximum scale
   Tuple get_max_scale() const;
   /// delete all saved couplings
   void clear();
   /// write couplings to file
   void write_to_file(const std::string&, bool overwrite = true) const;

private:
   std::vector<Tuple> couplings{}; ///< all couplings at all scales
   Data_getter data_getter{};      ///< returns parameters at given scale
   std::vector<std::string> parameter_names{}; ///< parameter names

   /// write line with parameter names
   void write_parameter_names_line(std::ostream&) const;
   /// write a comment line
   void write_comment_line(std::ostream&) const;
};

} // namespace flexiblesusy

#endif
