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

#include "coupling_monitor.hpp"

#include "error.hpp"
#include "logger.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace flexiblesusy {

namespace {

const int column_width{16}; ///< width of columns in output table

struct Scale_comp {
   bool operator() (const Coupling_monitor::Tuple& i,
                    const Coupling_monitor::Tuple& j) const {
      return i.first < j.first;
   }
};

} // anonymous namespace

Coupling_monitor::Coupling_monitor(
   const Data_getter& data_getter_,
   const std::vector<std::string>& parameter_names_)
   : data_getter(data_getter_), parameter_names(parameter_names_)
{
}

/**
 * Get the couplings at the largest scale
 *
 * @return a pair with the scale and a Eigen::ArrayXd which contains the
 * couplings at this scale
 */
Coupling_monitor::Tuple Coupling_monitor::get_max_scale() const
{
   if (couplings.empty()) {
      ERROR("Data container is empty!");
      return Tuple(0.0, Eigen::ArrayXd(1));
   }

   // find gauge couplings at the greatest scale
   auto maxScale
      = max_element(couplings.begin(), couplings.end(), Scale_comp());

   return *maxScale;
}

/**
 * Delete all internal couplings.
 */
void Coupling_monitor::clear()
{
   couplings.clear();
}

/**
 * write line with parameter names
 *
 * @param fout output stream
 */
void Coupling_monitor::write_parameter_names_line(std::ostream& fout) const
{
   if (!fout.good() || couplings.empty()) {
      return;
   }

   fout << std::left << std::setw(column_width) << "scale";

   for (const auto& p: parameter_names) {
      fout << std::left << std::setw(column_width) << p;
   }

   fout << '\n';
}

/**
 * write help line which describes the written data
 *
 * @param fout output stream
 */
void Coupling_monitor::write_comment_line(std::ostream& fout) const
{
   if (!fout.good() || couplings.empty())
      return;

   fout << std::left << std::setw(column_width) << "# [1] scale";

   for (std::size_t i = 0; i < parameter_names.size(); ++i) {
      fout << std::left << std::setw(column_width)
           << '[' + std::to_string(i+2) + "] " + parameter_names[i];
   }

   fout << '\n';
}

/**
 * Write all couplings to a text file.
 *
 * @param file_name name of file to write the data to
 * @param overwrite if true, file is overwritten, otherwise content is appended
 */
void Coupling_monitor::write_to_file(const std::string& file_name, bool overwrite) const
{
   if (couplings.empty())
      return;

   const std::ios_base::openmode openmode
      = (overwrite ? std::ios::out : std::ios::app);

   std::ofstream filestr(file_name, openmode);
   VERBOSE_MSG("Coupling_monitor<>::write_to_file: opening file: "
               << file_name);
   if (filestr.fail()) {
      ERROR("can't open file " << file_name
            << " for writing running couplings");
      return;
   }

   write_comment_line(filestr);
   write_parameter_names_line(filestr);

   // write data
   for (const auto& c: couplings) {
      if (!filestr.good()) {
         ERROR("file " << file_name << " is corrupted");
         break;
      }

      filestr << std::left << std::setw(column_width) << c.first;

      // write all gauge couplings in order
      for (int i = 0; i < c.second.size(); ++i) {
         filestr << std::left << std::setw(column_width) << c.second(i);
      }

      filestr << '\n';
   }

   filestr.close();
   VERBOSE_MSG("Coupling_monitor<>::write_to_file: file written: "
               << file_name);
}

/**
 * Add running couplings between scale q1 and q2.
 *
 * @param q1 scale to start at
 * @param q2 end scale
 * @param number_of_steps number of steps
 * @param include_endpoint include the endpoint q2 in the running
 *        (false by default)
 */
void Coupling_monitor::run(double q1, double q2,
                           int number_of_steps, bool include_endpoint)
{
   if (q1 <= 0.0 || q2 <= 0.0) {
      ERROR("negative scales are not allowed: q1=" << q1 << ", q2=" << q2);
      return;
   }

   if (number_of_steps < 1) {
      number_of_steps = 1;
   }

   // if the endpoint should be included, the scale loop must run from
   // (n == 0) to (n == number_of_steps); otherwise it runs from (n == 0) to (n
   // == number_of_steps - 1)
   const int endpoint_offset = include_endpoint ? 1 : 0;

   // run from q1 to q2
   for (int n = 0; n < number_of_steps + endpoint_offset; ++n) {
      const double lq1 = std::log(q1);
      const double lq2 = std::log(q2);
      const double scale = std::exp(lq1 + n * (lq2 - lq1) / number_of_steps);
      try {
         couplings.emplace_back(scale, data_getter(scale));
      } catch (const Error&) {
         ERROR("Coupling_monitor::run: run to scale "
               << scale << " failed");
         break;
      }
   }

   std::sort(couplings.begin(), couplings.end(), Scale_comp());
}

} // namespace flexiblesusy
