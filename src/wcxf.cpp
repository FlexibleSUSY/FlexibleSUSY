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

#include "wcxf.hpp"
#include "logger.hpp"

#include "json.hpp"

#include <fstream>
#include <iomanip>
#include <string>

namespace flexiblesusy
{
namespace wcxf
{

namespace
{

constexpr const char* EFT = "WET";
constexpr const char* BASIS = "flavio";

/// returns the content of the file if it holds WET Wilson coefficients
/// in the flavio basis at the given scale, an empty object otherwise
nlohmann::json read_same_scale(const std::string& file_name, double scale)
{
   std::ifstream ifs(file_name);
   if (!ifs.good()) {
      return {};
   }

   const auto j = nlohmann::json::parse(ifs, nullptr, false);

   if (j.is_discarded() || !j.is_object() || j.value("eft", "") != EFT ||
       j.value("basis", "") != BASIS || !j.contains("scale") ||
       !j["scale"].is_number() || j["scale"].get<double>() != scale ||
       !j.contains("values") || !j["values"].is_object()) {
      return {};
   }

   return j;
}

} // anonymous namespace

void update_json_flavio(const std::string& file_name, double scale,
                        const Values& values)
{
   auto j = read_same_scale(file_name, scale);

   if (j.empty()) {
      j["eft"] = EFT;
      j["basis"] = BASIS;
      j["scale"] = scale;
      j["values"] = nlohmann::json::object();
   }

   for (const auto& [name, value] : values) {
      j["values"][name] = {{"Re", value.real()}, {"Im", value.imag()}};
   }

   std::ofstream ofs(file_name);
   ofs << std::setw(4) << j << '\n';

   if (!ofs.good()) {
      ERROR("Cannot write Wilson coefficients to file " << file_name);
   }
}

} // namespace wcxf
} // namespace flexiblesusy
