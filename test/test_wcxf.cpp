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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_wcxf

#include <boost/test/unit_test.hpp>

#include "json.hpp"
#include "wcxf.hpp"

#include <cstdio>
#include <fstream>
#include <string>

using namespace flexiblesusy;

namespace
{

const std::string file_name = "test/test_wcxf.json";

nlohmann::json read(const std::string& name)
{
   std::ifstream ifs(name);
   return nlohmann::json::parse(ifs);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_write_new_file)
{
   (void)std::remove(file_name.c_str());

   wcxf::update_json_flavio(file_name, 100., {{"C9_bsmumu", {1., -2.}}});

   const auto j = read(file_name);
   BOOST_CHECK_EQUAL(j["eft"], "WET");
   BOOST_CHECK_EQUAL(j["basis"], "flavio");
   BOOST_CHECK_EQUAL(j["scale"].get<double>(), 100.);
   BOOST_CHECK_EQUAL(j["values"].size(), 1);
   BOOST_CHECK_EQUAL(j["values"]["C9_bsmumu"]["Re"].get<double>(), 1.);
   BOOST_CHECK_EQUAL(j["values"]["C9_bsmumu"]["Im"].get<double>(), -2.);
}

BOOST_AUTO_TEST_CASE(test_merge_same_scale)
{
   (void)std::remove(file_name.c_str());

   wcxf::update_json_flavio(file_name, 100.,
                            {{"C7_bs", {1., 0.}}, {"C9_bsmumu", {2., 0.}}});
   wcxf::update_json_flavio(
      file_name, 100., {{"C9_bsmumu", {3., 0.}}, {"C10_bsmumu", {4., 0.}}});

   const auto j = read(file_name);
   BOOST_CHECK_EQUAL(j["values"].size(), 3);
   BOOST_CHECK_EQUAL(j["values"]["C7_bs"]["Re"].get<double>(), 1.);
   BOOST_CHECK_EQUAL(j["values"]["C9_bsmumu"]["Re"].get<double>(), 3.);
   BOOST_CHECK_EQUAL(j["values"]["C10_bsmumu"]["Re"].get<double>(), 4.);
}

BOOST_AUTO_TEST_CASE(test_overwrite_other_scale)
{
   (void)std::remove(file_name.c_str());

   wcxf::update_json_flavio(file_name, 100., {{"C7_bs", {1., 0.}}});
   wcxf::update_json_flavio(file_name, 200., {{"C9_bsmumu", {2., 0.}}});

   const auto j = read(file_name);
   BOOST_CHECK_EQUAL(j["scale"].get<double>(), 200.);
   BOOST_CHECK_EQUAL(j["values"].size(), 1);
   BOOST_CHECK(!j["values"].contains("C7_bs"));
}

BOOST_AUTO_TEST_CASE(test_overwrite_invalid_file)
{
   {
      std::ofstream ofs(file_name);
      ofs << "no json {";
   }

   wcxf::update_json_flavio(file_name, 100., {{"C7_bs", {1., 0.}}});

   const auto j = read(file_name);
   BOOST_CHECK_EQUAL(j["values"].size(), 1);
   BOOST_CHECK_EQUAL(j["values"]["C7_bs"]["Re"].get<double>(), 1.);

   (void)std::remove(file_name.c_str());
}
