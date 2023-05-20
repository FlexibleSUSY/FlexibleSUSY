
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_slha_output

#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>

#include "slhaea.h"

int run_cmd(const std::string& cmd)
{
   if (!system(NULL)) {
      BOOST_TEST_MESSAGE("Error: command processor not available!");
      return -1;
   }

   BOOST_TEST_MESSAGE("Runnign command \"" << cmd << "\" ...");

   const int status = system(cmd.c_str());

   BOOST_TEST_MESSAGE("Command \"" << cmd << "\" returned with exit code " << status);

   return status;
}

int run_point(const std::string& slha_file,
              const std::string& flexiblesusy_output_file,
              const std::string& softsusy_output_file)
{
   int status;

   status = run_cmd("./models/CMSSM/run_CMSSM.x --slha-input-file=" +
                    slha_file + " --slha-output-file=" + flexiblesusy_output_file);

   if (status) {
      BOOST_TEST_MESSAGE("FlexibleSUSY failed with exit code " << status);
      return status;
   }

   status = run_cmd("./test/SOFTSUSY/run_softpoint.x leshouches < " +
                    slha_file + " > " + softsusy_output_file);

   if (status) {
      BOOST_TEST_MESSAGE("Softsusy failed with exit code " << status);
      return status;
   }

   return 0;
}

void compare_2component_block(const std::string& name,
                              const SLHAea::Coll& coll1,
                              const SLHAea::Coll& coll2,
                              const std::vector<int>& exclude)
{
   BOOST_REQUIRE(coll1.find(name) != coll1.end());
   BOOST_REQUIRE(coll2.find(name) != coll2.end());

   for (SLHAea::Block::const_iterator line = coll1.at(name).begin(),
           end = coll1.at(name).end(); line != end; ++line) {
      if (!line->is_data_line())
         continue;

      if (line->size() < 2)
         continue;

      const std::string index = (*line)[0];

      if (std::find(exclude.begin(), exclude.end(), SLHAea::to<int>(index))
          != exclude.end()) {
         BOOST_TEST_MESSAGE("  Note: skipping key " << index);
         continue;
      }

      const double value = std::fabs(SLHAea::to<double>((*line)[1]));

      // look for index in coll2
      SLHAea::Block::key_type key;
      key.push_back(index);
      SLHAea::Block::const_iterator line2 = coll2.at(name).find(key);
      if (line2 == coll2.at(name).end()) {
         BOOST_TEST_MESSAGE("Error: there is no line with key " << index << " in coll2");
         BOOST_CHECK(false);
         continue;
      }

      if (line2->size() < 2) {
         BOOST_TEST_MESSAGE("Error: line2 (key = " << index << ") contains no value");
         BOOST_CHECK(false);
         continue;
      }

      if (line->size() >= 3 && line2->size() >= 3)
         BOOST_TEST_MESSAGE("  comparing key " << index << " (" << (*line)[2]
                       << " vs " << (*line2)[2] << ")");

      const double value2 = std::fabs(SLHAea::to<double>((*line2)[1]));

      BOOST_CHECK_CLOSE_FRACTION(value, value2, 8.0e-3);
   }
}

void compare_block_gauge(const SLHAea::Coll& coll1, const SLHAea::Coll& coll2)
{
   BOOST_REQUIRE(coll1.find("gauge") != coll1.end());
   BOOST_REQUIRE(coll2.find("gauge") != coll2.end());

   // scale
   BOOST_CHECK_CLOSE_FRACTION(SLHAea::to<double>(coll1.at("gauge").at("Block").at(3)),
                              SLHAea::to<double>(coll2.at("gauge").at("Block").at(3)), 6e-4);

   BOOST_CHECK_CLOSE_FRACTION(SLHAea::to<double>(coll1.at("gauge").at("1").at(1)),
                              SLHAea::to<double>(coll2.at("gauge").at("1").at(1)), 1.0e-4);

   BOOST_CHECK_CLOSE_FRACTION(SLHAea::to<double>(coll1.at("gauge").at("2").at(1)),
                              SLHAea::to<double>(coll2.at("gauge").at("2").at(1)), 0.0003);

   BOOST_CHECK_CLOSE_FRACTION(SLHAea::to<double>(coll1.at("gauge").at("3").at(1)),
                              SLHAea::to<double>(coll2.at("gauge").at("3").at(1)), 3e-5);
}

void compare_slha_files(const std::string& file1, const std::string& file2)
{
   std::ifstream ifs1(file1);
   const SLHAea::Coll input1(ifs1);
   BOOST_REQUIRE(!input1.empty());

   std::ifstream ifs2(file2);
   const SLHAea::Coll input2(ifs2);
   BOOST_REQUIRE(!input2.empty());

   compare_block_gauge(input1, input2);

   // excluding sfermions because their ordering is different in
   // Softsusy
   int excluded_masses[] = {
      1000001, // Sd_1
      1000003, // Sd_2
      1000005, // Sd_3
      2000001, // Sd_4
      2000003, // Sd_5
      2000005, // Sd_6
      1000011, // Se_1
      1000013, // Se_2
      1000015, // Se_3
      2000011, // Se_4
      2000013, // Se_5
      2000015, // Se_6
      1000002, // Su_1
      1000004, // Su_2
      1000006, // Su_3
      2000002, // Su_4
      2000004, // Su_5
      2000006  // Su_6
   };

   const int number_of_excluded_masses
      = sizeof(excluded_masses) / sizeof(excluded_masses[0]);

   std::vector<int> excluded(excluded_masses, excluded_masses
                             + number_of_excluded_masses);

   compare_2component_block("mass", input1, input2, excluded);
}

BOOST_AUTO_TEST_CASE( test_slha_output )
{
   const std::string input_file("test/test_CMSSM_slha_output.in.spc");
   const std::string output_file("test/test_CMSSM_slha_output.out.spc");

   int status = run_point(input_file,
                          output_file + ".fs", output_file + ".ss");

   BOOST_REQUIRE_EQUAL(status, 0);

   compare_slha_files(output_file + ".fs", output_file + ".ss");
}
