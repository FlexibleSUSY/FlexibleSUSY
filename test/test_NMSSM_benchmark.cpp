
#include <iostream>
#include <fstream>
#include <string>

#include "test.hpp"
#include "run_cmd.hpp"
#include "slhaea.h"
#include "stopwatch.hpp"
#include "logger.hpp"

struct Data {
   Data() : number_of_valid_points(0), sum_of_times(0.) {}

   double time;
   int error;
   int number_of_valid_points;
   double sum_of_times;
};

bool is_valid_spectrum(const std::string& slha_file)
{
   std::ifstream ifs(slha_file);
   SLHAea::Coll coll(ifs);

   // find SPINFO block
   SLHAea::Coll::const_iterator block =
      coll.find(coll.cbegin(), coll.cend(), "SPINFO");

   if (block == coll.cend())
      throw std::string("Error: SPINFO block not found in file ") + slha_file;

   for (SLHAea::Block::const_iterator line = block->cbegin(),
           end = block->cend(); line != end; ++line) {
      if (line->is_data_line() && line->size() >= 2 &&
          (*line)[0] == "4" && (*line)[1] != "")
         return false;
   }

   return true;
}

void run_point(const std::string& slha_file, Data& fs_data, Data& ss_data)
{
   int status;
   flexiblesusy::Stopwatch stopwatch;

   const std::string slha_output_file("test/test_NMSSM_benchmark.out.spc");

   stopwatch.start();
   status = run_cmd("./models/NMSSM/run_NMSSM.x --slha-input-file=" +
                    slha_file + " --slha-output-file=" + slha_output_file +
                    " > /dev/null 2>&1");
   stopwatch.stop();

   fs_data.time = stopwatch.get_time_in_seconds();
   fs_data.sum_of_times += fs_data.time;
   fs_data.error = status;

   if (!fs_data.error) {
      // look for errors in the SLHA output file
      fs_data.error = !is_valid_spectrum(slha_output_file);
      if (!fs_data.error)
         fs_data.number_of_valid_points++;
   }

   stopwatch.start();
   status = run_cmd("./test/SOFTSUSY/run_softpoint.x leshouches < " +
                    slha_file + " > " + slha_output_file);
   stopwatch.stop();

   ss_data.time = stopwatch.get_time_in_seconds();
   ss_data.sum_of_times += ss_data.time;
   ss_data.error = status;

   if (!ss_data.error) {
      // look for errors in the SLHA output file
      ss_data.error = !is_valid_spectrum(slha_output_file);
      if (!ss_data.error)
         ss_data.number_of_valid_points++;
   }
}

SLHAea::Coll create_point(double tanBeta)
{
   std::ifstream ifs("test/test_NMSSM_benchmark.in.spc.in");
   SLHAea::Coll coll(ifs);
   SLHAea::Block minpar, extpar;

   const std::string minpar_str(
      "Block MINPAR\n"
      "   1   2.000000000e+02   # m0\n"
      "   2   5.000000000e+02   # m12\n"
      "   3   " + std::to_string(tanBeta) + "   # TanBeta\n"
      "   5  -5.000000000e+02   # A0\n");

   const std::string extpar_str(
      "Block EXTPAR\n"
      "  61   0.1               # LambdaInput\n");

   minpar.str(minpar_str);
   coll.push_back(minpar);

   extpar.str(extpar_str);
   coll.push_back(extpar);

   return coll;
}

void test_tanbeta_scan()
{
   const int num_points = 100;
   const double tanBeta_start = 2.;
   const double tanBeta_stop = 80.;
   const double tanBeta_step = (tanBeta_stop - tanBeta_start) / num_points;

   Data fs_data, ss_data;

   printf("%10s %30s %30s %20s\n", "tan(beta)",
          "Softsusy / s (status)",
          "FlexibleSUSY / s (status)",
          "FlexibleSUSY speed");

   for (int i = 0; i < num_points; i++) {
      const double tanBeta = tanBeta_start + i * tanBeta_step;
      const SLHAea::Coll coll(create_point(tanBeta));
      const std::string input_file("test/test_NMSSM_benchmark.in.spc");

      std::ofstream ofs(input_file);
      ofs << coll;
      ofs.close();

      run_point(input_file, fs_data, ss_data);

      printf("%10g %22g (%5d) %24g (%3d) %19.1f%%\n", tanBeta,
             ss_data.time, ss_data.error,
             fs_data.time, fs_data.error,
             100.*(fs_data.time - ss_data.time)/ss_data.time);
   }

   const double fs_average_time = fs_data.sum_of_times / fs_data.number_of_valid_points;
   const double ss_average_time = ss_data.sum_of_times / ss_data.number_of_valid_points;

   INFO("Summary: average times (in seconds) \n"
        "  FlexibleSUSY: " << fs_average_time <<
        " (" << fs_data.number_of_valid_points << "/" << num_points << " points)\n" <<
        "  Softsusy    : " << ss_average_time <<
        " (" << ss_data.number_of_valid_points << "/" << num_points << " points)\n" <<
        "  FlexibleSUSY speed: " << (100.*(fs_average_time - ss_average_time)/ss_average_time) << "%");

   TEST_GREATER(ss_average_time, fs_average_time);
}

int main()
{
   test_tanbeta_scan();

   return flexiblesusy::get_errors();
}
