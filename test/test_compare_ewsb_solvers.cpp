#define BOOST_TEST_IGNORE_NON_ZERO_CHILD_CODE

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_compare_ewsb_solvers

#include "run_cmd.hpp"
#include "slhaea.h"
#include "stopwatch.hpp"
#include <string>
#include <fstream>
#include <boost/test/unit_test.hpp>

struct Solvers {
   const std::string name;
   const std::string executable;
   double runtime;
   int status;
   int number_of_points;
   double runtime_sum;
};

struct Point {
   double m0, m12, tanBeta, a0;
   int signMu;
};

int run_point(Solvers solvers[], int number_of_solvers,
              const std::string& slha_file, const std::string& output_file)
{
   for (int i = 0; i < number_of_solvers; i++) {
      flexiblesusy::Stopwatch stopwatch;

      stopwatch.start();

      solvers[i].status = run_cmd(
         solvers[i].executable +
         " --slha-input-file=" + slha_file +
         " --slha-output-file=" + output_file + " > /dev/null 2>&1");

      stopwatch.stop();

      if (solvers[i].status == 0) {
         solvers[i].runtime = stopwatch.get_time_in_seconds();
         solvers[i].number_of_points++;
         solvers[i].runtime_sum += solvers[i].runtime;
      } else {
         solvers[i].runtime = 0;
      }
   }

   return 0;
}

SLHAea::Coll create_point(Point point, const std::string& input_file)
{
   std::ifstream ifs(input_file);

   if (!ifs.good()) {
      throw std::string("Error: could not open file ") + input_file;
   }

   SLHAea::Coll coll(ifs);
   SLHAea::Block minpar;

   const std::string str(
      "Block MINPAR\n"
      "   1   " + std::to_string(point.m0)      + "   # m0\n"
      "   2   " + std::to_string(point.m12)     + "   # m12\n"
      "   3   " + std::to_string(point.tanBeta) + "   # TanBeta\n"
      "   4   " + std::to_string(point.signMu)  + "   # sign(Mu)\n"
      "   5   " + std::to_string(point.a0)      + "   # A0\n");

   minpar.str(str);
   coll.push_back(minpar);

   return coll;
}

BOOST_AUTO_TEST_CASE( test_tanbeta_scan )
{
   Point point = {
      /* .m0      = */ 2000.,
      /* .m12     = */ 5000.,
      /* .tanBeta = */ 10.,
      /* .a0      = */ -3000.,
      /* .signMu  = */ 1
   };

   const int num_points = 20;
   const double tanBeta_start = 2.;
   const double tanBeta_stop = 50.;
   const double tanBeta_step = (tanBeta_stop - tanBeta_start) / num_points;

   Solvers solvers[] = {
      {"GSLHybrid"  , "./models/CMSSMGSLHybrid/run_CMSSMGSLHybrid.x"    , 0., 0},
      {"GSLHybridS" , "./models/CMSSMGSLHybridS/run_CMSSMGSLHybridS.x"  , 0., 0},
      {"GSLBroyden" , "./models/CMSSMGSLBroyden/run_CMSSMGSLBroyden.x"  , 0., 0},
      {"GSLNewton"  , "./models/CMSSMGSLNewton/run_CMSSMGSLNewton.x"    , 0., 0},
      {"FPIRelative", "./models/CMSSMFPIRelative/run_CMSSMFPIRelative.x", 0., 0},
      {"FPIAbsolute", "./models/CMSSMFPIAbsolute/run_CMSSMFPIAbsolute.x", 0., 0},
      {"FPITadpole" , "./models/CMSSMFPITadpole/run_CMSSMFPITadpole.x"  , 0., 0},
   };

   const int number_of_solvers = sizeof(solvers)/sizeof(*solvers);

   printf("m0 = %g, m12 = %g, signMu = %d, a0 = %g\n",
          point.m0, point.m12, point.signMu, point.a0);

   printf("%12s", "tan(beta)");
   for (int i = 0; i < number_of_solvers; i++)
      printf("%12s", solvers[i].name.c_str());
   printf("\n");

   for (int i = 0; i < num_points; i++) {
      point.tanBeta = tanBeta_start + i * tanBeta_step;
      const SLHAea::Coll coll(create_point(point, "test/test_compare_ewsb_solvers.in.spc.in"));
      const std::string input_file("test/test_compare_ewsb_solvers.in.spc");

      std::ofstream ofs(input_file);
      ofs << coll;
      ofs.close();

      run_point(solvers, number_of_solvers,
                input_file, "test/test_compare_ewsb_solvers.out.spc");

      printf("%12g", point.tanBeta);
      for (int i = 0; i < number_of_solvers; i++) {
         if (solvers[i].status == 0)
            printf("%12g", solvers[i].runtime);
         else
            printf("    %2s(%4d)", "- ", solvers[i].status);
      }
      printf("\n");
   }

   // print average times
      printf("%12s", "average");
   for (int i = 0; i < number_of_solvers; i++) {
      const int number_of_points = solvers[i].number_of_points;
      if (number_of_points > 0)
         printf("%12g", solvers[i].runtime_sum / number_of_points);
      else
         printf("%12g", "-");
   }
   printf("\n");
}
