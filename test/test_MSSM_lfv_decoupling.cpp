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
#define BOOST_TEST_MODULE test_MSSM_lfv_decoupling

#include <boost/test/unit_test.hpp>

#include "MSSM_observables.hpp"
#include "MSSM_slha_io.hpp"
#include "MSSM_two_scale_spectrum_generator.hpp"

#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "observables/l_to_l_conversion/settings.hpp"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"

#include <sstream>
#include <string>
#include <tuple>

using namespace flexiblesusy;

namespace
{

/// CMSSM-like point: m0 = m1/2 = M at the GUT scale, A0 = 0,
/// tan(beta) = 10, mu > 0, slepton flavour mixing (ml2)_12 = (me2)_12 = 0.1 M^2
std::string slha_input(double M)
{
   const double M2 = M * M;
   std::ostringstream s;
   s << "Block MODSEL\n"
        "    6   1\n"
        "Block FlexibleSUSY\n"
        "    0   1.0e-04\n"
        "    4   1\n"
        "    5   1\n"
        "    6   2\n"
        "    7   1\n"
        "   15   1\n"
        "   31   1\n"
        "Block MINPAR\n"
        "    3   10\n"
        "    4   1\n"
        "Block EXTPAR\n"
        "    0   1.9422e+16\n"
        "   21   "
     << M2
     << "\n"
        "   22   "
     << M2
     << "\n"
        "Block MSOFTIN\n"
        "    1   "
     << M
     << "\n"
        "    2   "
     << M
     << "\n"
        "    3   "
     << M << "\n";
   for (const char* block :
        {"MSQ2IN", "MSU2IN", "MSD2IN", "MSL2IN", "MSE2IN"}) {
      const bool sleptons =
         std::string(block) == "MSL2IN" || std::string(block) == "MSE2IN";
      s << "Block " << block << '\n';
      for (int i = 1; i <= 3; i++) {
         for (int k = 1; k <= 3; k++) {
            const double v =
               i == k ? M2 : (sleptons && i + k == 3 ? 0.1 * M2 : 0.);
            s << "  " << i << "  " << k << "  " << v << '\n';
         }
      }
   }
   return s.str();
}

/// returns BR(mu -> 3e) and CR(mu -> e, Al)
std::tuple<double, double> lfv_observables(double M)
{
   std::istringstream istr(slha_input(M));
   MSSM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   MSSM_input_parameters input;
   Spectrum_generator_settings settings;
   slha_io.fill(settings);
   slha_io.fill(qedqcd);
   slha_io.fill(input);
   settings.set(Spectrum_generator_settings::calculate_observables, 0);
   Loop_library::set(1);

   MSSM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   if (spectrum_generator.get_problems().have_problem()) {
      std::ostringstream msg;
      spectrum_generator.get_problems().print_problems(msg);
      BOOST_TEST_MESSAGE(msg.str());
   }
   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   const auto model = std::get<0>(spectrum_generator.get_models_slha());
   const auto observables = calculate_observables(
      model, qedqcd, LToLConversion_settings{}, Physical_input{}, settings);

   return {observables.Fe2to11bar1_All_1loop(0).real(),
           observables.Fe2Fe1inAl_All_1loop(0).real()};
}

} // anonymous namespace

// The new physics contributions to mu -> 3e and mu -> e conversion must
// decouple as 1/M^4. Until 2026 the Z penguins didn't decouple: from
// M = 3 TeV to 10 TeV, CR(mu -> e) stayed at 1.4e-12, BR(mu -> 3e) fell
// only by a factor 6.
BOOST_AUTO_TEST_CASE(test_decoupling)
{
   const auto [br3, cr3] = lfv_observables(3000.);
   const auto [br10, cr10] = lfv_observables(10000.);

   BOOST_TEST_MESSAGE("BR(mu -> 3e): " << br3 << " -> " << br10);
   BOOST_TEST_MESSAGE("CR(mu -> e, Al): " << cr3 << " -> " << cr10);

   BOOST_CHECK_GT(br3, 0.);
   BOOST_CHECK_GT(cr3, 0.);

   // (3/10)^4 = 0.0081
   BOOST_CHECK_LT(br10 / br3, 0.02);
   BOOST_CHECK_LT(cr10 / cr3, 0.02);
}
