// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"
#include "gm2calc/gm2_error.hpp"
#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_version.h"
#include "gm2calc/MSSMNoFV_onshell.hpp"

#include "gm2_1loop_helpers.hpp"
#include "gm2_2loop_helpers.hpp"
#include "gm2_config_options.hpp"
#include "gm2_log.hpp"
#include "gm2_slha_io.hpp"

#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <utility>

namespace {

/// reader function
using Reader = std::function<void(gm2calc::MSSMNoFV_onshell&,
                                  const gm2calc::GM2_slha_io& slha_io)>;

/// writer function
using Writer = std::function<void(const gm2calc::MSSMNoFV_onshell&,
                                  const gm2calc::Config_options& options,
                                  gm2calc::GM2_slha_io& slha_io)>;

/**
 * @class Gm2_cmd_line_options
 * @brief command line options for GM2Calc
 */
struct Gm2_cmd_line_options {
   enum E_input_type { SLHA, GM2Calc };

   std::string input_source; ///< input source (file name or `-' for stdin)
   E_input_type input_type{SLHA}; ///< input format (SLHA or GM2Calc)

   static bool starts_with(const std::string& str, const std::string& prefix) {
      return str.compare(0, prefix.size(), prefix) == 0;
   }
};

void print_usage(const char* program_name)
{
   std::cout <<
      "Usage: " << program_name << " [options]\n"
      "Options:\n"
      "  --slha-input-file=<source>      SLHA input source (file name or - for stdin)\n"
      "  --gm2calc-input-file=<source>   GM2Calc input source (file name or - for stdin)\n"
      "  --help,-h                       print this help message\n"
      "  --version,-v                    print version number"
             << std::endl;
}

/**
 * Parses command line options
 *
 * @param argc number of command line arguments
 * @param argv array of command line arguments
 *
 * @return object with extracted information
 */
Gm2_cmd_line_options get_cmd_line_options(int argc, const char* argv[])
{
   Gm2_cmd_line_options options;

   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);

      if (Gm2_cmd_line_options::starts_with(option, "--slha-input-file=")) {
         options.input_source = option.substr(18);
         options.input_type = Gm2_cmd_line_options::SLHA;
         continue;
      }

      if (Gm2_cmd_line_options::starts_with(option, "--gm2calc-input-file=")) {
         options.input_source = option.substr(21);
         options.input_type = Gm2_cmd_line_options::GM2Calc;
         continue;
      }

      if (option == "--help" || option == "-h") {
         print_usage(argv[0]);
         exit(EXIT_SUCCESS);
      }

      if (option == "--version" || option == "-v") {
         std::cout << GM2CALC_VERSION << '\n';
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }

   return options;
}

/**
 * Set the config options to default values, depending on the input
 * parameter set (chosen by the user).
 *
 * If SLHA input format has been selected, the default output format
 * will be SLHA format.  By default \f$a_\mu\f$ will be written to the
 * block GM2Calc[0].
 *
 * If GM2Calc input format has been chosen, the default values set in
 * \a Config_options are used.
 *
 * @param config_options configuration options
 * @param options command line options
 */
void set_to_default(gm2calc::Config_options& config_options,
                    const Gm2_cmd_line_options& options)
{
   switch (options.input_type) {
   case Gm2_cmd_line_options::SLHA:
      config_options.output_format = gm2calc::Config_options::GM2Calc;
      break;
   case Gm2_cmd_line_options::GM2Calc:
      config_options.output_format = gm2calc::Config_options::Detailed;
      break;
   default:
      throw gm2calc::ESetupError("Unknown input option");
      break;
   }
}

/**
 * Prints output if an error has occured.
 *
 * @param error error object
 * @param slha_io SLHA object
 * @param config_options configuration options
 */
void print_error(const gm2calc::Error& error,
                 gm2calc::GM2_slha_io& slha_io,
                 const gm2calc::Config_options& config_options)
{
   switch (config_options.output_format) {
   case gm2calc::Config_options::NMSSMTools:
   case gm2calc::Config_options::SPheno:
   case gm2calc::Config_options::GM2Calc:
      // print SPINFO block with error description
      slha_io.fill_block_entry("SPINFO", 1, "GM2Calc");
      slha_io.fill_block_entry("SPINFO", 2, GM2CALC_VERSION);
      slha_io.fill_block_entry("SPINFO", 4, error.what());
      slha_io.write_to_stream(std::cout);
      break;
   default:
      ERROR(error.what());
      break;
   }
}

/**
 * Calculates a_mu for a given set of configuration options (loop
 * order, tan(beta) resummation).
 *
 * @param model model (must be initialized)
 * @param options configuration options
 *
 * @return a_mu
 */
double calculate_amu(const gm2calc::MSSMNoFV_onshell& model,
                     const gm2calc::Config_options& options)
{
   double result = 0.0;

   if (options.tanb_resummation) {
      if (options.loop_order > 0) {
         result += gm2calc::calculate_amu_1loop(model);
      }
      if (options.loop_order > 1) {
         result += gm2calc::calculate_amu_2loop(model);
      }
   } else {
      // no tan(beta) resummation
      if (options.loop_order > 0) {
         result += gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model);
      }
      if (options.loop_order > 1) {
         result += gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model);
      }
   }

   return result;
}

/**
 * Calculates uncertainty of a_mu for a given set of configuration
 * options (loop order, tan(beta) resummation).
 *
 * @param model model (must be initialized)
 * @param options configuration options
 *
 * @return a_mu
 */
double calculate_uncertainty(const gm2calc::MSSMNoFV_onshell& model,
                             const gm2calc::Config_options& options)
{
   double result = std::numeric_limits<double>::signaling_NaN();

   switch (options.loop_order) {
   case 0:
      result = gm2calc::calculate_uncertainty_amu_0loop(model);
      break;
   case 1:
      result = gm2calc::calculate_uncertainty_amu_1loop(model);
      break;
   case 2:
      result = gm2calc::calculate_uncertainty_amu_2loop(model);
      break;
   default:
      ERROR("loop order > 2 not supported!");
      break;
   }

   return result;
}

/**
 * Reads parameters from SLHA i/o object (SLHA scheme) and initializes
 * model accordingly.
 *
 * @param model model to initialize
 * @param slha_io SLHA i/o object to read parameters from
 */
struct SLHA_reader {
   void operator()(gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::GM2_slha_io& slha_io)
   {
      slha_io.fill_slha(model);
      model.convert_to_onshell();
   }
};

/**
 * Reads parameters from SLHA i/o object in GM2Calc-specific input
 * scheme and initializes model accordingly.
 *
 * @param model model to initialize
 * @param slha_io SLHA i/o object to read parameters from
 */

struct GM2Calc_reader {
   void operator()(gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::GM2_slha_io& slha_io)
   {
      slha_io.fill_gm2calc(model);
      model.calculate_masses();
   }
};

/**
 * Prints a_mu (or the uncertainty) to stdout.
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
struct Minimal_writer {
   void operator()(const gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::Config_options& options,
                   gm2calc::GM2_slha_io&)
   {
      const auto value = options.calculate_uncertainty
                            ? calculate_uncertainty(model, options)
                            : calculate_amu(model, options);

      std::cout << boost::format("%.8e") % value << '\n';
   }
};

/**
 * Prints detailed a_mu calculation (1-loop w/ and w/o tan(beta)
 * resummation, 2-loop, and different contributions).
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
struct Detailed_writer {
   void operator()(const gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::Config_options&, gm2calc::GM2_slha_io&)
   {
#define FORMAT_AMU(amu) boost::format("% 14.8e") % (amu)
#define FORMAT_DEL(amu) boost::format("%14.8e") % (amu)
#define FORMAT_PCT(pct) boost::format("%2.1f") % (pct)

      const std::string error_str = model.get_problems().have_problem()
                                       ? model.get_problems().get_problems() +
                                            " (with tan(beta) resummation)\n\n"
                                       : "";

      const double amu_1l = gm2calc::calculate_amu_1loop(model);
      const double amu_2l_photonic_chipm = gm2calc::amu2LChipmPhotonic(model);
      const double amu_2l_photonic_chi0 = gm2calc::amu2LChi0Photonic(model);
      const double amu_2l_a_sfermion = gm2calc::amu2LaSferm(model);
      const double amu_2l_a_cha = gm2calc::amu2LaCha(model);
      const double amu_2l_ferm_sferm_approx = gm2calc::amu2LFSfapprox(model);
      const double amu_2l = gm2calc::calculate_amu_2loop(model);
      const double amu_2l_uncertainty =
         gm2calc::calculate_uncertainty_amu_2loop(model);
      const double tan_beta_cor = gm2calc::tan_beta_cor(model);

      // no tan(beta) resummation
      double amu_1l_non_tan_beta_resummed = 0.;
      double amu_2l_non_tan_beta_resummed = 0.;
      std::string error_str_non_tan_beta_resummation;

      try {
         // w/o tan(beta) resummation, allow throwing exceptions
         gm2calc::MSSMNoFV_onshell model_except(model);
         model_except.do_force_output(false);
         amu_1l_non_tan_beta_resummed =
            gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model_except);
         amu_2l_non_tan_beta_resummed =
            gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model_except);
      } catch (const gm2calc::Error& error) {
         error_str_non_tan_beta_resummation =
            " (" + std::string(error.what()) + ")";
         // try to redo calculation w/o throwing an exception
         gm2calc::MSSMNoFV_onshell model_no_except(model);
         model_no_except.do_force_output(true);
         amu_1l_non_tan_beta_resummed =
            gm2calc::calculate_amu_1loop_non_tan_beta_resummed(model_no_except);
         amu_2l_non_tan_beta_resummed =
            gm2calc::calculate_amu_2loop_non_tan_beta_resummed(model_no_except);
      }

      const double amu_2l_tanb_approx =
         (tan_beta_cor - 1.) * amu_1l_non_tan_beta_resummed;

      const double amu_best = amu_1l + amu_2l;

      std::cout
         << "====================================================================\n"
            "   amu (1-loop + 2-loop best) = "
         << FORMAT_AMU(amu_best) << " +- "
         << FORMAT_DEL(amu_2l_uncertainty) << '\n'
         << "====================================================================\n"
            "\n"
         << error_str
         << "==============================\n"
            "   amu (1-loop) corrections\n"
            "==============================\n"
            "\n"
            "full 1L with tan(beta) resummation:\n"
            "   chi^0     " << FORMAT_AMU(gm2calc::amu1LChi0(model)) << '\n'
         << "   chi^+-    " << FORMAT_AMU(gm2calc::amu1LChipm(model)) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_1l)
                            << " (" << FORMAT_PCT(100. * amu_1l / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "full 1L without tan(beta) resummation:\n"
            "             " << FORMAT_AMU(amu_1l_non_tan_beta_resummed)
         << error_str_non_tan_beta_resummation << '\n'
         << "\n"
            "1L approximation with tan(beta) resummation:\n"
            "   W-H-nu    " << FORMAT_AMU(gm2calc::amu1LWHnu(model) * tan_beta_cor) << '\n'
         << "   W-H-muL   " << FORMAT_AMU(gm2calc::amu1LWHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muL   " << FORMAT_AMU(gm2calc::amu1LBHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muR   " << FORMAT_AMU(gm2calc::amu1LBHmuR(model) * tan_beta_cor) << '\n'
         << "   B-muL-muR " << FORMAT_AMU(gm2calc::amu1LBmuLmuR(model) * tan_beta_cor) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(gm2calc::amu1Lapprox(model)) << '\n'
         << "\n"
            "==============================\n"
            "   amu (2-loop) corrections\n"
            "==============================\n"
            "\n"
            "2L best with tan(beta) resummation:\n"
            "             " << FORMAT_AMU(amu_2l) << " (" << FORMAT_PCT(100. * amu_2l / amu_best)
         << "% of full 1L + 2L result)\n"
            "\n"
            "2L best without tan(beta) resummation:\n"
            "             " << FORMAT_AMU(amu_2l_non_tan_beta_resummed)
         << error_str_non_tan_beta_resummation << '\n'
         << "\n"
            "photonic with tan(beta) resummation:\n"
            "   chi^0     " << FORMAT_AMU(amu_2l_photonic_chi0) << '\n'
         << "   chi^+-    " << FORMAT_AMU(amu_2l_photonic_chipm) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_2l_photonic_chipm + amu_2l_photonic_chi0)
                            << " (" << FORMAT_PCT(100. * (amu_2l_photonic_chipm + amu_2l_photonic_chi0) / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "fermion/sfermion approximation with tan(beta) resummation:\n"
            "   W-H-nu    " << FORMAT_AMU(gm2calc::amu2LWHnu(model) * tan_beta_cor) << '\n'
         << "   W-H-muL   " << FORMAT_AMU(gm2calc::amu2LWHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muL   " << FORMAT_AMU(gm2calc::amu2LBHmuL(model) * tan_beta_cor) << '\n'
         << "   B-H-muR   " << FORMAT_AMU(gm2calc::amu2LBHmuR(model) * tan_beta_cor) << '\n'
         << "   B-muL-muR " << FORMAT_AMU(gm2calc::amu2LBmuLmuR(model) * tan_beta_cor) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_2l_ferm_sferm_approx)
                            << " (" << FORMAT_PCT(100. * amu_2l_ferm_sferm_approx / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "2L(a) (1L insertions into 1L SM diagram) with tan(beta) "
            "resummation:\n"
            "   sfermion  " << FORMAT_AMU(amu_2l_a_sfermion) << '\n'
         << "   cha^+-    " << FORMAT_AMU(amu_2l_a_cha) << '\n'
         << "   -------------------------------\n"
            "   sum       " << FORMAT_AMU(amu_2l_a_sfermion + amu_2l_a_cha)
                            << " (" << FORMAT_PCT(100. * (amu_2l_a_sfermion + amu_2l_a_cha) / amu_best)
                            << "% of full 1L + 2L result)\n"
            "\n"
            "tan(beta) correction:\n"
            "   amu(1L) * (1 / (1 + Delta_mu) - 1) = " << FORMAT_AMU(amu_2l_tanb_approx)
                            << " (" << FORMAT_PCT(100. * amu_2l_tanb_approx / amu_1l_non_tan_beta_resummed)
         << "%)\n";

#undef FORMAT_AMU
#undef FORMAT_PCT
   }
};

/**
 * Calculates a_mu (and potentially also the uncertainty) and writes
 * it to the SLHA i/o object.
 *
 * @param model the model (must be initialized)
 * @param options calculation options
 * @param slha_io SLHA i/o object where results are stored
 */
struct SLHA_writer {
   /// SLHA entry (block name, key, value, comment)
   using SLHA_entry = std::tuple<std::string, int, double, std::string>;

   void operator()(const gm2calc::MSSMNoFV_onshell& model,
                   const gm2calc::Config_options& options,
                   gm2calc::GM2_slha_io& slha_io)
   {
      const SLHA_entry amu_entry = [&] {
         const auto amu = calculate_amu(model, options);
         const auto amu_comment = "Delta(g-2)_muon/2";

         switch (options.output_format) {
         case gm2calc::Config_options::NMSSMTools:
            return SLHA_entry{"LOWEN", 6, amu, amu_comment};
         case gm2calc::Config_options::SPheno:
            return SLHA_entry{"SPhenoLowEnergy", 21, amu, amu_comment};
         default:
            break;
         }

         return SLHA_entry{"GM2CalcOutput", 0, amu, amu_comment};
      }();

      set_SLHA_value(slha_io, amu_entry);

      if (options.calculate_uncertainty) {
         const auto damu = calculate_uncertainty(model, options);
         const SLHA_entry damu_entry{"GM2CalcOutput", 1, damu,
                                     "uncertainty of Delta(g-2)_muon/2"};
         set_SLHA_value(slha_io, damu_entry);
      }

      if (model.get_problems().have_warning()) {
         slha_io.fill_block_entry("SPINFO", 1, "GM2Calc");
         slha_io.fill_block_entry("SPINFO", 2, GM2CALC_VERSION);
         slha_io.fill_block_entry("SPINFO", 3, model.get_problems().get_warnings());
      }

      slha_io.write_to_stream(std::cout);
   }

private:
   /**
    * Sets a entry in a given SLHA block/key.
    *
    * @param slha_io SLHA input/output object
    * @param entry tuple defining the SLHA (block name, key, value, comment)
    */
   void set_SLHA_value(gm2calc::GM2_slha_io& slha_io, const SLHA_entry& entry)
   {
      slha_io.fill_block_entry(std::get<0>(entry), std::get<1>(entry),
                               std::get<2>(entry), std::get<3>(entry));
   }
};

/**
 * Class which handles input/output.
 */
class Setup
{
public:
   Setup(const gm2calc::Config_options& options_, const Reader& reader_,
         const Writer& writer_)
      : options(options_), reader(reader_), writer(writer_)
   {
      model.do_force_output(options_.force_output);
      model.set_verbose_output(options_.verbose_output);
   }

   /// returns whether the model discovered a problem
   bool have_problem() const { return model.get_problems().have_problem(); }

   /// read from SLHA i/o object and initialize model (via reader)
   void read(const gm2calc::GM2_slha_io& slha_io)
   {
      if (!reader) {
         throw gm2calc::ESetupError("No reader set");
      }

      reader(model, slha_io);

      if (options.verbose_output) {
         VERBOSE(model);
      }

      if (model.get_problems().have_problem() ||
          model.get_problems().have_warning()) {
         std::cerr << model.get_problems() << '\n';
      }
   }

   /// Output via the writer (potentially to SLHA i/o object)
   void write(gm2calc::GM2_slha_io& slha_io) const
   {
      if (!writer) {
         throw gm2calc::ESetupError("No writer set");
      }

      writer(model, options, slha_io);
   }

private:
   gm2calc::MSSMNoFV_onshell model;
   gm2calc::Config_options options;
   Reader reader{nullptr};
   Writer writer{nullptr};
};

/**
 * Returns properly configured (but not initialized) Setup object.
 *
 * @param input_type type of input (SLHA/GM2Calc)
 * @param options configuration options
 *
 * @return Setup object
 */
Setup make_setup(
   Gm2_cmd_line_options::E_input_type input_type,
   const gm2calc::Config_options& options)
{
   const Reader reader = [&] () -> Reader {
      switch (input_type) {
      case Gm2_cmd_line_options::SLHA:
         return SLHA_reader();
      case Gm2_cmd_line_options::GM2Calc:
         return GM2Calc_reader();
      }
      throw gm2calc::ESetupError("Unknown input type");
   }();

   const Writer writer = [&] () -> Writer {
      switch (options.output_format) {
      case gm2calc::Config_options::Minimal:
         return Minimal_writer();
      case gm2calc::Config_options::Detailed:
         return Detailed_writer();
      default:
         break;
      }
      return SLHA_writer();
   }();

   return Setup(options, reader, writer);
}

} // anonymous namespace

int main(int argc, const char* argv[])
{
   Gm2_cmd_line_options options(get_cmd_line_options(argc, argv));

   if (options.input_source.empty()) {
      ERROR("No input source given!\n"
            "   Please provide an SLHA input via the option --slha-input-file=\n"
            "   or a GM2Calc input via the option --gm2calc-input-file=");
      return EXIT_FAILURE;
   }

   gm2calc::GM2_slha_io slha_io;
   gm2calc::Config_options config_options;
   int exit_code = EXIT_SUCCESS;

   try {
      set_to_default(config_options, options);
      slha_io.read_from_source(options.input_source);
      slha_io.fill(config_options);

      Setup setup = make_setup(options.input_type, config_options);
      setup.read(slha_io);
      setup.write(slha_io);

      if (setup.have_problem()) {
         exit_code = EXIT_FAILURE;
      }
   } catch (const gm2calc::Error& error) {
      print_error(error, slha_io, config_options);
      exit_code = EXIT_FAILURE;
   }

   return exit_code;
}
