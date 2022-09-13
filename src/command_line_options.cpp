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

#include "command_line_options.hpp"
#include "array_view.hpp"
#include "build_info.hpp"
#include "logger.hpp"

#include <cstdio>

namespace flexiblesusy {

Command_line_options::Command_line_options(int argc, char* argv[])
   : Command_line_options(make_dynamic_array_view(&argv[0], argc))
{
}

Command_line_options::Command_line_options(const Dynamic_array_view<char*>& args)
{
   parse(args);
}

/**
 * Parse string of program options and store the given option values
 * in the member variables of this class.
 *
 * @param args command line arguments
 */
void Command_line_options::parse(const Dynamic_array_view<char*>& args)
{
   reset();

   if (args.empty())
      return;

   program = args[0];

   for (int i = 1; i < args.size(); ++i) {
      const std::string option(args[i]);
      if (starts_with(option,"--slha-input-file=")) {
         slha_input_file = option.substr(18);
      } else if (starts_with(option,"--slha-output-file=")) {
         slha_output_file = option.substr(19);
      } else if (starts_with(option,"--spectrum-output-file=")) {
         spectrum_file = option.substr(23);
      } else if (starts_with(option,"--database-output-file=")) {
         database_output_file = option.substr(23);
      } else if (starts_with(option,"--rgflow-output-file=")) {
         rgflow_file = option.substr(21);
      } else if (option == "--help" || option == "-h") {
         print_usage(std::cout);
         do_exit = true;
      } else if (option == "--build-info") {
         print_build_info(std::cout);
         do_exit = true;
      } else if (option == "--model-info") {
         do_print_model_info = true;
         do_exit = true;
      } else if (option == "--version" || option == "-v") {
         print_version(std::cout);
         do_exit = true;
      } else {
         ERROR("Unrecognized command line option: " << option);
         do_exit = true;
         exit_status = EXIT_FAILURE;
      }
   }
}

void Command_line_options::print_build_info(std::ostream& ostr) const
{
   flexiblesusy::print_all_info(ostr);
}

void Command_line_options::print_version(std::ostream& ostr) const
{
   ostr << "FlexibleSUSY ";
   print_flexiblesusy_version(ostr);
   ostr << std::endl;
}

void Command_line_options::print_usage(std::ostream& ostr) const
{
   ostr << "Usage: " << program << " [options]\n"
           "Options:\n"
           "  --slha-input-file=<filename>      SLHA input file (default: )\n"
           "                                    If <filename> is - then the"
                                                " SLHA input\n"
           "                                    is read from stdin.\n"
           "  --slha-output-file=<filename>     SLHA output file (default: -)\n"
           "                                    If <filename> is the empty string, then\n"
           "                                    no SLHA output is written.\n"
           "                                    If <filename> is - then the output is\n"
           "                                    printed to stdout.\n"
           "  --spectrum-output-file=<filename> file to write spectrum to\n"
           "  --database-output-file=<filename> SQLite database file to write\n"
           "                                    parameter point to\n"
           "  --rgflow-output-file=<filename>   file to write rgflow to\n"
           "  --build-info                      print build information\n"
           "  --model-info                      print model information\n"
           "  --help,-h                         print this help message\n"
           "  --version,-v                      print program version"
        << std::endl;
}

/**
 * Resets all command line options to their initial values.
 */
void Command_line_options::reset()
{
   *this = Command_line_options();
}

/**
 * Returns true if the string str starts with prefix, false otherwise.
 *
 * @param str string to search in
 * @param prefix string to search for
 *
 * @return true if the string str starts with prefix, false otherwise
 */
bool Command_line_options::starts_with(const std::string& str,
                                       const std::string& prefix)
{
   return !str.compare(0, prefix.size(), prefix);
}

/**
 * Extracts the parameter value from a command line option string of
 * the form --m0=125 .
 *
 * @param str full option string, including the parameter value (--m0=125)
 * @param prefix option string, without the parameter value (--m0=)
 * @param[out] parameter output parameter value (of type double)
 *
 * @return true, if str starts with prefix, false otherwise
 */
bool Command_line_options::get_parameter_value(const std::string& str,
                                               const std::string& prefix,
                                               double& parameter)
{
   if (starts_with(str, prefix)) {
      parameter = std::stod(str.substr(prefix.length()));
      return true;
   }
   return false;
}

/**
 * Extracts the parameter value from a command line option string of
 * the form --m0=125 .
 *
 * @param str full option string, including the parameter value (--m0=125)
 * @param prefix option string, without the parameter value (--m0=)
 * @param[out] parameter output parameter value (of type int)
 *
 * @return true, if str starts with prefix, false otherwise
 */
bool Command_line_options::get_parameter_value(const std::string& str,
                                               const std::string& prefix,
                                               int& parameter)
{
   if (starts_with(str, prefix)) {
      parameter = std::stoi(str.substr(prefix.length()));
      return true;
   }
   return false;
}

} // namespace flexiblesusy
