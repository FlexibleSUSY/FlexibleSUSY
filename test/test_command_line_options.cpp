#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_command_line_options

#include <boost/test/unit_test.hpp>

#include "command_line_options.hpp"

#include <array>
#include <cstdlib>
#include <span>
#include <string>

using namespace flexiblesusy;

namespace
{

/// returns the arguments as array of non-const char*, as in main()
template <std::size_t N>
std::array<char*, N> make_args(std::array<std::string, N>& strs)
{
   std::array<char*, N> args{};
   for (std::size_t i = 0; i < N; ++i) {
      args[i] = strs[i].data();
   }
   return args;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_no_arguments)
{
   const Command_line_options options(std::span<char* const>{});

   BOOST_CHECK(!options.must_exit());
   BOOST_CHECK_EQUAL(options.status(), EXIT_SUCCESS);
   BOOST_CHECK_EQUAL(options.get_program_name(), "");
   BOOST_CHECK_EQUAL(options.get_slha_output_file(), "-");
}

BOOST_AUTO_TEST_CASE(test_program_name_only)
{
   std::array<std::string, 1> strs{"run.x"};
   auto args = make_args(strs);
   const Command_line_options options(static_cast<int>(args.size()), args.data());

   BOOST_CHECK(!options.must_exit());
   BOOST_CHECK_EQUAL(options.status(), EXIT_SUCCESS);
   BOOST_CHECK_EQUAL(options.get_program_name(), "run.x");
   BOOST_CHECK_EQUAL(options.get_slha_input_file(), "");
}

BOOST_AUTO_TEST_CASE(test_file_options)
{
   std::array<std::string, 5> strs{
      "run.x", "--slha-input-file=in.slha", "--slha-output-file=out.slha",
      "--spectrum-output-file=spec.dat", "--rgflow-output-file=rgflow.dat"};
   auto args = make_args(strs);
   const Command_line_options options(static_cast<int>(args.size()), args.data());

   BOOST_CHECK(!options.must_exit());
   BOOST_CHECK_EQUAL(options.status(), EXIT_SUCCESS);
   BOOST_CHECK_EQUAL(options.get_program_name(), "run.x");
   BOOST_CHECK_EQUAL(options.get_slha_input_file(), "in.slha");
   BOOST_CHECK_EQUAL(options.get_slha_output_file(), "out.slha");
   BOOST_CHECK_EQUAL(options.get_spectrum_file(), "spec.dat");
   BOOST_CHECK_EQUAL(options.get_rgflow_file(), "rgflow.dat");
}

BOOST_AUTO_TEST_CASE(test_unknown_option)
{
   std::array<std::string, 3> strs{"run.x", "--slha-input-file=in.slha",
                                   "--unknown"};
   auto args = make_args(strs);
   const Command_line_options options{std::span<char* const>(args)};

   BOOST_CHECK(options.must_exit());
   BOOST_CHECK_EQUAL(options.status(), EXIT_FAILURE);
   BOOST_CHECK_EQUAL(options.get_slha_input_file(), "in.slha");
}

BOOST_AUTO_TEST_CASE(test_parse_resets_options)
{
   std::array<std::string, 2> strs1{"run.x", "--unknown"};
   auto args1 = make_args(strs1);
   Command_line_options options{std::span<char* const>(args1)};

   BOOST_CHECK(options.must_exit());

   std::array<std::string, 2> strs2{"run2.x", "--model-info"};
   auto args2 = make_args(strs2);
   options.parse(args2);

   BOOST_CHECK(options.must_exit());
   BOOST_CHECK(options.must_print_model_info());
   BOOST_CHECK_EQUAL(options.status(), EXIT_SUCCESS);
   BOOST_CHECK_EQUAL(options.get_program_name(), "run2.x");
}

BOOST_AUTO_TEST_CASE(test_get_parameter_value)
{
   double m0 = 0.;
   int sign_mu = 0;

   BOOST_CHECK(
      Command_line_options::get_parameter_value("--m0=125.5", "--m0=", m0));
   BOOST_CHECK_EQUAL(m0, 125.5);
   BOOST_CHECK(
      !Command_line_options::get_parameter_value("--m12=1", "--m0=", m0));
   BOOST_CHECK_EQUAL(m0, 125.5);
   BOOST_CHECK(Command_line_options::get_parameter_value("--SignMu=-1",
                                                         "--SignMu=", sign_mu));
   BOOST_CHECK_EQUAL(sign_mu, -1);
}
