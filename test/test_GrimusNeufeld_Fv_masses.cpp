
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_GrimusNeufeld_Fv_masses

#include <boost/test/unit_test.hpp>

#include "physical_input.hpp"

#include "GrimusNeufeld_two_scale_spectrum_generator.hpp"
#include "GrimusNeufeld_two_scale_model.hpp"
#include "GrimusNeufeld_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_GrimusNeufeld_Fv_masses )
{

  char const * const slha_input = R"(
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # solver (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   1                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   0                    # beta-functions loop order
    7   0                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   1                    # calculate observables (a_muon, ...)
   16   1                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   23   1                    # calculate BSM pole masses
   24   124111421            # individual threshold correction loop orders
   31   0                    # loop library (0 = softsusy)
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279160000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.776990000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.385               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583715e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR
   1   0.258766              # Lambda1IN
   2   0.051388              # Lambda2IN
   3   0.008974              # Lambda3IN
   4   0.430832              # Lambda4IN
   5   0.000030              # Lambda5IN
   9   14197.502739          # M222IN
Block EXTPAR
   0   110                   # Qin
Block Theta122313IN
   1   0.59                  # Theta122313(1)
   2   0.84                  # Theta122313(2)
   3   0.15                  # Theta122313(3)
Block deltaCP
   4.5                       # deltaCP
Block deltaM2
   0                         # deltaM2
Block Inverted
   0                         # 0 - NO, 1 - IO
Block MnuIN
   2   8.6948e-12            # MnuIN(2)
   3   5.1240e-11            # MnuIN(3)
   4   1.00e-02              # MnuIN(4)
Block ROPhiIN
   1 0.7                     # r
   2 1.0                     # omega22
   3 0.0                     # phiR
)";

   std::stringstream istr(slha_input);

   GrimusNeufeld_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   GrimusNeufeld_input_parameters input;
   Spectrum_generator_settings settings;

   // extract the input parameters from spectrum string
   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(physical_input);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);
   GrimusNeufeld_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   auto m = std::get<0>(spectrum_generator.get_models_slha());

   // MFv(0) == 0 both at tree and 1-loop level
   BOOST_CHECK_SMALL(m.get_physical().MFv(0), 1e-16);
   BOOST_CHECK_SMALL(m.get_MFv(0), 1e-16);

   // MFv(1) == 0 at the tree-level but massive at 1-loop
   BOOST_CHECK_SMALL(m.get_MFv(1), 1e-16);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MFv(1), 8.6949585400840875e-12, 4e-16);

   BOOST_CHECK_CLOSE_FRACTION(m.get_MFv(2), 2.3901681143242786e-11, 1e-16);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MFv(2), 5.1239706327974192e-11, 3e-16);

   // MFv(3) == MNuIn(4) (exactly) at the tree level
   BOOST_CHECK_CLOSE_FRACTION(m.get_MFv(3), 0.01, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MFv(3), 0.010000000879295179, 1e-16);
}
