
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_THDMII_FlexibleDecay

#include <boost/test/unit_test.hpp>

#include "THDMII_two_scale_spectrum_generator.hpp"
#include "THDMII_two_scale_model.hpp"
#include "decays/THDMII_decays.hpp"
#include "THDMII_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_THDMII_FlexibleDecay )
{

  char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.4.1
     5   THDMII
     9   4.14.3
Block MODSEL                 # Select model
#   12    1000                # parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   1                    # calculate SM pole masses
    4   0                    # pole mass loop order
    5   0                    # EWSB loop order
    6   0                    # beta-functions loop order
    7   0                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   24   0                    # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   -1                   # 0(Softsusy),1(Collier),2(Looptools),3(fflite)
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733000000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.404               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR
     1     1.30331451E+00   # Lambda1IN
     2     1.20621999E-01   # Lambda2IN
     3     6.02834555E+00   # Lambda3IN
     4    -3.22346966E+00   # Lambda4IN
     5    -1.49148565E+00   # Lambda5IN
     6     0.00000000E+00   # Lambda6IN
     7     0.00000000E+00   # Lambda7IN
     8     1.58000000E+04   # M122IN
    10     1.00000000E+01   # TanBeta
Block EXTPAR
     0     1.25000000E+02   # Qin
Block gauge Q= 1.25000000E+02
     1     3.58002663E-01   # g1 * 0.7745966692414834
     2     6.48438444E-01   # g2
     3     1.21565011E+00   # g3
Block Yu Q= 1.25000000E+02
  1  1     8.02484356E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.66968592E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     1.00034899E+00   # Yu(3,3)
Block Yd Q= 1.25000000E+02
  1  1     1.59566777E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     3.49367259E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     1.66494343E-01   # Yd(3,3)
Block Ye Q= 1.25000000E+02
  1  1     2.94966667E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     6.09897463E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     1.02574734E-01   # Ye(3,3)
Block HMIX Q= 1.25000000E+02
    31     1.30331451E+00   # Lambda1
    32     1.20621999E-01   # Lambda2
    33     6.02834555E+00   # Lambda3
    34    -3.22346966E+00   # Lambda4
    35    -1.49148565E+00   # Lambda5
    36     0.00000000E+00   # Lambda6
    37     0.00000000E+00   # Lambda7
    20     1.17800313E+05   # M112
    21    -6.05437397E+03   # M222
    22     1.58000000E+04   # M122
   102     2.44997710E+01   # v1
   103     2.44997710E+02   # v2
Block MASS
        24     8.04040000E+01   # VWm
        37     5.50000000E+02   # Hm(2)
        25     1.25000000E+02   # hh(1)
        35     4.00000000E+02   # hh(2)
        36     5.00000000E+02   # Ah(2)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
         1     2.76432753E-03   # Fd(1)
         3     6.05242238E-02   # Fd(2)
         5     2.88434037E+00   # Fd(3)
         2     1.39022022E-03   # Fu(1)
         4     6.35734708E-01   # Fu(2)
         6     1.73300000E+02   # Fu(3)
        11     5.10998902E-04   # Fe(1)
        13     1.05658357E-01   # Fe(2)
        15     1.77700000E+00   # Fe(3)
        22     0.00000000E+00   # VP
        23     9.11876000E+01   # VZ
Block PSEUDOSCALARMIX
  1  1     9.95037190E-02   # ZA(1,1)
  1  2     9.95037190E-01   # ZA(1,2)
  2  1     9.95037190E-01   # ZA(2,1)
  2  2    -9.95037190E-02   # ZA(2,2)
Block SCALARMIX
  1  1     5.49159256E-02   # ZH(1,1)
  1  2     9.98490982E-01   # ZH(1,2)
  2  1     9.98490982E-01   # ZH(2,1)
  2  2    -5.49159256E-02   # ZH(2,2)
Block CHARGEMIX
  1  1     9.95037190E-02   # ZP(1,1)
  1  2     9.95037190E-01   # ZP(1,2)
  2  1     9.95037190E-01   # ZP(2,1)
  2  2    -9.95037190E-02   # ZP(2,2)
Block UULMIX
  1  1     1.00000000E+00   # Re(Vu(1,1))
  1  2     0.00000000E+00   # Re(Vu(1,2))
  1  3     0.00000000E+00   # Re(Vu(1,3))
  2  1     0.00000000E+00   # Re(Vu(2,1))
  2  2     1.00000000E+00   # Re(Vu(2,2))
  2  3     0.00000000E+00   # Re(Vu(2,3))
  3  1     0.00000000E+00   # Re(Vu(3,1))
  3  2     0.00000000E+00   # Re(Vu(3,2))
  3  3     1.00000000E+00   # Re(Vu(3,3))
Block UDLMIX
  1  1     1.00000000E+00   # Re(Vd(1,1))
  1  2     0.00000000E+00   # Re(Vd(1,2))
  1  3     0.00000000E+00   # Re(Vd(1,3))
  2  1     0.00000000E+00   # Re(Vd(2,1))
  2  2     1.00000000E+00   # Re(Vd(2,2))
  2  3     0.00000000E+00   # Re(Vd(2,3))
  3  1     0.00000000E+00   # Re(Vd(3,1))
  3  2     0.00000000E+00   # Re(Vd(3,2))
  3  3     1.00000000E+00   # Re(Vd(3,3))
Block UURMIX
  1  1     1.00000000E+00   # Re(Uu(1,1))
  1  2     0.00000000E+00   # Re(Uu(1,2))
  1  3     0.00000000E+00   # Re(Uu(1,3))
  2  1     0.00000000E+00   # Re(Uu(2,1))
  2  2     1.00000000E+00   # Re(Uu(2,2))
  2  3     0.00000000E+00   # Re(Uu(2,3))
  3  1     0.00000000E+00   # Re(Uu(3,1))
  3  2     0.00000000E+00   # Re(Uu(3,2))
  3  3     1.00000000E+00   # Re(Uu(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(Ud(1,1))
  1  2     0.00000000E+00   # Re(Ud(1,2))
  1  3     0.00000000E+00   # Re(Ud(1,3))
  2  1     0.00000000E+00   # Re(Ud(2,1))
  2  2     1.00000000E+00   # Re(Ud(2,2))
  2  3     0.00000000E+00   # Re(Ud(2,3))
  3  1     0.00000000E+00   # Re(Ud(3,1))
  3  2     0.00000000E+00   # Re(Ud(3,2))
  3  3     1.00000000E+00   # Re(Ud(3,3))
Block UELMIX
  1  1     1.00000000E+00   # Re(Ve(1,1))
  1  2     0.00000000E+00   # Re(Ve(1,2))
  1  3     0.00000000E+00   # Re(Ve(1,3))
  2  1     0.00000000E+00   # Re(Ve(2,1))
  2  2     1.00000000E+00   # Re(Ve(2,2))
  2  3     0.00000000E+00   # Re(Ve(2,3))
  3  1     0.00000000E+00   # Re(Ve(3,1))
  3  2     0.00000000E+00   # Re(Ve(3,2))
  3  3     1.00000000E+00   # Re(Ve(3,3))
Block UERMIX
  1  1     1.00000000E+00   # Re(Ue(1,1))
  1  2     0.00000000E+00   # Re(Ue(1,2))
  1  3     0.00000000E+00   # Re(Ue(1,3))
  2  1     0.00000000E+00   # Re(Ue(2,1))
  2  2     1.00000000E+00   # Re(Ue(2,2))
  2  3     0.00000000E+00   # Re(Ue(2,3))
  3  1     0.00000000E+00   # Re(Ue(3,1))
  3  2     0.00000000E+00   # Re(Ue(3,2))
  3  3     1.00000000E+00   # Re(Ue(3,3))
)";

   std::stringstream istr(slha_input);

   THDMII_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   THDMII_input_parameters input;
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

   THDMII_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   THDMII_slha m = std::get<0>(spectrum_generator.get_models_slha());

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   THDMII_decays decays_with_HO(m, qedqcd, physical_input, SM_higher_order_corrections::enable);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.00074272223389794053, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              0.00011992760375249642, 2e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              7.9645091090513334e-05, 1e-15);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
   //                            0.0007420072130724369, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
                              0.00081324389707053922, 1e-4);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
   //                           8.3396512232490906e-05, 2e-14);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
                              0.00010470911466556374, 1e-4);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VGVG(&m, 0), 0.00035337029100761187, 4e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 1.0119606786130193e-05, 3e-12);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVZ(&m, 0), 6.4088569360381209e-06, 8e-12);

   // Ah -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00078773836543896236, 3e-14);
   // Ah -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VPVP(&m, 1), 4.0996746335596648e-06, 3e-12);
   // Ah -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VPVZ(&m, 1), 5.8406939485474199e-07, 8e-12);

   // Ah -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFdFd(&m, 0, 2, 2),
                              0.0018671335954316421, 2e-15);
   // Ah -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFuFu(&m, 0, 1, 1),
                              9.3482696650389066e-05, 2e-16);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   THDMII_decays decays_without_HO(m, qedqcd, physical_input, SM_higher_order_corrections::disable);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.00062168170474045936, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              0.00010014749564849829, 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              7.8811709644693243e-05, 1e-15);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VZVZ(&m, 0),
   //                           8.3396512232490906e-05, 2e-14);
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VZVZ(&m, 0),
                              0.00010470911466556374, 1e-3);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VGVG(&m, 0), 0.00020488527576778421, 4e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVP(&m, 0), 9.9357276197179806e-06, 3e-12);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVZ(&m, 0), 6.3818951172675953e-06, 8e-12);

   // Ah -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00078773836543896236, 4e-13);
   // Ah -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Ah_to_VPVP(&m, 1), 2.9199349207376965e-06, 3e-12);
   // Ah -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Ah_to_VPVZ(&m, 1), 5.8406939485474199e-07, 8e-12);
}
