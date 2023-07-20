
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_BLSMlightZp_amm

#include <boost/test/unit_test.hpp>

#include "BLSMlightZp_amm.hpp"
#include "BLSMlightZp_two_scale_spectrum_generator.hpp"
#include "BLSMlightZp_two_scale_model.hpp"
#include "BLSMlightZp_slha_io.hpp"
#include "cxx_qft/BLSMlightZp_qft.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_BLSMlightZp_amm )
{

  char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.7.1
     5   BLSMlightZp
     9   4.15.1
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   1                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   0                    # loop library (0 = softsusy)
   32   2                    # loop level to calculate AMM
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
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR
     1    -5.73224594E-01   # Lambda1INPUT
     2    -3.37519068E-01   # Lambda2INPUT
     3     7.19718182E-01   # Lambda3INPUT
    10     5.13036082E-01   # g1pINPUT
    11     5.07775177E-01   # g1p1INPUT
    12     8.44949257E-01   # g11pINPUT
    20     5.81797727E-02   # vXinput
Block YvIN
  1  1     0.00000000E+00   # YvInput(1,1)
  1  2     0.00000000E+00   # YvInput(1,2)
  1  3     0.00000000E+00   # YvInput(1,3)
  2  1     0.00000000E+00   # YvInput(2,1)
  2  2     0.00000000E+00   # YvInput(2,2)
  2  3     0.00000000E+00   # YvInput(2,3)
  3  1     0.00000000E+00   # YvInput(3,1)
  3  2     0.00000000E+00   # YvInput(3,2)
  3  3     0.00000000E+00   # YvInput(3,3)
Block YXIN
  1  1     3.20000000E-01   # YxInput(1,1)
  1  2     0.00000000E+00   # YxInput(1,2)
  1  3     0.00000000E+00   # YxInput(1,3)
  2  1     0.00000000E+00   # YxInput(2,1)
  2  2     3.20000000E-01   # YxInput(2,2)
  2  3     0.00000000E+00   # YxInput(2,3)
  3  1     0.00000000E+00   # YxInput(3,1)
  3  2     0.00000000E+00   # YxInput(3,2)
  3  3     3.20000000E-01   # YxInput(3,3)
Block gauge Q= 9.11876000E+01
     1     3.81222554E-01   # g1 * 0.7745966692414834
    11     8.44949257E-01   # g11p
    10     5.07775177E-01   # g1p1
     4     6.28338311E-01   # g1p * 1.224744871391589
     2     5.49315571E-01   # g2
     3     1.21088651E+00   # g3
Block Yu Q= 9.11876000E+01
  1  1     6.72312149E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.07442059E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.29487430E-01   # Yu(3,3)
Block Yd Q= 9.11876000E+01
  1  1     1.33683207E-05   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     2.92695865E-04   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     1.39519332E-02   # Yd(3,3)
Block Ye Q= 9.11876000E+01
  1  1     2.06688605E-06   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     4.32568503E-04   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     7.58885263E-03   # Ye(3,3)
Block HMIX Q= 9.11876000E+01
     3     2.92433848E+02   # vH
Block BL Q= 9.11876000E+01
    43     5.81797727E-02   # vX
     1    -5.73224594E-01   # L1
     2    -3.37519068E-01   # L2
     3     7.19718182E-01   # L3
    10     3.00105695E+04   # MuP
    11    -5.04275482E+04   # mu2
Block YX Q= 9.11876000E+01
  1  1     3.20000000E-01   # Yx(1,1)
  1  2     0.00000000E+00   # Yx(1,2)
  1  3     0.00000000E+00   # Yx(1,3)
  2  1     0.00000000E+00   # Yx(2,1)
  2  2     3.20000000E-01   # Yx(2,2)
  2  3     0.00000000E+00   # Yx(2,3)
  3  1     0.00000000E+00   # Yx(3,1)
  3  2     0.00000000E+00   # Yx(3,2)
  3  3     3.20000000E-01   # Yx(3,3)
Block Yv Q= 9.11876000E+01
  1  1     0.00000000E+00   # Yv(1,1)
  1  2     0.00000000E+00   # Yv(1,2)
  1  3     0.00000000E+00   # Yv(1,3)
  2  1     0.00000000E+00   # Yv(2,1)
  2  2     0.00000000E+00   # Yv(2,2)
  2  3     0.00000000E+00   # Yv(2,3)
  3  1     0.00000000E+00   # Yv(3,1)
  3  2     0.00000000E+00   # Yv(3,2)
  3  3     0.00000000E+00   # Yv(3,3)
Block MASS
        24     8.03850000E+01   # VWm
        25     3.19306239E-02   # hh(1)
        35     3.24775536E+02   # hh(2)
        22     0.00000000E+00   # VP
        31     2.77814658E-02   # VZp
        23     1.53867421E+02   # VZ
        21     0.00000000E+00   # VG
         1     5.04825539E-03   # Fd(1)
         3     9.66182962E-02   # Fd(2)
         5     3.76960800E+00   # Fd(3)
         2     2.69515856E-03   # Fu(1)
         4     9.29763888E-01   # Fu(2)
         6     1.73230885E+02   # Fu(3)
        11     5.31475215E-04   # Fe(1)
        13     1.08393319E-01   # Fe(2)
        15     1.78521611E+00   # Fe(3)
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
   8810012     2.35650723E-02   # Fv(4)
   8810014     2.35650723E-02   # Fv(5)
   8810016     2.35650723E-02   # Fv(6)
Block PSEUDOSCALARMIX
  1  1     0.00000000E+00   # ZA(1,1)
  1  2     0.00000000E+00   # ZA(1,2)
  2  1     0.00000000E+00   # ZA(2,1)
  2  2     0.00000000E+00   # ZA(2,2)
Block SCALARMIX
  1  1     1.05832257E-04   # ZH(1,1)
  1  2     9.99999994E-01   # ZH(1,2)
  2  1     9.99999994E-01   # ZH(2,1)
  2  2    -1.05832257E-04   # ZH(2,2)
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
Block UVMIX
  1  1     1.00000000E+00   # Re(ZM(1,1))
  1  2     0.00000000E+00   # Re(ZM(1,2))
  1  3     0.00000000E+00   # Re(ZM(1,3))
  1  4     0.00000000E+00   # Re(ZM(1,4))
  1  5     0.00000000E+00   # Re(ZM(1,5))
  1  6     0.00000000E+00   # Re(ZM(1,6))
  2  1     0.00000000E+00   # Re(ZM(2,1))
  2  2     1.00000000E+00   # Re(ZM(2,2))
  2  3     0.00000000E+00   # Re(ZM(2,3))
  2  4     0.00000000E+00   # Re(ZM(2,4))
  2  5     0.00000000E+00   # Re(ZM(2,5))
  2  6     0.00000000E+00   # Re(ZM(2,6))
  3  1     0.00000000E+00   # Re(ZM(3,1))
  3  2     0.00000000E+00   # Re(ZM(3,2))
  3  3     1.00000000E+00   # Re(ZM(3,3))
  3  4     0.00000000E+00   # Re(ZM(3,4))
  3  5     0.00000000E+00   # Re(ZM(3,5))
  3  6     0.00000000E+00   # Re(ZM(3,6))
  4  1     0.00000000E+00   # Re(ZM(4,1))
  4  2     0.00000000E+00   # Re(ZM(4,2))
  4  3     0.00000000E+00   # Re(ZM(4,3))
  4  4     1.00000000E+00   # Re(ZM(4,4))
  4  5     0.00000000E+00   # Re(ZM(4,5))
  4  6     0.00000000E+00   # Re(ZM(4,6))
  5  1     0.00000000E+00   # Re(ZM(5,1))
  5  2     0.00000000E+00   # Re(ZM(5,2))
  5  3     0.00000000E+00   # Re(ZM(5,3))
  5  4     0.00000000E+00   # Re(ZM(5,4))
  5  5     1.00000000E+00   # Re(ZM(5,5))
  5  6     0.00000000E+00   # Re(ZM(5,6))
  6  1     0.00000000E+00   # Re(ZM(6,1))
  6  2     0.00000000E+00   # Re(ZM(6,2))
  6  3     0.00000000E+00   # Re(ZM(6,3))
  6  4     0.00000000E+00   # Re(ZM(6,4))
  6  5     0.00000000E+00   # Re(ZM(6,5))
  6  6     1.00000000E+00   # Re(ZM(6,6))
)";

   std::stringstream istr(slha_input);

   BLSMlightZp_slha_io slha_io;
   slha_io.read_from_stream(istr);

   // extract the input parameters
   softsusy::QedQcd qedqcd;
   BLSMlightZp_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(qedqcd);
      slha_io.fill(input);
      slha_io.fill(settings);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   BLSMlightZp_slha m(input);
   slha_io.fill(m);
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   using BLSMlightZp_cxx_diagrams::fields::Fe;

   // 1L
   settings.set(Spectrum_generator_settings::calculate_amm, 1.0);

   auto ae = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, 5.1118603233954821e-07, 1e-7);
   auto amu = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, 0.0026497259684861866, 1e-7);
   auto atau = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, 0.0035964457880142101, 1e-7);

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 1.5);

   ae = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, 5.1118603233954821e-07, 1e-7);
   amu = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, 0.0026497259684861866, 1e-7);
   atau = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, 0.0035964457880142101, 1e-7);

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 2.0);
   ae = BLSMlightZp_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   std::cout << ae << std::endl;
}
