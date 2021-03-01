
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMCPV_FlexibleDecay

#include <boost/test/unit_test.hpp>

#include "MSSMCPV_two_scale_spectrum_generator.hpp"
#include "MSSMCPV_two_scale_model.hpp"
#include "decays/MSSMCPV_decays.hpp"
#include "MSSMCPV_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MSSMCPV_FlexibleDecay )
{

  char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.6.0
     5   MSSMCPV
     9   4.14.3
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
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
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
     3     5.00000000E+00   # TanBeta
Block EXTPAR
     0     2.00000000E+03   # MSUSY
     1     2.00000000E+03   # ReM1Input
     2     2.00000000E+03   # ReM2Input
     3     2.00000000E+03   # ReM3Input
    23     2.00000000E+03   # ReMuInput
    24     4.00000000E+06   # mA2Input
   100     1.00000000E-02   # etaInput
Block IMEXTPAR
     1     1.00000000E+02   # ImM1Input
     2     1.00000000E+02   # ImM2Input
     3     1.00000000E+02   # ImM3Input
    23     1.00000000E+02   # ImMuInput
Block IMADIN
  1  1     0.00000000E+00   # ImAdInput(1,1)
  1  2     0.00000000E+00   # ImAdInput(1,2)
  1  3     0.00000000E+00   # ImAdInput(1,3)
  2  1     0.00000000E+00   # ImAdInput(2,1)
  2  2     0.00000000E+00   # ImAdInput(2,2)
  2  3     0.00000000E+00   # ImAdInput(2,3)
  3  1     0.00000000E+00   # ImAdInput(3,1)
  3  2     0.00000000E+00   # ImAdInput(3,2)
  3  3     0.00000000E+00   # ImAdInput(3,3)
Block IMAEIN
  1  1     0.00000000E+00   # ImAeInput(1,1)
  1  2     0.00000000E+00   # ImAeInput(1,2)
  1  3     0.00000000E+00   # ImAeInput(1,3)
  2  1     0.00000000E+00   # ImAeInput(2,1)
  2  2     0.00000000E+00   # ImAeInput(2,2)
  2  3     0.00000000E+00   # ImAeInput(2,3)
  3  1     0.00000000E+00   # ImAeInput(3,1)
  3  2     0.00000000E+00   # ImAeInput(3,2)
  3  3     0.00000000E+00   # ImAeInput(3,3)
Block IMAUIN
  1  1     0.00000000E+00   # ImAuInput(1,1)
  1  2     0.00000000E+00   # ImAuInput(1,2)
  1  3     0.00000000E+00   # ImAuInput(1,3)
  2  1     0.00000000E+00   # ImAuInput(2,1)
  2  2     0.00000000E+00   # ImAuInput(2,2)
  2  3     0.00000000E+00   # ImAuInput(2,3)
  3  1     0.00000000E+00   # ImAuInput(3,1)
  3  2     0.00000000E+00   # ImAuInput(3,2)
  3  3     0.00000000E+00   # ImAuInput(3,3)
Block MSD2IN
  1  1     4.00000000E+06   # md2Input(1,1)
  1  2     0.00000000E+00   # md2Input(1,2)
  1  3     0.00000000E+00   # md2Input(1,3)
  2  1     0.00000000E+00   # md2Input(2,1)
  2  2     4.00000000E+06   # md2Input(2,2)
  2  3     0.00000000E+00   # md2Input(2,3)
  3  1     0.00000000E+00   # md2Input(3,1)
  3  2     0.00000000E+00   # md2Input(3,2)
  3  3     4.00000000E+06   # md2Input(3,3)
Block MSE2IN
  1  1     4.00000000E+06   # me2Input(1,1)
  1  2     0.00000000E+00   # me2Input(1,2)
  1  3     0.00000000E+00   # me2Input(1,3)
  2  1     0.00000000E+00   # me2Input(2,1)
  2  2     4.00000000E+06   # me2Input(2,2)
  2  3     0.00000000E+00   # me2Input(2,3)
  3  1     0.00000000E+00   # me2Input(3,1)
  3  2     0.00000000E+00   # me2Input(3,2)
  3  3     4.00000000E+06   # me2Input(3,3)
Block MSL2IN
  1  1     4.00000000E+06   # ml2Input(1,1)
  1  2     0.00000000E+00   # ml2Input(1,2)
  1  3     0.00000000E+00   # ml2Input(1,3)
  2  1     0.00000000E+00   # ml2Input(2,1)
  2  2     4.00000000E+06   # ml2Input(2,2)
  2  3     0.00000000E+00   # ml2Input(2,3)
  3  1     0.00000000E+00   # ml2Input(3,1)
  3  2     0.00000000E+00   # ml2Input(3,2)
  3  3     4.00000000E+06   # ml2Input(3,3)
Block MSQ2IN
  1  1     4.00000000E+06   # mq2Input(1,1)
  1  2     0.00000000E+00   # mq2Input(1,2)
  1  3     0.00000000E+00   # mq2Input(1,3)
  2  1     0.00000000E+00   # mq2Input(2,1)
  2  2     4.00000000E+06   # mq2Input(2,2)
  2  3     0.00000000E+00   # mq2Input(2,3)
  3  1     0.00000000E+00   # mq2Input(3,1)
  3  2     0.00000000E+00   # mq2Input(3,2)
  3  3     4.00000000E+06   # mq2Input(3,3)
Block MSU2IN
  1  1     4.00000000E+06   # mu2Input(1,1)
  1  2     0.00000000E+00   # mu2Input(1,2)
  1  3     0.00000000E+00   # mu2Input(1,3)
  2  1     0.00000000E+00   # mu2Input(2,1)
  2  2     4.00000000E+06   # mu2Input(2,2)
  2  3     0.00000000E+00   # mu2Input(2,3)
  3  1     0.00000000E+00   # mu2Input(3,1)
  3  2     0.00000000E+00   # mu2Input(3,2)
  3  3     4.00000000E+06   # mu2Input(3,3)
Block ADIN
  1  1     1.00000000E+04   # ReAdInput(1,1)
  1  2     0.00000000E+00   # ReAdInput(1,2)
  1  3     0.00000000E+00   # ReAdInput(1,3)
  2  1     0.00000000E+00   # ReAdInput(2,1)
  2  2     1.00000000E+04   # ReAdInput(2,2)
  2  3     0.00000000E+00   # ReAdInput(2,3)
  3  1     0.00000000E+00   # ReAdInput(3,1)
  3  2     0.00000000E+00   # ReAdInput(3,2)
  3  3     1.00000000E+04   # ReAdInput(3,3)
Block AEIN
  1  1     1.00000000E+04   # ReAeInput(1,1)
  1  2     0.00000000E+00   # ReAeInput(1,2)
  1  3     0.00000000E+00   # ReAeInput(1,3)
  2  1     0.00000000E+00   # ReAeInput(2,1)
  2  2     1.00000000E+04   # ReAeInput(2,2)
  2  3     0.00000000E+00   # ReAeInput(2,3)
  3  1     0.00000000E+00   # ReAeInput(3,1)
  3  2     0.00000000E+00   # ReAeInput(3,2)
  3  3     1.00000000E+04   # ReAeInput(3,3)
Block AUIN
  1  1     4.00000000E+02   # ReAuInput(1,1)
  1  2     0.00000000E+00   # ReAuInput(1,2)
  1  3     0.00000000E+00   # ReAuInput(1,3)
  2  1     0.00000000E+00   # ReAuInput(2,1)
  2  2     4.00000000E+02   # ReAuInput(2,2)
  2  3     0.00000000E+00   # ReAuInput(2,3)
  3  1     0.00000000E+00   # ReAuInput(3,1)
  3  2     0.00000000E+00   # ReAuInput(3,2)
  3  3     4.00000000E+02   # ReAuInput(3,3)
Block gauge Q= 2.00000000E+03
     1     3.63590842E-01   # g1 * 0.7745966692414834
     2     6.36432161E-01   # g2
     3     1.02558883E+00   # g3
Block Yu Q= 2.00000000E+03
  1  1     7.22296753E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.30299799E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.54220076E-01   # Yu(3,3)
Block IMYu Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(Yu(1,1))
  1  2     0.00000000E+00   # Im(Yu(1,2))
  1  3     0.00000000E+00   # Im(Yu(1,3))
  2  1     0.00000000E+00   # Im(Yu(2,1))
  2  2     0.00000000E+00   # Im(Yu(2,2))
  2  3     0.00000000E+00   # Im(Yu(2,3))
  3  1     0.00000000E+00   # Im(Yu(3,1))
  3  2     0.00000000E+00   # Im(Yu(3,2))
  3  3     0.00000000E+00   # Im(Yu(3,3))
Block Yd Q= 2.00000000E+03
  1  1     6.87050263E-05   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     1.50427905E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     6.78851670E-02   # Yd(3,3)
Block IMYd Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(Yd(1,1))
  1  2     0.00000000E+00   # Im(Yd(1,2))
  1  3     0.00000000E+00   # Im(Yd(1,3))
  2  1     0.00000000E+00   # Im(Yd(2,1))
  2  2     0.00000000E+00   # Im(Yd(2,2))
  2  3     0.00000000E+00   # Im(Yd(2,3))
  3  1     0.00000000E+00   # Im(Yd(3,1))
  3  2     0.00000000E+00   # Im(Yd(3,2))
  3  3     0.00000000E+00   # Im(Yd(3,3))
Block Ye Q= 2.00000000E+03
  1  1     1.44626694E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     2.99042112E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     5.02942376E-02   # Ye(3,3)
Block IMYe Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(Ye(1,1))
  1  2     0.00000000E+00   # Im(Ye(1,2))
  1  3     0.00000000E+00   # Im(Ye(1,3))
  2  1     0.00000000E+00   # Im(Ye(2,1))
  2  2     0.00000000E+00   # Im(Ye(2,2))
  2  3     0.00000000E+00   # Im(Ye(2,3))
  3  1     0.00000000E+00   # Im(Ye(3,1))
  3  2     0.00000000E+00   # Im(Ye(3,2))
  3  3     0.00000000E+00   # Im(Ye(3,3))
Block Te Q= 2.00000000E+03
  1  1     1.44626927E-01   # Re(TYe(1,1))
  1  2     0.00000000E+00   # Re(TYe(1,2))
  1  3     0.00000000E+00   # Re(TYe(1,3))
  2  1     0.00000000E+00   # Re(TYe(2,1))
  2  2     2.99042592E+01   # Re(TYe(2,2))
  2  3     0.00000000E+00   # Re(TYe(2,3))
  3  1     0.00000000E+00   # Re(TYe(3,1))
  3  2     0.00000000E+00   # Re(TYe(3,2))
  3  3     5.02943166E+02   # Re(TYe(3,3))
Block IMTe Q= 2.00000000E+03
  1  1    -1.28160489E-10   # Im(TYe(1,1))
  1  2     0.00000000E+00   # Im(TYe(1,2))
  1  3     0.00000000E+00   # Im(TYe(1,3))
  2  1     0.00000000E+00   # Im(TYe(2,1))
  2  2    -2.64973639E-08   # Im(TYe(2,2))
  2  3     0.00000000E+00   # Im(TYe(2,3))
  3  1     0.00000000E+00   # Im(TYe(3,1))
  3  2     0.00000000E+00   # Im(TYe(3,2))
  3  3    -4.35626607E-07   # Im(TYe(3,3))
Block Td Q= 2.00000000E+03
  1  1     6.87048270E-01   # Re(TYd(1,1))
  1  2     0.00000000E+00   # Re(TYd(1,2))
  1  3     0.00000000E+00   # Re(TYd(1,3))
  2  1     0.00000000E+00   # Re(TYd(2,1))
  2  2     1.50427469E+01   # Re(TYd(2,2))
  2  3     0.00000000E+00   # Re(TYd(2,3))
  3  1     0.00000000E+00   # Re(TYd(3,1))
  3  2     0.00000000E+00   # Re(TYd(3,2))
  3  3     6.78846562E+02   # Re(TYd(3,3))
Block IMTd Q= 2.00000000E+03
  1  1     8.17085641E-09   # Im(TYd(1,1))
  1  2     0.00000000E+00   # Im(TYd(1,2))
  1  3     0.00000000E+00   # Im(TYd(1,3))
  2  1     0.00000000E+00   # Im(TYd(2,1))
  2  2     1.78898783E-07   # Im(TYd(2,2))
  2  3     0.00000000E+00   # Im(TYd(2,3))
  3  1     0.00000000E+00   # Im(TYd(3,1))
  3  2     0.00000000E+00   # Im(TYd(3,2))
  3  3     1.78391086E-05   # Im(TYd(3,3))
Block Tu Q= 2.00000000E+03
  1  1     2.88919526E-03   # Re(TYu(1,1))
  1  2     0.00000000E+00   # Re(TYu(1,2))
  1  3     0.00000000E+00   # Re(TYu(1,3))
  2  1     0.00000000E+00   # Re(TYu(2,1))
  2  2     1.32120297E+00   # Re(TYu(2,2))
  2  3     0.00000000E+00   # Re(TYu(2,3))
  3  1     0.00000000E+00   # Re(TYu(3,1))
  3  2     0.00000000E+00   # Re(TYu(3,2))
  3  3     3.41688289E+02   # Re(TYu(3,3))
Block IMTu Q= 2.00000000E+03
  1  1     8.45729492E-10   # Im(TYu(1,1))
  1  2     0.00000000E+00   # Im(TYu(1,2))
  1  3     0.00000000E+00   # Im(TYu(1,3))
  2  1     0.00000000E+00   # Im(TYu(2,1))
  2  2     3.86744318E-07   # Im(TYu(2,2))
  2  3     0.00000000E+00   # Im(TYu(2,3))
  3  1     0.00000000E+00   # Im(TYu(3,1))
  3  2     0.00000000E+00   # Im(TYu(3,2))
  3  3     1.82838625E-06   # Im(TYu(3,3))
Block MSQ2 Q= 2.00000000E+03
  1  1     3.99999902E+06   # Re(mq2(1,1))
  1  2     0.00000000E+00   # Re(mq2(1,2))
  1  3     0.00000000E+00   # Re(mq2(1,3))
  2  1     0.00000000E+00   # Re(mq2(2,1))
  2  2     3.99999902E+06   # Re(mq2(2,2))
  2  3     0.00000000E+00   # Re(mq2(2,3))
  3  1     0.00000000E+00   # Re(mq2(3,1))
  3  2     0.00000000E+00   # Re(mq2(3,2))
  3  3     3.99999867E+06   # Re(mq2(3,3))
Block IMMSQ2 Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(mq2(1,1))
  1  2     0.00000000E+00   # Im(mq2(1,2))
  1  3     0.00000000E+00   # Im(mq2(1,3))
  2  1     0.00000000E+00   # Im(mq2(2,1))
  2  2     0.00000000E+00   # Im(mq2(2,2))
  2  3     0.00000000E+00   # Im(mq2(2,3))
  3  1     0.00000000E+00   # Im(mq2(3,1))
  3  2     0.00000000E+00   # Im(mq2(3,2))
  3  3     0.00000000E+00   # Im(mq2(3,3))
Block MSE2 Q= 2.00000000E+03
  1  1     3.99999997E+06   # Re(me2(1,1))
  1  2     0.00000000E+00   # Re(me2(1,2))
  1  3     0.00000000E+00   # Re(me2(1,3))
  2  1     0.00000000E+00   # Re(me2(2,1))
  2  2     3.99999997E+06   # Re(me2(2,2))
  2  3     0.00000000E+00   # Re(me2(2,3))
  3  1     0.00000000E+00   # Re(me2(3,1))
  3  2     0.00000000E+00   # Re(me2(3,2))
  3  3     3.99999997E+06   # Re(me2(3,3))
Block IMMSE2 Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(me2(1,1))
  1  2     0.00000000E+00   # Im(me2(1,2))
  1  3     0.00000000E+00   # Im(me2(1,3))
  2  1     0.00000000E+00   # Im(me2(2,1))
  2  2     0.00000000E+00   # Im(me2(2,2))
  2  3     0.00000000E+00   # Im(me2(2,3))
  3  1     0.00000000E+00   # Im(me2(3,1))
  3  2     0.00000000E+00   # Im(me2(3,2))
  3  3     0.00000000E+00   # Im(me2(3,3))
Block MSL2 Q= 2.00000000E+03
  1  1     3.99999998E+06   # Re(ml2(1,1))
  1  2     0.00000000E+00   # Re(ml2(1,2))
  1  3     0.00000000E+00   # Re(ml2(1,3))
  2  1     0.00000000E+00   # Re(ml2(2,1))
  2  2     3.99999998E+06   # Re(ml2(2,2))
  2  3     0.00000000E+00   # Re(ml2(2,3))
  3  1     0.00000000E+00   # Re(ml2(3,1))
  3  2     0.00000000E+00   # Re(ml2(3,2))
  3  3     3.99999997E+06   # Re(ml2(3,3))
Block IMMSL2 Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(ml2(1,1))
  1  2     0.00000000E+00   # Im(ml2(1,2))
  1  3     0.00000000E+00   # Im(ml2(1,3))
  2  1     0.00000000E+00   # Im(ml2(2,1))
  2  2     0.00000000E+00   # Im(ml2(2,2))
  2  3     0.00000000E+00   # Im(ml2(2,3))
  3  1     0.00000000E+00   # Im(ml2(3,1))
  3  2     0.00000000E+00   # Im(ml2(3,2))
  3  3     0.00000000E+00   # Im(ml2(3,3))
Block MSU2 Q= 2.00000000E+03
  1  1     3.99999902E+06   # Re(mu2(1,1))
  1  2     0.00000000E+00   # Re(mu2(1,2))
  1  3     0.00000000E+00   # Re(mu2(1,3))
  2  1     0.00000000E+00   # Re(mu2(2,1))
  2  2     3.99999902E+06   # Re(mu2(2,2))
  2  3     0.00000000E+00   # Re(mu2(2,3))
  3  1     0.00000000E+00   # Re(mu2(3,1))
  3  2     0.00000000E+00   # Re(mu2(3,2))
  3  3     3.99999822E+06   # Re(mu2(3,3))
Block IMMSU2 Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(mu2(1,1))
  1  2     0.00000000E+00   # Im(mu2(1,2))
  1  3     0.00000000E+00   # Im(mu2(1,3))
  2  1     0.00000000E+00   # Im(mu2(2,1))
  2  2     0.00000000E+00   # Im(mu2(2,2))
  2  3     0.00000000E+00   # Im(mu2(2,3))
  3  1     0.00000000E+00   # Im(mu2(3,1))
  3  2     0.00000000E+00   # Im(mu2(3,2))
  3  3     0.00000000E+00   # Im(mu2(3,3))
Block MSD2 Q= 2.00000000E+03
  1  1     3.99999901E+06   # Re(md2(1,1))
  1  2     0.00000000E+00   # Re(md2(1,2))
  1  3     0.00000000E+00   # Re(md2(1,3))
  2  1     0.00000000E+00   # Re(md2(2,1))
  2  2     3.99999901E+06   # Re(md2(2,2))
  2  3     0.00000000E+00   # Re(md2(2,3))
  3  1     0.00000000E+00   # Re(md2(3,1))
  3  2     0.00000000E+00   # Re(md2(3,2))
  3  3     3.99999909E+06   # Re(md2(3,3))
Block IMMSD2 Q= 2.00000000E+03
  1  1     0.00000000E+00   # Im(md2(1,1))
  1  2     0.00000000E+00   # Im(md2(1,2))
  1  3     0.00000000E+00   # Im(md2(1,3))
  2  1     0.00000000E+00   # Im(md2(2,1))
  2  2     0.00000000E+00   # Im(md2(2,2))
  2  3     0.00000000E+00   # Im(md2(2,3))
  3  1     0.00000000E+00   # Im(md2(3,1))
  3  2     0.00000000E+00   # Im(md2(3,2))
  3  3     0.00000000E+00   # Im(md2(3,3))
Block Phases Q= 2.00000000E+03
     1     9.99688036E-01   # Re(PhaseGlu)
Block IMPhases Q= 2.00000000E+03
     1     2.49766003E-02   # Im(PhaseGlu)
Block MASS
   1000021     2.12248650E+03   # Glu
        24     8.03577013E+01   # VWm
   1000024     1.94149080E+03   # Cha(1)
   1000037     2.06815666E+03   # Cha(2)
        37     2.01410464E+03   # Hpm(2)
   1000012     2.01444544E+03   # Sv(1)
   1000014     2.01453694E+03   # Sv(2)
   1000016     2.01453726E+03   # Sv(3)
   1000022    -1.92660987E+03   # Chi(1)
   1000023    -1.98974263E+03   # Chi(2)
   1000025    -2.00050983E+03   # Chi(3)
   1000035    -2.07532614E+03   # Chi(4)
        25     1.05649119E+02   # hh(2)
        35     2.01213397E+03   # hh(3)
        36     2.01249453E+03   # hh(4)
   1000001     2.07062699E+03   # Sd(1)
   1000003     2.07100950E+03   # Sd(2)
   1000005     2.07100970E+03   # Sd(3)
   2000001     2.08385152E+03   # Sd(4)
   2000003     2.08482184E+03   # Sd(5)
   2000005     2.08489187E+03   # Sd(6)
   1000011     2.00757312E+03   # Se(1)
   1000013     2.00777308E+03   # Se(2)
   1000015     2.00777413E+03   # Se(3)
   2000011     2.01624667E+03   # Se(4)
   2000013     2.01632849E+03   # Se(5)
   2000015     2.01633906E+03   # Se(6)
   1000002     2.07187112E+03   # Su(1)
   1000004     2.07187113E+03   # Su(2)
   1000006     2.07242490E+03   # Su(3)
   2000002     2.08366996E+03   # Su(4)
   2000004     2.08444844E+03   # Su(5)
   2000006     2.09071039E+03   # Su(6)
Block UMIX
  1  1    -6.74308528E-01   # Re(UM(1,1))
  1  2     7.37548381E-01   # Re(UM(1,2))
  2  1     7.38430942E-01   # Re(UM(2,1))
  2  2     6.73554514E-01   # Re(UM(2,2))
Block VMIX
  1  1    -6.89511561E-01   # Re(UP(1,1))
  1  2     7.23719240E-01   # Re(UP(1,2))
  2  1     7.23107866E-01   # Re(UP(2,1))
  2  2     6.90094531E-01   # Re(UP(2,2))
Block DSQMIX
  1  1     0.00000000E+00   # Re(ZD(1,1))
  1  2    -1.02049981E-16   # Re(ZD(1,2))
  1  3    -1.64836623E-03   # Re(ZD(1,3))
  1  4    -1.50618287E-33   # Re(ZD(1,4))
  1  5    -3.96738074E-15   # Re(ZD(1,5))
  1  6    -8.16343871E-01   # Re(ZD(1,6))
  2  1     0.00000000E+00   # Re(ZD(2,1))
  2  2     3.97181164E-04   # Re(ZD(2,2))
  2  3     5.67887989E-15   # Re(ZD(2,3))
  2  4    -4.74338409E-20   # Re(ZD(2,4))
  2  5     1.54110923E-01   # Re(ZD(2,5))
  2  6    -6.20537550E-15   # Re(ZD(2,6))
  3  1    -1.97836577E-05   # Re(ZD(3,1))
  3  2     0.00000000E+00   # Re(ZD(3,2))
  3  3    -3.54618571E-35   # Re(ZD(3,3))
  3  4     2.52921286E-01   # Re(ZD(3,4))
  3  5     5.14079088E-20   # Re(ZD(3,5))
  3  6    -1.37363525E-33   # Re(ZD(3,6))
  4  1     0.00000000E+00   # Re(ZD(4,1))
  4  2     5.90085390E-14   # Re(ZD(4,2))
  4  3    -9.58023202E-01   # Re(ZD(4,3))
  4  4    -3.51971449E-32   # Re(ZD(4,4))
  4  5    -2.20215547E-14   # Re(ZD(4,5))
  4  6    -5.64034393E-03   # Re(ZD(4,6))
  5  1     0.00000000E+00   # Re(ZD(5,1))
  5  2     9.18747130E-01   # Re(ZD(5,2))
  5  3     6.14772771E-14   # Re(ZD(5,3))
  5  4     1.99493487E-23   # Re(ZD(5,4))
  5  5    -6.48147500E-05   # Re(ZD(5,5))
  5  6     2.72965740E-16   # Re(ZD(5,6))
  6  1     1.00000000E+00   # Re(ZD(6,1))
  6  2     0.00000000E+00   # Re(ZD(6,2))
  6  3    -7.01565243E-40   # Re(ZD(6,3))
  6  4     5.00370814E-06   # Re(ZD(6,4))
  6  5     1.01703647E-24   # Re(ZD(6,5))
  6  6    -2.71755296E-38   # Re(ZD(6,6))
Block SELMIX
  1  1    -0.00000000E+00   # Re(ZE(1,1))
  1  2     3.13284975E-16   # Re(ZE(1,2))
  1  3    -3.39214507E-02   # Re(ZE(1,3))
  1  4     1.03896833E-31   # Re(ZE(1,4))
  1  5    -4.13604096E-15   # Re(ZE(1,5))
  1  6     5.14601937E-01   # Re(ZE(1,6))
  2  1     0.00000000E+00   # Re(ZE(2,1))
  2  2     2.03843013E-03   # Re(ZE(2,2))
  2  3     2.44285156E-15   # Re(ZE(2,3))
  2  4     6.81532055E-19   # Re(ZE(2,4))
  2  5    -5.88927093E-01   # Re(ZE(2,5))
  2  6     8.72426912E-15   # Re(ZE(2,6))
  3  1     9.88700435E-06   # Re(ZE(3,1))
  3  2     0.00000000E+00   # Re(ZE(3,2))
  3  3    -3.23176259E-34   # Re(ZE(3,3))
  3  4    -5.26454855E-01   # Re(ZE(3,4))
  3  5    -1.58979968E-22   # Re(ZE(3,5))
  3  6    -1.79844093E-36   # Re(ZE(3,6))
  4  1    -0.00000000E+00   # Re(ZE(4,1))
  4  2     9.11266516E-13   # Re(ZE(4,2))
  4  3    -9.99325981E-01   # Re(ZE(4,3))
  4  4     3.04350187E-28   # Re(ZE(4,4))
  4  5     3.55937106E-15   # Re(ZE(4,5))
  4  6    -1.82567993E-02   # Re(ZE(4,6))
  5  1     0.00000000E+00   # Re(ZE(5,1))
  5  2     9.97163661E-01   # Re(ZE(5,2))
  5  3     9.13460204E-13   # Re(ZE(5,3))
  5  4     3.33066210E-16   # Re(ZE(5,4))
  5  5     1.20390140E-03   # Re(ZE(5,5))
  5  6     1.65052435E-14   # Re(ZE(5,6))
  6  1     1.00000000E+00   # Re(ZE(6,1))
  6  2     0.00000000E+00   # Re(ZE(6,2))
  6  3     3.19524507E-39   # Re(ZE(6,3))
  6  4     5.20506144E-06   # Re(ZE(6,4))
  6  5     1.57183563E-27   # Re(ZE(6,5))
  6  6     1.77811933E-41   # Re(ZE(6,6))
Block SCALARMIX
  1  1     6.58317063E-07   # ZH(1,1)
  1  2     5.75447562E-07   # ZH(1,2)
  1  3    -2.04921897E-01   # ZH(1,3)
  1  4     9.78778328E-01   # ZH(1,4)
  2  1    -2.05941537E-01   # ZH(2,1)
  2  2    -9.78564297E-01   # ZH(2,2)
  2  3     2.56774341E-06   # ZH(2,3)
  2  4     1.25143159E-06   # ZH(2,4)
  3  1    -3.20729293E-03   # ZH(3,1)
  3  2     6.72153232E-04   # ZH(3,2)
  3  3    -9.78773073E-01   # ZH(3,3)
  3  4    -2.04920795E-01   # ZH(3,4)
  4  1    -9.78559041E-01   # ZH(4,1)
  4  2     2.05940440E-01   # ZH(4,2)
  4  3     3.20731619E-03   # ZH(4,3)
  4  4     6.72036756E-04   # ZH(4,4)
Block NMIX
  1  1    -1.17411817E-02   # Re(ZN(1,1))
  1  2     1.27757982E-02   # Re(ZN(1,2))
  1  3    -3.33170047E-02   # Re(ZN(1,3))
  1  4    -8.21772589E-03   # Re(ZN(1,4))
  2  1     2.09697462E-02   # Re(ZN(2,1))
  2  2     1.20739572E-02   # Re(ZN(2,2))
  2  3    -1.23544567E-02   # Re(ZN(2,3))
  2  4    -7.22039804E-03   # Re(ZN(2,4))
  3  1     6.01397605E-03   # Re(ZN(3,1))
  3  2    -1.08795375E-02   # Re(ZN(3,2))
  3  3    -7.05302993E-01   # Re(ZN(3,3))
  3  4    -7.07073094E-01   # Re(ZN(3,4))
  4  1     6.95368868E-03   # Re(ZN(4,1))
  4  2    -1.50631891E-02   # Re(ZN(4,2))
  4  3    -3.14432341E-02   # Re(ZN(4,3))
  4  4    -8.63878097E-03   # Re(ZN(4,4))
Block CHARGEMIX
  1  1     1.98170113E-01   # Re(ZP(1,1))
  1  2    -9.80167081E-01   # Re(ZP(1,2))
  2  1     9.80167642E-01   # Re(ZP(2,1))
  2  2     1.98169999E-01   # Re(ZP(2,2))
Block USQMIX
  1  1    -0.00000000E+00   # Re(ZU(1,1))
  1  2     2.02870323E-03   # Re(ZU(1,2))
  1  3    -7.95342182E-14   # Re(ZU(1,3))
  1  4     2.42381861E-19   # Re(ZU(1,4))
  1  5     9.88815806E-01   # Re(ZU(1,5))
  1  6    -1.66170608E-13   # Re(ZU(1,6))
  2  1    -4.43695922E-06   # Re(ZU(2,1))
  2  2     0.00000000E+00   # Re(ZU(2,2))
  2  3     2.31285481E-33   # Re(ZU(2,3))
  2  4    -9.90641516E-01   # Re(ZU(2,4))
  2  5     6.03639393E-19   # Re(ZU(2,5))
  2  6     1.37788954E-31   # Re(ZU(2,6))
  3  1     0.00000000E+00   # Re(ZU(3,1))
  3  2    -1.11984484E-15   # Re(ZU(3,2))
  3  3     2.21063653E-01   # Re(ZU(3,3))
  3  4     4.66267870E-31   # Re(ZU(3,4))
  3  5     1.38948153E-13   # Re(ZU(3,5))
  3  6     6.81758753E-01   # Re(ZU(3,6))
  4  1     1.00000000E+00   # Re(ZU(4,1))
  4  2     0.00000000E+00   # Re(ZU(4,2))
  4  3     1.02620424E-38   # Re(ZU(4,3))
  4  4    -4.39543601E-06   # Re(ZU(4,4))
  4  5     2.67832337E-24   # Re(ZU(4,5))
  4  6     6.11363969E-37   # Re(ZU(4,6))
  5  1     0.00000000E+00   # Re(ZU(5,1))
  5  2     9.99923977E-01   # Re(ZU(5,2))
  5  3    -5.10272788E-15   # Re(ZU(5,3))
  5  4     2.22044570E-16   # Re(ZU(5,4))
  5  5    -2.00636162E-03   # Re(ZU(5,5))
  5  6     3.49098316E-15   # Re(ZU(5,6))
  6  1     0.00000000E+00   # Re(ZU(6,1))
  6  2    -7.24045556E-15   # Re(ZU(6,2))
  6  3    -7.61359167E-01   # Re(ZU(6,3))
  6  4    -1.53685523E-30   # Re(ZU(6,4))
  6  5    -1.30781855E-14   # Re(ZU(6,5))
  6  6     3.23925874E-01   # Re(ZU(6,6))
Block SNUMIX
  1  1     0.00000000E+00   # Re(ZV(1,1))
  1  2     4.93507117E-16   # Re(ZV(1,2))
  1  3    -8.13200476E-01   # Re(ZV(1,3))
  2  1     0.00000000E+00   # Re(ZV(2,1))
  2  2    -3.25932406E-01   # Re(ZV(2,2))
  2  3    -6.47877747E-16   # Re(ZV(2,3))
  3  1     1.00000000E+00   # Re(ZV(3,1))
  3  2     0.00000000E+00   # Re(ZV(3,2))
  3  3     0.00000000E+00   # Re(ZV(3,3))
Block FlexibleSUSYOutput
     0     0.00000000E+00   # HighScale
     1     2.00000000E+03   # SUSYScale
     2     9.11876000E+01   # LowScale
Block FlexibleSUSYLowEnergy Q= 2.00000000E+03
    21     0.00000000E+00   # Delta(g-2)_muon/2 FlexibleSUSY
    23     0.00000000E+00   # electric dipole moment of Fe(0) [1/GeV]
    24     0.00000000E+00   # electric dipole moment of Fe(1) [1/GeV]
    25     0.00000000E+00   # electric dipole moment of Fe(2) [1/GeV]
    26     0.00000000E+00   # BR(Fe1 -> Fe0 VP)
Block EFFHIGGSCOUPLINGS
       25       22       22     0.00000000E+00   # Abs(effective H-Photon-Photon coupling)
       35       22       22     0.00000000E+00   # Abs(effective H-Photon-Photon coupling)
       36       22       22     0.00000000E+00   # Abs(effective H-Photon-Photon coupling)
       25       21       21     0.00000000E+00   # Abs(effective H-Gluon-Gluon coupling)
       35       21       21     0.00000000E+00   # Abs(effective H-Gluon-Gluon coupling)
       36       21       21     0.00000000E+00   # Abs(effective H-Gluon-Gluon coupling)
Block ALPHA
          -1.36337055E+00   # ArcSin(Pole(ZH(2,2)))
Block HMIX Q= 2.00000000E+03
     1     1.99999979E+03   # Re(Mu)
     2     4.77520537E+00   # vu/vd
     3     2.44022303E+02   # Sqrt(Sqr(vd) + Sqr(vu))
   101     8.02468438E+05   # Re(BMu)
   102     5.00169780E+01   # vd
   103     2.38841341E+02   # vu
Block ImHMIX Q= 2.00000000E+03
     1     9.99999897E+01   # Im(Mu)
   101    -8.04093338E+03   # Im(BMu)
Block Au Q= 2.00000000E+03
  1  1     4.00001141E+02   # Re(TYu(1,1)/Yu(1,1))
  2  2     4.00001141E+02   # Re(TYu(2,2)/Yu(2,2))
  3  3     4.00000303E+02   # Re(TYu(3,3)/Yu(3,3))
Block Ad Q= 2.00000000E+03
  1  1     9.99997100E+03   # Re(TYd(1,1)/Yd(1,1))
  2  2     9.99997100E+03   # Re(TYd(2,2)/Yd(2,2))
  3  3     9.99992475E+03   # Re(TYd(3,3)/Yd(3,3))
Block Ae Q= 2.00000000E+03
  1  1     1.00000160E+04   # Re(TYe(1,1)/Ye(1,1))
  2  2     1.00000160E+04   # Re(TYe(2,2)/Ye(2,2))
  3  3     1.00000157E+04   # Re(TYe(3,3)/Ye(3,3))
Block ImAu Q= 2.00000000E+03
  1  1     1.17088923E-04   # Im(TYu(1,1)/Yu(1,1))
  2  2     1.17088875E-04   # Im(TYu(2,2)/Yu(2,2))
  3  3     2.14041592E-06   # Im(TYu(3,3)/Yu(3,3))
Block ImAd Q= 2.00000000E+03
  1  1     1.18926618E-04   # Im(TYd(1,1)/Yd(1,1))
  2  2     1.18926593E-04   # Im(TYd(2,2)/Yd(2,2))
  3  3     2.62783601E-04   # Im(TYd(3,3)/Yd(3,3))
Block ImAe Q= 2.00000000E+03
  1  1    -8.86146846E-06   # Im(TYe(1,1)/Ye(1,1))
  2  2    -8.86074664E-06   # Im(TYe(2,2)/Ye(2,2))
  3  3    -8.66156099E-06   # Im(TYe(3,3)/Ye(3,3))
Block MSOFT Q= 2.00000000E+03
     1     2.00000004E+03   # Re(MassB)
     2     2.00000002E+03   # Re(MassWB)
     3     1.99999977E+03   # Re(MassG)
    21    -1.88645947E+05   # mHd2
    22    -3.77261113E+06   # mHu2
    31     1.99999999E+03   # SignedAbsSqrt(Re(ml2(1,1)))
    32     1.99999999E+03   # SignedAbsSqrt(Re(ml2(2,2)))
    33     1.99999999E+03   # SignedAbsSqrt(Re(ml2(3,3)))
    34     1.99999999E+03   # SignedAbsSqrt(Re(me2(1,1)))
    35     1.99999999E+03   # SignedAbsSqrt(Re(me2(2,2)))
    36     1.99999999E+03   # SignedAbsSqrt(Re(me2(3,3)))
    41     1.99999975E+03   # SignedAbsSqrt(Re(mq2(1,1)))
    42     1.99999975E+03   # SignedAbsSqrt(Re(mq2(2,2)))
    43     1.99999967E+03   # SignedAbsSqrt(Re(mq2(3,3)))
    44     1.99999975E+03   # SignedAbsSqrt(Re(mu2(1,1)))
    45     1.99999975E+03   # SignedAbsSqrt(Re(mu2(2,2)))
    46     1.99999955E+03   # SignedAbsSqrt(Re(mu2(3,3)))
    47     1.99999975E+03   # SignedAbsSqrt(Re(md2(1,1)))
    48     1.99999975E+03   # SignedAbsSqrt(Re(md2(2,2)))
    49     1.99999977E+03   # SignedAbsSqrt(Re(md2(3,3)))
Block ImMSOFT Q= 2.00000000E+03
     1     1.00000002E+02   # Im(MassB)
     2     1.00000001E+02   # Im(MassWB)
     3     9.99999883E+01   # Im(MassG)
    31     0.00000000E+00   # SignedAbsSqrt(Im(ml2(1,1)))
    32     0.00000000E+00   # SignedAbsSqrt(Im(ml2(2,2)))
    33     0.00000000E+00   # SignedAbsSqrt(Im(ml2(3,3)))
    34     0.00000000E+00   # SignedAbsSqrt(Im(me2(1,1)))
    35     0.00000000E+00   # SignedAbsSqrt(Im(me2(2,2)))
    36     0.00000000E+00   # SignedAbsSqrt(Im(me2(3,3)))
    41     0.00000000E+00   # SignedAbsSqrt(Im(mq2(1,1)))
    42     0.00000000E+00   # SignedAbsSqrt(Im(mq2(2,2)))
    43     0.00000000E+00   # SignedAbsSqrt(Im(mq2(3,3)))
    44     0.00000000E+00   # SignedAbsSqrt(Im(mu2(1,1)))
    45     0.00000000E+00   # SignedAbsSqrt(Im(mu2(2,2)))
    46     0.00000000E+00   # SignedAbsSqrt(Im(mu2(3,3)))
    47     0.00000000E+00   # SignedAbsSqrt(Im(md2(1,1)))
    48     0.00000000E+00   # SignedAbsSqrt(Im(md2(2,2)))
    49     0.00000000E+00   # SignedAbsSqrt(Im(md2(3,3)))
)";

   std::stringstream istr(slha_input);

   MSSMCPV_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   MSSMCPV_input_parameters input;
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

   MSSMCPV_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   MSSMCPV_slha m = std::get<0>(spectrum_generator.get_models_slha());

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   MSSMCPV_decays decays_with_HO(m, qedqcd, physical_input, SM_higher_order_corrections::enable);

   // scalar Higgs

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFdFd(&m, 1, 2, 2),
                              0.0015435164852836149, 3e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFuFu(&m, 1, 1, 1),
                              8.3735251701183208e-05, 2e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFeFe(&m, 1, 2, 2),
                              0.00022695651730241407, 3e-13);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 1),
   //                            4.271655828335071e-05, 2e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 1),
                              6.0919503405648103e-05, 1e-3);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 1),
   //                            1.3778143534859871e-06, 4e-11);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 1),
                              5.7304088223134617e-06, 1e-3);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VGVG(&m, 1), 0.00017296000561389487, 2e-10);
   // h -> gamma gamma
   // without 2L QCD for squark
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 6.3284545616000571e-06, 4e-11);
   // with 2L QCD for squark
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 1), 6.2715279664521928e-06, 4e-10);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVZ(&m, 1), 5.7227729785272779e-07, 7e-11);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   MSSMCPV_decays decays_without_HO(m, qedqcd, physical_input, SM_higher_order_corrections::disable);

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFdFd(&m, 1, 2, 2),
                              0.0012302546720318283, 2e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFuFu(&m, 1, 1, 1),
                              6.5889191532059478e-05, 2e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFeFe(&m, 1, 2, 2),
                              0.00022454719334137702, 3e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVP(&m, 1), 6.1880194357760738e-06, 4e-10);

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VGVG(&m, 1), 6.1972398675528688e-05, 2e-10);

   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVZ(&m, 1), 5.7059182799004216e-07, 7e-11);
}
