
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
     2   2.6.1
     5   MSSMCPV
     9   4.14.5
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
  1  1     7.22296754E-06   # Yu(1,1)
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
  1  1     1.44626926E-01   # Re(TYe(1,1))
  1  2     0.00000000E+00   # Re(TYe(1,2))
  1  3     0.00000000E+00   # Re(TYe(1,3))
  2  1     0.00000000E+00   # Re(TYe(2,1))
  2  2     2.99042592E+01   # Re(TYe(2,2))
  2  3     0.00000000E+00   # Re(TYe(2,3))
  3  1     0.00000000E+00   # Re(TYe(3,1))
  3  2     0.00000000E+00   # Re(TYe(3,2))
  3  3     5.02943166E+02   # Re(TYe(3,3))
Block IMTe Q= 2.00000000E+03
  1  1    -1.28097413E-10   # Im(TYe(1,1))
  1  2     0.00000000E+00   # Im(TYe(1,2))
  1  3     0.00000000E+00   # Im(TYe(1,3))
  2  1     0.00000000E+00   # Im(TYe(2,1))
  2  2    -2.64843217E-08   # Im(TYe(2,2))
  2  3     0.00000000E+00   # Im(TYe(2,3))
  3  1     0.00000000E+00   # Im(TYe(3,1))
  3  2     0.00000000E+00   # Im(TYe(3,2))
  3  3    -4.35407343E-07   # Im(TYe(3,3))
Block Td Q= 2.00000000E+03
  1  1     6.87048270E-01   # Re(TYd(1,1))
  1  2     0.00000000E+00   # Re(TYd(1,2))
  1  3     0.00000000E+00   # Re(TYd(1,3))
  2  1     0.00000000E+00   # Re(TYd(2,1))
  2  2     1.50427469E+01   # Re(TYd(2,2))
  2  3     0.00000000E+00   # Re(TYd(2,3))
  3  1     0.00000000E+00   # Re(TYd(3,1))
  3  2     0.00000000E+00   # Re(TYd(3,2))
  3  3     6.78846561E+02   # Re(TYd(3,3))
Block IMTd Q= 2.00000000E+03
  1  1     8.17236531E-09   # Im(TYd(1,1))
  1  2     0.00000000E+00   # Im(TYd(1,2))
  1  3     0.00000000E+00   # Im(TYd(1,3))
  2  1     0.00000000E+00   # Im(TYd(2,1))
  2  2     1.78931820E-07   # Im(TYd(2,2))
  2  3     0.00000000E+00   # Im(TYd(2,3))
  3  1     0.00000000E+00   # Im(TYd(3,1))
  3  2     0.00000000E+00   # Im(TYd(3,2))
  3  3     1.78404600E-05   # Im(TYd(3,3))
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
  1  1     8.45875713E-10   # Im(TYu(1,1))
  1  2     0.00000000E+00   # Im(TYu(1,2))
  1  3     0.00000000E+00   # Im(TYu(1,3))
  2  1     0.00000000E+00   # Im(TYu(2,1))
  2  2     3.86811184E-07   # Im(TYu(2,2))
  2  3     0.00000000E+00   # Im(TYu(2,3))
  3  1     0.00000000E+00   # Im(TYu(3,1))
  3  2     0.00000000E+00   # Im(TYu(3,2))
  3  3     1.84365506E-06   # Im(TYu(3,3))
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
  1  1    -7.27595761E-12   # Im(mu2(1,1))
  1  2     0.00000000E+00   # Im(mu2(1,2))
  1  3     0.00000000E+00   # Im(mu2(1,3))
  2  1     0.00000000E+00   # Im(mu2(2,1))
  2  2     0.00000000E+00   # Im(mu2(2,2))
  2  3     0.00000000E+00   # Im(mu2(2,3))
  3  1     0.00000000E+00   # Im(mu2(3,1))
  3  2     0.00000000E+00   # Im(mu2(3,2))
  3  3     7.27595761E-12   # Im(mu2(3,3))
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
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.30033519E-04   # Fe(1)
        13     1.07458698E-01   # Fe(2)
        15     1.78826361E+00   # Fe(3)
         1     4.25777269E-03   # Fd(1)
         3     8.44043803E-02   # Fd(2)
         5     3.35408940E+00   # Fd(3)
         2     2.21823311E-03   # Fu(1)
         4     8.26552842E-01   # Fu(2)
         6     1.74023906E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.10363042E+01   # VZ
Block UMIX
  1  1    -6.74308528E-01   # Re(UM(1,1))
  1  2     7.37548380E-01   # Re(UM(1,2))
  2  1     7.38430942E-01   # Re(UM(2,1))
  2  2     6.73554514E-01   # Re(UM(2,2))
Block VMIX
  1  1    -6.89511561E-01   # Re(UP(1,1))
  1  2     7.23719240E-01   # Re(UP(1,2))
  2  1     7.23107865E-01   # Re(UP(2,1))
  2  2     6.90094531E-01   # Re(UP(2,2))
Block DSQMIX
  1  1     0.00000000E+00   # Re(ZD(1,1))
  1  2     4.05876750E-16   # Re(ZD(1,2))
  1  3     2.25193872E-02   # Re(ZD(1,3))
  1  4     4.52983839E-32   # Re(ZD(1,4))
  1  5     1.60605380E-14   # Re(ZD(1,5))
  1  6    -5.79521697E-01   # Re(ZD(1,6))
  2  1     0.00000000E+00   # Re(ZD(2,1))
  2  2     4.29985522E-04   # Re(ZD(2,2))
  2  3     9.80325892E-17   # Re(ZD(2,3))
  2  4     4.94933664E-20   # Re(ZD(2,4))
  2  5    -3.67898261E-01   # Re(ZD(2,5))
  2  6     1.32024654E-15   # Re(ZD(2,6))
  3  1     1.97836576E-05   # Re(ZD(3,1))
  3  2     0.00000000E+00   # Re(ZD(3,2))
  3  3    -2.44243972E-34   # Re(ZD(3,3))
  3  4    -2.52921273E-01   # Re(ZD(3,4))
  3  5     1.85262697E-20   # Re(ZD(3,5))
  3  6    -1.26186452E-34   # Re(ZD(3,6))
  4  1    -0.00000000E+00   # Re(ZD(4,1))
  4  2     5.94241871E-15   # Re(ZD(4,2))
  4  3    -2.60399345E-01   # Re(ZD(4,3))
  4  4     6.88934165E-31   # Re(ZD(4,4))
  4  5    -1.09481417E-15   # Re(ZD(4,5))
  4  6    -2.17123218E-02   # Re(ZD(4,6))
  5  1     0.00000000E+00   # Re(ZD(5,1))
  5  2     9.92682840E-01   # Re(ZD(5,2))
  5  3     1.54570710E-15   # Re(ZD(5,3))
  5  4     1.11022291E-16   # Re(ZD(5,4))
  5  5     1.59356966E-04   # Re(ZD(5,5))
  5  6     3.64350311E-16   # Re(ZD(5,6))
  6  1     1.00000000E+00   # Re(ZD(6,1))
  6  2     0.00000000E+00   # Re(ZD(6,2))
  6  3     4.83203912E-39   # Re(ZD(6,3))
  6  4     5.00370787E-06   # Re(ZD(6,4))
  6  5    -3.66517376E-25   # Re(ZD(6,5))
  6  6     2.49642955E-39   # Re(ZD(6,6))
Block SELMIX
  1  1     0.00000000E+00   # Re(ZE(1,1))
  1  2    -1.91618813E-16   # Re(ZE(1,2))
  1  3     2.53949812E-02   # Re(ZE(1,3))
  1  4     7.99510941E-29   # Re(ZE(1,4))
  1  5     3.96104865E-15   # Re(ZE(1,5))
  1  6    -9.57207848E-01   # Re(ZE(1,6))
  2  1     0.00000000E+00   # Re(ZE(2,1))
  2  2    -1.99734605E-03   # Re(ZE(2,2))
  2  3     9.88180987E-15   # Re(ZE(2,3))
  2  4    -3.85177746E-18   # Re(ZE(2,4))
  2  5     3.33349838E-01   # Re(ZE(2,5))
  2  6     2.79952017E-14   # Re(ZE(2,6))
  3  1    -9.88700429E-06   # Re(ZE(3,1))
  3  2    -1.55431223E-15   # Re(ZE(3,2))
  3  3    -1.83840514E-28   # Re(ZE(3,3))
  3  4     5.26454847E-01   # Re(ZE(3,4))
  3  5     2.34083863E-17   # Re(ZE(3,5))
  3  6     1.38970100E-28   # Re(ZE(3,6))
  4  1     0.00000000E+00   # Re(ZE(4,1))
  4  2    -4.43153669E-13   # Re(ZE(4,2))
  4  3     5.01349084E-01   # Re(ZE(4,3))
  4  4     3.94626733E-30   # Re(ZE(4,4))
  4  5    -8.03083558E-15   # Re(ZE(4,5))
  4  6     3.39104666E-02   # Re(ZE(4,6))
  5  1     0.00000000E+00   # Re(ZE(5,1))
  5  2    -9.86503896E-01   # Re(ZE(5,2))
  5  3    -2.28203261E-13   # Re(ZE(5,3))
  5  4    -3.88571854E-16   # Re(ZE(5,4))
  5  5    -7.77112837E-04   # Re(ZE(5,5))
  5  6    -1.50698730E-14   # Re(ZE(5,6))
  6  1     1.00000000E+00   # Re(ZE(6,1))
  6  2    -1.53674917E-20   # Re(ZE(6,2))
  6  3    -1.81763195E-33   # Re(ZE(6,3))
  6  4     5.20506133E-06   # Re(ZE(6,4))
  6  5     2.31438815E-22   # Re(ZE(6,5))
  6  6     1.37399797E-33   # Re(ZE(6,6))
Block SCALARMIX
  1  1     6.58317063E-07   # ZH(1,1)
  1  2     5.75447560E-07   # ZH(1,2)
  1  3    -2.04921897E-01   # ZH(1,3)
  1  4     9.78778328E-01   # ZH(1,4)
  2  1    -2.05941537E-01   # ZH(2,1)
  2  2    -9.78564297E-01   # ZH(2,2)
  2  3     2.56774341E-06   # ZH(2,3)
  2  4     1.25143159E-06   # ZH(2,4)
  3  1    -3.20729294E-03   # ZH(3,1)
  3  2     6.72153232E-04   # ZH(3,2)
  3  3    -9.78773073E-01   # ZH(3,3)
  3  4    -2.04920795E-01   # ZH(3,4)
  4  1    -9.78559041E-01   # ZH(4,1)
  4  2     2.05940440E-01   # ZH(4,2)
  4  3     3.20731619E-03   # ZH(4,3)
  4  4     6.72036757E-04   # ZH(4,4)
Block NMIX
  1  1    -1.17411817E-02   # Re(ZN(1,1))
  1  2     1.27757982E-02   # Re(ZN(1,2))
  1  3    -3.33170047E-02   # Re(ZN(1,3))
  1  4    -8.21772588E-03   # Re(ZN(1,4))
  2  1     2.09697461E-02   # Re(ZN(2,1))
  2  2     1.20739572E-02   # Re(ZN(2,2))
  2  3    -1.23544567E-02   # Re(ZN(2,3))
  2  4    -7.22039797E-03   # Re(ZN(2,4))
  3  1     6.01397606E-03   # Re(ZN(3,1))
  3  2    -1.08795375E-02   # Re(ZN(3,2))
  3  3    -7.05302993E-01   # Re(ZN(3,3))
  3  4    -7.07073094E-01   # Re(ZN(3,4))
  4  1     6.95368867E-03   # Re(ZN(4,1))
  4  2    -1.50631891E-02   # Re(ZN(4,2))
  4  3    -3.14432341E-02   # Re(ZN(4,3))
  4  4    -8.63878097E-03   # Re(ZN(4,4))
Block CHARGEMIX
  1  1     1.98170113E-01   # Re(ZP(1,1))
  1  2    -9.80167081E-01   # Re(ZP(1,2))
  2  1     9.80167642E-01   # Re(ZP(2,1))
  2  2     1.98169999E-01   # Re(ZP(2,2))
Block USQMIX
  1  1     0.00000000E+00   # Re(ZU(1,1))
  1  2     1.99245866E-03   # Re(ZU(1,2))
  1  3     3.32241102E-14   # Re(ZU(1,3))
  1  4    -4.50484651E-19   # Re(ZU(1,4))
  1  5     9.98603891E-01   # Re(ZU(1,5))
  1  6     5.93054460E-14   # Re(ZU(1,6))
  2  1     4.43695921E-06   # Re(ZU(2,1))
  2  2     0.00000000E+00   # Re(ZU(2,2))
  2  3    -1.39141046E-33   # Re(ZU(2,3))
  2  4     9.90641516E-01   # Re(ZU(2,4))
  2  5    -1.69357107E-23   # Re(ZU(2,5))
  2  6    -3.47197908E-33   # Re(ZU(2,6))
  3  1     0.00000000E+00   # Re(ZU(3,1))
  3  2     4.35986043E-15   # Re(ZU(3,2))
  3  3    -3.60472795E-01   # Re(ZU(3,3))
  3  4    -9.85964623E-31   # Re(ZU(3,4))
  3  5     7.65773006E-14   # Re(ZU(3,5))
  3  6    -9.14458679E-01   # Re(ZU(3,6))
  4  1     1.00000000E+00   # Re(ZU(4,1))
  4  2     0.00000000E+00   # Re(ZU(4,2))
  4  3     6.17363147E-39   # Re(ZU(4,3))
  4  4    -4.39543600E-06   # Re(ZU(4,4))
  4  5     7.51430575E-29   # Re(ZU(4,5))
  4  6     1.54050296E-38   # Re(ZU(4,6))
  5  1     0.00000000E+00   # Re(ZU(5,1))
  5  2     9.82052145E-01   # Re(ZU(5,2))
  5  3    -1.66400697E-14   # Re(ZU(5,3))
  5  4    -2.22044148E-16   # Re(ZU(5,4))
  5  5    -2.02604005E-03   # Re(ZU(5,5))
  5  6     1.09391609E-14   # Re(ZU(5,6))
  6  1     0.00000000E+00   # Re(ZU(6,1))
  6  2    -1.91544388E-14   # Re(ZU(6,2))
  6  3    -9.32648116E-01   # Re(ZU(6,3))
  6  4     4.32828091E-30   # Re(ZU(6,4))
  6  5     5.44694874E-15   # Re(ZU(6,5))
  6  6     3.56501913E-01   # Re(ZU(6,6))
Block SNUMIX
  1  1     0.00000000E+00   # Re(ZV(1,1))
  1  2     3.92143553E-16   # Re(ZV(1,2))
  1  3    -3.20556426E-01   # Re(ZV(1,3))
  2  1     0.00000000E+00   # Re(ZV(2,1))
  2  2    -9.74685428E-01   # Re(ZV(2,2))
  2  3    -3.94618972E-16   # Re(ZV(2,3))
  3  1     1.00000000E+00   # Re(ZV(3,1))
  3  2     0.00000000E+00   # Re(ZV(3,2))
  3  3     0.00000000E+00   # Re(ZV(3,3))
Block UELMIX
  1  1     9.99999532E-01   # Re(ZEL(1,1))
  1  2     0.00000000E+00   # Re(ZEL(1,2))
  1  3     0.00000000E+00   # Re(ZEL(1,3))
  2  1     0.00000000E+00   # Re(ZEL(2,1))
  2  2     9.99999511E-01   # Re(ZEL(2,2))
  2  3     0.00000000E+00   # Re(ZEL(2,3))
  3  1     0.00000000E+00   # Re(ZEL(3,1))
  3  2     0.00000000E+00   # Re(ZEL(3,2))
  3  3     9.99999501E-01   # Re(ZEL(3,3))
Block UERMIX
  1  1     1.00000000E+00   # Re(ZER(1,1))
  1  2     0.00000000E+00   # Re(ZER(1,2))
  1  3     0.00000000E+00   # Re(ZER(1,3))
  2  1     0.00000000E+00   # Re(ZER(2,1))
  2  2     1.00000000E+00   # Re(ZER(2,2))
  2  3     0.00000000E+00   # Re(ZER(2,3))
  3  1     0.00000000E+00   # Re(ZER(3,1))
  3  2     0.00000000E+00   # Re(ZER(3,2))
  3  3     1.00000000E+00   # Re(ZER(3,3))
Block UDLMIX
  1  1     9.99999628E-01   # Re(ZDL(1,1))
  1  2     0.00000000E+00   # Re(ZDL(1,2))
  1  3     0.00000000E+00   # Re(ZDL(1,3))
  2  1     0.00000000E+00   # Re(ZDL(2,1))
  2  2    -4.55435498E-02   # Re(ZDL(2,2))
  2  3    -2.17623843E-17   # Re(ZDL(2,3))
  3  1     0.00000000E+00   # Re(ZDL(3,1))
  3  2     1.58539148E-19   # Re(ZDL(3,2))
  3  3     9.99999284E-01   # Re(ZDL(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(ZDR(1,1))
  1  2     0.00000000E+00   # Re(ZDR(1,2))
  1  3     0.00000000E+00   # Re(ZDR(1,3))
  2  1     0.00000000E+00   # Re(ZDR(2,1))
  2  2    -4.65090294E-02   # Re(ZDR(2,2))
  2  3    -9.50114445E-16   # Re(ZDR(2,3))
  3  1     0.00000000E+00   # Re(ZDR(3,1))
  3  2    -4.41889006E-17   # Re(ZDR(3,2))
  3  3     1.00000000E+00   # Re(ZDR(3,3))
Block UULMIX
  1  1     9.99999999E-01   # Re(ZUL(1,1))
  1  2     0.00000000E+00   # Re(ZUL(1,2))
  1  3     0.00000000E+00   # Re(ZUL(1,3))
  2  1     0.00000000E+00   # Re(ZUL(2,1))
  2  2     9.99999999E-01   # Re(ZUL(2,2))
  2  3     0.00000000E+00   # Re(ZUL(2,3))
  3  1     0.00000000E+00   # Re(ZUL(3,1))
  3  2     0.00000000E+00   # Re(ZUL(3,2))
  3  3     9.99999999E-01   # Re(ZUL(3,3))
Block UURMIX
  1  1     1.00000000E+00   # Re(ZUR(1,1))
  1  2     0.00000000E+00   # Re(ZUR(1,2))
  1  3     0.00000000E+00   # Re(ZUR(1,3))
  2  1     0.00000000E+00   # Re(ZUR(2,1))
  2  2     1.00000000E+00   # Re(ZUR(2,2))
  2  3     0.00000000E+00   # Re(ZUR(2,3))
  3  1     0.00000000E+00   # Re(ZUR(3,1))
  3  2     0.00000000E+00   # Re(ZUR(3,2))
  3  3     1.00000000E+00   # Re(ZUR(3,3))
Block FlexibleSUSYOutput
     0     0.00000000E+00   # HighScale
     1     2.00000000E+03   # SUSYScale
     2     9.11876000E+01   # LowScale
Block ALPHA
          -1.36337055E+00   # ArcSin(Pole(ZH(2,2)))
Block HMIX Q= 2.00000000E+03
     1     1.99999979E+03   # Re(Mu)
     2     4.77520537E+00   # vu/vd
     3     2.44022303E+02   # Sqrt(Sqr(vd) + Sqr(vu))
   101     8.02468438E+05   # Re(BMu)
   102     5.00169779E+01   # vd
   103     2.38841341E+02   # vu
Block ImHMIX Q= 2.00000000E+03
     1     9.99999897E+01   # Im(Mu)
   101    -8.04093338E+03   # Im(BMu)
Block Au Q= 2.00000000E+03
  1  1     4.00001141E+02   # Re(TYu(1,1)/Yu(1,1))
  2  2     4.00001141E+02   # Re(TYu(2,2)/Yu(2,2))
  3  3     4.00000303E+02   # Re(TYu(3,3)/Yu(3,3))
Block Ad Q= 2.00000000E+03
  1  1     9.99997099E+03   # Re(TYd(1,1)/Yd(1,1))
  2  2     9.99997099E+03   # Re(TYd(2,2)/Yd(2,2))
  3  3     9.99992474E+03   # Re(TYd(3,3)/Yd(3,3))
Block Ae Q= 2.00000000E+03
  1  1     1.00000160E+04   # Re(TYe(1,1)/Ye(1,1))
  2  2     1.00000160E+04   # Re(TYe(2,2)/Ye(2,2))
  3  3     1.00000157E+04   # Re(TYe(3,3)/Ye(3,3))
Block ImAu Q= 2.00000000E+03
  1  1     1.17109167E-04   # Im(TYu(1,1)/Yu(1,1))
  2  2     1.17109119E-04   # Im(TYu(2,2)/Yu(2,2))
  3  3     2.15829048E-06   # Im(TYu(3,3)/Yu(3,3))
Block ImAd Q= 2.00000000E+03
  1  1     1.18948580E-04   # Im(TYd(1,1)/Yd(1,1))
  2  2     1.18948555E-04   # Im(TYd(2,2)/Yd(2,2))
  3  3     2.62803507E-04   # Im(TYd(3,3)/Yd(3,3))
Block ImAe Q= 2.00000000E+03
  1  1    -8.85710713E-06   # Im(TYe(1,1)/Ye(1,1))
  2  2    -8.85638533E-06   # Im(TYe(2,2)/Ye(2,2))
  3  3    -8.65720136E-06   # Im(TYe(3,3)/Ye(3,3))
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
   FlexibleDecay_settings flexibledecay_settings;

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

   MSSMCPV_slha m(input);
   slha_io.fill(m);
   m.calculate_DRbar_masses();

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   MSSMCPV_decays decays_with_HO(m, qedqcd, physical_input, flexibledecay_settings);

   // scalar Higgs

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFdFd(&m, 1, 2, 2),
                              0.0021079024633380125, 6e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFuFu(&m, 1, 1, 1),
                              0.00010447088544366981, 2e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFeFe(&m, 1, 2, 2),
                              0.00022563024219708657, 3e-13);
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
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VGVG(&m, 1), 0.0002318886239533966, 2e-10);
   // h -> gamma gamma
   // without 2L QCD for squark
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 6.3284545616000571e-06, 4e-11);
   // with 2L QCD for squark
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 1), 5.3089300093115946e-06, 7e-10);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVZ(&m, 1), 5.1994372474535712e-07, 7e-11);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   MSSMCPV_decays decays_without_HO(m, qedqcd, physical_input, flexibledecay_settings);

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFdFd(&m, 1, 2, 2),
                              0.0012241960498053724, 2e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFuFu(&m, 1, 1, 1),
                              6.5917338328992629e-05, 2e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFeFe(&m, 1, 2, 2),
                              0.00022323499746669591, 3e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVP(&m, 1), 5.2382389517397406e-06, 7e-10);

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VGVG(&m, 1), 6.1972398675528688e-05, 2e-10);

   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVZ(&m, 1), 5.1994372474535712e-07, 7e-11);
}
