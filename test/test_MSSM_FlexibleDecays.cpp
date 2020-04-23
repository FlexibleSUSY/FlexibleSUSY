
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_FlexibleDecays

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_spectrum_generator.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSM_decays.hpp"
#include "MSSM_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MSSM_FlexibleDecays )
{

  char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.4.1
     5   MSSM
     9   4.14.3
Block MODSEL                 # Select model
    6   0                    # flavour violation
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
    9   80.404               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
# Block VCKMIN                 # CKM matrix input (Wolfenstein parameters)
#     1   2.272000000e-01      # lambda(MZ) SM DR-bar
#     2   8.180000000e-01      # A(MZ) SM DR-bar
#     3   2.210000000e-01      # rhobar(MZ) SM DR-bar
#     4   3.400000000e-01      # etabar(MZ) SM DR-bar
# Block UPMNSIN                # PMNS matrix input
#     1   5.837630000e-01      # theta_12
#     2   7.695840000e-01      # theta_23
#     3   1.549480000e-01      # theta_13
#     4   0.000000000e+00      # delta
#     5   0.000000000e+00      # alpha_1
#     6   0.000000000e+00      # alpha_2
Block MINPAR
     3     1.00000000E+01   # TanBeta
     4     1.00000000E+00   # SignMu
Block EXTPAR
     0     1.94220000E+16   # Qin
    21     1.56250000E+04   # mHd2IN
    22     1.56250000E+04   # mHu2IN
Block ADIN
  1  1     0.00000000E+00   # Adij(1,1)
  1  2     0.00000000E+00   # Adij(1,2)
  1  3     0.00000000E+00   # Adij(1,3)
  2  1     0.00000000E+00   # Adij(2,1)
  2  2     0.00000000E+00   # Adij(2,2)
  2  3     0.00000000E+00   # Adij(2,3)
  3  1     0.00000000E+00   # Adij(3,1)
  3  2     0.00000000E+00   # Adij(3,2)
  3  3     0.00000000E+00   # Adij(3,3)
Block AEIN
  1  1     0.00000000E+00   # Aeij(1,1)
  1  2     0.00000000E+00   # Aeij(1,2)
  1  3     0.00000000E+00   # Aeij(1,3)
  2  1     0.00000000E+00   # Aeij(2,1)
  2  2     0.00000000E+00   # Aeij(2,2)
  2  3     0.00000000E+00   # Aeij(2,3)
  3  1     0.00000000E+00   # Aeij(3,1)
  3  2     0.00000000E+00   # Aeij(3,2)
  3  3     0.00000000E+00   # Aeij(3,3)
Block AUIN
  1  1     0.00000000E+00   # Auij(1,1)
  1  2     0.00000000E+00   # Auij(1,2)
  1  3     0.00000000E+00   # Auij(1,3)
  2  1     0.00000000E+00   # Auij(2,1)
  2  2     0.00000000E+00   # Auij(2,2)
  2  3     0.00000000E+00   # Auij(2,3)
  3  1     0.00000000E+00   # Auij(3,1)
  3  2     0.00000000E+00   # Auij(3,2)
  3  3     0.00000000E+00   # Auij(3,3)
Block MSOFTIN
     1     5.00000000E+02   # MassBInput
     3     5.00000000E+02   # MassGInput
     2     5.00000000E+02   # MassWBInput
Block MSD2IN
  1  1     1.56250000E+04   # md2Input(1,1)
  1  2     0.00000000E+00   # md2Input(1,2)
  1  3     0.00000000E+00   # md2Input(1,3)
  2  1     0.00000000E+00   # md2Input(2,1)
  2  2     1.56250000E+04   # md2Input(2,2)
  2  3     0.00000000E+00   # md2Input(2,3)
  3  1     0.00000000E+00   # md2Input(3,1)
  3  2     0.00000000E+00   # md2Input(3,2)
  3  3     1.56250000E+04   # md2Input(3,3)
Block MSE2IN
  1  1     1.56250000E+04   # me2Input(1,1)
  1  2     0.00000000E+00   # me2Input(1,2)
  1  3     0.00000000E+00   # me2Input(1,3)
  2  1     0.00000000E+00   # me2Input(2,1)
  2  2     1.56250000E+04   # me2Input(2,2)
  2  3     0.00000000E+00   # me2Input(2,3)
  3  1     0.00000000E+00   # me2Input(3,1)
  3  2     0.00000000E+00   # me2Input(3,2)
  3  3     1.56250000E+04   # me2Input(3,3)
Block MSL2IN
  1  1     1.56250000E+04   # ml2Input(1,1)
  1  2     0.00000000E+00   # ml2Input(1,2)
  1  3     0.00000000E+00   # ml2Input(1,3)
  2  1     0.00000000E+00   # ml2Input(2,1)
  2  2     1.56250000E+04   # ml2Input(2,2)
  2  3     0.00000000E+00   # ml2Input(2,3)
  3  1     0.00000000E+00   # ml2Input(3,1)
  3  2     0.00000000E+00   # ml2Input(3,2)
  3  3     1.56250000E+04   # ml2Input(3,3)
Block MSQ2IN
  1  1     1.56250000E+04   # mq2Input(1,1)
  1  2     0.00000000E+00   # mq2Input(1,2)
  1  3     0.00000000E+00   # mq2Input(1,3)
  2  1     0.00000000E+00   # mq2Input(2,1)
  2  2     1.56250000E+04   # mq2Input(2,2)
  2  3     0.00000000E+00   # mq2Input(2,3)
  3  1     0.00000000E+00   # mq2Input(3,1)
  3  2     0.00000000E+00   # mq2Input(3,2)
  3  3     1.56250000E+04   # mq2Input(3,3)
Block MSU2IN
  1  1     1.56250000E+04   # mu2Input(1,1)
  1  2     0.00000000E+00   # mu2Input(1,2)
  1  3     0.00000000E+00   # mu2Input(1,3)
  2  1     0.00000000E+00   # mu2Input(2,1)
  2  2     1.56250000E+04   # mu2Input(2,2)
  2  3     0.00000000E+00   # mu2Input(2,3)
  3  1     0.00000000E+00   # mu2Input(3,1)
  3  2     0.00000000E+00   # mu2Input(3,2)
  3  3     1.56250000E+04   # mu2Input(3,3)
Block gauge Q= 8.61574711E+02
     1     3.59999781E-01   # g1 * 0.7745966692414834
     2     6.56898251E-01   # g2
     3     1.06298750E+00   # g3
Block Yu Q= 8.61574711E+02
  1  1     7.44114581E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.40276871E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.75768601E-01   # Yu(3,3)
Block Yd Q= 8.61574711E+02
  1  1     1.43200353E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     3.13533605E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     1.35247807E-01   # Yd(3,3)
Block Ye Q= 8.61574711E+02
  1  1     2.94098526E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     6.08102487E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     1.02276019E-01   # Ye(3,3)
Block Te Q= 8.61574711E+02
  1  1    -9.05425923E-03   # TYe(1,1)
  1  2     0.00000000E+00   # TYe(1,2)
  1  3     0.00000000E+00   # TYe(1,3)
  2  1     0.00000000E+00   # TYe(2,1)
  2  2    -1.87209650E+00   # TYe(2,2)
  2  3     0.00000000E+00   # TYe(2,3)
  3  1     0.00000000E+00   # TYe(3,1)
  3  2     0.00000000E+00   # TYe(3,2)
  3  3    -3.13103118E+01   # TYe(3,3)
Block Td Q= 8.61574711E+02
  1  1    -1.99663480E-01   # TYd(1,1)
  1  2     0.00000000E+00   # TYd(1,2)
  1  3     0.00000000E+00   # TYd(1,3)
  2  1     0.00000000E+00   # TYd(2,1)
  2  2    -4.37156706E+00   # TYd(2,2)
  2  3     0.00000000E+00   # TYd(2,3)
  3  1     0.00000000E+00   # TYd(3,1)
  3  2     0.00000000E+00   # TYd(3,2)
  3  3    -1.75898126E+02   # TYd(3,3)
Block Tu Q= 8.61574711E+02
  1  1    -8.41105668E-03   # TYu(1,1)
  1  2     0.00000000E+00   # TYu(1,2)
  1  3     0.00000000E+00   # TYu(1,3)
  2  1     0.00000000E+00   # TYu(2,1)
  2  2    -3.84628253E+00   # TYu(2,2)
  2  3     0.00000000E+00   # TYu(2,3)
  3  1     0.00000000E+00   # TYu(3,1)
  3  2     0.00000000E+00   # TYu(3,2)
  3  3    -7.55562692E+02   # TYu(3,3)
Block MSQ2 Q= 8.61574711E+02
  1  1     1.01698410E+06   # mq2(1,1)
  1  2     0.00000000E+00   # mq2(1,2)
  1  3     0.00000000E+00   # mq2(1,3)
  2  1     0.00000000E+00   # mq2(2,1)
  2  2     1.01697877E+06   # mq2(2,2)
  2  3     0.00000000E+00   # mq2(2,3)
  3  1     0.00000000E+00   # mq2(3,1)
  3  2     0.00000000E+00   # mq2(3,2)
  3  3     8.60663456E+05   # mq2(3,3)
Block MSE2 Q= 8.61574711E+02
  1  1     4.92190886E+04   # me2(1,1)
  1  2     0.00000000E+00   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     0.00000000E+00   # me2(2,1)
  2  2     4.92136724E+04   # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     4.76891146E+04   # me2(3,3)
Block MSL2 Q= 8.61574711E+02
  1  1     1.28303328E+05   # ml2(1,1)
  1  2     0.00000000E+00   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     0.00000000E+00   # ml2(2,1)
  2  2     1.28300676E+05   # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     1.27554323E+05   # ml2(3,3)
Block MSU2 Q= 8.61574711E+02
  1  1     9.37379062E+05   # mu2(1,1)
  1  2     0.00000000E+00   # mu2(1,2)
  1  3     0.00000000E+00   # mu2(1,3)
  2  1     0.00000000E+00   # mu2(2,1)
  2  2     9.37373619E+05   # mu2(2,2)
  2  3     0.00000000E+00   # mu2(2,3)
  3  1     0.00000000E+00   # mu2(3,1)
  3  2     0.00000000E+00   # mu2(3,2)
  3  3     6.27037222E+05   # mu2(3,3)
Block MSD2 Q= 8.61574711E+02
  1  1     9.28209501E+05   # md2(1,1)
  1  2     0.00000000E+00   # md2(1,2)
  1  3     0.00000000E+00   # md2(1,3)
  2  1     0.00000000E+00   # md2(2,1)
  2  2     9.28204149E+05   # md2(2,2)
  2  3     0.00000000E+00   # md2(2,3)
  3  1     0.00000000E+00   # md2(3,1)
  3  2     0.00000000E+00   # md2(3,2)
  3  3     9.18700843E+05   # md2(3,3)
Block Phases Q= 8.61574711E+02
     1     1.00000000E+00   # Re(PhaseGlu)
Block IMPhases Q= 8.61574711E+02
     1     0.00000000E+00   # Im(PhaseGlu)
Block MASS
   1000021     1.14482639E+03   # Glu
        24     8.22179853E+01   # VWm
   1000024     3.80849556E+02   # Cha(1)
   1000037     6.51784974E+02   # Cha(2)
        25     1.15461248E+02   # hh(1)
        35     7.22610281E+02   # hh(2)
        37     7.27129000E+02   # Hpm(2)
        36     7.22328661E+02   # Ah(2)
   1000012     3.55652239E+02   # Sv(1)
   1000014     3.56861268E+02   # Sv(2)
   1000016     3.56865558E+02   # Sv(3)
   1000022     2.07672439E+02   # Chi(1)
   1000023     3.80829020E+02   # Chi(2)
   1000025    -6.38064535E+02   # Chi(3)
   1000035     6.51307635E+02   # Chi(4)
   1000001     9.54750500E+02   # Sd(1)
   1000003     9.94898036E+02   # Sd(2)
   1000005     9.98047186E+02   # Sd(3)
   2000001     9.98051060E+02   # Sd(4)
   2000003     1.04512659E+03   # Sd(5)
   2000005     1.04512853E+03   # Sd(6)
   1000011     2.22269063E+02   # Se(1)
   1000013     2.29495276E+02   # Se(2)
   1000015     2.29520911E+02   # Se(3)
   2000011     3.65774038E+02   # Se(4)
   2000013     3.65777881E+02   # Se(5)
   2000015     3.66795565E+02   # Se(6)
   1000002     7.90438045E+02   # Su(1)
   1000004     9.98575595E+02   # Su(2)
   1000006     1.00154427E+03   # Su(3)
   2000002     1.00233322E+03   # Su(4)
   2000004     1.04221767E+03   # Su(5)
   2000006     1.04221811E+03   # Su(6)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.30386450E-04   # Fe(1)
        13     1.07519240E-01   # Fe(2)
        15     1.78907515E+00   # Fe(3)
         1     4.61547151E-03   # Fd(1)
         3     9.12933235E-02   # Fd(2)
         5     3.39876655E+00   # Fd(3)
         2     2.32347760E-03   # Fu(1)
         4     8.54837592E-01   # Fu(2)
         6     1.73817518E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.11338975E+01   # VZ
Block UMIX
  1  1     9.60383246E-01   # Re(UM(1,1))
  1  2    -2.78682654E-01   # Re(UM(1,2))
  2  1     2.78682654E-01   # Re(UM(2,1))
  2  2     9.60383246E-01   # Re(UM(2,2))
Block VMIX
  1  1     9.82927175E-01   # Re(UP(1,1))
  1  2    -1.83995021E-01   # Re(UP(1,2))
  2  1     1.83995021E-01   # Re(UP(2,1))
  2  2     9.82927175E-01   # Re(UP(2,2))
Block PSEUDOSCALARMIX
  1  1    -9.90160581E-02   # ZA(1,1)
  1  2     9.95085836E-01   # ZA(1,2)
  2  1     9.95085836E-01   # ZA(2,1)
  2  2     9.90160581E-02   # ZA(2,2)
Block DSQMIX
  1  1    -0.00000000E+00   # ZD(1,1)
  1  2    -0.00000000E+00   # ZD(1,2)
  1  3    -9.74286165E-01   # ZD(1,3)
  1  4    -0.00000000E+00   # ZD(1,4)
  1  5    -0.00000000E+00   # ZD(1,5)
  1  6    -2.25314157E-01   # ZD(1,6)
  2  1     0.00000000E+00   # ZD(2,1)
  2  2     0.00000000E+00   # ZD(2,2)
  2  3     2.25314157E-01   # ZD(2,3)
  2  4     0.00000000E+00   # ZD(2,4)
  2  5     0.00000000E+00   # ZD(2,5)
  2  6    -9.74286165E-01   # ZD(2,6)
  3  1    -0.00000000E+00   # ZD(3,1)
  3  2    -4.48067024E-03   # ZD(3,2)
  3  3    -0.00000000E+00   # ZD(3,3)
  3  4    -0.00000000E+00   # ZD(3,4)
  3  5    -9.99989962E-01   # ZD(3,5)
  3  6    -0.00000000E+00   # ZD(3,6)
  4  1     2.04652197E-04   # ZD(4,1)
  4  2     0.00000000E+00   # ZD(4,2)
  4  3     0.00000000E+00   # ZD(4,3)
  4  4     9.99999979E-01   # ZD(4,4)
  4  5     0.00000000E+00   # ZD(4,5)
  4  6     0.00000000E+00   # ZD(4,6)
  5  1     0.00000000E+00   # ZD(5,1)
  5  2    -9.99989962E-01   # ZD(5,2)
  5  3     0.00000000E+00   # ZD(5,3)
  5  4     0.00000000E+00   # ZD(5,4)
  5  5     4.48067024E-03   # ZD(5,5)
  5  6     0.00000000E+00   # ZD(5,6)
  6  1     9.99999979E-01   # ZD(6,1)
  6  2     0.00000000E+00   # ZD(6,2)
  6  3     0.00000000E+00   # ZD(6,3)
  6  4    -2.04652197E-04   # ZD(6,4)
  6  5     0.00000000E+00   # ZD(6,5)
  6  6     0.00000000E+00   # ZD(6,6)
Block SELMIX
  1  1     0.00000000E+00   # ZE(1,1)
  1  2     2.71107806E-17   # ZE(1,2)
  1  3     1.37829507E-01   # ZE(1,3)
  1  4     0.00000000E+00   # ZE(1,4)
  1  5     3.31164198E-16   # ZE(1,5)
  1  6     9.90455969E-01   # ZE(1,6)
  2  1     0.00000000E+00   # ZE(2,1)
  2  2    -8.52595046E-03   # ZE(2,2)
  2  3     3.05760718E-17   # ZE(2,3)
  2  4     0.00000000E+00   # ZE(2,4)
  2  5    -9.99963653E-01   # ZE(2,5)
  2  6     3.38375884E-16   # ZE(2,6)
  3  1     4.12403764E-05   # ZE(3,1)
  3  2     0.00000000E+00   # ZE(3,2)
  3  3     0.00000000E+00   # ZE(3,3)
  3  4     9.99999999E-01   # ZE(3,4)
  3  5     0.00000000E+00   # ZE(3,5)
  3  6     0.00000000E+00   # ZE(3,6)
  4  1     9.99999999E-01   # ZE(4,1)
  4  2     0.00000000E+00   # ZE(4,2)
  4  3     0.00000000E+00   # ZE(4,3)
  4  4    -4.12403764E-05   # ZE(4,4)
  4  5     0.00000000E+00   # ZE(4,5)
  4  6     0.00000000E+00   # ZE(4,6)
  5  1     0.00000000E+00   # ZE(5,1)
  5  2     9.99963653E-01   # ZE(5,2)
  5  3     1.87938838E-14   # ZE(5,3)
  5  4     0.00000000E+00   # ZE(5,4)
  5  5    -8.52595046E-03   # ZE(5,5)
  5  6    -2.63976402E-15   # ZE(5,6)
  6  1     0.00000000E+00   # ZE(6,1)
  6  2    -1.89778014E-14   # ZE(6,2)
  6  3     9.90455969E-01   # ZE(6,3)
  6  4     0.00000000E+00   # ZE(6,4)
  6  5     1.43100164E-16   # ZE(6,5)
  6  6    -1.37829507E-01   # ZE(6,6)
Block SCALARMIX
  1  1     1.06626992E-01   # ZH(1,1)
  1  2     9.94299092E-01   # ZH(1,2)
  2  1     9.94299092E-01   # ZH(2,1)
  2  2    -1.06626992E-01   # ZH(2,2)
Block NMIX
  1  1    -9.95937884E-01   # Re(ZN(1,1))
  1  2     1.82403903E-02   # Re(ZN(1,2))
  1  3    -8.13208590E-02   # Re(ZN(1,3))
  1  4     3.40871879E-02   # Re(ZN(1,4))
  2  1     3.80941901E-02   # Re(ZN(2,1))
  2  2     9.71669835E-01   # Re(ZN(2,2))
  2  3    -1.94844064E-01   # Re(ZN(2,3))
  2  4     1.28227750E-01   # Re(ZN(2,4))
  3  1     3.23009163E-02   # Re(ZN(3,1))
  3  2    -4.88516996E-02   # Re(ZN(3,2))
  3  3    -7.03454555E-01   # Re(ZN(3,3))
  3  4    -7.08323268E-01   # Re(ZN(3,4))
  4  1     7.49213842E-02   # Re(ZN(4,1))
  4  2    -2.30517965E-01   # Re(ZN(4,2))
  4  3    -6.78656318E-01   # Re(ZN(4,3))
  4  4     6.93306466E-01   # Re(ZN(4,4))
Block CHARGEMIX
  1  1    -9.96359870E-02   # ZP(1,1)
  1  2     9.95023955E-01   # ZP(1,2)
  2  1     9.95023955E-01   # ZP(2,1)
  2  2     9.96359870E-02   # ZP(2,2)
Block USQMIX
  1  1     0.00000000E+00   # ZU(1,1)
  1  2     0.00000000E+00   # ZU(1,2)
  1  3     4.19318558E-01   # ZU(1,3)
  1  4     0.00000000E+00   # ZU(1,4)
  1  5     0.00000000E+00   # ZU(1,5)
  1  6     9.07839164E-01   # ZU(1,6)
  2  1    -0.00000000E+00   # ZU(2,1)
  2  2    -8.95797956E-03   # ZU(2,2)
  2  3    -0.00000000E+00   # ZU(2,3)
  2  4    -0.00000000E+00   # ZU(2,4)
  2  5    -9.99959876E-01   # ZU(2,5)
  2  6    -0.00000000E+00   # ZU(2,6)
  3  1     1.95917144E-05   # ZU(3,1)
  3  2     0.00000000E+00   # ZU(3,2)
  3  3     0.00000000E+00   # ZU(3,3)
  3  4     1.00000000E+00   # ZU(3,4)
  3  5     0.00000000E+00   # ZU(3,5)
  3  6     0.00000000E+00   # ZU(3,6)
  4  1     0.00000000E+00   # ZU(4,1)
  4  2     0.00000000E+00   # ZU(4,2)
  4  3     9.07839164E-01   # ZU(4,3)
  4  4     0.00000000E+00   # ZU(4,4)
  4  5     0.00000000E+00   # ZU(4,5)
  4  6    -4.19318558E-01   # ZU(4,6)
  5  1     1.00000000E+00   # ZU(5,1)
  5  2     0.00000000E+00   # ZU(5,2)
  5  3     0.00000000E+00   # ZU(5,3)
  5  4    -1.95917144E-05   # ZU(5,4)
  5  5     0.00000000E+00   # ZU(5,5)
  5  6     0.00000000E+00   # ZU(5,6)
  6  1     0.00000000E+00   # ZU(6,1)
  6  2    -9.99959876E-01   # ZU(6,2)
  6  3     0.00000000E+00   # ZU(6,3)
  6  4     0.00000000E+00   # ZU(6,4)
  6  5     8.95797956E-03   # ZU(6,5)
  6  6     0.00000000E+00   # ZU(6,6)
Block SNUMIX
  1  1     0.00000000E+00   # ZV(1,1)
  1  2     0.00000000E+00   # ZV(1,2)
  1  3     1.00000000E+00   # ZV(1,3)
  2  1     0.00000000E+00   # ZV(2,1)
  2  2     1.00000000E+00   # ZV(2,2)
  2  3     0.00000000E+00   # ZV(2,3)
  3  1     1.00000000E+00   # ZV(3,1)
  3  2     0.00000000E+00   # ZV(3,2)
  3  3     0.00000000E+00   # ZV(3,3)
Block UELMIX
  1  1     1.00000000E+00   # Re(ZEL(1,1))
  1  2     0.00000000E+00   # Re(ZEL(1,2))
  1  3     0.00000000E+00   # Re(ZEL(1,3))
  2  1     0.00000000E+00   # Re(ZEL(2,1))
  2  2     1.00000000E+00   # Re(ZEL(2,2))
  2  3     0.00000000E+00   # Re(ZEL(2,3))
  3  1     0.00000000E+00   # Re(ZEL(3,1))
  3  2     0.00000000E+00   # Re(ZEL(3,2))
  3  3     1.00000000E+00   # Re(ZEL(3,3))
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
  1  1     1.00000000E+00   # Re(ZDL(1,1))
  1  2     0.00000000E+00   # Re(ZDL(1,2))
  1  3     0.00000000E+00   # Re(ZDL(1,3))
  2  1     0.00000000E+00   # Re(ZDL(2,1))
  2  2     1.00000000E+00   # Re(ZDL(2,2))
  2  3     0.00000000E+00   # Re(ZDL(2,3))
  3  1     0.00000000E+00   # Re(ZDL(3,1))
  3  2     0.00000000E+00   # Re(ZDL(3,2))
  3  3     1.00000000E+00   # Re(ZDL(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(ZDR(1,1))
  1  2     0.00000000E+00   # Re(ZDR(1,2))
  1  3     0.00000000E+00   # Re(ZDR(1,3))
  2  1     0.00000000E+00   # Re(ZDR(2,1))
  2  2     1.00000000E+00   # Re(ZDR(2,2))
  2  3     0.00000000E+00   # Re(ZDR(2,3))
  3  1     0.00000000E+00   # Re(ZDR(3,1))
  3  2     0.00000000E+00   # Re(ZDR(3,2))
  3  3     1.00000000E+00   # Re(ZDR(3,3))
Block UULMIX
  1  1     1.00000000E+00   # Re(ZUL(1,1))
  1  2     0.00000000E+00   # Re(ZUL(1,2))
  1  3     0.00000000E+00   # Re(ZUL(1,3))
  2  1     0.00000000E+00   # Re(ZUL(2,1))
  2  2     1.00000000E+00   # Re(ZUL(2,2))
  2  3     0.00000000E+00   # Re(ZUL(2,3))
  3  1     0.00000000E+00   # Re(ZUL(3,1))
  3  2     0.00000000E+00   # Re(ZUL(3,2))
  3  3     1.00000000E+00   # Re(ZUL(3,3))
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
     0     1.94220000E+16   # HighScale
     1     8.61574711E+02   # SUSYScale
     2     9.11876000E+01   # LowScale
Block ALPHA
          -1.06830078E-01   # ArcSin(Pole(ZH(2,2)))
Block HMIX Q= 8.61574711E+02
     1     6.32596004E+02   # Mu
     2     9.66521530E+00   # vu/vd
     3     2.40393898E+02   # Sqrt(Sqr(vd) + Sqr(vu))
     4     5.41191164E+05   # Sqr(MAh(2))
   101     5.54006492E+04   # BMu
   102     2.47400032E+01   # vd
   103     2.39117457E+02   # vu
Block Au Q= 8.61574711E+02
  1  1    -1.13034429E+03   # TYu(1,1)/Yu(1,1)
  2  2    -1.13033910E+03   # TYu(2,2)/Yu(2,2)
  3  3    -8.62742386E+02   # TYu(3,3)/Yu(3,3)
Block Ad Q= 8.61574711E+02
  1  1    -1.39429461E+03   # TYd(1,1)/Yd(1,1)
  2  2    -1.39428979E+03   # TYd(2,2)/Yd(2,2)
  3  3    -1.30056176E+03   # TYd(3,3)/Yd(3,3)
Block Ae Q= 8.61574711E+02
  1  1    -3.07864828E+02   # TYe(1,1)/Ye(1,1)
  2  2    -3.07858714E+02   # TYe(2,2)/Ye(2,2)
  3  3    -3.06135419E+02   # TYe(3,3)/Ye(3,3)
Block MSOFT Q= 8.61574711E+02
     1     2.12674076E+02   # MassB
     2     3.82588792E+02   # MassWB
     3     1.11445312E+03   # MassG
    21     1.12097875E+05   # mHd2
    22    -3.87890853E+05   # mHu2
    31     3.58194540E+02   # SignedAbsSqrt(ml2(1,1))
    32     3.58190838E+02   # SignedAbsSqrt(ml2(2,2))
    33     3.57147481E+02   # SignedAbsSqrt(ml2(3,3))
    34     2.21853755E+02   # SignedAbsSqrt(me2(1,1))
    35     2.21841548E+02   # SignedAbsSqrt(me2(2,2))
    36     2.18378375E+02   # SignedAbsSqrt(me2(3,3))
    41     1.00845629E+03   # SignedAbsSqrt(mq2(1,1))
    42     1.00845366E+03   # SignedAbsSqrt(mq2(2,2))
    43     9.27719492E+02   # SignedAbsSqrt(mq2(3,3))
    44     9.68183382E+02   # SignedAbsSqrt(mu2(1,1))
    45     9.68180571E+02   # SignedAbsSqrt(mu2(2,2))
    46     7.91856819E+02   # SignedAbsSqrt(mu2(3,3))
    47     9.63436298E+02   # SignedAbsSqrt(md2(1,1))
    48     9.63433521E+02   # SignedAbsSqrt(md2(2,2))
    49     9.58488833E+02   # SignedAbsSqrt(md2(3,3))
)";

   std::stringstream istr(slha_input);

   MSSM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   // extract the input parameters
   softsusy::QedQcd qedqcd;
   MSSM_input_parameters input;
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::calculate_sm_masses, 0);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 0);

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   MSSM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   MSSM_slha<MSSM<Two_scale>> m = std::get<0>(spectrum_generator.get_models_slha());

   MSSM_decays decays(m, qedqcd, input, HigherOrderSMCorrections::enable);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0016298093732297074, 6e-14);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.0002600448772469762, 6e-14);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_conjVWmVWm(&m, 0),
                              0.00022767000236269188, 8e-13);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVZ(&m, 0),
                              1.6116438957094875e-05, 2e-12);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(&m, 0), 0.00015442075031863592, 3e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(&m, 0), 6.5333738573774704e-06, 4e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVZ(&m, 0), 3.2717635916980618e-06, 1e-12);

}
