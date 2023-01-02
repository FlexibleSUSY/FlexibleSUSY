
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_FlexibleDecay

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_spectrum_generator.hpp"
#include "MSSM_two_scale_model.hpp"
#include "decays/MSSM_decays.hpp"
#include "MSSM_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MSSM_FlexibleDecay )
{

  char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.6.1
     5   MSSM
     9   4.14.5
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
   31   -1                   # 0(Softsusy),1(Collier),2(Looptools),3(fflite)
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
Block gauge Q= 8.64566212E+02
     1     3.62397260E-01   # g1 * 0.7745966692414834
     2     6.43021599E-01   # g2
     3     1.06289850E+00   # g3
Block Yu Q= 8.64566212E+02
  1  1     7.32865780E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.35132890E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.61499146E-01   # Yu(3,3)
Block Yd Q= 8.64566212E+02
  1  1     1.41176492E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     3.09102407E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     1.33305769E-01   # Yd(3,3)
Block Ye Q= 8.64566212E+02
  1  1     2.89548332E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     5.98694134E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     1.00693901E-01   # Ye(3,3)
Block Te Q= 8.64566212E+02
  1  1    -8.66772958E-03   # TYe(1,1)
  1  2     0.00000000E+00   # TYe(1,2)
  1  3     0.00000000E+00   # TYe(1,3)
  2  1     0.00000000E+00   # TYe(2,1)
  2  2    -1.79217692E+00   # TYe(2,2)
  2  3     0.00000000E+00   # TYe(2,3)
  3  1     0.00000000E+00   # TYe(3,1)
  3  2     0.00000000E+00   # TYe(3,2)
  3  3    -2.99779990E+01   # TYe(3,3)
Block Td Q= 8.64566212E+02
  1  1    -1.95555678E-01   # TYd(1,1)
  1  2     0.00000000E+00   # TYd(1,2)
  1  3     0.00000000E+00   # TYd(1,3)
  2  1     0.00000000E+00   # TYd(2,1)
  2  2    -4.28162839E+00   # TYd(2,2)
  2  3     0.00000000E+00   # TYd(2,3)
  3  1     0.00000000E+00   # TYd(3,1)
  3  2     0.00000000E+00   # TYd(3,2)
  3  3    -1.72594209E+02   # TYd(3,3)
Block Tu Q= 8.64566212E+02
  1  1    -8.29130231E-03   # TYu(1,1)
  1  2     0.00000000E+00   # TYu(1,2)
  1  3     0.00000000E+00   # TYu(1,3)
  2  1     0.00000000E+00   # TYu(2,1)
  2  2    -3.79152071E+00   # TYu(2,2)
  2  3     0.00000000E+00   # TYu(2,3)
  3  1     0.00000000E+00   # TYu(3,1)
  3  2     0.00000000E+00   # TYu(3,2)
  3  3    -7.52334904E+02   # TYu(3,3)
Block MSQ2 Q= 8.64566212E+02
  1  1     1.01404039E+06   # mq2(1,1)
  1  2     0.00000000E+00   # mq2(1,2)
  1  3     0.00000000E+00   # mq2(1,3)
  2  1     0.00000000E+00   # mq2(2,1)
  2  2     1.01403524E+06   # mq2(2,2)
  2  3     0.00000000E+00   # mq2(2,3)
  3  1     0.00000000E+00   # mq2(3,1)
  3  2     0.00000000E+00   # mq2(3,2)
  3  3     8.61830295E+05   # mq2(3,3)
Block MSE2 Q= 8.64566212E+02
  1  1     4.93996314E+04   # me2(1,1)
  1  2     0.00000000E+00   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     0.00000000E+00   # me2(2,1)
  2  2     4.93944907E+04   # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     4.79473518E+04   # me2(3,3)
Block MSL2 Q= 8.64566212E+02
  1  1     1.24998047E+05   # ml2(1,1)
  1  2     0.00000000E+00   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     0.00000000E+00   # ml2(2,1)
  2  2     1.24995528E+05   # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     1.24286462E+05   # ml2(3,3)
Block MSU2 Q= 8.64566212E+02
  1  1     9.37329592E+05   # mu2(1,1)
  1  2     0.00000000E+00   # mu2(1,2)
  1  3     0.00000000E+00   # mu2(1,3)
  2  1     0.00000000E+00   # mu2(2,1)
  2  2     9.37324316E+05   # mu2(2,2)
  2  3     0.00000000E+00   # mu2(2,3)
  3  1     0.00000000E+00   # mu2(3,1)
  3  2     0.00000000E+00   # mu2(3,2)
  3  3     6.35411190E+05   # mu2(3,3)
Block MSD2 Q= 8.64566212E+02
  1  1     9.28110603E+05   # md2(1,1)
  1  2     0.00000000E+00   # md2(1,2)
  1  3     0.00000000E+00   # md2(1,3)
  2  1     0.00000000E+00   # md2(2,1)
  2  2     9.28105434E+05   # md2(2,2)
  2  3     0.00000000E+00   # md2(2,3)
  3  1     0.00000000E+00   # md2(3,1)
  3  2     0.00000000E+00   # md2(3,2)
  3  3     9.18917534E+05   # md2(3,3)
Block Phases Q= 8.64566212E+02
     1     1.00000000E+00   # Re(PhaseGlu)
Block IMPhases Q= 8.64566212E+02
     1     0.00000000E+00   # Im(PhaseGlu)
Block MASS
   1000021     1.14507765E+03   # Glu
        24     8.03633615E+01   # VWm
   1000024     3.84946652E+02   # Cha(1)
   1000037     6.43567570E+02   # Cha(2)
        25     1.14844829E+02   # hh(1)
        35     7.12647702E+02   # hh(2)
        37     7.17159962E+02   # Hpm(2)
        36     7.12376885E+02   # Ah(2)
   1000012     3.50946348E+02   # Sv(1)
   1000014     3.52110103E+02   # Sv(2)
   1000016     3.52114232E+02   # Sv(3)
   1000022     2.03765032E+02   # Chi(1)
   1000023     3.84941967E+02   # Chi(2)
   1000025    -6.29278072E+02   # Chi(3)
   1000035     6.43254469E+02   # Chi(4)
   1000001     9.55616899E+02   # Sd(1)
   1000003     9.95172642E+02   # Sd(2)
   1000005     9.98169155E+02   # Sd(3)
   2000001     9.98172941E+02   # Sd(4)
   2000003     1.04379200E+03   # Sd(5)
   2000005     1.04379385E+03   # Sd(6)
   1000011     2.22944690E+02   # Se(1)
   1000013     2.30022228E+02   # Se(2)
   1000015     2.30047372E+02   # Se(3)
   2000011     3.61034664E+02   # Se(4)
   2000013     3.61038910E+02   # Se(5)
   2000015     3.62164628E+02   # Se(6)
   1000002     7.94247721E+02   # Su(1)
   1000004     1.00025984E+03   # Su(2)
   1000006     1.00165304E+03   # Su(3)
   2000002     1.00272702E+03   # Su(4)
   2000004     1.04091773E+03   # Su(5)
   2000006     1.04091837E+03   # Su(6)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.30293274E-04   # Fe(1)
        13     1.07502302E-01   # Fe(2)
        15     1.78881716E+00   # Fe(3)
         1     4.61593304E-03   # Fd(1)
         3     9.13029979E-02   # Fd(2)
         5     3.39894600E+00   # Fd(3)
         2     2.32430773E-03   # Fu(1)
         4     8.55208202E-01   # Fu(2)
         6     1.73822611E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.11298854E+01   # VZ
Block UMIX
  1  1     9.57429771E-01   # Re(UM(1,1))
  1  2    -2.88666301E-01   # Re(UM(1,2))
  2  1     2.88666301E-01   # Re(UM(2,1))
  2  2     9.57429771E-01   # Re(UM(2,2))
Block VMIX
  1  1     9.80920736E-01   # Re(UP(1,1))
  1  2    -1.94408100E-01   # Re(UP(1,2))
  2  1     1.94408100E-01   # Re(UP(2,1))
  2  2     9.80920736E-01   # Re(UP(2,2))
Block PSEUDOSCALARMIX
  1  1    -9.89804863E-02   # ZA(1,1)
  1  2     9.95089375E-01   # ZA(1,2)
  2  1     9.95089375E-01   # ZA(2,1)
  2  2     9.89804863E-02   # ZA(2,2)
Block DSQMIX
  1  1     0.00000000E+00   # ZD(1,1)
  1  2    -7.29213802E-16   # ZD(1,2)
  1  3     9.74106342E-01   # ZD(1,3)
  1  4     0.00000000E+00   # ZD(1,4)
  1  5    -7.61221493E-16   # ZD(1,5)
  1  6     2.26090322E-01   # ZD(1,6)
  2  1     0.00000000E+00   # ZD(2,1)
  2  2     3.07818909E-16   # ZD(2,2)
  2  3    -2.26090322E-01   # ZD(2,3)
  2  4     0.00000000E+00   # ZD(2,4)
  2  5     2.40876201E-15   # ZD(2,5)
  2  6     9.74106342E-01   # ZD(2,6)
  3  1     0.00000000E+00   # ZD(3,1)
  3  2    -4.57005211E-03   # ZD(3,2)
  3  3    -1.26686878E-15   # ZD(3,3)
  3  4     0.00000000E+00   # ZD(3,4)
  3  5    -9.99989557E-01   # ZD(3,5)
  3  6     2.16934167E-15   # ZD(3,6)
  4  1     2.08734887E-04   # ZD(4,1)
  4  2     0.00000000E+00   # ZD(4,2)
  4  3     0.00000000E+00   # ZD(4,3)
  4  4     9.99999978E-01   # ZD(4,4)
  4  5     0.00000000E+00   # ZD(4,5)
  4  6     0.00000000E+00   # ZD(4,6)
  5  1     0.00000000E+00   # ZD(5,1)
  5  2     9.99989557E-01   # ZD(5,2)
  5  3     7.74145093E-16   # ZD(5,3)
  5  4     0.00000000E+00   # ZD(5,4)
  5  5    -4.57005211E-03   # ZD(5,5)
  5  6    -1.25067470E-16   # ZD(5,6)
  6  1     9.99999978E-01   # ZD(6,1)
  6  2     0.00000000E+00   # ZD(6,2)
  6  3     0.00000000E+00   # ZD(6,3)
  6  4    -2.08734887E-04   # ZD(6,4)
  6  5     0.00000000E+00   # ZD(6,5)
  6  6     0.00000000E+00   # ZD(6,6)
Block SELMIX
  1  1     0.00000000E+00   # ZE(1,1)
  1  2     0.00000000E+00   # ZE(1,2)
  1  3     1.42036469E-01   # ZE(1,3)
  1  4     0.00000000E+00   # ZE(1,4)
  1  5     0.00000000E+00   # ZE(1,5)
  1  6     9.89861425E-01   # ZE(1,6)
  2  1     0.00000000E+00   # ZE(2,1)
  2  2     8.80185979E-03   # ZE(2,2)
  2  3     0.00000000E+00   # ZE(2,3)
  2  4     0.00000000E+00   # ZE(2,4)
  2  5     9.99961263E-01   # ZE(2,5)
  2  6     0.00000000E+00   # ZE(2,6)
  3  1     4.25752572E-05   # ZE(3,1)
  3  2     0.00000000E+00   # ZE(3,2)
  3  3     0.00000000E+00   # ZE(3,3)
  3  4     9.99999999E-01   # ZE(3,4)
  3  5     0.00000000E+00   # ZE(3,5)
  3  6     0.00000000E+00   # ZE(3,6)
  4  1     9.99999999E-01   # ZE(4,1)
  4  2     0.00000000E+00   # ZE(4,2)
  4  3     0.00000000E+00   # ZE(4,3)
  4  4    -4.25752572E-05   # ZE(4,4)
  4  5     0.00000000E+00   # ZE(4,5)
  4  6     0.00000000E+00   # ZE(4,6)
  5  1     0.00000000E+00   # ZE(5,1)
  5  2     9.99961263E-01   # ZE(5,2)
  5  3     0.00000000E+00   # ZE(5,3)
  5  4     0.00000000E+00   # ZE(5,4)
  5  5    -8.80185979E-03   # ZE(5,5)
  5  6     0.00000000E+00   # ZE(5,6)
  6  1     0.00000000E+00   # ZE(6,1)
  6  2     0.00000000E+00   # ZE(6,2)
  6  3     9.89861425E-01   # ZE(6,3)
  6  4     0.00000000E+00   # ZE(6,4)
  6  5     0.00000000E+00   # ZE(6,5)
  6  6    -1.42036469E-01   # ZE(6,6)
Block SCALARMIX
  1  1     1.06631109E-01   # ZH(1,1)
  1  2     9.94298651E-01   # ZH(1,2)
  2  1     9.94298651E-01   # ZH(2,1)
  2  2    -1.06631109E-01   # ZH(2,2)
Block NMIX
  1  1    -9.95692479E-01   # Re(ZN(1,1))
  1  2     1.78138151E-02   # Re(ZN(1,2))
  1  3    -8.39671875E-02   # Re(ZN(1,3))
  1  4     3.50523347E-02   # Re(ZN(1,4))
  2  1     3.91216656E-02   # Re(ZN(2,1))
  2  2     9.69238435E-01   # Re(ZN(2,2))
  2  3    -2.01735414E-01   # Re(ZN(2,3))
  2  4     1.35459125E-01   # Re(ZN(2,4))
  3  1     3.35071390E-02   # Re(ZN(3,1))
  3  2    -4.87542661E-02   # Re(ZN(3,2))
  3  3    -7.03377609E-01   # Re(ZN(3,3))
  3  4    -7.08350360E-01   # Re(ZN(3,4))
  4  1     7.70925008E-02   # Re(ZN(4,1))
  4  2    -2.40587917E-01   # Re(ZN(4,2))
  4  3    -6.76396535E-01   # Re(ZN(4,3))
  4  4     6.91853978E-01   # Re(ZN(4,4))
Block CHARGEMIX
  1  1    -9.95719443E-02   # ZP(1,1)
  1  2     9.95030365E-01   # ZP(1,2)
  2  1     9.95030365E-01   # ZP(2,1)
  2  2     9.95719443E-02   # ZP(2,2)
Block USQMIX
  1  1     0.00000000E+00   # ZU(1,1)
  1  2     0.00000000E+00   # ZU(1,2)
  1  3     4.27858606E-01   # ZU(1,3)
  1  4     0.00000000E+00   # ZU(1,4)
  1  5     0.00000000E+00   # ZU(1,5)
  1  6     9.03845680E-01   # ZU(1,6)
  2  1     0.00000000E+00   # ZU(2,1)
  2  2    -9.27975914E-03   # ZU(2,2)
  2  3     0.00000000E+00   # ZU(2,3)
  2  4     0.00000000E+00   # ZU(2,4)
  2  5    -9.99956942E-01   # ZU(2,5)
  2  6     0.00000000E+00   # ZU(2,6)
  3  1     2.02956466E-05   # ZU(3,1)
  3  2     0.00000000E+00   # ZU(3,2)
  3  3     0.00000000E+00   # ZU(3,3)
  3  4     1.00000000E+00   # ZU(3,4)
  3  5     0.00000000E+00   # ZU(3,5)
  3  6     0.00000000E+00   # ZU(3,6)
  4  1     0.00000000E+00   # ZU(4,1)
  4  2     0.00000000E+00   # ZU(4,2)
  4  3     9.03845680E-01   # ZU(4,3)
  4  4     0.00000000E+00   # ZU(4,4)
  4  5     0.00000000E+00   # ZU(4,5)
  4  6    -4.27858606E-01   # ZU(4,6)
  5  1     1.00000000E+00   # ZU(5,1)
  5  2     0.00000000E+00   # ZU(5,2)
  5  3     0.00000000E+00   # ZU(5,3)
  5  4    -2.02956466E-05   # ZU(5,4)
  5  5     0.00000000E+00   # ZU(5,5)
  5  6     0.00000000E+00   # ZU(5,6)
  6  1     0.00000000E+00   # ZU(6,1)
  6  2    -9.99956942E-01   # ZU(6,2)
  6  3     0.00000000E+00   # ZU(6,3)
  6  4     0.00000000E+00   # ZU(6,4)
  6  5     9.27975914E-03   # ZU(6,5)
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
     1     8.64566212E+02   # SUSYScale
     2     9.11876000E+01   # LowScale
Block ALPHA
          -1.06834219E-01   # ArcSin(Pole(ZH(2,2)))
Block HMIX Q= 8.64566212E+02
     1     6.23808476E+02   # Mu
     2     9.67481745E+00   # vu/vd
     3     2.44130978E+02   # Sqrt(Sqr(vd) + Sqr(vu))
     4     5.26379373E+05   # Sqr(MAh(2))
   101     5.38320475E+04   # BMu
   102     2.50999305E+01   # vd
   103     2.42837245E+02   # vu
Block Au Q= 8.64566212E+02
  1  1    -1.13135345E+03   # TYu(1,1)/Yu(1,1)
  2  2    -1.13134844E+03   # TYu(2,2)/Yu(2,2)
  3  3    -8.73285723E+02   # TYu(3,3)/Yu(3,3)
Block Ad Q= 8.64566212E+02
  1  1    -1.38518584E+03   # TYd(1,1)/Yd(1,1)
  2  2    -1.38518119E+03   # TYd(2,2)/Yd(2,2)
  3  3    -1.29472423E+03   # TYd(3,3)/Yd(3,3)
Block Ae Q= 8.64566212E+02
  1  1    -2.99353463E+02   # TYe(1,1)/Ye(1,1)
  2  2    -2.99347667E+02   # TYe(2,2)/Ye(2,2)
  3  3    -2.97714149E+02   # TYe(3,3)/Ye(3,3)
Block MSOFT Q= 8.64566212E+02
     1     2.08850011E+02   # MassB
     2     3.87896014E+02   # MassWB
     3     1.11451862E+03   # MassG
    21     1.09331762E+05   # mHd2
    22    -3.77162157E+05   # mHu2
    31     3.53550629E+02   # SignedAbsSqrt(ml2(1,1))
    32     3.53547066E+02   # SignedAbsSqrt(ml2(2,2))
    33     3.52542851E+02   # SignedAbsSqrt(ml2(3,3))
    34     2.22260279E+02   # SignedAbsSqrt(me2(1,1))
    35     2.22248714E+02   # SignedAbsSqrt(me2(2,2))
    36     2.18968838E+02   # SignedAbsSqrt(me2(3,3))
    41     1.00699572E+03   # SignedAbsSqrt(mq2(1,1))
    42     1.00699317E+03   # SignedAbsSqrt(mq2(2,2))
    43     9.28348154E+02   # SignedAbsSqrt(mq2(3,3))
    44     9.68157834E+02   # SignedAbsSqrt(mu2(1,1))
    45     9.68155109E+02   # SignedAbsSqrt(mu2(2,2))
    46     7.97126835E+02   # SignedAbsSqrt(mu2(3,3))
    47     9.63384971E+02   # SignedAbsSqrt(md2(1,1))
    48     9.63382289E+02   # SignedAbsSqrt(md2(2,2))
    49     9.58601864E+02   # SignedAbsSqrt(md2(3,3))
)";

   std::stringstream istr(slha_input);

   MSSM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   MSSM_input_parameters input;
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

   MSSM_slha m(input);
   slha_io.fill(m);
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   MSSM_decays decays_with_HO(m, qedqcd, physical_input, flexibledecay_settings);

   // scalar Higgs

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0024454872837716721, 7e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              0.00011321683135547809, 3e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.00026234514953420564, 3e-13);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
   //                            0.0001976368796373175, 2e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
                              0.00023570018177709328, 1e-3);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
   //                            1.4826613977728886e-05, 3e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
                              2.4998529690114355e-05, 1e-3);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VGVG(&m, 0), 0.00029107803462142189, 2e-10);
   // h -> gamma gamma
   // without 2L QCD for squark
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 6.3284545616000571e-06, 4e-11);
   // with 2L QCD for squark
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 7.2710487609165856e-06, 4e-11);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVZ(&m, 0), 2.2198087946884208e-06, 2e-10);

   // pseudoscalar Higgs
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00043378061307989104, 4e-13);

   // ------------ tree-level decays ------------

   // Ah -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFdFd(&m, 1, 2, 2),
                              0.88658115467704091, 4e-13);
   // Ah -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFeFe(&m, 1, 2, 2),
                              0.14370135383550148, 4e-13);
   // Ah -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFuFu(&m, 1, 1, 1),
                              6.1500292100561543e-06, 3e-13);
   // Ah -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_conjVWmVWm(&m, 1),
                              3.926795985114403e-05, 2e-12);
   // Ah -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VZVZ(&m, 1),
                              1.0275363930925979e-05, 2e-12);

   // ------------ loop-induces decays ------------

   // Ah -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00043378061307989104, 4e-13);
   // Ah -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VPVP(&m, 1), 4.1843914309546067e-06, 4e-12);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   MSSM_decays decays_without_HO(m, qedqcd, physical_input, flexibledecay_settings);

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0013683069345404055, 3e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              7.6128491585600595e-05, 3e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.00025956002881119691, 3e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVP(&m, 0), 7.1507061728988442e-06, 5e-11);

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VGVG(&m, 0), 9.3812245440699979e-05, 2e-10);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVZ(&m, 0), 2.2198087946884208e-06, 2e-10);

   // Ah -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00043378061307989104, 4e-13);
   // Ah -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Ah_to_VPVP(&m, 1), 4.0630569784682124e-06, 4e-12);
   // Ah -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Ah_to_VPVZ(&m, 1), 9.3179926421364644e-06, 2e-11);

   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVZ(&m, 0), 2.2198087946884208e-06, 2e-10);

   // Sd5 -> Chi d
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Sd_to_ChiFd(&m, 5, 0, 0),
                              0.16874114148504477, 9e-14);
   // Su5 -> Cha d
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Su_to_barChaFd(&m, 4, 0, 0),
                              6.1539820064913817, 8e-14);

   // Sv -> Chi Fv
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Sv_to_FvChi(&m, 2, 0, 0),
                              0.21565996414399113, 7e-15);

   // Se -> Chi Fv
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Se_to_ChiFe(&m, 2, 0, 0),
                              0.055308969643458182, 2e-13);

   // hh(2) -> hh(1) hh(1)
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_hhhh(&m, 1, 0, 0),
                              0.0055271603394791615, 5e-13);
}
