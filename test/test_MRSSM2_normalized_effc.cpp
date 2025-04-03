
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MRSSM2_normalized_effc

#include <boost/test/unit_test.hpp>

#include "test_MRSSM2.hpp"
#include "MRSSM2_two_scale_model.hpp"
#include "MRSSM2_two_scale_spectrum_generator.hpp"
#include "decays/MRSSM2_decays.hpp"
#include "MRSSM2_slha_io.hpp"
#include "decays/experimental_constraints.hpp"

using namespace flexiblesusy;

char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.6.1
     5   MRSSM2
     9   4.14.5
Block MODSEL                 # Select model
    6   1                    # flavour violation
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
   15   1                    # calculate all observables
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
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.385               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block VCKMIN                 # CKM matrix input (Wolfenstein parameters)
    1   2.272000000e-01      # lambda(MZ) SM DR-bar
    2   8.180000000e-01      # A(MZ) SM DR-bar
    3   2.210000000e-01      # rhobar(MZ) SM DR-bar
    4   3.400000000e-01      # etabar(MZ) SM DR-bar
# Block UPMNSIN                # PMNS matrix input
#     1   5.837630000e-01      # theta_12
#     2   7.695840000e-01      # theta_23
#     3   1.549480000e-01      # theta_13
#     4   0.000000000e+00      # delta
#     5   0.000000000e+00      # alpha_1
#     6   0.000000000e+00      # alpha_2
Block FlexibleDecay
   0   1                    # calculate decays (0 = no, 1 = yes)
   1   1e-5                 # minimum BR to print
   2   4                    # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
   3   1                    # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
   4   2                    # off-shell decays into VV pair
Block MINPAR
     3     5.00000000E+00   # TanBeta
Block EXTPAR
     0     1.00000000E+03   # Ms
Block HMIXIN
   101     1.00000000E+03   # BMuInput
   301    -1.00000000E-01   # LamSDInput
   302    -1.00000000E-01   # LamSUInput
   303    -1.00000000E-01   # LamTDInput
   304    -1.00000000E-01   # LamTUInput
   201     1.00000000E+03   # MuDInput
   202     1.00000000E+03   # MuUInput
Block MSD2IN
  1  1     1.00000000E+06   # md2Input(1,1)
  1  2     0.00000000E+00   # md2Input(1,2)
  1  3     0.00000000E+00   # md2Input(1,3)
  2  1     0.00000000E+00   # md2Input(2,1)
  2  2     1.00000000E+06   # md2Input(2,2)
  2  3     0.00000000E+00   # md2Input(2,3)
  3  1     0.00000000E+00   # md2Input(3,1)
  3  2     0.00000000E+00   # md2Input(3,2)
  3  3     1.00000000E+06   # md2Input(3,3)
Block MSOFTIN
   300     1.00000000E+03   # MDBSInput
   302     1.00000000E+03   # MDGocInput
   301     1.00000000E+03   # MDWBTInput
   111     1.00000000E+06   # moc2Input
    50     1.00000000E+06   # mRd2Input
    51     1.00000000E+06   # mRu2Input
   110     1.00000000E+06   # mT2Input
Block MSE2IN
  1  1     1.00000000E+06   # me2Input(1,1)
  1  2     0.00000000E+00   # me2Input(1,2)
  1  3     0.00000000E+00   # me2Input(1,3)
  2  1     0.00000000E+00   # me2Input(2,1)
  2  2     1.00000000E+06   # me2Input(2,2)
  2  3     0.00000000E+00   # me2Input(2,3)
  3  1     0.00000000E+00   # me2Input(3,1)
  3  2     0.00000000E+00   # me2Input(3,2)
  3  3     1.00000000E+06   # me2Input(3,3)
Block MSL2IN
  1  1     1.00000000E+06   # ml2Input(1,1)
  1  2     1.00000000E+03   # ml2Input(1,2)
  1  3     0.00000000E+00   # ml2Input(1,3)
  2  1     1.00000000E+03   # ml2Input(2,1)
  2  2     1.00000000E+06   # ml2Input(2,2)
  2  3     0.00000000E+00   # ml2Input(2,3)
  3  1     0.00000000E+00   # ml2Input(3,1)
  3  2     0.00000000E+00   # ml2Input(3,2)
  3  3     1.00000000E+06   # ml2Input(3,3)
Block MSQ2IN
  1  1     1.00000000E+06   # mq2Input(1,1)
  1  2     0.00000000E+00   # mq2Input(1,2)
  1  3     0.00000000E+00   # mq2Input(1,3)
  2  1     0.00000000E+00   # mq2Input(2,1)
  2  2     1.00000000E+06   # mq2Input(2,2)
  2  3     0.00000000E+00   # mq2Input(2,3)
  3  1     0.00000000E+00   # mq2Input(3,1)
  3  2     0.00000000E+00   # mq2Input(3,2)
  3  3     1.00000000E+06   # mq2Input(3,3)
Block NMSSMRUNIN
    10     1.00000000E+06   # mS2Input
Block MSU2IN
  1  1     1.00000000E+06   # mu2Input(1,1)
  1  2     0.00000000E+00   # mu2Input(1,2)
  1  3     0.00000000E+00   # mu2Input(1,3)
  2  1     0.00000000E+00   # mu2Input(2,1)
  2  2     1.00000000E+06   # mu2Input(2,2)
  2  3     0.00000000E+00   # mu2Input(2,3)
  3  1     0.00000000E+00   # mu2Input(3,1)
  3  2     0.00000000E+00   # mu2Input(3,2)
  3  3     1.00000000E+06   # mu2Input(3,3)
Block gauge Q= 1.00000000E+03
     1     3.62109856E-01   # g1 * 0.7745966692414834
     2     6.41073466E-01   # g2
     3     1.06728947E+00   # g3
Block Yu Q= 1.00000000E+03
  1  1     7.37742725E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.35845591E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.99634427E-01   # Yu(3,3)
Block Yd Q= 1.00000000E+03
  1  1     7.06528699E-05   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     1.54695813E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     7.12921031E-02   # Yd(3,3)
Block Ye Q= 1.00000000E+03
  1  1     1.43788747E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     2.97309565E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     5.00057951E-02   # Ye(3,3)
Block MSQ2 Q= 1.00000000E+03
  1  1     1.00000003E+06   # mq2(1,1)
  1  2    -1.15271744E-05   # mq2(1,2)
  1  3     2.54981370E-04   # mq2(1,3)
  2  1    -1.15271744E-05   # mq2(2,1)
  2  2     1.00000003E+06   # mq2(2,2)
  2  3    -1.41896179E-03   # mq2(2,3)
  3  1     2.54981370E-04   # mq2(3,1)
  3  2    -1.41896179E-03   # mq2(3,2)
  3  3     1.00000007E+06   # mq2(3,3)
Block MSE2 Q= 1.00000000E+03
  1  1     1.00000000E+06   # me2(1,1)
  1  2     1.72071753E-11   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     1.72071753E-11   # me2(2,1)
  2  2     1.00000000E+06   # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     1.00000000E+06   # me2(3,3)
Block MSL2 Q= 1.00000000E+03
  1  1     1.00000000E+06   # ml2(1,1)
  1  2     1.00000000E+03   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     1.00000000E+03   # ml2(2,1)
  2  2     1.00000000E+06   # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     1.00000000E+06   # ml2(3,3)
Block MSU2 Q= 1.00000000E+03
  1  1     1.00000003E+06   # mu2(1,1)
  1  2     9.84986832E-15   # mu2(1,2)
  1  3     9.51220997E-13   # mu2(1,3)
  2  1     9.84988187E-15   # mu2(2,1)
  2  2     1.00000003E+06   # mu2(2,2)
  2  3     7.89637698E-09   # mu2(2,3)
  3  1     9.51220553E-13   # mu2(3,1)
  3  2     7.89637601E-09   # mu2(3,2)
  3  3     1.00000010E+06   # mu2(3,3)
Block MSD2 Q= 1.00000000E+03
  1  1     1.00000003E+06   # md2(1,1)
  1  2     1.00014441E-12   # md2(1,2)
  1  3    -1.48252011E-09   # md2(1,3)
  2  1     1.00008890E-12   # md2(2,1)
  2  2     1.00000003E+06   # md2(2,2)
  2  3     1.80649921E-07   # md2(2,3)
  3  1    -1.48252016E-09   # md2(3,1)
  3  2     1.80649932E-07   # md2(3,2)
  3  3     1.00000004E+06   # md2(3,3)
Block MSOFT Q= 1.00000000E+03
    21    -1.00598832E+06   # mHd2
    22    -9.86792837E+05   # mHu2
   110     1.00000000E+06   # mT2
   111     1.00000008E+06   # moc2
   300     1.00000001E+03   # MDBS
   301     1.00000000E+03   # MDWBT
   302     9.99999889E+02   # MDGoc
    50     1.00000000E+06   # mRd2
    51     1.00000000E+06   # mRu2
Block NMSSMRUN Q= 1.00000000E+03
    10     1.00000000E+06   # mS2
     5    -1.76605894E-01   # vS
Block MASS
        24     8.04262848E+01   # VWm
       404     1.42821630E+03   # SRdp
       403     1.42667975E+03   # SRum
   1000021     1.12251521E+03   # Glu
   3000022     1.13838851E+03   # sigmaO
   3000021     2.22702985E+03   # phiO
   2000024     9.67010105E+02   # Cha2(1)
   2000037     1.06053307E+03   # Cha2(2)
   1000024     1.00458574E+03   # Cha1(1)
   1000037     1.02950396E+03   # Cha1(2)
       401     1.42706132E+03   # Rh(1)
       402     1.42833286E+03   # Rh(2)
   1000012     1.01017224E+03   # Sv(1)
   1000014     1.01066934E+03   # Sv(2)
   1000016     1.01116257E+03   # Sv(3)
   1000022     9.73658992E+02   # Chi(1)
   1000023     9.98787128E+02   # Chi(2)
   1000025     1.00905248E+03   # Chi(3)
   1000035     1.05027881E+03   # Chi(4)
        37     8.98947832E+01   # Hpm(2)
        47     1.03476501E+03   # Hpm(3)
        57     2.22384542E+03   # Hpm(4)
        25     6.57268566E+01   # hh(1)
        35     9.70046477E+01   # hh(2)
        45     2.22348245E+03   # hh(3)
        55     2.22441207E+03   # hh(4)
        36     7.20555296E+01   # Ah(2)
        46     1.00029610E+03   # Ah(3)
        56     1.03459487E+03   # Ah(4)
   1000001     1.05594527E+03   # Sd(1)
   1000003     1.05594527E+03   # Sd(2)
   1000005     1.05595379E+03   # Sd(3)
   2000001     1.06579881E+03   # Sd(4)
   2000003     1.06579883E+03   # Sd(5)
   2000005     1.06674918E+03   # Sd(6)
   1000011     1.00558845E+03   # Se(1)
   1000013     1.00559638E+03   # Se(2)
   1000015     1.00560468E+03   # Se(3)
   2000011     1.01145577E+03   # Se(4)
   2000013     1.01194197E+03   # Se(5)
   2000015     1.01242562E+03   # Se(6)
   1000002     1.05664271E+03   # Su(1)
   1000004     1.05664292E+03   # Su(2)
   1000006     1.06505407E+03   # Su(3)
   2000002     1.06505429E+03   # Su(4)
   2000004     1.06948947E+03   # Su(5)
   2000006     1.07808171E+03   # Su(6)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.29077915E-04   # Fe(1)
        13     1.07290815E-01   # Fe(2)
        15     1.78570576E+00   # Fe(3)
         1     4.44464601E-03   # Fd(1)
         3     8.75553618E-02   # Fd(2)
         5     3.46472467E+00   # Fd(3)
         2     2.30337494E-03   # Fu(1)
         4     8.42093264E-01   # Fu(2)
         6     1.77826061E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.11732679E+01   # VZ
Block U1MIX
  1  1    -4.15955693E-01   # Re(UM1(1,1))
  1  2     9.09384881E-01   # Re(UM1(1,2))
  2  1     9.09384881E-01   # Re(UM1(2,1))
  2  2     4.15955693E-01   # Re(UM1(2,2))
Block U2MIX
  1  1    -6.06588090E-01   # Re(UM2(1,1))
  1  2    -7.95016282E-01   # Re(UM2(1,2))
  2  1     7.95016282E-01   # Re(UM2(2,1))
  2  2    -6.06588090E-01   # Re(UM2(2,2))
Block V1MIX
  1  1    -4.03914007E-01   # Re(UP1(1,1))
  1  2     9.14796958E-01   # Re(UP1(1,2))
  2  1     9.14796958E-01   # Re(UP1(2,1))
  2  2     4.03914007E-01   # Re(UP1(2,2))
Block V2MIX
  1  1    -6.55340052E-01   # Re(UP2(1,1))
  1  2     7.55333976E-01   # Re(UP2(1,2))
  2  1     7.55333976E-01   # Re(UP2(2,1))
  2  2     6.55340052E-01   # Re(UP2(2,2))
Block PSEUDOSCALARMIX
  1  1     2.38556140E-01   # ZA(1,1)
  1  2    -9.71128708E-01   # ZA(1,2)
  1  3     1.70186260E-07   # ZA(1,3)
  1  4    -1.35534914E-06   # ZA(1,4)
  2  1     9.71128708E-01   # ZA(2,1)
  2  2     2.38556140E-01   # ZA(2,2)
  2  3    -8.12517432E-08   # ZA(2,3)
  2  4     1.06754027E-07   # ZA(2,4)
  3  1     3.93425582E-08   # ZA(3,1)
  3  2     1.78325395E-07   # ZA(3,2)
  3  3     9.99988876E-01   # ZA(3,3)
  3  4     4.71676592E-03   # ZA(3,4)
  4  1     2.19471830E-07   # ZA(4,1)
  4  2    -1.34254134E-06   # ZA(4,2)
  4  3    -4.71676592E-03   # ZA(4,3)
  4  4     9.99988876E-01   # ZA(4,4)
Block DSQMIX
  1  1     0.00000000E+00   # ZD(1,1)
  1  2     0.00000000E+00   # ZD(1,2)
  1  3     0.00000000E+00   # ZD(1,3)
  1  4     1.00000000E+00   # ZD(1,4)
  1  5     4.47765580E-06   # ZD(1,5)
  1  6    -2.01293688E-06   # ZD(1,6)
  2  1     0.00000000E+00   # ZD(2,1)
  2  2     0.00000000E+00   # ZD(2,2)
  2  3     0.00000000E+00   # ZD(2,3)
  2  4     4.47716156E-06   # ZD(2,4)
  2  5    -9.99999970E-01   # ZD(2,5)
  2  6    -2.45461952E-04   # ZD(2,6)
  3  1    -0.00000000E+00   # ZD(3,1)
  3  2    -0.00000000E+00   # ZD(3,2)
  3  3    -0.00000000E+00   # ZD(3,3)
  3  4    -2.01403591E-06   # ZD(3,4)
  3  5     2.45461942E-04   # ZD(3,5)
  3  6    -9.99999970E-01   # ZD(3,6)
  4  1     9.84473512E-01   # ZD(4,1)
  4  2     1.75533191E-01   # ZD(4,2)
  4  3    -5.71083316E-05   # ZD(4,3)
  4  4    -0.00000000E+00   # ZD(4,4)
  4  5    -0.00000000E+00   # ZD(4,5)
  4  6    -0.00000000E+00   # ZD(4,6)
  5  1     1.75376050E-01   # ZD(5,1)
  5  2    -9.83605846E-01   # ZD(5,2)
  5  3    -4.19854904E-02   # ZD(5,3)
  5  4     0.00000000E+00   # ZD(5,4)
  5  5     0.00000000E+00   # ZD(5,5)
  5  6     0.00000000E+00   # ZD(5,6)
  6  1     7.42601920E-03   # ZD(6,1)
  6  2    -4.13235877E-02   # ZD(6,2)
  6  3     9.99118219E-01   # ZD(6,3)
  6  4     0.00000000E+00   # ZD(6,4)
  6  5     0.00000000E+00   # ZD(6,5)
  6  6     0.00000000E+00   # ZD(6,6)
Block SELMIX
  1  1     0.00000000E+00   # ZE(1,1)
  1  2     0.00000000E+00   # ZE(1,2)
  1  3     0.00000000E+00   # ZE(1,3)
  1  4    -1.00000000E+00   # ZE(1,4)
  1  5     5.46975111E-06   # ZE(1,5)
  1  6     0.00000000E+00   # ZE(1,6)
  2  1    -0.00000000E+00   # ZE(2,1)
  2  2    -0.00000000E+00   # ZE(2,2)
  2  3    -0.00000000E+00   # ZE(2,3)
  2  4    -5.46975111E-06   # ZE(2,4)
  2  5    -1.00000000E+00   # ZE(2,5)
  2  6    -0.00000000E+00   # ZE(2,6)
  3  1     0.00000000E+00   # ZE(3,1)
  3  2     0.00000000E+00   # ZE(3,2)
  3  3     0.00000000E+00   # ZE(3,3)
  3  4     0.00000000E+00   # ZE(3,4)
  3  5     0.00000000E+00   # ZE(3,5)
  3  6     1.00000000E+00   # ZE(3,6)
  4  1    -7.07111125E-01   # ZE(4,1)
  4  2     7.07102437E-01   # ZE(4,2)
  4  3     0.00000000E+00   # ZE(4,3)
  4  4     0.00000000E+00   # ZE(4,4)
  4  5     0.00000000E+00   # ZE(4,5)
  4  6     0.00000000E+00   # ZE(4,6)
  5  1     0.00000000E+00   # ZE(5,1)
  5  2     0.00000000E+00   # ZE(5,2)
  5  3     1.00000000E+00   # ZE(5,3)
  5  4     0.00000000E+00   # ZE(5,4)
  5  5     0.00000000E+00   # ZE(5,5)
  5  6     0.00000000E+00   # ZE(5,6)
  6  1    -7.07102437E-01   # ZE(6,1)
  6  2    -7.07111125E-01   # ZE(6,2)
  6  3    -0.00000000E+00   # ZE(6,3)
  6  4    -0.00000000E+00   # ZE(6,4)
  6  5    -0.00000000E+00   # ZE(6,5)
  6  6    -0.00000000E+00   # ZE(6,6)
Block SCALARMIX
  1  1     9.30263639E-01   # ZH(1,1)
  1  2     3.66864112E-01   # ZH(1,2)
  1  3     7.96058458E-04   # ZH(1,3)
  1  4     4.43294946E-03   # ZH(1,4)
  2  1     3.66818780E-01   # ZH(2,1)
  2  2    -9.29862128E-01   # ZH(2,2)
  2  3     1.16085940E-02   # ZH(2,3)
  2  4    -2.58001274E-02   # ZH(2,4)
  3  1     5.00665039E-03   # ZH(3,1)
  3  2    -2.48905349E-02   # ZH(3,2)
  3  3     6.51034366E-02   # ZH(3,3)
  3  4     9.97555481E-01   # ZH(3,4)
  4  1    -5.33642783E-03   # ZH(4,1)
  4  2     1.21494049E-02   # ZH(4,2)
  4  3     9.97810678E-01   # ZH(4,3)
  4  4    -6.47901622E-02   # ZH(4,4)
Block RHMIX
  1  1     3.97946315E-03   # ZHR(1,1)
  1  2     9.99992082E-01   # ZHR(1,2)
  2  1     9.99992082E-01   # ZHR(2,1)
  2  2    -3.97946315E-03   # ZHR(2,2)
Block N1MIX
  1  1     4.61503745E-01   # Re(ZN1(1,1))
  1  2    -4.61411757E-01   # Re(ZN1(1,2))
  1  3    -1.75189425E-01   # Re(ZN1(1,3))
  1  4    -7.37171723E-01   # Re(ZN1(1,4))
  2  1     8.11929040E-01   # Re(ZN1(2,1))
  2  2     3.81709442E-01   # Re(ZN1(2,2))
  2  3    -2.85269999E-01   # Re(ZN1(2,3))
  2  4     3.37179721E-01   # Re(ZN1(2,4))
  3  1    -3.06754856E-01   # Re(ZN1(3,1))
  3  2    -1.54090932E-01   # Re(ZN1(3,2))
  3  3    -9.30795097E-01   # Re(ZN1(3,3))
  3  4     1.25610232E-01   # Re(ZN1(3,4))
  4  1    -1.83540149E-01   # Re(ZN1(4,1))
  4  2     7.85909077E-01   # Re(ZN1(4,2))
  4  3    -1.46799796E-01   # Re(ZN1(4,3))
  4  4    -5.71935098E-01   # Re(ZN1(4,4))
Block N2MIX
  1  1     4.83784620E-01   # Re(ZN2(1,1))
  1  2    -4.94862959E-01   # Re(ZN2(1,2))
  1  3     1.69624354E-01   # Re(ZN2(1,3))
  1  4    -7.01634286E-01   # Re(ZN2(1,4))
  2  1     8.02344380E-01   # Re(ZN2(2,1))
  2  2     3.94369709E-01   # Re(ZN2(2,2))
  2  3     2.86586803E-01   # Re(ZN2(2,3))
  2  4     3.44360325E-01   # Re(ZN2(2,4))
  3  1    -3.07978309E-01   # Re(ZN2(3,1))
  3  2    -1.57060096E-01   # Re(ZN2(3,2))
  3  3     9.30203210E-01   # Re(ZN2(3,3))
  3  4     1.23302374E-01   # Re(ZN2(3,4))
  4  1    -1.65364141E-01   # Re(ZN2(4,1))
  4  2     7.58231700E-01   # Re(ZN2(4,2))
  4  3     1.54329423E-01   # Re(ZN2(4,3))
  4  4    -6.11491471E-01   # Re(ZN2(4,4))
Block CHARGEMIX
  1  1     1.55583881E-01   # ZP(1,1)
  1  2    -9.87520172E-01   # ZP(1,2)
  1  3     1.72844949E-02   # ZP(1,3)
  1  4     1.72861782E-02   # ZP(1,4)
  2  1    -9.87808234E-01   # ZP(2,1)
  2  2    -1.55400109E-01   # ZP(2,2)
  2  3     6.54590778E-03   # ZP(2,3)
  2  4     6.54592230E-03   # ZP(2,4)
  3  1     7.21943130E-06   # ZP(3,1)
  3  2     3.65888162E-05   # ZP(3,2)
  3  3    -7.06127856E-01   # ZP(3,3)
  3  4     7.08084352E-01   # ZP(3,4)
  4  1     5.34300704E-03   # ZP(4,1)
  4  2     2.55873917E-02   # ZP(4,2)
  4  3     7.07843095E-01   # ZP(4,3)
  4  4     7.05885890E-01   # ZP(4,4)
Block USQMIX
  1  1     0.00000000E+00   # ZU(1,1)
  1  2     0.00000000E+00   # ZU(1,2)
  1  3     0.00000000E+00   # ZU(1,3)
  1  4     1.00000000E+00   # ZU(1,4)
  1  5    -5.32754896E-11   # ZU(1,5)
  1  6    -2.22465617E-12   # ZU(1,6)
  2  1     0.00000000E+00   # ZU(2,1)
  2  2     0.00000000E+00   # ZU(2,2)
  2  3     0.00000000E+00   # ZU(2,3)
  2  4    -5.32754896E-11   # ZU(2,4)
  2  5    -1.00000000E+00   # ZU(2,5)
  2  6     1.78687179E-08   # ZU(2,6)
  3  1     9.75032192E-01   # ZU(3,1)
  3  2     2.22054998E-01   # ZU(3,2)
  3  3     1.95009600E-03   # ZU(3,3)
  3  4    -0.00000000E+00   # ZU(3,4)
  3  5    -0.00000000E+00   # ZU(3,5)
  3  6    -0.00000000E+00   # ZU(3,6)
  4  1     2.21937700E-01   # ZU(4,1)
  4  2    -9.74146366E-01   # ZU(4,2)
  4  3    -4.22198313E-02   # ZU(4,3)
  4  4     0.00000000E+00   # ZU(4,4)
  4  5     0.00000000E+00   # ZU(4,5)
  4  6     0.00000000E+00   # ZU(4,6)
  5  1     0.00000000E+00   # ZU(5,1)
  5  2     0.00000000E+00   # ZU(5,2)
  5  3     0.00000000E+00   # ZU(5,3)
  5  4    -2.22465712E-12   # ZU(5,4)
  5  5    -1.78687181E-08   # ZU(5,5)
  5  6    -1.00000000E+00   # ZU(5,6)
  6  1     7.47544562E-03   # ZU(6,1)
  6  2    -4.15984945E-02   # ZU(6,2)
  6  3     9.99106442E-01   # ZU(6,3)
  6  4     0.00000000E+00   # ZU(6,4)
  6  5     0.00000000E+00   # ZU(6,5)
  6  6     0.00000000E+00   # ZU(6,6)
Block SNUMIX
  1  1    -7.07109115E-01   # ZV(1,1)
  1  2     7.07104448E-01   # ZV(1,2)
  1  3     0.00000000E+00   # ZV(1,3)
  2  1     0.00000000E+00   # ZV(2,1)
  2  2     0.00000000E+00   # ZV(2,2)
  2  3     1.00000000E+00   # ZV(2,3)
  3  1    -7.07104448E-01   # ZV(3,1)
  3  2    -7.07109115E-01   # ZV(3,2)
  3  3    -0.00000000E+00   # ZV(3,3)
Block UELMIX
  1  1     1.00000000E+00   # Re(ZEL(1,1))
  1  2     5.57315388E-07   # Re(ZEL(1,2))
  1  3     0.00000000E+00   # Re(ZEL(1,3))
  2  1    -5.57315388E-07   # Re(ZEL(2,1))
  2  2     1.00000000E+00   # Re(ZEL(2,2))
  2  3     0.00000000E+00   # Re(ZEL(2,3))
  3  1     0.00000000E+00   # Re(ZEL(3,1))
  3  2     0.00000000E+00   # Re(ZEL(3,2))
  3  3     1.00000000E+00   # Re(ZEL(3,3))
Block UERMIX
  1  1     1.00000000E+00   # Re(ZER(1,1))
  1  2     5.44762710E-09   # Re(ZER(1,2))
  1  3     0.00000000E+00   # Re(ZER(1,3))
  2  1    -5.44762710E-09   # Re(ZER(2,1))
  2  2     1.00000000E+00   # Re(ZER(2,2))
  2  3     0.00000000E+00   # Re(ZER(2,3))
  3  1     0.00000000E+00   # Re(ZER(3,1))
  3  2     0.00000000E+00   # Re(ZER(3,2))
  3  3     1.00000000E+00   # Re(ZER(3,3))
Block UDLMIX
  1  1     9.99999998E-01   # Re(ZDL(1,1))
  1  2     2.96158613E-06   # Re(ZDL(1,2))
  1  3    -6.65065731E-05   # Re(ZDL(1,3))
  2  1    -2.93695099E-06   # Re(ZDL(2,1))
  2  2     9.99999931E-01   # Re(ZDL(2,2))
  2  3     3.70413735E-04   # Re(ZDL(2,3))
  3  1     6.65076655E-05   # Re(ZDL(3,1))
  3  2    -3.70413539E-04   # Re(ZDL(3,2))
  3  3     9.99999929E-01   # Re(ZDL(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(ZDR(1,1))
  1  2     2.63518091E-07   # Re(ZDR(1,2))
  1  3    -1.25014992E-07   # Re(ZDR(1,3))
  2  1    -2.63516129E-07   # Re(ZDR(2,1))
  2  2     1.00000000E+00   # Re(ZDR(2,2))
  2  3     1.56985285E-05   # Re(ZDR(2,3))
  3  1     1.25019129E-07   # Re(ZDR(3,1))
  3  2    -1.56985285E-05   # Re(ZDR(3,2))
  3  3     1.00000000E+00   # Re(ZDR(3,3))
Block UULMIX
  1  1     9.73845644E-01   # Re(ZUL(1,1))
  1  2     2.27200871E-01   # Re(ZUL(1,2))
  1  3     2.10368507E-03   # Re(ZUL(1,3))
  2  1    -2.27095337E-01   # Re(ZUL(2,1))
  2  2     9.73014965E-01   # Re(ZUL(2,2))
  2  3     4.08605764E-02   # Re(ZUL(2,3))
  3  1     7.23664147E-03   # Re(ZUL(3,1))
  3  2    -4.02696314E-02   # Re(ZUL(3,2))
  3  3     9.99162643E-01   # Re(ZUL(3,3))
Block UURMIX
  1  1     1.00000000E+00   # Re(ZUR(1,1))
  1  2    -7.49078919E-09   # Re(ZUR(1,2))
  1  3    -1.29544109E-09   # Re(ZUR(1,3))
  2  1     7.49077654E-09   # Re(ZUR(2,1))
  2  2     1.00000000E+00   # Re(ZUR(2,2))
  2  3    -9.76722652E-06   # Re(ZUR(2,3))
  3  1     1.29551426E-09   # Re(ZUR(3,1))
  3  2     9.76722652E-06   # Re(ZUR(3,2))
  3  3     1.00000000E+00   # Re(ZUR(3,3))
Block VCKM Q= 1.00000000E+03
  1  1     9.73846499E-01   # Re(CKM)(1,1)
  1  2     2.27196801E-01   # Re(CKM)(1,2)
  1  3     2.14700952E-03   # Re(CKM)(1,3)
  2  1    -2.27086872E-01   # Re(CKM)(2,1)
  2  2     9.72981216E-01   # Re(CKM)(2,2)
  2  3     4.17025829E-02   # Re(CKM)(2,3)
  3  1     7.38569349E-03   # Re(CKM)(3,1)
  3  2    -4.10994721E-02   # Re(CKM)(3,2)
  3  3     9.99127762E-01   # Re(CKM)(3,3)
Block IMVCKM Q= 1.00000000E+03
  1  1     0.00000000E+00   # Im(CKM)(1,1)
  1  2     0.00000000E+00   # Im(CKM)(1,2)
  1  3     0.00000000E+00   # Im(CKM)(1,3)
  2  1     0.00000000E+00   # Im(CKM)(2,1)
  2  2    -0.00000000E+00   # Im(CKM)(2,2)
  2  3     0.00000000E+00   # Im(CKM)(2,3)
  3  1     0.00000000E+00   # Im(CKM)(3,1)
  3  2     0.00000000E+00   # Im(CKM)(3,2)
  3  3     0.00000000E+00   # Im(CKM)(3,3)
)";

BOOST_AUTO_TEST_CASE( test_MRSSM2_normalized_effective_couplings )
{

   std::stringstream istr(slha_input);

   MRSSM2_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   MRSSM2_input_parameters input;
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

   MRSSM2_slha m(input);
   slha_io.fill(m);
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   flexibledecay_settings.set(FlexibleDecay_settings::calculate_normalized_effc, 1);

   MRSSM2_decays decays = MRSSM2_decays(m, qedqcd, physical_input, flexibledecay_settings);
   decays.calculate_decays();
   const auto effc = get_normalized_effective_couplings(decays.get_neutral_higgs_effc(), physical_input, qedqcd, settings, flexibledecay_settings);

   // tolerance in %
   BOOST_CHECK_CLOSE_FRACTION(effc[1].gg.second, 1.0263117658410017,      5e-12);
   BOOST_CHECK_CLOSE_FRACTION(effc[1].gamgam.second, 0.35684597371984528, 2e-11);
   BOOST_CHECK_CLOSE_FRACTION(effc[1].Zgam.second, 0.72007052033522212,   2e-10);

   BOOST_CHECK_CLOSE_FRACTION(effc[1].ZZ.second, 0.84090671645096526, 6e-5);
   BOOST_CHECK_CLOSE_FRACTION(effc[1].WW.second, 0.83500663372868467, 2e-6);

   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].ee.second),     1.8693282675754028, 3e-16);
   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].mumu.second),   1.8693282402317875, 2e-16);
   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].tautau.second), 1.8693212737802281, 2e-16);

   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].bb.second), 1.8564986659704481,  4e-16);
   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].cc.second), 0.94769505194041959, 2e-16);
   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].ss.second), 1.828375808475458,   2e-16);
   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].dd.second), 1.7919566274448113,  2e-16);
   BOOST_CHECK_CLOSE_FRACTION(std::real(effc[1].uu.second), 0.94768964478700091, 3e-16);

   // Ah
   // -----------------------------------------------------------------
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].bb.second), 4.9471107239345873,  0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].cc.second), 0.24445945353375145, 0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].ss.second), 4.8845714244777474,  0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].dd.second), 4.8233392849452121,  0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].uu.second), 0.24724669904197341, 0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].ee.second), 4.9521274449818486, 0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].mumu.second), 4.9521492793168838, 0.);
   BOOST_CHECK_CLOSE_FRACTION(std::imag(effc[2].tautau.second), 4.9581898246025222, 0.);

   // pure pseudoscalar has no renormalizable VV couplings and only gamma5 couplings to fermions
   BOOST_REQUIRE_EQUAL(effc[2].ZZ.second,                0.);
   BOOST_REQUIRE_EQUAL(effc[2].WW.second,                0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].bb.second),     0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].cc.second),     0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].ss.second),     0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].dd.second),     0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].uu.second),     0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].ee.second),     0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].mumu.second),   0.);
   BOOST_REQUIRE_EQUAL(std::real(effc[2].tautau.second), 0.);
}
