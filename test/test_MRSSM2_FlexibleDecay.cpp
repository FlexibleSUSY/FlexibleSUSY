
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MRSSM2_FlexibleDecay

#include <boost/test/unit_test.hpp>

#include "test_MRSSM2.hpp"
#include "MRSSM2_two_scale_model.hpp"
#include "MRSSM2_two_scale_spectrum_generator.hpp"
#include "decays/MRSSM2_decays.hpp"
#include "MRSSM2_slha_io.hpp"
#include "MRSSM2_mass_eigenstates_running.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MRSSM2_FlexibleDecay )
{

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
     3     4.00000000E+01   # TanBeta
Block EXTPAR
     0     1.00000000E+03   # Ms
Block HMIXIN
   101     4.00000000E+04   # BMuInput
   301     1.50000000E-01   # LamSDInput
   302    -1.50000000E-01   # LamSUInput
   303    -1.00000000E+00   # LamTDInput
   304    -1.15000000E+00   # LamTUInput
   201     4.00000000E+02   # MuDInput
   202     4.00000000E+02   # MuUInput
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
   300     2.50000000E+02   # MDBSInput
   302     1.00000000E+03   # MDGocInput
   301     5.00000000E+02   # MDWBTInput
   111     1.00000000E+06   # moc2Input
    50     4.90000000E+05   # mRd2Input
    51     1.00000000E+06   # mRu2Input
   110     9.00000000E+06   # mT2Input
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
Block NMSSMRUNIN
    10     4.00000000E+06   # mS2Input
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
     1     3.62350606E-01   # g1 * 0.7745966692414834
     2     6.42118963E-01   # g2
     3     1.06215840E+00   # g3
Block Yu Q= 1.00000000E+03
  1  1     7.32341299E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.33390202E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.83286581E-01   # Yu(3,3)
Block Yd Q= 1.00000000E+03
  1  1     5.69014312E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     1.24587493E-02   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     5.61642799E-01   # Yd(3,3)
Block Ye Q= 1.00000000E+03
  1  1     1.15519032E-04   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     2.38856226E-02   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     4.01531161E-01   # Ye(3,3)
Block HMIX Q= 1.00000000E+03
     1     0.00000000E+00   # Mu
   101     3.99999974E+04   # BMu
   102     6.17979538E+00   # vd
   103     2.41154149E+02   # vu
   310    -3.57900109E-01   # vT
   201     4.00000090E+02   # MuD
   202     3.99999878E+02   # MuU
   203     0.00000000E+00   # BMuD
   204     0.00000000E+00   # BMuU
   301     1.50000034E-01   # LamSD
   302    -1.49999954E-01   # LamSU
   303    -1.00000021E+00   # LamTD
   304    -1.14999963E+00   # LamTU
Block MSQ2 Q= 1.00000000E+03
  1  1     4.00000005E+06   # mq2(1,1)
  1  2     2.16076037E-04   # mq2(1,2)
  1  3    -5.26363540E-03   # mq2(1,3)
  2  1     2.16076037E-04   # mq2(2,1)
  2  2     4.00000005E+06   # mq2(2,2)
  2  3     2.92912673E-02   # mq2(2,3)
  3  1    -5.26363540E-03   # mq2(3,1)
  3  2     2.92912673E-02   # mq2(3,2)
  3  3     4.00000021E+06   # mq2(3,3)
Block MSE2 Q= 1.00000000E+03
  1  1     9.99999999E+05   # me2(1,1)
  1  2     0.00000000E+00   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     0.00000000E+00   # me2(2,1)
  2  2     9.99999998E+05   # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     9.99999888E+05   # me2(3,3)
Block MSL2 Q= 1.00000000E+03
  1  1     4.00000001E+06   # ml2(1,1)
  1  2     0.00000000E+00   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     0.00000000E+00   # ml2(2,1)
  2  2     4.00000001E+06   # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     3.99999995E+06   # ml2(3,3)
Block MSU2 Q= 1.00000000E+03
  1  1     1.00000004E+06   # mu2(1,1)
  1  2     8.13498574E-14   # mu2(1,2)
  1  3     6.90953136E-10   # mu2(1,3)
  2  1     8.13481227E-14   # mu2(2,1)
  2  2     1.00000004E+06   # mu2(2,2)
  2  3     6.11378928E-06   # mu2(2,3)
  3  1     6.90953133E-10   # mu2(3,1)
  3  2     6.11378928E-06   # mu2(3,2)
  3  3     9.99998614E+05   # mu2(3,3)
Block MSD2 Q= 1.00000000E+03
  1  1     1.00000004E+06   # md2(1,1)
  1  2     4.71937489E-11   # md2(1,2)
  1  3    -4.08043448E-07   # md2(1,3)
  2  1     4.71937489E-11   # md2(2,1)
  2  2     1.00000004E+06   # md2(2,2)
  2  3     4.97361335E-05   # md2(2,3)
  3  1    -4.08043448E-07   # md2(3,1)
  3  2     4.97361335E-05   # md2(3,2)
  3  3     1.00000177E+06   # md2(3,3)
Block MSOFT Q= 1.00000000E+03
    21     1.28465220E+06   # mHd2
    22    -3.33824153E+05   # mHu2
   110     8.99999999E+06   # mT2
   111     1.00000010E+06   # moc2
   300     2.50000000E+02   # MDBS
   301     5.00000001E+02   # MDWBT
   302     9.99999929E+02   # MDGoc
    50     4.90000092E+05   # mRd2
    51     9.99999885E+05   # mRu2
Block NMSSMRUN Q= 1.00000000E+03
    10     4.00000000E+06   # mS2
     5    -6.58611175E-02   # vS
Block MASS
        24     8.03946575E+01   # VWm
       404     8.91726546E+02   # SRdp
       403     1.17471104E+03   # SRum
   1000021     1.15530397E+03   # Glu
   3000022     1.15881785E+03   # sigmaO
   3000021     2.29249134E+03   # phiO
   2000024     4.25299413E+02   # Cha2(1)
   2000037     5.66773575E+02   # Cha2(2)
   1000024     4.09265442E+02   # Cha1(1)
   1000037     5.27489928E+02   # Cha1(2)
       401     8.95815263E+02   # Rh(1)
       402     1.16262647E+03   # Rh(2)
   1000012     2.00175391E+03   # Sv(1)
   1000014     2.00297474E+03   # Sv(2)
   1000016     2.00297908E+03   # Sv(3)
   1000022     2.50941894E+02   # Chi(1)
   1000023     4.09056207E+02   # Chi(2)
   1000025     4.22075658E+02   # Chi(3)
   1000035     5.47585018E+02   # Chi(4)
        37     1.23229115E+03   # Hpm(2)
        47     3.01676952E+03   # Hpm(3)
        57     3.17984669E+03   # Hpm(4)
        25     1.27060566E+02   # hh(1)
        35     1.22926626E+03   # hh(2)
        45     2.06032688E+03   # hh(3)
        55     3.17958470E+03   # hh(4)
        36     1.22925494E+03   # Ah(2)
        46     1.99954615E+03   # Ah(3)
        56     3.01657098E+03   # Ah(4)
   1000001     1.05730617E+03   # Sd(1)
   1000003     1.05919387E+03   # Sd(2)
   1000005     1.05919484E+03   # Sd(3)
   2000001     2.03399075E+03   # Sd(4)
   2000003     2.04190092E+03   # Sd(5)
   2000005     2.04190519E+03   # Sd(6)
   1000011     1.00136991E+03   # Se(1)
   1000013     1.00226372E+03   # Se(2)
   1000015     1.00226689E+03   # Se(3)
   2000011     2.00369715E+03   # Se(4)
   2000013     2.00492394E+03   # Se(5)
   2000015     2.00492830E+03   # Se(6)
   1000002     1.05839522E+03   # Su(1)
   1000004     1.05839538E+03   # Su(2)
   1000006     1.07014623E+03   # Su(3)
   2000002     2.03841573E+03   # Su(4)
   2000004     2.04048569E+03   # Su(5)
   2000006     2.04048811E+03   # Su(6)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.29111495E-04   # Fe(1)
        13     1.07287537E-01   # Fe(2)
        15     1.78530042E+00   # Fe(3)
         1     4.42593001E-03   # Fd(1)
         3     8.75538751E-02   # Fd(2)
         5     3.45866793E+00   # Fd(3)
         2     2.29131381E-03   # Fu(1)
         4     8.45162005E-01   # Fu(2)
         6     1.77307012E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.11536697E+01   # VZ
Block U1MIX
  1  1     2.74744790E-03   # Re(UM1(1,1))
  1  2     9.99996226E-01   # Re(UM1(1,2))
  2  1     9.99996226E-01   # Re(UM1(2,1))
  2  2    -2.74744790E-03   # Re(UM1(2,2))
Block U2MIX
  1  1     4.68016439E-01   # Re(UM2(1,1))
  1  2    -8.83719759E-01   # Re(UM2(1,2))
  2  1     8.83719759E-01   # Re(UM2(2,1))
  2  2     4.68016439E-01   # Re(UM2(2,2))
Block V1MIX
  1  1     1.09030525E-02   # Re(UP1(1,1))
  1  2     9.99940560E-01   # Re(UP1(1,2))
  2  1     9.99940560E-01   # Re(UP1(2,1))
  2  2    -1.09030525E-02   # Re(UP1(2,2))
Block V2MIX
  1  1     1.71904445E-01   # Re(UP2(1,1))
  1  2     9.85113629E-01   # Re(UP2(1,2))
  2  1     9.85113629E-01   # Re(UP2(2,1))
  2  2    -1.71904445E-01   # Re(UP2(2,2))
Block PSEUDOSCALARMIX
  1  1    -2.50514967E-02   # ZA(1,1)
  1  2     9.99686153E-01   # ZA(1,2)
  1  3    -6.92103957E-05   # ZA(1,3)
  1  4     1.11524894E-04   # ZA(1,4)
  2  1     9.99686162E-01   # ZA(2,1)
  2  2     2.50514967E-02   # ZA(2,2)
  2  3    -1.05059746E-06   # ZA(2,3)
  2  4     1.12046293E-06   # ZA(2,4)
  3  1     6.82326059E-07   # ZA(3,1)
  3  2    -6.91330221E-05   # ZA(3,2)
  3  3    -9.99999728E-01   # ZA(3,3)
  3  4    -7.34890714E-04   # ZA(3,4)
  4  1     1.67425612E-06   # ZA(4,1)
  4  2    -1.11568798E-04   # ZA(4,2)
  4  3    -7.34882997E-04   # ZA(4,3)
  4  4     9.99999724E-01   # ZA(4,4)
Block DSQMIX
  1  1    -0.00000000E+00   # ZD(1,1)
  1  2    -0.00000000E+00   # ZD(1,2)
  1  3    -0.00000000E+00   # ZD(1,3)
  1  4     1.48140447E-07   # ZD(1,4)
  1  5    -1.80584972E-05   # ZD(1,5)
  1  6    -1.00000000E+00   # ZD(1,6)
  2  1     0.00000000E+00   # ZD(2,1)
  2  2     0.00000000E+00   # ZD(2,2)
  2  3     0.00000000E+00   # ZD(2,3)
  2  4    -2.77307016E-07   # ZD(2,4)
  2  5    -1.00000000E+00   # ZD(2,5)
  2  6     1.80584971E-05   # ZD(2,6)
  3  1     0.00000000E+00   # ZD(3,1)
  3  2     0.00000000E+00   # ZD(3,2)
  3  3     0.00000000E+00   # ZD(3,3)
  3  4     1.00000000E+00   # ZD(3,4)
  3  5    -2.77304341E-07   # ZD(3,5)
  3  6     1.48145455E-07   # ZD(3,6)
  4  1     4.88721699E-03   # ZD(4,1)
  4  2    -2.72006501E-02   # ZD(4,2)
  4  3     9.99618047E-01   # ZD(4,3)
  4  4     0.00000000E+00   # ZD(4,4)
  4  5     0.00000000E+00   # ZD(4,5)
  4  6     0.00000000E+00   # ZD(4,6)
  5  1    -1.25596414E-01   # ZD(5,1)
  5  2     9.91697447E-01   # ZD(5,2)
  5  3     2.75991738E-02   # ZD(5,3)
  5  4     0.00000000E+00   # ZD(5,4)
  5  5     0.00000000E+00   # ZD(5,5)
  5  6     0.00000000E+00   # ZD(5,6)
  6  1    -9.92069381E-01   # ZD(6,1)
  6  2    -1.25683326E-01   # ZD(6,2)
  6  3     1.43033649E-03   # ZD(6,3)
  6  4    -0.00000000E+00   # ZD(6,4)
  6  5    -0.00000000E+00   # ZD(6,5)
  6  6    -0.00000000E+00   # ZD(6,6)
Block SELMIX
  1  1     0.00000000E+00   # ZE(1,1)
  1  2     0.00000000E+00   # ZE(1,2)
  1  3     0.00000000E+00   # ZE(1,3)
  1  4     0.00000000E+00   # ZE(1,4)
  1  5     0.00000000E+00   # ZE(1,5)
  1  6     1.00000000E+00   # ZE(1,6)
  2  1     0.00000000E+00   # ZE(2,1)
  2  2     0.00000000E+00   # ZE(2,2)
  2  3     0.00000000E+00   # ZE(2,3)
  2  4     0.00000000E+00   # ZE(2,4)
  2  5     1.00000000E+00   # ZE(2,5)
  2  6     0.00000000E+00   # ZE(2,6)
  3  1     0.00000000E+00   # ZE(3,1)
  3  2     0.00000000E+00   # ZE(3,2)
  3  3     0.00000000E+00   # ZE(3,3)
  3  4     1.00000000E+00   # ZE(3,4)
  3  5     0.00000000E+00   # ZE(3,5)
  3  6     0.00000000E+00   # ZE(3,6)
  4  1     0.00000000E+00   # ZE(4,1)
  4  2     0.00000000E+00   # ZE(4,2)
  4  3     1.00000000E+00   # ZE(4,3)
  4  4     0.00000000E+00   # ZE(4,4)
  4  5     0.00000000E+00   # ZE(4,5)
  4  6     0.00000000E+00   # ZE(4,6)
  5  1     0.00000000E+00   # ZE(5,1)
  5  2     1.00000000E+00   # ZE(5,2)
  5  3     0.00000000E+00   # ZE(5,3)
  5  4     0.00000000E+00   # ZE(5,4)
  5  5     0.00000000E+00   # ZE(5,5)
  5  6     0.00000000E+00   # ZE(5,6)
  6  1     1.00000000E+00   # ZE(6,1)
  6  2     0.00000000E+00   # ZE(6,2)
  6  3     0.00000000E+00   # ZE(6,3)
  6  4     0.00000000E+00   # ZE(6,4)
  6  5     0.00000000E+00   # ZE(6,5)
  6  6     0.00000000E+00   # ZE(6,6)
Block SCALARMIX
  1  1    -2.60123622E-02   # ZH(1,1)
  1  2    -9.99657149E-01   # ZH(1,2)
  1  3     1.78509883E-04   # ZH(1,3)
  1  4     2.98483135E-03   # ZH(1,4)
  2  1     9.99661619E-01   # ZH(2,1)
  2  2    -2.60120618E-02   # ZH(2,2)
  2  3     1.67182039E-05   # ZH(2,3)
  2  4     1.38556929E-04   # ZH(2,4)
  3  1     1.21066976E-05   # ZH(3,1)
  3  2    -1.80729770E-04   # ZH(3,2)
  3  3    -9.99999793E-01   # ZH(3,3)
  3  4    -6.17467555E-04   # ZH(3,4)
  4  1     6.08603375E-05   # ZH(4,1)
  4  2    -2.98731446E-03   # ZH(4,2)
  4  3     6.18005442E-04   # ZH(4,3)
  4  4    -9.99995345E-01   # ZH(4,4)
Block RHMIX
  1  1    -9.99999639E-01   # ZHR(1,1)
  1  2     8.49807450E-04   # ZHR(1,2)
  2  1    -8.49807450E-04   # ZHR(2,1)
  2  2    -9.99999639E-01   # ZHR(2,2)
Block N1MIX
  1  1     9.93812659E-01   # Re(ZN1(1,1))
  1  2     3.54137057E-02   # Re(ZN1(1,2))
  1  3    -2.66214317E-03   # Re(ZN1(1,3))
  1  4    -1.05238683E-01   # Re(ZN1(1,4))
  2  1    -9.37734293E-04   # Re(ZN1(2,1))
  2  2     1.22447426E-03   # Re(ZN1(2,2))
  2  3    -9.99856852E-01   # Re(ZN1(2,3))
  2  4     1.68492517E-02   # Re(ZN1(2,4))
  3  1     1.10956357E-01   # Re(ZN1(3,1))
  3  2    -3.58545649E-01   # Re(ZN1(3,2))
  3  3     1.50744962E-02   # Re(ZN1(3,3))
  3  4     9.26772067E-01   # Re(ZN1(3,4))
  4  1     4.91990158E-03   # Re(ZN1(4,1))
  4  2     9.32839422E-01   # Re(ZN1(4,2))
  4  3     7.20753241E-03   # Re(ZN1(4,3))
  4  4     3.60186699E-01   # Re(ZN1(4,4))
Block N2MIX
  1  1     9.99845443E-01   # Re(ZN2(1,1))
  1  2     1.58613853E-02   # Re(ZN2(1,2))
  1  3     1.82055639E-04   # Re(ZN2(1,3))
  1  4    -7.58118239E-03   # Re(ZN2(1,4))
  2  1    -4.87682769E-05   # Re(ZN2(2,1))
  2  2    -3.97310844E-04   # Re(ZN2(2,2))
  2  3     9.99859667E-01   # Re(ZN2(2,3))
  2  4     1.67477109E-02   # Re(ZN2(2,4))
  3  1     9.74898766E-03   # Re(ZN2(3,1))
  3  2    -1.41269095E-01   # Re(ZN2(3,2))
  3  3    -1.66345939E-02   # Re(ZN2(3,3))
  3  4     9.89783456E-01   # Re(ZN2(3,4))
  4  1    -1.46303068E-02   # Re(ZN2(4,1))
  4  2     9.89844079E-01   # Re(ZN2(4,2))
  4  3    -1.97565115E-03   # Re(ZN2(4,3))
  4  4     1.41388647E-01   # Re(ZN2(4,4))
Block CHARGEMIX
  1  1     2.51125867E-02   # ZP(1,1)
  1  2    -9.99680917E-01   # ZP(1,2)
  1  3    -2.00213201E-03   # ZP(1,3)
  1  4    -1.84765551E-03   # ZP(1,4)
  2  1    -9.99684628E-01   # ZP(2,1)
  2  2    -2.51123438E-02   # ZP(2,2)
  2  3    -8.80116470E-05   # ZP(2,3)
  2  4    -8.65276923E-05   # ZP(2,4)
  3  1     9.56888105E-07   # ZP(3,1)
  3  2    -1.45712048E-04   # ZP(3,2)
  3  3     7.16516053E-01   # ZP(3,3)
  3  4    -6.97570588E-01   # ZP(3,4)
  4  1     5.50352333E-05   # ZP(4,1)
  4  2     2.72274398E-03   # ZP(4,2)
  4  3    -6.97567724E-01   # ZP(4,3)
  4  4    -7.16513680E-01   # ZP(4,4)
Block USQMIX
  1  1     0.00000000E+00   # ZU(1,1)
  1  2     0.00000000E+00   # ZU(1,2)
  1  3     0.00000000E+00   # ZU(1,3)
  1  4     1.00000000E+00   # ZU(1,4)
  1  5    -4.41889299E-09   # ZU(1,5)
  1  6    -1.79858505E-10   # ZU(1,6)
  2  1     0.00000000E+00   # ZU(2,1)
  2  2     0.00000000E+00   # ZU(2,2)
  2  3     0.00000000E+00   # ZU(2,3)
  2  4    -4.41889270E-09   # ZU(2,4)
  2  5    -1.00000000E+00   # ZU(2,5)
  2  6     1.59127507E-06   # ZU(2,6)
  3  1     0.00000000E+00   # ZU(3,1)
  3  2     0.00000000E+00   # ZU(3,2)
  3  3     0.00000000E+00   # ZU(3,3)
  3  4    -1.79865536E-10   # ZU(3,4)
  3  5    -1.59127507E-06   # ZU(3,5)
  3  6    -1.00000000E+00   # ZU(3,6)
  4  1    -1.56974454E-01   # ZU(4,1)
  4  2     9.82832740E-01   # ZU(4,2)
  4  3     9.69475387E-02   # ZU(4,3)
  4  4     0.00000000E+00   # ZU(4,4)
  4  5     0.00000000E+00   # ZU(4,5)
  4  6     0.00000000E+00   # ZU(4,6)
  5  1    -9.87453546E-01   # ZU(5,1)
  5  2    -1.57898504E-01   # ZU(5,2)
  5  3     1.88595291E-03   # ZU(5,3)
  5  4     0.00000000E+00   # ZU(5,4)
  5  5     0.00000000E+00   # ZU(5,5)
  5  6     0.00000000E+00   # ZU(5,6)
  6  1     1.71614476E-02   # ZU(6,1)
  6  2    -9.54351444E-02   # ZU(6,2)
  6  3     9.95287706E-01   # ZU(6,3)
  6  4     0.00000000E+00   # ZU(6,4)
  6  5     0.00000000E+00   # ZU(6,5)
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
  1  1     9.99999988E-01   # Re(ZDL(1,1))
  1  2     6.22824967E-06   # Re(ZDL(1,2))
  1  3    -1.57234209E-04   # Re(ZDL(1,3))
  2  1    -6.09051785E-06   # Re(ZDL(2,1))
  2  2     9.99999616E-01   # Re(ZDL(2,2))
  2  3     8.75951247E-04   # Re(ZDL(2,3))
  3  1     1.57239604E-04   # Re(ZDL(3,1))
  3  2    -8.75950279E-04   # Re(ZDL(3,2))
  3  3     9.99999604E-01   # Re(ZDL(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(ZDR(1,1))
  1  2     5.75816956E-07   # Re(ZDR(1,2))
  1  3    -3.41712787E-07   # Re(ZDR(1,3))
  2  1    -5.75803149E-07   # Re(ZDR(2,1))
  2  2     9.99999999E-01   # Re(ZDR(2,2))
  2  3     4.04034945E-05   # Re(ZDR(2,3))
  3  1     3.41736052E-07   # Re(ZDR(3,1))
  3  2    -4.04034943E-05   # Re(ZDR(3,2))
  3  3     9.99999999E-01   # Re(ZDR(3,3))
Block UULMIX
  1  1     9.73845795E-01   # Re(ZUL(1,1))
  1  2     2.27200301E-01   # Re(ZUL(1,2))
  1  3     2.09546246E-03   # Re(ZUL(1,3))
  2  1    -2.27095589E-01   # Re(ZUL(2,1))
  2  2     9.73021583E-01   # Re(ZUL(2,2))
  2  3     4.07012648E-02   # Re(ZUL(2,3))
  3  1     7.20840943E-03   # Re(ZUL(3,1))
  3  2    -4.01126259E-02   # Re(ZUL(3,2))
  3  3     9.99169163E-01   # Re(ZUL(3,3))
Block UURMIX
  1  1     1.00000000E+00   # Re(ZUR(1,1))
  1  2    -1.14917713E-08   # Re(ZUR(1,2))
  1  3    -1.42990213E-09   # Re(ZUR(1,3))
  2  1     1.14917555E-08   # Re(ZUR(2,1))
  2  2     1.00000000E+00   # Re(ZUR(2,2))
  2  3    -1.10574958E-05   # Re(ZUR(2,3))
  3  1     1.43002920E-09   # Re(ZUR(3,1))
  3  2     1.10574958E-05   # Re(ZUR(3,2))
  3  3     1.00000000E+00   # Re(ZUR(3,3))
Block VCKM Q= 1.00000000E+03
  1  1     9.73846592E-01   # Re(CKM)(1,1)
  1  2     2.27196483E-01   # Re(CKM)(1,2)
  1  3     2.13864228E-03   # Re(CKM)(1,3)
  2  1    -2.27087408E-01   # Re(CKM)(2,1)
  2  2     9.72988026E-01   # Re(CKM)(2,2)
  2  3     4.15404624E-02   # Re(CKM)(2,3)
  3  1     7.35697364E-03   # Re(CKM)(3,1)
  3  2    -4.09396965E-02   # Re(CKM)(3,2)
  3  3     9.99134534E-01   # Re(CKM)(3,3)
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

   MRSSM2_mass_eigenstates_running m(input);
   slha_io.fill(m);
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   MRSSM2_decays decays_with_HO(m, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0025634108458618405, 5e-12);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              0.00012289267709167259, 2e-13);
   // QED corrections
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
   //                            2.6059181498481999E-003, 5e-15);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.0002754515192858279, 5e-12);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
   //                           0.00066154345019159267, 5e-11);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
                              0.00094083425095728521, 1e-3);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
   //                            7.5383132433569488e-05, 9e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
                              0.00012722814312214181, 1e-3);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VGVG(&m, 0), 0.00040230168141814043, 7e-11);
   // h -> gamma gamma
   // without 2-loop QCD corrections to squark loop
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 8.3519576334971031e-06, 4e-11);
   // with 2-loop QCD corrections to squark loop
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 1.1323516625659973e-05, 4e-11);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVZ(&m, 0), 7.8704345707650875e-06, 5e-11);

   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00029585462066712443, 2e-13);

   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFdFd(&m, 1, 2, 2),
                              27.422744591984291, 5e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_barFuFu(&m, 1, 1, 1),
                              6.5374236623950305e-07, 2e-13);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   MRSSM2_decays decays_without_HO(m, qedqcd, physical_input, flexibledecay_settings);

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0016025829292552107, 2e-13);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              8.4238943320917841e-05, 2e-13);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.00027252978357707359, 5e-12);

   // ------------ loop-induces decays ------------

   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVP(&m, 0), 1.1149676028035005e-05, 4e-11);

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VGVG(&m, 0), 0.00013277227074162971, 8e-11);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVZ(&m, 0), 7.8704345707650875e-06, 5e-11);

   // scalar sgluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_phiO_to_VGVG(&m), 0.057227807158571162, 3e-11);
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_phiO_to_barFuFu(&m, 2, 2), 9.300625612620111e-05, 2e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_phiO_to_SuconjSu(&m, 0, 0), 7.5178588442628946, 1e-16);

   // pseudoscalar sgluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_sigmaO_to_VGVG(&m), 0, 1e-16);
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_sigmaO_to_barFuFu(&m, 2, 2), 9.057611234446413e-05, 2e-14);

   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_Su_to_GluFu(&m, 3, 0), 56.21133546983561, 1e-16);
}

