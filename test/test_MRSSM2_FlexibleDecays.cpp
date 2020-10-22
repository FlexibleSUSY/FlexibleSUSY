
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MRSSM2_FlexibleDecays

#include <boost/test/unit_test.hpp>

#include "test_MRSSM2.hpp"
#include "MRSSM2_two_scale_model.hpp"
#include "MRSSM2_two_scale_spectrum_generator.hpp"
#include "decays/MRSSM2_decays.hpp"
#include "MRSSM2_slha_io.hpp"

using namespace flexiblesusy;

/* BMP3 of arXiv:1410.4791
   considered in Adam Buchner's thesis */

BOOST_AUTO_TEST_CASE( test_MRSSM2_FlexibleDecays )
{

   char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.4.1
     5   MRSSM2
     9   4.13.0
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
  1  1     6.25000000E+06   # md2Input(1,1)
  1  2     0.00000000E+00   # md2Input(1,2)
  1  3     0.00000000E+00   # md2Input(1,3)
  2  1     0.00000000E+00   # md2Input(2,1)
  2  2     6.25000000E+06   # md2Input(2,2)
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
  1  1     1.00000000E+06   # ml2Input(1,1)
  1  2     0.00000000E+00   # ml2Input(1,2)
  1  3     0.00000000E+00   # ml2Input(1,3)
  2  1     0.00000000E+00   # ml2Input(2,1)
  2  2     1.00000000E+06   # ml2Input(2,2)
  2  3     0.00000000E+00   # ml2Input(2,3)
  3  1     0.00000000E+00   # ml2Input(3,1)
  3  2     0.00000000E+00   # ml2Input(3,2)
  3  3     1.00000000E+06   # ml2Input(3,3)
Block MSQ2IN
  1  1     6.25000000E+06   # mq2Input(1,1)
  1  2     0.00000000E+00   # mq2Input(1,2)
  1  3     0.00000000E+00   # mq2Input(1,3)
  2  1     0.00000000E+00   # mq2Input(2,1)
  2  2     6.25000000E+06   # mq2Input(2,2)
  2  3     0.00000000E+00   # mq2Input(2,3)
  3  1     0.00000000E+00   # mq2Input(3,1)
  3  2     0.00000000E+00   # mq2Input(3,2)
  3  3     1.00000000E+06   # mq2Input(3,3)
Block NMSSMRUNIN
    10     4.00000000E+06   # mS2Input
Block MSU2IN
  1  1     6.25000000E+06   # mu2Input(1,1)
  1  2     0.00000000E+00   # mu2Input(1,2)
  1  3     0.00000000E+00   # mu2Input(1,3)
  2  1     0.00000000E+00   # mu2Input(2,1)
  2  2     6.25000000E+06   # mu2Input(2,2)
  2  3     0.00000000E+00   # mu2Input(2,3)
  3  1     0.00000000E+00   # mu2Input(3,1)
  3  2     0.00000000E+00   # mu2Input(3,2)
  3  3     1.00000000E+06   # mu2Input(3,3)
Block gauge Q= 1.00000000E+03
     1     3.62120030E-01   # g1 * 0.7745966692414834
     2     6.43012598E-01   # g2
     3     1.05807099E+00   # g3
Block Yu Q= 1.00000000E+03
  1  1     7.33065972E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.33720117E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.89862945E-01   # Yu(3,3)
Block Yd Q= 1.00000000E+03
  1  1     5.69356767E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     1.24662505E-02   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     5.63512800E-01   # Yd(3,3)
Block Ye Q= 1.00000000E+03
  1  1     1.15181876E-04   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     2.38159700E-02   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     4.00639040E-01   # Ye(3,3)
Block HMIX Q= 1.00000000E+03
     1     0.00000000E+00   # Mu
   101     3.99999831E+04   # BMu
   102     6.17997917E+00   # vd
   103     2.41067730E+02   # vu
   310    -3.44417030E-01   # vT
   201     4.00000047E+02   # MuD
   202     3.99999778E+02   # MuU
   203     0.00000000E+00   # BMuD
   204     0.00000000E+00   # BMuU
   301     1.50000018E-01   # LamSD
   302    -1.49999917E-01   # LamSU
   303    -1.00000010E+00   # LamTD
   304    -1.14999934E+00   # LamTU
Block MSQ2 Q= 1.00000000E+03
  1  1     6.25000008E+06   # mq2(1,1)
  1  2     6.46548635E-04   # mq2(1,2)
  1  3    -1.59044113E-02   # mq2(1,3)
  2  1     6.46548635E-04   # mq2(2,1)
  2  2     6.25000007E+06   # mq2(2,2)
  2  3     8.85814761E-02   # mq2(2,3)
  3  1    -1.59044113E-02   # mq2(3,1)
  3  2     8.85814761E-02   # mq2(3,2)
  3  3     9.99999938E+05   # mq2(3,3)
Block MSE2 Q= 1.00000000E+03
  1  1     9.99999996E+05   # me2(1,1)
  1  2     0.00000000E+00   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     0.00000000E+00   # me2(2,1)
  2  2     9.99999995E+05   # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     9.99999790E+05   # me2(3,3)
Block MSL2 Q= 1.00000000E+03
  1  1     1.00000001E+06   # ml2(1,1)
  1  2     0.00000000E+00   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     0.00000000E+00   # ml2(2,1)
  2  2     1.00000001E+06   # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     9.99999903E+05   # ml2(3,3)
Block MSU2 Q= 1.00000000E+03
  1  1     6.25000008E+06   # mu2(1,1)
  1  2     5.96176233E-12   # mu2(1,2)
  1  3     6.08888029E-08   # mu2(1,3)
  2  1     5.96176927E-12   # mu2(2,1)
  2  2     6.25000008E+06   # mu2(2,2)
  2  3     5.38336802E-04   # mu2(2,3)
  3  1     6.08888029E-08   # mu2(3,1)
  3  2     5.38336802E-04   # mu2(3,2)
  3  3     9.99999268E+05   # mu2(3,3)
Block MSD2 Q= 1.00000000E+03
  1  1     6.25000007E+06   # md2(1,1)
  1  2    -7.86388732E-11   # md2(1,2)
  1  3    -2.11796663E-05   # md2(1,3)
  2  1    -7.86388732E-11   # md2(2,1)
  2  2     6.25000007E+06   # md2(2,2)
  2  3     2.58275994E-03   # md2(2,3)
  3  1    -2.11796663E-05   # md2(3,1)
  3  2     2.58275994E-03   # md2(3,2)
  3  3     1.00000059E+06   # md2(3,3)
Block MSOFT Q= 1.00000000E+03
    21     1.30202577E+06   # mHd2
    22    -2.95541171E+05   # mHu2
   110     8.99999992E+06   # mT2
   111     1.00000017E+06   # moc2
   300     2.50000001E+02   # MDBS
   301     5.00000000E+02   # MDWBT
   302     9.99999929E+02   # MDGoc
    50     4.90000039E+05   # mRd2
    51     9.99999833E+05   # mRu2
Block NMSSMRUN Q= 1.00000000E+03
    10     4.00000000E+06   # mS2
     5    -6.31078301E-02   # vS
Block MASS
        24     8.04935130E+01   # VWm
       404     8.91170992E+02   # SRdp
       403     1.17375227E+03   # SRum
   1000021     1.16847311E+03   # Glu
   3000022     1.15761125E+03   # sigmaO
   3000021     2.32557770E+03   # phiO
   2000024     4.23953307E+02   # Cha2(1)
   2000037     5.64621108E+02   # Cha2(2)
   1000024     4.08278270E+02   # Cha1(1)
   1000037     5.25549522E+02   # Cha1(2)
       401     8.95570270E+02   # Rh(1)
       402     1.16211833E+03   # Rh(2)
   1000012     9.99942682E+02   # Sv(1)
   1000014     1.00175553E+03   # Sv(2)
   1000016     1.00176197E+03   # Sv(3)
   1000022     2.51353830E+02   # Chi(1)
   1000023     4.08048927E+02   # Chi(2)
   1000025     4.20363857E+02   # Chi(3)
   1000035     5.45509248E+02   # Chi(4)
        37     1.23200002E+03   # Hpm(2)
        47     3.01680652E+03   # Hpm(3)
        57     3.17924784E+03   # Hpm(4)
        25     1.24775982E+02   # hh(1)
        35     1.22898763E+03   # hh(2)
        45     2.06068208E+03   # hh(3)
        55     3.17898558E+03   # hh(4)
        36     1.22897644E+03   # Ah(2)
        46     1.99954600E+03   # Ah(3)
        56     3.01660750E+03   # Ah(4)
   1000001     1.05165217E+03   # Sd(1)
   1000003     1.05404121E+03   # Sd(2)
   1000005     2.53612871E+03   # Sd(3)
   2000001     2.53613011E+03   # Sd(4)
   2000003     2.54251640E+03   # Sd(5)
   2000005     2.54252062E+03   # Sd(6)
   1000011     9.98655341E+02   # Se(1)
   1000013     1.00224777E+03   # Se(2)
   1000015     1.00226050E+03   # Se(3)
   2000011     1.00349567E+03   # Se(4)
   2000013     1.00531846E+03   # Se(5)
   2000015     1.00532492E+03   # Se(6)
   1000002     1.05710992E+03   # Su(1)
   1000004     1.06420319E+03   # Su(2)
   1000006     2.53648695E+03   # Su(3)
   2000002     2.53648715E+03   # Su(4)
   2000004     2.54139480E+03   # Su(5)
   2000006     2.54139875E+03   # Su(6)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.29033878E-04   # Fe(1)
        13     1.07278707E-01   # Fe(2)
        15     1.78529243E+00   # Fe(3)
         1     4.43123064E-03   # Fd(1)
         3     8.77332488E-02   # Fd(2)
         5     3.45634038E+00   # Fd(3)
         2     2.29355408E-03   # Fu(1)
         4     8.47491831E-01   # Fu(2)
         6     1.77434227E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.11601029E+01   # VZ
Block U1MIX
  1  1     2.75021091E-03   # Re(UM1(1,1))
  1  2     9.99996218E-01   # Re(UM1(1,2))
  2  1     9.99996218E-01   # Re(UM1(2,1))
  2  2    -2.75021091E-03   # Re(UM1(2,2))
Block U2MIX
  1  1     4.64965885E-01   # Re(UM2(1,1))
  1  2    -8.85328598E-01   # Re(UM2(1,2))
  2  1     8.85328598E-01   # Re(UM2(2,1))
  2  2     4.64965885E-01   # Re(UM2(2,2))
Block V1MIX
  1  1     1.08747242E-02   # Re(UP1(1,1))
  1  2     9.99940868E-01   # Re(UP1(1,2))
  2  1     9.99940868E-01   # Re(UP1(2,1))
  2  2    -1.08747242E-02   # Re(UP1(2,2))
Block V2MIX
  1  1     1.66010343E-01   # Re(UP2(1,1))
  1  2     9.86124011E-01   # Re(UP2(1,2))
  2  1     9.86124011E-01   # Re(UP2(2,1))
  2  2    -1.66010343E-01   # Re(UP2(2,2))
Block PSEUDOSCALARMIX
  1  1    -2.50525910E-02   # ZA(1,1)
  1  2     9.99686126E-01   # ZA(1,2)
  1  3    -6.90965974E-05   # ZA(1,3)
  1  4     1.11592680E-04   # ZA(1,4)
  2  1     9.99686135E-01   # ZA(2,1)
  2  2     2.50525911E-02   # ZA(2,2)
  2  3    -1.04941819E-06   # ZA(2,3)
  2  4     1.12201699E-06   # ZA(2,4)
  3  1     6.80730364E-07   # ZA(3,1)
  3  2    -6.90192312E-05   # ZA(3,2)
  3  3    -9.99999728E-01   # ZA(3,3)
  3  4    -7.34427429E-04   # ZA(3,4)
  4  1     1.67452135E-06   # ZA(4,1)
  4  2    -1.11636483E-04   # ZA(4,2)
  4  3    -7.34419720E-04   # ZA(4,3)
  4  4     9.99999724E-01   # ZA(4,4)
Block DSQMIX
  1  1     0.00000000E+00   # ZD(1,1)
  1  2     0.00000000E+00   # ZD(1,2)
  1  3     0.00000000E+00   # ZD(1,3)
  1  4     1.75745044E-07   # ZD(1,4)
  1  5    -2.14231084E-05   # ZD(1,5)
  1  6     1.00000000E+00   # ZD(1,6)
  2  1     1.06097390E-04   # ZD(2,1)
  2  2    -5.90875452E-04   # ZD(2,2)
  2  3     9.99999820E-01   # ZD(2,3)
  2  4     0.00000000E+00   # ZD(2,4)
  2  5     0.00000000E+00   # ZD(2,5)
  2  6     0.00000000E+00   # ZD(2,6)
  3  1    -1.97712110E-01   # ZD(3,1)
  3  2     9.80259946E-01   # ZD(3,2)
  3  3     6.00188385E-04   # ZD(3,3)
  3  4     0.00000000E+00   # ZD(3,4)
  3  5     0.00000000E+00   # ZD(3,5)
  3  6     0.00000000E+00   # ZD(3,6)
  4  1    -9.80260124E-01   # ZD(4,1)
  4  2    -1.97712138E-01   # ZD(4,2)
  4  3    -1.28202108E-05   # ZD(4,3)
  4  4    -0.00000000E+00   # ZD(4,4)
  4  5    -0.00000000E+00   # ZD(4,5)
  4  6    -0.00000000E+00   # ZD(4,6)
  5  1    -0.00000000E+00   # ZD(5,1)
  5  2    -0.00000000E+00   # ZD(5,2)
  5  3    -0.00000000E+00   # ZD(5,3)
  5  4    -1.00000000E+00   # ZD(5,4)
  5  5    -5.78329000E-07   # ZD(5,5)
  5  6     1.75732655E-07   # ZD(5,6)
  6  1     0.00000000E+00   # ZD(6,1)
  6  2     0.00000000E+00   # ZD(6,2)
  6  3     0.00000000E+00   # ZD(6,3)
  6  4     5.78325235E-07   # ZD(6,4)
  6  5    -1.00000000E+00   # ZD(6,5)
  6  6    -2.14231085E-05   # ZD(6,6)
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
  1  1    -2.60150592E-02   # ZH(1,1)
  1  2    -9.99657406E-01   # ZH(1,2)
  1  3     1.56979855E-04   # ZH(1,3)
  1  4     2.87442402E-03   # ZH(1,4)
  2  1     9.99661549E-01   # ZH(2,1)
  2  2    -2.60147782E-02   # ZH(2,2)
  2  3     1.66647476E-05   # ZH(2,3)
  2  4     1.34305770E-04   # ZH(2,4)
  3  1     1.26120035E-05   # ZH(3,1)
  3  2    -1.59136493E-04   # ZH(3,2)
  3  3    -9.99999797E-01   # ZH(3,3)
  3  4    -6.17184141E-04   # ZH(3,4)
  4  1     5.94744765E-05   # ZH(4,1)
  4  2    -2.87684744E-03   # ZH(4,2)
  4  3     6.17640155E-04   # ZH(4,3)
  4  4    -9.99995669E-01   # ZH(4,4)
Block RHMIX
  1  1    -9.99999638E-01   # ZHR(1,1)
  1  2     8.50495379E-04   # ZHR(1,2)
  2  1    -8.50495379E-04   # ZHR(2,1)
  2  2    -9.99999638E-01   # ZHR(2,2)
Block N1MIX
  1  1     9.93910506E-01   # Re(ZN1(1,1))
  1  2     3.53299104E-02   # Re(ZN1(1,2))
  1  3    -2.67585189E-03   # Re(ZN1(1,3))
  1  4    -1.04338600E-01   # Re(ZN1(1,4))
  2  1    -8.11104643E-04   # Re(ZN1(2,1))
  2  2     7.75211685E-04   # Re(ZN1(2,2))
  2  3    -9.99834143E-01   # Re(ZN1(2,3))
  2  4     1.81776526E-02   # Re(ZN1(2,4))
  3  1     1.10102244E-01   # Re(ZN1(3,1))
  3  2    -3.55647997E-01   # Re(ZN1(3,2))
  3  3     1.65059610E-02   # Re(ZN1(3,3))
  3  4     9.27965275E-01   # Re(ZN1(3,4))
  4  1     4.32945548E-03   # Re(ZN1(4,1))
  4  2     9.33951658E-01   # Re(ZN1(4,2))
  4  3     7.21657556E-03   # Re(ZN1(4,3))
  4  4     3.57300263E-01   # Re(ZN1(4,4))
Block N2MIX
  1  1     9.99846857E-01   # Re(ZN2(1,1))
  1  2     1.59398355E-02   # Re(ZN2(1,2))
  1  3     1.94223332E-04   # Re(ZN2(1,3))
  1  4    -7.22129694E-03   # Re(ZN2(1,4))
  2  1    -5.47010467E-05   # Re(ZN2(2,1))
  2  2    -5.30320738E-04   # Re(ZN2(2,2))
  2  3     9.99835186E-01   # Re(ZN2(2,3))
  2  4     1.81470763E-02   # Re(ZN2(2,4))
  3  1     9.33126946E-03   # Re(ZN2(3,1))
  3  2    -1.36389336E-01   # Re(ZN2(3,2))
  3  3    -1.80485104E-02   # Re(ZN2(3,3))
  3  4     9.90446933E-01   # Re(ZN2(3,4))
  4  1    -1.48049865E-02   # Re(ZN2(4,1))
  4  2     9.90526925E-01   # Re(ZN2(4,2))
  4  3    -1.95298770E-03   # Re(ZN2(4,3))
  4  4     1.36504245E-01   # Re(ZN2(4,4))
Block CHARGEMIX
  1  1     2.51132077E-02   # ZP(1,1)
  1  2    -9.99681236E-01   # ZP(1,2)
  1  3    -1.91329206E-03   # ZP(1,3)
  1  4    -1.75875021E-03   # ZP(1,4)
  2  1    -9.99684612E-01   # ZP(2,1)
  2  2    -2.51129839E-02   # ZP(2,2)
  2  3    -8.47447660E-05   # ZP(2,3)
  2  4    -8.32594559E-05   # ZP(2,4)
  3  1     9.76024908E-07   # ZP(3,1)
  3  2    -1.44110284E-04   # ZP(3,2)
  3  3     7.16525947E-01   # ZP(3,3)
  3  4    -6.97560425E-01   # ZP(3,4)
  4  1     5.35704486E-05   # ZP(4,1)
  4  2     2.59699135E-03   # ZP(4,2)
  4  3    -6.97557811E-01   # ZP(4,3)
  4  4    -7.16523798E-01   # ZP(4,4)
Block USQMIX
  1  1     0.00000000E+00   # ZU(1,1)
  1  2     0.00000000E+00   # ZU(1,2)
  1  3     0.00000000E+00   # ZU(1,3)
  1  4     3.91334119E-10   # ZU(1,4)
  1  5     3.46105492E-06   # ZU(1,5)
  1  6     1.00000000E+00   # ZU(1,6)
  2  1     6.88673900E-05   # ZU(2,1)
  2  2    -3.83703467E-04   # ZU(2,2)
  2  3     9.99999924E-01   # ZU(2,3)
  2  4     0.00000000E+00   # ZU(2,4)
  2  5     0.00000000E+00   # ZU(2,5)
  2  6     0.00000000E+00   # ZU(2,6)
  3  1     9.87267869E-01   # ZU(3,1)
  3  2     1.59066508E-01   # ZU(3,2)
  3  3    -6.95619121E-06   # ZU(3,3)
  3  4    -0.00000000E+00   # ZU(3,4)
  3  5    -0.00000000E+00   # ZU(3,5)
  3  6    -0.00000000E+00   # ZU(3,6)
  4  1     1.59066493E-01   # ZU(4,1)
  4  2    -9.87267795E-01   # ZU(4,2)
  4  3    -3.89772599E-04   # ZU(4,3)
  4  4     0.00000000E+00   # ZU(4,4)
  4  5     0.00000000E+00   # ZU(4,5)
  4  6     0.00000000E+00   # ZU(4,6)
  5  1     0.00000000E+00   # ZU(5,1)
  5  2     0.00000000E+00   # ZU(5,2)
  5  3     0.00000000E+00   # ZU(5,3)
  5  4    -1.00000000E+00   # ZU(5,4)
  5  5    -1.95514893E-07   # ZU(5,5)
  5  6     3.92010806E-10   # ZU(5,6)
  6  1     0.00000000E+00   # ZU(6,1)
  6  2     0.00000000E+00   # ZU(6,2)
  6  3     0.00000000E+00   # ZU(6,3)
  6  4     1.95514895E-07   # ZU(6,4)
  6  5    -1.00000000E+00   # ZU(6,5)
  6  6     3.46105492E-06   # ZU(6,6)
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
  1  1     9.99999987E-01   # Re(ZDL(1,1))
  1  2     6.30436594E-06   # Re(ZDL(1,2))
  1  3    -1.59439405E-04   # Re(ZDL(1,3))
  2  1    -6.16274368E-06   # Re(ZDL(2,1))
  2  2     9.99999605E-01   # Re(ZDL(2,2))
  2  3     8.88236153E-04   # Re(ZDL(2,3))
  3  1     1.59444941E-04   # Re(ZDL(3,1))
  3  2    -8.88235159E-04   # Re(ZDL(3,2))
  3  3     9.99999593E-01   # Re(ZDL(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(ZDR(1,1))
  1  2     5.82600624E-07   # Re(ZDR(1,2))
  1  3    -3.46055997E-07   # Re(ZDR(1,3))
  2  1    -5.82586462E-07   # Re(ZDR(2,1))
  2  2     9.99999999E-01   # Re(ZDR(2,2))
  2  3     4.09220373E-05   # Re(ZDR(2,3))
  3  1     3.46079838E-07   # Re(ZDR(3,1))
  3  2    -4.09220371E-05   # Re(ZDR(3,2))
  3  3     9.99999999E-01   # Re(ZDR(3,3))
Block UULMIX
  1  1     9.73845815E-01   # Re(ZUL(1,1))
  1  2     2.27200171E-01   # Re(ZUL(1,2))
  1  3     2.10036017E-03   # Re(ZUL(1,3))
  2  1    -2.27094969E-01   # Re(ZUL(2,1))
  2  2     9.73017747E-01   # Re(ZUL(2,2))
  2  3     4.07963018E-02   # Re(ZUL(2,3))
  3  1     7.22523903E-03   # Re(ZUL(3,1))
  3  2    -4.02062890E-02   # Re(ZUL(3,2))
  3  3     9.99165277E-01   # Re(ZUL(3,3))
Block UURMIX
  1  1     1.00000000E+00   # Re(ZUR(1,1))
  1  2    -1.11482059E-08   # Re(ZUR(1,2))
  1  3    -1.26334348E-09   # Re(ZUR(1,3))
  2  1     1.11481937E-08   # Re(ZUR(2,1))
  2  2     1.00000000E+00   # Re(ZUR(2,2))
  2  3    -9.71199909E-06   # Re(ZUR(2,3))
  3  1     1.26345175E-09   # Re(ZUR(3,1))
  3  2     9.71199909E-06   # Re(ZUR(3,2))
  3  3     1.00000000E+00   # Re(ZUR(3,3))
Block VCKM Q= 1.00000000E+03
  1  1     9.73846600E-01   # Re(CKM)(1,1)
  1  2     2.27196452E-01   # Re(CKM)(1,2)
  1  3     2.13826417E-03   # Re(CKM)(1,3)
  2  1    -2.27087415E-01   # Re(CKM)(2,1)
  2  2     9.72988337E-01   # Re(CKM)(2,2)
  2  3     4.15331325E-02   # Re(CKM)(2,3)
  3  1     7.35567425E-03   # Re(CKM)(3,1)
  3  2    -4.09324727E-02   # Re(CKM)(3,2)
  3  3     9.99134839E-01   # Re(CKM)(3,3)
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
Block FlexibleSUSYLowEnergy Q= 1.00000000E+03
    21     0.00000000E+00   # Delta(g-2)_muon/2 FlexibleSUSY
    23     0.00000000E+00   # electric dipole moment of Fe(0) [1/GeV]
    24     0.00000000E+00   # electric dipole moment of Fe(1) [1/GeV]
    25     0.00000000E+00   # electric dipole moment of Fe(2) [1/GeV]
    26     0.00000000E+00   # BR(Fe1 -> Fe0 VP)
)";

   std::stringstream istr(slha_input);

   MRSSM2_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   MRSSM2_input_parameters input;
   Spectrum_generator_settings settings;

   // extract the input parameters from spectrum string
   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   MRSSM2_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   MRSSM2_slha m = std::get<0>(spectrum_generator.get_models_slha());

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   MRSSM2_decays decays_with_HO(m, qedqcd, SM_higher_order_corrections::enable);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0018866791466703291, 5e-12);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              9.9580460952386956e-05, 4e-14);
   // QED corrections
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
   //                            2.6059181498481999E-003, 5e-15);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.00026946060398832237, 5e-12);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
   //                           0.00066154345019159267, 5e-11);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_conjVWmVWm(&m, 0),
                              0.00073836837044127768, 2e-4);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
   //                            7.5383132433569488e-05, 9e-12);
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VZVZ(&m, 0),
                              9.638232014475222e-05, 1e-4);

   // ------------ loop-induces decays ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VGVG(&m, 0), 0.00020682880248722835, 7e-11);
   // h -> gamma gamma
   // without 2-loop QCD corrections to squark loop
   // BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 8.3519576334971031e-06, 4e-11);
   // with 2-loop QCD corrections to squark loop
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVP(&m, 0), 1.2159638473406977e-05, 4e-11);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_hh_to_VPVZ(&m, 0), 6.9440628969348722e-06, 5e-11);

   BOOST_CHECK_CLOSE_FRACTION(decays_with_HO.partial_width_Ah_to_VGVG(&m, 1), 0.00029502623532270532, 8e-14);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   MRSSM2_decays decays_without_HO(m, qedqcd, SM_higher_order_corrections::disable);

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFdFd(&m, 0, 2, 2),
                              0.0015852320624501718, 4e-14);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFuFu(&m, 0, 1, 1),
                              8.2876032145842477e-05, 3e-14);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_barFeFe(&m, 0, 2, 2),
                              0.00026660324954596258, 5e-12);

   // ------------ loop-induces decays ------------

   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VPVP(&m, 0), 1.1943830812277772e-05, 4e-11);

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_without_HO.partial_width_hh_to_VGVG(&m, 0), 0.00012423136936565911, 7e-11);
}
