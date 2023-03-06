#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_unitarity

#include <boost/test/unit_test.hpp>

#include "CMSSM_slha_io.hpp"
#include "lowe.h"
#include "error.hpp"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "CMSSM_input_parameters.hpp"
#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_model_slha.hpp"
#include "CMSSM_slha_io.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "CMSSM_unitarity.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_FlexibleDecay )
{

   char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.7.1
     5   CMSSM
     9   4.15.1
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
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
     1     9.00000000E+02   # m0
     2     1.20000000E+03   # m12
     3     2.00000000E+01   # TanBeta
     4     1.00000000E+00   # SignMu
     5    -3.00000000E+03   # Azero
Block gauge Q= 1.69653267E+03
     1     3.63323632E-01   # g1 * 0.7745966692414834
     2     6.37468910E-01   # g2
     3     1.02444596E+00   # g3
Block Yu Q= 1.69653267E+03
  1  1     7.15847023E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.27350609E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.29803780E-01   # Yu(3,3)
Block Yd Q= 1.69653267E+03
  1  1     2.74872963E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     6.01828620E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     2.44611829E-01   # Yd(3,3)
Block Ye Q= 1.69653267E+03
  1  1     5.86902388E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     1.21352568E-02   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     2.03920370E-01   # Ye(3,3)
Block Te Q= 1.69653267E+03
  1  1    -2.07352238E-01   # TYe(1,1)
  1  2     0.00000000E+00   # TYe(1,2)
  1  3     0.00000000E+00   # TYe(1,3)
  2  1     0.00000000E+00   # TYe(2,1)
  2  2    -4.28686768E+01   # TYe(2,2)
  2  3     0.00000000E+00   # TYe(2,3)
  3  1     0.00000000E+00   # TYe(3,1)
  3  2     0.00000000E+00   # TYe(3,2)
  3  3    -6.96225137E+02   # TYe(3,3)
Block Td Q= 1.69653267E+03
  1  1    -1.61582711E+00   # TYd(1,1)
  1  2     0.00000000E+00   # TYd(1,2)
  1  3     0.00000000E+00   # TYd(1,3)
  2  1     0.00000000E+00   # TYd(2,1)
  2  2    -3.53776793E+01   # TYd(2,2)
  2  3     0.00000000E+00   # TYd(2,3)
  3  1     0.00000000E+00   # TYd(3,1)
  3  2     0.00000000E+00   # TYd(3,2)
  3  3    -1.28836088E+03   # TYd(3,3)
Block Tu Q= 1.69653267E+03
  1  1    -3.26797947E-02   # TYu(1,1)
  1  2     0.00000000E+00   # TYu(1,2)
  1  3     0.00000000E+00   # TYu(1,3)
  2  1     0.00000000E+00   # TYu(2,1)
  2  2    -1.49440255E+01   # TYu(2,2)
  2  3     0.00000000E+00   # TYu(2,3)
  3  1     0.00000000E+00   # TYu(3,1)
  3  2     0.00000000E+00   # TYu(3,2)
  3  3    -2.54167207E+03   # TYu(3,3)
Block MSQ2 Q= 1.69653267E+03
  1  1     5.87420303E+06   # mq2(1,1)
  1  2     0.00000000E+00   # mq2(1,2)
  1  3     0.00000000E+00   # mq2(1,3)
  2  1     0.00000000E+00   # mq2(2,1)
  2  2     5.87400263E+06   # mq2(2,2)
  2  3     0.00000000E+00   # mq2(2,3)
  3  1     0.00000000E+00   # mq2(3,1)
  3  2     0.00000000E+00   # mq2(3,2)
  3  3     3.99172485E+06   # mq2(3,3)
Block MSE2 Q= 1.69653267E+03
  1  1     1.00109925E+06   # me2(1,1)
  1  2     0.00000000E+00   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     0.00000000E+00   # me2(2,1)
  2  2     9.99900678E+05   # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     6.72283276E+05   # me2(3,3)
Block MSL2 Q= 1.69653267E+03
  1  1     1.41541159E+06   # ml2(1,1)
  1  2     0.00000000E+00   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     0.00000000E+00   # ml2(2,1)
  2  2     1.41481774E+06   # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     1.25262789E+06   # ml2(3,3)
Block MSU2 Q= 1.69653267E+03
  1  1     5.44637416E+06   # mu2(1,1)
  1  2     0.00000000E+00   # mu2(1,2)
  1  3     0.00000000E+00   # mu2(1,3)
  2  1     0.00000000E+00   # mu2(2,1)
  2  2     5.44629860E+06   # mu2(2,2)
  2  3     0.00000000E+00   # mu2(2,3)
  3  1     0.00000000E+00   # mu2(3,1)
  3  2     0.00000000E+00   # mu2(3,2)
  3  3     2.09801547E+06   # mu2(3,3)
Block MSD2 Q= 1.69653267E+03
  1  1     5.39397190E+06   # md2(1,1)
  1  2     0.00000000E+00   # md2(1,2)
  1  3     0.00000000E+00   # md2(1,3)
  2  1     0.00000000E+00   # md2(2,1)
  2  2     5.39364200E+06   # md2(2,2)
  2  3     0.00000000E+00   # md2(2,3)
  3  1     0.00000000E+00   # md2(3,1)
  3  2     0.00000000E+00   # md2(3,2)
  3  3     4.90425043E+06   # md2(3,3)
Block Phases Q= 1.69653267E+03
     1     1.00000000E+00   # Re(PhaseGlu)
Block IMPhases Q= 1.69653267E+03
     1     0.00000000E+00   # Im(PhaseGlu)
Block MASS
   1000021     2.59926541E+03   # Glu
        24     8.03552325E+01   # VWm
   1000024     9.82646494E+02   # Cha(1)
   1000037     2.00134175E+03   # Cha(2)
        25     1.24953773E+02   # hh(1)
        35     2.12829258E+03   # hh(2)
        37     2.13036038E+03   # Hpm(2)
        36     2.12844754E+03   # Ah(2)
   1000012     1.12323577E+03   # Sv(1)
   1000014     1.19517456E+03   # Sv(2)
   1000016     1.19542931E+03   # Sv(3)
   1000022     5.19516874E+02   # Chi(1)
   1000023     9.82484228E+02   # Chi(2)
   1000025    -1.99818327E+03   # Chi(3)
   1000035     2.00040731E+03   # Chi(4)
   1000001     2.03470543E+03   # Sd(1)
   1000003     2.27294662E+03   # Sd(2)
   1000005     2.37930063E+03   # Sd(3)
   2000001     2.37937291E+03   # Sd(4)
   2000003     2.48396840E+03   # Sd(5)
   2000005     2.48400742E+03   # Sd(6)
   1000011     8.16470659E+02   # Se(1)
   1000013     1.00467447E+03   # Se(2)
   1000015     1.00530529E+03   # Se(3)
   2000011     1.13061565E+03   # Se(4)
   2000013     1.19806402E+03   # Se(5)
   2000015     1.19829755E+03   # Se(6)
   1000002     1.45861649E+03   # Su(1)
   1000004     2.06630028E+03   # Su(2)
   1000006     2.38988052E+03   # Su(3)
   2000002     2.38990051E+03   # Su(4)
   2000004     2.48284741E+03   # Su(5)
   2000006     2.48288543E+03   # Su(6)
Block UMIX
  1  1     9.97258631E-01   # Re(UM(1,1))
  1  2    -7.39947425E-02   # Re(UM(1,2))
  2  1     7.39947425E-02   # Re(UM(2,1))
  2  2     9.97258631E-01   # Re(UM(2,2))
Block VMIX
  1  1     9.99231772E-01   # Re(UP(1,1))
  1  2    -3.91901127E-02   # Re(UP(1,2))
  2  1     3.91901127E-02   # Re(UP(2,1))
  2  2     9.99231772E-01   # Re(UP(2,2))
Block PSEUDOSCALARMIX
  1  1    -5.55207927E-02   # ZA(1,1)
  1  2     9.98457531E-01   # ZA(1,2)
  2  1     9.98457531E-01   # ZA(2,1)
  2  2     5.55207927E-02   # ZA(2,2)
Block DSQMIX
  1  1     0.00000000E+00   # ZD(1,1)
  1  2     3.52200044E-16   # ZD(1,2)
  1  3    -9.95819693E-01   # ZD(1,3)
  1  4     0.00000000E+00   # ZD(1,4)
  1  5    -4.39823147E-16   # ZD(1,5)
  1  6    -9.13407825E-02   # ZD(1,6)
  2  1     0.00000000E+00   # ZD(2,1)
  2  2    -6.46680472E-17   # ZD(2,2)
  2  3     9.13407825E-02   # ZD(2,3)
  2  4     0.00000000E+00   # ZD(2,4)
  2  5     1.53271070E-16   # ZD(2,5)
  2  6    -9.95819693E-01   # ZD(2,6)
  3  1    -0.00000000E+00   # ZD(3,1)
  3  2     4.90894302E-03   # ZD(3,2)
  3  3    -5.73300987E-16   # ZD(3,3)
  3  4    -0.00000000E+00   # ZD(3,4)
  3  5     9.99987951E-01   # ZD(3,5)
  3  6     1.12947362E-16   # ZD(3,6)
  4  1     2.24272845E-04   # ZD(4,1)
  4  2     0.00000000E+00   # ZD(4,2)
  4  3     0.00000000E+00   # ZD(4,3)
  4  4     9.99999975E-01   # ZD(4,4)
  4  5     0.00000000E+00   # ZD(4,5)
  4  6     0.00000000E+00   # ZD(4,6)
  5  1     0.00000000E+00   # ZD(5,1)
  5  2    -9.99987951E-01   # ZD(5,2)
  5  3    -3.59453202E-16   # ZD(5,3)
  5  4     0.00000000E+00   # ZD(5,4)
  5  5     4.90894302E-03   # ZD(5,5)
  5  6     3.27823345E-17   # ZD(5,6)
  6  1     9.99999975E-01   # ZD(6,1)
  6  2     0.00000000E+00   # ZD(6,2)
  6  3     0.00000000E+00   # ZD(6,3)
  6  4    -2.24272845E-04   # ZD(6,4)
  6  5     0.00000000E+00   # ZD(6,5)
  6  6     0.00000000E+00   # ZD(6,6)
Block SELMIX
  1  1    -0.00000000E+00   # ZE(1,1)
  1  2     2.63909013E-17   # ZE(1,2)
  1  3    -1.26313485E-01   # ZE(1,3)
  1  4    -0.00000000E+00   # ZE(1,4)
  1  5     3.78644394E-17   # ZE(1,5)
  1  6    -9.91990375E-01   # ZE(1,6)
  2  1     0.00000000E+00   # ZE(2,1)
  2  2    -1.07602944E-02   # ZE(2,2)
  2  3     1.27559077E-16   # ZE(2,3)
  2  4     0.00000000E+00   # ZE(2,4)
  2  5    -9.99942106E-01   # ZE(2,5)
  2  6    -4.24138649E-17   # ZE(2,6)
  3  1     5.21251682E-05   # ZE(3,1)
  3  2     0.00000000E+00   # ZE(3,2)
  3  3     0.00000000E+00   # ZE(3,3)
  3  4     9.99999999E-01   # ZE(3,4)
  3  5     0.00000000E+00   # ZE(3,5)
  3  6     0.00000000E+00   # ZE(3,6)
  4  1    -0.00000000E+00   # ZE(4,1)
  4  2    -1.00482522E-15   # ZE(4,2)
  4  3     9.91990375E-01   # ZE(4,3)
  4  4    -0.00000000E+00   # ZE(4,4)
  4  5     1.46982243E-16   # ZE(4,5)
  4  6    -1.26313485E-01   # ZE(4,6)
  5  1     0.00000000E+00   # ZE(5,1)
  5  2    -9.99942106E-01   # ZE(5,2)
  5  3    -1.00154103E-15   # ZE(5,3)
  5  4     0.00000000E+00   # ZE(5,4)
  5  5     1.07602944E-02   # ZE(5,5)
  5  6     1.01205699E-16   # ZE(5,6)
  6  1     9.99999999E-01   # ZE(6,1)
  6  2     0.00000000E+00   # ZE(6,2)
  6  3     0.00000000E+00   # ZE(6,3)
  6  4    -5.21251682E-05   # ZE(6,4)
  6  5     0.00000000E+00   # ZE(6,5)
  6  6     0.00000000E+00   # ZE(6,6)
Block SCALARMIX
  1  1     5.19950069E-02   # ZH(1,1)
  1  2     9.98647345E-01   # ZH(1,2)
  2  1     9.98647345E-01   # ZH(2,1)
  2  2    -5.19950069E-02   # ZH(2,2)
Block NMIX
  1  1    -9.99669997E-01   # Re(ZN(1,1))
  1  2     1.59533855E-03   # Re(ZN(1,2))
  1  3    -2.44551895E-02   # Re(ZN(1,3))
  1  4     7.70033491E-03   # Re(ZN(1,4))
  2  1     3.07991927E-03   # Re(ZN(2,1))
  2  2     9.98257674E-01   # Re(ZN(2,2))
  2  3    -5.21267359E-02   # Re(ZN(2,3))
  2  4     2.74760422E-02   # Re(ZN(2,4))
  3  1     1.18114404E-02   # Re(ZN(3,1))
  3  2    -1.74702875E-02   # Re(ZN(3,2))
  3  3    -7.06656457E-01   # Re(ZN(3,3))
  3  4    -7.07242483E-01   # Re(ZN(3,4))
  4  1     2.26031089E-02   # Re(ZN(4,1))
  4  2    -5.63370158E-02   # Re(ZN(4,2))
  4  3    -7.05210180E-01   # Re(ZN(4,3))
  4  4     7.06394962E-01   # Re(ZN(4,4))
Block CHARGEMIX
  1  1    -5.55277431E-02   # ZP(1,1)
  1  2     9.98457145E-01   # ZP(1,2)
  2  1     9.98457145E-01   # ZP(2,1)
  2  2     5.55277431E-02   # ZP(2,2)
Block USQMIX
  1  1     0.00000000E+00   # ZU(1,1)
  1  2     0.00000000E+00   # ZU(1,2)
  1  3     2.21507342E-01   # ZU(1,3)
  1  4     0.00000000E+00   # ZU(1,4)
  1  5     0.00000000E+00   # ZU(1,5)
  1  6     9.75158704E-01   # ZU(1,6)
  2  1     0.00000000E+00   # ZU(2,1)
  2  2     0.00000000E+00   # ZU(2,2)
  2  3     9.75158704E-01   # ZU(2,3)
  2  4     0.00000000E+00   # ZU(2,4)
  2  5     0.00000000E+00   # ZU(2,5)
  2  6    -2.21507342E-01   # ZU(2,6)
  3  1     0.00000000E+00   # ZU(3,1)
  3  2    -6.05414320E-03   # ZU(3,2)
  3  3     0.00000000E+00   # ZU(3,3)
  3  4     0.00000000E+00   # ZU(3,4)
  3  5    -9.99981674E-01   # ZU(3,5)
  3  6     0.00000000E+00   # ZU(3,6)
  4  1     1.32363003E-05   # ZU(4,1)
  4  2     0.00000000E+00   # ZU(4,2)
  4  3     0.00000000E+00   # ZU(4,3)
  4  4     1.00000000E+00   # ZU(4,4)
  4  5     0.00000000E+00   # ZU(4,5)
  4  6     0.00000000E+00   # ZU(4,6)
  5  1     0.00000000E+00   # ZU(5,1)
  5  2    -9.99981674E-01   # ZU(5,2)
  5  3     0.00000000E+00   # ZU(5,3)
  5  4     0.00000000E+00   # ZU(5,4)
  5  5     6.05414320E-03   # ZU(5,5)
  5  6     0.00000000E+00   # ZU(5,6)
  6  1     1.00000000E+00   # ZU(6,1)
  6  2     0.00000000E+00   # ZU(6,2)
  6  3     0.00000000E+00   # ZU(6,3)
  6  4    -1.32363003E-05   # ZU(6,4)
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
Block FlexibleSUSYOutput
     0     1.43018463E+16   # HighScale
     1     1.69653267E+03   # SUSYScale
     2     9.11876000E+01   # LowScale
Block FlexibleSUSYLowEnergy Q= 1.69653267E+03
    20     3.51727141E-15   # Delta(g-2)/2 of Fe(1) (calculated with FlexibleSUSY)
    21     1.59810752E-10   # Delta(g-2)/2 of Fe(2) (calculated with FlexibleSUSY)
    22     5.52622749E-08   # Delta(g-2)/2 of Fe(3) (calculated with FlexibleSUSY)
Block ALPHA
          -5.20184634E-02   # ArcSin(Pole(ZH(2,2)))
Block HMIX Q= 1.69653267E+03
     1     2.00326649E+03   # Mu
     2     1.92778330E+01   # vu/vd
     3     2.43366355E+02   # Sqrt(Sqr(vd) + Sqr(vu))
     4     4.71979738E+06   # Sqr(MAh(2))
   101     2.44173262E+05   # BMu
   102     1.26072047E+01   # vd
   103     2.43039588E+02   # vu
Block Au Q= 1.69653267E+03
  1  1    -4.56519251E+03   # TYu(1,1)/Yu(1,1)
  2  2    -4.56514364E+03   # TYu(2,2)/Yu(2,2)
  3  3    -3.06297963E+03   # TYu(3,3)/Yu(3,3)
Block Ad Q= 1.69653267E+03
  1  1    -5.87845053E+03   # TYd(1,1)/Yd(1,1)
  2  2    -5.87836439E+03   # TYd(2,2)/Yd(2,2)
  3  3    -5.26696067E+03   # TYd(3,3)/Yd(3,3)
Block Ae Q= 1.69653267E+03
  1  1    -3.53299359E+03   # TYe(1,1)/Ye(1,1)
  2  2    -3.53257268E+03   # TYe(2,2)/Ye(2,2)
  3  3    -3.41420103E+03   # TYe(3,3)/Ye(3,3)
Block MSOFT Q= 1.69653267E+03
     1     5.23791145E+02   # MassB
     2     9.52920422E+02   # MassWB
     3     2.54069787E+03   # MassG
    21     4.79434582E+05   # mHd2
    22    -3.97106861E+06   # mHu2
    31     1.18971072E+03   # SignedAbsSqrt(ml2(1,1))
    32     1.18946111E+03   # SignedAbsSqrt(ml2(2,2))
    33     1.11920860E+03   # SignedAbsSqrt(ml2(3,3))
    34     1.00054947E+03   # SignedAbsSqrt(me2(1,1))
    35     9.99950338E+02   # SignedAbsSqrt(me2(2,2))
    36     8.19928824E+02   # SignedAbsSqrt(me2(3,3))
    41     2.42367552E+03   # SignedAbsSqrt(mq2(1,1))
    42     2.42363418E+03   # SignedAbsSqrt(mq2(2,2))
    43     1.99793014E+03   # SignedAbsSqrt(mq2(3,3))
    44     2.33374681E+03   # SignedAbsSqrt(mu2(1,1))
    45     2.33373062E+03   # SignedAbsSqrt(mu2(2,2))
    46     1.44845278E+03   # SignedAbsSqrt(mu2(3,3))
    47     2.32249260E+03   # SignedAbsSqrt(md2(1,1))
    48     2.32242158E+03   # SignedAbsSqrt(md2(2,2))
    49     2.21455423E+03   # SignedAbsSqrt(md2(3,3))
)";

   std::stringstream istr(slha_input);

   CMSSM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   CMSSM_input_parameters input;
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

   CMSSM_slha m(input);
   slha_io.fill(m);
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   static constexpr int size = 36;

   Eigen::SparseMatrix<double> ref_matrix(size, size);
   ref_matrix.insert(0, 0)   = -0.0046229614555045135;
   ref_matrix.insert(0, 8)   = -0.0015409312401939398;
   ref_matrix.insert(0, 18)  = -0.0021792849295675151;
   ref_matrix.insert(0, 24)  = -0.0015391607421558148;
   ref_matrix.insert(0, 29)  =  0.002191044670651851;
   ref_matrix.insert(1, 1)   = -0.0030818624803878791;
   ref_matrix.insert(2, 5)   = -0.0030819743036696748;
   ref_matrix.insert(3, 6)   = -0.0021767019962289917;
   ref_matrix.insert(4, 7)   =  0.0030986050890011384;
   ref_matrix.insert(8, 8)   = -0.0046226259917450515;
   ref_matrix.insert(8, 18)  = -0.00217920594520562;
   ref_matrix.insert(8, 24)  = -0.001539104896783423;
   ref_matrix.insert(8, 29)  =  0.0021909651729722815;
   ref_matrix.insert(9, 12)  = -0.0030818626029138676;
   ref_matrix.insert(9, 28)  = -0.0028428770609625234;
   ref_matrix.insert(10, 13) = -0.0021766230189459593;
   ref_matrix.insert(10, 32) = -0.0028428770609625234;
   ref_matrix.insert(11, 14) =  0.0030984926623045142;
   ref_matrix.insert(11, 19) = -0.0028428770609625234;
   ref_matrix.insert(12, 25) = -0.0028428770609625234;
   ref_matrix.insert(13, 17) = -0.0028428770609625234;
   ref_matrix.insert(14, 23) = -0.0028428770609625234;
   ref_matrix.insert(15, 30) = -0.0030819743036696761;
   ref_matrix.insert(16, 31) =  0.0030986050740594282;
   ref_matrix.insert(17, 32) = -0.00092412316568786881;
   ref_matrix.insert(18, 18) = -0.0061639486073393496;
   ref_matrix.insert(18, 24) =  0.0030986050740594282;
   ref_matrix.insert(18, 29) = -0.00092412316568786881;
   ref_matrix.insert(19, 23) =  0.0030986050740594282;
   ref_matrix.insert(20, 27) = -0.00092412316568786881;
   ref_matrix.insert(21, 33) = -0.0043768955684380072;
   ref_matrix.insert(22, 34) = -0.0040419594819556121;
   ref_matrix.insert(24, 24) = -0.0087537911368760126;
   ref_matrix.insert(24, 29) = -0.0040419594819556121;
   ref_matrix.insert(25, 28) = -0.0040419594819556121;
   ref_matrix.insert(26, 35) = -0.0031153256164892489;
   ref_matrix.insert(29, 29) = -0.006230651232978497;

   auto unitarityMatrix = CMSSM_unitarity::max_scattering_eigenvalue_infinite_s(m);
   for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
         BOOST_CHECK_CLOSE_FRACTION(unitarityMatrix.second.coeff(i,j), ref_matrix.coeff(i, j), 1e-16);
      }
   }
}
