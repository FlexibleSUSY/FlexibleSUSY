// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCPV_edm

#include <boost/test/unit_test.hpp>

#include "CMSSMCPV_two_scale_spectrum_generator.hpp"
#include "CMSSMCPV_slha_io.hpp"
#include "CMSSMCPV_edm.hpp"
#include "cxx_qft/CMSSMCPV_qft.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_edm, * boost::unit_test::tolerance(1e-5) )
{
   char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.7.1
     5   CMSSMCPV
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
    6   2                    # beta-functions loop order
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
     1     5.00000000E+02   # m0
     2     5.00000000E+02   # m12
     3     1.00000000E+01   # TanBeta
     4     1.00000000E+00   # CosPhiMu
     5     0.00000000E+00   # Azero
   100     1.00000000E+00   # PhaseMu
Block EXTPAR
   100     1.00000000E-01   # etaInput
Block IMMINPAR
     2     0.00000000E+00   # Imm12
     4     0.00000000E+00   # SinPhiMu
     5     0.00000000E+00   # ImAzero
Block gauge Q= 9.34088556E+02
     1     3.62245784E-01   # g1 * 0.7745966692414834
     2     6.42568447E-01   # g2
     3     1.05968121E+00   # g3
Block Yu Q= 9.34088556E+02
  1  1     7.30612124E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.34102321E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.59720076E-01   # Yu(3,3)
Block IMYu Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(Yu(1,1))
  1  2     0.00000000E+00   # Im(Yu(1,2))
  1  3     0.00000000E+00   # Im(Yu(1,3))
  2  1     0.00000000E+00   # Im(Yu(2,1))
  2  2     0.00000000E+00   # Im(Yu(2,2))
  2  3     0.00000000E+00   # Im(Yu(2,3))
  3  1     0.00000000E+00   # Im(Yu(3,1))
  3  2     0.00000000E+00   # Im(Yu(3,2))
  3  3     0.00000000E+00   # Im(Yu(3,3))
Block Yd Q= 9.34088556E+02
  1  1     1.40592524E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     3.07823830E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     1.33251137E-01   # Yd(3,3)
Block IMYd Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(Yd(1,1))
  1  2     0.00000000E+00   # Im(Yd(1,2))
  1  3     0.00000000E+00   # Im(Yd(1,3))
  2  1     0.00000000E+00   # Im(Yd(2,1))
  2  2     0.00000000E+00   # Im(Yd(2,2))
  2  3     0.00000000E+00   # Im(Yd(2,3))
  3  1     0.00000000E+00   # Im(Yd(3,1))
  3  2     0.00000000E+00   # Im(Yd(3,2))
  3  3     0.00000000E+00   # Im(Yd(3,3))
Block Ye Q= 9.34088556E+02
  1  1     2.88537489E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     5.96604084E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     1.00345029E-01   # Ye(3,3)
Block IMYe Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(Ye(1,1))
  1  2     0.00000000E+00   # Im(Ye(1,2))
  1  3     0.00000000E+00   # Im(Ye(1,3))
  2  1     0.00000000E+00   # Im(Ye(2,1))
  2  2     0.00000000E+00   # Im(Ye(2,2))
  2  3     0.00000000E+00   # Im(Ye(2,3))
  3  1     0.00000000E+00   # Im(Ye(3,1))
  3  2     0.00000000E+00   # Im(Ye(3,2))
  3  3     0.00000000E+00   # Im(Ye(3,3))
Block Te Q= 9.34088556E+02
  1  1    -8.63764502E-03   # Re(TYe(1,1))
  1  2     0.00000000E+00   # Re(TYe(1,2))
  1  3     0.00000000E+00   # Re(TYe(1,3))
  2  1     0.00000000E+00   # Re(TYe(2,1))
  2  2    -1.78595687E+00   # Re(TYe(2,2))
  2  3     0.00000000E+00   # Re(TYe(2,3))
  3  1     0.00000000E+00   # Re(TYe(3,1))
  3  2     0.00000000E+00   # Re(TYe(3,2))
  3  3    -2.98756908E+01   # Re(TYe(3,3))
Block IMTe Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(TYe(1,1))
  1  2     0.00000000E+00   # Im(TYe(1,2))
  1  3     0.00000000E+00   # Im(TYe(1,3))
  2  1     0.00000000E+00   # Im(TYe(2,1))
  2  2     0.00000000E+00   # Im(TYe(2,2))
  2  3     0.00000000E+00   # Im(TYe(2,3))
  3  1     0.00000000E+00   # Im(TYe(3,1))
  3  2     0.00000000E+00   # Im(TYe(3,2))
  3  3     0.00000000E+00   # Im(TYe(3,3))
Block Td Q= 9.34088556E+02
  1  1    -1.95939081E-01   # Re(TYd(1,1))
  1  2     0.00000000E+00   # Re(TYd(1,2))
  1  3     0.00000000E+00   # Re(TYd(1,3))
  2  1     0.00000000E+00   # Re(TYd(2,1))
  2  2    -4.29002298E+00   # Re(TYd(2,2))
  2  3     0.00000000E+00   # Re(TYd(2,3))
  3  1     0.00000000E+00   # Re(TYd(3,1))
  3  2     0.00000000E+00   # Re(TYd(3,2))
  3  3    -1.73579027E+02   # Re(TYd(3,3))
Block IMTd Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(TYd(1,1))
  1  2     0.00000000E+00   # Im(TYd(1,2))
  1  3     0.00000000E+00   # Im(TYd(1,3))
  2  1     0.00000000E+00   # Im(TYd(2,1))
  2  2     0.00000000E+00   # Im(TYd(2,2))
  2  3     0.00000000E+00   # Im(TYd(2,3))
  3  1     0.00000000E+00   # Im(TYd(3,1))
  3  2     0.00000000E+00   # Im(TYd(3,2))
  3  3     0.00000000E+00   # Im(TYd(3,3))
Block Tu Q= 9.34088556E+02
  1  1    -8.32505011E-03   # Re(TYu(1,1))
  1  2     0.00000000E+00   # Re(TYu(1,2))
  1  3     0.00000000E+00   # Re(TYu(1,3))
  2  1     0.00000000E+00   # Re(TYu(2,1))
  2  2    -3.80695332E+00   # Re(TYu(2,2))
  2  3     0.00000000E+00   # Re(TYu(2,3))
  3  1     0.00000000E+00   # Re(TYu(3,1))
  3  2     0.00000000E+00   # Re(TYu(3,2))
  3  3    -7.56353190E+02   # Re(TYu(3,3))
Block IMTu Q= 9.34088556E+02
  1  1     1.08420217E-19   # Im(TYu(1,1))
  1  2     0.00000000E+00   # Im(TYu(1,2))
  1  3     0.00000000E+00   # Im(TYu(1,3))
  2  1     0.00000000E+00   # Im(TYu(2,1))
  2  2    -5.55111512E-17   # Im(TYu(2,2))
  2  3     0.00000000E+00   # Im(TYu(2,3))
  3  1     0.00000000E+00   # Im(TYu(3,1))
  3  2     0.00000000E+00   # Im(TYu(3,2))
  3  3     4.26325641E-14   # Im(TYu(3,3))
Block MSQ2 Q= 9.34088556E+02
  1  1     1.25740796E+06   # Re(mq2(1,1))
  1  2     0.00000000E+00   # Re(mq2(1,2))
  1  3     0.00000000E+00   # Re(mq2(1,3))
  2  1     0.00000000E+00   # Re(mq2(2,1))
  2  2     1.25740060E+06   # Re(mq2(2,2))
  2  3     0.00000000E+00   # Re(mq2(2,3))
  3  1     0.00000000E+00   # Re(mq2(3,1))
  3  2     0.00000000E+00   # Re(mq2(3,2))
  3  3     1.02671403E+06   # Re(mq2(3,3))
Block IMMSQ2 Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(mq2(1,1))
  1  2     0.00000000E+00   # Im(mq2(1,2))
  1  3     0.00000000E+00   # Im(mq2(1,3))
  2  1     0.00000000E+00   # Im(mq2(2,1))
  2  2     0.00000000E+00   # Im(mq2(2,2))
  2  3     0.00000000E+00   # Im(mq2(2,3))
  3  1     0.00000000E+00   # Im(mq2(3,1))
  3  2     0.00000000E+00   # Im(mq2(3,2))
  3  3     0.00000000E+00   # Im(mq2(3,3))
Block MSE2 Q= 9.34088556E+02
  1  1     2.82655212E+05   # Re(me2(1,1))
  1  2     0.00000000E+00   # Re(me2(1,2))
  1  3     0.00000000E+00   # Re(me2(1,3))
  2  1     0.00000000E+00   # Re(me2(2,1))
  2  2     2.82636031E+05   # Re(me2(2,2))
  2  3     0.00000000E+00   # Re(me2(2,3))
  3  1     0.00000000E+00   # Re(me2(3,1))
  3  2     0.00000000E+00   # Re(me2(3,2))
  3  3     2.77231254E+05   # Re(me2(3,3))
Block IMMSE2 Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(me2(1,1))
  1  2     0.00000000E+00   # Im(me2(1,2))
  1  3     0.00000000E+00   # Im(me2(1,3))
  2  1     0.00000000E+00   # Im(me2(2,1))
  2  2     0.00000000E+00   # Im(me2(2,2))
  2  3     0.00000000E+00   # Im(me2(2,3))
  3  1     0.00000000E+00   # Im(me2(3,1))
  3  2     0.00000000E+00   # Im(me2(3,2))
  3  3     0.00000000E+00   # Im(me2(3,3))
Block MSL2 Q= 9.34088556E+02
  1  1     3.57226602E+05   # Re(ml2(1,1))
  1  2     0.00000000E+00   # Re(ml2(1,2))
  1  3     0.00000000E+00   # Re(ml2(1,3))
  2  1     0.00000000E+00   # Re(ml2(2,1))
  2  2     3.57217095E+05   # Re(ml2(2,2))
  2  3     0.00000000E+00   # Re(ml2(2,3))
  3  1     0.00000000E+00   # Re(ml2(3,1))
  3  2     0.00000000E+00   # Re(ml2(3,2))
  3  3     3.54538716E+05   # Re(ml2(3,3))
Block IMMSL2 Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(ml2(1,1))
  1  2     0.00000000E+00   # Im(ml2(1,2))
  1  3     0.00000000E+00   # Im(ml2(1,3))
  2  1     0.00000000E+00   # Im(ml2(2,1))
  2  2     0.00000000E+00   # Im(ml2(2,2))
  2  3     0.00000000E+00   # Im(ml2(2,3))
  3  1     0.00000000E+00   # Im(ml2(3,1))
  3  2     0.00000000E+00   # Im(ml2(3,2))
  3  3     0.00000000E+00   # Im(ml2(3,3))
Block MSU2 Q= 9.34088556E+02
  1  1     1.18342122E+06   # Re(mu2(1,1))
  1  2     0.00000000E+00   # Re(mu2(1,2))
  1  3     0.00000000E+00   # Re(mu2(1,3))
  2  1     0.00000000E+00   # Re(mu2(2,1))
  2  2     1.18341357E+06   # Re(mu2(2,2))
  2  3     0.00000000E+00   # Re(mu2(2,3))
  3  1     0.00000000E+00   # Re(mu2(3,1))
  3  2     0.00000000E+00   # Re(mu2(3,2))
  3  3     7.26030648E+05   # Re(mu2(3,3))
Block IMMSU2 Q= 9.34088556E+02
  1  1    -1.45519152E-11   # Im(mu2(1,1))
  1  2     0.00000000E+00   # Im(mu2(1,2))
  1  3     0.00000000E+00   # Im(mu2(1,3))
  2  1     0.00000000E+00   # Im(mu2(2,1))
  2  2     1.45519152E-11   # Im(mu2(2,2))
  2  3     0.00000000E+00   # Im(mu2(2,3))
  3  1     0.00000000E+00   # Im(mu2(3,1))
  3  2     0.00000000E+00   # Im(mu2(3,2))
  3  3     0.00000000E+00   # Im(mu2(3,3))
Block MSD2 Q= 9.34088556E+02
  1  1     1.17470089E+06   # Re(md2(1,1))
  1  2     0.00000000E+00   # Re(md2(1,2))
  1  3     0.00000000E+00   # Re(md2(1,3))
  2  1     0.00000000E+00   # Re(md2(2,1))
  2  2     1.17469365E+06   # Re(md2(2,2))
  2  3     0.00000000E+00   # Re(md2(2,3))
  3  1     0.00000000E+00   # Re(md2(3,1))
  3  2     0.00000000E+00   # Re(md2(3,2))
  3  3     1.16158079E+06   # Re(md2(3,3))
Block IMMSD2 Q= 9.34088556E+02
  1  1     0.00000000E+00   # Im(md2(1,1))
  1  2     0.00000000E+00   # Im(md2(1,2))
  1  3     0.00000000E+00   # Im(md2(1,3))
  2  1     0.00000000E+00   # Im(md2(2,1))
  2  2     0.00000000E+00   # Im(md2(2,2))
  2  3     0.00000000E+00   # Im(md2(2,3))
  3  1     0.00000000E+00   # Im(md2(3,1))
  3  2     0.00000000E+00   # Im(md2(3,2))
  3  3     0.00000000E+00   # Im(md2(3,3))
Block Phases Q= 9.34088556E+02
     1     1.00000000E+00   # Re(PhaseGlu)
Block IMPhases Q= 9.34088556E+02
     1     0.00000000E+00   # Im(PhaseGlu)
Block MASS
   1000021     1.16773028E+03   # Glu
        24     8.03783208E+01   # VWm
   1000024     3.89029195E+02   # Cha(1)
   1000037     6.49497652E+02   # Cha(2)
        37     8.64998421E+02   # Hpm(2)
   1000012     5.94173855E+02   # Sv(1)
   1000014     5.96523130E+02   # Sv(2)
   1000016     5.96531451E+02   # Sv(3)
   1000022    -2.05656303E+02   # Chi(1)
   1000023    -3.89020005E+02   # Chi(2)
   1000025    -6.35232661E+02   # Chi(3)
   1000035    -6.49244258E+02   # Chi(4)
        25     1.10668434E+02   # hh(2)
        35     8.60958860E+02   # hh(3)
        36     8.61163364E+02   # hh(4)
   1000001     1.04318879E+03   # Sd(1)
   1000003     1.11399514E+03   # Sd(2)
   1000005     1.11940984E+03   # Sd(3)
   2000001     1.11941416E+03   # Sd(4)
   2000003     1.15854461E+03   # Sd(5)
   2000005     1.15854719E+03   # Sd(6)
   1000011     5.28103166E+02   # Se(1)
   1000013     5.34933792E+02   # Se(2)
   1000015     5.34958213E+02   # Se(3)
   2000011     6.01060815E+02   # Se(4)
   2000013     6.02008352E+02   # Se(5)
   2000015     6.02011419E+02   # Se(6)
   1000002     8.59787187E+02   # Su(1)
   1000004     1.07955473E+03   # Su(2)
   1000006     1.12233122E+03   # Su(3)
   2000002     1.12233786E+03   # Su(4)
   2000004     1.15600996E+03   # Su(5)
   2000006     1.15601014E+03   # Su(6)
        21     0.00000000E+00   # VG
        12     0.00000000E+00   # Fv(1)
        14     0.00000000E+00   # Fv(2)
        16     0.00000000E+00   # Fv(3)
        11     5.29985893E-04   # Fe(1)
        13     1.07449437E-01   # Fe(2)
        15     1.78802360E+00   # Fe(3)
         1     4.57455916E-03   # Fd(1)
         3     9.04955381E-02   # Fd(2)
         5     3.39061995E+00   # Fd(3)
         2     2.31276738E-03   # Fu(1)
         4     8.52001277E-01   # Fu(2)
         6     1.73801276E+02   # Fu(3)
        22     0.00000000E+00   # VP
        23     9.11240165E+01   # VZ
Block UMIX
  1  1     9.54717937E-01   # Re(UM(1,1))
  1  2    -2.86777509E-01   # Re(UM(1,2))
  2  1     2.85627366E-01   # Re(UM(2,1))
  2  2     9.57995189E-01   # Re(UM(2,2))
Block VMIX
  1  1     9.77756882E-01   # Re(UP(1,1))
  1  2    -1.93244089E-01   # Re(UP(1,2))
  2  1     1.92575652E-01   # Re(UP(2,1))
  2  2     9.81150713E-01   # Re(UP(2,2))
Block DSQMIX
  1  1    -0.00000000E+00   # Re(ZD(1,1))
  1  2    -1.30944174E-15   # Re(ZD(1,2))
  1  3     9.92782094E-01   # Re(ZD(1,3))
  1  4    -1.58606855E-31   # Re(ZD(1,4))
  1  5    -2.56258734E-15   # Re(ZD(1,5))
  1  6     1.12402263E-01   # Re(ZD(1,6))
  2  1     0.00000000E+00   # Re(ZD(2,1))
  2  2     4.05874839E-16   # Re(ZD(2,2))
  2  3    -1.13170939E-01   # Re(ZD(2,3))
  2  4     4.68074937E-32   # Re(ZD(2,4))
  2  5     4.54701848E-15   # Re(ZD(2,5))
  2  6     9.86857095E-01   # Re(ZD(2,6))
  3  1     0.00000000E+00   # Re(ZD(3,1))
  3  2    -4.76688854E-03   # Re(ZD(3,2))
  3  3    -3.60930281E-15   # Re(ZD(3,3))
  3  4    -5.30219302E-19   # Re(ZD(3,4))
  3  5    -9.96733657E-01   # Re(ZD(3,5))
  3  6     4.31826953E-15   # Re(ZD(3,6))
  4  1    -2.17725619E-04   # Re(ZD(4,1))
  4  2     0.00000000E+00   # Re(ZD(4,2))
  4  3     0.00000000E+00   # Re(ZD(4,3))
  4  4    -9.96687456E-01   # Re(ZD(4,4))
  4  5     9.87959964E-22   # Re(ZD(4,5))
  4  6     0.00000000E+00   # Re(ZD(4,6))
  5  1     0.00000000E+00   # Re(ZD(5,1))
  5  2     9.99988386E-01   # Re(ZD(5,2))
  5  3     1.44126259E-15   # Re(ZD(5,3))
  5  4     1.11021036E-16   # Re(ZD(5,4))
  5  5    -4.75137343E-03   # Re(ZD(5,5))
  5  6    -2.34864571E-16   # Re(ZD(5,6))
  6  1     9.99999976E-01   # Re(ZD(6,1))
  6  2     0.00000000E+00   # Re(ZD(6,2))
  6  3     0.00000000E+00   # Re(ZD(6,3))
  6  4    -2.17004398E-04   # Re(ZD(6,4))
  6  5     2.15104200E-25   # Re(ZD(6,5))
  6  6     0.00000000E+00   # Re(ZD(6,6))
Block SELMIX
  1  1     0.00000000E+00   # Re(ZE(1,1))
  1  2    -7.67736965E-17   # Re(ZE(1,2))
  1  3     1.42015083E-01   # Re(ZE(1,3))
  1  4     1.71660484E-32   # Re(ZE(1,4))
  1  5     5.37643874E-16   # Re(ZE(1,5))
  1  6     9.82576706E-01   # Re(ZE(1,6))
  2  1     0.00000000E+00   # Re(ZE(2,1))
  2  2    -9.03456807E-03   # Re(ZE(2,2))
  2  3    -4.17217213E-16   # Re(ZE(2,3))
  2  4     2.00678726E-18   # Re(ZE(2,4))
  2  5    -9.92152962E-01   # Re(ZE(2,5))
  2  6     6.36461714E-16   # Re(ZE(2,6))
  3  1     4.37248200E-05   # Re(ZE(3,1))
  3  2     0.00000000E+00   # Re(ZE(3,2))
  3  3     1.35781131E-35   # Re(ZE(3,3))
  3  4     9.95473412E-01   # Re(ZE(3,4))
  3  5    -1.83184323E-22   # Re(ZE(3,5))
  3  6     2.38441537E-37   # Re(ZE(3,6))
  4  1    -0.00000000E+00   # Re(ZE(4,1))
  4  2    -5.69053404E-14   # Re(ZE(4,2))
  4  3     9.89617908E-01   # Re(ZE(4,3))
  4  4     1.26417254E-29   # Re(ZE(4,4))
  4  5     7.21571754E-17   # Re(ZE(4,5))
  4  6    -1.41088927E-01   # Re(ZE(4,6))
  5  1     0.00000000E+00   # Re(ZE(5,1))
  5  2    -9.99513712E-01   # Re(ZE(5,2))
  5  3    -5.63560741E-14   # Re(ZE(5,3))
  5  4     2.22035536E-16   # Re(ZE(5,4))
  5  5     8.96803452E-03   # Re(ZE(5,5))
  5  6     7.92397293E-15   # Re(ZE(5,6))
  6  1     9.99999999E-01   # Re(ZE(6,1))
  6  2     0.00000000E+00   # Re(ZE(6,2))
  6  3    -5.93700552E-40   # Re(ZE(6,3))
  6  4    -4.35268958E-05   # Re(ZE(6,4))
  6  5     8.00970157E-27   # Re(ZE(6,5))
  6  6    -1.04258133E-41   # Re(ZE(6,6))
Block SCALARMIX
  1  1     7.73380792E-07   # ZH(1,1)
  1  2    -2.28377554E-06   # ZH(1,2)
  1  3     1.02827415E-01   # ZH(1,3)
  1  4    -9.94699212E-01   # ZH(1,4)
  2  1     1.05482902E-01   # ZH(2,1)
  2  2     9.94421116E-01   # ZH(2,2)
  2  3     2.93910980E-05   # ZH(2,3)
  2  4     8.37192245E-07   # ZH(2,4)
  3  1    -1.07881558E-02   # ZH(3,1)
  3  2     1.11486602E-03   # ZH(3,2)
  3  3     9.94640710E-01   # ZH(3,3)
  3  4     1.02821356E-01   # ZH(3,4)
  4  1     9.94362596E-01   # ZH(4,1)
  4  2    -1.05477014E-01   # ZH(4,2)
  4  3     1.07879753E-02   # ZH(4,3)
  4  4     1.11622640E-03   # ZH(4,4)
Block NMIX
  1  1     1.36642976E-04   # Re(ZN(1,1))
  1  2    -1.78346681E-04   # Re(ZN(1,2))
  1  3     7.74387293E-03   # Re(ZN(1,3))
  1  4     2.50283059E-03   # Re(ZN(1,4))
  2  1    -2.70818204E-04   # Re(ZN(2,1))
  2  2     2.48603194E-04   # Re(ZN(2,2))
  2  3    -1.79113487E-02   # Re(ZN(2,3))
  2  4    -1.10669609E-02   # Re(ZN(2,4))
  3  1     3.32790406E-02   # Re(ZN(3,1))
  3  2    -4.83666676E-02   # Re(ZN(3,2))
  3  3    -7.00646295E-01   # Re(ZN(3,3))
  3  4    -7.05449929E-01   # Re(ZN(3,4))
  4  1     5.66015588E-04   # Re(ZN(4,1))
  4  2    -6.67361123E-04   # Re(ZN(4,2))
  4  3     6.06131280E-02   # Re(ZN(4,3))
  4  4     6.22970451E-02   # Re(ZN(4,4))
Block CHARGEMIX
  1  1     9.92305414E-02   # Re(ZP(1,1))
  1  2    -9.95064389E-01   # Re(ZP(1,2))
  2  1     9.95064470E-01   # Re(ZP(2,1))
  2  2     9.92305333E-02   # Re(ZP(2,2))
Block USQMIX
  1  1     0.00000000E+00   # Re(ZU(1,1))
  1  2     1.09768609E-16   # Re(ZU(1,2))
  1  3    -3.67048906E-01   # Re(ZU(1,3))
  1  4    -7.10146060E-33   # Re(ZU(1,4))
  1  5    -1.58553955E-16   # Re(ZU(1,5))
  1  6    -9.30106591E-01   # Re(ZU(1,6))
  2  1     0.00000000E+00   # Re(ZU(2,1))
  2  2    -1.31561962E-15   # Re(ZU(2,2))
  2  3     9.30172245E-01   # Re(ZU(2,3))
  2  4     1.57330598E-31   # Re(ZU(2,4))
  2  5     1.33179361E-15   # Re(ZU(2,5))
  2  6    -3.67022374E-01   # Re(ZU(2,6))
  3  1     0.00000000E+00   # Re(ZU(3,1))
  3  2     9.77143400E-03   # Re(ZU(3,2))
  3  3    -1.39710574E-15   # Re(ZU(3,3))
  3  4    -1.08506394E-18   # Re(ZU(3,4))
  3  5     9.99691881E-01   # Re(ZU(3,5))
  3  6     4.19872647E-16   # Re(ZU(3,6))
  4  1    -2.13746161E-05   # Re(ZU(4,1))
  4  2     0.00000000E+00   # Re(ZU(4,2))
  4  3     0.00000000E+00   # Re(ZU(4,3))
  4  4    -9.99985979E-01   # Re(ZU(4,4))
  4  5    -5.02023179E-23   # Re(ZU(4,5))
  4  6     0.00000000E+00   # Re(ZU(4,6))
  5  1     0.00000000E+00   # Re(ZU(5,1))
  5  2    -9.99798692E-01   # Re(ZU(5,2))
  5  3    -1.43500367E-15   # Re(ZU(5,3))
  5  4     1.11017000E-16   # Re(ZU(5,4))
  5  5     9.77039009E-03   # Re(ZU(5,5))
  5  6     3.87572313E-16   # Re(ZU(5,6))
  6  1     1.00000000E+00   # Re(ZU(6,1))
  6  2     0.00000000E+00   # Re(ZU(6,2))
  6  3     0.00000000E+00   # Re(ZU(6,3))
  6  4    -2.13743164E-05   # Re(ZU(6,4))
  6  5    -1.07305527E-27   # Re(ZU(6,5))
  6  6     0.00000000E+00   # Re(ZU(6,6))
Block SNUMIX
  1  1     0.00000000E+00   # Re(ZV(1,1))
  1  2    -7.67408795E-17   # Re(ZV(1,2))
  1  3    -9.86916981E-01   # Re(ZV(1,3))
  2  1     0.00000000E+00   # Re(ZV(2,1))
  2  2    -9.99997391E-01   # Re(ZV(2,2))
  2  3     7.89370401E-17   # Re(ZV(2,3))
  3  1     1.00000000E+00   # Re(ZV(3,1))
  3  2     0.00000000E+00   # Re(ZV(3,2))
  3  3     0.00000000E+00   # Re(ZV(3,3))
Block UELMIX
  1  1     9.99999007E-01   # Re(ZEL(1,1))
  1  2     0.00000000E+00   # Re(ZEL(1,2))
  1  3     0.00000000E+00   # Re(ZEL(1,3))
  2  1     0.00000000E+00   # Re(ZEL(2,1))
  2  2     9.99998964E-01   # Re(ZEL(2,2))
  2  3     0.00000000E+00   # Re(ZEL(2,3))
  3  1     0.00000000E+00   # Re(ZEL(3,1))
  3  2     0.00000000E+00   # Re(ZEL(3,2))
  3  3     9.99998937E-01   # Re(ZEL(3,3))
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
  1  1     9.99997310E-01   # Re(ZDL(1,1))
  1  2     0.00000000E+00   # Re(ZDL(1,2))
  1  3     0.00000000E+00   # Re(ZDL(1,3))
  2  1    -0.00000000E+00   # Re(ZDL(2,1))
  2  2    -9.99989450E-01   # Re(ZDL(2,2))
  2  3    -4.31329332E-17   # Re(ZDL(2,3))
  3  1     0.00000000E+00   # Re(ZDL(3,1))
  3  2    -4.31325543E-17   # Re(ZDL(3,2))
  3  3     9.99998021E-01   # Re(ZDL(3,3))
Block UDRMIX
  1  1     1.00000000E+00   # Re(ZDR(1,1))
  1  2     0.00000000E+00   # Re(ZDR(1,2))
  1  3     0.00000000E+00   # Re(ZDR(1,3))
  2  1     0.00000000E+00   # Re(ZDR(2,1))
  2  2    -9.99998022E-01   # Re(ZDR(2,2))
  2  3    -1.61184970E-15   # Re(ZDR(2,3))
  3  1     0.00000000E+00   # Re(ZDR(3,1))
  3  2    -1.61184651E-15   # Re(ZDR(3,2))
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
)";

   std::stringstream istr(slha_input);

   CMSSMCPV_slha_io slha_io;
   slha_io.read_from_stream(istr);

   // extract the input parameters
   softsusy::QedQcd qedqcd;
   CMSSMCPV_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   settings.set(Spectrum_generator_settings::calculate_sm_masses, 0);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 0);
   CMSSMCPV_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   auto m = std::get<0>(spectrum_generator.get_models_slha());

   using CMSSMCPV_cxx_diagrams::fields::Fe;

   const double eps = 1e-11;

   auto de = CMSSMCPV_edm::calculate_edm<Fe>(m, qedqcd, 0);
   BOOST_CHECK_CLOSE_FRACTION(de, -3.2103089248010137e-13, 9e-13);

   auto dmu = CMSSMCPV_edm::calculate_edm<Fe>(m, qedqcd, 1);
   BOOST_CHECK_CLOSE_FRACTION(dmu, -6.6379955917989794e-11, 6e-13);

   auto dtau = CMSSMCPV_edm::calculate_edm<Fe>(m, qedqcd, 2);
   BOOST_CHECK_CLOSE_FRACTION(dtau, -1.1209214170045263e-09, 2e-12);
}
