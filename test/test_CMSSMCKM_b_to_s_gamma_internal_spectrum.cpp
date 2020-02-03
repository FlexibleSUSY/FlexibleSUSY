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
#define BOOST_TEST_MODULE test_CMSSMCKM_b_to_s_gamma

#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include <complex>

#include "lowe.h"
#include "test_complex_equality.hpp"

#include "test_CMSSMCKM.hpp"

#include "wrappers.hpp"
#include "CMSSMCKM_b_to_s_gamma.hpp"
#include "CMSSMCKM_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_b_to_s_gamma, * boost::unit_test::tolerance(1e-13) )
{
   char const * const slha_input = R"(
Block SPINFO
     1   FlexibleSUSY
     2   2.3.0
     5   CMSSMCKM
     9   4.14.1
Block MODSEL                 # Select model
    6   1                    # flavour violation
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
   17   160                    # pole mass renormalization scale (0 = SUSY scale)
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
Block VCKMIN                 # CKM matrix input (Wolfenstein parameters)
    1   2.272000000e-01      # lambda(MZ) SM DR-bar
    2   8.180000000e-01      # A(MZ) SM DR-bar
    3   2.210000000e-01      # rhobar(MZ) SM DR-bar
    4   3.400000000e-01      # etabar(MZ) SM DR-bar
Block MINPAR
     3     1.00000000E+01   # TanBeta
     4                  1   # SignMu
Block EXTPAR
    21     1.56250000E+04   # mHd2IN
    22     1.56250000E+04   # mHu2IN
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
  2  3     1.00000000E+05   # md2Input(2,3)
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
  2  3     1.00000000E+05   # mq2Input(2,3)
  3  1     0.00000000E+00   # mq2Input(3,1)
  3  2     0.00000000E+00   # mq2Input(3,2)
  3  3     1.56250000E+04   # mq2Input(3,3)
Block MSU2IN
  1  1     1.56250000E+04   # mu2Input(1,1)
  1  2     0.00000000E+00   # mu2Input(1,2)
  1  3     0.00000000E+00   # mu2Input(1,3)
  2  1     0.00000000E+00   # mu2Input(2,1)
  2  2     1.56250000E+04   # mu2Input(2,2)
  2  3     1.00000000E+05   # mu2Input(2,3)
  3  1     0.00000000E+00   # mu2Input(3,1)
  3  2     0.00000000E+00   # mu2Input(3,2)
  3  3     1.56250000E+04   # mu2Input(3,3)
Block TdIN
  1  1     0.00000000E+00   # TYdInput(1,1)
  1  2     0.00000000E+00   # TYdInput(1,2)
  1  3     0.00000000E+00   # TYdInput(1,3)
  2  1     0.00000000E+00   # TYdInput(2,1)
  2  2     0.00000000E+00   # TYdInput(2,2)
  2  3     0.00000000E+00   # TYdInput(2,3)
  3  1     0.00000000E+00   # TYdInput(3,1)
  3  2     0.00000000E+00   # TYdInput(3,2)
  3  3     0.00000000E+00   # TYdInput(3,3)
Block TeIN
  1  1     0.00000000E+00   # TYeInput(1,1)
  1  2     0.00000000E+00   # TYeInput(1,2)
  1  3     0.00000000E+00   # TYeInput(1,3)
  2  1     0.00000000E+00   # TYeInput(2,1)
  2  2     0.00000000E+00   # TYeInput(2,2)
  2  3     0.00000000E+00   # TYeInput(2,3)
  3  1     0.00000000E+00   # TYeInput(3,1)
  3  2     0.00000000E+00   # TYeInput(3,2)
  3  3     0.00000000E+00   # TYeInput(3,3)
Block TuIN
  1  1     0.00000000E+00   # TYuInput(1,1)
  1  2     0.00000000E+00   # TYuInput(1,2)
  1  3     0.00000000E+00   # TYuInput(1,3)
  2  1     0.00000000E+00   # TYuInput(2,1)
  2  2     0.00000000E+00   # TYuInput(2,2)
  2  3     0.00000000E+00   # TYuInput(2,3)
  3  1     0.00000000E+00   # TYuInput(3,1)
  3  2     0.00000000E+00   # TYuInput(3,2)
  3  3     0.00000000E+00   # TYuInput(3,3)
Block gauge Q= 1.60000000E+02
     1     3.56809477E-01   # g1 * 0.7745966692414834
     2     6.39545857E-01   # g2
     3     1.10162444E+00   # g3
Block Yu Q= 1.60000000E+02
  1  1     7.74816144E-06   # Yu(1,1)
  1  2     0.00000000E+00   # Yu(1,2)
  1  3     0.00000000E+00   # Yu(1,3)
  2  1     0.00000000E+00   # Yu(2,1)
  2  2     3.54316285E-03   # Yu(2,2)
  2  3     0.00000000E+00   # Yu(2,3)
  3  1     0.00000000E+00   # Yu(3,1)
  3  2     0.00000000E+00   # Yu(3,2)
  3  3     8.88812731E-01   # Yu(3,3)
Block IMYu Q= 1.60000000E+02
  1  1     0.00000000E+00   # Im(Yu(1,1))
  1  2     0.00000000E+00   # Im(Yu(1,2))
  1  3     0.00000000E+00   # Im(Yu(1,3))
  2  1     0.00000000E+00   # Im(Yu(2,1))
  2  2     0.00000000E+00   # Im(Yu(2,2))
  2  3     0.00000000E+00   # Im(Yu(2,3))
  3  1     0.00000000E+00   # Im(Yu(3,1))
  3  2     0.00000000E+00   # Im(Yu(3,2))
  3  3     0.00000000E+00   # Im(Yu(3,3))
Block Yd Q= 1.60000000E+02
  1  1     1.52806543E-04   # Yd(1,1)
  1  2     0.00000000E+00   # Yd(1,2)
  1  3     0.00000000E+00   # Yd(1,3)
  2  1     0.00000000E+00   # Yd(2,1)
  2  2     3.34567502E-03   # Yd(2,2)
  2  3     0.00000000E+00   # Yd(2,3)
  3  1     0.00000000E+00   # Yd(3,1)
  3  2     0.00000000E+00   # Yd(3,2)
  3  3     1.43080939E-01   # Yd(3,3)
Block IMYd Q= 1.60000000E+02
  1  1     0.00000000E+00   # Im(Yd(1,1))
  1  2     0.00000000E+00   # Im(Yd(1,2))
  1  3     0.00000000E+00   # Im(Yd(1,3))
  2  1     0.00000000E+00   # Im(Yd(2,1))
  2  2     0.00000000E+00   # Im(Yd(2,2))
  2  3     0.00000000E+00   # Im(Yd(2,3))
  3  1     0.00000000E+00   # Im(Yd(3,1))
  3  2     0.00000000E+00   # Im(Yd(3,2))
  3  3     0.00000000E+00   # Im(Yd(3,3))
Block Ye Q= 1.60000000E+02
  1  1     2.94339955E-05   # Ye(1,1)
  1  2     0.00000000E+00   # Ye(1,2)
  1  3     0.00000000E+00   # Ye(1,3)
  2  1     0.00000000E+00   # Ye(2,1)
  2  2     6.08600968E-03   # Ye(2,2)
  2  3     0.00000000E+00   # Ye(2,3)
  3  1     0.00000000E+00   # Ye(3,1)
  3  2     0.00000000E+00   # Ye(3,2)
  3  3     1.02325851E-01   # Ye(3,3)
Block IMYe Q= 1.60000000E+02
  1  1     0.00000000E+00   # Im(Ye(1,1))
  1  2     0.00000000E+00   # Im(Ye(1,2))
  1  3     0.00000000E+00   # Im(Ye(1,3))
  2  1     0.00000000E+00   # Im(Ye(2,1))
  2  2     0.00000000E+00   # Im(Ye(2,2))
  2  3     0.00000000E+00   # Im(Ye(2,3))
  3  1     0.00000000E+00   # Im(Ye(3,1))
  3  2     0.00000000E+00   # Im(Ye(3,2))
  3  3     0.00000000E+00   # Im(Ye(3,3))
Block Te Q= 1.60000000E+02
  1  1    -9.10779374E-03   # Re(TYe(1,1))
  1  2     0.00000000E+00   # Re(TYe(1,2))
  1  3     0.00000000E+00   # Re(TYe(1,3))
  2  1     0.00000000E+00   # Re(TYe(2,1))
  2  2    -1.88316099E+00   # Re(TYe(2,2))
  2  3     0.00000000E+00   # Re(TYe(2,3))
  3  1     0.00000000E+00   # Re(TYe(3,1))
  3  2     0.00000000E+00   # Re(TYe(3,2))
  3  3    -3.14741145E+01   # Re(TYe(3,3))
Block IMTe Q= 1.60000000E+02
  1  1     0.00000000E+00   # Im(TYe(1,1))
  1  2     0.00000000E+00   # Im(TYe(1,2))
  1  3     0.00000000E+00   # Im(TYe(1,3))
  2  1     0.00000000E+00   # Im(TYe(2,1))
  2  2     0.00000000E+00   # Im(TYe(2,2))
  2  3     0.00000000E+00   # Im(TYe(2,3))
  3  1     0.00000000E+00   # Im(TYe(3,1))
  3  2     0.00000000E+00   # Im(TYe(3,2))
  3  3     0.00000000E+00   # Im(TYe(3,3))
Block Td Q= 1.60000000E+02
  1  1    -2.38738505E-01   # Re(TYd(1,1))
  1  2    -4.73458875E-06   # Re(TYd(1,2))
  1  3     1.13037494E-04   # Re(TYd(1,3))
  2  1    -1.03196673E-04   # Re(TYd(2,1))
  2  2    -5.22657173E+00   # Re(TYd(2,2))
  2  3    -1.37731471E-02   # Re(TYd(2,3))
  3  1     1.05884237E-01   # Re(TYd(3,1))
  3  2    -5.89216998E-01   # Re(TYd(3,2))
  3  3    -2.08209330E+02   # Re(TYd(3,3))
Block IMTd Q= 1.60000000E+02
  1  1     8.55068136E-14   # Im(TYd(1,1))
  1  2    -2.16194965E-06   # Im(TYd(1,2))
  1  3     4.92508887E-05   # Im(TYd(1,3))
  2  1     4.66783136E-05   # Im(TYd(2,1))
  2  2    -3.69146241E-13   # Im(TYd(2,2))
  2  3     2.51614115E-04   # Im(TYd(2,3))
  3  1    -4.61388659E-02   # Im(TYd(3,1))
  3  2    -1.07642236E-02   # Im(TYd(3,2))
  3  3     1.74996183E-14   # Im(TYd(3,3))
Block Tu Q= 1.60000000E+02
  1  1    -9.77367229E-03   # Re(TYu(1,1))
  1  2     4.37387654E-09   # Re(TYu(1,2))
  1  3     4.40194499E-08   # Re(TYu(1,3))
  2  1     2.00013906E-06   # Re(TYu(2,1))
  2  2    -4.46937194E+00   # Re(TYu(2,2))
  2  3     3.91702400E-04   # Re(TYu(2,3))
  3  1     5.11544026E-03   # Re(TYu(3,1))
  3  2     9.95405044E-02   # Re(TYu(3,2))
  3  3    -8.50321580E+02   # Re(TYu(3,3))
Block IMTu Q= 1.60000000E+02
  1  1     9.30342798E-11   # Im(TYu(1,1))
  1  2    -2.65897782E-09   # Im(TYu(1,2))
  1  3    -6.79796980E-08   # Im(TYu(1,3))
  2  1     1.32699334E-06   # Im(TYu(2,1))
  2  2     1.49071554E-11   # Im(TYu(2,2))
  2  3     4.00034092E-09   # Im(TYu(2,3))
  3  1     7.89895344E-03   # Im(TYu(3,1))
  3  2    -9.54814497E-07   # Im(TYu(3,2))
  3  3     6.45528886E-14   # Im(TYu(3,3))
Block MSQ2 Q= 1.60000000E+02
  1  1     1.21997890E+06   # Re(mq2(1,1))
  1  2     5.64663948E+01   # Re(mq2(1,2))
  1  3    -1.35311883E+03   # Re(mq2(1,3))
  2  1     1.19021422E+02   # Re(mq2(2,1))
  2  2     1.21928784E+06   # Re(mq2(2,2))
  2  3    -8.13189866E+04   # Re(mq2(2,3))
  3  1    -1.36318982E+03   # Re(mq2(3,1))
  3  2     7.58523820E+03   # Re(mq2(3,2))
  3  3     1.03064794E+06   # Re(mq2(3,3))
Block IMMSQ2 Q= 1.60000000E+02
  1  1     9.57269332E-02   # Im(mq2(1,1))
  1  2     2.57392147E+01   # Im(mq2(1,2))
  1  3    -5.94890233E+02   # Im(mq2(1,3))
  2  1    -6.95642072E+01   # Im(mq2(2,1))
  2  2     7.17960065E+01   # Im(mq2(2,2))
  2  3     1.81697219E+04   # Im(mq2(2,3))
  3  1     5.93753361E+02   # Im(mq2(3,1))
  3  2     1.39771577E+02   # Im(mq2(3,2))
  3  3     3.89673323E+01   # Im(mq2(3,3))
Block MSE2 Q= 1.60000000E+02
  1  1     4.97104540E+04   # Re(me2(1,1))
  1  2     0.00000000E+00   # Re(me2(1,2))
  1  3     0.00000000E+00   # Re(me2(1,3))
  2  1     0.00000000E+00   # Re(me2(2,1))
  2  2     4.97046983E+04   # Re(me2(2,2))
  2  3     0.00000000E+00   # Re(me2(2,3))
  3  1     0.00000000E+00   # Re(me2(3,1))
  3  2     0.00000000E+00   # Re(me2(3,2))
  3  3     4.80857958E+04   # Re(me2(3,3))
Block IMMSE2 Q= 1.60000000E+02
  1  1     8.39650732E-01   # Im(me2(1,1))
  1  2     0.00000000E+00   # Im(me2(1,2))
  1  3     0.00000000E+00   # Im(me2(1,3))
  2  1     0.00000000E+00   # Im(me2(2,1))
  2  2     8.39646606E-01   # Im(me2(2,2))
  2  3     0.00000000E+00   # Im(me2(2,3))
  3  1     0.00000000E+00   # Im(me2(3,1))
  3  2     0.00000000E+00   # Im(me2(3,2))
  3  3     8.38484270E-01   # Im(me2(3,3))
Block MSL2 Q= 1.60000000E+02
  1  1     1.29097540E+05   # Re(ml2(1,1))
  1  2     0.00000000E+00   # Re(ml2(1,2))
  1  3     0.00000000E+00   # Re(ml2(1,3))
  2  1     0.00000000E+00   # Re(ml2(2,1))
  2  2     1.29094716E+05   # Re(ml2(2,2))
  2  3     0.00000000E+00   # Re(ml2(2,3))
  3  1     0.00000000E+00   # Re(ml2(3,1))
  3  2     0.00000000E+00   # Re(ml2(3,2))
  3  3     1.28300461E+05   # Re(ml2(3,3))
Block IMMSL2 Q= 1.60000000E+02
  1  1    -4.19825366E-01   # Im(ml2(1,1))
  1  2     0.00000000E+00   # Im(ml2(1,2))
  1  3     0.00000000E+00   # Im(ml2(1,3))
  2  1     0.00000000E+00   # Im(ml2(2,1))
  2  2    -4.19827419E-01   # Im(ml2(2,2))
  2  3     0.00000000E+00   # Im(ml2(2,3))
  3  1     0.00000000E+00   # Im(ml2(3,1))
  3  2     0.00000000E+00   # Im(ml2(3,2))
  3  3    -4.20405468E-01   # Im(ml2(3,3))
Block MSU2 Q= 1.60000000E+02
  1  1     1.14238297E+06   # Re(mu2(1,1))
  1  2     9.61225305E-06   # Re(mu2(1,2))
  1  3     5.64003237E-02   # Re(mu2(1,3))
  2  1     1.95511554E-06   # Re(mu2(2,1))
  2  2     1.14237598E+06   # Re(mu2(2,2))
  2  3     8.25273179E+04   # Re(mu2(2,3))
  3  1    -1.66829140E-05   # Re(mu2(3,1))
  3  2    -9.02262411E-02   # Re(mu2(3,2))
  3  3     7.66143489E+05   # Re(mu2(3,3))
Block IMMSU2 Q= 1.60000000E+02
  1  1    -5.59767153E-01   # Im(mu2(1,1))
  1  2    -1.98246332E-06   # Im(mu2(1,2))
  1  3    -1.16349431E-02   # Im(mu2(1,3))
  2  1     1.53181781E-06   # Im(mu2(2,1))
  2  2    -5.63643335E-01   # Im(mu2(2,2))
  2  3    -2.23579811E+01   # Im(mu2(2,3))
  3  1    -1.14165731E-05   # Im(mu2(3,1))
  3  2     3.74815439E-02   # Im(mu2(3,2))
  3  3     2.23197482E+02   # Im(mu2(3,3))
Block MSD2 Q= 1.60000000E+02
  1  1     1.13351654E+06   # Re(md2(1,1))
  1  2    -4.96207169E-06   # Re(md2(1,2))
  1  3    -3.18634946E-03   # Re(md2(1,3))
  2  1     4.11027218E-03   # Re(md2(2,1))
  2  2     1.13350866E+06   # Re(md2(2,2))
  2  3    -9.76065176E+04   # Re(md2(2,3))
  3  1     5.06325551E-03   # Re(md2(3,1))
  3  2    -6.16910440E-01   # Re(md2(3,2))
  3  3     1.12077464E+06   # Re(md2(3,3))
Block IMMSD2 Q= 1.60000000E+02
  1  1     2.79883576E-01   # Im(md2(1,1))
  1  2    -2.28048996E-06   # Im(md2(1,2))
  1  3     1.53179193E-04   # Im(md2(1,3))
  2  1    -2.89660811E-03   # Im(md2(2,1))
  2  2     3.83405948E-01   # Im(md2(2,2))
  2  3     2.01115310E+04   # Im(md2(2,3))
  3  1    -2.20517757E-03   # Im(md2(3,1))
  3  2    -1.13516475E-02   # Im(md2(3,2))
  3  3     1.89404158E-01   # Im(md2(3,3))
Block MSOFT Q= 1.60000000E+02
    21     1.07757752E+05   # mHd2
    22    -4.94282066E+05   # mHu2
     1     2.02252990E+02   # MassB
     2     3.82513966E+02   # MassWB
     3     1.19605665E+03   # MassG
Block Phases Q= 1.60000000E+02
     1     1.00000000E+00   # Re(PhaseGlu)
Block IMPhases Q= 1.60000000E+02
     1     0.00000000E+00   # Im(PhaseGlu)
Block MASS
   1000021     1.13943763E+03   # Glu
        24     8.04246449E+01   # VWm
   1000024     3.84347841E+02   # Cha(1)
   1000037     6.50743163E+02   # Cha(2)
        25     1.16106308E+02   # hh(1)
        35     7.17938170E+02   # hh(2)
        37     7.22631825E+02   # Hpm(2)
        36     7.17662200E+02   # Ah(2)
   1000012     3.51521068E+02   # Sv(1)
   1000014     3.52595542E+02   # Sv(2)
   1000016     3.52599356E+02   # Sv(3)
   1000022    -2.03668068E+02   # Chi(1)
   1000023    -3.84334442E+02   # Chi(2)
   1000025    -6.36960733E+02   # Chi(3)
   1000035    -6.50381584E+02   # Chi(4)
   1000001     9.41923225E+02   # Sd(1)
   1000003     9.55876265E+02   # Sd(2)
   1000005     1.00628253E+03   # Sd(3)
   2000001     1.05010719E+03   # Sd(4)
   2000003     1.05290483E+03   # Sd(5)
   2000005     1.06702422E+03   # Sd(6)
   1000011     2.21376343E+02   # Se(1)
   1000013     2.28253618E+02   # Se(2)
   1000015     2.28278063E+02   # Se(3)
   2000011     3.61475345E+02   # Se(4)
   2000013     3.61479933E+02   # Se(5)
   2000015     3.62704146E+02   # Se(6)
   1000002     7.92855780E+02   # Su(1)
   1000004     9.75957490E+02   # Su(2)
   1000006     1.00982711E+03   # Su(3)
   2000002     1.02333967E+03   # Su(4)
   2000004     1.04727616E+03   # Su(5)
   2000006     1.07217895E+03   # Su(6)
Block UMIX
  1  1     9.59728905E-01   # Re(UM(1,1))
  1  2    -2.80927800E-01   # Re(UM(1,2))
  2  1     2.80927800E-01   # Re(UM(2,1))
  2  2     9.59728905E-01   # Re(UM(2,2))
Block VMIX
  1  1     9.82333253E-01   # Re(UP(1,1))
  1  2    -1.87140003E-01   # Re(UP(1,2))
  2  1     1.87140003E-01   # Re(UP(2,1))
  2  2     9.82333253E-01   # Re(UP(2,2))
Block PSEUDOSCALARMIX
  1  1    -1.01086769E-01   # ZA(1,1)
  1  2     9.94877613E-01   # ZA(1,2)
  2  1     9.94877613E-01   # ZA(2,1)
  2  2     1.01086769E-01   # ZA(2,2)
Block DSQMIX
  1  1     4.68992183E-03   # Re(ZD(1,1))
  1  2     2.49744658E-01   # Re(ZD(1,2))
  1  3     6.95188146E-01   # Re(ZD(1,3))
  1  4     6.68515143E-07   # Re(ZD(1,4))
  1  5     2.97913312E-01   # Re(ZD(1,5))
  1  6     3.99626350E-01   # Re(ZD(1,6))
  2  1     3.62139889E-03   # Re(ZD(2,1))
  2  2     1.89657653E-01   # Re(ZD(2,2))
  2  3     4.71568070E-01   # Re(ZD(2,3))
  2  4     6.55432861E-07   # Re(ZD(2,4))
  2  5    -4.83321960E-01   # Re(ZD(2,5))
  2  6    -5.23517898E-01   # Re(ZD(2,6))
  3  1     2.15429181E-04   # Re(ZD(3,1))
  3  2    -1.42361290E-06   # Re(ZD(3,2))
  3  3    -1.49890103E-06   # Re(ZD(3,3))
  3  4     9.99999977E-01   # Re(ZD(3,4))
  3  5     2.37153981E-07   # Re(ZD(3,5))
  3  6     6.20175819E-09   # Re(ZD(3,6))
  4  1     9.99846824E-01   # Re(ZD(4,1))
  4  2    -1.43926670E-02   # Re(ZD(4,2))
  4  3     7.91713201E-04   # Re(ZD(4,3))
  4  4    -2.15424920E-04   # Re(ZD(4,4))
  4  5    -7.22476837E-04   # Re(ZD(4,5))
  4  6     7.36572252E-04   # Re(ZD(4,6))
  5  1    -4.27252456E-03   # Re(ZD(5,1))
  5  2    -1.54711872E-01   # Re(ZD(5,2))
  5  3     2.15607145E-02   # Re(ZD(5,3))
  5  4     8.46880217E-07   # Re(ZD(5,4))
  5  5    -6.07196216E-01   # Re(ZD(5,5))
  5  6     6.33763023E-01   # Re(ZD(5,6))
  6  1     1.59033577E-02   # Re(ZD(6,1))
  6  2     7.46467745E-01   # Re(ZD(6,2))
  6  3    -3.56376760E-01   # Re(ZD(6,3))
  6  4    -2.46137229E-06   # Re(ZD(6,4))
  6  5    -9.55007867E-02   # Re(ZD(6,5))
  6  6     1.25316898E-01   # Re(ZD(6,6))
Block SELMIX
  1  1    -0.00000000E+00   # Re(ZE(1,1))
  1  2     3.45384142E-17   # Re(ZE(1,2))
  1  3     1.41529899E-01   # Re(ZE(1,3))
  1  4    -0.00000000E+00   # Re(ZE(1,4))
  1  5     1.06928764E-15   # Re(ZE(1,5))
  1  6     9.89933981E-01   # Re(ZE(1,6))
  2  1     0.00000000E+00   # Re(ZE(2,1))
  2  2     8.76118474E-03   # Re(ZE(2,2))
  2  3     1.74080561E-16   # Re(ZE(2,3))
  2  4     0.00000000E+00   # Re(ZE(2,4))
  2  5     9.99961620E-01   # Re(ZE(2,5))
  2  6    -1.14829523E-15   # Re(ZE(2,6))
  3  1     4.23783376E-05   # Re(ZE(3,1))
  3  2     0.00000000E+00   # Re(ZE(3,2))
  3  3     0.00000000E+00   # Re(ZE(3,3))
  3  4     9.99999999E-01   # Re(ZE(3,4))
  3  5     0.00000000E+00   # Re(ZE(3,5))
  3  6     0.00000000E+00   # Re(ZE(3,6))
  4  1     9.99999999E-01   # Re(ZE(4,1))
  4  2     0.00000000E+00   # Re(ZE(4,2))
  4  3     0.00000000E+00   # Re(ZE(4,3))
  4  4    -4.23783376E-05   # Re(ZE(4,4))
  4  5     0.00000000E+00   # Re(ZE(4,5))
  4  6     0.00000000E+00   # Re(ZE(4,6))
  5  1     0.00000000E+00   # Re(ZE(5,1))
  5  2    -9.99961620E-01   # Re(ZE(5,2))
  5  3    -1.55613581E-14   # Re(ZE(5,3))
  5  4     0.00000000E+00   # Re(ZE(5,4))
  5  5     8.76118474E-03   # Re(ZE(5,5))
  5  6     2.24984044E-15   # Re(ZE(5,6))
  6  1     0.00000000E+00   # Re(ZE(6,1))
  6  2     1.57254670E-14   # Re(ZE(6,2))
  6  3    -9.89933981E-01   # Re(ZE(6,3))
  6  4     0.00000000E+00   # Re(ZE(6,4))
  6  5     4.33561647E-16   # Re(ZE(6,5))
  6  6     1.41529899E-01   # Re(ZE(6,6))
Block SCALARMIX
  1  1     1.04010140E-01   # ZH(1,1)
  1  2     9.94576237E-01   # ZH(1,2)
  2  1     9.94576237E-01   # ZH(2,1)
  2  2    -1.04010140E-01   # ZH(2,2)
Block NMIX
  1  1    -5.43246223E-11   # Re(ZN(1,1))
  1  2    -1.21262761E-10   # Re(ZN(1,2))
  1  3    -2.88926809E-10   # Re(ZN(1,3))
  1  4    -3.20672554E-10   # Re(ZN(1,4))
  2  1    -1.77648990E-10   # Re(ZN(2,1))
  2  2     6.70322561E-11   # Re(ZN(2,2))
  2  3     1.39384259E-09   # Re(ZN(2,3))
  2  4     1.22796249E-09   # Re(ZN(2,4))
  3  1    -3.32005621E-02   # Re(ZN(3,1))
  3  2     4.84157265E-02   # Re(ZN(3,2))
  3  3     7.03431028E-01   # Re(ZN(3,3))
  3  4     7.08334969E-01   # Re(ZN(3,4))
  4  1     2.24445103E-09   # Re(ZN(4,1))
  4  2    -3.30634894E-09   # Re(ZN(4,2))
  4  3    -4.84065842E-08   # Re(ZN(4,3))
  4  4    -4.98677466E-08   # Re(ZN(4,4))
Block CHARGEMIX
  1  1     1.00358870E-01   # Re(ZP(1,1))
  1  2    -9.94951304E-01   # Re(ZP(1,2))
  2  1     9.94951304E-01   # Re(ZP(2,1))
  2  2     1.00358870E-01   # Re(ZP(2,2))
Block USQMIX
  1  1    -3.28054384E-03   # Re(ZU(1,1))
  1  2    -5.43572804E-02   # Re(ZU(1,2))
  1  3    -3.80791616E-01   # Re(ZU(1,3))
  1  4     2.54316304E-08   # Re(ZU(1,4))
  1  5     1.79627839E-01   # Re(ZU(1,5))
  1  6    -8.17717829E-01   # Re(ZU(1,6))
  2  1     3.10387849E-03   # Re(ZU(2,1))
  2  2     3.53798660E-01   # Re(ZU(2,2))
  2  3     6.51286643E-01   # Re(ZU(2,3))
  2  4     2.02252294E-06   # Re(ZU(2,4))
  2  5     3.22206588E-01   # Re(ZU(2,5))
  2  6    -2.64027639E-01   # Re(ZU(2,6))
  3  1    -2.03966973E-05   # Re(ZU(3,1))
  3  2    -1.84101034E-06   # Re(ZU(3,2))
  3  3     2.44655574E-06   # Re(ZU(3,3))
  3  4    -9.99999496E-01   # Re(ZU(3,4))
  3  5     4.54017496E-06   # Re(ZU(3,5))
  3  6     2.46239293E-08   # Re(ZU(3,6))
  4  1     1.30720020E-03   # Re(ZU(4,1))
  4  2    -4.30249732E-03   # Re(ZU(4,2))
  4  3    -2.87311832E-02   # Re(ZU(4,3))
  4  4    -5.92070565E-08   # Re(ZU(4,4))
  4  5     1.82026991E-01   # Re(ZU(4,5))
  4  6     5.87687128E-02   # Re(ZU(4,6))
  5  1    -9.99957829E-01   # Re(ZU(5,1))
  5  2     6.69450787E-03   # Re(ZU(5,2))
  5  3    -3.83644967E-04   # Re(ZU(5,3))
  5  4     2.03882619E-05   # Re(ZU(5,4))
  5  5     1.56520766E-03   # Re(ZU(5,5))
  5  6     3.36521173E-03   # Re(ZU(5,6))
  6  1     7.88893627E-03   # Re(ZU(6,1))
  6  2     6.87466608E-01   # Re(ZU(6,2))
  6  3    -4.58463524E-01   # Re(ZU(6,3))
  6  4    -1.95232044E-06   # Re(ZU(6,4))
  6  5     1.16160519E-01   # Re(ZU(6,5))
  6  6     1.80657800E-01   # Re(ZU(6,6))
Block SNUMIX
  1  1     0.00000000E+00   # Re(ZV(1,1))
  1  2     0.00000000E+00   # Re(ZV(1,2))
  1  3     1.00000000E+00   # Re(ZV(1,3))
  2  1     0.00000000E+00   # Re(ZV(2,1))
  2  2     1.00000000E+00   # Re(ZV(2,2))
  2  3     0.00000000E+00   # Re(ZV(2,3))
  3  1     1.00000000E+00   # Re(ZV(3,1))
  3  2     0.00000000E+00   # Re(ZV(3,2))
  3  3     0.00000000E+00   # Re(ZV(3,3))
Block VCKM Q= 1.60000000E+02
  1  1     9.73840542E-01   # Re(CKM)(1,1)
  1  2     2.27197593E-01   # Re(CKM)(1,2)
  1  3     2.16779611E-03   # Re(CKM)(1,3)
  2  1    -2.27086799E-01   # Re(CKM)(2,1)
  2  2     9.72963908E-01   # Re(CKM)(2,2)
  2  3     4.21046284E-02   # Re(CKM)(2,3)
  3  1     7.45698952E-03   # Re(CKM)(3,1)
  3  2    -4.14959293E-02   # Re(CKM)(3,2)
  3  3     9.99105274E-01   # Re(CKM)(3,3)
Block IMVCKM Q= 1.60000000E+02
  1  1     8.35089783E-18   # Im(CKM)(1,1)
  1  2    -8.74312702E-18   # Im(CKM)(1,2)
  1  3    -3.33952067E-03   # Im(CKM)(1,3)
  2  1    -1.36933185E-04   # Im(CKM)(2,1)
  2  2    -3.19465956E-05   # Im(CKM)(2,2)
  2  3     2.70863600E-19   # Im(CKM)(2,3)
  3  1    -3.24930233E-03   # Im(CKM)(3,1)
  3  2    -7.58064217E-04   # Im(CKM)(3,2)
  3  3    -3.73012844E-17   # Im(CKM)(3,3)
Block ALPHA
          -1.04198591E-01   # ArcSin(Pole(ZH(2,2)))
Block HMIX Q= 1.60000000E+02
     1     6.23246434E+02   # Mu
     2     9.91512175E+00   # vu/vd
     3     2.48511560E+02   # Sqrt(Sqr(vd) + Sqr(vu))
     4     7.98685571E+05   # Sqr(MAh(2))
   101     7.97411484E+04   # BMu
   102     2.49373840E+01   # vd
   103     2.47257199E+02   # vu
Block FlexibleSUSYLowEnergy Q= 1.60000000E+02
    27     4.27183066E-02   # Re(C7) for b -> s gamma
)";

   std::stringstream istr(slha_input);

   CMSSMCKM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   slha_io.fill(qedqcd);

   CMSSMCKM_slha<CMSSMCKM<Two_scale>> model;
   slha_io.fill(model);
   model.calculate_spectrum();

   const auto reference_value = CMSSMCKM_b_to_s_gamma::calculate_b_to_s_gamma(model, qedqcd);
   constexpr std::complex<double> C7NP  {-0.043845625628194593, -0.00080098722119433636};
   constexpr std::complex<double> C7pNP {-0.011077563909982211, -0.00020236881118087242};
   constexpr std::complex<double> C8NP  {-0.032281142988385345, -0.00058972320838812729};
   constexpr std::complex<double> C8pNP {-0.024669363784006193, -0.00045066856414695082};

   TEST_COMPLEX_EQUALITY(C7NP, reference_value[0]);
   TEST_COMPLEX_EQUALITY(C7pNP, reference_value[1]);
   TEST_COMPLEX_EQUALITY(C8NP, reference_value[2]);
   TEST_COMPLEX_EQUALITY(C8pNP, reference_value[3]);
}
