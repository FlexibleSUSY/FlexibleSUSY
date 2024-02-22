(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

Needs["TestSuite`", "TestSuite.m"];
Needs["ReadSLHA`", "ReadSLHA.m"];

slha = "
    1 1 1  10.1
    2 2 2  20.2
    3 3 3  30.3
Block MODSEL                 # Select model
#   12    1000                # parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate all observables
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
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166370000e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
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
Block MINPAR                 # Input parameters
    1   0.192                # Lambda(Qin)
Block EXTPAR                 # Input parameters
    0   1000                 # input scale Qin
    1   173.34               # scale QEWSB
Block EXTPAR
    0   1.1
    1   2.2
Block V
    1   10.1
    2   20.2
    3   a
Block  V2
    1   1.1
    2   2.2
Block
    1   3.1
    2   4.2
Block M
    1 1  10.1
    1 2  1
    2 2  20.2
    3 3  30.3
Block T3
    1 1 1  10.1
    2 2 2  20.2
    3 3 3  30.3
Block T4
    1 1 1 1  10.1
    2 2 2 2  20.2
    3 3 3 3  30.3
";

pars = {
    {Qin  , {0}            , {EXTPAR, 0}},
    {QEWSB, {0}            , {EXTPAR, 1}},
    {V    , {3}            , V},
    {V2   , {2}            , V2},
    {M    , {3, 3}         , M},
    {T3   , {3, 3, 3}      , T3},
    {T4   , {3, 3, 3, 3}   , T4},
    {Ts   , {0}            , {T4, 2, 2, 2, 2}}
};

values = ReadSLHAString[slha, pars];

TestEquality[Qin /. values, 1.1];
TestEquality[QEWSB /. values, 2.2];
TestEquality[V /. values, {10.1, 20.2, 0}];
TestEquality[V2 /. values, {1.1, 2.2}];
TestEquality[M /. values, {{10.1, 1, 0}, {0, 20.2, 0}, {0, 0, 30.3}}];
TestEquality[T3 /. values, {{{10.1, 0, 0}, {0, 0, 0}, {0, 0, 0}},
                            {{0, 0, 0}, {0, 20.2, 0}, {0, 0, 0}},
                            {{0, 0, 0}, {0, 0, 0}, {0, 0, 30.3}}}];

t4 = Table[0, {i,1,3}, {j,1,3}, {k,1,3}, {l,1,3}];
t4[[1,1,1,1]] = 10.1;
t4[[2,2,2,2]] = 20.2;
t4[[3,3,3,3]] = 30.3;

TestEquality[T4 /. values, t4];
TestEquality[Ts /. values, 20.2];

Print[""];

PrintTestSummary[];
