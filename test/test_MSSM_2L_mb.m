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

b2l = << "meta/MSSM/bquark_2loop_sqcd_decoupling.m";
<< "meta/MSSM/twoloopbubble.m";

colorCA = 3; colorCF = 4/3; GS = g3;
zt2 = Zeta[2]; zt3 = Zeta[3];
MGl = mgl; MT = mt; MB = mb;
mmgl = mgl^2; mmt = mt^2;
mmst1 = mst1^2; mmst2 = mst2^2;
mmsb1 = msb1^2; mmsb2 = msb2^2;
mmsusy = msusy^2; mmu = scale^2;
s2t = 2 mt Xt/(mmst1 - mmst2);
s2b = 2 mb Xb/(mmsb1 - mmsb2);

Delta[mass1_, mass2_, mass3_, x_] := (mass1^2 + mass2^2 + mass3^2 - 2*(mass1*mass2 + mass1*mass3 + mass2*mass3))^x;
fin[0, mass1_, mass2_]            := Fin20[mass1, mass2, mmu];
fin[mass1_, mass2_, mass3_]       := Fin3[mass1, mass2, mass3, mmu];

MakePoint[MG_, Mst1_, Mst2_, Msb1_, Msb2_, MS_, xt_, xb_, Q_] :=
    {
        g3 -> Sqrt[4 Pi 0.1184],
        msusy -> MS,
        mb -> 3.2,
        mt -> 173.34,
        Xt -> xt,
        Xb -> xb,
        scale -> Q,
        mgl -> MG,
        mst1 -> Mst1, mst2 -> Mst2,
        msb1 -> Msb1, msb2 -> Msb2
    }

(* 2-loop SQCD correction to mb in the MSSM at non-universal SUSY
   masses; the reference values are those of the C++ implementation
   (test/test_mssm_twoloop_mb.cpp) *)
points = {
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,               0, 0, 1000],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,               0, 0, 1100],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,  Sqrt[6] * 1300, 0, 1000],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,  Sqrt[6] * 1300, 0, 1100],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400, -Sqrt[6] * 1300, 0, 1000],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400, -Sqrt[6] * 1300, 0, 1100],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,               0,  Sqrt[6], 1000],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,               0, -Sqrt[6], 1100],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,  Sqrt[6] * 1300,  Sqrt[6], 1000],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400,  Sqrt[6] * 1300, -Sqrt[6], 1100],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400, -Sqrt[6] * 1300,  Sqrt[6], 1000],
    MakePoint[1000, 1200, 1300, 1500, 1600, 1400, -Sqrt[6] * 1300, -Sqrt[6], 1100],
    MakePoint[1600, 1500, 1400, 1300, 1200, 1100,               0, -Sqrt[6], 1000],
    MakePoint[1600, 1500, 1400, 1300, 1200, 1100,  Sqrt[6] * 1300, -Sqrt[6], 1000],
    MakePoint[1600, 1500, 1400, 1300, 1200, 1100, -Sqrt[6] * 1300, -Sqrt[6], 1000],
    MakePoint[1600, 1500, 1400, 1300, 1200, 1100,               0,        0, 1000],
    MakePoint[1600, 1500, 1400, 1300, 1200, 1100,  Sqrt[6] * 1300,  Sqrt[6], 1000]
};

reference = {
    0.0060994559262305161,
    0.0055186465519328747,
    0.006097076094553221,
    0.0055162667202560654,
    0.0061018357579078102,
    0.0055210263836096849,
    0.006097451701735793,
    0.0055206961724585694,
    0.0060950747439344347,
    0.0055183134669058234,
    0.006099828659537153,
    0.0055230788780113154,
    0.0062822480184172172,
    0.0062798983371975597,
    0.0062845976996368739,
    0.0062791938133536125,
    0.0062737898965503939
};

Print["Testing 2L SQCD correction to mb at non-universal SUSY masses ..."];

MapThread[TestCloseRel[Re @ N[b2l / (4 Pi)^4 //. #1], #2, 1*^-8]&, {points, reference}];

Print[];
PrintTestSummary[];
