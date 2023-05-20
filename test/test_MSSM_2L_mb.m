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

Print[];
Print["# non-universal SUSY masses: "];
Print[Re @ N[b2l / (4 Pi)^4 //. #]]& /@ {
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

Print[];
PrintTestSummary[];
