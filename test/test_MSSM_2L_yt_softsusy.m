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

colorCA = 3; colorCF = 4/3; Tf = 1/2; GS = g3; zt3 = zeta[3]; zt2 = Zeta[2]; log2 = Log[2];
cst = s2t/2/snt; mmt = mt^2; mmst1 = mst1^2; mmst2 = mst2^2; mmsusy = msusy^2;
mmgl = mgl^2; MGl = mgl; MT = mt; mmu = scale^2; mu = scale;

fin[0, mm1_, mm2_] := fpart[0, Sqrt[mm1], Sqrt[mm2], scale];
fpart[0, mst2, mst1, scale] = fpart[0, mst1, mst2, scale];
fpart[0, mgl, mst1, scale] = fpart[0, mst1, mgl, scale];
fpart[0, mgl, mst2, scale] = fpart[0, mst2, mgl, scale];
fpart[0, msusy, mgl, scale] = fpart[0, mgl, msusy, scale];
fpart[0, msusy, mst1, scale] = fpart[0, mst1, msusy, scale];
fpart[0, msusy, mst2, scale] = fpart[0, mst2, msusy, scale];

exFS = Get[FileNameJoin[{"meta","MSSM","tquark_2loop_strong.m"}]] / g3^4;
exSS = Get[FileNameJoin[{"meta","MSSM","dmtas2.m"}]] / g3^4;

diff = Collect[
    Refine[exFS - exSS, 
           mst1 > 0 && mst2 > 0 && mt > 0 && mgl > 0 && msusy > 0 && scale > 0
    ],
    {Log[__], s2t}, Together
];

TestEquality[diff, 0];

PrintTestSummary[];
