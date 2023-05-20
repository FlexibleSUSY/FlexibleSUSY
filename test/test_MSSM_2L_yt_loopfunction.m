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

(* from twoloopbubble.cpp, used in SOFTSUSY *)
fpart[m1_, m2_, m3_, scale_] := Module[
   {Y = m2^2/m3^2, m3d = m3/scale},
   (
      2*Y*
        Log[Y]*(3 - 4*Log[m3d]) - (1 + Y)*(7 + 8*Log[m3d]^2 - 
          12*Log[m3d] + Pi^2/6) - 
       Y*Log[Y]^2 + (1 - Y)*(2*Log[Y]*Log[1 - Y] + 2*Li2[Y] - Pi^2/3)
      ) m3^2/2 /. Li2[a_] :> PolyLog[2, a]
   ];

(* from test2t.m, used in FlexibleSUSY *)
loopFunctions = {
    Hmine[mm1_, mm2_] :> 2*PolyLog[2, 1 - mm1/mm2] + 1/2*Log[mm1/mm2]^2, 

    fin[mm1_, mm2_, mmu_] :> 
    1/2*(-(mm1 + mm2)*(7 + Zeta[2]) + 
       6*(mm1*Log[mm1/mmu] + mm2*Log[mm2/mmu]) - 
       2*(mm1*Log[mm1/mmu]^2 + mm2*Log[mm2/mmu]^2) + 
       1/2*(mm1 + mm2)*Log[mm1/mm2]^2 + (mm1 - mm2)*Hmine[mm1, mm2])
};

simp = {mm1 -> m1^2, mm2 -> m2^2, mmu -> mu^2};

diff = FullSimplify[
    fpart[0, m1, m2, mu] - fin[mm1, mm2, mmu] //. loopFunctions //. simp,
    m1 > 0 && m2 > 0 && mu > 0
];

TestEquality[diff, 0];

PrintTestSummary[];
