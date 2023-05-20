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
Get[FileNameJoin[{"meta", "TwoLoopMSSM.m"}]];

t2l    = Get[FileNameJoin[{"meta", "MSSM", "tquark_2loop_strong.m"}]];
t2lqcd = Get[FileNameJoin[{"meta", "MSSM", "tquark_2loop_qcd.m"}]] /. GSY -> GS;
t2lDeg = GetDeltaMPoleOverMRunningMSSMSQCDDRbar2LUniversalMSUSY[];

colorCA = 3; colorCF = 4/3; Tf = 1/2; GS = g3;
MGl = mgl; MT = mt; SX = 2 mt Xt; s2t = SX / (mmst1 - mmst2);
fin[0, args__] := fin[args, mmu];

CollectTerms[expr_] := Collect[expr /. zt2 -> Zeta[2], {Log[__]}, Together];

loopFunctions = {
    Hmine[mm1_,mm2_] :> 2 * PolyLog[2, 1-mm1/mm2] + 1/2 * Log[mm1/mm2]^2,
    fin[mm1_,mm2_,mmu_] :>
       1/2 * ( - (mm1 + mm2) * ( 7 + Zeta[2] )
               + 6 * (mm1 * Log[mm1/mmu] + mm2 * Log[mm2/mmu])
               - 2 * (mm1 * Log[mm1/mmu]^2 + mm2 * Log[mm2/mmu]^2 )
               +1/2 * (mm1 + mm2) * Log[mm1/mm2]^2 + (mm1-mm2)*Hmine[mm1,mm2] )
};

(* ******* calculate limit ******* *)

t2lLimit = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. {
    mmst1  -> mmsusy + x * dst1,
    mmst2  -> mmsusy + x * dst2,
    mmgl   -> mmsusy + x * dst3
    }, {x,0,0}]];

(* ******* check difference ******* *)

diff = 
 Collect[(t2lLimit - t2lDeg (4 Pi)^4) //. {
     mmgl -> mgl^2, mgl -> MSUSY, mmu -> Q^2, mmsusy -> MSUSY^2,
     mmt -> mt^2
   }, {Xt, g3}, 
   Simplify[#, MSUSY > 0 && Q > 0 && mt > 0] &];

(* diff still contains power-suppressed terms O(mt^2/MSUSY^2) *)

diff = Normal[Series[diff /. mt -> x MSUSY, {x,0,0}]];

TestEquality[diff, 0];

PrintTestSummary[];
