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
Needs["ThreeLoopSM`", "ThreeLoopSM.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

Print["testing BetaSM[] ..."];

TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`hyperchargeCoupling]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`leftCoupling       ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`strongCoupling     ]], 5];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`UpYukawa           ]], 4];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`DownYukawa         ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[SARAH`ElectronYukawa     ]], 3];
TestEquality[Length[ThreeLoopSM`BetaSM[\[Lambda]                ]], 4];

Print["checking BetaSM[] and SM[] ..."];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];
SARAH`Start["SM"];
SARAH`CalcRGEs[ReadLists -> True, TwoLoop -> True,
               NoMatrixMultiplication -> False];

FindBetaFunction[lst_List, c_] :=
    Cases[lst, {c | c[__], beta__} :> {beta}][[1]];

TestBetaEquality[lst_, c_, loop_] :=
    Module[{sa, sm},
           sm = ThreeLoopSM`BetaSM[c][[loop]];
           sa = FindBetaFunction[lst, c][[loop]] /. {
               SARAH`MatMul[a__][_,_] :> Times[a],
               SARAH`MatMul[a__]    :> Times[a],
               SARAH`trace[a__]     :> Times[a],
               SARAH`trace[a_]      :> a
           } //. {
               SARAH`Adj[a_]        :> a,
               SARAH`UpYukawa       -> gt,
               SARAH`DownYukawa     -> gb,
               SARAH`ElectronYukawa -> g\[Tau],
               a_[Susyno`LieGroups`i1,SARAH`i2] :> a
           } /. {
               gt -> SARAH`UpYukawa[3,3],
               gb -> SARAH`DownYukawa[3,3],
               g\[Tau] -> SARAH`ElectronYukawa[3,3]
           };
           TestEquality[Expand[sm], Expand[sa]]
          ];

For[l = 1, l <= 2, l++,
    TestBetaEquality[SARAH`BetaGauge, SARAH`hyperchargeCoupling, l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`leftCoupling       , l];
    TestBetaEquality[SARAH`BetaGauge, SARAH`strongCoupling     , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`UpYukawa           , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`DownYukawa         , l];
    TestBetaEquality[SARAH`BetaYijk , SARAH`ElectronYukawa     , l];
    TestBetaEquality[SARAH`BetaLijkl, \[Lambda]                , l];
    TestBetaEquality[SARAH`BetaBij  , mu2                      , l];
   ];

PrintTestSummary[];
