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

(* This test compares the implementation of arxiv:0901.2065 with the
   implementation of arxiv:hep-ph/9307201 and arxiv:1508.00576 *)

Needs["TestSuite`", "TestSuite.m"];
Needs["THDMThresholds1L`", FileNameJoin[{Directory[], "meta", "THDM", "Thresholds_1L_full.m"}]];

Print["testing THDM loop order ..."];

TestEquality[GetTHDMThresholds1L[loopOrder -> {0,0}],
             Table[0, {i,1,7}]];

Print["testing THDM tree-level ..."];

(* tree-level lambda couplings in convention of
   arxiv:hep-ph/9307201 and arxiv:1508.00576 *)
lambdaTree = {
    1/4 (gY^2 + g2^2),
    1/4 (gY^2 + g2^2),
    1/4 (-gY^2 + g2^2),
    -1/2 g2^2,
    0, 0, 0
};

TestEquality[Simplify[lambdaTree - GetTHDMThresholds1L[loopOrder -> {1,0}]],
             Table[0, {i,1,7}]];

Print["testing THDM flags ..."];

zeroFlags = { flagSferm -> 0, flagIno -> 0, flagdg -> 0 };

TestEquality[Simplify[lambdaTree - GetTHDMThresholds1L[flags -> zeroFlags]],
             Table[0, {i,1,7}]];

approx = {
    Yu[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Yd[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Ye[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Tu[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Td[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Te[i_, j_] :> 0 /; i < 3 || j < 3 || i != j,
    Yu[3, 3] -> ht,
    Yd[3, 3] -> hb,
    Ye[3, 3] -> htau,
    Tu[3, 3] -> At ht,
    Td[3, 3] -> Ab hb,
    Te[3, 3] -> Atau htau,
    mse[_] :> MSUSY,
    msu[_] :> MSUSY,
    msd[_] :> MSUSY,
    msq[_] :> MSUSY,
    msl[_] :> MSUSY,
    M1 -> Mu,
    M2 -> Mu,
    Abs[p_] :> p,
    Conjugate[p_] :> p,
    Re[p_] :> p
};

result = GetTHDMThresholds1L[];

TestEquality[FreeQ[result, Undef] && FreeQ[result, Null], True];

result = result //. approx //. GetTHDMThresholds1LLoopFunctions[] /. Q -> MSUSY // Expand;

(* collect only terms involving Yukawa couplings *)
result = Simplify[Plus @@ Select[List @@ Expand[#],
                                 Function[x,!FreeQ[x,ht] || !FreeQ[x,hb] || !FreeQ[x,htau]]]]& /@ result

Print["testing THDM threshold corrections[] ..."];

Get["model_files/THDMIIMSSMBC/FlexibleSUSY.m.in"];

renameRules = {
    AtauInput -> Atau,
    AtInput -> At,
    AbInput -> Ab,
    Yu[3,3] -> ht,
    Yd[3,3] -> hb,
    Ye[3,3] -> htau,
    MuInput -> Mu,
    \[Mu] -> Mu
};

gRules = {
    g1 -> gY Sqrt[5/3],
    GUTNormalization[g1] -> Sqrt[3/5]
};

deltaLambdaTh = {
    deltaLambda1th1L, deltaLambda2th1L, deltaLambda3th1L, 
    deltaLambda4th1L, deltaLambda5th1L, deltaLambda6th1L, 
    deltaLambda7th1L
} //. renameRules //. gRules //. approx;

deltaLambdaPhi = {
    deltaLambda1Phi1L, deltaLambda2Phi1L, deltaLambda3Phi1L, 
    deltaLambda4Phi1L, deltaLambda5Phi1L, deltaLambda6Phi1L, 
    deltaLambda7Phi1L
} //. renameRules //. gRules //. approx;

deltaLambdaFull = deltaLambdaTh + deltaLambdaPhi;

lamThDiff  = Simplify[result - deltaLambdaFull];

TestEquality[lamThDiff, Table[0, {i,1,7}]];

Print["testing MS-DR conversion ..."];

(* tag MS-DR conversion terms and filter them out *)
conv = Coefficient[
    GetTHDMThresholds1L[loopOrder -> {0, k (4 Pi)^2},
                        flags -> Join[{flagMSDRg2 -> tag, flagMSDRlam -> tag},
                                      GetTHDMThresholds1LFlags[]]],
    tag
] /. gRules;

(* [1805.00867] Eqs.(104) *)
lit = {
    -1/12 k (7 g2^4 + 6 g2^2 gY^2 + 3 gY^4),
    -1/12 k (7 g2^4 + 6 g2^2 gY^2 + 3 gY^4),
    -1/12 k (7 g2^4 - 6 g2^2 gY^2 + 3 gY^4),
    -1/3  k g2^2 (g2^2 + 3 gY^2),
    0, 0, 0
};

TestEquality[conv - lit // Simplify, Table[0, {i,1,7}]];

PrintTestSummary[];
