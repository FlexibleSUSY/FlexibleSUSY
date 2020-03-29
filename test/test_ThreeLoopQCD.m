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

Needs["SARAH`"];
Needs["TestSuite`", "TestSuite.m"];
Needs["TwoLoopQCD`", "TwoLoopQCD.m"];
Needs["ThreeLoopQCD`", "ThreeLoopQCD.m"];
Needs["ThreeLoopSM`", "ThreeLoopSM.m"];

FlexibleSUSY`$flexiblesusyMetaDir = FileNameJoin[{Directory[], "meta"}];

Start["SM"];

Print["Testing 1L arxiv:hep-ph/9912391 vs. arxiv:hep-ph/9803493 ..."];

m1 = M/(1 + 
    h GetDeltaMPoleOverMRunningQCDOneLoopMSbar[TopQuark, Q]);

m1 = Simplify[Normal[Series[m1, {h, 0, 1}]] /. h -> 1];

m2 = M GetMTopMSbarOverMTopPole[{1, 1, 0, 0}];

TestEquality[Simplify[m1 - m2], 0];

Print["Testing 2L arxiv:hep-ph/9912391 vs. arxiv:hep-ph/9803493 ..."];

m1 = M/(1 + 
    h GetDeltaMPoleOverMRunningQCDOneLoopMSbar[TopQuark, Q] + 
    h^2 GetDeltaMPoleOverMRunningQCDTwoLoopMSbar[TopQuark, Q]);

m1 = Simplify[Normal[Series[m1, {h, 0, 2}]] /. h -> 1];

m2 = M GetMTopMSbarOverMTopPole[{1, 1, 1, 0}];

TestEquality[Simplify[(m1 - m2) /. Log[x_/Q^2] :> -Log[Q^2/x]], 0];

Print["Testing 2L arxiv:hep-ph/9803493 vs. arxiv:hep-ph/9912391 ..."];

M1 = m (1 +
        h   GetDeltaMPoleOverMRunningQCDOneLoopMSbar[TopQuark, Q] +
        h^2 GetDeltaMPoleOverMRunningQCDTwoLoopMSbar[TopQuark, Q]);

M2 = m / (GetMTopMSbarOverMTopPole[{1, 0, 0  , 0}] +
          GetMTopMSbarOverMTopPole[{0, h, 0  , 0}] +
          GetMTopMSbarOverMTopPole[{0, 0, h^2, 0}])

M2 = Simplify[Normal[Series[M2, {h, 0, 2}]]];

TestEquality[Simplify[(M1 - M2) /. h -> 1 /. Log[x_/Q^2] :> -Log[Q^2/x]], 0];

Print["Testing 3L renormalization scale invariance ..."];

gRules = {
    g3 -> Sqrt[as[Q] 4 Pi],
    Yu[__] -> 0,
    Yd[__] -> 0,
    g1 -> 0,
    g2 -> 0,
    FlexibleSUSY`M[Fu] -> m[Q]
};

ass = {Q > 0, FlexibleSUSY`M[Fu] > 0, m[Q] > 0};

(* beta functions of g3 *)
betag3 = ThreeLoopSM`BetaSM[SARAH`strongCoupling];
betag3[[1]] *= 1/(4 Pi)^2;
betag3[[2]] *= 1/(4 Pi)^4;
betag3[[3]] *= 1/(4 Pi)^6;

(* beta functions of alpha_S *)
betaAlphaS = Simplify[(g3/(2 Pi) betag3) /. gRules];

(* beta function of MS-bar up-quark masses *)
betam = -2 as[Q] m[Q] {
     1/Pi,
     as[Q] (202/3 - (20 Nf)/9)/(16 Pi^2),
     as[Q]^2 (1249 - ((2216 Nf)/27 + 160 Zeta[3] Nf/3) -
         140 Nf^2/81)/(64 Pi^3),
     as[Q]^3 ((1379027 + 6048*Pi^4 - 2277312*Zeta[3] +
               561600*Zeta[5])/162)/(4 Pi)^4
} /. Nf -> 6;

(* define derivative of alpha_S *)
Derivative[1][as][Q] = (
    h^1 betaAlphaS[[1]] +
    h^2 betaAlphaS[[2]] +
    h^3 betaAlphaS[[3]]
)/Q;

(* define derivative of MS-bar quark masses *)
Derivative[1][m][Q] = (
    h^1 betam[[1]] +
    h^2 betam[[2]] +
    h^3 betam[[3]] +
    h^4 betam[[4]]
)/Q;

(* pole mass up to 3L order *)
M = Simplify[
    m[Q] GetMTopPoleOverMTopMSbar[{1,h,h^2,h^3}] /. gRules
    , ass];

deriv = D[M, Q];

deriv = FullSimplify[
    Normal[Series[deriv, {h, 0, 3}]],
    ass
];

TestEquality[deriv, 0];

Print["Testing cached 3-loop expression Mt/mt ..."];

Get[FileNameJoin[{FlexibleSUSY`$flexiblesusyMetaDir, "SM", "mf_3loop_qcd.m"}]];

TestEquality[
    FullSimplify[
        (MfOvermf /. { L -> Log[mf^2/Q^2] }) -
        (GetMTopPoleOverMTopMSbar[{1,h^1,h^2,h^3}, TopQuark, Q, NH, NL] //. {
            h -> k (4Pi)^2, FlexibleSUSY`M[Fu] -> mf,
            Log[Q^2/mf^2] -> -Log[mf^2/Q^2]
         })
    ]
    , 0
];

Print["Testing cached 3-loop expression mt/Mt ..."];

TestEquality[
    FullSimplify[
        (mfOverMf /. { L -> Log[mf^2/Q^2] }) -
        (GetMTopMSbarOverMTopPole[{1,h^1,h^2,h^3}, TopQuark, Q, NH, NL] //. {
            h -> k (4Pi)^2, FlexibleSUSY`M[Fu] -> mf,
            Log[Q^2/mf^2] -> -Log[mf^2/Q^2]
         })
    ]
    , 0
];

Print["Testing cached 4-loop expression Mt/mt ..."];

Get[FileNameJoin[{FlexibleSUSY`$flexiblesusyMetaDir, "SM", "mt_4loop_qcd.m"}]];

TestEquality[
    FullSimplify[
        Normal[Series[MtOvermt - MfOvermf, {k,0,3}]] //. {
            NH -> 1,
            NL -> 5
        }
    ]
    , 0
];

Print["Testing cached 4-loop expression mt/Mt ..."];

TestEquality[
    FullSimplify[
        Normal[Series[mtOverMt - mfOverMf, {k,0,3}]] //. {
            NH -> 1,
            NL -> 5
        }
    ]
    , 0
];

Print["Testing 4L renormalization scale invariance ..."];

(* pole mass up to 4L order *)
M = Simplify[m[Q] MtOvermt /. gRules //. {
        k -> h/(4Pi)^2,
        L -> Log[m[Q]^2/Q^2]
    }, ass];

deriv = D[M, Q];

deriv = FullSimplify[
    Normal[Series[deriv, {h, 0, 4}]],
    ass
];

TestEquality[deriv, 0];


PrintTestSummary[];
