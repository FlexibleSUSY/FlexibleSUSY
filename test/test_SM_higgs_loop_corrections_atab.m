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
SMdir = FileNameJoin[{Directory[], "meta", "SM"}];

yt = gt;
yb = gb;

DMh20L = \[Lambda] v^2;

(* 1-loop contribution, derived from the eff. potential *)
DMh21L = (
    + 2 yt^2 v^2 (\[Lambda] (2+3Lt)/2 - 3yt^2 Lt)
    + 2 yb^2 v^2 (\[Lambda] (2+3Lb)/2 - 3yb^2 Lb)
) /. { Lt -> Log[t/Q^2], Lb-> Log[b/Q^2] } /.
     { b -> mb^2, t -> mt^2 } /.
     { mt -> v yt/Sqrt[2], mb -> v yb/Sqrt[2] };

(* 2-loop contributions, derived from the eff. potential, which are to be tested *)
Mh22lEffPot = mt^6/v^4 (
    - 8 (6+\[Pi]^2 + 3 Log[mt^2/Q^2] (-7 + 3 Log[mt^2/Q^2]))
    + 24 (\[Pi]^2 + Log[mt^2/Q^2] + 3 Log[mt^2/Q^2]^2) (mb/mt)^2
    + 12 (-15 - 2 \[Pi]^2 + 2 Log[mt^2/Q^2] + 6 Log[mt^2/Q^2]^2
          + 8 Log[mb/mt] (2 + 3 Log[mt^2/Q^2])) (mb/mt)^4
    + 4/3 (49 + 6 \[Pi]^2 - 144 Log[mb/mt]^2 + 24 Log[mb/mt] (5 - 9 Log[mt^2/Q^2])
           + 18 (7 - 3 Log[mt^2/Q^2]) Log[mt^2/Q^2]) (mb/mt)^6
);

(* 2-loop contributions which come from momentum iteration *)
Mh22lpIter = 1/v^4 (
    + mt^2 2 (2 + 3 Log[mt^2/Q^2])
    + mb^2 2 (2 + 3 Log[mb^2/Q^2])
) (-24 mt^4 Log[mt^2/Q^2] - 24 mb^4 Log[mb^2/Q^2]);

DMh22L = Mh22lpIter + Mh22lEffPot;

Mh22L = DMh20L + h DMh21L + h^2 DMh22L;

toRunningPars = {
    g1 -> 0,
    g2 -> 0,
    g3 -> 0,
    g\[Tau] -> 0,
    gb -> gb[Q],
    gt -> gt[Q],
    \[Lambda] -> \[Lambda][Q],
    v -> v[Q],
    m2 -> m2[Q],
    mt -> gt[Q] v[Q] / Sqrt[2],
    mb -> gb[Q] v[Q] / Sqrt[2]
};

simp = {
    g3[Q] -> 0,
    gb[Q] -> mb/v Sqrt[2],
    gt[Q] -> mt/v Sqrt[2],
    \[Lambda][Q] -> 0,
    v[Q] -> v,
    m2[Q] -> m2
};

ass = { mt > 0, mb > 0, Q > 0, gt > 0, gb > 0, v > 0 };

(* beta functions *)
betagb = Get[FileNameJoin[{SMdir, "beta_gb.m"}]] /. toRunningPars;
betagt = Get[FileNameJoin[{SMdir, "beta_gt.m"}]] /. toRunningPars;
betala = Get[FileNameJoin[{SMdir, "beta_lambda.m"}]] /. toRunningPars;
betav  = Get[FileNameJoin[{SMdir, "beta_v.m"}]] /. toRunningPars;
betam2 = Get[FileNameJoin[{SMdir, "beta_m2.m"}]] /. toRunningPars;

(* define derivative of gb *)
Derivative[1][gb][Q] = (
    h^1 betagb[[1]] +
    h^2 betagb[[2]]
)/Q;

(* define derivative of gt *)
Derivative[1][gt][Q] = (
    h^1 betagt[[1]] +
    h^2 betagt[[2]]
)/Q;

(* define derivative of lambda *)
Derivative[1][\[Lambda]][Q] = (
    h^1 betala[[1]] +
    h^2 betala[[2]]
)/Q;

(* define derivative of VEV *)
Derivative[1][v][Q] = (
    h^1 betav[[1]] +
    h^2 betav[[2]]
)/Q;

(* define derivative of m^2 *)
Derivative[1][m2][Q] = (
    h^1 betam2[[1]] +
    h^2 betam2[[2]]
)/Q;

deriv = D[Mh22L /. toRunningPars, Q] //. simp;

deriv = Collect[
    Normal[Series[deriv, {h, 0, 2}]],
    {h, gt}, FullSimplify[#, ass]&
];


Print["Testing 0L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 0], 0];

Print["Testing 1L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 1], 0];

Print["Testing 2L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 2], 0];

PrintTestSummary[];
