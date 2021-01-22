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
Needs["LoopFunctions`", "LoopFunctions.m"];

assumptions = m > 0 && m1 > 0 && m2 > 0 && p > 0 && Q > 0 && m2 != p;

RuleToExpansion[a_ -> b_, order_:0] := {a, b, 0};

RunTests[] := (
    Print["running all tests with $BPMZSign = ", $BPMZSign];
    Print["testing p -> 0 limit ..."];

    loopFunctions = {
        {A0[m,Q], A0[m,Q]},
        {B0[p,m,m,Q], B0[0,m,m,Q]},
        {B1[p,m,m,Q], B1[0,m,m,Q]},
        {B00[p,m,m,Q], B00[0,m,m,Q]},
        {B22[p,m,m,Q], B22[0,m,m,Q]},
        {B22tilde[p,m,m,Q], B22tilde[0,m,m,Q]},
        {F[p,m,m,Q], F[0,m,m,Q]},
        {G[p,m,m,Q], G[0,m,m,Q]},
        {H[p,m,m,Q], H[0,m,m,Q]}
        (* {C0[p1,p2,m1,m2,m3,Q], C0[0,0,m1,m2,m3,Q]}, *)
        (* {D0[p1,p2,p3,m1,m2,m3,m4,Q], D0[0,0,0,m1,m2,m3,m4,Q]}, *)
        (* {D27[p1,p2,p3,m1,m2,m3,m4,Q], D27[0,0,0,m1,m2,m3,m4,Q]} *)
    };

    For[i = 1, i <= Length[loopFunctions], i++,
        lffull = loopFunctions[[i,1]];
        lfzero = loopFunctions[[i,2]];
        Print["   testing ", lffull, " ..."];
        expr1 = Normal @ Series[lffull /. LFFull[], {p, 0, 0},
                                Assumptions :> assumptions];
        expr2 = lfzero /. LFZeroMomentum[];
        TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
       ];

    Print["testing m -> 0 limits ..."];

    loopFunctions = {
        {B0[p,m1,m2,Q], m1 -> 0},
        {B0[p,m1,m2,Q], m2 -> 0},
        {B0[0,m1,m2,Q], m1 -> 0},
        {B0[0,m1,m2,Q], m2 -> 0},
        {B1[p,m1,m2,Q], m1 -> 0},
        {B1[p,m1,m2,Q], m2 -> 0},
        {B1[0,m1,m2,Q], m1 -> 0},
        {B1[0,m1,m2,Q], m2 -> 0},
        {B00[p,m1,m2,Q], m1 -> 0},
        {B00[p,m1,m2,Q], m2 -> 0},
        {B00[0,m1,m2,Q], m1 -> 0},
        {B00[0,m1,m2,Q], m2 -> 0},
        {B11[p,m1,m2,Q], m1 -> 0},
        {B11[p,m1,m2,Q], m2 -> 0},
        {B22[p,m1,m2,Q], m1 -> 0},
        {B22[p,m1,m2,Q], m2 -> 0},
        {B22[0,m1,m2,Q], m1 -> 0},
        {B22[0,m1,m2,Q], m2 -> 0},
        {B22tilde[p,m1,m2,Q], m1 -> 0},
        {B22tilde[p,m1,m2,Q], m2 -> 0},
        {B22tilde[0,m1,m2,Q], m1 -> 0},
        {B22tilde[0,m1,m2,Q], m2 -> 0},
        {F[p,m1,m2,Q], m1 -> 0},
        {F[p,m1,m2,Q], m2 -> 0},
        {F[0,m1,m2,Q], m1 -> 0},
        {F[0,m1,m2,Q], m2 -> 0},
        {G[p,m1,m2,Q], m1 -> 0},
        {G[p,m1,m2,Q], m2 -> 0},
        {G[0,m1,m2,Q], m1 -> 0},
        {G[0,m1,m2,Q], m2 -> 0},
        {H[p,m1,m2,Q], m1 -> 0},
        {H[p,m1,m2,Q], m2 -> 0},
        {H[0,m1,m2,Q], m1 -> 0},
        {H[0,m1,m2,Q], m2 -> 0}
    };

    For[i = 1, i <= Length[loopFunctions], i++,
        lf = First[loopFunctions[[i]]];
        limit = Drop[loopFunctions[[i]], 1];
        Print["   testing ", lf, " in the limit ", limit," ..."];
        expr1 = Normal @ Series[lf /. LFFull[], Sequence @@ (RuleToExpansion /@ limit), Assumptions :> assumptions];
        expr2 = (lf /. limit) /. LFFull[];
        TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
       ];

    (* Print["testing B0[] function ..."]; *)

    (* TestEquality[Simplify[ *)
    (*     (B0[1, 2, 3, 4] /. LFFull[]) - *)
    (*     LoopFunctions`Private`B0integral[1, 2, 3, 4]], *)
    (*              0 *)
    (* ]; *)

    Print["testing divergences ..."];

    loopFunctions = {
        A0[m,Q],
        B0[p,m1,m2,Q],
        B1[p,m1,m2,Q],
        B00[p,m1,m2,Q],
        B11[p,m1,m2,Q],
        B22[p,m1,m2,Q],
        B22tilde[p,m1,m2,Q],
        F[p,m1,m2,Q],
        G[p,m1,m2,Q],
        H[p,m1,m2,Q]
    };

    For[i = 1, i <= Length[loopFunctions], i++,
        Print["   testing ", loopFunctions[[i]], " ..."];
        expr1 = 1/2 Q D[loopFunctions[[i]] /. LFFull[], Q];
        expr2 = Coefficient[loopFunctions[[i]] /. LFDivergence[], Delta];
        TestEquality[FullSimplify[expr1 - expr2, assumptions], 0];
       ];

    loopFunctions = {
        C0[p1,p2,m1,m2,m3,Q],
        D0[p1,p2,p3,m1,m2,m3,m4,Q],
        D27[p1,p2,p3,m1,m2,m3,m4,Q]
    };

    For[i = 1, i <= Length[loopFunctions], i++,
        Print["   testing ", loopFunctions[[i]], " ..."];
        expr2 = loopFunctions[[i]] /. LFDivergence[];
        TestEquality[expr2, 0];
       ];

    Print["testing logarithms ..."];

    loopFunctions = {
        A0[m, mu],
        B0[p, m1, m2, mu],
        B1[p, m1, m2, mu],
        B00[p, m1, m2, mu],
        B11[p, m1, m2, mu],
        B22[p, m1, m2, mu],
        B22tilde[p, m1, m2, mu],
        C0[p1, p2, m1, m2, m3, mu],
        D0[p1, p2, p3, m1, m2, m3, m4, mu],
        D27[p1, p2, p3, m1, m2, m3, m4, mu],
        F[p, m1, m2, mu],
        G[p, m1, m2, mu],
        H[p, m1, m2, mu]
    };

    divs = loopFunctions /. LFDivergence[] /. Delta -> Log[mu^2];
    logs = loopFunctions /. LFScaleDependence[];

    divs = (mu D[#, mu]) & /@ divs;
    logs = (mu D[#, mu]) & /@ logs;

    TestEquality[Simplify[divs - logs], Table[0,{i,1,Length[loopFunctions]}]];
);

$BPMZSign = 1;
RunTests[];

$BPMZSign = -1;
RunTests[];

PrintTestSummary[];
