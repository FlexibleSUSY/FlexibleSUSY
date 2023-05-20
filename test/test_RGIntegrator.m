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
Needs["RGIntegrator`", "RGIntegrator.m"];

On[Assert];

betas = {
    {},
    {{g}},
    {{g}, {l}},
    {{g, h a g}},
    {{g, h a g, h^2 b g^2}},
    {{g, h a g, h^2 b g^2 l^2},
     {l, h (c l + d g), h^2 (e l^2 + f g^2 + k l g)}}
};

Print["Testing running up and back down ..."];

For[i = 1, i <= Length[betas], i++,
    For[p = 1, p <= Length[betas[[i]]], p++,
        par = First[betas[[i, p]]];
        lo = Length[betas[[i, p]]] - 1;
        Print["  Running ", par, " using ", InputForm[betas[[i]]], " (",lo,"L)"];
        runningPar = FullSimplify[
            Normal@Series[(par[Q0] /. RGIntegrate[betas[[i]], Q0, Q1] /. 
                           RGIntegrate[betas[[i]], Q1, Q0]),
                          {h, 0, lo}], Assumptions :> Q0 > 0 && Q1 > 0];
        TestEquality[runningPar, par[Q0]];
       ];
   ];

betas = {
    {{g}},
    {{g}, {l}},
    {{g, h a g}},
    {{g, h a g, h^2 b g^2}},
    {{g, h a g, h^2 b g^2 l^2},
     {l, h (c l + d g), h^2 (e l^2 + f g^2 + k l g)}}
};

Print["Testing running up and back down with higher loop order ..."];

For[i = 1, i <= Length[betas], i++,
    For[p = 1, p <= Length[betas[[i]]], p++,
        par = First[betas[[i, p]]];
        lo = Length[betas[[i, p]]] - 1 + 1;
        Print["  Running ", par, " using ", InputForm[betas[[i]]], " (",lo,"L)"];
        runningPar = FullSimplify[
            Normal@Series[(par[Q0] /. RGIntegrate[betas[[i]], Q0, Q1, loopOrder -> lo] /. 
                           RGIntegrate[betas[[i]], Q1, Q0, loopOrder -> lo]),
                          {h, 0, lo}], Assumptions :> Q0 > 0 && Q1 > 0];
        TestEquality[runningPar, par[Q0]];
       ];
   ];

Print["Testing scale dependence of solution ..."];

For[i = 1, i <= Length[betas], i++,
    pars = #[[1]]& /@ betas[[i]];
    parsQ0 = #[[1]][Q0]& /@ betas[[i]];
    sols = parsQ0 /. RGIntegrate[betas[[i]], Q0, Q];
    deriv = (Q0 D[#, Q0]& /@ sols) /. Q -> Q0;
    beta = Total /@ (Drop[#,1]& /@ (betas[[i]] /. (Rule[#,#[Q0]]& /@ pars)));

    Print["  Testing ", InputForm[beta]];

    TestEquality[Simplify[deriv - beta], Table[0,{i,1,Length[betas[[i]]]}]];
   ];

PrintTestSummary[];
