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
Needs["Utils`", "Utils.m"];
Needs["CCompilerDriver`"];

Get[FileNameJoin[{Directory[], "meta", "TwoLoopMSSM.m"}]];

Print["Comparing numerically with Pietro Slavich's routines ... "];

points = {
   {mt -> 175, M3 -> 1000, mst1 -> 1001, mst2 -> 2001   , sinTheta -> 0.2, Q -> 900, Mu -> 100, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1},
   {mt -> 175, M3 -> 2000, mst1 -> 1001, mst2 -> 2001   , sinTheta -> 0.2, Q -> 900, Mu -> 100, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1},
   {mt -> 175, M3 -> 2000, mst1 -> 2001, mst2 -> 1001   , sinTheta -> 0.2, Q -> 900, Mu -> 100, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1},
   {mt -> 175, M3 -> 2000, mst1 -> 1000, mst2 -> 2000.01, sinTheta -> 0  , Q -> 900, Mu -> 100, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1},
   {mt -> 175, M3 -> 2000, mst1 -> 1000, mst2 -> 2000.01, sinTheta -> 0  , Q -> 900, Mu ->   0, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1},
   {mt -> 175, M3 -> 2000, mst1 -> 1000, mst2 -> 1000.01, sinTheta -> 0  , Q -> 900, Mu -> 100, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1},
   {mt -> 175, M3 -> 2000, mst1 -> 1000, mst2 -> 1000.01, sinTheta -> 0  , Q -> 900, Mu ->   0, TanBeta -> 10, v -> 245, g3 -> 0.118, signMu -> 1}
};

randomPoints = BlockRandom[
    SeedRandom[1];
    {mt -> RandomReal[{100,200}],
     M3 -> RandomReal[{100,10000}],
     mst1 -> RandomReal[{100,10000}],
     mst2 -> RandomReal[{100,10000}],
     sinTheta -> RandomReal[{-1,1}],
     Q -> RandomReal[{10,10000}],
     Mu -> RandomReal[{0,10000}],
     TanBeta -> RandomReal[{1,100}],
     v -> RandomReal[{240,250}],
     g3 -> RandomReal[{0.1,0.2}],
     signMu -> RandomChoice[{-1,1}]}& /@ Table[i, {i,1,100}]
];

points = Join[points, randomPoints];

CalculatePointFromAnalyticExpr[point_] :=
    Module[{s2t, yt, at, deltaMh, deltaMa, pars},
           s2t = Sin[2 ArcSin[sinTheta /. point]];
           yt = (Sqrt[2] mt/(v Sin[ArcTan[TanBeta]])) /. point;
           at = ((mst1^2 - mst2^2) s2t/(2 mt) - signMu Mu/TanBeta) /. point;
           pars = Join[
               { sin2Theta -> s2t, ht -> yt, At -> at },
               point
           ];
           deltaMh = Simplify @ Re @ N[
               GetMSSMCPEvenHiggsLoopMassMatrix[
                   loopOrder -> {0,0,1}, parameters -> pars, corrections -> {1,0}]];
           deltaMa = Simplify @ Re @ N[
               GetMSSMCPOddHiggsLoopMass[
                   loopOrder -> {0,0,1}, parameters -> pars, corrections -> {1,0}]];
           {deltaMh[[1, 1]], deltaMh[[2, 2]], deltaMh[[1, 2]], deltaMa}
          ];

CalculatePointNumerical[point_] :=
    Module[{progr, exec, fortComp},
           progr = "
#include <stdio.h>
#include <math.h>
#include \"mssm_twoloophiggs.h\"

double sqr(double x) { return x*x; }

int calc_Sij(double* S11, double* S22, double* S12)
{
   int OS = 0;
   double mt2 = sqr(" <> ToString[mt /. point] <> ");
   double mg = " <> ToString[M3 /. point] <> ";
   double mst12 = sqr(" <> ToString[mst1 /. point] <> ");
   double mst22 = sqr(" <> ToString[mst2 /. point] <> ");
   double st = " <> ToString[sinTheta /. point] <> ";
   double ct = sqrt(1. - sqr(st));
   double q2 = sqr(" <> ToString[Q /. point] <> ");
   double mu = " <> ToString[signMu Mu /. point] <> ";
   double tb = " <> ToString[TanBeta /. point] <> ";
   double v2 = sqr(" <> ToString[v /. point] <> ");
   double g3 = " <> ToString[g3 /. point] <> ";

   return dszhiggs_(
      &mt2, &mg, &mst12, &mst22, &st, &ct,
      &q2, &mu, &tb, &v2, &g3, &OS, S11, S22, S12);
}

int calc_A(double* A)
{
   double mt2 = sqr(" <> ToString[mt /. point] <> ");
   double mg = " <> ToString[M3 /. point] <> ";
   double mst12 = sqr(" <> ToString[mst1 /. point] <> ");
   double mst22 = sqr(" <> ToString[mst2 /. point] <> ");
   double st = " <> ToString[sinTheta /. point] <> ";
   double ct = sqrt(1. - sqr(st));
   double q2 = sqr(" <> ToString[Q /. point] <> ");
   double mu = " <> ToString[signMu Mu /. point] <> ";
   double tb = " <> ToString[TanBeta /. point] <> ";
   double v2 = sqr(" <> ToString[v /. point] <> ");
   double g3 = " <> ToString[g3 /. point] <> ";

   return dszodd_(
      &mt2, &mg, &mst12, &mst22, &st, &ct,
      &q2, &mu, &tb, &v2, &g3, A);
}

int main(){
   double S11 = 0., S22 = 0., S12 = 0., A = 0.;
   calc_Sij(&S11, &S22, &S12);
   calc_A(&A);
   printf(\"{%g, %g, %g, %g}\\n\", S11, S22, S12, A);
   return 0;
}
";
           fortComp =
              Select[
                 StringSplit@Import[File["flexiblesusy-config"]],
                 StringMatchQ[#, "fc=" ~~ __] &
              ];
           If[Length[fortComp]===1,
              fortComp = Last@StringSplit[First@fortComp, "="],
              Print["Could not identify a Fortran compiler"]; Quit[1];
           ];
           fortComp = StringReplace[fortComp, "\"" -> ""];
           exec = CreateExecutable[
               progr, "exec",
               "CompileOptions" -> StringJoin[Riffle[
                   {"-I" <> FileNameJoin[{Directory[], "src"}],
                    "-I" <> FileNameJoin[{Directory[], "model_specific", "MSSM_higgs"}]},
                   " "]],
               "Libraries" -> {FileNameJoin[{Directory[], "model_specific", "MSSM_higgs", "libmodel_specific_MSSM_higgs.a"}],
                               FileNameJoin[{Directory[], "src", "libflexisusy.a"}],
                               Switch[fortComp,
                                  "gfortran", "gfortran",
                                  "ifort",    "imf",
                                  _, Print["Unknown Fortran compiler " <> fortComp]; Quit[1]
                               ]}
               (* , "ShellOutputFunction"->Print, "ShellCommandFunction"->Print *)
           ];
           If[exec === $Failed,
              Print["Error: cannot create executable"];
              Return[{0,0,0,0}];
             ];
           ToExpression[
               StringReplace[
                   Import["!" <> QuoteFile[exec], "String"], {"e+" :> "*^", "e-" :> "*^-"}
               ]
           ]
          ];

For[i = 1, i <= Length[points], i++,
    pAnalytic = CalculatePointFromAnalyticExpr[points[[i]]];
    pNumeric = CalculatePointNumerical[points[[i]]];

    max = N[10^(-3)];
    relDiff = Max[MaxRelDiff /@ Zip[pAnalytic, pNumeric]];
    test = relDiff < max;

    Print["Running point ", i];
    Print["analytic: ", pAnalytic];
    Print["numeric : ", pNumeric];
    Print[""];

    If[!test,
       Print["Error: rel. difference = ", InputForm[relDiff], " > ", InputForm[max]];
       Print["Point was: ", points[[i]]];
      ];

    TestEquality[test, True];
   ];

Print["done"];
Print[""];

Print["Testing limits ... "];

points = {
    {Mu -> 0, At -> 0},
    {Mu -> 0, At -> 0, sin2Theta -> 0},
    {Mu -> 0, At -> 0, sin2Theta -> 0, mst1 -> mst2},
    {sin2Theta -> 0},
    {sin2Theta -> 0, mst1 -> mst2}
};

For[i = 1, i <= Length[points], i++,
    ex = GetMSSMCPEvenHiggsLoopMassMatrix[loopOrder -> {0, 0, 1}, parameters -> points[[i]], corrections -> {1,0}];
    re = ReplaceStopMasses[parameters -> points[[i]]];
    TestEquality[ex[[1,1]] =!= Indeterminate, True];
    TestEquality[ex[[1,2]] =!= Indeterminate, True];
    TestEquality[ex[[2,1]] =!= Indeterminate, True];
    TestEquality[ex[[2,2]] =!= Indeterminate, True];
    TestEquality[(mst1 /. re) =!= Indeterminate, True];
    TestEquality[(mst2 /. re) =!= Indeterminate, True];
    TestEquality[(sin2Theta /. re) =!= Indeterminate, True];

    ex = GetMSSMCPOddHiggsLoopMass[loopOrder -> {0, 0, 1}, parameters -> points[[i]], corrections -> {1,0}];
    TestEquality[ex =!= Indeterminate, True];
   ];

Print["done"];
Print[""];

PrintTestSummary[];
