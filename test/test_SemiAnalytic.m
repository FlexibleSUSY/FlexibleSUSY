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
Needs["SemiAnalytic`", "SemiAnalytic.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["MSSM"];
SARAH`CalcRGEs[];

betas = { SARAH`BetaWijkl, SARAH`BetaYijk, SARAH`BetaMuij,
          SARAH`BetaLi, SARAH`BetaGauge, SARAH`BetaVEV,
          SARAH`BetaQijkl, SARAH`BetaTijk, SARAH`BetaBij,
          SARAH`BetaLSi, SARAH`Betam2ij, SARAH`BetaMi,
          SARAH`BetaDGi, SARAH`BetaLijkl };

If[Head[#] === Symbol && !ValueQ[#], Set[#,{}]]& /@ betas;
If[!ValueQ[SARAH`Gij] || Head[SARAH`Gij] =!= List,
   SARAH`Gij = {};
  ];

susyBetaFunctions = { SARAH`BetaLijkl, SARAH`BetaWijkl,
                      SARAH`BetaYijk, SARAH`BetaMuij,
                      SARAH`BetaLi, SARAH`BetaGauge,
                      SARAH`BetaVEV };

susyBreakingBetaFunctions = { SARAH`BetaQijkl, SARAH`BetaTijk ,
                              SARAH`BetaBij, SARAH`BetaLSi,
                              SARAH`Betam2ij, SARAH`BetaMi,
                              SARAH`BetaDGi };

allModelParameters = ((#[[1]])& /@ Join[Join @@ susyBetaFunctions,
                                        Join @@ susyBreakingBetaFunctions]) /.
                               a_[Susyno`LieGroups`i1] :> a /.
                               a_[Susyno`LieGroups`i1,SARAH`i2] :> a;

Parameters`SetModelParameters[allModelParameters];

allInputParameters = {{m0, {"MINPAR", 1}, CConversion`ScalarType[CConversion`realScalarCType]},
                      {m12, {"MINPAR", 2}, CConversion`ScalarType[CConversion`realScalarCType]},
                      {TanBeta, {"MINPAR", 3}, CConversion`ScalarType[CConversion`realScalarCType]},
                      {Azero, {"MINPAR", 4}, CConversion`ScalarType[CConversion`realScalarCType]}};

Parameters`SetInputParameters[allInputParameters];

SortSolutionBasis[SemiAnalytic`SemiAnalyticSolution[par_, basis_List]] :=
    SemiAnalytic`SemiAnalyticSolution[par, Sort[basis]];

Print["testing IsAllowedSemiAnalyticParameters ..."];

TestEquality[SemiAnalytic`IsAllowedSemiAnalyticParameter[T[Yu]], True];
TestEquality[SemiAnalytic`IsAllowedSemiAnalyticParameter[T[Yu][i1,i2]], True];
TestEquality[SemiAnalytic`IsAllowedSemiAnalyticParameter[T[Yu][1,2]], True];
TestEquality[SemiAnalytic`IsAllowedSemiAnalyticParameter[Yu], False];
TestEquality[SemiAnalytic`IsAllowedSemiAnalyticParameter[Yu[i1,i2]], False];
TestEquality[SemiAnalytic`IsAllowedSemiAnalyticParameter[vu], False];

Print["testing SetSemiAnalyticParameters ..."];

semianalyticPars = Sort[{ T[Yu], T[Yd], T[Ye], B[\[Mu]], mHu2, mHd2, mq2, mu2,
                          md2, ml2, me2, MassB, MassWB, MassG }];

SemiAnalytic`SetSemiAnalyticParameters[Parameters`GetModelParameters[]];

TestEquality[semianalyticPars, Sort[SemiAnalytic`GetSemiAnalyticParameters[]]];

Print["testing IsSemiAnalyticParameter ..."];

TestEquality[SemiAnalytic`IsSemiAnalyticParameter[T[Yu]], True];
TestEquality[SemiAnalytic`IsSemiAnalyticParameter[T[Yu][i1,i2]], True];
TestEquality[SemiAnalytic`IsSemiAnalyticParameter[T[Yu][3,2]], True];
TestEquality[SemiAnalytic`IsSemiAnalyticParameter[Yu], False];

Print["testing GetBoundaryValueSubstitutions ..."];

boundaryCondition = { {mHd2, m0^2}, {mHu2, m0^2}, FlexibleSUSY`FSSolveEWSBFor[{\[Mu], B[\[Mu]]}] };
expected = { Rule[mHd2, m0^2], Rule[mHu2, m0^2], Rule[B[\[Mu]], BMuEWSBSol] };
TestEquality[SemiAnalytic`Private`GetBoundaryValueSubstitutions[boundaryCondition], expected];

boundaryCondition = { {mHd2, m0^2}, {mHu2, m0^2}, FlexibleSUSY`FSSolveEWSBFor[{\[Mu], B[\[Mu]]}], {mHu2, \[Mu]^2} };
expected = { Rule[mHd2, m0^2], Rule[mHu2, \[Mu]^2], Rule[B[\[Mu]], BMuEWSBSol] };
TestEquality[SemiAnalytic`Private`GetBoundaryValueSubstitutions[boundaryCondition], expected];

boundaryCondition = { {mHd2, m0^2}, {mHu2, m0^2}, FlexibleSUSY`FSSolveEWSBFor[{\[Mu], B[\[Mu]]}],
                      {mHu2, \[Mu]^2}, FlexibleSUSY`FSMinimize[{T[Yu][1,1], mq2[1,1]}, {T[Yu][1,1]^2 + mq2[1,1]^2}] };
expected = { Rule[mHd2, m0^2], Rule[mHu2, \[Mu]^2], Rule[B[\[Mu]], BMuEWSBSol],
             Rule[T[Yu][1,1], TYu11MinSol], Rule[mq2[1,1], mq211MinSol] };
TestEquality[SemiAnalytic`Private`GetBoundaryValueSubstitutions[boundaryCondition], expected];

boundaryCondition = { {mHd2, m0^2}, {mHu2, m0^2}, {mq2, CConversion`UNITMATRIX[3] m0^2} };
expected = { Rule[mHd2, m0^2], Rule[mHu2, m0^2], Rule[mq2[1,1], m0^2], Rule[mq2[1,2], 0], Rule[mq2[1,3], 0],
             Rule[mq2[2,1], 0], Rule[mq2[2,2], m0^2], Rule[mq2[2,3], 0], Rule[mq2[3,1], 0], Rule[mq2[3,2], 0],
             Rule[mq2[3,3], m0^2] };
TestEquality[SemiAnalytic`Private`GetBoundaryValueSubstitutions[boundaryCondition], expected];

Print["testing ReplaceImplicitConstraints ..."];
boundaryCondition = { {mHd2, m0^2}, {mHu2, m0^2}, FlexibleSUSY`FSSolveEWSBFor[{\[Mu], B[\[Mu]]}],
                      {MassB, m12}, FlexibleSUSY`FSMinimize[{T[Yu][1,1], mq2[1,1]},{T[Yu][1,1]^2 + mq2[1,1]^2}],
                      {T[Yu][2,2], Azero} };
expected = { {mHd2, m0^2}, {mHu2, m0^2}, {B[\[Mu]], BMuEWSBSol}, {MassB, m12}, {T[Yu][1,1], TYu11MinSol},
             {mq2[1,1], mq211MinSol}, {T[Yu][2,2], Azero} };
TestEquality[SemiAnalytic`Private`ReplaceImplicitConstraints[boundaryCondition], expected];

Print["testing GetSolutionBasis ..."];

boundaryValues = { T[Yu][1,1] -> Yu[1,1] Azero, T[Yu][1,2] -> Yu[1,2] Azero,
                   T[Yu][1,3] -> Yu[1,3] Azero, T[Yu][2,1] -> Yu[2,1] Azero,
                   T[Yu][2,2] -> Yu[2,2] Azero, T[Yu][2,3] -> Yu[2,3] Azero,
                   T[Yu][3,1] -> Yu[3,1] Azero, T[Yu][3,2] -> Yu[3,2] Azero,
                   T[Yu][3,3] -> Yu[3,3] Azero, T[Yd][1,1] -> Yd[1,1] Azero,
                   T[Yd][1,2] -> Yd[1,2] Azero, T[Yd][1,3] -> Yd[1,3] Azero,
                   T[Yd][2,1] -> Yd[2,1] Azero, T[Yd][2,2] -> Yd[2,2] Azero,
                   T[Yd][2,3] -> Yd[2,3] Azero, T[Yd][3,1] -> Yd[3,1] Azero,
                   T[Yd][3,2] -> Yd[3,2] Azero, T[Yd][3,3] -> Yd[3,3] Azero,
                   T[Ye][1,1] -> Ye[1,1] Azero, T[Ye][1,2] -> Ye[1,2] Azero,
                   T[Ye][1,3] -> Ye[1,3] Azero, T[Ye][2,1] -> Ye[2,1] Azero,
                   T[Ye][2,2] -> Ye[2,2] Azero, T[Ye][2,3] -> Ye[2,3] Azero,
                   T[Ye][3,1] -> Ye[3,1] Azero, T[Ye][3,2] -> Ye[3,2] Azero,
                   T[Ye][3,3] -> Ye[3,3] Azero, MassB -> m12, MassWB -> m12,
                   MassG -> m12 };
dimOnePars = { T[Yu], T[Yd], T[Ye], MassB, MassWB, MassG };
expected = Sort[{m12, Azero}];
TestEquality[Sort[SemiAnalytic`Private`GetSolutionBasis[dimOnePars, boundaryValues]], expected];

boundaryValues = { mHd2 -> m0^2 + Azero^2 + m12^2 + Azero m12,
                   mHu2 -> m0^2 + Azero^2 + m12^2 + Azero m12 };
dimTwoPars = { mHd2, mHu2 };
expected = Sort[{m0^2, Azero^2, m12^2, Azero m12}];
TestEquality[Sort[SemiAnalytic`Private`GetSolutionBasis[dimTwoPars, boundaryValues]], expected];

Print["testing GetSemiAnalyticSolutions ..."];

SemiAnalytic`SetSemiAnalyticParameters[{T[Yu], T[Yd], T[Ye], MassB, MassWB, MassG}];
boundaryCondition = {{T[Yu], Azero*Yu}, {T[Yd], Azero*Yd}, {T[Ye], Azero*Ye},
                     {MassB, m12}, {MassWB, m12}, {MassG, m12}};
expected = Sort[SortSolutionBasis[#]& /@ {SemiAnalytic`SemiAnalyticSolution[MassB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassWB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassG, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yu], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yd], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Ye], {Azero, m12}]}];
TestEquality[Sort[(SortSolutionBasis[#])& /@ #]& @ SemiAnalytic`GetSemiAnalyticSolutions[boundaryCondition], expected];

boundaryCondition = {{T[Yu], T[Yd]}, {T[Yd], Azero*Yd}, {T[Ye], Azero*Ye},
                     {MassB, m12}, {MassWB, m12}, {MassG, m12}};
expected = Sort[SortSolutionBasis[#]& /@ {SemiAnalytic`SemiAnalyticSolution[MassB, {m12, Azero, T[Yd][0,0], T[Yd][0,1], T[Yd][0,2],
                                                                                    T[Yd][1,0], T[Yd][1,1], T[Yd][1,2],
                                                                                    T[Yd][2,0], T[Yd][2,1], T[Yd][2,2]}],
                                          SemiAnalytic`SemiAnalyticSolution[MassWB, {m12, Azero, T[Yd][0,0], T[Yd][0,1], T[Yd][0,2],
                                                                                     T[Yd][1,0], T[Yd][1,1], T[Yd][1,2],
                                                                                     T[Yd][2,0], T[Yd][2,1], T[Yd][2,2]}],
                                          SemiAnalytic`SemiAnalyticSolution[MassG, {m12, Azero, T[Yd][0,0], T[Yd][0,1], T[Yd][0,2],
                                                                                    T[Yd][1,0], T[Yd][1,1], T[Yd][1,2],
                                                                                    T[Yd][2,0], T[Yd][2,1], T[Yd][2,2]}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yu], {m12, Azero, T[Yd][0,0], T[Yd][0,1], T[Yd][0,2],
                                                                                    T[Yd][1,0], T[Yd][1,1], T[Yd][1,2],
                                                                                    T[Yd][2,0], T[Yd][2,1], T[Yd][2,2]}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yd], {m12, Azero, T[Yd][0,0], T[Yd][0,1], T[Yd][0,2],
                                                                                    T[Yd][1,0], T[Yd][1,1], T[Yd][1,2],
                                                                                    T[Yd][2,0], T[Yd][2,1], T[Yd][2,2]}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Ye], {m12, Azero, T[Yd][0,0], T[Yd][0,1], T[Yd][0,2],
                                                                                    T[Yd][1,0], T[Yd][1,1], T[Yd][1,2],
                                                                                    T[Yd][2,0], T[Yd][2,1], T[Yd][2,2]}]}];
TestEquality[Sort[(SortSolutionBasis[#])& /@ #]& @ SemiAnalytic`GetSemiAnalyticSolutions[boundaryCondition], expected];

SemiAnalytic`SetSemiAnalyticParameters[{T[Yu], T[Yd], T[Ye], MassB, MassWB, MassG, mHd2, mHu2}];
boundaryCondition = {{T[Yu], Azero*Yu}, {T[Yd], Azero*Yd}, {T[Ye], Azero*Ye},
                     {MassB, m12}, {MassWB, m12}, {MassG, m12}, {mHd2, m0^2}, {mHu2, m0^2}};
expected = Sort[SortSolutionBasis[#]& /@ {SemiAnalytic`SemiAnalyticSolution[MassB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassWB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassG, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yu], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yd], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Ye], {Azero, m12}],
                                          SemiAnalytic`SemiAnalyticSolution[mHd2, {m0^2, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[mHu2, {m0^2, m12^2, m12 Azero, Azero^2}]}];
TestEquality[Sort[(SortSolutionBasis[#])& /@ #]& @ SemiAnalytic`GetSemiAnalyticSolutions[boundaryCondition], expected];

SemiAnalytic`SetSemiAnalyticParameters[{T[Yu], T[Yd], T[Ye], MassB, MassWB, MassG, mHd2, mHu2, B[\[Mu]]}];
boundaryCondition = {{T[Yu], Azero*Yu}, {T[Yd], Azero*Yd}, {T[Ye], Azero*Ye},
                     {MassB, m12}, {MassWB, m12}, {MassG, m12}, {mHd2, m0^2}, {mHu2, m0^2},
                     {B[\[Mu]], \[Mu]*B0}};
expected = Sort[SortSolutionBasis[#]& /@ {SemiAnalytic`SemiAnalyticSolution[MassB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassWB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassG, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yu], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yd], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Ye], {Azero, m12}],
                                          SemiAnalytic`SemiAnalyticSolution[mHd2, {m0^2, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[mHu2, {m0^2, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[B[\[Mu]], {\[Mu] B0, \[Mu] Azero, \[Mu] m12}]}];
TestEquality[Sort[(SortSolutionBasis[#])& /@ #]& @ SemiAnalytic`GetSemiAnalyticSolutions[boundaryCondition], expected];

boundaryCondition = {{T[Yu], Azero*Yu}, {T[Yd], Azero*Yd}, {T[Ye], Azero*Ye},
                     {MassB, m12}, {MassWB, m12}, {MassG, m12}, FlexibleSUSY`FSSolveEWSBFor[{mHd2, mHu2}],
                     {B[\[Mu]], \[Mu]*B0}};
expected = Sort[SortSolutionBasis[#]& /@ {SemiAnalytic`SemiAnalyticSolution[MassB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassWB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassG, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yu], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yd], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Ye], {Azero, m12}],
                                          SemiAnalytic`SemiAnalyticSolution[mHd2, {mHu2EWSBSol, mHd2EWSBSol, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[mHu2, {mHu2EWSBSol, mHd2EWSBSol, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[B[\[Mu]], {\[Mu] B0, \[Mu] Azero, \[Mu] m12}]}];
TestEquality[Sort[(SortSolutionBasis[#])& /@ #]& @ SemiAnalytic`GetSemiAnalyticSolutions[boundaryCondition], expected];

boundaryCondition = {{T[Yu], Azero*Yu}, {T[Yd], Azero*Yd}, {T[Ye], Azero*Ye},
                     {MassB, m12}, {MassWB, m12}, {MassG, m12}, FlexibleSUSY`FSFindRoot[{mHd2, mHu2}, {mHd2, mHu2}],
                     {B[\[Mu]], \[Mu]*B0}};
expected = Sort[SortSolutionBasis[#]& /@ {SemiAnalytic`SemiAnalyticSolution[MassB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassWB, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[MassG, {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yu], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Yd], {m12, Azero}],
                                          SemiAnalytic`SemiAnalyticSolution[T[Ye], {Azero, m12}],
                                          SemiAnalytic`SemiAnalyticSolution[mHd2, {mHu2RootSol, mHd2RootSol, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[mHu2, {mHu2RootSol, mHd2RootSol, m12^2, m12 Azero, Azero^2}],
                                          SemiAnalytic`SemiAnalyticSolution[B[\[Mu]], {\[Mu] B0, \[Mu] Azero, \[Mu] m12}]}];
TestEquality[Sort[(SortSolutionBasis[#])& /@ #]& @ SemiAnalytic`GetSemiAnalyticSolutions[boundaryCondition], expected];

DeleteDirectory[SARAH`SARAH[OutputDirectory], DeleteContents -> True];

PrintTestSummary[];
