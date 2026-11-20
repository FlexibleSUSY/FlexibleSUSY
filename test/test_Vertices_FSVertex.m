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

(* Tests Vertices`FSVertex[], Vertices`IsNonZeroVertex[] and
   Vertices`LoadVerticesIfNecessary[] against SARAH`Vertex[] in the SM. *)

Needs["TestSuite`", "TestSuite.m"];
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{workingDirectory, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["SM"];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

(* compares the coefficients of two vertices, which may differ in form *)
SameVertexQ[v1_, v2_] :=
    v1[[1]] === v2[[1]] &&
    v1[[2 ;;, 2]] === v2[[2 ;;, 2]] &&
    Simplify[(v1[[2 ;;, 1]] - v2[[2 ;;, 1]]) /.
                SARAH`sum[i_, a_, b_, x_] :> Sum[x, {i, a, b}]] === Table[0, {Length[v1] - 1}];

StripVertexFields[v_] := Vertices`StripFieldIndices /@ v[[1]];

Print["testing FSVertex[] without the vertex lists ..."];

TestEquality[Head[SARAH`VertexList3] === Symbol || SARAH`VertexList3 === {}, True];
TestEquality[SameVertexQ[Vertices`FSVertex[{bar[Fu], Fu, hh}], SARAH`Vertex[{bar[Fu], Fu, hh}]], True];

Print["writing and reading SARAH's vertex lists ..."];

SA`CurrentStates = SARAH`EWSB;
SARAH`InitVertexCalculation[SARAH`EWSB, False];
SARAH`partDefinition = SARAH`ParticleDefinitions[SARAH`EWSB];
SARAH`Particles[SARAH`Current] = SARAH`Particles[SARAH`EWSB];
SARAH`MakeVertexList[SARAH`EWSB, SARAH`VerticesForLoops -> True];

SARAH`VertexList3 = {};
SARAH`VertexList4 = {};
Vertices`LoadVerticesIfNecessary[];

TestEquality[ListQ[SARAH`VertexList3] && Length[SARAH`VertexList3] > 0, True];
TestEquality[ListQ[SARAH`VertexList4] && Length[SARAH`VertexList4] > 0, True];

Print["testing FSVertex[] without dependences ..."];

(* the entries of the vertex lists are SARAH`Vertex[] without dependences *)
TestEquality[
    Select[Join[SARAH`VertexList3, SARAH`VertexList4],
           !SameVertexQ[Vertices`FSVertex[StripVertexFields[#]], SARAH`Vertex[StripVertexFields[#]]] &],
    {}];

(* the entries of the vertex lists are in SARAH's order *)
TestEquality[
    Select[StripVertexFields /@ Join[SARAH`VertexList3, SARAH`VertexList4],
           Vertices`SortFieldsInCp[#] =!= # &],
    {}];

Print["testing FSVertex[] with dependences ..."];

TestEquality[
    Select[Join[SARAH`VertexList3, SARAH`VertexList4],
           Vertices`FSVertex[StripVertexFields[#], True] =!=
              SARAH`Vertex[StripVertexFields[#], SARAH`UseDependences -> True] &],
    {}];

Print["testing IsNonZeroVertex[] ..."];

(* vertices that vanish only with the dependences (through
   g1 Cos[ThetaW] = g2 Sin[ThetaW]); they are not in the vertex lists *)
zeroWithDependences = {
    {hh, VP, VP}, {hh, VP, VZ}, {Ah, hh, VP}, {bar[Fv], Fv, VP},
    {hh, hh, VP, VP}, {hh, hh, VP, VZ}, {Ah, Ah, VP, VP}, {Ah, Ah, VP, VZ}
};

TestEquality[Vertices`IsNonZeroVertex[#], True]& /@ zeroWithDependences;
TestEquality[Vertices`IsNonZeroVertex[#, False], True]& /@ zeroWithDependences;
TestEquality[Vertices`IsNonZeroVertex[#, True], False]& /@ zeroWithDependences;
TestEquality[MemberQ[Vertices`FSVertex[#, True][[2 ;;, 1]], Except[0]], False]& /@ zeroWithDependences;

(* in the vertex lists, but zero with dependences: SARAH doesn't
   simplify ghost vertices before checking them, so it counts as
   non-zero *)
TestEquality[Vertices`IsNonZeroVertex[{bar[gZ], gP, hh}, True], True];
TestEquality[MemberQ[Vertices`FSVertex[{bar[gZ], gP, hh}, True][[2 ;;, 1]], Except[0]], False];

(* fields in any order *)
TestEquality[Vertices`IsNonZeroVertex[{VP, hh, VP}, True], False];
TestEquality[Vertices`IsNonZeroVertex[{VP, hh, VP}, False], True];
TestEquality[Vertices`IsNonZeroVertex[{hh, Fu, bar[Fu]}, True], True];
TestEquality[Vertices`IsNonZeroVertex[{VWp, VP, conj[VWp]}, True], True];
TestEquality[Vertices`IsNonZeroVertex[{VP, conj[VWp], VWp, VZ}, True], True];

(* charge-violating vertices are zero without asking SARAH *)
numberOfFSVertexResults[] := Length[DownValues[Vertices`FSVertex]];
Module[{n = numberOfFSVertexResults[]},
   TestEquality[Vertices`IsNonZeroVertex[{VZ, bar[Fd], Fu}, #], False]& /@ {False, True};
   TestEquality[Vertices`IsNonZeroVertex[{VWp, Hp, hh}, #], False]& /@ {False, True};
   TestEquality[Vertices`IsNonZeroVertex[{VP, Hp, Hp, hh}, #], False]& /@ {False, True};
   TestEquality[numberOfFSVertexResults[], n];
];

(* the charge check never rejects a vertex of SARAH's lists, also with
   conjugated fields *)
TestEquality[
    Select[StripVertexFields /@ Join[SARAH`VertexList3, SARAH`VertexList4],
           Vertices`Private`ViolatesElectricChargeQ],
    {}];
TestEquality[Vertices`IsNonZeroVertex[{conj[VWp], Hp, hh}, True], True];

(* vertices that vanish without dependences *)
TestEquality[Vertices`IsNonZeroVertex[{hh, VG, VG}, #], False]& /@ {False, True};
TestEquality[Vertices`IsNonZeroVertex[{bar[Fe], Fu, hh}, #], False]& /@ {False, True};

PrintTestSummary[];
