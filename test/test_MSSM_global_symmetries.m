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

(* Tests Vertices`GetExactGlobalSymmetries[] and
   Vertices`ViolatesGlobalSymmetryQ[] in the MSSM, where R-parity is
   exact. *)

Needs["TestSuite`", "TestSuite.m"];
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{workingDirectory, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["MSSM"];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

SA`CurrentStates = SARAH`EWSB;
SARAH`InitVertexCalculation[SARAH`EWSB, False];
SARAH`partDefinition = SARAH`ParticleDefinitions[SARAH`EWSB];
SARAH`Particles[SARAH`Current] = SARAH`Particles[SARAH`EWSB];
SARAH`MakeVertexList[SARAH`EWSB, SARAH`VerticesForLoops -> True];
SARAH`VertexList3 = {};
SARAH`VertexList4 = {};

Print["testing GetExactGlobalSymmetries[] ..."];

TestEquality[Vertices`GetExactGlobalSymmetries[], {RParity}];

(* no vertex of SARAH's lists is rejected *)
TestEquality[
    Select[Join[SARAH`VertexList3, SARAH`VertexList4],
           Vertices`ViolatesGlobalSymmetryQ[#[[1]]]&],
    {}];

Print["testing ViolatesGlobalSymmetryQ[] ..."];

(* R-parity violating vertices; their electric charges add up to 0 *)
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{Su, Sv, conj[Su]}], True];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{Ah, hh, Sv}], True];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{bar[Fe], Fe, Sv}], True];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{VZ, Chi, Fv}], True];
(* with indices *)
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{Ah[{gt1}], hh[{gt2}], Sv[{gt3}]}], True];

(* R-parity conserving vertices *)
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{hh, Su, conj[Su]}], False];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{bar[Chi], Fe, conj[Se]}], False];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{bar[Cha], Cha, VZ}], False];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{Glu, Glu, VG}], False];

(* a decay into an R-parity violating final state, squark -> quark neutrino *)
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{conj[Su], Fu, Fv}], True];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{conj[Su], Fu, Chi}], False];

Print["testing IsNonZeroVertex[] ..."];

(* violating vertices are zero without asking SARAH *)
numberOfFSVertexResults[] := Length[DownValues[Vertices`FSVertex]];
Module[{n = numberOfFSVertexResults[]},
   TestEquality[Vertices`IsNonZeroVertex[{Ah, hh, Sv}, #], False]& /@ {False, True};
   TestEquality[Vertices`IsNonZeroVertex[{Su, Sv, conj[Su]}, #], False]& /@ {False, True};
   TestEquality[numberOfFSVertexResults[], n];
];
TestEquality[Vertices`IsNonZeroVertex[{hh, Su, conj[Su]}, True], True];

PrintTestSummary[];
