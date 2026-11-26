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

(* Tests the selection of 2-body final states in Decays` in the MSSM. *)

Needs["TestSuite`", "TestSuite.m"];
Needs["Decays`", "Decays.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{workingDirectory, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["MSSM"];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

Print["testing IsScalarToScalarAndMasslessVector[] ..."];

isForbidden = Decays`Private`IsScalarToScalarAndMasslessVector;

(* scalar -> scalar + photon or gluon vanishes at every order *)
TestEquality[isForbidden[hh, {Ah, VP}], True];
TestEquality[isForbidden[Ah, {VP, hh}], True];
TestEquality[isForbidden[Sv, {Sv, VP}], True];
TestEquality[isForbidden[Su, {Su, VG}], True];
TestEquality[isForbidden[Su, {VP, Su}], True];
TestEquality[isForbidden[Hpm, {Hpm, VP}], True];
TestEquality[isForbidden[conj[Se], {conj[Se], VP}], True];

(* massive vectors have a longitudinal polarization *)
TestEquality[isForbidden[hh, {Ah, VZ}], False];
TestEquality[isForbidden[Hpm, {hh, VWm}], False];
TestEquality[isForbidden[Su, {Sd, conj[VWm]}], False];

(* other final states with a massless vector *)
TestEquality[isForbidden[hh, {VP, VP}], False];
TestEquality[isForbidden[hh, {VP, VZ}], False];
TestEquality[isForbidden[Hpm, {VWm, VP}], False];
TestEquality[isForbidden[hh, {VG, VG}], False];

(* the initial particle isn't a scalar *)
TestEquality[isForbidden[Chi, {Chi, VP}], False];
TestEquality[isForbidden[VZ, {hh, VP}], False];

(* no vector in the final state *)
TestEquality[isForbidden[hh, {Ah, Ah}], False];
TestEquality[isForbidden[Su, {Fu, Chi}], False];

(* IsCPForbiddenScalarDecay[] checks the tree-level vertex, which needs
   SARAH's vertex lists *)
SA`CurrentStates = SARAH`EWSB;
SARAH`InitVertexCalculation[SARAH`EWSB, False];
SARAH`partDefinition = SARAH`ParticleDefinitions[SARAH`EWSB];
SARAH`Particles[SARAH`Current] = SARAH`Particles[SARAH`EWSB];
SARAH`MakeVertexList[SARAH`EWSB, SARAH`VerticesForLoops -> True];
SARAH`VertexList3 = {};
SARAH`VertexList4 = {};

Print["testing IsCPForbiddenScalarDecay[] ..."];

TestEquality[Decays`Private`GetDefiniteCP[hh], 1];
TestEquality[Decays`Private`GetDefiniteCP[Ah], -1];
TestEquality[Decays`Private`GetDefiniteCP[Hpm], 0];
TestEquality[Decays`Private`GetDefiniteCP[Sv], 0];

isCPForbidden = Decays`Private`IsCPForbiddenScalarDecay;

(* CP(S) != CP(S1) CP(S2) *)
TestEquality[isCPForbidden[hh, {hh, Ah}], True];
TestEquality[isCPForbidden[hh, {Ah, hh}], True];
TestEquality[isCPForbidden[Ah, {hh, hh}], True];
TestEquality[isCPForbidden[Ah, {Ah, Ah}], True];

(* CP(S) = CP(S1) for S -> S1 Z *)
TestEquality[isCPForbidden[hh, {hh, VZ}], True];
TestEquality[isCPForbidden[Ah, {VZ, Ah}], True];

(* CP-conserving decays *)
TestEquality[isCPForbidden[hh, {hh, hh}], False];
TestEquality[isCPForbidden[hh, {Ah, Ah}], False];
TestEquality[isCPForbidden[Ah, {hh, Ah}], False];
TestEquality[isCPForbidden[hh, {Ah, VZ}], False];
TestEquality[isCPForbidden[Ah, {hh, VZ}], False];

(* particles without a definite CP property, other final states *)
TestEquality[isCPForbidden[Hpm, {Hpm, hh}], False];
TestEquality[isCPForbidden[Hpm, {hh, VWm}], False];
TestEquality[isCPForbidden[hh, {Sv, conj[Sv]}], False];
TestEquality[isCPForbidden[Ah, {VZ, VZ}], False];
TestEquality[isCPForbidden[Ah, {VP, VP}], False];
TestEquality[isCPForbidden[Ah, {bar[Fd], Fd}], False];
TestEquality[isCPForbidden[Chi, {Chi, hh}], False];

PrintTestSummary[];
