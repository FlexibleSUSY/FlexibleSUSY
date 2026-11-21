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

(* Tests the global charges of a Z[3] symmetry in the NMSSM: the
   charges are phases, and an antiparticle has the inverse phase.  The
   Z[3] isn't exact: the charginos and neutralinos mix states of
   different charge. *)

Needs["TestSuite`", "TestSuite.m"];
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{workingDirectory, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["NMSSM"];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

charge[f_] := Vertices`Private`GlobalCharge[f, Z3];
conserves[fields_] := Vertices`Private`ConservesGlobalSymmetryQ[fields, Z3];

TestEquality[Vertices`Private`GlobalSymmetryGroup[Z3], "ZN"];
TestEquality[Vertices`Private`GlobalSymmetryGroup[RParity], "ZN"];
TestEquality[Simplify[charge[Sd]^3], 1];
TestEquality[Simplify[charge[conj[Sd]] charge[Sd]], 1];
TestEquality[Simplify[charge[Sd[{gt1, ct1}]] - charge[Sd]], 0];

TestEquality[conserves[{Sd, conj[Sd]}], True];
TestEquality[conserves[{Sd, Sd, Sd}], True];
TestEquality[conserves[{Sd, Sd}], False];
TestEquality[conserves[{Sd, Sd, conj[Sd]}], False];
(* unknown charges never reject a vertex *)
TestEquality[conserves[{Chi, Chi, Sd}], True];

TestEquality[Vertices`Private`ExactGlobalSymmetryQ[Z3], False];

PrintTestSummary[];
