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

(* Tests Vertices`GetExactGlobalSymmetries[] in the MRSSM.  The
   R-symmetry is used as broken: the model file contains the
   superpotential term mu Hu Hd, which violates it, so that SARAH
   generates the corresponding B-term; mu is set to zero later.  As
   the code generation can't assume that, vertices violating the
   R-symmetry must not be skipped. *)

Needs["TestSuite`", "TestSuite.m"];
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{workingDirectory, "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Off[Superpotential::ViolationGlobal];
Start["MRSSM"];
On[Superpotential::ViolationGlobal];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

TestEquality[SARAH`Global[[All, 2]], {RSymmetry}];

(* rejected by the superpotential term, before the vertex lists are needed *)
TestEquality[Vertices`Private`ModelTermsConserveGlobalSymmetryQ[RSymmetry], False];
TestEquality[Vertices`GetExactGlobalSymmetries[], {}];

(* nothing is skipped, although the charges are known *)
TestEquality[Vertices`Private`GlobalCharge[Rh, RSymmetry], 2];
TestEquality[Vertices`ViolatesGlobalSymmetryQ[{Rh, hh, hh}], False];

PrintTestSummary[];
