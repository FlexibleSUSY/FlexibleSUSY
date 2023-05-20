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

Get["utils/load-FlexibleSUSY.m"];

Needs["TestSuite`", "TestSuite.m"];
Needs["CXXDiagrams`", "CXXDiagrams.m"];

Start["E6SSM"];

{susyBetaFunctions, susyBreakingBetaFunctions} =
  FlexibleSUSY`ReadSARAHBetaFunctions[];
allParameters =
  FlexibleSUSY`SetupModelParameters[susyBetaFunctions,
   susyBreakingBetaFunctions];

{massMatrices, Lat$massMatrices} =
  FlexibleSUSY`SetupMassMatrices[allParameters];
FlexibleSUSY`SetupOutputParameters[massMatrices];

(* ChiralVertex with left = -right *)
TestEquality[
   Vertex[{FSI, FSI, VZ}][[2, 1]] + Vertex[{FSI, FSI, VZ}][[3, 1]],
   0
];

TestEquality[
   StringContainsQ[
      CXXDiagrams`CreateVertices[{{FSI, FSI, VZ}}][[1,2]],
      "const auto gN = MODELPARAMETER(gN);"
   ], True
];

PrintTestSummary[];

