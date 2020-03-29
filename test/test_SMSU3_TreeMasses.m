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

AppendTo[$Path, FileNameJoin[{Directory[], "meta"}]];

Needs["SARAH`"];
Needs["FlexibleSUSY`", "FlexibleSUSY.m"];
Needs["TestSuite`", "TestSuite.m"];
Needs["TreeMasses`", "TreeMasses.m"];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{Directory[], "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};

Start["SMSU3"];

Print["testing IsMassless[] ..."];

TestEquality[TreeMasses`IsMassless[gG], True];
TestEquality[TreeMasses`IsMassless[VG], True];
TestEquality[TreeMasses`IsMassless[gP], True];
TestEquality[TreeMasses`IsMassless[VP], True];
TestEquality[TreeMasses`IsMassless[gZ], False];
TestEquality[TreeMasses`IsMassless[VZ], False];
TestEquality[TreeMasses`IsMassless[gWp], False];
TestEquality[TreeMasses`IsMassless[VWp], False];
TestEquality[TreeMasses`IsMassless[gCo2], True];
TestEquality[TreeMasses`IsMassless[VCo2], True];

PrintTestSummary[];
