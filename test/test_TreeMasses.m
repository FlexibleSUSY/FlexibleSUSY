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
Needs["TreeMasses`", "TreeMasses.m"];

Print["testing IsSymmetric[] ..."];

TestEquality[TreeMasses`Private`IsSymmetric[{}], True];
TestEquality[TreeMasses`Private`IsSymmetric[{a}], False];
TestEquality[TreeMasses`Private`IsSymmetric[{{a}}], True];
TestEquality[TreeMasses`Private`IsSymmetric[{{a,b},{b,a}}], True];
TestEquality[TreeMasses`Private`IsSymmetric[{{a,b},{c,a}}], False];

PrintTestSummary[];
