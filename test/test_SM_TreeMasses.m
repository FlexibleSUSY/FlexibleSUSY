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
Start["SM"];
Needs["TestSuite`", "TestSuite.m"];

Print["testing IsSMParticleElementwise[] ..."];

TestEquality[IsSMParticleElementwise[Fu], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fd], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fe], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fv], {True, True, True}];
TestEquality[IsSMParticleElementwise[VP], {True}];
TestEquality[IsSMParticleElementwise[VZ], {True}];
TestEquality[IsSMParticleElementwise[VWp], {True}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> True], {True}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> False], {True}];
TestEquality[IsSMParticleElementwise[Ah], {True}];
TestEquality[IsSMParticleElementwise[Hp], {True}];

PrintTestSummary[];
