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
Start["munuSSM"];
Needs["TestSuite`", "TestSuite.m"];

Print["testing IsSMParticleElementwise[] ..."];

TestEquality[IsSMParticleElementwise[Fu], {True, True, True}];
TestEquality[IsSMParticleElementwise[Fd], {True, True, True}];
TestEquality[IsSMParticleElementwise[Chi], {True, True, True, False, False, False, False, False}];
TestEquality[IsSMParticleElementwise[Cha], {True, True, True, False, False}];
TestEquality[IsSMParticleElementwise[Glu], {False}];
TestEquality[IsSMParticleElementwise[VP], {True}];
TestEquality[IsSMParticleElementwise[VZ], {True}];
TestEquality[IsSMParticleElementwise[VWm], {True}];
TestEquality[IsSMParticleElementwise[Su], {False, False, False, False, False, False}];
TestEquality[IsSMParticleElementwise[Sd], {False, False, False, False, False, False}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> True], {True, False, False, False, False, False}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> False], {False, False, False, False, False, False}];
TestEquality[IsSMParticleElementwise[Ah], {True, False, False, False, False, False}];
TestEquality[IsSMParticleElementwise[Hpm], {True, False, False, False, False, False, False, False}];

PrintTestSummary[];
