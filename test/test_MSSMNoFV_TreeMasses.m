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
Start["MSSMNoFV"];
Needs["TestSuite`", "TestSuite.m"];

Print["testing IsSMParticleElementwise[] ..."];

TestEquality[IsSMParticleElementwise[Fu], {True}];
TestEquality[IsSMParticleElementwise[Fd], {True}];
TestEquality[IsSMParticleElementwise[Fs], {True}];
TestEquality[IsSMParticleElementwise[Fc], {True}];
TestEquality[IsSMParticleElementwise[Ft], {True}];
TestEquality[IsSMParticleElementwise[Fb], {True}];
TestEquality[IsSMParticleElementwise[Fe], {True}];
TestEquality[IsSMParticleElementwise[Fm], {True}];
TestEquality[IsSMParticleElementwise[Ftau], {True}];
TestEquality[IsSMParticleElementwise[Chi], {False, False, False, False}];
TestEquality[IsSMParticleElementwise[Cha], {False, False}];
TestEquality[IsSMParticleElementwise[Glu], {False}];
TestEquality[IsSMParticleElementwise[VP], {True}];
TestEquality[IsSMParticleElementwise[VZ], {True}];
TestEquality[IsSMParticleElementwise[VWm], {True}];
TestEquality[IsSMParticleElementwise[Su], {False, False}];
TestEquality[IsSMParticleElementwise[Sd], {False, False}];
TestEquality[IsSMParticleElementwise[Ss], {False, False}];
TestEquality[IsSMParticleElementwise[Sc], {False, False}];
TestEquality[IsSMParticleElementwise[St], {False, False}];
TestEquality[IsSMParticleElementwise[Sb], {False, False}];
TestEquality[IsSMParticleElementwise[Se], {False, False}];
TestEquality[IsSMParticleElementwise[Sm], {False, False}];
TestEquality[IsSMParticleElementwise[Stau], {False, False}];
TestEquality[IsSMParticleElementwise[SveL], {False}];
TestEquality[IsSMParticleElementwise[SvmL], {False}];
TestEquality[IsSMParticleElementwise[SvtL], {False}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> True], {True, False}];
TestEquality[IsSMParticleElementwise[hh, IncludeHiggs -> False], {False, False}];
TestEquality[IsSMParticleElementwise[Ah], {True, False}];
TestEquality[IsSMParticleElementwise[Hpm], {True, False}];

PrintTestSummary[];
