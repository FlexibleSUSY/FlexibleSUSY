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

Print["testing FSAntiField[] ..."];

TestEquality[FSAntiField[hh], hh];
TestEquality[FSAntiField[Ah], Ah];
TestEquality[FSAntiField[Hp], conj[Hp]];
TestEquality[FSAntiField[VP], VP];
TestEquality[FSAntiField[VZ], VZ];
TestEquality[FSAntiField[VWp], conj[VWp]];
TestEquality[FSAntiField[gP], conj[gP]];
TestEquality[FSAntiField[gZ], conj[gZ]];
TestEquality[FSAntiField[gWp], conj[gWp]];
TestEquality[FSAntiField[gWpC], conj[gWpC]];
TestEquality[FSAntiField[Fv], bar[Fv]];
TestEquality[FSAntiField[Fe], bar[Fe]];
TestEquality[FSAntiField[Fu], bar[Fu]];
TestEquality[FSAntiField[Fd], bar[Fd]];

Print["testing GetSMParticles[] ..."];

smParticles = {VG, gG, Hp, Fv, Ah, hh, VP, VZ, gP, gZ, gWp, gWpC, Fd, Fu, Fe, VWp};

TestEquality[Sort[TreeMasses`GetSMParticles[]], Sort[smParticles]];
TestEquality[Sort[TreeMasses`GetParticles[]]  , Sort[smParticles]];

Print["testing IsParticle[] ..."];

TestEquality[TreeMasses`IsParticle /@ smParticles, Array[True&, Length[smParticles]]];

conjSMParticles = {VG, conj[gG], conj[Hp], bar[Fv], Ah, hh, VP, VZ,
                   conj[gP], conj[gZ], conj[gWp], conj[gWpC], bar[Fd],
                   bar[Fu], bar[Fe], conj[VWp]};

TestEquality[Sort[FSAntiField /@ smParticles], Sort[conjSMParticles]];

TestEquality[TreeMasses`IsParticle /@ FSAntiField /@ smParticles,
             Array[True&, Length[conjSMParticles]]];

Print["testing IsScalar[] ..."];

TestEquality[IsScalar[hh], True];
TestEquality[IsScalar[Ah], True];
TestEquality[IsScalar[Hp], True];
TestEquality[IsScalar[FSAntiField[Hp]], True];
TestEquality[IsScalar[VP], False];
TestEquality[IsScalar[VZ], False];
TestEquality[IsScalar[VWp], False];
TestEquality[IsScalar[FSAntiField[VWp]], False];
TestEquality[IsScalar[gP], False];
TestEquality[IsScalar[gZ], False];
TestEquality[IsScalar[gWp], False];
TestEquality[IsScalar[gWpC], False];
TestEquality[IsScalar[FSAntiField[gP]], False];
TestEquality[IsScalar[FSAntiField[gZ]], False];
TestEquality[IsScalar[FSAntiField[gWp]], False];
TestEquality[IsScalar[FSAntiField[gWpC]], False];
TestEquality[IsScalar[Fv], False];
TestEquality[IsScalar[Fe], False];
TestEquality[IsScalar[Fu], False];
TestEquality[IsScalar[Fd], False];
TestEquality[IsScalar[FSAntiField[Fv]], False];
TestEquality[IsScalar[FSAntiField[Fe]], False];
TestEquality[IsScalar[FSAntiField[Fu]], False];
TestEquality[IsScalar[FSAntiField[Fd]], False];

Print["testing IsRealScalar[] ..."];

TestEquality[IsRealScalar[hh], True];
TestEquality[IsRealScalar[Ah], True];
TestEquality[IsRealScalar[Hp], False];
TestEquality[IsRealScalar[FSAntiField[Hp]], False];

Print["testing IsComplexScalar[] ..."];

TestEquality[IsComplexScalar[hh], False];
TestEquality[IsComplexScalar[Ah], False];
TestEquality[IsComplexScalar[Hp], True];
TestEquality[IsComplexScalar[FSAntiField[Hp]], True];

Print["testing IsFermion[] ..."];

TestEquality[IsFermion[hh], False];
TestEquality[IsFermion[Ah], False];
TestEquality[IsFermion[Hp], False];
TestEquality[IsFermion[FSAntiField[Hp]], False];
TestEquality[IsFermion[VP], False];
TestEquality[IsFermion[VZ], False];
TestEquality[IsFermion[VWp], False];
TestEquality[IsFermion[FSAntiField[VWp]], False];
TestEquality[IsFermion[gP], False];
TestEquality[IsFermion[gZ], False];
TestEquality[IsFermion[gWp], False];
TestEquality[IsFermion[gWpC], False];
TestEquality[IsFermion[FSAntiField[gP]], False];
TestEquality[IsFermion[FSAntiField[gZ]], False];
TestEquality[IsFermion[FSAntiField[gWp]], False];
TestEquality[IsFermion[FSAntiField[gWpC]], False];
TestEquality[IsFermion[Fv], True];
TestEquality[IsFermion[Fe], True];
TestEquality[IsFermion[Fu], True];
TestEquality[IsFermion[Fd], True];
TestEquality[IsFermion[FSAntiField[Fv]], True];
TestEquality[IsFermion[FSAntiField[Fe]], True];
TestEquality[IsFermion[FSAntiField[Fu]], True];
TestEquality[IsFermion[FSAntiField[Fd]], True];

Print["testing IsMajoranaFermion[] ..."];

TestEquality[IsMajoranaFermion[Fv], False];
TestEquality[IsMajoranaFermion[Fe], False];
TestEquality[IsMajoranaFermion[Fu], False];
TestEquality[IsMajoranaFermion[Fd], False];
TestEquality[IsMajoranaFermion[FSAntiField[Fv]], False];
TestEquality[IsMajoranaFermion[FSAntiField[Fe]], False];
TestEquality[IsMajoranaFermion[FSAntiField[Fu]], False];
TestEquality[IsMajoranaFermion[FSAntiField[Fd]], False];

Print["testing IsDiracFermion[] ..."];

TestEquality[IsDiracFermion[Fv], True];
TestEquality[IsDiracFermion[Fe], True];
TestEquality[IsDiracFermion[Fu], True];
TestEquality[IsDiracFermion[Fd], True];
TestEquality[IsDiracFermion[FSAntiField[Fv]], True];
TestEquality[IsDiracFermion[FSAntiField[Fe]], True];
TestEquality[IsDiracFermion[FSAntiField[Fu]], True];
TestEquality[IsDiracFermion[FSAntiField[Fd]], True];

Print["testing IsVector[] ..."];

TestEquality[IsVector[hh], False];
TestEquality[IsVector[Ah], False];
TestEquality[IsVector[Hp], False];
TestEquality[IsVector[FSAntiField[Hp]], False];
TestEquality[IsVector[VP], True];
TestEquality[IsVector[VZ], True];
TestEquality[IsVector[VWp], True];
TestEquality[IsVector[FSAntiField[VWp]], True];
TestEquality[IsVector[gP], False];
TestEquality[IsVector[gZ], False];
TestEquality[IsVector[gWp], False];
TestEquality[IsVector[gWpC], False];
TestEquality[IsVector[FSAntiField[gP]], False];
TestEquality[IsVector[FSAntiField[gZ]], False];
TestEquality[IsVector[FSAntiField[gWp]], False];
TestEquality[IsVector[FSAntiField[gWpC]], False];
TestEquality[IsVector[Fv], False];
TestEquality[IsVector[Fe], False];
TestEquality[IsVector[Fu], False];
TestEquality[IsVector[Fd], False];
TestEquality[IsVector[FSAntiField[Fv]], False];
TestEquality[IsVector[FSAntiField[Fe]], False];
TestEquality[IsVector[FSAntiField[Fu]], False];
TestEquality[IsVector[FSAntiField[Fd]], False];

Print["testing IsRealVector[] ..."];

TestEquality[IsRealVector[VP], True];
TestEquality[IsRealVector[VZ], True];
TestEquality[IsRealVector[VWp], False];
TestEquality[IsRealVector[FSAntiField[VWp]], False];

Print["testing IsComplexVector[] ..."];

TestEquality[IsComplexVector[VP], False];
TestEquality[IsComplexVector[VZ], False];
TestEquality[IsComplexVector[VWp], True];
TestEquality[IsComplexVector[FSAntiField[VWp]], True];

Print["testing IsGhost[] ..."];

TestEquality[IsGhost[hh], False];
TestEquality[IsGhost[Ah], False];
TestEquality[IsGhost[Hp], False];
TestEquality[IsGhost[FSAntiField[Hp]], False];
TestEquality[IsGhost[VP], False];
TestEquality[IsGhost[VZ], False];
TestEquality[IsGhost[VWp], False];
TestEquality[IsGhost[FSAntiField[VWp]], False];
TestEquality[IsGhost[gP], True];
TestEquality[IsGhost[gZ], True];
TestEquality[IsGhost[gWp], True];
TestEquality[IsGhost[gWpC], True];
TestEquality[IsGhost[FSAntiField[gP]], True];
TestEquality[IsGhost[FSAntiField[gZ]], True];
TestEquality[IsGhost[FSAntiField[gWp]], True];
TestEquality[IsGhost[FSAntiField[gWpC]], True];
TestEquality[IsGhost[Fv], False];
TestEquality[IsGhost[Fe], False];
TestEquality[IsGhost[Fu], False];
TestEquality[IsGhost[Fd], False];
TestEquality[IsGhost[FSAntiField[Fv]], False];
TestEquality[IsGhost[FSAntiField[Fe]], False];
TestEquality[IsGhost[FSAntiField[Fu]], False];
TestEquality[IsGhost[FSAntiField[Fd]], False];

Print["testing IsGoldstone[] ..."];

TestEquality[IsGoldstone[hh], False];
TestEquality[IsGoldstone[Ah], True];
TestEquality[IsGoldstone[Hp], True];
TestEquality[IsGoldstone[FSAntiField[Hp]], True];
TestEquality[IsGoldstone[VP], False];
TestEquality[IsGoldstone[VZ], False];
TestEquality[IsGoldstone[VWp], False];
TestEquality[IsGoldstone[FSAntiField[VWp]], False];
TestEquality[IsGoldstone[gP], False];
TestEquality[IsGoldstone[gZ], False];
TestEquality[IsGoldstone[gWp], False];
TestEquality[IsGoldstone[gWpC], False];
TestEquality[IsGoldstone[FSAntiField[gP]], False];
TestEquality[IsGoldstone[FSAntiField[gZ]], False];
TestEquality[IsGoldstone[FSAntiField[gWp]], False];
TestEquality[IsGoldstone[FSAntiField[gWpC]], False];
TestEquality[IsGoldstone[Fv], False];
TestEquality[IsGoldstone[Fe], False];
TestEquality[IsGoldstone[Fu], False];
TestEquality[IsGoldstone[Fd], False];
TestEquality[IsGoldstone[FSAntiField[Fv]], False];
TestEquality[IsGoldstone[FSAntiField[Fe]], False];
TestEquality[IsGoldstone[FSAntiField[Fu]], False];
TestEquality[IsGoldstone[FSAntiField[Fd]], False];

Print["testing IsAuxiliary[] ..."];

TestEquality[IsAuxiliary[hh], False];
TestEquality[IsAuxiliary[Ah], False];
TestEquality[IsAuxiliary[Hp], False];
TestEquality[IsAuxiliary[FSAntiField[Hp]], False];
TestEquality[IsAuxiliary[VP], False];
TestEquality[IsAuxiliary[VZ], False];
TestEquality[IsAuxiliary[VWp], False];
TestEquality[IsAuxiliary[FSAntiField[VWp]], False];
TestEquality[IsAuxiliary[gP], False];
TestEquality[IsAuxiliary[gZ], False];
TestEquality[IsAuxiliary[gWp], False];
TestEquality[IsAuxiliary[gWpC], False];
TestEquality[IsAuxiliary[FSAntiField[gP]], False];
TestEquality[IsAuxiliary[FSAntiField[gZ]], False];
TestEquality[IsAuxiliary[FSAntiField[gWp]], False];
TestEquality[IsAuxiliary[FSAntiField[gWpC]], False];
TestEquality[IsAuxiliary[Fv], False];
TestEquality[IsAuxiliary[Fe], False];
TestEquality[IsAuxiliary[Fu], False];
TestEquality[IsAuxiliary[Fd], False];
TestEquality[IsAuxiliary[FSAntiField[Fv]], False];
TestEquality[IsAuxiliary[FSAntiField[Fe]], False];
TestEquality[IsAuxiliary[FSAntiField[Fu]], False];
TestEquality[IsAuxiliary[FSAntiField[Fd]], False];

Print["testing IsSMParticle[] ..."];

TestEquality[TreeMasses`IsSMParticle /@ TreeMasses`GetSMParticles[],
             Array[True&, Length[TreeMasses`GetSMParticles[]]]];

TestEquality[TreeMasses`IsSMParticle /@ FSAntiField /@ TreeMasses`GetSMParticles[],
             Array[True&, Length[TreeMasses`GetSMParticles[]]]];

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

Print["testing IsSMGoldstone[] ..."];

TestEquality[IsSMGoldstone[hh], False];
TestEquality[IsSMGoldstone[Ah], True];
TestEquality[IsSMGoldstone[Hp], True];
TestEquality[IsSMGoldstone[FSAntiField[Hp]], True];
TestEquality[IsSMGoldstone[VP], False];
TestEquality[IsSMGoldstone[VZ], False];
TestEquality[IsSMGoldstone[VWp], False];
TestEquality[IsSMGoldstone[FSAntiField[VWp]], False];
TestEquality[IsSMGoldstone[gP], False];
TestEquality[IsSMGoldstone[gZ], False];
TestEquality[IsSMGoldstone[gWp], False];
TestEquality[IsSMGoldstone[gWpC], False];
TestEquality[IsSMGoldstone[FSAntiField[gP]], False];
TestEquality[IsSMGoldstone[FSAntiField[gZ]], False];
TestEquality[IsSMGoldstone[FSAntiField[gWp]], False];
TestEquality[IsSMGoldstone[FSAntiField[gWpC]], False];
TestEquality[IsSMGoldstone[Fv], False];
TestEquality[IsSMGoldstone[Fe], False];
TestEquality[IsSMGoldstone[Fu], False];
TestEquality[IsSMGoldstone[Fd], False];
TestEquality[IsSMGoldstone[FSAntiField[Fv]], False];
TestEquality[IsSMGoldstone[FSAntiField[Fe]], False];
TestEquality[IsSMGoldstone[FSAntiField[Fu]], False];
TestEquality[IsSMGoldstone[FSAntiField[Fd]], False];

Print["testing IsSMHiggs[] ..."];

TestEquality[IsSMHiggs[hh], True];
TestEquality[IsSMHiggs[Ah], False];
TestEquality[IsSMHiggs[Hp], False];
TestEquality[IsSMHiggs[FSAntiField[Hp]], False];
TestEquality[IsSMHiggs[VP], False];
TestEquality[IsSMHiggs[VZ], False];
TestEquality[IsSMHiggs[VWp], False];
TestEquality[IsSMHiggs[FSAntiField[VWp]], False];
TestEquality[IsSMHiggs[gP], False];
TestEquality[IsSMHiggs[gZ], False];
TestEquality[IsSMHiggs[gWp], False];
TestEquality[IsSMHiggs[gWpC], False];
TestEquality[IsSMHiggs[FSAntiField[gP]], False];
TestEquality[IsSMHiggs[FSAntiField[gZ]], False];
TestEquality[IsSMHiggs[FSAntiField[gWp]], False];
TestEquality[IsSMHiggs[FSAntiField[gWpC]], False];
TestEquality[IsSMHiggs[Fv], False];
TestEquality[IsSMHiggs[Fe], False];
TestEquality[IsSMHiggs[Fu], False];
TestEquality[IsSMHiggs[Fd], False];
TestEquality[IsSMHiggs[FSAntiField[Fv]], False];
TestEquality[IsSMHiggs[FSAntiField[Fe]], False];
TestEquality[IsSMHiggs[FSAntiField[Fu]], False];
TestEquality[IsSMHiggs[FSAntiField[Fd]], False];

Print["testing IsMassless[] ..."];

TestEquality[IsMassless[hh], False];
TestEquality[IsMassless[Ah], False];
TestEquality[IsMassless[Hp], False];
TestEquality[IsMassless[FSAntiField[Hp]], False];
TestEquality[IsMassless[VP], True];
TestEquality[IsMassless[VZ], False];
TestEquality[IsMassless[VWp], False];
TestEquality[IsMassless[FSAntiField[VWp]], False];
TestEquality[IsMassless[gP], True];
TestEquality[IsMassless[gZ], False];
TestEquality[IsMassless[gWp], False];
TestEquality[IsMassless[gWpC], False];
TestEquality[IsMassless[FSAntiField[gP]], True];
TestEquality[IsMassless[FSAntiField[gZ]], False];
TestEquality[IsMassless[FSAntiField[gWp]], False];
TestEquality[IsMassless[FSAntiField[gWpC]], False];
TestEquality[IsMassless[Fv], True];
TestEquality[IsMassless[Fe], False];
TestEquality[IsMassless[Fu], False];
TestEquality[IsMassless[Fd], False];
TestEquality[IsMassless[FSAntiField[Fv]], True];
TestEquality[IsMassless[FSAntiField[Fe]], False];
TestEquality[IsMassless[FSAntiField[Fu]], False];
TestEquality[IsMassless[FSAntiField[Fd]], False];

Print["testing IsQuark[] ..."];

TestEquality[IsQuark[hh], False];
TestEquality[IsQuark[Ah], False];
TestEquality[IsQuark[Hp], False];
TestEquality[IsQuark[FSAntiField[Hp]], False];
TestEquality[IsQuark[VP], False];
TestEquality[IsQuark[VZ], False];
TestEquality[IsQuark[VWp], False];
TestEquality[IsQuark[FSAntiField[VWp]], False];
TestEquality[IsQuark[gP], False];
TestEquality[IsQuark[gZ], False];
TestEquality[IsQuark[gWp], False];
TestEquality[IsQuark[gWpC], False];
TestEquality[IsQuark[FSAntiField[gP]], False];
TestEquality[IsQuark[FSAntiField[gZ]], False];
TestEquality[IsQuark[FSAntiField[gWp]], False];
TestEquality[IsQuark[FSAntiField[gWpC]], False];
TestEquality[IsQuark[Fv], False];
TestEquality[IsQuark[Fe], False];
TestEquality[IsQuark[Fu], True];
TestEquality[IsQuark[Fd], True];
TestEquality[IsQuark[FSAntiField[Fv]], False];
TestEquality[IsQuark[FSAntiField[Fe]], False];
TestEquality[IsQuark[FSAntiField[Fu]], True];
TestEquality[IsQuark[FSAntiField[Fd]], True];

Print["testing IsLepton[] ..."];

TestEquality[IsLepton[hh], False];
TestEquality[IsLepton[Ah], False];
TestEquality[IsLepton[Hp], False];
TestEquality[IsLepton[FSAntiField[Hp]], False];
TestEquality[IsLepton[VP], False];
TestEquality[IsLepton[VZ], False];
TestEquality[IsLepton[VWp], False];
TestEquality[IsLepton[FSAntiField[VWp]], False];
TestEquality[IsLepton[gP], False];
TestEquality[IsLepton[gZ], False];
TestEquality[IsLepton[gWp], False];
TestEquality[IsLepton[gWpC], False];
TestEquality[IsLepton[FSAntiField[gP]], False];
TestEquality[IsLepton[FSAntiField[gZ]], False];
TestEquality[IsLepton[FSAntiField[gWp]], False];
TestEquality[IsLepton[FSAntiField[gWpC]], False];
TestEquality[IsLepton[Fv], True];
TestEquality[IsLepton[Fe], True];
TestEquality[IsLepton[Fu], False];
TestEquality[IsLepton[Fd], False];
TestEquality[IsLepton[FSAntiField[Fv]], True];
TestEquality[IsLepton[FSAntiField[Fe]], True];
TestEquality[IsLepton[FSAntiField[Fu]], False];
TestEquality[IsLepton[FSAntiField[Fd]], False];

Print["testing IsSMChargedLepton[] ..."];

TestEquality[IsSMChargedLepton[hh], False];
TestEquality[IsSMChargedLepton[Ah], False];
TestEquality[IsSMChargedLepton[Hp], False];
TestEquality[IsSMChargedLepton[FSAntiField[Hp]], False];
TestEquality[IsSMChargedLepton[VP], False];
TestEquality[IsSMChargedLepton[VZ], False];
TestEquality[IsSMChargedLepton[VWp], False];
TestEquality[IsSMChargedLepton[FSAntiField[VWp]], False];
TestEquality[IsSMChargedLepton[gP], False];
TestEquality[IsSMChargedLepton[gZ], False];
TestEquality[IsSMChargedLepton[gWp], False];
TestEquality[IsSMChargedLepton[gWpC], False];
TestEquality[IsSMChargedLepton[FSAntiField[gP]], False];
TestEquality[IsSMChargedLepton[FSAntiField[gZ]], False];
TestEquality[IsSMChargedLepton[FSAntiField[gWp]], False];
TestEquality[IsSMChargedLepton[FSAntiField[gWpC]], False];
TestEquality[IsSMChargedLepton[Fv], False];
TestEquality[IsSMChargedLepton[Fe], True];
TestEquality[IsSMChargedLepton[Fu], False];
TestEquality[IsSMChargedLepton[Fd], False];
TestEquality[IsSMChargedLepton[FSAntiField[Fv]], False];
TestEquality[IsSMChargedLepton[FSAntiField[Fe]], True];
TestEquality[IsSMChargedLepton[FSAntiField[Fu]], False];
TestEquality[IsSMChargedLepton[FSAntiField[Fd]], False];

Print["testing IsSMNeutralLepton[] ..."];

TestEquality[IsSMNeutralLepton[hh], False];
TestEquality[IsSMNeutralLepton[Ah], False];
TestEquality[IsSMNeutralLepton[Hp], False];
TestEquality[IsSMNeutralLepton[FSAntiField[Hp]], False];
TestEquality[IsSMNeutralLepton[VP], False];
TestEquality[IsSMNeutralLepton[VZ], False];
TestEquality[IsSMNeutralLepton[VWp], False];
TestEquality[IsSMNeutralLepton[FSAntiField[VWp]], False];
TestEquality[IsSMNeutralLepton[gP], False];
TestEquality[IsSMNeutralLepton[gZ], False];
TestEquality[IsSMNeutralLepton[gWp], False];
TestEquality[IsSMNeutralLepton[gWpC], False];
TestEquality[IsSMNeutralLepton[FSAntiField[gP]], False];
TestEquality[IsSMNeutralLepton[FSAntiField[gZ]], False];
TestEquality[IsSMNeutralLepton[FSAntiField[gWp]], False];
TestEquality[IsSMNeutralLepton[FSAntiField[gWpC]], False];
TestEquality[IsSMNeutralLepton[Fv], True];
TestEquality[IsSMNeutralLepton[Fe], False];
TestEquality[IsSMNeutralLepton[Fu], False];
TestEquality[IsSMNeutralLepton[Fd], False];
TestEquality[IsSMNeutralLepton[FSAntiField[Fv]], True];
TestEquality[IsSMNeutralLepton[FSAntiField[Fe]], False];
TestEquality[IsSMNeutralLepton[FSAntiField[Fu]], False];
TestEquality[IsSMNeutralLepton[FSAntiField[Fd]], False];

Print["testing IsSMQuark[] ..."];

TestEquality[IsSMQuark[hh], False];
TestEquality[IsSMQuark[Ah], False];
TestEquality[IsSMQuark[Hp], False];
TestEquality[IsSMQuark[FSAntiField[Hp]], False];
TestEquality[IsSMQuark[VP], False];
TestEquality[IsSMQuark[VZ], False];
TestEquality[IsSMQuark[VWp], False];
TestEquality[IsSMQuark[FSAntiField[VWp]], False];
TestEquality[IsSMQuark[gP], False];
TestEquality[IsSMQuark[gZ], False];
TestEquality[IsSMQuark[gWp], False];
TestEquality[IsSMQuark[gWpC], False];
TestEquality[IsSMQuark[FSAntiField[gP]], False];
TestEquality[IsSMQuark[FSAntiField[gZ]], False];
TestEquality[IsSMQuark[FSAntiField[gWp]], False];
TestEquality[IsSMQuark[FSAntiField[gWpC]], False];
TestEquality[IsSMQuark[Fv], False];
TestEquality[IsSMQuark[Fe], False];
TestEquality[IsSMQuark[Fu], True];
TestEquality[IsSMQuark[Fd], True];
TestEquality[IsSMQuark[FSAntiField[Fv]], False];
TestEquality[IsSMQuark[FSAntiField[Fe]], False];
TestEquality[IsSMQuark[FSAntiField[Fu]], True];
TestEquality[IsSMQuark[FSAntiField[Fd]], True];

Print["testing IsSMUpQuark[] ..."];

TestEquality[IsSMUpQuark[hh], False];
TestEquality[IsSMUpQuark[Ah], False];
TestEquality[IsSMUpQuark[Hp], False];
TestEquality[IsSMUpQuark[FSAntiField[Hp]], False];
TestEquality[IsSMUpQuark[VP], False];
TestEquality[IsSMUpQuark[VZ], False];
TestEquality[IsSMUpQuark[VWp], False];
TestEquality[IsSMUpQuark[FSAntiField[VWp]], False];
TestEquality[IsSMUpQuark[gP], False];
TestEquality[IsSMUpQuark[gZ], False];
TestEquality[IsSMUpQuark[gWp], False];
TestEquality[IsSMUpQuark[gWpC], False];
TestEquality[IsSMUpQuark[FSAntiField[gP]], False];
TestEquality[IsSMUpQuark[FSAntiField[gZ]], False];
TestEquality[IsSMUpQuark[FSAntiField[gWp]], False];
TestEquality[IsSMUpQuark[FSAntiField[gWpC]], False];
TestEquality[IsSMUpQuark[Fv], False];
TestEquality[IsSMUpQuark[Fe], False];
TestEquality[IsSMUpQuark[Fu], True];
TestEquality[IsSMUpQuark[Fd], False];
TestEquality[IsSMUpQuark[FSAntiField[Fv]], False];
TestEquality[IsSMUpQuark[FSAntiField[Fe]], False];
TestEquality[IsSMUpQuark[FSAntiField[Fu]], True];
TestEquality[IsSMUpQuark[FSAntiField[Fd]], False];

Print["testing IsSMDownQuark[] ..."];

TestEquality[IsSMDownQuark[hh], False];
TestEquality[IsSMDownQuark[Ah], False];
TestEquality[IsSMDownQuark[Hp], False];
TestEquality[IsSMDownQuark[FSAntiField[Hp]], False];
TestEquality[IsSMDownQuark[VP], False];
TestEquality[IsSMDownQuark[VZ], False];
TestEquality[IsSMDownQuark[VWp], False];
TestEquality[IsSMDownQuark[FSAntiField[VWp]], False];
TestEquality[IsSMDownQuark[gP], False];
TestEquality[IsSMDownQuark[gZ], False];
TestEquality[IsSMDownQuark[gWp], False];
TestEquality[IsSMDownQuark[gWpC], False];
TestEquality[IsSMDownQuark[FSAntiField[gP]], False];
TestEquality[IsSMDownQuark[FSAntiField[gZ]], False];
TestEquality[IsSMDownQuark[FSAntiField[gWp]], False];
TestEquality[IsSMDownQuark[FSAntiField[gWpC]], False];
TestEquality[IsSMDownQuark[Fv], False];
TestEquality[IsSMDownQuark[Fe], False];
TestEquality[IsSMDownQuark[Fu], False];
TestEquality[IsSMDownQuark[Fd], True];
TestEquality[IsSMDownQuark[FSAntiField[Fv]], False];
TestEquality[IsSMDownQuark[FSAntiField[Fe]], False];
TestEquality[IsSMDownQuark[FSAntiField[Fu]], False];
TestEquality[IsSMDownQuark[FSAntiField[Fd]], True];

PrintTestSummary[];
