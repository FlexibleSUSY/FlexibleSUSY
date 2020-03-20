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
Needs["TreeMasses`", "TreeMasses.m"];

Start["MRSSM"];

Print[""];

Print["Testing getters for particle collection..."];

TestEquality[
   TreeMasses`GetSusyParticles[],
   {Glu, SRdp, SRum, sigmaO, phiO, Sd, Sv, Su, Se, hh, Ah, Rh, Hpm, Chi, Cha1, Cha2}];

TestEquality[
   TreeMasses`GetColoredParticles[],
   {VG, gG, Glu, sigmaO, phiO, Sd, Su, Fd, Fu}
];

TestEquality[
   TreeMasses`GetVectorBosons[], {VG, VP, VZ, VWm}
];

Print["Testing getters for individual particles..."];

TestEquality[TreeMasses`GetPhoton[], VP];
TestEquality[TreeMasses`GetGluon[], VG];
TestEquality[TreeMasses`GetZBoson[], VZ];
TestEquality[TreeMasses`GetWBoson[], VWm];
TestEquality[TreeMasses`GetHiggsBoson[], hh];
TestEquality[TreeMasses`GetChargedHiggsBoson[], Hpm];
TestEquality[TreeMasses`GetPseudoscalarHiggsBoson[], Ah];

Print["Testing particle properties..."];

TestEquality[TreeMasses`IsSMParticle[#], True]& /@ {Fe, Fd, Fu, Fv, gG, VG, gP, VP, VZ, gZ, VWm, gWm, gWmC};
TestEquality[TreeMasses`IsSMParticle[#], False]& /@ {Glu, SRdp, SRum, sigmaO, phiO, Sd, Sv, Su, Se, hh, Ah, Rh, Hpm, Chi, Cha1, Cha2};

TestEquality[TreeMasses`IsMassless[#], True]& /@ {Fv, gG, VG, gP, VP};
TestEquality[TreeMasses`IsMassless[#], False]& /@ {
   gZ, VZ, gWm, gWmC, VWm, SRdp, SRum, sigmaO, phiO, Sd, Sv, Su, Se, hh, Ah, Rh, Hpm, Fe, Fd, Fu, Cha1, Cha2, Chi, Glu};

TestEquality[TreeMasses`IsScalar[#], True]& /@ {SRdp, SRum, sigmaO, phiO, Sd, Sv, Su, Se, hh, Ah, Rh, Hpm};
TestEquality[TreeMasses`IsFermion[#], True]& /@ {Fe, Fd, Fu, Fv, Cha1, Cha2, Chi, Glu};
TestEquality[TreeMasses`IsVector[#], True]& /@ {VG, VP, VZ, VWm};
TestEquality[TreeMasses`IsGhost[#], True]& /@ {gWm, gWmC, gP, gZ, gG};
TestEquality[TreeMasses`IsGoldstone[#], True]& /@ {Ah[{1}], Hpm[{1}]};

TestEquality[TreeMasses`IsElectricallyCharged[#], True]& /@ {SRdp, SRum, Sd, Su, Se, Fe, Fd, Fu, Hpm, Cha1, Cha2, VWm, gWm, gWmC};
TestEquality[TreeMasses`ColorChargedQ[#], True]& /@ {Fd, Fu, VG, gG, Su, Sd, Glu, sigmaO, phiO};

TestEquality[SA`CPViolationHiggsSector, False];

Print[""];

PrintTestSummary[];
