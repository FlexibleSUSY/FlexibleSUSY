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
Needs["ReadSLHA`", "ReadSLHA.m"];

Get["models/CMSSMCPV/CMSSMCPV_librarylink.m"];

Print["Comparing SLHA output to LibraryLink output"];

settings = {
    precisionGoal -> 0.0001,
    maxIterations -> 0,
    calculateStandardModelMasses -> 0,
    poleMassLoopOrder -> 2,
    ewsbLoopOrder -> 2,
    betaFunctionLoopOrder -> 3,
    thresholdCorrectionsLoopOrder -> 2,
    higgs2loopCorrectionAtAs -> 1,
    higgs2loopCorrectionAbAs -> 1,
    higgs2loopCorrectionAtAt -> 1,
    higgs2loopCorrectionAtauAtau -> 1,
    forceOutput -> 0,
    topPoleQCDCorrections -> 1,
    betaZeroThreshold -> 1.*10^-11,
    forcePositiveMasses -> 0,
    poleMassScale -> 0.,
    parameterOutputScale -> 0.
};

smInputs = {
    alphaEmMZ -> 1./127.934,
    GF -> 0.0000116637,
    alphaSMZ -> 0.1176,
    MZ -> 91.1876,
    mbmb -> 4.2,
    Mt -> 173.3,
    Mtau -> 1.777,
    Mv3 -> 0,
    MW -> 80.404,
    Me -> 0.000510998902,
    Mv1 -> 0.,
    Mm -> 0.105658357,
    Mv2 -> 0.,
    md2GeV -> 0.00475,
    mu2GeV -> 0.0024,
    ms2GeV -> 0.104,
    mcmc -> 1.27
};

handle = FSCMSSMCPVOpenHandle[
    fsSettings -> settings,
    fsSMParameters -> smInputs,
    fsModelParameters -> {
        m0 -> 125, m12 -> 500, TanBeta -> 10, CosPhiMu -> 1, Azero -> 0,
        Imm12 -> 10, SinPhiMu -> 0, ImAzero -> 10, etaInput -> 0.1
    }
];

specML = CMSSMCPV /. FSCMSSMCPVCalculateSpectrum[handle];
obsML  = CMSSMCPV /. FSCMSSMCPVCalculateObservables[handle];
probML = FSCMSSMCPVGetProblems[handle];
warnML = FSCMSSMCPVGetWarnings[handle];

FSCMSSMCPVCloseHandle[handle];

inputFile = "test/test_CMSSMCPV_librarylink.in.spc";
outputFile = "test/test_CMSSMCPV_librarylink.out.spc";
cmd = "models/CMSSMCPV/run_CMSSMCPV.x --slha-input-file=" <> inputFile <>
    " --slha-output-file=" <> outputFile;

Run["rm -f " <> outputFile];
Run[cmd];
slhaStr = Import[outputFile, "String"];

parameters = {
    {g1, {0}   , {GAUGE, 1}},
    {g2, {0}   , {GAUGE, 2}},
    {g3, {0}   , {GAUGE, 3}},
    {Yu, {3, 3}, Yu},
    {Yd, {3, 3}, Yd},
    {Ye, {3, 3}, Ye},
    {Tu, {3, 3}, Tu},
    {Td, {3, 3}, Td},
    {Te, {3, 3}, Te},
    {ImTu, {3, 3}, ImTu},
    {ImTd, {3, 3}, ImTd},
    {ImTe, {3, 3}, ImTe},
    {mq2, {3, 3}, MSQ2},
    {ml2, {3, 3}, MSL2},
    {mu2, {3, 3}, MSU2},
    {md2, {3, 3}, MSD2},
    {me2, {3, 3}, MSE2},
    {Immq2, {3, 3}, ImMSQ2},
    {Imml2, {3, 3}, ImMSL2},
    {Immu2, {3, 3}, ImMSU2},
    {Immd2, {3, 3}, ImMSD2},
    {Imme2, {3, 3}, ImMSE2},
    {mHu2, {0}, {MSOFT, 22}},
    {mHd2, {0}, {MSOFT, 21}},
    {vu, {0}, {HMIX, 103}},
    {vd, {0}, {HMIX, 102}},
    {Mu, {0}, {HMIX, 1}},
    {BMu, {0}, {HMIX, 101}},
    {ImMu, {0}, {ImHMIX, 1}},
    {ImBMu, {0}, {ImHMIX, 101}},
    {CpHPP1, {0}, {EFFHIGGSCOUPLINGS, 25, 22, 22}},
    {CpHPP2, {0}, {EFFHIGGSCOUPLINGS, 35, 22, 22}},
    {CpHPP3, {0}, {EFFHIGGSCOUPLINGS, 36, 22, 22}},
    {CpHGG1, {0}, {EFFHIGGSCOUPLINGS, 25, 21, 21}},
    {CpHGG2, {0}, {EFFHIGGSCOUPLINGS, 35, 21, 21}},
    {CpHGG3, {0}, {EFFHIGGSCOUPLINGS, 36, 21, 21}},
    {aMuon, {0}, {FlexibleSUSYLowEnergy, 21}},
    {EDM[Fe[1]], {0}, {FlexibleSUSYLowEnergy, 23}},
    {EDM[Fe[2]], {0}, {FlexibleSUSYLowEnergy, 24}},
    {EDM[Fe[3]], {0}, {FlexibleSUSYLowEnergy, 25}}
};

slhaData = ReadSLHAString[slhaStr, parameters];

delta = 1*^-4;

TestCloseRel[g1 * Sqrt[5/3] /. slhaData, g1 /. specML, delta];
TestCloseRel[g2   /. slhaData, g2 /. specML, delta];
TestCloseRel[g3   /. slhaData, g3 /. specML, delta];
TestCloseRel[Abs[Yu] /. slhaData, Abs[Yu] /. specML, delta];
TestCloseRel[Abs[Yd] /. slhaData, Abs[Yd] /. specML, delta];
TestCloseRel[Abs[Ye] /. slhaData, Abs[Ye] /. specML, delta];
TestCloseRel[Abs[Tu + I ImTu] /. slhaData, Abs[T[Yu]] /. specML, delta];
TestCloseRel[Abs[Td + I ImTd] /. slhaData, Abs[T[Yd]] /. specML, delta];
TestCloseRel[Abs[Te + I ImTe] /. slhaData, Abs[T[Ye]] /. specML, delta];
TestCloseRel[Abs[mq2 + I Immq2] /. slhaData, Abs[mq2] /. specML, delta];
TestCloseRel[Abs[ml2 + I Imml2] /. slhaData, Abs[ml2] /. specML, delta];
TestCloseRel[Abs[mu2 + I Immu2] /. slhaData, Abs[mu2] /. specML, delta];
TestCloseRel[Abs[md2 + I Immd2] /. slhaData, Abs[md2] /. specML, delta];
TestCloseRel[Abs[me2 + I Imme2] /. slhaData, Abs[me2] /. specML, delta];
TestCloseRel[mHu2 /. slhaData, mHu2 /. specML, delta];
TestCloseRel[mHd2 /. slhaData, mHd2 /. specML, delta];
TestCloseRel[vu /. slhaData, vu /. specML, delta];
TestCloseRel[vd /. slhaData, vd /. specML, delta];
TestCloseRel[Mu + I ImMu /. slhaData, \[Mu] /. specML, delta];
TestCloseRel[BMu + I ImBMu /. slhaData, B[\[Mu]] /. specML, delta];

delta = 1*^-5;

TestCloseRel[aMuon  /. slhaData, FlexibleSUSYObservable`aMuon /. obsML, delta];
TestCloseRel[EDM[Fe[1]] /. slhaData, FlexibleSUSYObservable`EDM[Fe[1]] /. obsML, delta];
TestCloseRel[EDM[Fe[2]] /. slhaData, FlexibleSUSYObservable`EDM[Fe[2]] /. obsML, delta];
TestCloseRel[EDM[Fe[3]] /. slhaData, FlexibleSUSYObservable`EDM[Fe[3]] /. obsML, delta];

TestEquality[CMSSMCPV /. probML, {}];
TestEquality[CMSSMCPV /. warnML, {}];

PrintTestSummary[];
