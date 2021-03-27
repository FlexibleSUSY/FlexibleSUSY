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

Get["models/CMSSM/CMSSM_librarylink.m"];

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

handle = FSCMSSMOpenHandle[
    fsSettings -> settings,
    fsSMParameters -> smInputs,
    fsModelParameters -> {
        m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
];

specML = CMSSM /. FSCMSSMCalculateSpectrum[handle];
obsML  = CMSSM /. FSCMSSMCalculateObservables[handle];
probML = FSCMSSMGetProblems[handle];
warnML = FSCMSSMGetWarnings[handle];

FSCMSSMCloseHandle[handle];

inputFile = "test/test_CMSSM_librarylink.in.spc";
outputFile = "test/test_CMSSM_librarylink.out.spc";
cmd = "models/CMSSM/run_CMSSM.x --slha-input-file=" <> inputFile <>
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
    {mq2, {3, 3}, MSQ2},
    {ml2, {3, 3}, MSL2},
    {mu2, {3, 3}, MSU2},
    {md2, {3, 3}, MSD2},
    {me2, {3, 3}, MSE2},
    {mHu2, {0}, {MSOFT, 22}},
    {mHd2, {0}, {MSOFT, 21}},
    {vu, {0}, {HMIX, 103}},
    {vd, {0}, {HMIX, 102}},
    {Mu, {0}, {HMIX, 1}},
    {BMu, {0}, {HMIX, 101}},
    {CpHPP1, {0}, {EFFHIGGSCOUPLINGS, 25, 22, 22}},
    {CpHPP2, {0}, {EFFHIGGSCOUPLINGS, 35, 22, 22}},
    {CpHGG1, {0}, {EFFHIGGSCOUPLINGS, 25, 21, 21}},
    {CpHGG2, {0}, {EFFHIGGSCOUPLINGS, 35, 21, 21}},
    {CpAPP, {0}, {EFFHIGGSCOUPLINGS, 36, 22, 22}},
    {CpAGG, {0}, {EFFHIGGSCOUPLINGS, 36, 21, 21}},
    {aMuon, {0}, {FlexibleSUSYLowEnergy, 21}}
};

slhaData = ReadSLHAString[slhaStr, parameters];

delta = 1*^-8;

TestCloseRel[g1 * Sqrt[5/3] /. slhaData, g1 /. specML, delta];
TestCloseRel[g2 /. slhaData, g2 /. specML, delta];
TestCloseRel[g3 /. slhaData, g3 /. specML, delta];
TestCloseRel[Yu /. slhaData, Yu /. specML, delta];
TestCloseRel[Yd /. slhaData, Yd /. specML, delta];
TestCloseRel[Ye /. slhaData, Ye /. specML, delta];
TestCloseRel[Tu /. slhaData, T[Yu] /. specML, delta];
TestCloseRel[Td /. slhaData, T[Yd] /. specML, delta];
TestCloseRel[Te /. slhaData, T[Ye] /. specML, delta];
TestCloseRel[mq2 /. slhaData, mq2 /. specML, delta];
TestCloseRel[ml2 /. slhaData, ml2 /. specML, delta];
TestCloseRel[mu2 /. slhaData, mu2 /. specML, delta];
TestCloseRel[md2 /. slhaData, md2 /. specML, delta];
TestCloseRel[me2 /. slhaData, me2 /. specML, delta];
TestCloseRel[mHu2 /. slhaData, mHu2 /. specML, delta];
TestCloseRel[mHd2 /. slhaData, mHd2 /. specML, delta];
TestCloseRel[vu /. slhaData, vu /. specML, delta];
TestCloseRel[vd /. slhaData, vd /. specML, delta];
TestCloseRel[Mu /. slhaData, \[Mu] /. specML, delta];
TestCloseRel[BMu /. slhaData, B[\[Mu]] /. specML, delta];

delta = 1*^-6;

TestCloseRel[aMuon  /. slhaData, FlexibleSUSYObservable`aMuon /. obsML, delta];

TestEquality[CMSSM /. probML, {}];
TestEquality[CMSSM /. warnML, {}];

Print["Check re-calculation of spectrum yields the same"];

CalcMh[TB_] :=
    Module[{spec, handle},
           handle = FSCMSSMOpenHandle[
               fsSettings -> settings,
               fsSMParameters -> smInputs,
               fsModelParameters -> {
                   m0 -> 125, m12 -> 500, TanBeta -> TB, SignMu -> 1, Azero -> 0 }
           ];
           spec = CMSSM /. FSCMSSMCalculateSpectrum[handle];
           FSCMSSMCloseHandle[handle];
           If[spec === $Failed, 0,
              (Pole[M[hh]] /. spec)[[1]]
             ]
          ];

Mhh1 = CalcMh[45];
Mhh2 = CalcMh[45];

TestEquality[Mhh1, Mhh2];

Print["Testing Set[] function"];

CalcMh[handle_] :=
    Module[{spec},
           spec = CMSSM /. FSCMSSMCalculateSpectrum[handle];
           If[spec === $Failed, 0,
              (Pole[M[hh]] /. spec)[[1]]
             ]
          ];

handle = FSCMSSMOpenHandle[
    fsSettings -> settings,
    fsSMParameters -> smInputs,
    fsModelParameters -> {
        m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
];

Mhh1 = CalcMh[handle];

FSCMSSMSet[handle, TanBeta -> 10];

Mhh2 = CalcMh[handle];

FSCMSSMSet[handle, fsModelParameters -> { TanBeta -> 10 }];

Mhh3 = CalcMh[handle];

TestEquality[Mhh1, Mhh2];
TestEquality[Mhh1, Mhh3];

FSCMSSMSet[handle, TanBeta -> 20];

Mhh4 = CalcMh[handle];

TestNonEquality[Mhh1, Mhh4];

PrintTestSummary[];
