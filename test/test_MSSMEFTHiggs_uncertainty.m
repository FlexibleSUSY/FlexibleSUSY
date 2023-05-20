Needs["TestSuite`", "TestSuite.m"];

Get["models/MSSMEFTHiggs/MSSMEFTHiggs_librarylink.m"];
Get["model_files/MSSMEFTHiggs/MSSMEFTHiggs_uncertainty_estimate.m"];

Mtpole = 173.34;

settings = {
    precisionGoal -> 1.*^-5,
    maxIterations -> 100,
    calculateStandardModelMasses -> 1,
    calculateBSMMasses -> 0,
    betaFunctionLoopOrder -> 3,
    poleMassLoopOrder -> 3,
    ewsbLoopOrder -> 3,
    thresholdCorrectionsLoopOrder -> 3,
    thresholdCorrections -> 123111321
};

smpars = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.166378700*^-5,  (* SMINPUTS[2] *)
    alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.18,           (* SMINPUTS[5] *)
    Mt -> Mtpole,           (* SMINPUTS[6] *)
    Mtau -> 1.777,          (* SMINPUTS[7] *)
    Mv3 -> 0,               (* SMINPUTS[8] *)
    MW -> 80.385,           (* SMINPUTS[9] *)
    Me -> 0.000510998902,   (* SMINPUTS[11] *)
    Mv1 -> 0,               (* SMINPUTS[12] *)
    Mm -> 0.1056583715,     (* SMINPUTS[13] *)
    Mv2 -> 0,               (* SMINPUTS[14] *)
    md2GeV -> 0.00475,      (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
    ms2GeV -> 0.104,        (* SMINPUTS[23] *)
    mcmc -> 1.27,           (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

MSSMEFTHiggsCalcMhDMh[MS_, TB_, Xtt_] :=
    CalcMSSMEFTHiggsDMh[
        fsSettings -> settings,
        fsSMParameters -> smpars,
        fsModelParameters -> {
            MSUSY   -> MS,
            M1Input -> MS,
            M2Input -> MS,
            M3Input -> MS,
            MuInput -> MS,
            mAInput -> MS,
            TanBeta -> TB,
            mq2Input -> MS^2 IdentityMatrix[3],
            mu2Input -> MS^2 IdentityMatrix[3],
            md2Input -> MS^2 IdentityMatrix[3],
            ml2Input -> MS^2 IdentityMatrix[3],
            me2Input -> MS^2 IdentityMatrix[3],
            AuInput -> {{MS/TB, 0    , 0},
                        {0    , MS/TB, 0},
                        {0    , 0    , MS/TB + Xtt MS}},
            AdInput -> MS TB IdentityMatrix[3],
            AeInput -> MS TB IdentityMatrix[3]
        }
   ];

MSSMEFTHiggsCalcMh[MS_, TB_, Xtt_] :=
    Module[{handle, spec},
           handle = FSMSSMEFTHiggsOpenHandle[
               fsSettings -> settings,
               fsSMParameters -> smpars,
               fsModelParameters -> {
                   MSUSY   -> MS,
                   M1Input -> MS,
                   M2Input -> MS,
                   M3Input -> MS,
                   MuInput -> MS,
                   mAInput -> MS,
                   TanBeta -> TB,
                   mq2Input -> MS^2 IdentityMatrix[3],
                   mu2Input -> MS^2 IdentityMatrix[3],
                   md2Input -> MS^2 IdentityMatrix[3],
                   ml2Input -> MS^2 IdentityMatrix[3],
                   me2Input -> MS^2 IdentityMatrix[3],
                   AuInput -> {{MS/TB, 0    , 0},
                               {0    , MS/TB, 0},
                               {0    , 0    , MS/TB + Xtt MS}},
                   AdInput -> MS TB IdentityMatrix[3],
                   AeInput -> MS TB IdentityMatrix[3]
               }
           ];
           spec = FSMSSMEFTHiggsCalculateSpectrum[handle];
           FSMSSMEFTHiggsCloseHandle[handle];
           If[spec === $Failed, $Failed,
              Pole[M[hh]] /. (MSSMEFTHiggs /. spec)]
          ];

Xtt = 0;
TBX = 5;
MSX = 5000;

Mh1 = MSSMEFTHiggsCalcMh[MSX, TBX, Xtt][[1]];
Mh2 = MSSMEFTHiggsCalcMhDMh[MSX, TBX, Xtt][[1]];

TestEquality[Mh1, Mh2];

PrintTestSummary[];
