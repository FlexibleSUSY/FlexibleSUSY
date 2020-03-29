Get["models/MSSMEFTHiggs/MSSMEFTHiggs_librarylink.m"];
Get["model_files/MSSMEFTHiggs/MSSMEFTHiggs_uncertainty_estimate.m"];

Get["models/HSSUSY/HSSUSY_librarylink.m"];
Get["model_files/HSSUSY/HSSUSY_uncertainty_estimate.m"];

Get["models/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_librarylink.m"];
Get["model_files/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_uncertainty_estimate.m"];

settings = {
    precisionGoal -> 1.*^-5,
    maxIterations -> 1000,
    poleMassLoopOrder -> 4,
    ewsbLoopOrder -> 4,
    thresholdCorrectionsLoopOrder -> 4,
    thresholdCorrections -> 124111421
};

smpars = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.166378700*^-5,  (* SMINPUTS[2] *)
    alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.18,           (* SMINPUTS[5] *)
    Mt -> 173.34,           (* SMINPUTS[6] *)
    Mtau -> 1.77699,        (* SMINPUTS[7] *)
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

MSSMEFTHiggsCalcMh[MS_, TB_, Xtt_] :=
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

HSSUSYCalcMh[MS_, TB_, Xtt_] :=
    CalcHSSUSYDMh[
        fsSettings -> settings,
        fsSMParameters -> smpars,
        fsModelParameters -> {
            TanBeta -> TB,
            MEWSB -> 173.34,
            MSUSY -> MS,
            M1Input -> MS,
            M2Input -> MS,
            M3Input -> MS,
            MuInput -> MS,
            mAInput -> MS,
            AtInput -> (Xtt + 1/TB) MS,
            msq2 -> MS^2 IdentityMatrix[3],
            msu2 -> MS^2 IdentityMatrix[3],
            msd2 -> MS^2 IdentityMatrix[3],
            msl2 -> MS^2 IdentityMatrix[3],
            mse2 -> MS^2 IdentityMatrix[3],
            LambdaLoopOrder -> 2,
            TwoLoopAtAs -> 1,
            TwoLoopAbAs -> 1,
            TwoLoopAtAb -> 1,
            TwoLoopAtauAtau -> 1,
            TwoLoopAtAt -> 1
        }
   ];

NUHMSSMNoFVCalcMh[MS_, TB_, Xtt_] :=
    CalcNUHMSSMNoFVHimalayaDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            maxIterations -> 1000,
            poleMassLoopOrder -> 2,
            ewsbLoopOrder -> 2,
            thresholdCorrectionsLoopOrder -> 2,
            thresholdCorrections -> 121111121
        },
        fsSMParameters -> smpars,
        fsModelParameters -> {
            TanBeta -> TB,
            Qin -> MS,
            M1 -> MS,
            M2 -> MS,
            M3 -> MS,
            AtIN -> MS/TB + Xtt MS,
            AbIN -> MS TB,
            AtauIN -> MS TB,
            AcIN -> MS/TB,
            AsIN -> MS TB,
            AmuonIN -> MS TB,
            AuIN -> MS/TB,
            AdIN -> MS TB,
            AeIN -> MS TB,
            MuIN -> MS,
            mA2IN -> MS^2,
            ml11IN -> MS,
            ml22IN -> MS,
            ml33IN -> MS,
            me11IN -> MS,
            me22IN -> MS,
            me33IN -> MS,
            mq11IN -> MS,
            mq22IN -> MS,
            mq33IN -> MS,
            mu11IN -> MS,
            mu22IN -> MS,
            mu33IN -> MS,
            md11IN -> MS,
            md22IN -> MS,
            md33IN -> MS
        }
   ];

NUHMSSMNoFVHimalayaCalcMh[MS_, TB_, Xtt_] :=
    CalcNUHMSSMNoFVHimalayaDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            maxIterations -> 1000,
            poleMassLoopOrder -> 3,
            ewsbLoopOrder -> 3,
            thresholdCorrectionsLoopOrder -> 2,
            thresholdCorrections -> 122111121
        },
        fsSMParameters -> smpars,
        fsModelParameters -> {
            TanBeta -> TB,
            Qin -> MS,
            M1 -> MS,
            M2 -> MS,
            M3 -> MS,
            AtIN -> MS/TB + Xtt MS,
            AbIN -> MS TB,
            AtauIN -> MS TB,
            AcIN -> MS/TB,
            AsIN -> MS TB,
            AmuonIN -> MS TB,
            AuIN -> MS/TB,
            AdIN -> MS TB,
            AeIN -> MS TB,
            MuIN -> MS,
            mA2IN -> MS^2,
            ml11IN -> MS,
            ml22IN -> MS,
            ml33IN -> MS,
            me11IN -> MS,
            me22IN -> MS,
            me33IN -> MS,
            mq11IN -> MS,
            mq22IN -> MS,
            mq33IN -> MS,
            mu11IN -> MS,
            mu22IN -> MS,
            mu33IN -> MS,
            md11IN -> MS,
            md22IN -> MS,
            md33IN -> MS
        }
   ];

(*** Xt = 0 ***)

TB  = 5;
Xtt = 0;

data = Map[
    { N[#],
      Sequence @@ NUHMSSMNoFVCalcMh[#, TB, Xtt],
      Sequence @@ NUHMSSMNoFVHimalayaCalcMh[#, TB, Xtt],
      Sequence @@ MSSMEFTHiggsCalcMh[#, TB, Xtt],
      Sequence @@ HSSUSYCalcMh[#, TB, Xtt]
    }&,
    LogRange[300, 5 10^4, 60]
];

Export["test/test_Mh_uncertainties_TB-5_Xt-0.dat", data];

(*** Xt = -2 ***)

Xtt = -2;

data = Map[
    { N[#],
      Sequence @@ NUHMSSMNoFVCalcMh[#, TB, Xtt],
      Sequence @@ NUHMSSMNoFVHimalayaCalcMh[#, TB, Xtt],
      Sequence @@ MSSMEFTHiggsCalcMh[#, TB, Xtt],
      Sequence @@ HSSUSYCalcMh[#, TB, Xtt]
    }&,
    LogRange[300, 5 10^4, 60]
];

Export["test/test_Mh_uncertainties_TB-5_Xt--2.dat", data];
