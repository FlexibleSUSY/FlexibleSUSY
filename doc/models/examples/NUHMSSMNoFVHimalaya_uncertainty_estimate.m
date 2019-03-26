Get["models/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_librarylink.m"];
Get["model_files/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_uncertainty_estimate.m"];

settings = {
    precisionGoal -> 1.*^-5,
    maxIterations -> 1000,
    betaFunctionLoopOrder -> 3,
    poleMassLoopOrder -> 3,
    ewsbLoopOrder -> 3,
    thresholdCorrectionsLoopOrder -> 2,
    thresholdCorrections -> 122111121
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

NUHMSSMNoFVHimalayaCalcMh[MS_, TB_, Xtt_] :=
    CalcNUHMSSMNoFVHimalayaDMh[
        fsSettings -> settings,
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

LinearRange[start_, stop_, steps_] :=
    Range[start, stop, (stop - start)/steps];

Xtt = -2;
TB  = 5;

data = ParallelMap[
    { N[#], Sequence @@ NUHMSSMNoFVHimalayaCalcMh[#, TB, Xtt] }&,
    LinearRange[500, 10^4, 100]
];

MhMin[{MS_, Mh_, DMh_}]  := {MS, Mh - DMh};
MhMax[{MS_, Mh_, DMh_}]  := {MS, Mh + DMh};
MhBest[{MS_, Mh_, DMh_}] := {MS, Mh};

dataMhMin  = MhMin  /@ data;
dataMhMax  = MhMax  /@ data;
dataMhBest = MhBest /@ data;

plot2 = ListLinePlot[dataMhBest,
                     PlotStyle -> {Red, Thick}];

plot1 = ListLinePlot[{dataMhMax, dataMhMin},
                     PlotStyle -> LightGray,
                     Filling -> {1 -> {{2}, LightGray}},
                     PlotRange -> All];

plot = Show[{plot1, plot2},
            BaseStyle -> {FontSize -> 16, FontFamily -> "Helvetica"},
            PlotLabel -> Style["\*SubscriptBox[X, t] = -2 \*SubscriptBox[M, S], tan\[Beta] = 5"],
            PlotRange -> Automatic,
            Axes -> False, Frame -> True,
            FrameLabel -> {Style["\*SubscriptBox[M, S] / GeV"],
                           Style["\*SubscriptBox[M, h] / GeV"]}];

Export["NUHMSSMNoFVHimalaya_Mh_MS.png", plot, ImageSize -> 600];
