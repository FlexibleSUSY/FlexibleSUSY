Get["models/HSSUSY/HSSUSY_librarylink.m"];
Get["model_files/HSSUSY/HSSUSY_uncertainty_estimate.m"];

settings = {
    precisionGoal -> 1.*^-5,
    poleMassLoopOrder -> 4,
    ewsbLoopOrder -> 4,
    betaFunctionLoopOrder -> 5,
    thresholdCorrectionsLoopOrder -> 3,
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
            AbInput -> 0,
            AtauInput -> 0,
            msq2 -> MS^2 IdentityMatrix[3],
            msu2 -> MS^2 IdentityMatrix[3],
            msd2 -> MS^2 IdentityMatrix[3],
            msl2 -> MS^2 IdentityMatrix[3],
            mse2 -> MS^2 IdentityMatrix[3],
            LambdaLoopOrder -> 3,
            TwoLoopAtAs -> 1,
            TwoLoopAbAs -> 1,
            TwoLoopAtAb -> 1,
            TwoLoopAtauAtau -> 1,
            TwoLoopAtAt -> 1,
            ThreeLoopAtAsAs -> 1
        }
   ];

LinearRange[start_, stop_, steps_] :=
    Range[start, stop, (stop - start)/steps];

Xtt = Sqrt[6];
TB  = 20;

data = ParallelMap[
    { N[#], Sequence @@ HSSUSYCalcMh[#, TB, Xtt] }&,
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
            PlotLabel -> Style["\*SubscriptBox[X, t] = 2.44949 \*SubscriptBox[M, S], tan\[Beta] = 20"],
            PlotRange -> Automatic,
            Axes -> False, Frame -> True,
            FrameLabel -> {Style["\*SubscriptBox[M, S] / GeV"],
                           Style["\*SubscriptBox[M, h] / GeV"]}];

Export["HSSUSY_Mh_MS.png", plot, ImageSize -> 600];
