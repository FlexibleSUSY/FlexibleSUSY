Get["models/HSSUSY/HSSUSY_librarylink.m"];

(* generate logarithmically spaced range [start, stop] *)
LogRange[start_, stop_, steps_] :=
    Exp /@ Range[Log[start], Log[stop],
                 (Log[stop] - Log[start])/steps];

(* generate logarithmically spaced range [2 Q, Q / 2] *)
GenerateScales[Q_] := LogRange[Q/2, 2 Q, 10];

CalcMh[MS_, TB_, Xt_, ytLoops_, asLoops_, Qpole_] :=
    Module[{handle, spec},
    handle = FSHSSUSYOpenHandle[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            calculateStandardModelMasses -> 1,
            poleMassLoopOrder -> 2,
            ewsbLoopOrder -> 2,
            betaFunctionLoopOrder -> 3,
            thresholdCorrectionsLoopOrder -> 3,
            poleMassScale -> Qpole,
            thresholdCorrections -> 120111021 +
                ytLoops * 10^6 + asLoops * 10^2
        },
        fsModelParameters -> {
            TanBeta -> TB,
            MEWSB -> 173.34,
            MSUSY -> MS,
            M1Input -> MS,
            M2Input -> MS,
            M3Input -> MS,
            MuInput -> MS,
            mAInput -> MS,
            AtInput -> (Xt + 1/TB) * MS,
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
    spec = FSHSSUSYCalculateSpectrum[handle];
    FSHSSUSYCloseHandle[handle];
    If[spec =!= $Failed, Pole[M[hh]] /. (HSSUSY /. spec), 0]
];

(* calculate Higgs mass with uncertainty estimate *)
CalcDMh[MS_, TB_, Xt_] :=
    Module[{Mh, MhYt3L, MhAs3L, varyQpole, DMh},
           Mh     = CalcMh[MS, TB, Xt, 2, 2, 0];
           MhYt3L = CalcMh[MS, TB, Xt, 3, 2, 0];
           MhAs3L = CalcMh[MS, TB, Xt, 2, 3, 0];
           varyQpole = CalcMh[MS, TB, Xt, 2, 2, #]& /@
                       GenerateScales[173.34];
           (* combine uncertainty estimates *)
           DMh = Max[Abs[Max[varyQpole] - Mh],
                     Abs[Min[varyQpole] - Mh]] +
                 Abs[Mh - MhYt3L] + Abs[Mh - MhAs3L];
           { Mh, DMh }
          ];

LaunchKernels[];
DistributeDefinitions[CalcDMh];

data = {
    ParallelMap[{#, Sequence @@ CalcDMh[1000 , 5, #]}&, Range[-3.5, 3.5, 0.1]],
    ParallelMap[{#, Sequence @@ CalcDMh[2000 , 5, #]}&, Range[-3.5, 3.5, 0.1]],
    ParallelMap[{#, Sequence @@ CalcDMh[10000, 5, #]}&, Range[-3.5, 3.5, 0.1]]
};

Export["HSSUSY_uncertainty_Mh_Xt_MS-1000.dat", data[[1]]];
Export["HSSUSY_uncertainty_Mh_Xt_MS-2000.dat", data[[2]]];
Export["HSSUSY_uncertainty_Mh_Xt_MS-10000.dat", data[[3]]];
