Get["models/HSSUSY/HSSUSY_librarylink.m"];

CalcMh[TB_, Xtt_, MS_] := Module[{handle, spec},
    handle = FSHSSUSYOpenHandle[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            calculateStandardModelMasses -> 1,
            poleMassLoopOrder -> 2,
            ewsbLoopOrder -> 2,
            betaFunctionLoopOrder -> 3,
            thresholdCorrectionsLoopOrder -> 2,
            poleMassScale -> 173.34,
            parameterOutputScale -> 173.34
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
            AtInput -> (Xtt + 1/TB) * MS,
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

LaunchKernels[];
DistributeDefinitions[CalcMh];

data = {
    ParallelMap[{#, CalcMh[5, #, 1000 ]}&, Range[-3.5, 3.5, 0.1]],
    ParallelMap[{#, CalcMh[5, #, 2000 ]}&, Range[-3.5, 3.5, 0.1]],
    ParallelMap[{#, CalcMh[5, #, 10000]}&, Range[-3.5, 3.5, 0.1]]
};

Export["HSSUSY_Mh_Xt_MS-1000.dat", data[[1]]];
Export["HSSUSY_Mh_Xt_MS-2000.dat", data[[2]]];
Export["HSSUSY_Mh_Xt_MS-10000.dat", data[[3]]];
