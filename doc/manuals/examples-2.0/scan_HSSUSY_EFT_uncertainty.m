Get["models/HSSUSY/HSSUSY_librarylink.m"];

(* generate logarithmically spaced range [start, stop] *)
LogRange[start_, stop_, steps_] :=
    Exp /@ Range[Log[start], Log[stop],
                 (Log[stop] - Log[start])/steps];

CalcMh[MS_, TB_, Xt_, deltaEFT_] :=
    Module[{handle, spec},
    handle = FSHSSUSYOpenHandle[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            calculateStandardModelMasses -> 1,
            poleMassLoopOrder -> 2,
            ewsbLoopOrder -> 2,
            betaFunctionLoopOrder -> 3,
            thresholdCorrectionsLoopOrder -> 3,
            thresholdCorrections -> 122111221
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
            TwoLoopAtAt -> 1,
            DeltaEFT -> deltaEFT
        }
    ];
    spec = FSHSSUSYCalculateSpectrum[handle];
    FSHSSUSYCloseHandle[handle];
    If[spec =!= $Failed, Pole[M[hh]] /. (HSSUSY /. spec), 0]
];

(* calculate Higgs mass with uncertainty estimate *)
CalcDMh[MS_, TB_, Xt_] :=
    Module[{Mh, MhEFT},
           Mh    = CalcMh[MS, TB, Xt, 0];
           MhEFT = CalcMh[MS, TB, Xt, 1];
           { Mh, Abs[Mh - MhEFT] }
          ];

LaunchKernels[];
DistributeDefinitions[CalcDMh];

data = ParallelMap[{#, Sequence @@ CalcDMh[#, 5, Sqrt[6]]}&,
                   LogRange[173.34, 10^5, 60]];

Export["HSSUSY_EFT_uncertainty.dat", data];
