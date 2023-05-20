deltaLambdaQmatch1L = With[{
    k = 1/(4 Pi)^2
    },
    k Get["meta/SM/HSSUSY_scale_variation_1L.m"] /. {
        ml2[a__] :> msl2[a],
        me2[a__] :> mse2[a],
        mq2[a__] :> msq2[a],
        mu2[a__] :> msu2[a],
        md2[a__] :> msd2[a],
        MS -> MSUSY,
        Q -> Qmatch,
        gt -> Yu[3,3],
        gb -> Yd[3,3],
        gtau -> Ye[3,3],
        tb -> TanBeta,
        mA -> mAInput,
        M1 -> M1Input,
        M2 -> M2Input,
        M3 -> M3Input,
        \[Mu] -> MuInput
    }
];

deltaLambdaQmatch2L = With[{
    k = 1/(4 Pi)^2
    },
    k^2 Get["meta/SM/HSSUSY_scale_variation_2L.m"] /. {
        Derivative[1][TCF[1]] -> TCD1F[1],
        Derivative[1][TCF[2]] -> TCD1F[2],
        Derivative[1][TCF[3]] -> TCD1F[3],
        Derivative[1][TCF[4]] -> TCD1F[4],
        Derivative[1][TCF[5]] -> TCD1F[5],
        Derivative[1][TCf[1]] -> TCD1f[1],
        Derivative[1][TCf[2]] -> TCD1f[2],
        Derivative[1][TCf[3]] -> TCD1f[3],
        Derivative[1][TCf[4]] -> TCD1f[4],
        Derivative[1][TCf0]   -> TCD1f0,
        Derivative[1][TCg0]   -> TCD1g0,
        Derivative[1,0][TCf[i_]] :> TCD10f[i],
        Derivative[0,1][TCf[i_]] :> TCD01f[i],
        ml2[a__] :> msl2[a],
        me2[a__] :> mse2[a],
        mq2[a__] :> msq2[a],
        mu2[a__] :> msu2[a],
        md2[a__] :> msd2[a],
        MS -> MSUSY,
        Q -> Qmatch,
        gt -> Yu[3,3],
        gb -> Yd[3,3],
        gtau -> Ye[3,3],
        tb -> TanBeta,
        mA -> mAInput,
        M1 -> M1Input,
        M2 -> M2Input,
        M3 -> M3Input,
        \[Mu] -> MuInput
    }
];

deltaLambdaQmatch3L = With[{
    k = 1/(4 Pi)^2
    },
    k^3 Get["meta/SM/HSSUSY_scale_variation_3L.m"] /. {
        Derivative[1][TCF[1]] -> TCD1F[1],
        Derivative[1][TCF[2]] -> TCD1F[2],
        Derivative[2][TCF[1]] -> TCD2F[1],
        Derivative[2][TCF[2]] -> TCD2F[2],
        ml2[a__] :> msl2[a],
        me2[a__] :> mse2[a],
        mq2[a__] :> msq2[a],
        mu2[a__] :> msu2[a],
        md2[a__] :> msd2[a],
        MS -> MSUSY,
        Q -> Qmatch,
        gt -> Yu[3,3],
        tb -> TanBeta,
        mA -> mAInput,
        M3 -> M3Input,
        \[Mu] -> MuInput
    }
];
