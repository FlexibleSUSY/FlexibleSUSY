Get["models/MSSMEFTHiggs/MSSMEFTHiggs_librarylink.m"];

Mtpole = 173.34;

(* generate logarithmically spaced range [start, stop] *)
LogRange[start_, stop_, steps_] :=
    Exp /@ Range[Log[start], Log[stop],
                 (Log[stop] - Log[start])/steps];

(* generate logarithmically spaced range [2 Q, Q / 2] *)
GenerateScales[Q_] := LogRange[Q/2, 2 Q, 10];

(* run MSSMEFTHiggs spectrum generator *)
RunMSSMEFTHiggs[MS_, TB_, Xt_, Qpole_, Qmatch_] :=
    Module[{handle, spectrum},
           handle = FSMSSMEFTHiggsOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,
                   maxIterations -> 10000,
                   poleMassLoopOrder -> 2,
                   ewsbLoopOrder -> 2,
                   betaFunctionLoopOrder -> 3,
                   thresholdCorrectionsLoopOrder -> 2,
                   poleMassScale -> 0,
                   eftPoleMassScale -> Qpole,
                   eftMatchingScale -> Qmatch,
                   eftMatchingLoopOrderUp -> 1,
                   eftMatchingLoopOrderDown -> 1,
                   calculateBSMMasses -> 0
               },
               fsSMParameters -> {
                   Mt -> Mtpole
               },
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
                               {0    , 0    , MS/TB + Xt MS}},
                   AdInput -> MS TB IdentityMatrix[3],
                   AeInput -> MS TB IdentityMatrix[3]
               }
           ];
           spectrum = FSMSSMEFTHiggsCalculateSpectrum[handle];
           FSMSSMEFTHiggsCloseHandle[handle];
           spectrum
          ];

(* extract lightest Higgs pole mass Pole[M[hh]] from spectrum *)
RunMSSMEFTHiggsMh[pars__] :=
    (Pole[M[hh]] /. (MSSMEFTHiggs /. RunMSSMEFTHiggs[pars]))[[1]];

(* calculate Higgs mass and perform scale variation *)
RunMSSMEFTHiggsUncertainty[MS_, TB_, Xt_] :=
    Module[{MhMean, DMh, varyQpole, varyQmatch},
           MhMean = RunMSSMEFTHiggsMh[MS, TB, Xt, 0, 0];
           varyQpole  = RunMSSMEFTHiggsMh[MS, TB, Xt, #, 0]& /@ GenerateScales[Mtpole];
           varyQmatch = RunMSSMEFTHiggsMh[MS, TB, Xt, 0, #]& /@ GenerateScales[MS];
           (* combine uncertainty estimates *)
           DMh = Max[Abs[Max[varyQpole] - MhMean],
                     Abs[Min[varyQpole] - MhMean]] +
                 Max[Abs[Max[varyQmatch] - MhMean],
                     Abs[Min[varyQmatch] - MhMean]];
           { MhMean, DMh }
          ];

{Mh, DMh} = RunMSSMEFTHiggsUncertainty[2500, 20, Sqrt[6]];

Print["Mh = (", Mh, " +- ", DMh, ") GeV"];
