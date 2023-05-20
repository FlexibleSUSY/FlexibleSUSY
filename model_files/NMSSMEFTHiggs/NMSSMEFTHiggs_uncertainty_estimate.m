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

(**
  Functions to perform an uncertainty estimate of NMSSMEFTHiggs.
 *)

CalcNMSSMEFTHiggsDMh::usage="\
The function takes the parameter point as input (same syntax as
FSNMSSMEFTHiggsOpenHandle[]) and returns the 2-component list { Mh, DMh }
where Mh is the Higgs mass and DMh is an uncertainty estimate of
missing 2-loop corrections.

The uncertainty estimate takes two sources into account:

 * SM uncertainty: Missing higher order corrections from the
   extraction of the running SM parameters and to the calculation of
   the Higgs pole mass.

 * SUSY uncertainty: Missing 3-loop contributions to the quartic Higgs
   coupling \[Lambda] from SUSY particles.

An EFT uncertainty does not exist in in the FlexibleEFTHiggs approach.

Example: Peform a parameter scan over the SUSY scale in the interval
[1000, 10^10] GeV for tan(beta) = 20 and Xt/MS = Sqrt[6].

Get[\"models/NMSSMEFTHiggs/NMSSMEFTHiggs_librarylink.m\"];
Get[\"model_files/NMSSMEFTHiggs/NMSSMEFTHiggs_uncertainty_estimate.m\"];

CalcTlambda[mA2_, TB_, lam_, kap_, muEff_] :=
    Module[{vS = muEff Sqrt[2] / lam},
           (mA2 Sqrt[2] TB / (vS (TB^2 + 1)) - muEff kap)
           ];

CalcTkappa[m2_, lam_, muEff_] :=
    Module[{vS = muEff Sqrt[2] / lam},
           - Sqrt[2] m2 / vS
           ];

CalcMh[MS_, TB_, Xtt_] :=
    CalcNMSSMEFTHiggsDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            thresholdCorrectionsLoopOrder -> 2,
            thresholdCorrections -> 122111121
        },
        fsModelParameters -> {
            MSUSY   -> MS,
            M1Input -> MS,
            M2Input -> MS,
            M3Input -> MS,
            MuInput -> MS,
            TanBeta -> TB,
            LambdaInput -> lam,
            KappaInput -> kap,
            ALambdaInput -> CalcTlambda[MS^2, TB, lam, kap, MS] / lam,
            AKappaInput -> CalcTkappa[MS^2, lam, MS] / kap,
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

LaunchKernels[];
DistributeDefinitions[CalcMh];

data = ParallelMap[
    { N[#], CalcMh[#, 20, Sqrt[6]] }&,
    LogRange[10^3, 10^10, 50]
];
";

(* get digit of [num] at position [pos] *)
GetDigit[num_, pos_, base_:10] :=
    IntegerPart[Mod[num / base^pos, base]];

(* set digit of [num] at position [pos] to [val] *)
SetDigit[num_, pos_, val_, base_:10] :=
    num + (val - GetDigit[num,pos,base]) base^pos;

(* generate logarithmically spaced range [start, stop] *)
LogRange[start_, stop_, steps_] :=
    Exp /@ Range[Log[start], Log[stop], (Log[stop] - Log[start])/steps];

(* calculate Higgs mass *)
CalcNMSSMEFTHiggsMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcNMSSMEFTHiggsMh[a, Sequence @@ s, r];

CalcNMSSMEFTHiggsMh[ytLoops_?NumericQ, Qpole_?NumericQ, Qm_?NumericQ, args__] :=
    Module[{handle, spec, tc},
           tc = thresholdCorrections /. { args };
           tc = If[IntegerQ[tc], tc,
                   thresholdCorrections /. Options[FSNMSSMEFTHiggsOpenHandle]];
           handle = FSNMSSMEFTHiggsOpenHandle[args];
           FSNMSSMEFTHiggsSet[handle,
               fsSettings -> {
                   calculateStandardModelMasses -> 0,
                   calculateBSMMasses -> 0,
                   thresholdCorrectionsLoopOrder -> 3,
                   eftPoleMassScale -> Qpole,
                   eftMatchingScale -> Qm,
                   thresholdCorrections -> SetDigit[tc, 6, ytLoops]
               }
           ];
           spec = FSNMSSMEFTHiggsCalculateSpectrum[handle];
           FSNMSSMEFTHiggsCloseHandle[handle];
           If[spec === $Failed, $Failed,
              (Pole[M[hh]] /. (NMSSMEFTHiggs /. spec))[[1]]]
          ];

(* calculate Higgs mass and uncertainty estimate *)
CalcNMSSMEFTHiggsDMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcNMSSMEFTHiggsDMh[a, Sequence @@ s, r];

CalcNMSSMEFTHiggsDMh[args__] :=
    Module[{Mh, MhYt3L, varyQpole, varyQmatch,
            DMhSM, DMhSUSY,
            MS = MSUSY /. { args }, Mlow = Mt /. { args }},
           If[!NumericQ[Mlow],
              Mlow = Mt /. Options[FSNMSSMEFTHiggsOpenHandle]
             ];
           Mh         = CalcNMSSMEFTHiggsMh[2, 0, 0, args];
           If[Mh === $Failed, Return[{ $Failed, $Failed }]];
           MhYt3L     = CalcNMSSMEFTHiggsMh[3, 0, 0, args];
           If[MhYt3L === $Failed, Return[{ Mh, $Failed }]];
           varyQpole  = CalcNMSSMEFTHiggsMh[2, #, 0, args]& /@
                        LogRange[Mlow/2, 2 Mlow, 10];
           If[MemberQ[varyQpole, $Failed], Return[{ Mh, $Failed }]];
           varyQmatch = CalcNMSSMEFTHiggsMh[2, 0, #, args]& /@
                        LogRange[MS/2, 2 MS, 10];
           If[MemberQ[varyQmatch, $Failed], Return[{ Mh, $Failed }]];
           (* combine uncertainty estimates *)
           DMhSM   = Max[Abs[Max[varyQpole] - Mh],
                         Abs[Min[varyQpole] - Mh]] +
                     Abs[Mh - MhYt3L];
           DMhSUSY = Max[Abs[Max[varyQmatch] - Mh],
                         Abs[Min[varyQmatch] - Mh]];
           { Mh, DMhSM + DMhSUSY }
          ];
