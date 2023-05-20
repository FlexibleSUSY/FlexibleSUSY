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
  Functions to perform an uncertainty estimate of HSSUSY.
 *)

CalcHSSUSYDMh::usage="\
The function takes the parameter point as input (same syntax as
FSHSSUSYOpenHandle[]) and returns the 2-component list { Mh, DMh }
where Mh is the Higgs mass and DMh is an uncertainty estimate of
missing higher-order corrections.

The uncertainty estimate takes three sources into account:

 * SM uncertainty: Missing higher order corrections from the
   extraction of the running SM parameters and to the calculation of
   the Higgs pole mass.

 * EFT uncertainty: Missing terms of the order O(v^2/MSUSY^2).

 * SUSY uncertainty: Missing 3-loop contributions to the quartic Higgs
   coupling \[Lambda] from SUSY particles.

Important note: The uncertainty estimate assumes that all 2-loop
threshold corrections at the SUSY scale are enabled:

  fsModelParameters -> {
      LambdaLoopOrder -> 2,
      TwoLoopAtAs -> 1,
      TwoLoopAbAs -> 1,
      TwoLoopAtAb -> 1,
      TwoLoopAtauAtau -> 1,
      TwoLoopAtAt -> 1
  }

Example: Peform a parameter scan over the SUSY scale in the interval
[1000, 10^10] GeV for tan(beta) = 20 and Xt/MS = Sqrt[6].

Get[\"models/HSSUSY/HSSUSY_librarylink.m\"];
Get[\"model_files/HSSUSY/HSSUSY_uncertainty_estimate.m\"];

CalcMh[MS_, TB_, Xtt_] :=
    CalcHSSUSYDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            thresholdCorrectionsLoopOrder -> 4,
            thresholdCorrections -> 124111421
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
CalcHSSUSYMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcHSSUSYMh[a, Sequence @@ s, r];

CalcHSSUSYMh[ytLoops_?NumericQ, Qpole_?NumericQ, Qm_?NumericQ, eft_?NumericQ, ytMSSM_?NumericQ, args__] :=
    Module[{handle, spec, tc},
           tc = thresholdCorrections /. { args };
           tc = If[IntegerQ[tc], tc,
                   thresholdCorrections /. Options[FSHSSUSYOpenHandle]];
           If[ytLoops >= 0, tc = SetDigit[tc, 6, ytLoops]];
           handle = FSHSSUSYOpenHandle[args];
           FSHSSUSYSet[handle,
               fsSettings -> {
                   calculateStandardModelMasses -> 1,
                   thresholdCorrectionsLoopOrder -> 4,
                   poleMassScale -> Qpole,
                   thresholdCorrections -> tc
               },
               fsModelParameters -> {
                   DeltaEFT -> eft,
                   DeltaYt -> ytMSSM,
                   Qmatch -> Qm
               }
           ];
           spec = FSHSSUSYCalculateSpectrum[handle];
           FSHSSUSYCloseHandle[handle];
           If[spec === $Failed, $Failed,
              Pole[M[hh]] /. (HSSUSY /. spec)]
          ];

(* calculate Higgs mass and uncertainty estimate *)
CalcHSSUSYDMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcHSSUSYDMh[a, Sequence @@ s, r];

CalcHSSUSYDMh[args__] :=
    Module[{Mh0, Mh, MhYt3L, MhEFT, MhYtMSSM, varyQpole, varyQmatch,
            DMhSM, DMhEFT, DMhSUSY,
            MS = MSUSY /. { args }, Mlow = MEWSB /. { args }},
           Mh0        = CalcHSSUSYMh[-1, 0, 0, 0, 0, args];
           If[Mh0 === $Failed, Return[{$Failed, $Failed}]];
           Mh         = CalcHSSUSYMh[3, 0, 0, 0, 0, args];
           If[Mh === $Failed, Return[{Mh0, $Failed}]];
           MhYt3L     = CalcHSSUSYMh[4, 0, 0, 0, 0, args];
           If[MhYt3L === $Failed, Return[{Mh0, $Failed}]];
           MhEFT      = CalcHSSUSYMh[3, 0, 0, 1, 0, args];
           If[MhEFT === $Failed, Return[{Mh0, $Failed}]];
           MhYtMSSM   = CalcHSSUSYMh[3, 0, 0, 0, 1, args];
           If[MhYtMSSM === $Failed, Return[{Mh0, $Failed}]];
           varyQpole  = CalcHSSUSYMh[3, #, 0, 0, 0, args]& /@
                        LogRange[Mlow/2, 2 Mlow, 10];
           varyQmatch = CalcHSSUSYMh[3, 0, #, 0, 0, args]& /@
                        LogRange[MS/2, 2 MS, 10];
           varyQpole  = Select[varyQpole , NumericQ];
           varyQmatch = Select[varyQmatch, NumericQ];
           If[varyQmatch === {} || varyQpole === {}, Return[{Mh0, $Failed}]];
           (* combine uncertainty estimates *)
           DMhSM   = Max[Abs[Max[varyQpole] - Mh],
                         Abs[Min[varyQpole] - Mh]] +
                     Abs[Mh - MhYt3L];
           DMhEFT  = Abs[Mh - MhEFT];
           DMhSUSY = Max[Abs[Max[varyQmatch] - Mh],
                         Abs[Min[varyQmatch] - Mh]] +
                     Abs[Mh - MhYtMSSM];
           { Mh0, DMhSM + DMhEFT + DMhSUSY }
          ];
