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
  Functions to perform an uncertainty estimate of MRSSMEFTHiggs.
 *)

CalcMRSSMEFTHiggsDMh::usage="\
The function takes the parameter point as input (same syntax as
FSMRSSMEFTHiggsOpenHandle[]) and returns the 2-component list { Mh, DMh }
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
[100, 10^4] GeV for tan(beta) = 10 and lambda_i = -0.5.

Get[\"models/MRSSMEFTHiggs/MRSSMEFTHiggs_librarylink.m\"];
Get[\"model_files/MRSSMEFTHiggs/MRSSMEFTHiggs_uncertainty_estimate.m\"];

CalcMh[MSUSY_, TB_, lam_] :=
    CalcMRSSMEFTHiggsDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            thresholdCorrectionsLoopOrder -> 2,
            thresholdCorrections -> 122111121
        },
        fsModelParameters -> {
            TanBeta -> TB,
            MS -> MSUSY,
            LamTDInput -> lam,
            LamTUInput -> lam,
            LamSDInput -> lam,
            LamSUInput -> lam,
            MuDInput -> -MSUSY,
            MuUInput -> MSUSY,
            BMuInput -> MSUSY^2,
            mq2Input -> scaleFactor MSUSY^2 IdentityMatrix[3],
            mu2Input -> scaleFactor MSUSY^2 IdentityMatrix[3],
            ml2Input -> MSUSY^2 IdentityMatrix[3],
            md2Input -> MSUSY^2 IdentityMatrix[3],
            me2Input -> MSUSY^2 IdentityMatrix[3],
            mS2Input -> scaleFactor MSUSY^2,
            moc2Input -> scaleFactor MSUSY^2,
            mT2Input -> MSUSY^2,
            mRd2Input -> MSUSY^2,
            mRu2Input -> MSUSY^2,
            MDBSInput -> MSUSY,
            MDWBTInput -> MSUSY,
            MDGocInput -> MSUSY
        }
   ];

LaunchKernels[];
DistributeDefinitions[CalcMh];

data = ParallelMap[
    { N[#], CalcMh[#, 10, -0.5] }&,
    LogRange[10^2, 10^4, 50]
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
CalcMRSSMEFTHiggsMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcMRSSMEFTHiggsMh[a, Sequence @@ s, r];

CalcMRSSMEFTHiggsMh[ytLoops_?NumericQ, Qpole_?NumericQ, Qm_?NumericQ, args__] :=
    Module[{handle, spec, tc},
           tc = thresholdCorrections /. { args };
           tc = If[IntegerQ[tc], tc,
                   thresholdCorrections /. Options[FSMRSSMEFTHiggsOpenHandle]];
           handle = FSMRSSMEFTHiggsOpenHandle[args];
           FSMRSSMEFTHiggsSet[handle,
               fsSettings -> {
                   calculateStandardModelMasses -> 0,
                   calculateBSMMasses -> 0,
                   thresholdCorrectionsLoopOrder -> 3,
                   eftPoleMassScale -> Qpole,
                   eftMatchingScale -> Qm,
                   thresholdCorrections -> SetDigit[tc, 6, ytLoops]
               }
           ];
           spec = FSMRSSMEFTHiggsCalculateSpectrum[handle];
           FSMRSSMEFTHiggsCloseHandle[handle];
           If[spec === $Failed, $Failed,
              (Pole[M[hh]] /. (MRSSMEFTHiggs /. spec))[[1]]]
          ];

(* calculate Higgs mass and uncertainty estimate *)
CalcMRSSMEFTHiggsDMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcMRSSMEFTHiggsDMh[a, Sequence @@ s, r];

CalcMRSSMEFTHiggsDMh[args__] :=
    Module[{Mh, MhYt3L, varyQpole, varyQmatch,
            DMhSM, DMhSUSY,
            MSUSY = MS /. { args }, Mlow = Mt /. { args }},
           If[!NumericQ[Mlow],
              Mlow = Mt /. Options[FSMRSSMEFTHiggsOpenHandle]
             ];
           Mh         = CalcMRSSMEFTHiggsMh[2, 0, 0, args];
           If[Mh === $Failed, Return[{ $Failed, $Failed }]];
           MhYt3L     = CalcMRSSMEFTHiggsMh[3, 0, 0, args];
           If[MhYt3L === $Failed, Return[{ Mh, $Failed }]];
           varyQpole  = CalcMRSSMEFTHiggsMh[2, #, 0, args]& /@
                        LogRange[Mlow/2, 2 Mlow, 10];
           If[MemberQ[varyQpole, $Failed], Return[{ Mh, $Failed }]];
           varyQmatch = CalcMRSSMEFTHiggsMh[2, 0, #, args]& /@
                        LogRange[MSUSY/2, 2 MSUSY, 10];
           If[MemberQ[varyQmatch, $Failed], Return[{ Mh, $Failed }]];
           (* combine uncertainty estimates *)
           DMhSM   = Max[Abs[Max[varyQpole] - Mh],
                         Abs[Min[varyQpole] - Mh]] +
                     Abs[Mh - MhYt3L];
           DMhSUSY = Max[Abs[Max[varyQmatch] - Mh],
                         Abs[Min[varyQmatch] - Mh]];
           { Mh, DMhSM + DMhSUSY }
          ];
