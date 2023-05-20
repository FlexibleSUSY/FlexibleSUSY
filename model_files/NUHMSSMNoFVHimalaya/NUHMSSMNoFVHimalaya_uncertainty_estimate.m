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
  Functions to perform an uncertainty estimate of NUHMSSMNoFVHimalaya.
 *)

CalcNUHMSSMNoFVHimalayaDMh::usage="\
The function takes the parameter point as input (same syntax as
FSNUHMSSMNoFVHimalayaOpenHandle[]) and returns the 2-component list { Mh, DMh }
where Mh is the Higgs mass and DMh is an uncertainty estimate of
missing 3-loop and 4-loop corrections.

Example: Peform a parameter scan over the SUSY scale in the interval
[1000, 10^10] GeV for tan(beta) = 20 and Xt/MS = Sqrt[6].

Get[\"models/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_librarylink.m\"];
Get[\"model_files/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_uncertainty_estimate.m\"];

CalcMh[MS_, TB_, Xtt_] :=
    CalcNUHMSSMNoFVHimalayaDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            poleMassLoopOrder -> 3,
            ewsbLoopOrder -> 3,
            thresholdCorrectionsLoopOrder -> 2,
            thresholdCorrections -> 122111121
        },
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
CalcNUHMSSMNoFVHimalayaMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcNUHMSSMNoFVHimalayaMh[a, Sequence @@ s, r];

CalcNUHMSSMNoFVHimalayaMh[asLoops_?NumericQ, ytLoops_?NumericQ, Qpole_?NumericQ, args__] :=
    Module[{handle, spec, tc},
           tc = thresholdCorrections /. { args };
           tc = If[IntegerQ[tc], tc,
                   thresholdCorrections /. Options[FSNUHMSSMNoFVHimalayaOpenHandle]];
           If[asLoops >= 0, tc = SetDigit[tc, 2, asLoops]];
           If[ytLoops >= 0, tc = SetDigit[tc, 6, ytLoops]];
           handle = FSNUHMSSMNoFVHimalayaOpenHandle[args];
           FSNUHMSSMNoFVHimalayaSet[handle,
               fsSettings -> {
                   calculateStandardModelMasses -> 0,
                   calculateBSMMasses -> 1,
                   thresholdCorrectionsLoopOrder -> 3,
                   poleMassScale -> Qpole,
                   thresholdCorrections -> tc
               }
           ];
           spec = FSNUHMSSMNoFVHimalayaCalculateSpectrum[handle];
           FSNUHMSSMNoFVHimalayaCloseHandle[handle];
           If[spec === $Failed, $Failed,
              (Pole[M[hh]] /. (NUHMSSMNoFVHimalaya /. spec))[[1]]]
          ];

(* calculate Higgs mass and uncertainty estimate *)
CalcNUHMSSMNoFVHimalayaDMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcNUHMSSMNoFVHimalayaDMh[a, Sequence @@ s, r];

CalcNUHMSSMNoFVHimalayaDMh[args__] :=
    Module[{Mh0, Mh, MhAs, MhYt, varyQpole, DMh,
            mhLoops = poleMassLoopOrder /. { args }, ytLoops, asLoops,
            MS = Qin /. { args }},
           ytLoops    = Max[mhLoops - 1, 1];
           asLoops    = Max[mhLoops - 2, 1];
           Mh0        = CalcNUHMSSMNoFVHimalayaMh[-1         , -1         , 0, args];
           If[Mh0 === $Failed, Return[{ $Failed, $Failed }]];
           Mh         = CalcNUHMSSMNoFVHimalayaMh[asLoops    , ytLoops    , 0, args];
           If[Mh === $Failed, Return[{ Mh0, $Failed }]];
           MhAs       = CalcNUHMSSMNoFVHimalayaMh[asLoops + 1, ytLoops    , 0, args];
           If[MhAs === $Failed, Return[{ Mh0, $Failed }]];
           MhYt       = CalcNUHMSSMNoFVHimalayaMh[asLoops    , ytLoops + 1, 0, args];
           If[MhYt === $Failed, Return[{ Mh0, $Failed }]];
           varyQpole  = CalcNUHMSSMNoFVHimalayaMh[asLoops    , ytLoops    , #, args]& /@
                        LogRange[MS/2, 2 MS, 10];
           If[MemberQ[varyQpole, $Failed], Return[{ Mh0, $Failed }]];
           (* combine uncertainty estimates *)
           DMh   = Max[Abs[Max[varyQpole] - Mh],
                       Abs[Min[varyQpole] - Mh]] +
                   Abs[Mh - MhAs] + Abs[Mh - MhYt];
           { Mh0, DMh }
          ];
