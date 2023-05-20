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
  Functions to perform an uncertainty estimate of NUHNMSSM.
 *)

CalcNUHNMSSMDMh::usage="\
The function takes the parameter point as input (same syntax as
FSNUHNMSSMOpenHandle[]) and returns the 2-component list { Mh, DMh }
where Mh is the Higgs mass and DMh is an uncertainty estimate of
missing 3-loop and 4-loop corrections.

Example: Peform a parameter scan over the SUSY scale in the interval
[1000, 10^10] GeV for tan(beta) = 20 and Xt/MS = Sqrt[6].

Get[\"models/NUHNMSSM/NUHNMSSM_librarylink.m\"];
Get[\"model_files/NUHNMSSM/NUHNMSSM_uncertainty_estimate.m\"];

CalcTlambda[mA2_, TB_, lam_, kap_, muEff_] :=
    Module[{vS = muEff Sqrt[2] / lam},
           (mA2 Sqrt[2] TB / (vS (TB^2 + 1)) - muEff kap)
           ];

CalcTkappa[m2_, lam_, muEff_] :=
    Module[{vS = muEff Sqrt[2] / lam},
           - Sqrt[2] m2 / vS
           ];

CalcMh[MS_, TB_, Xtt_] :=
    CalcNUHNMSSMDMh[
        fsSettings -> {
            precisionGoal -> 1.*^-5,
            poleMassLoopOrder -> 3,
            ewsbLoopOrder -> 3,
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
CalcNUHNMSSMMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcNUHNMSSMMh[a, Sequence @@ s, r];

CalcNUHNMSSMMh[asLoops_?NumericQ, ytLoops_?NumericQ, Qpole_?NumericQ, args__] :=
    Module[{handle, spec, tc},
           tc = thresholdCorrections /. { args };
           tc = If[IntegerQ[tc], tc,
                   thresholdCorrections /. Options[FSNUHNMSSMOpenHandle]];
           handle = FSNUHNMSSMOpenHandle[args];
           FSNUHNMSSMSet[handle,
               fsSettings -> {
                   calculateStandardModelMasses -> 0,
                   calculateBSMMasses -> 1,
                   thresholdCorrectionsLoopOrder -> 3,
                   poleMassScale -> Qpole,
                   thresholdCorrections -> SetDigit[SetDigit[tc, 2, asLoops], 6, ytLoops]
               }
           ];
           spec = FSNUHNMSSMCalculateSpectrum[handle];
           FSNUHNMSSMCloseHandle[handle];
           If[spec === $Failed, $Failed,
              (Pole[M[hh]] /. (NUHNMSSM /. spec))[[1]]]
          ];

(* calculate Higgs mass and uncertainty estimate *)
CalcNUHNMSSMDMh[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    CalcNUHNMSSMDMh[a, Sequence @@ s, r];

CalcNUHNMSSMDMh[args__] :=
    Module[{Mh, MhAs, MhYt, varyQpole, DMh,
            mhLoops = poleMassLoopOrder /. { args }, ytLoops, asLoops,
            MS = MSUSY /. { args }},
           ytLoops    = Max[mhLoops - 1, 1];
           asLoops    = Max[mhLoops - 2, 1];
           Mh         = CalcNUHNMSSMMh[asLoops    , ytLoops    , 0, args];
           MhAs       = CalcNUHNMSSMMh[asLoops + 1, ytLoops    , 0, args];
           MhYt       = CalcNUHNMSSMMh[asLoops    , ytLoops + 1, 0, args];
           varyQpole  = CalcNUHNMSSMMh[asLoops    , ytLoops    , #, args]& /@
                        LogRange[MS/2, 2 MS, 10];
           If[Mh === $Failed, Return[{ $Failed, $Failed }]];
           If[MhAs === $Failed || MhYt === $Failed || MemberQ[varyQpole, $Failed],
              Return[{ Mh, $Failed }]];
           (* combine uncertainty estimates *)
           DMh   = Max[Abs[Max[varyQpole] - Mh],
                       Abs[Min[varyQpole] - Mh]] +
                   Abs[Mh - MhAs] + Abs[Mh - MhYt];
           { Mh, DMh }
          ];
