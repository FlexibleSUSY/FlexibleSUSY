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

BeginPackage["TwoLoopSM`"];
EndPackage[];

GetSMHiggsMass::usage = "Returns the loop corrections to the Higgs
 mass in the Standard Model.  The loop-corrections are taken from
 arxiv:1205.6497 and arxiv:1504.05200.

Note: The return value contains the contributions from tadpole
 diagrams.

Usage: GetSMHiggsMass[loopOrder -> {1,1,1}, corrections -> {1,1}, simplifications -> { p -> 0 }]

Parameters:

- loopOrder (optional): List of factors multiplied by each loop order.
  #1: tree-level
  #2: 1-loop level
  #3: 2-loop level
  (default: {1,1,1})

- corrections (optional): List of factors multiplied by each 2-loop
  correction. (default: {1, 1})
  #1: alpha_t * alpha_s (arxiv:1205.6497, Eq. (20))
  #2: alpha_t^2         (arxiv:1504.05200, Eq. (20))

- simplifications (optional): List of replacement rules to simplify
  the result (default: { p -> 0, B0[0,m12_,m22_,Q2_] :>
     B0zero[Sqrt[m12],Sqrt[m22],Sqrt[Q2]] }).
  Example: simplifications -> {} (* no simplifications *)
";

(* MS-bar parameters *)
{ p, yt, mt, g3, Q };

(* loop functions *)
{ B0 };

(* options *)
{ loopOrder, corrections };

Begin["TwoLoopSM`Private`"];

(* B0 for p = 0 *)
B0zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1 - m2],
          Log[mu^2/m2^2],
          PossibleZeroQ[m1],
          1 + Log[mu^2/m2^2],
          PossibleZeroQ[m2],
          1 + Log[mu^2/m1^2],
          True,
          1 + (m1^2 Log[mu^2/m1^2] - m2^2 Log[mu^2/m2^2])/(m1^2 - m2^2)
   ];

Options[GetSMHiggsMass] = {
    loopOrder -> {1,1,1},
    corrections -> {1,1},
    simplifications -> { p -> 0,
                         B0[0,m12_,m22_,Q2_] :> B0zero[Sqrt[m12],Sqrt[m22],Sqrt[Q2]] }
};

GetSMHiggsMass1LAlphaT[] :=
    Module[{h = 1/(4 Pi)^2},
           h 3 yt^2 (4 mt^2 - p^2) B0[p^2,mt^2,mt^2,Q^2]
          ];

GetSMHiggsMass2LAlphaTAlphaS[] :=
    Module[{LogT = Log[mt^2/Q^2], h = 1/(4 Pi)^2},
           h^2 2 mt^2 16 g3^2 yt^2 (3 LogT^2 + LogT)
          ];

GetSMHiggsMass2LAlphaTAlphaT[] :=
    Module[{LogT = Log[mt^2/Q^2], h = 1/(4 Pi)^2},
           h^2 2 mt^2 (-3 yt^4 (3 LogT^2 - 7 LogT + 2 + Pi^2/3))
          ];

GetSMHiggsMass[OptionsPattern[]] :=
    (
        OptionValue[loopOrder][[1]] 0 +
        OptionValue[loopOrder][[2]] GetSMHiggsMass1LAlphaT[] +
        OptionValue[loopOrder][[3]] OptionValue[corrections][[1]] GetSMHiggsMass2LAlphaTAlphaS[] +
        OptionValue[loopOrder][[3]] OptionValue[corrections][[2]] GetSMHiggsMass2LAlphaTAlphaT[]
    ) //. OptionValue[simplifications];

End[];
