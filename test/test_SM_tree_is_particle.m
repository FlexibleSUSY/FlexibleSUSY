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

Needs["TestSuite`", "TestSuite.m"];
Block[{Print}, (* We do not test SARAH, so let us reduce output. *)
   Needs["SARAH`"];
   Start@"SM";
];

FlexibleSUSY`FSEnableColors = False; (* Desired by Utils`FSFancyWarning *)
FlexibleSUSY`FSEigenstates = EWSB; (* Required by TreeMasses`Private`IsOfType *)
Needs["TreeMasses`", "TreeMasses.m"];

check[function_, particle_] :=
Module[{typeQ, message, stream},
   Print["\ncheck"];
   stream = OpenWrite[];
   AppendTo[$Output, stream];
      typeQ = function@particle;
   $Output = Most@$Output;
      message = ReadString[stream[[1]]];
   Close[stream];
   {typeQ, MatchQ[message, _String]}
];

With[{BAR = SARAH`bar, CONJ = Susyno`LieGroups`conj},
   Utils`FSFancyLine[];
   TestEquality[{True, True}, check[TreeMasses`IsScalar, BAR@SARAH`hh]];
   (*            ^~~~ Is it a scalar?
    *                  ^~~~ Should a warning be created?
    *)
   TestEquality[{True, False}, check[TreeMasses`IsScalar, CONJ@SARAH`hh]];
   TestEquality[{True, False}, check[TreeMasses`IsScalar, SARAH`hh]];

   TestEquality[{True, True}, check[TreeMasses`IsFermion, CONJ@Global`Fe]];
   TestEquality[{True, False}, check[TreeMasses`IsFermion, BAR@Global`Fe]];
   TestEquality[{True, False}, check[TreeMasses`IsFermion, Global`Fe]];

   TestEquality[{True, True}, check[TreeMasses`IsVector, BAR@SARAH`VP]];
   TestEquality[{True, False}, check[TreeMasses`IsVector, CONJ@SARAH`VP]];
   TestEquality[{True, False}, check[TreeMasses`IsVector, SARAH`VP]];

   TestEquality[{True, False}, check[TreeMasses`IsGhost, CONJ@Global`gG]];
   TestEquality[{True, False}, check[TreeMasses`IsGhost, BAR@Global`gG]];
   TestEquality[{True, False}, check[TreeMasses`IsGhost, Global`gG]];

   TestEquality[{True, True}, check[TreeMasses`IsGoldstone, BAR@SARAH`Ah]];
   TestEquality[{True, False}, check[TreeMasses`IsGoldstone, CONJ@SARAH`Ah]];
   TestEquality[{True, False}, check[TreeMasses`IsGoldstone, SARAH`Ah]];
];

PrintTestSummary[];
