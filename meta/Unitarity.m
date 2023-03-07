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

BeginPackage["Unitarity`", {"SARAH`", "Parameters`", "CConversion`", "TreeMasses`", "TextFormatting`", "Vertices`", "CXXDiagrams`"}];

GetScatteringMatrix::usage = "";

Begin["`Private`"];

InfiniteS[a0Input_, generationSizes_, FSScatteringPairs_] := Module[{params = Parameters`FindAllParametersClassified[a0Input], paramsCPP, mixingCPP, a0, decrementIndices},
   (* CPP definitions of parameters present in the expression *)
   paramsCPP =
      StringJoin[
         ("const auto " <> ToString@CConversion`ToValidCSymbol[#] <>
            " = model.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "();\n")& /@ Join[Parameters`FSModelParameters /. params, Parameters`FSOutputParameters /. params]
      ];
   decrementIndices = (#[i_, j_] /; Or@@(IntegerQ/@{i,j}) :> #@@(If[IntegerQ@#, #-1, #]& /@ {i, j}))& /@ Select[Join[Parameters`FSModelParameters /. params, Parameters`FSOutputParameters /. params], Length[GetParameterDimensions[#]]===2&];

   (* replace input parameters with their FS names *)
   a0 = a0Input /. Thread[(Parameters`FSModelParameters /. params) -> CConversion`ToValidCSymbol /@ (Parameters`FSModelParameters /. params)];
   a0 = a0 /. Susyno`LieGroups`conj -> Conj;
   a0 = a0 //. SARAH`sum[idx_, start_, stop_, exp_] :> FlexibleSUSY`SUM[idx, start-1, stop-1, exp];
   (* only for test
   a0 = a0 /. Delta[c1_?SarahColorIndexQ, c2_?SarahColorIndexQ] :> 1;*)
   resultInfinite = "";
   For[i=1, i<=Length[a0], i++,
      For[j=i, j<=Length[a0[[i]]], j++,
         If[a0[[i,j]] =!= 0,
            resultInfinite = resultInfinite <> "\n// " <> ToString[FSScatteringPairs[[i]]] <> "->" <> ToString[FSScatteringPairs[[j]]] <> "\n" <>
                     If[Count[generationSizes[[i,j]], x_/;x>1] === 0, "{\n", ""] <>
                        If[generationSizes[[i,j,1]] > 1,
                           "for (int in1=0; in1<" <> CXXNameOfField[First@FSScatteringPairs[[i]] /. Susyno`LieGroups`conj -> Identity] <> "::numberOfGenerations; ++in1) {\n",
                           ""
                        ] <>
                        If[generationSizes[[i,j,2]] > 1,
                           Nest[TextFormatting`IndentText, "for (int in2=0; in2<" <> CXXNameOfField[Last@FSScatteringPairs[[i]] /. Susyno`LieGroups`conj -> Identity] <> "::numberOfGenerations; ++in2) {\n", If[generationSizes[[i,j,1]] > 1, 1, 0]],
                           ""
                        ] <>
                        If[generationSizes[[i,j,3]] > 1,
                           Nest[TextFormatting`IndentText, "for (int out1=0; out1<" <> CXXNameOfField[First@FSScatteringPairs[[j]] /. Susyno`LieGroups`conj -> Identity] <> "::numberOfGenerations; ++out1) {\n", If[generationSizes[[i,j,1]] > 1, 1, 0] +If[generationSizes[[i,j,2]] > 1, 1, 0] ],
                           ""
                        ] <>
                        If[generationSizes[[i,j,4]] > 1,
                           Nest[TextFormatting`IndentText, "for (int out2=0; out2<" <> CXXNameOfField[Last@FSScatteringPairs[[j]] /. Susyno`LieGroups`conj -> Identity] <> "::numberOfGenerations; ++out2) {\n", If[generationSizes[[i,j,1]] > 1, 1, 0] +If[generationSizes[[i,j,2]] > 1, 1, 0] + If[generationSizes[[i,j,3]] > 1, 1, 0]],
                           ""
                        ] <>
                        Nest[
                           TextFormatting`IndentText,
                              "double temp = std::real(" <> ToString@CForm[FullSimplify[a0[[i,j]]] /. decrementIndices] <> ");\n" <>
                              "if (std::abs(matrix.coeff(" <> ToString[i-1] <> ", " <> ToString[j-1] <> ")) < std::abs(temp)) {\n" <>
                                 TextFormatting`IndentText["matrix.coeffRef(" <> ToString[i-1] <> ", " <> ToString[j-1] <> ") = temp;\n"] <>
                              "}\n",
                           If[MemberQ[{0, 1}, Count[generationSizes[[i,j]], x_/;x>1]], 1, Count[generationSizes[[i,j]], x_/;x>1]]
                        ] <>
                        If[generationSizes[[i,j,4]] > 1,
                           Nest[TextFormatting`IndentText, "}\n", If[generationSizes[[i,j,1]] > 1, 1, 0] +If[generationSizes[[i,j,2]] > 1, 1, 0] + If[generationSizes[[i,j,3]] > 1, 1, 0]],
                           ""
                        ] <>
                        If[generationSizes[[i,j,3]] > 1,
                           Nest[TextFormatting`IndentText, "}\n", If[generationSizes[[i,j,1]] > 1, 1, 0] +If[generationSizes[[i,j,2]] > 1, 1, 0] ],
                           ""
                        ] <>
                        If[generationSizes[[i,j,2]] > 1,
                           Nest[TextFormatting`IndentText, "}\n", If[generationSizes[[i,j,1]] > 1, 1, 0]],
                           ""
                        ] <>
                        If[Count[generationSizes[[i,j]], x_/;x>1] ===0 || generationSizes[[i,j, 1]] > 1,
                           "}\n",
                           ""
                        ]
         ]
      ]
   ];
   TextFormatting`IndentText[paramsCPP <> "\n" <> resultInfinite]
];

GetScatteringMatrix[] := Module[{generationSizes, a0, a0InfiniteS, FSScatteringPairs},
   InitUnitarity[];

   (* only color neutral final states *)
   FSScatteringPairs = Select[scatteringPairs, (TreeMasses`GetColorRepresentation /@ #) == {S,S}&];

   a0 = Outer[GetScatteringDiagrams[#1 -> #2]&, FSScatteringPairs, FSScatteringPairs, 1];
   (* obviously 's' lives in Susyno`LieGroups` *)
   a0InfiniteS = Map[Limit[# /. qChan -> 1 /. sChan|tChan|uChan -> 0, Susyno`LieGroups`s->Infinity]&, a0, {2}];

   generationSizes = Table[{i, j}, {i, Length[FSScatteringPairs]}, {j, Length[FSScatteringPairs]}];
   generationSizes = Apply[Join[TreeMasses`GetDimension /@ FSScatteringPairs[[#1]], TreeMasses`GetDimension /@ FSScatteringPairs[[#2]]]&, generationSizes, {2}];

   {Length[FSScatteringPairs], InfiniteS[a0InfiniteS, generationSizes, FSScatteringPairs]}
];

End[];
EndPackage[];

