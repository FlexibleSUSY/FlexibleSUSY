(* ::Package:: *)

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

BeginPackage["FFVFormFactors`",
   {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`", "Utils`"}
];

FFVFormFactorsCreateInterfaceFunction::usage = "";
FFVGraphs::usage = "Returns adjacency matrix, which is wrapped in a List
for a triangle diagram which can be seen by typing a command:
AdjacencyGraph[FFVFormFactors`FFVGraphs[][[1]], VertexLabels -> \"Name\"]

why such a weird numbering?
because some procedure uses this order for fermion number flow.";
FFVContributingDiagramsForGraph::usage = "";

Begin["Private`"];

vertexCorrectionGraph = {
   {0, 1, 0, 0, 0, 0},
   {1, 0, 1, 0, 1, 0},
   {0, 1, 0, 0, 1, 1},
   {0, 0, 0, 0, 1, 0},
   {0, 1, 1, 1, 0, 0},
   {0, 0, 1, 0, 0, 0}
};
contributingGraphs = {vertexCorrectionGraph};
FFVGraphs[] := contributingGraphs;

EmitterL[diagram_] := diagram[[3,2]];
EmitterR[diagram_] := diagram[[2,2]];
Spectator[diagram_] := diagram[[2,3]];

FFVContributingDiagramsForGraph[graph_, Fj_ -> {Fi_, V_}] :=
   Module[{diagrams},
      diagrams =
         CXXDiagrams`FeynmanDiagramsOfType[
            graph,
            {4 -> Fj, 1 -> CXXDiagrams`LorentzConjugate[Fi], 6 -> CXXDiagrams`LorentzConjugate[V]}
         ];

      (* for now we only support selected classes of diagrams *)
      Select[diagrams, IsDiagramSupported[graph, #]&]
   ];

IsDiagramSupported[graph_, diagram_] :=
   Module[{photonEmitter, exchangeParticle, photonEmitterAfter},

      photonEmitter = diagram[[3,2]]; (* Edge between vertices ? and ? (3rd edge of vertex 4) *)
      photonEmitterAfter = diagram[[2,2]];
      exchangeParticle = diagram[[2,3]]; (* Edge between vertices ? and ? (2nd edge of vertex 4) *)

      If[TreeMasses`IsFermion[photonEmitter] &&
         TreeMasses`IsFermion[photonEmitterAfter]
         && TreeMasses`IsScalar[exchangeParticle],
         Return[True]
      ];

      If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter] &&TreeMasses`IsScalar[photonEmitterAfter],
         Return[True]
      ];

      Print["Warning: Diagram with internal particles of type ",
         StringJoin @@ (ToString /@ SARAH`getType /@ {EmitterL[diagram], EmitterR[diagram], Spectator[diagram]})];
      Print["         is currently not supported."];
      Print["         Discarding diagram with particles ",
         {EmitterL[diagram], EmitterR[diagram], Spectator[diagram]}, "."];
      Return[False];
   ];

FFVFormFactorsCreateInterfaceFunction[Fj_ -> {Fi_, V_}, topologies_, diagrams_] :=
   Module[{prototype, definition,
           numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[Fj],
           numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[Fi],
           numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[V], temp = {}},

      Utils`AssertWithMessage[Length[topologies] === Length[diagrams],
         "Length of diagrams should be the same as length of topologies"];

      For[i = 1, i <= Length[topologies], i++,
         For[j = 1, j <= Length[diagrams[[i]]], j++,
                 temp = temp <> CreateCall[Fj, Fi, V, topologies[[i]], diagrams[[i,j]]];
            ];
         ];
      prototype =
         "std::valarray<std::complex<double>> calculate_" <> CXXNameOfField[Fj] <>
            "_" <> CXXNameOfField[Fi] <> "_" <> CXXNameOfField[V] <> "_form_factors (\n" <>
            IndentText[
               If[TreeMasses`GetDimension[Fj] =!= 1,
                  "int generationIndex1, ",
                  " "
               ] <>
               If[TreeMasses`GetDimension[Fi] =!= 1,
                  "int generationIndex2,\n",
                  " "
               ] <>
               "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, bool discard_SM_contributions)"
            ];

      definition =
         prototype <> " {\n\n" <>
            IndentText[
               "context_base context {model};\n" <>
               "std::array<int, " <> ToString @ numberOfIndices1 <> "> indices1 = {" <>
                       If[TreeMasses`GetDimension[Fj] =!= 1,
                          "generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ", "],
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[Fi] =!= 1,
                          "generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ", "] <> " ",
                             ""]
                         ] <> "};\n\n" <>

               "std::valarray<std::complex<double>> val {0.0, 0.0, 0.0, 0.0};\n\n" <>

               temp <>

               "\n" <>
               "return val;"

            ] <> "\n}\n\n";

      {prototype <> ";", definition}
   ];

CreateCall[Fj_, Fi_, V_, topology_, diagram_] :=
   "val += std::complex<double> " <>
      ToString @ N @ FSReIm @ CXXDiagrams`ExtractColourFactor @ CXXDiagrams`ColorFactorForDiagram[topology, diagram] <> " * FFV_" <>
      StringJoin @@ (ToString /@ SARAH`getType /@ {EmitterL[diagram], EmitterR[diagram], Spectator[diagram]}) <> "<" <>
         StringJoin @ Riffle[CXXDiagrams`CXXNameOfField /@ {Fj, Fi, V, EmitterL[diagram], EmitterR[diagram], Spectator[diagram]}, ","]  <>
                     ">::value(indices1, indices2, context, discard_SM_contributions);\n";

End[];
EndPackage[];
