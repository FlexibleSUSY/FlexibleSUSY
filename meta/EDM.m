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

BeginPackage["EDM`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`"}];

EDMCreateInterfaceFunctionForField::usage="";
EDMContributingDiagramsForFieldAndGraph::usage="";
EDMContributingGraphs::usage="";

Begin["`Private`"];

(* The graphs that contribute to the EDM are precisely those with three
   external lines given by the field in question, its Lorentz conjugate
   and a photon.
   They are given as a List of undirected adjacency matrices where
    1 is the field itself
    2 is its Lorentz conjugate
    3 is the photon
   and all other indices unspecified. *)
vertexCorrectionGraph = {{0,0,0,1,0,0},
                         {0,0,0,0,1,0},
                         {0,0,0,0,0,1},
                         {1,0,0,0,1,1},
                         {0,1,0,1,0,1},
                         {0,0,1,1,1,0}};
contributingGraphs = {vertexCorrectionGraph};

EDMContributingGraphs[] := contributingGraphs

EDMContributingDiagramsForFieldAndGraph[field_,graph_] :=
  Module[{diagrams},
    diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph,
         {1 -> field, 2 -> CXXDiagrams`LorentzConjugate[field], 3 -> TreeMasses`GetPhoton[]}];

    Select[diagrams,IsDiagramSupported[field,graph,#] &]
 ]

IsDiagramSupported[field_,vertexCorrectionGraph,diagram_] :=
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
       Return[False]];
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[True]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[True]];
    
    Return[False];
  ]

EDMCreateInterfaceFunctionForField[field_,gTaggedDiagrams_List] :=
	Module[{prototype,definition,numberOfIndices = CXXDiagrams`NumberOfFieldIndices[field]},
		prototype = "double calculate_edm_" <> CXXNameOfField[field] <>
			"(" <> If[TreeMasses`GetDimension[field] =!= 1,
				" int generationIndex, ",
				" "] <>
			"const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model );";
                 
    definition = "double calculate_edm_" <> CXXNameOfField[field] <>
			"(" <> If[TreeMasses`GetDimension[field] =!= 1,
				" int generationIndex, ",
				" "] <>
			"const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model )\n" <>
			"{\n" <> IndentText[
				"context_base context{ model };\n" <>
				"std::array<int, " <> ToString @ numberOfIndices <> "> indices = {" <>
					If[TreeMasses`GetDimension[field] =!= 1,
						" generationIndex" <>
						If[numberOfIndices =!= 1,
							StringJoin @ Table[", 0", {numberOfIndices-1}],
							""] <> " ",
						If[numberOfIndices =!= 0,
							StringJoin @ Riffle[Table[" 0", {numberOfIndices}], ","] <> " ",
							""]
					] <> "};\n\n" <>

				"double val = 0.0;\n\n" <>

				"using namespace " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields;\n\n" <>

				StringJoin @ Riffle[Module[{graph = #[[1]], diagrams = #[[2]]},
					StringJoin @ Riffle[Module[{diagram = #, indexedDiagram},
						indexedDiagram = CXXDiagrams`IndexDiagramFromGraph[diagram, graph];
						
						"val += " <> 
						ToString @ ProjectColourFactor[field,
							CXXDiagrams`ColourFactorForIndexedDiagramFromGraph[indexedDiagram, graph]] <>
						" * " <> 
						CXXEvaluatorForFieldAndDiagramFromGraph[field, diagram, graph] <>
						"::value(indices, context);"
					] & /@ diagrams, "\n"]
				] & /@ gTaggedDiagrams, "\n"] <> "\n\n" <>

				"return val;"
			] <> "\n}";
    
    {prototype, definition}
  ];

(* Because of colour conservation, we know that the resulting colour
 * structure is either trivial or a ``SARAH`Delta[]``.
 * A more general approach should use Coefficient[] instead.
 *)
ProjectColourFactor[field_, colourFactor_] := colourFactor /;
	!TreeMasses`ColorChargedQ[field]
	
ProjectColourFactor[field_, colourFactor_ * SARAH`Delta[cIndex1, cIndex2]] := 
	colourFactor /;
	TreeMasses`ColorChargedQ[field] && 
		Vertices`SarahColorIndexQ[cIndex1] && Vertices`SarahColorIndexQ[cIndex2]

CXXEvaluatorForFieldAndDiagramFromGraph[field_,diagram_,vertexCorrectionGraph] := 
  Module[{photonEmitter,exchangeParticle},
    photonEmitter = diagram[[4,3]]; (* Edge between vertices 4 and 6 (3rd edge of vertex 4) *)
    exchangeParticle = diagram[[4,2]]; (* Edge between vertices 4 and 5 (2nd edge of vertex 4) *)
    
    If[diagram[[6]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
       Return["(unknown diagram)"]];
    If[TreeMasses`IsFermion[photonEmitter] && TreeMasses`IsScalar[exchangeParticle],
       Return[CXXEvaluatorFS[field,photonEmitter,exchangeParticle]]];
    If[TreeMasses`IsFermion[exchangeParticle] && TreeMasses`IsScalar[photonEmitter],
       Return[CXXEvaluatorSF[field,photonEmitter,exchangeParticle]]];
    
    Return["(unknown diagram)"];
  ]

CXXEvaluatorFS[field_,photonEmitter_,exchangeParticle_] :=
  "EDMVertexCorrectionFS<" <> CXXDiagrams`CXXNameOfField[field] <> ", " <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">"

CXXEvaluatorSF[field_,photonEmitter_,exchangeParticle_] :=
  "EDMVertexCorrectionSF<" <> CXXDiagrams`CXXNameOfField[field] <> ", " <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">"

End[];
EndPackage[];
