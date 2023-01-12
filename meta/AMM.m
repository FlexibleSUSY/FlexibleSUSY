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

BeginPackage["AMM`", {"SARAH`", "CXXDiagrams`", "TextFormatting`", "TreeMasses`", "LoopMasses`"}];

AMMGetMSUSY::usage="";
AMMContributingGraphs::usage="";
AMMContributingDiagramsForGraph::usage="";
CXXEvaluatorForDiagramFromGraph::usage="";
ForwardDeclaration::usage = "";

Begin["`Private`"];

barZeeGraph = {{0,0,0,1,0,0,0,0},
               {0,0,0,0,1,0,0,0},
               {0,0,0,0,0,0,0,1},
               {1,0,0,0,1,1,0,0},
               {0,1,0,1,0,0,1,0},
               {0,0,0,1,0,0,1,1},
               {0,0,0,0,1,1,0,1},
               {0,0,1,0,0,1,1,0}};

contributingGraphs = {barZeeGraph};

AMMContributingGraphs[] := contributingGraphs;

AMMContributingDiagramsForGraph[graph_, field_] :=
  Module[{diagrams},
    diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph,
         {1 -> field, 2 -> SARAH`AntiField[field], 3 -> GetPhoton[]}];

    Select[diagrams,IsDiagramSupported[graph,#] &]
 ]

IsDiagramSupported[barZeeGraph, diagram_] :=
  Module[{photonEmitter,exchangeParticle1,exchangeParticle2},

    photonEmitter = diagram[[8,3]];
    exchangeParticle1 = diagram[[4,3]];
    exchangeParticle2 = diagram[[5,3]];
    If[diagram[[8]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
        Return[False]];
    If[exchangeParticle1 == TreeMasses`GetPhoton[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsFermion[photonEmitter] && photonEmitter == diagram[[7,2]],
        Return[True]];
    If[exchangeParticle1 == TreeMasses`GetZBoson[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsFermion[photonEmitter] && photonEmitter == diagram[[7,2]],
        Return[True]];
    If[exchangeParticle1 == TreeMasses`GetPhoton[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsScalar[photonEmitter] && photonEmitter == diagram[[7,2]],
        Return[True]];
    If[exchangeParticle1 == TreeMasses`GetPhoton[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsVector[photonEmitter] && photonEmitter == diagram[[7,2]],
        Return[True]];

    Return[False];
  ];

CXXEvaluatorForDiagramFromGraph[diagram_, barZeeGraph] :=
  Module[{photonEmitter,exchangeParticle1,exchangeParticle2},

    photonEmitter = diagram[[8,3]];
    exchangeParticle1 = diagram[[4,3]];
    exchangeParticle2 = diagram[[5,3]];
    If[diagram[[8]] =!= {TreeMasses`GetPhoton[],CXXDiagrams`LorentzConjugate[photonEmitter],photonEmitter},
        Return["(unknown diagram)"]];
    If[exchangeParticle1 == TreeMasses`GetPhoton[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsFermion[photonEmitter],
        Return[CXXEvaluatorBZFL[photonEmitter,exchangeParticle2]]];
    If[exchangeParticle1 == TreeMasses`GetZBoson[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsFermion[photonEmitter],
        Return[CXXEvaluatorBZFLZ[photonEmitter,exchangeParticle2]]];
    If[exchangeParticle1 == TreeMasses`GetPhoton[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsScalar[photonEmitter],
        Return[CXXEvaluatorBZSL[photonEmitter,exchangeParticle2]]];
    If[exchangeParticle1 == TreeMasses`GetPhoton[] && TreeMasses`IsScalar[exchangeParticle2] && TreeMasses`IsVector[photonEmitter],
        Return[CXXEvaluatorBZVL[photonEmitter,exchangeParticle2]]];

    Return["(unknown diagram)"];
  ];

CXXEvaluatorBZFL[photonEmitter_,exchangeParticle_] :=
  "AMMBarZeeFermionLoop<Lepton," <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">";

CXXEvaluatorBZFLZ[photonEmitter_,exchangeParticle_] :=
  "AMMBarZeeFermionLoopZ<Lepton," <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">";

CXXEvaluatorBZSL[photonEmitter_,exchangeParticle_] :=
  "AMMBarZeeScalarLoop<Lepton," <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">";

CXXEvaluatorBZVL[photonEmitter_,exchangeParticle_] :=
  "AMMBarZeeVectorLoop<Lepton," <>
  CXXDiagrams`CXXNameOfField[photonEmitter] <> ", " <>
  CXXDiagrams`CXXNameOfField[exchangeParticle] <> ">";

GetMinMass[particle_] :=
    Module[{dim = TreeMasses`GetDimension[particle],
            mStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]],
            tail},
           If[dim == 1,
              "model.get_" <> mStr <> "()",
              tail = ToString[GetDimension[particle] - GetDimensionStartSkippingGoldstones[particle] + 1];
              "model.get_" <> mStr <> "().tail<" <> tail <> ">().minCoeff()"
             ]
          ];

AMMGetMSUSY[] :=
    Module[{susyParticles},
           susyParticles = Select[TreeMasses`GetSusyParticles[], IsElectricallyCharged];
           If[susyParticles === {},
              "return 0.;",
              "return Min(" <>
                 StringRiffle[GetMinMass /@ susyParticles, ", "] <>
              ");"
             ]
          ];

ForwardDeclaration[field_Symbol, functionName_String] :=
"template double " <> FlexibleSUSY`FSModelName <> "_amm::" <> functionName <> "<" <>
   CXXDiagrams`CXXNameOfField[field, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <>
   ">(" <>
   "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates&, const softsusy::QedQcd&, const Spectrum_generator_settings&" <> If[TreeMasses`GetDimension[field] =!= 1, ", int", ""] <> ");"

End[];
EndPackage[];

