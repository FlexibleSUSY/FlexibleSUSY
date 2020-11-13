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

(* @todo: rename diagram to insertion *)
BeginPackage["Decays`",
   {"SARAH`", "CConversion`", "CXXDiagrams`", "TreeMasses`", "TextFormatting`", "Utils`", "Vertices`"}
];

FSParticleDecay::usage="head used for storing details of an particle decay,
in the format
   FSParticleDecay[decaying particle, {final state particles}, {{loopOrder, {{topology, {diagrams}}, {topology, {diagrams}}}}, ...}]
where the topology is encoded as the relevant adjacency matrix.
";

IsSupportedDecayParticle::usage="returns True if decays for the given particle are supported.";

CreateCompleteParticleList::usage="";
GetDecaysForParticle::usage = "Creates 'objects' FSParticleDecay";
GetVerticesForDecays::usage="gets required vertices for a list of decays";

CreateSMParticleAliases::usage="creates aliases for SM particles present in model.";
CreateBSMParticleAliasList::usage="";

CallDecaysCalculationFunctions::usage="creates calls to functions calculating
decays of the given particles.";
CreateDecaysCalculationPrototypes::usage="creates prototypes for convenience
functions calculating all decays for each particle.";
CreateDecaysCalculationFunctions::usage="creates definitions for convenience
functions calculating all decays for each particle.";
CreatePartialWidthCalculationPrototypes::usage="create prototypes for
functions computing partial widths of all decays.";
CreatePartialWidthCalculationFunctions::usage="creates definitions for
functions computing partial widths of all decays.";
CreateDecaysGetterFunctions::usage="create getters for specific particle decays";
CreateDecayTableGetterPrototypes::usage="create getter prototypes for C++ decay table.";
CreateDecayTableGetterFunctions::usage="create getter definitions for C++ decay table.";
CreateDecayTableInitialization::usage="create C++ initializer for decay table.";
CreateTotalAmplitudeSpecializations::usage="creates specialized functions for higher-order decays.";
CreatePartialWidthSpecializations::usage="creates specialized functions for particular decays.";

Begin["`Private`"];

FieldPositionInVertex[field_, vertex_] := Module[{pos},

   pos = Position[vertex, field, {1}];
   Utils`AssertWithMessage[MatchQ[pos, {{_Integer}}],
      "Error! Could't find the location of field " <> ToString@field <> " in vertex " <> ToString@vertex <>
         ". Got " <> ToString@pos <> "."
   ];

   pos /. {{n_}} :> n
];

GetInitialState[FSParticleDecay[particle_, finalState_List, diagrams_List]] := particle;

GetFinalState[FSParticleDecay[particle_, finalState_List, diagrams_List]] := finalState;

(*
  Returns a list of the form

     {{integer loop order, {{topology one, {diagram one, diagram two, ...}}, ...}},
      {integer loop order, {{topology one, {diagram one, diagram two, ...}}, ...}}}

  where each sublist contains the topologies and corresponding diagrams contributing
  to the decay at that particular loop order.
*)
GetDecayTopologiesAndDiagrams[FSParticleDecay[particle_, finalState_List, diagrams_List]] := diagrams;

(*
  Returns a list of the form

     {{integer loop order, {diagram one, diagram two, ...}},
      {integer loop order, {diagram one, diagram two, ...}}, ...}

  where each sublist contains the diagrams contributing to the decay at that
  particular loop-order.
*)
GetDecayDiagramsOnly[FSParticleDecay[particle_, finalState_List, diagrams_List]] :=
    {#[[1]], Flatten[Last /@ #[[2]], 1]}& /@ diagrams;

(*
   Returns a list of the form

      {{topology one, {diagram one, diagram two, ...}},
       {topology two, {diagram one, diagram two, ...}}, ...}

   containing all topologies and diagrams contributing to the decay
   at the given loop order.
*)
GetDecayTopologiesAndDiagramsAtLoopOrder[decay_FSParticleDecay, loopOrder_Integer] :=
    Module[{toposAndDiags},
           toposAndDiags = Select[GetDecayTopologiesAndDiagrams[decay], (First[#] == loopOrder)&];
           If[toposAndDiags =!= {},
              toposAndDiags = Last[Flatten[toposAndDiags, 1]];
             ];
           toposAndDiags
          ];

(*
  Returns a list of the form

     {diagram one, diagram two, ...}

  containing all diagrams contributing to the decay at the given loop order.
*)
GetDecayDiagramsOnlyAtLoopOrder[decay_FSParticleDecay, loopOrder_Integer] :=
    Module[{diags},
           diags = Select[GetDecayDiagramsOnly[decay], (First[#] == loopOrder)&];
           If[diags =!= {},
              diags = Last[Flatten[diags, 1]];
             ];
           diags
          ];

GetPossibleDecayTopologies[nProducts_, nLoops_] :=
    (
     Print["Error: decay topology with ", nProducts, " particles and ", nLoops, " loops not supported"];
     Quit[1]
    );

(* tree-level two-body decay, with
   vertex 1 = incoming state
   vertex 2 = outgoing state
   vertex 3 = outgoing state
   vertex 4 = internal vertex *)
GetPossibleDecayTopologies[2, 0] := {
   {{0,0,0,1},
    {0,0,0,1},
    {0,0,0,1},
    {1,1,1,0}}
};

$translationFile = Once[Get["meta/generic_loop_decay_diagram_classes.m"]];
SetAtrributes[$translationFile, {Protected, Locked}];

(* get topologies as generated by FeynArts
   Some diagrams are 0 because of spin algebra yet they exist in the file - we filter them out *)
GetPossibleDecayTopologies[2, 1] :=
    DeleteDuplicates[#[[2]] & /@ Select[$translationFile, (#[[-3]] =!= {})&]];

IsTreeLevelDecayTopology[t_] := MemberQ[GetPossibleDecayTopologies[2, 0], t];
IsTreeLevelTwoBodyDecayTopology[t_] := IsTreeLevelDecayTopology[t];

IsOneLoopDecayTopology[t_] := MemberQ[GetPossibleDecayTopologies[2, 1], t];
IsOneLoopTwoBodyDecayTopology[t_] := IsOneLoopDecayTopology[t];

GetPossibleDecayTopologies[nProducts_] := Join @@ (GetPossibleDecayTopologies[nProducts, #]& /@ {0, 1});

GetTreeLevelDecayTopologyName[nFinalParticles_Integer] := "e" <> ToString[nFinalParticles + 1] <> "_l0_t1";

GetOneLoopDecayTopologyName[nFinalParticles_Integer, topology_] :=
    Module[{pos},
           pos = Position[GetPossibleDecayTopologies[nFinalParticles, 1], topology];
           If[pos === {},
              Print["Error: unknown one-loop topology."];
              Quit[1];
             ];
           pos = First[First[pos]];
           "e" <> ToString[nFinalParticles + 1] <> "_l1_t" <> ToString[pos]
          ];

GetDecayTopologyName[t_] :=
    Which[MemberQ[GetPossibleDecayTopologies[2, 0], t],
          GetTreeLevelDecayTopologyName[2],
          MemberQ[GetPossibleDecayTopologies[2, 1], t],
          GetOneLoopDecayTopologyName[2, t],
          True,
          Print["Error: unknown decay topology provided"];
          Quit[1];
         ];

CreateCompleteParticleList[particles_List] := DeleteDuplicates[Join[particles, SARAH`AntiField /@ particles]];

GenericScalarName[] := "scalar";
GenericVectorName[] := "vector";
GenericFermionName[] := "fermion";
GenericGhostName[] := "ghost";

SimplifiedName[Susyno`LieGroups`conj[particle_]] :=
    Susyno`LieGroups`conj[SimplifiedName[particle]];
SimplifiedName[SARAH`bar[particle_]] :=
    SARAH`bar[SimplifiedName[particle]];

SimplifiedName[particle_?TreeMasses`IsSMLepton] := "lep";
SimplifiedName[particle_?TreeMasses`IsSMDownQuark] := "dq";
SimplifiedName[particle_?TreeMasses`IsSMUpQuark] := "uq";
SimplifiedName[particle_ /; TreeMasses`GetHiggsBoson[] =!= Null && particle === TreeMasses`GetHiggsBoson[]] := "H";
SimplifiedName[particle_ /; TreeMasses`GetPseudoscalarHiggsBoson[] =!= Null && particle === TreeMasses`GetPseudoscalarHiggsBoson[]] := "AH";
SimplifiedName[particle_ /; TreeMasses`GetWBoson[] =!= Null && particle === TreeMasses`GetWBoson[]] := "W";
SimplifiedName[particle_ /; TreeMasses`GetZBoson[] =!= Null && particle === TreeMasses`GetZBoson[]] := "Z";
SimplifiedName[particle_ /; TreeMasses`GetPhoton[] =!= Null && particle === TreeMasses`GetPhoton[]] := "A";
SimplifiedName[particle_ /; TreeMasses`GetGluon[] =!= Null && particle === TreeMasses`GetGluon[]] := "G";
SimplifiedName[particle_] := particle;

CreateParticleAlias[particle_, namespace_String] :=
    "using " <> SimplifiedName[particle] <> " = " <>
    CXXDiagrams`CXXNameOfField[particle, prefixNamespace -> namespace] <> ";";

CreateParticleAliases[particles_, namespace_:""] :=
    Utils`StringJoinWithSeparator[CreateParticleAlias[#, namespace]& /@ particles, "\n"];

CreateSMParticleAliases[namespace_:""] :=
    Module[{smParticlesToAlias},
           smParticlesToAlias = Select[{TreeMasses`GetHiggsBoson[],
                                        TreeMasses`GetPseudoscalarHiggsBoson[],
                                        TreeMasses`GetWBoson[], TreeMasses`GetZBoson[],
                                        TreeMasses`GetGluon[], TreeMasses`GetPhoton[],
                                        TreeMasses`GetDownLepton[1] /. field_[generation_] :> field,
                                        TreeMasses`GetUpQuark[1] /. field_[generation_] :> field,
                                        TreeMasses`GetDownQuark[1] /.field_[generation_] :> field
                                       }, (# =!= Null)&];
           CreateParticleAliases[smParticlesToAlias, namespace]
          ];

CreateBSMParticleAliasList[namespace_:""] :=
   Module[{bsmForZdecay, bsmForWdecay},
      bsmForZdecay = Select[Prepend[#, TreeMasses`GetZBoson[]]& /@ DeleteDuplicates@Sort@Tuples[Join[TreeMasses`GetSusyParticles[], SARAH`AntiField /@ TreeMasses`GetSusyParticles[]], 2], IsPossibleNonZeroVertex[#, True]&];
      bsmForZdecay =
         Join[
            bsmForZdecay,
            Select[Prepend[#, TreeMasses`GetZBoson[]]& /@ DeleteDuplicates@Sort@Tuples[{Join[TreeMasses`GetSusyParticles[], SARAH`AntiField /@ TreeMasses`GetSusyParticles[]], Join[TreeMasses`GetSMParticles[], SARAH`AntiField /@ TreeMasses`GetSMParticles[]]}], IsPossibleNonZeroVertex[#, True]&]
         ];
      bsmForWdecay = Select[Prepend[#, TreeMasses`GetWBoson[]]& /@ DeleteDuplicates@Sort@Tuples[Join[TreeMasses`GetSusyParticles[], SARAH`AntiField /@ TreeMasses`GetSusyParticles[]], 2], IsPossibleNonZeroVertex[#, True]&];
      bsmForWdecay =
         Join[
            bsmForWdecay,
            Select[Prepend[#, TreeMasses`GetWBoson[]]& /@ DeleteDuplicates@Sort@Tuples[{Join[TreeMasses`GetSusyParticles[], SARAH`AntiField /@ TreeMasses`GetSusyParticles[]], Join[TreeMasses`GetSMParticles[], SARAH`AntiField /@ TreeMasses`GetSMParticles[]]}], IsPossibleNonZeroVertex[#, True]&]
         ];
      {
         Join[bsmForZdecay, bsmForWdecay],
      "// List of potential Z boson decay products excluding pure SM decays\n" <>
      "typedef boost::mpl::list<\n" <>
         TextFormatting`IndentText@StringJoin@Riffle[("boost::mpl::list<" <> CXXDiagrams`CXXNameOfField[#1, prefixNamespace -> namespace] <> ", " <> CXXDiagrams`CXXNameOfField[#2, prefixNamespace -> namespace] <> ">")& @@@ (Drop[#, {1}]& /@ bsmForZdecay), ",\n"] <> If[Length@bsmForZdecay > 0, "\n", ""] <> "> BSMForZdecay;\n\n" <>
      "// List of potential W boson decay products excluding pure SM decays\n" <>
      "typedef boost::mpl::list<\n" <>
         TextFormatting`IndentText@StringJoin@Riffle[("boost::mpl::list<" <> CXXDiagrams`CXXNameOfField[#1, prefixNamespace -> namespace] <> ", " <> CXXDiagrams`CXXNameOfField[#2, prefixNamespace -> namespace] <> ">")& @@@ (Drop[#, {1}]& /@ bsmForWdecay), ",\n"] <> If[Length@bsmForWdecay > 0, "\n", ""] <> "> BSMForWdecay;"
      }
   ];

GetGenericTypeName[p_?TreeMasses`IsScalar] := GenericScalarName[];
GetGenericTypeName[p_?TreeMasses`IsVector] := GenericVectorName[];
GetGenericTypeName[p_?TreeMasses`IsFermion] := GenericFermionName[];
GetGenericTypeName[p_?TreeMasses`IsGhost] := GenericGhostName[];

GetFeynArtsTypeName[p_?TreeMasses`IsScalar] := S;
GetFeynArtsTypeName[p_?TreeMasses`IsVector] := V;
GetFeynArtsTypeName[p_?TreeMasses`IsFermion] := F;
GetFeynArtsTypeName[p_?TreeMasses`IsGhost] := U;

IsSupportedDecayParticle[particle_] :=
    (
     !TreeMasses`IsGhost[particle] &&
     !TreeMasses`IsVector[particle] &&
     !TreeMasses`IsMassless[particle] &&
     TreeMasses`GetDimensionWithoutGoldstones[particle] > 0
    );

(* returns False if state consists only of Goldstones or ghosts, for any set of generation indices *)
IsPhysicalFinalState[finalState_List] :=
    Module[{goldstones},
           goldstones = Select[finalState, TreeMasses`IsGoldstone];
           isAlwaysGoldstone = (TreeMasses`GetDimensionWithoutGoldstones[#] == 0)& /@ goldstones;
           (Select[finalState, TreeMasses`IsGhost] === {}) && !(Or @@ isAlwaysGoldstone)
          ];

IsElectricChargeConservingDecay[initialParticle_, finalState_List] :=
    Module[{chargeSum},
           chargeSum = Simplify[Plus @@ (Join[{-TreeMasses`GetElectricCharge[initialParticle]}, TreeMasses`GetElectricCharge /@ finalState])];
           PossibleZeroQ[chargeSum]
          ];

(* @todo handle more than 2 particles in final state and non-SM color representations *)
IsColorInvariantDecay[initialParticle_, finalState_List] :=
    Module[{initialStateRep, finalStateReps, result = True},
           If[Length[finalState] == 2,
              initialStateRep = TreeMasses`GetColorRepresentation[initialParticle];
              finalStateReps = Sort[TreeMasses`GetColorRepresentation /@ finalState];
              Switch[initialStateRep,
                     S, result = ((finalStateReps === {S, S}) ||
                                  (finalStateReps === {T, T}) ||
                                  (finalStateReps === {-T, -T}) ||
                                  (finalStateReps === Sort[{-T, T}]) ||
                                  (finalStateReps === {O, O}));,
                     T|-T, result = ((finalStateReps === Sort[{T, S}]) ||
                                     (finalStateReps === Sort[{-T, S}]) ||
                                     (finalStateReps === Sort[{O, T}]) ||
                                     (finalStateReps === Sort[{O, -T}]));,
                     O, result = ((finalStateReps === Sort[{O, S}]) ||
                                  (finalStateReps === {T, T}) ||
                                  (finalStateReps === {-T, -T}) ||
                                  (finalStateReps === Sort[{-T, T}])
                                  (* uncomment to enable O -> OO decays once
                                     the handling of multiple color structures
                                     is introduced
                                  || (finalStateReps === Sort[{O, O}])*));,
                     _, result = True; (* unhandled case *)
                    ];
             ];
           result
          ];

(* don't generate decays like Fe3 -> Fe1 VP
   @todo: should we generate them? *)
IsSelfDecay[initialParticle_, finalState_List] :=
   (MemberQ[finalState, initialParticle] && (MemberQ[finalState, TreeMasses`GetPhoton[]] || MemberQ[finalState, TreeMasses`GetGluon[]]));

FinalStateContainsInitialState[initialParticle_, finalState_List] :=
    Module[{containsInitialMultiplet, dim},
           containsInitialMultiplet = !FreeQ[finalState, initialParticle];
           dim = TreeMasses`GetDimension[initialParticle];
           containsInitialMultiplet && dim == 1
          ];

IsPossibleNonZeroVertex[fields_List, useDependences_:False] :=
    Module[{numFields},
           numFields = Length[fields];
           Vertices`IsNonZeroVertex[fields, Vertices`GetCachedVertices[numFields, useDependences], useDependences]
         ];

IsPossibleNonZeroDiagram[diagram_, useDependences_:False] :=
    Module[{vertices},
           vertices = CXXDiagrams`VerticesForDiagram[diagram];
           And @@ (IsPossibleNonZeroVertex[#, useDependences]& /@ vertices)
          ];

IsPossibleTreeLevelDecay[decay_FSParticleDecay, useDependences_:False] :=
    Module[{treeLevelDiags = GetDecayDiagramsOnlyAtLoopOrder[decay, 0]},
           treeLevelDiags =!= {} && (And @@ (IsPossibleNonZeroDiagram[#, useDependences]& /@ treeLevelDiags))
          ];

IsPossibleOneLoopDecay[decay_FSParticleDecay] :=
    GetDecayDiagramsOnlyAtLoopOrder[decay, 1] =!= {};

ContainsOnlySupportedVertices[diagram_] :=
    Module[{vertices, vertexTypes, unsupportedVertices},
           vertices = CXXDiagrams`VerticesForDiagram[diagram];
           (* This function indirectly calls a GaugeStructureOfVertexLorentzPart.
              That's why we don't make GaugeStructureOfVertexLorentzPart function
                 to fail on unsuported color structure.
              We generated tree and 1-loop amplitude all at once and then select the
              lowest order one. For the process we considered so far the amplitudes with
              unsuported vertices don't contribitu in the end. *)
           vertexTypes = CXXDiagrams`VertexTypeForFields /@ vertices;
           unsupportedVertices = Complement[vertexTypes, CXXDiagrams`VertexTypes[]];
           If[unsupportedVertices =!= {},
              MapIndexed[(If[!MemberQ[CXXDiagrams`VertexTypes[], vertexTypes[[First[#2]]]],
                             Print["Warning: vertex with fields ", #1, " is not currently supported."];
                             Print["    Diagrams involving this vertex will be discarded."];
                            ];)&, vertices];
             ];
           unsupportedVertices === {}
          ];

IsSupportedDiagram[diagram_] := ContainsOnlySupportedVertices[diagram];

GetContributingDiagramsForDecayGraph[initialField_, finalFields_List, graph_] :=
   Module[
      {
         externalFields = Join[{1 -> initialField}, MapIndexed[(First[#2] + 1 -> #1)&, finalFields]],
         diagrams
      },
      (* vertices in diagrams are not SortCp'ed *)
      diagrams =
         CXXDiagrams`FeynmanDiagramsOfType[
            graph,
            externalFields,
            (* One loop decay topologies T2, T3 & T5 contain an A0 bubble on external leg.
                  With below argument set to True, charged particles are inserted twice in
               such bubble - once as particle and once as antiparticle. *)
            If[IsOneLoopDecayTopology[graph], !MemberQ[{"T2","T3","T5"}, FeynArtsTopologyName[graph]], True]
         ];
      Select[diagrams, IsPossibleNonZeroDiagram[#, True]&]
   ];

(* returns list of {{number of loops, {{topology, list of insertions}}}, ...} *)
GetContributingGraphsForDecay[initialParticle_, finalParticles_List, maxLoops_Integer] :=
    Module[{nFinalParticles = Length[finalParticles], topologies, diagrams},
           (* Until we include non-renormalizable operators, existence
              of tree-level amplitude for h->Zgamma or h->gamma gamma
              would violate Ward identity (as at the level of dim 4
              operators there can be only terms proportional to gmunu
              and not to p2mu p1nu). So even if SARAH claims that there
              is such a coupling, it must be 0. So for those processes
              we start generating amplitudes from the 1-loop level. *)
           minLoops =
              If[MemberQ[{TreeMasses`GetHiggsBoson[], TreeMasses`GetPseudoscalarHiggsBoson}, initialParticle] &&
                 (Sort@finalParticles === Sort[{TreeMasses`GetPhoton[], TreeMasses`GetPhoton[]}] ||
                 Sort@finalParticles === Sort[{TreeMasses`GetPhoton[], TreeMasses`GetZBoson[]}]),
                 1,
                 0
              ];
           topologies = Join[Table[{i, GetPossibleDecayTopologies[nFinalParticles, i]}, {i, minLoops, maxLoops}]];
           diagrams = {#[[1]], {#, GetContributingDiagramsForDecayGraph[initialParticle, SARAH`AntiField /@ finalParticles, #]}& /@ #[[2]]}&
                      /@ topologies;
           diagrams = {#[[1]], With[{toposAndDiags = #[[2]]}, Select[toposAndDiags, #[[2]] =!= {}&]]}& /@ diagrams;
           diagrams = DeleteCases[diagrams, {_Integer, {}}];
           {#[[1]], With[{toposAndDiags = #[[2]]}, {#[[1]], Select[#[[2]], IsSupportedDiagram]}& /@ toposAndDiags]}& /@ diagrams
          ];

GetContributingGraphsForDecay[initialParticle_, finalParticles_List] :=
    GetContributingGraphsForDecay[initialParticle, finalParticles, 1];

(* defines a fixed ordering for final state particles  *)
(* @todo decide on what this ordering actually will be *)
OrderFinalState[initialParticle_?TreeMasses`IsScalar, finalParticles_List] :=
    Module[{orderedFinalState},
           orderedFinalState = First[Vertices`SortCp[SARAH`Cp[Join[{initialParticle}, finalParticles]]]];
           orderedFinalState = Drop[orderedFinalState, First[Position[orderedFinalState, initialParticle]]];
           If[Length[orderedFinalState] === 2,
              (* re-order to SSV *)
              If[TreeMasses`IsVector[orderedFinalState[[1]]] && TreeMasses`IsScalar[orderedFinalState[[2]]],
                 orderedFinalState = Reverse[orderedFinalState]
                ];
              (* re-order to S bar[F] F *)
              If[TreeMasses`IsFermion[orderedFinalState[[1]]] && TreeMasses`IsFermion[orderedFinalState[[2]]],
                 If[Head[orderedFinalState[[2]]] === SARAH`bar && !Head[orderedFinalState[[1]]] === SARAH`bar,
                    orderedFinalState = Reverse[orderedFinalState];
                   ];
                ];
              If[TreeMasses`IsVector[orderedFinalState[[1]]] && TreeMasses`IsVector[orderedFinalState[[2]]],
                 If[Head[orderedFinalState[[2]]] === Susyno`LieGroups`conj && !Head[orderedFinalState[[1]]] === Susyno`LieGroups`conj,
                    orderedFinalState = Reverse[orderedFinalState];
                 ];
              ];
           ];
            orderedFinalState
          ];

OrderFinalState[initialParticle_?TreeMasses`IsFermion, finalParticles_List] :=
   Module[{orderedFinalState},
      orderedFinalState = First[Vertices`SortCp[SARAH`Cp[Join[{initialParticle}, finalParticles]]]];
      orderedFinalState = Drop[orderedFinalState, First[Position[orderedFinalState, initialParticle]]];
      If[Length[orderedFinalState] === 2,
         (* re-order to FFV *)
         If[TreeMasses`IsVector[orderedFinalState[[1]]] && TreeMasses`IsFermion[orderedFinalState[[2]]],
            orderedFinalState = Reverse[orderedFinalState]
         ];
         (* re-order to FFS *)
         If[TreeMasses`IsScalar[orderedFinalState[[1]]] && TreeMasses`IsFermion[orderedFinalState[[2]]],
            orderedFinalState = Reverse[orderedFinalState]
         ];
      ];
      orderedFinalState
   ];

OrderFinalState[initialParticle_, finalParticles_List] :=
    Module[{orderedFinalState},
           orderedFinalState = First[Vertices`SortCp[SARAH`Cp[Join[{initialParticle}, finalParticles]]]];
           Drop[orderedFinalState, First[Position[orderedFinalState, initialParticle]]]
          ];

GetDecaysForParticle[particle_, Infinity, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    );

GetDecaysForParticle[particle_, {Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    );

GetDecaysForParticle[particle_, {n_, Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    );

GetDecaysForParticle[particle_, {Infinity, n_}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    );

GetDecaysForParticle[particle_, {Infinity, Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    );

GetDecaysForParticle[particle_, maxNumberOfProducts_Integer /; maxNumberOfProducts >= 2,
                     allowedFinalStateParticles_List] :=
    GetDecaysForParticle[particle, {2, maxNumberOfProducts}, allowedFinalStateParticles];

GetDecaysForParticle[particle_, {minNumberOfProducts_Integer /; minNumberOfProducts >= 2,
                                 maxNumberOfProducts_Integer /; maxNumberOfProducts >= 2},
                                allowedFinalStateParticles_List] :=
    Module[{finalStateSizes},
           finalStateSizes = Table[{i}, {i, minNumberOfProducts, maxNumberOfProducts}];
           Flatten[GetDecaysForParticle[particle, #, allowedFinalStateParticles]& /@ finalStateSizes, 1]
          ];

(* Deletes diagrams from FSParticleDecay which are 0 because of color algebra.
   If all diagrams for a given topology are deleted, the topology is also removed *)
DeleteDiagramsVanishingDueToColor[decay_FSParticleDecay /; GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 0] === {} ]  :=
    Module[{wrong, topoAndinsertions, emptyTopos = {}},

        topoAndinsertions = GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1];

        For[i = 1, i <= Length@topoAndinsertions, i++,
           wrong = {};

          (* loop over diagrams *)
          For[j=1, j <= Length[topoAndinsertions[[i,2]]], j++,
            If[CXXDiagrams`ColorFactorForDiagram[topoAndinsertions[[i,1]], topoAndinsertions[[i,2, j]]] === 0,
              AppendTo[wrong, {j}];
            ]
          ];

          topoAndinsertions[[i,2]] = Delete[topoAndinsertions[[i,2]], wrong];
          If[Length[topoAndinsertions[[i,2]]]===0, AppendTo[emptyTopos, {i}]];
      ];

      topoAndinsertions = Delete[topoAndinsertions, emptyTopos];

      FSParticleDecay@@MapAt[
          (# /. {{1, x_}} :> {{1, topoAndinsertions}})&,
          List@@(decay),
          {3}
      ]
    ];
DeleteDiagramsVanishingDueToColor[decay_FSParticleDecay/; GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 0] =!= {}] := decay;

(* outputs list of FSParticleDecay objects *)
GetDecaysForParticle[particle_, {exactNumberOfProducts_Integer}, allowedFinalStateParticles_List] :=
    Module[{genericFinalStates,
            isPossibleDecay, concreteFinalStates, decays, temp},

           Utils`AssertWithMessage[exactNumberOfProducts === 2,
              "decays with " <> ToString@exactNumberOfProducts <>
                    " final particles are not currently supported."
           ];

           (* list with {S, S}, {V,V}, etc *)
           genericFinalStates = GetAllowedGenericFinalStates[particle, exactNumberOfProducts];

           (* @todo checks on colour and Lorentz structure *)
           isPossibleDecay[finalState_] := (IsPhysicalFinalState[finalState] &&
                                            IsElectricChargeConservingDecay[particle, finalState] &&
                                            IsColorInvariantDecay[particle, finalState] &&
                                            !FinalStateContainsInitialState[particle, finalState] &&
                                               !IsSelfDecay[particle, finalState]);
           concreteFinalStates = Join @@ (GetParticleCombinationsOfType[#, allowedFinalStateParticles, isPossibleDecay]& /@ genericFinalStates);
           concreteFinalStates = OrderFinalState[particle, #] & /@ concreteFinalStates;

           Print[""];
           Print["Creating amplitudes for ", particle, " decays..."];
           (*ParallelEvaluate[Get["/home/wojciech/HEP-software/mathematica/SARAH-4.14.3/SARAH.m"]]
           ParallelEvaluate[Start["MSSM"]];*)
           decays =
              Map[
                 (
                    (* @todo StringPadRigh was introduced only in 10.1 *)
                    WriteString["stdout", StringPadRight["   - Creating amplitude for " <> ToString@particle <> " -> " <> ToString@#, 64, "."]];
                    temp = FSParticleDecay[particle, #, GetContributingGraphsForDecay[particle, #]];
                    Print[" Done."];
                    temp
                 )&,
                 concreteFinalStates(*,
                 DistributedContexts -> All*)
              ];

           decays = Select[decays, GetDecayTopologiesAndDiagrams[#] =!= {}&];
           DeleteDiagramsVanishingDueToColor /@ decays
    ];

GetDecaysForParticle[particle_, n_, allowedFinalStateParticles_List] :=
    (
     Print["Error: invalid number of final state particles: ", n];
     Quit[1];
    );

GatherParticlesByType[particles_List] :=
    Module[{areSameType},
           areSameType[p1_, p2_] := Or @@ ((#[p1] && #[p2])& /@ { TreeMasses`IsScalar,
                                                                  TreeMasses`IsVector,
                                                                  TreeMasses`IsFermion,
                                                                  TreeMasses`IsGhost });
           Gather[particles, areSameType]
          ];

(* returns a list of lists of the form {{particle type 1, {particles of type 1}}, {particle type 2, {particles of type 2}}, ...} *)
GetClassifiedParticleLists[particles_List] :=
    Module[{classified, foundTypes},
           classified = {GetGenericTypeName[First[#]], #}& /@ GatherParticlesByType[particles];
           foundTypes = First /@ classified;
           If[Length[Union[foundTypes]] != Length[foundTypes],
              Print["Error: particles incorrectly classified: ", classified];
              Quit[1];
             ];
           classified
          ];

BaseMulticombination[k_] := Table[1, {i, 1, k}];

(* see, e.g., http://www.martinbroadhurst.com/multicombinations.html *)
NextMulticombination[n_, combination_] :=
    Module[{k = Length[combination], incrementable, pos, val},
           incrementable = Position[combination, x_ /; x < n];
           If[Length[incrementable] == 0,
              {},
              pos = First[Last[incrementable]];
              val = combination[[pos]] + 1;
              Join[Take[combination, pos - 1], {val}, Table[val, {i, 1, k - pos}]]
             ]
          ];

NextMulticombinationsList[setSizes_List, combinations_List] :=
    Module[{numCombinations = Length[combinations], next},
           next = combinations;
           For[i = numCombinations, i > 0, i--,
               nextCombination = NextMulticombination[setSizes[[i]], combinations[[i]]];
               If[nextCombination =!= {},
                  next[[i]] = nextCombination;
                  Return[next];,
                  next[[i]] = BaseMulticombination[Length[next[[i]]]];
                 ];
              ];
           {}
          ];

GetParticleCombinationsOfType[genericState_List, particles_List, isValidTest_:Function[True]] :=
    Module[{genericTypeCounts, classifiedParticles, indexLists, candidate, combinations},
           genericTypeCounts = {#, Count[genericState, #]}& /@ DeleteDuplicates[genericState];
           classifiedParticles = Select[GetClassifiedParticleLists[particles], MemberQ[genericState, First[#]]&];
           genericTypeCounts = genericTypeCounts[[Flatten[Position[First /@ classifiedParticles, First[#]]& /@ genericTypeCounts]]];
           indexLists = BaseMulticombination[Last[#]]& /@ genericTypeCounts;
           combinations = Reap[
               While[indexLists =!= {},
                     candidate = Flatten[MapIndexed[With[{pos = First[#2], indices = #1},
                                                    Last[classifiedParticles[[pos]]][[#]]& /@ indices]&, indexLists]];
                     If[isValidTest[candidate],
                        Sow[candidate];
                       ];
                     indexLists = NextMulticombinationsList[Length[Last[#]]& /@ classifiedParticles, indexLists];
                    ];
               ];
           Flatten[Last[combinations], 1]
          ];

GetAllowedGenericFinalStates[particle_ /; (TreeMasses`IsScalar[particle] || TreeMasses`IsVector[particle]),
                             n_Integer] :=
    Switch[n,
           2, {{GenericScalarName[], GenericScalarName[]},
               {GenericScalarName[], GenericVectorName[]},
               {GenericVectorName[], GenericVectorName[]},
               {GenericFermionName[], GenericFermionName[]}},
           _, Print["Error: cannot determine allowed generic final states for n = ", n]; Quit[1];
          ];

GetAllowedGenericFinalStates[particle_?TreeMasses`IsFermion, n_Integer] :=
    Switch[n,
           2, {{GenericScalarName[], GenericFermionName[]},
               {GenericVectorName[], GenericFermionName[]}},
           _, Print["Error: cannot determine allowed generic final states for n = ", n]; Quit[1];
          ];

GetVerticesForDecay[decay_FSParticleDecay] :=
    Module[{diagrams = GetDecayDiagramsOnly[decay]},
           DeleteDuplicates[Join @@ (Flatten[(CXXDiagrams`VerticesForDiagram /@ Last[#]), 1]& /@ diagrams)]
          ];

GetVerticesForDecays[particleDecays_List] :=
    DeleteDuplicates[Flatten[GetVerticesForDecay /@ particleDecays, 1]];

LoopOverIndex[loopBody_String, index_, start_, stop_, type_:CConversion`ScalarType[CConversion`integerScalarCType]] :=
    Module[{idxStr, startStr, stopStr},
           idxStr = CConversion`ToValidCSymbolString[index];
           startStr = CConversion`ToValidCSymbolString[start];
           stopStr = CConversion`ToValidCSymbolString[stop];
           "for (" <> CConversion`CreateCType[type] <> " " <> idxStr <>
           " = " <> startStr <> "; " <> idxStr <> " < " <> stopStr <> "; ++" <>
           idxStr <> ") {\n" <> TextFormatting`IndentText[loopBody] <> "\n}"
          ];

(* generates a loop over the given indices, in the form
   {{idx1, start, stop}, {idx2, start, stop}, ...} with
   the first list entry being the innermost loop *)
LoopOverIndexCollection[loopBody_String, indices_List] :=
    "\n" <> Fold[LoopOverIndex[#1, Sequence @@ #2]&, loopBody, indices] <> "\n";

CreateGenericGetPartialWidthFunctionName[] := "get_partial_width";

CreateSpecializedPartialWidthCalculationName[initialState_, finalState_List, fieldsNamespace_] :=
    CreateGenericGetPartialWidthFunctionName[] <> "<" <>
    CXXDiagrams`CXXNameOfField[initialState, prefixNamespace -> fieldsNamespace] <> "," <>
    Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#, prefixNamespace -> fieldsNamespace]& /@ finalState, ","] <> " >";
CreateSpecializedPartialWidthCalculationName[initialState_, finalState_List] :=
    CreateGenericGetPartialWidthFunctionName[] <> "<" <>
        CXXDiagrams`CXXNameOfField[initialState] <> ", " <>
        Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#]& /@ finalState, ", "] <> ">";

CreatePartialWidthCalculationName[decay_FSParticleDecay, scope_:""] :=
    Module[{initialState, initialStateName,
            finalState, finalStateName},
           initialState = GetInitialState[decay];
           initialStateName = CConversion`ToValidCSymbolString[initialState];
           finalState = GetFinalState[decay];
           finalStateName = StringJoin[CConversion`ToValidCSymbolString /@ finalState];
           scope <> If[scope != "", "::", ""] <> "partial_width_" <> initialStateName <> "_to_" <> finalStateName
          ];

CreatePartialWidthCalculationPrototype[decay_FSParticleDecay] :=
    Module[{returnType = "", functionName = "", functionArgs = "",
            initialStateDim, finalStateDims},
           initialStateDim = TreeMasses`GetDimension[GetInitialState[decay]];
           finalStateDims = TreeMasses`GetDimension /@ GetFinalState[decay];
           functionArgs = FlexibleSUSY`FSModelName <> "_mass_eigenstates_interface*" <>
                          If[initialStateDim > 1, ", int", ""] <>
                          StringJoin[If[# > 1, ", int", ""]& /@ finalStateDims];
           functionName = CreatePartialWidthCalculationName[decay];
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           returnType <> " " <> functionName <> "(" <> functionArgs <> ") const;"
          ];

CreatePartialWidthCalculationFunction[decay_FSParticleDecay, fieldsNamespace_] :=
    Module[{returnType = "", functionName = "", functionArgs = "",
            initialState = GetInitialState[decay], initialStateDim,
            finalState = GetFinalState[decay], finalStateDims, setFieldIndices, body = ""},
           initialStateDim = TreeMasses`GetDimension[GetInitialState[decay]];
           finalStateDims = TreeMasses`GetDimension /@ GetFinalState[decay];
           functionArgs = FlexibleSUSY`FSModelName <> "_mass_eigenstates_interface* model" <>
                          If[initialStateDim > 1, ", int gI1", ""] <>
                          StringJoin[MapIndexed[If[#1 > 1, ", int gO" <> ToString[First[#2]], ""]&, finalStateDims]];
           functionName = CreatePartialWidthCalculationName[decay, "CLASSNAME"];
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           setFieldIndices[field_, indicesName_, indexVal_] :=
               Module[{dim, numIndices, result = ""},
                      dim = TreeMasses`GetDimension[field];
                      numIndices = CXXDiagrams`NumberOfFieldIndices[field];
                      result = "const typename field_indices<" <> CXXDiagrams`CXXNameOfField[field] <> ">::type " <> indicesName;
                      If[numIndices == 0 || dim <= 1,
                         result = result <> " {};\n";,
                         result = result <> " {{" <> ToString[indexVal] <>
                                  StringJoin[Table[", 0", {i, 1, numIndices - 1}]] <> "}};\n";
                        ];
                      result
                     ];
           body = "context_base context {model};\n" <>
                  StringJoin[setFieldIndices[#[[1]], #[[2]], #[[3]]]& /@
                                 Join[{{initialState, "in_indices", If[initialStateDim > 1, "gI1", ""]}},
                                      MapIndexed[{#1, "out_" <> ToString[First[#2]] <> "_indices",
                                                  If[finalStateDims[[First[#2]]] > 1, "gO" <> ToString[First[#2]], ""]}&, finalState]]];
           body = body <> "\nreturn " <> CreateSpecializedPartialWidthCalculationName[initialState, finalState] <>
                  "(context, in_indices" <> StringJoin[Table[", out_" <> ToString[i] <> "_indices", {i, 1, Length[finalState]}]] <> ");\n";
           returnType <> " " <> functionName <> "(" <> functionArgs <> ") const\n{\n" <>
                  TextFormatting`IndentText[body] <> "}\n"
          ];

CreatePartialWidthCalculationPrototypes[particleDecays_List] :=
    Module[{allDecays},
           allDecays = Flatten[Last /@ particleDecays];
           Utils`StringJoinWithSeparator[CreatePartialWidthCalculationPrototype /@ allDecays, "\n"]
          ];

CreatePartialWidthCalculationFunctions[particleDecays_List, fieldsNamespace_] :=
    Module[{allDecays},
           allDecays = Flatten[Last /@ particleDecays];
           Utils`StringJoinWithSeparator[CreatePartialWidthCalculationFunction[#, fieldsNamespace]& /@ allDecays, "\n"]
          ];

CallPDGCodeGetter[SARAH`bar[particle_], args__] :=
    "-" <> CallPDGCodeGetter[particle, args];

CallPDGCodeGetter[Susyno`LieGroups`conj[particle_], args__] :=
   "-" <> CallPDGCodeGetter[particle, args];

CallPDGCodeGetter[particle_, idx_String, namespace_] :=
    Module[{dim = TreeMasses`GetDimension[particle], particleStr, result = ""},
           particleStr = namespace <> If[namespace != "", "::", ""] <> CConversion`ToValidCSymbolString[particle];
           result = namespace <> If[namespace != "", "::", ""] <> "get_pdg_code_for_particle(" <>
                    particleStr;
           If[dim > 1,
              result = result <> ", " <> idx;
             ];
           result <> ")"
          ];

CallPartialWidthCalculation[decay_FSParticleDecay] :=
    Module[{initialState = GetInitialState[decay], initialStateDim,
            finalState = GetFinalState[decay], finalStateDims, functionArgs = "",
            pdgsList = "", loopIndices, body = "", finalStateStarts, skip},
           initialStateDim = TreeMasses`GetDimension[initialState];
           finalStateDims = TreeMasses`GetDimension /@ finalState;
           finalStateStarts = TreeMasses`GetDimensionStartSkippingGoldstones /@ finalState;
           functionArgs = "dec_model.get()" <> If[initialStateDim > 1, ", gI1", ""] <>
                          MapIndexed[If[#1 > 1, ", gO" <> ToString[First[#2]], ""]&, finalStateDims];
           pdgsList = MapIndexed[With[{idx = First[#2]},
                                      CallPDGCodeGetter[#1, If[finalStateDims[[idx]] > 1, "gO" <> ToString[idx], ""], FlexibleSUSY`FSModelName <> "_info"]]&,
                                 finalState];
           pdgsList = "{" <> Utils`StringJoinWithSeparator[pdgsList, ", "] <> "}";

           CheckOffShellDecay[x_, y_] := (
           initialState === x && (
                 finalState === {y, SARAH`AntiField[y]} ||
                     finalState === {SARAH`AntiField[y], y}

               ));
           (* skip decays of a particle into itself *)
           skip =
              MapIndexed[
                 With[{idx = First[#2]},
                    If[
                       initialState === #1, "if (gI1 == gO" <> ToString[idx] <> ") {\n" <> TextFormatting`IndentText["continue;\n}\n"], ""
                    ]]&, finalState
              ] <>
              If[
                CheckOffShellDecay[TreeMasses`GetHiggsBoson[], TreeMasses`GetWBoson[]] || CheckOffShellDecay[TreeMasses`GetHiggsBoson[], TreeMasses`GetZBoson[]], "",

              "if (context.physical_mass<" <> CXXNameOfField[initialState] <> ">(std::array<int, " <> If[initialStateDim > 1, "1", "0"] <> ">{" <> If[initialStateDim > 1, "gI1", ""] <> "}) < " <>
               StringJoin @ Riffle[
                  MapIndexed[
                     With[{idx = First[#2]},
                        "context.physical_mass<" <> CXXNameOfField[#1] <> ">(" <>
                                      "std::array<int, " <> If[finalStateDims[[idx]] > 1, "1", "0"] <> "> {" <> If[finalStateDims[[idx]] > 1, "gO" <> ToString[idx], ""] <> "})"
                     ]&,
                     finalState
                  ], " + "
              ] <> ") {\n" <> TextFormatting`IndentText["continue;\n"] <> "}\n"
                ];
           (* call decay *)
                      body = "decays.set_decay(" <> CreatePartialWidthCalculationName[decay] <> "(" <> functionArgs <> "), " <> pdgsList <>
               ", create_process_string<" <> CXXNameOfField[initialState] <> ", " <> StringRiffle[CXXNameOfField/@finalState, ", "] <> ">(" <>
                  If[initialStateDim > 1, "{gI1}", "{}"] <> "," <>
                  StringJoin @ Riffle[
                  MapIndexed[
                     With[{idx = First[#2]},
                        If[finalStateDims[[idx]] > 1, "{gO" <> ToString[idx] <> "}", "{}"]
                     ]&,
                     finalState
                  ], ","] <>
                  ")" <>
           ");\n";

           loopIndices = Reverse[Select[MapIndexed[With[{idx = First[#2]},
                                                        If[#1 > 1,
                                                           {"gO" <> ToString[idx], Utils`MathIndexToCPP@finalStateStarts[[idx]], #1},
                                                           {}
                                                          ]
                                                       ]&, finalStateDims], (# =!= {})&]];
           If[loopIndices =!= {},
              body = LoopOverIndexCollection[skip <> body, loopIndices],
             body = If[
               CheckOffShellDecay[TreeMasses`GetHiggsBoson[], TreeMasses`GetWBoson[]] || CheckOffShellDecay[TreeMasses`GetHiggsBoson[], TreeMasses`GetZBoson[]],
               "",



               "\nif (context.physical_mass<" <> CXXNameOfField[initialState] <>
                 ">(std::array<int, " <> If[initialStateDim > 1, "1", "0"] <> ">{" <> If[initialStateDim > 1, "gI1", ""] <> "}) > " <>
                 StringJoin @ Riffle[
                    MapIndexed[
                       With[{idx = First[#2]},
                          "context.physical_mass<" <> CXXNameOfField[#1] <> ">(" <>
                             "std::array<int, " <> If[finalStateDims[[idx]] > 1, "1", "0"] <> "> {" <> If[finalStateDims[[idx]] > 1, "gO" <> ToString[idx], ""] <> "})"
                       ]&,
                       finalState
                    ], " + "
                 ] <> ") {\n"] <> TextFormatting`IndentText[body] <>
                 If[CheckOffShellDecay[TreeMasses`GetHiggsBoson[], TreeMasses`GetWBoson[]] || CheckOffShellDecay[TreeMasses`GetHiggsBoson[], TreeMasses`GetZBoson[]],
                   "",

                 "}\n"
                   ]
             ]
          ];

CreateDecaysCalculationFunctionName[particle_, scope_:""] :=
    scope <> If[scope =!= "", "::", ""] <> "calculate_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysCalculationPrototype[particle_] :=
    "void " <> CreateDecaysCalculationFunctionName[particle] <> "();";

CreateDecaysCalculationPrototypes[decayLists_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationPrototype[First[#]]& /@ decayLists, "\n"];

CreateLocalScope[body_] := "{\n" <> TextFormatting`IndentText[body] <> "}\n";

CreateDecaysCalculationFunction[decaysList_] :=
    Module[{particle = First[decaysList], particleDim, particleStart,
            decayChannels = Last[decaysList],
            runToScale = "", body = ""},

           particleDim = TreeMasses`GetDimension[particle];
           particleStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle];

           runToScale =
              "auto decay_mass = PHYSICAL(" <>
                 CConversion`ToValidCSymbolString[TreeMasses`GetMass[particle]] <> ");\n" <>
                 "if (decay_mass" <> If[particleDim > 1, "(gI1)", ""] <> " > qedqcd.displayPoleMZ()) {\n" <>
                 TextFormatting`IndentText[
                    "model.run_to(decay_mass" <> If[particleDim > 1, "(gI1)", ""] <>  ");\n" <>
                    "model.solve_ewsb();\n"
                 ] <> "}\n";

           body = StringJoin[CallPartialWidthCalculation /@ decayChannels];
           body = "\nauto& decays = decay_table.get_" <> CConversion`ToValidCSymbolString[particle] <>
                  "_decays(" <> If[particleDim > 1, "gI1", ""] <> ");\n" <> body;
           body =
              "const int flag = 1;\n" <>
              "auto dec_model = [&] () -> std::unique_ptr<" <> FlexibleSUSY`FSModelName <> "_mass_eigenstates_interface> {\n" <>
              TextFormatting`IndentText[
                 "switch (flag) {\n" <> TextFormatting`IndentText[
                 "case 1: {\n" <>
                 TextFormatting`IndentText[
                    "auto dm = std::make_unique<" <> FlexibleSUSY`FSModelName <> "_mass_eigenstates_decoupling_scheme>(model.get_input());\n" <>
                    "// fill_from BSM model has to be called before fill_from SM\n" <>
                    "// both calls are required\n" <>
                    "dm->fill_from(model);\n" <>
                    "standard_model::Standard_model sm{};\n" <>
                    "sm.initialise_from_input(qedqcd);\n" <>
                    "// set loop level for RGE running to match RGE setting\n" <>
                    "// of BSM model\n" <>
                    "sm.set_loops(model.get_loops());\n" <>
                    "if (run_to_decay_particle_scale) {\n" <>
                    TextFormatting`IndentText[
                       "auto decay_mass = PHYSICAL(" <>
                          CConversion`ToValidCSymbolString[TreeMasses`GetMass[particle]] <> ");\n" <>
                       "if (decay_mass" <> If[particleDim > 1, "(gI1)", ""] <> " > qedqcd.displayPoleMZ()) {\n" <>
                       TextFormatting`IndentText[
                          "sm.run_to(decay_mass" <> If[particleDim > 1, "(gI1)", ""] <>  ");\n"
                       ] <> "}\n"
                    ] <>
                    "}\n" <>
                    "sm.solve_ewsb_tree_level();\n" <>
                    "sm.calculate_DRbar_masses();\n" <>
                    "dm->fill_from(sm);\n" <>
                    "return dm;\n" <>
                    "break;\n"
                 ] <>
                 "}\n" <>
                 "case 2:\n" <>
                 TextFormatting`IndentText[
                    "return std::make_unique<" <> FlexibleSUSY`FSModelName <> "_mass_eigenstates>(model);\n" <>
                    "break;\n"
                  ] <>
                  "default:\n" <>
                  TextFormatting`IndentText[
                     "throw SetupError(\"flag value is not supported\");\n"
                 ]
                 ] <>
                 "}\n"
                 ] <>
                 "}();\n\n" <>
              "context_base context {dec_model.get()};\n" <>
              body;

           body = "\nif (run_to_decay_particle_scale) {\n" <>
                  TextFormatting`IndentText[runToScale] <> "}\n\n" <> body;
           If[particleDim > 1,
              body = LoopOverIndexCollection[body, {{"gI1", Utils`MathIndexToCPP@particleStart, particleDim}}] <> "\n";
             ];
           "void " <> CreateDecaysCalculationFunctionName[particle, "CLASSNAME"] <>
           "()\n{\n"
           <> TextFormatting`IndentText[body] <> "\n}\n"
          ];

CreateDecaysCalculationFunctions[particleDecays_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationFunction /@ particleDecays, "\n"];

CallDecaysFunction[particle_, arg_:"model", obj_:""] :=
    obj <> CreateDecaysCalculationFunctionName[particle] <> "();\n";

CallThreadedDecaysFunction[particle_, ptr_:"this", pool_:"tp"] :=
    pool <> ".run_task([" <> ptr <> "] () { " <>
    If[ptr === "this", "", ptr <> "->"] <>
    CreateDecaysCalculationFunctionName[particle] <> "(); });\n";

CallDecaysCalculationFunctions[particles_List, enableDecaysThreads_] :=
    Module[{result = ""},
           If[enableDecaysThreads,
              result = "Thread_pool tp(std::min(std::thread::hardware_concurrency(), " <>
                       ToString[Length[particles]] <> "u));\n\n" <>
                       StringJoin[CallThreadedDecaysFunction /@ particles];
              ,
              result = StringJoin[CallDecaysFunction /@ particles];
             ];
           result
          ];

GetDecayAmplitudeType[initialParticle_?TreeMasses`IsScalar, finalState_List] :=
    Module[{vertexType},
           vertexType = CXXDiagrams`VertexTypeForFields[Join[{initialParticle}, finalState]];
           Switch[vertexType,
                  ScalarVertex, "Decay_amplitude_SSS",
                  ChiralVertex, "Decay_amplitude_SFF",
                  MomentumDifferenceVertex, "Decay_amplitude_SSV",
                  InverseMetricVertex, "Decay_amplitude_SVV",
                  _, Print["Warning: decay ", initialParticle, " -> ", finalState, " is not supported."];
                     "Unknown_amplitude_type"
                 ]
          ];

GetDecayAmplitudeType[initialParticle_?TreeMasses`IsFermion, finalState_List] :=
    Module[{vertexType},
           vertexType = CXXDiagrams`VertexTypeForFields[Join[{initialParticle}, finalState]];
           Switch[{vertexType, GetFeynArtsTypeName@Last@finalState},
                  {ChiralVertex, S}, "Decay_amplitude_FFS",
                  {ChiralVertex, V}, "Decay_amplitude_FFV",
                  _, Print["Warning: decay ", initialParticle, " -> ", finalState, " is not supported."];
                     "Unknown_amplitude_type"
                 ]
          ];

GetDecayAmplitudeType[decay_FSParticleDecay] :=
    GetDecayAmplitudeType[GetInitialState[decay], GetFinalState[decay]];

CreateFieldIndices[particle_String] :=
    "typename " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::field_indices<" <> particle <> " >::type";

CreateFieldIndices[particle_, fieldsNamespace_] :=
    CreateFieldIndices[CXXDiagrams`CXXNameOfField[particle, prefixNamespace -> fieldsNamespace]];

CreateTotalAmplitudeFunctionName[] := "calculate_amplitude";

CreateTotalAmplitudeSpecializationDecl[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            returnType = "", fieldsNamespace, fieldsList, templatePars = "", args = ""},
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           fieldsNamespace = FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields";
           fieldsList = Join[{initialParticle}, finalState];
           templatePars = "<" <> Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#, prefixNamespace -> modelName <> "_cxx_diagrams::fields"]& /@
                                                               fieldsList, ", "] <> ">";
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[#, fieldsNamespace] <> "&")& /@ fieldsList, ", "];
           "template<>\n" <> returnType <> " " <> modelName <> "_decays::" <>
           CreateTotalAmplitudeFunctionName[] <> templatePars <> "(" <> args <> ") const;\n"
          ];

FillSSSDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{assignments = ""},
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay]] <>
                         ">(idx_1);\n";
           assignments = assignments <> structName <> ".m_scalar_1 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[First[GetFinalState[decay]]] <>
                         ">(idx_2);\n";
           assignments = assignments <> structName <> ".m_scalar_2 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[Last[GetFinalState[decay]]] <>
                         ">(idx_3);\n";
           assignments
          ];

FillSFFDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{assignments = ""},
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay]] <>
                         ">(idx_1);\n";
           assignments = assignments <> structName <> ".m_fermion_1 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[First[GetFinalState[decay]]] <>
                         ">(idx_2);\n";
           assignments = assignments <> structName <> ".m_fermion_2 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[Last[GetFinalState[decay]]] <>
                         ">(idx_3);\n";
           assignments
          ];

FillSSVDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{finalState, scalar, scalarPos, vector, vectorPos, assignments = ""},
           finalState = GetFinalState[decay];
           scalar = First[Select[finalState, TreeMasses`IsScalar]];
           scalarPos = First[First[Position[finalState, scalar]]];
           vector = First[Select[finalState, TreeMasses`IsVector]];
           vectorPos = First[First[Position[finalState, vector]]];
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay]] <>
                         ">(idx_1);\n";
           assignments = assignments <> structName <> ".m_scalar = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[scalar] <>
                         ">(idx_" <> ToString[scalarPos + 1] <> ");\n";
           assignments = assignments <> structName <> ".m_vector = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[vector] <>
                         ">(idx_" <> ToString[vectorPos + 1] <> ");\n";
           assignments
          ];

(* replace 'physical_mass' with 'mass' to check Ward identity *)
FillSVVDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{assignments = ""},
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay]] <>
                         ">(idx_1);\n";
           assignments = assignments <> structName <> ".m_vector_1 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[First[GetFinalState[decay]]] <>
                         ">(idx_2);\n";
           assignments = assignments <> structName <> ".m_vector_2 = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[Last[GetFinalState[decay]]] <>
                         ">(idx_3);\n";
           assignments
          ];

FillFFSDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{finalState, scalar, scalarPos, fermion, fermionPos, assignments = ""},
           finalState = GetFinalState[decay];
           fermion = First[Select[finalState, TreeMasses`IsFermion]];
           fermionPos = First[First[Position[finalState, fermion]]];
           scalar = First[Select[finalState, TreeMasses`IsScalar]];
           scalarPos = First[First[Position[finalState, scalar]]];
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay]] <>
                         ">(idx_1);\n";
           assignments = assignments <> structName <> ".m_fermion = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[fermion] <>
                         ">(idx_" <> ToString[fermionPos + 1] <> ");\n";
           assignments = assignments <> structName <> ".m_scalar = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[scalar] <>
                         ">(idx_" <> ToString[scalarPos + 1] <> ");\n";
           assignments
          ];

FillFFVDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{finalState, fermion, fermionPos, vector, vectorPos, assignments = ""},
           finalState = GetFinalState[decay];
           fermion = First[Select[finalState, TreeMasses`IsFermion]];
           fermionPos = First[First[Position[finalState, fermion]]];
           vector = First[Select[finalState, TreeMasses`IsVector]];
           vectorPos = First[First[Position[finalState, vector]]];
           assignments = assignments <> structName <> ".m_decay = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[GetInitialState[decay]] <>
                         ">(idx_1);\n";
           assignments = assignments <> structName <> ".m_fermion = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[fermion] <>
                         ">(idx_" <> ToString[fermionPos + 1] <> ");\n";
           assignments = assignments <> structName <> ".m_vector = " <> paramsStruct <> ".physical_mass<" <>
                         CXXDiagrams`CXXNameOfField[vector] <>
                         ">(idx_" <> ToString[vectorPos + 1] <> ");\n";
           assignments
          ];

FillDecayAmplitudeMasses[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Switch[GetDecayAmplitudeType[decay],
           "Decay_amplitude_SSS", FillSSSDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SFF", FillSFFDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SSV", FillSSVDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SVV", FillSVVDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFS", FillFFSDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFV", FillFFVDecayAmplitudeMasses[decay, modelName, structName, paramsStruct],
           _, ""
          ];
FillMasses[decay_FSParticleDecay] :=
   Switch[GetDecayAmplitudeType[decay],
      "Decay_amplitude_SSS", "result.m_decay, result.m_scalar_1, result.m_scalar_2",
      "Decay_amplitude_SFF", "result.m_decay, result.m_fermion_1, result.m_fermion_2",
      "Decay_amplitude_SSV", "result.m_decay, result.m_scalar, result.m_vector",
      "Decay_amplitude_SVV", "result.m_decay, result.m_vector_1, result.m_vector_2",
      "Decay_amplitude_FFS", "result.m_decay, result.m_fermion, result.m_scalar",
      "Decay_amplitude_FFV", "result.m_decay, result.m_fermion, result.m_vector",
      _, ""
   ];

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SSS",
                              structName_] :=
    structName <> ".form_factor = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SSV",
                              structName_] :=
    structName <> ".form_factor = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SVV",
                              structName_] :=
    structName <> ".form_factor_g = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_11 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_12 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_21 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_22 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_eps = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_SFF",
                              structName_] :=
    structName <> ".form_factor_left = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_right = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_FFS",
                              structName_] :=
    structName <> ".form_factor_left = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_right = std::complex<double>(0., 0.);\n";

ZeroDecayAmplitudeFormFactors[decay_FSParticleDecay /; GetDecayAmplitudeType[decay] == "Decay_amplitude_FFV",
                              structName_] :=
    structName <> ".form_factor_gam_left = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_gam_right = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_p_1 = std::complex<double>(0., 0.);\n" <>
    structName <> ".form_factor_p_2 = std::complex<double>(0., 0.);\n";

GetTreeLevelTwoBodyDecayVertex[decay_FSParticleDecay] :=
    Module[{treeLevelDiags, vertices = {}},
           treeLevelDiags = GetDecayDiagramsOnlyAtLoopOrder[decay, 0];
           If[treeLevelDiags =!= {},
              vertices = Flatten[CXXDiagrams`VerticesForDiagram /@ treeLevelDiags, 1];
             ];
           vertices
          ];

EvaluateTreeLevelTwoBodyDecayVertex[decay_FSParticleDecay, modelName_, indicesName_, paramsStruct_, resultName_:"vertex"] :=
    Module[{vertexFields, templatePars, sortedVertexFields},
           vertexFields = GetTreeLevelTwoBodyDecayVertex[decay];
           If[Length[vertexFields] > 1,
              Print["Error: more than a single vertex in tree-level decays."];
              Quit[1];
             ];
           If[vertexFields =!= {},
              vertexFields = First[vertexFields];
              sortedVertexFields = SortFieldsInCp[vertexFields];
              templatePars = "<" <>
                              Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField[#]&
                                                            /@ sortedVertexFields, ", "] <> " >";
              "const auto " <> resultName <> " =  Vertex" <> templatePars <> "::evaluate(" <>
              indicesName <> ", " <> paramsStruct <> ");\n",
              ""
             ]
          ];

FillSSSTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor += vertex.value();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillSFFTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           sortedFieldsList = SortFieldsInCp[fieldsList];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[
                        Permute[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], FindPermutation[fieldsList, sortedFieldsList]],
                        ", "
                     ] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_left += vertex.left();\n" <>
                         structName <> ".form_factor_right += vertex.right();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillSSVTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor += vertex.value(0, 1);\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillSVVTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsList, fieldsNamespace, indices, vertex, assignments, sortedFieldsList},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           sortedFieldsList = SortFieldsInCp[fieldsList];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[
                        Permute[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], FindPermutation[fieldsList, sortedFieldsList]],
                        ", "
                     ] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_g += vertex.value();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillFFSTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsList, fieldsNamespace, indices, vertex, assignments, sortedFieldsList},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           sortedFieldsList = SortFieldsInCp[fieldsList];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[
                        Permute[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], FindPermutation[fieldsList, sortedFieldsList]],
                        ", "
                     ] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_left += vertex.left();\n" <>
                         structName <> ".form_factor_right += vertex.right();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillFFVTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{fieldsList, fieldsNamespace, indices, vertex, assignments},
           fieldsList = Join[{GetInitialState[decay]}, GetFinalState[decay]];
           fieldsNamespace = modelName <> "_cxx_diagrams::fields";
           indices = "const auto indices = concatenate(" <>
                     Utils`StringJoinWithSeparator[Table["idx_" <> ToString[i], {i, 1, Length[fieldsList]}], ", "] <> ");\n";
           vertex = EvaluateTreeLevelTwoBodyDecayVertex[decay, modelName, "indices", paramsStruct];
           assignments = structName <> ".form_factor_gam_left += vertex.left();\n" <>
                         structName <> ".form_factor_gam_right += vertex.right();\n";
           "// tree-level amplitude\n" <> indices <> vertex <> "\n" <> assignments
          ];

FillTreeLevelDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Switch[GetDecayAmplitudeType[decay],
           "Decay_amplitude_SSS", FillSSSTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SFF", FillSFFTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SSV", FillSSVTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_SVV", FillSVVTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFS", FillFFSTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           "Decay_amplitude_FFV", FillFFVTreeLevelDecayAmplitudeFormFactors[decay, modelName, structName, paramsStruct],
           _, ""
          ];

(* the order of lists a and b might be different *)
Compare[aIn_List, bIn_List] := Module[{
  a = SortBy[aIn, First],
  b = SortBy[bIn, First]
},

  Utils`AssertWithMessage[Length[a] === Length[b], "Compare called with lists of unequal length " <> ToString@a <> ToString@b];
   For[i = 1, i <= Length[a], i++,
     If[(
       (Rule @@ MapAt[Sort, List @@ (a[[i]]), 1]) =!= (Rule @@ MapAt[Sort, List @@(b[[i]]), 1])), Return[False]]
   ];

   True
];

(* returns an element of the list from
   meta/generic_loop_decay_diagram_classes.m
   corresponding to a given insertion of fields *)
GenericTranslationForInsertion[topology_, insertion_] := Module[{
   genericDiagramsWithCorrectTopology,
   genericFieldTypesOnEdges,
   res
},

   If[$translationFile === $Failed,
      Print["Getting translation file failed"];Quit[1]
   ];

   genericDiagramsWithCorrectTopology = Select[$translationFile, MemberQ[#, topology]&];
   Utils`AssertWithMessage[genericDiagramsWithCorrectTopology =!= {},
      "Can't find topology " <> ToString@topology <> " for insertion " <> ToString@insertion
   ];

   genericFieldTypesOnEdges = (#[[1]] -> GetFeynArtsTypeName[#[[2]]])& /@ InsertionsOnEdgesForDiagram[topology, insertion];

   (* find the generic-level diagram in diagramsWithCorrectTopology that corresponds
      to the given generic insertion *)
   res =
      Select[
         genericDiagramsWithCorrectTopology,
         Compare[#[[3]]/.#[[4]], genericFieldTypesOnEdges]&
      ];

   (* there should be one-to-one correspondence *)
   Utils`AssertWithMessage[Length[res] === 1 || Plus@@Flatten[(Length/@Take[#,{-3}])& /@ res,1] === 0,
      "Couldn't identify the C++ translation for a needed diagram"
   ];

   First@res
];

(* ContractionsBetweenVerticesForDiagramFromGraph returns indices as array,
   even for external particles which are in diagrams not a member of array *)
FieldFromDiagram[diagram_, {i_, j_}] :=
   If[!ListQ[diagram[[i]]],
      diagram[[i]],
      diagram[[i,j]]
   ];

(* Returns a list of particles on edges of a diagram in form of list of elements as
   {vertex1, vertex2} -> concrete particle connecting vertex1 and 2 (e.g. Fe)
   This is the core function that allows matching between diagram at mathematica and C++
   level. *)
(* @todo: not clear about the conjugation of particle *)
InsertionsOnEdgesForDiagram[topology_, insertion_] :=
   Module[{sortedVertexCombinations, connectedVertices, res},

      sortedVertexCombinations =
         Select[
            Tuples[
               Range[CXXDiagrams`NumberOfExternalParticlesInTopology[topology] + CXXDiagrams`NumberOfPropagatorsInTopology[topology]],
               2
            ],
            OrderedQ
         ];
      (* remove not connected pairs *)
      connectedVertices = Select[sortedVertexCombinations,
         CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[Sequence@@#, insertion, topology] =!= {}&
      ];

      res = Flatten[
            Switch[Length[CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, insertion, topology]],
               (* vertices connected by 1 particle *)
               1, {{#1, #2} -> FieldFromDiagram[insertion, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, insertion, topology][[1,1]]}]},
               (* vertices connected by 2 particle *)
               2, {
                  {#1, #2} -> FieldFromDiagram[insertion, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, insertion, topology][[1,1]]}],
                  {#1, #2} -> FieldFromDiagram[insertion, {#1, CXXDiagrams`ContractionsBetweenVerticesForDiagramFromGraph[#1, #2, insertion, topology][[2,1]]}]
               },
               (* error *)
               _, Print["Only 1 or 2 connections between vertices are supported"]; Quit[1];
            ]& @@@ connectedVertices,
         1
      ];

      res
];

(* returns vertex needed by FA generated routines, but in terms of SortCp'ed vertex *)
ConvertCouplingToCPP[Decays`Private`FACp[particles__][lor_], fieldAssociation_, vertices_, indices_] :=
   Module[{vertexEdges, res, pos, globalMinus = 1, vertex, lorSorted, temp},

   vertexEdges = (List[particles] /. Index[Generic, n_] :> n);
   pos = First@First@Position[vertices, vertexEdges];
   lorSorted =
      lor /. Index[Generic, n_] :> n;
   (* we use field position in vertex in Mom and not field itself as there
      are vertices like {hh, hh, VZ} in the MSSMCPV which would give
      Mom[hh] - Mom[hh] which is 0 *)
   lorSorted = lorSorted /. Mom[-f_[n_]] :> Mom[First@First@Position[vertexEdges, -f[n]]] /. Mom[f_[n_]] :> Mom[First@First@Position[vertexEdges, f[n]]];
   (* SortCp requires that FFV vertex has only PL|PR as lorent structure,
      not LorentzProduct[gamma[lt], PL|PR] *)
   lorSorted = lorSorted /. LorentzProduct[gamma[_], lr_:PL|PR] :> lr;
   vertex =
      SARAH`Cp[Sequence@@(vertexEdges /. -f_[n_] /; !IsParticle[f] :> SARAH`AntiField[Field[n] /. fieldAssociation] /. f_[n_] /; !IsParticle[f] && f =!= Susyno`LieGroups`conj && f =!= SARAH`bar :> (Field[n] /. fieldAssociation))][lorSorted];
   temp =
      SortCp[vertex];
   temp = temp /. -SARAH`Cp[y___][x___] :> (globalMinus*=-1;{{y}, x}) /. SARAH`Cp[y___][x___] :> {{y}, x};
   vertex = First@temp;
   lorSorted = Last@temp;

   res =
      Replace[lorSorted, {
         (* pure only scalar vertices *)
         1 -> "value()",

         (* fermion vertices *)
         PL -> "left()",
         PR -> "right()",
         LorentzProduct[_, PL] -> "left()",
         LorentzProduct[_, PR] -> "right()",

         (* @todo: rules below need checking! *)
         (* momentum difference vertices *)
         Mom[i1_Integer] - Mom[i2_Integer] :>
            "value(" <> ToString[i1-1] <> ", " <> ToString[i2-1] <> ")",

        (* metric tensor vertices *)
        (* there's a global - sign in those expression, that's why permutation with even
           signature is replaced with odd_permutation *)
        g[lt1, lt2] (-Mom[i1_] + Mom[i2_])
            + g[lt1, lt3] (Mom[i1_] - Mom[i3_])
            + g[lt2, lt3] (-Mom[i2_] + Mom[i3_]) :> (
               Switch[Signature[{i1, i2, i3}],
                  1, "value(TripleVectorVertex::odd_permutation {})",
                 -1, "value(TripleVectorVertex::even_permutation {})",
                 _,  (Print["Can't find TripleVectorVertex permutation"]; Quit[1])
               ]

        ),

         g[l1_, l2_] g[l3_, l4_] :>
            Switch[{l1, l2, l3, l4},
               {lt1, lt2, lt3, lt4}, "value1()",
               {lt1, lt3, lt2, lt4}, "value2()",
               {lt1, lt4, lt2, lt3}, "value3()",
               _, (Print["Unknown lorentz index combination " <> ToString@{l1, l2, l3, l4} <> " in 4-vector vertex."]; Quit[1])
            ],

         g[lt1_, lt2_] -> "value()",

         (* momentum vertices *)
         (* apparently FeynArts writes all ghost-ghost-vector vertices as proportional to momentum
            of bared ghost. This is opposite to Sarah where all such vertices are written using
            the momentum of non-bared ghost *)
         Mom[i_Integer] :> (
            globalMinus *= -1;
            With[{unbarredGhosts = Select[Delete[vertex, i], (IsGhost[#] && Head[#]=!=SARAH`bar)&]},
               If[Head[vertex[[i]]]===SARAH`bar && Length@unbarredGhosts =!= 1,
                  Print["Error! Couldn't identify ghost in vertex"];
                  Quit[1],
                  "value(" <> ToString@Utils`MathIndexToCPP@FieldPositionInVertex[If[Head[vertex[[i]]]=!=SARAH`bar, vertex[[i]], First@unbarredGhosts], vertex] <> ")"
               ]
            ]
         ),

        (* error *)
        lorSorted :> (Print["Unhandled lorentz structure " <> ToString@lorSorted]; Quit[1])
      }
   ];

   If[globalMinus === -1, "-", ""] <>
   "1.0i*vertex" <> ToString@indices[[pos]] <> "Val." <> res
];

(*
   GetFieldsAssociations is the core function that matches concrete Feynman diagram generated
   by FS onto a generic diagram from FeynArts/FormCalc.

   concrete fields on edges, e.g.
   {{1, 4} -> hh, {2, 5} -> SRdp, {3, 6} -> VWm, {4, 5} -> SRdp, {4, 6} -> Hpm, {5, 6} -> VZ}
   come from FS

   fieldNumberOnEdgeBetweenVertices, e.g.
   {{1, 4} -> Field[1], {2, 5} -> Field[2], {3, 6} -> Field[3], {4, 5} -> Field[4], {4, 6} -> Field[5], {5, 6} -> Field[6]}
   come from FeynArts/FormCalc translation files

   fieldTypes, e. g.
   {Field[1] -> S, Field[2] -> S, Field[3] -> V, Field[4] -> S, Field[5] -> S, Field[6] -> V}

   GetFieldsAssociations returns translation in the form
   {Field[1] -> hh, Field[2] -> SRdp, Field[3] -> VWm, Field[4] -> SRdp, Field[5] -> Hpm, Field[6] -> VZ}

   @todo: There is an ambiguity whether Field[n] -> somefield or Field[n] -> AntiField[somefield]
*)
GetFieldsAssociations[concreteFieldOnEdgeBetweenVertices_, fieldNumberOnEdgeBetweenVertices_, fieldTypes_, diagram_, verticesInFieldTypesForFACp_] :=
   Module[{temp = {}, concreteFieldOnEdgeBetweenVerticesLocal, verticesForFACp, temp2},

      concreteFieldOnEdgeBetweenVerticesLocal = concreteFieldOnEdgeBetweenVertices;
      temp = SortBy[Reverse /@ fieldNumberOnEdgeBetweenVertices, Last];

      For[i = 1, i <= Length[temp], i++,

         If[GetFeynArtsTypeName[concreteFieldOnEdgeBetweenVerticesLocal[[i,2]]] === (temp[[i,1]] /. fieldTypes),

            temp[[i]] = temp[[i]] /. concreteFieldOnEdgeBetweenVerticesLocal[[i]];
            temp[[i]] = temp[[i]] /. (Reverse@concreteFieldOnEdgeBetweenVerticesLocal[[i,1]] -> SARAH`AntiField[concreteFieldOnEdgeBetweenVerticesLocal[[i,2]]]),

            concreteFieldOnEdgeBetweenVerticesLocal[[{i, i+1}]] = concreteFieldOnEdgeBetweenVerticesLocal[[{i+1, i}]];
            i=i-1;
          ]
      ];

      temp
];

(* map topology to FeynArts name *)
FeynArtsTopologyName[topology_] :=
   Switch[topology,
      {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 0,
        0, 1, 1}, {0, 1, 0, 1, 0, 1}, {0, 0, 1, 1, 1, 0}}, "T1",
      {{0, 0, 0, 1, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {1, 1, 0, 0,
        1}, {0, 0, 1, 1, 1}}, "T2",
      {{0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {0, 0, 0, 1, 0}, {1, 0, 1, 0,
        1}, {0, 1, 0, 1, 1}}, "T3",
      {{0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {0, 0, 0, 0, 1}, {1, 0, 0, 0,
        2}, {0, 1, 1, 2, 0}}, "T4",
      {{0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {0, 0, 0, 0, 1}, {1, 0, 0, 1,
        1}, {0, 1, 1, 1, 0}}, "T5",
      {{0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {0, 0, 0, 1, 0}, {1, 0, 1, 0,
        2}, {0, 1, 0, 2, 0}}, "T6",
      {{0, 0, 0, 1, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {1, 1, 0, 0,
        2}, {0, 0, 1, 2, 0}}, "T7",
      {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {1, 1, 0,
        0, 0, 1}, {0, 0, 1, 0, 0, 2}, {0, 0, 0, 1, 2, 0}}, "T8",
      {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 1, 0, 0}, {1, 0, 1,
        0, 0, 1}, {0, 1, 0, 0, 0, 2}, {0, 0, 0, 1, 2, 0}}, "T9",
      {{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 1, 0}, {1, 0, 0,
        0, 0, 2}, {0, 1, 1, 0, 0, 1}, {0, 0, 0, 2, 1, 0}}, "T10",
      _, Print["Error: Cannot map topology to FeynArts name"]; Quit[1]
   ];

WrapCodeInLoop[indices_, code_] :=
   (
      "\n// loops over vertices' indices\n" <>
      StringJoin@@(MapIndexed[
      Nest[
         TextFormatting`IndentText,
         "for (const auto& index" <> ToString@#1 <> ": index_range<vertex" <> ToString@#1 <> ">()) {\n" <>
            If[First@#2===Length[indices], "\n" <> TextFormatting`IndentText[code], ""],
         First@#2 -1
      ]&,
      indices
   ]) <>
     StringJoin@@ (Nest[
         TextFormatting`IndentText,
         "}\n",
         #-1
         ]& /@ Reverse@Range@Length[indices]));

WrapCodeInLoopOverInternalVertices[decay_, topology_, diagram_] :=
   Module[{(*vertices, *)indices, cppVertices,
      mass = {}, translation, fieldAssociation,
      externalEdges,
     (*verticesInFieldTypes, *)matchExternalFieldIndicesCode, matchInternalFieldIndicesCode = "", functionBody = "",
     verticesInFieldTypesForFACp, verticesForFACp, colorFac = "colorFac", symmetryFac = "symmetryFac", FinitePart = False,
     whereToConj, verticesForFACp2
   },

      translation = GenericTranslationForInsertion[topology, diagram];
      If[Length[translation]===7 && translation[[-3]] === {},
        Return[{FinitePart, ""}]
      ];

      (* the vertices positions FACp are not ordered according to numbers in element 2 & 3 of translation *)
      (* vertex in terms of field types (S,V,F,...) and indices 1, 2 *)
      (*      verticesInFieldTypes =*)
      (*          Select[translation[[-2]], ListQ] /. Susyno`LieGroups`conj[f_] :> -f /. Field[n_] :> (Field[n]/.translation[[4]])[n];*)
      verticesInFieldTypesForFACp =
         List @@@ (DeleteDuplicates[translation[[-3]] /. Index[Generic, n_Integer] :> n /. Decays`Private`FACp[x__][__] :> FACp[x]]);

      (* {Field[1] -> concrete field, ...} *)
      fieldAssociation =
          GetFieldsAssociations[
            InsertionsOnEdgesForDiagram[topology, diagram],
            translation[[3]],
            translation[[4]],
            diagram,
            verticesInFieldTypesForFACp
          ];

      (* vertices in an orientation as required by Cp *)
      verticesForFACp = verticesInFieldTypesForFACp /. (fieldAssociation /. ((#1 -> #2@@#1)& @@@ translation[[4]])) /. - e_ :> SARAH`AntiField[e];

      (* extra conjugations to bring verticesForFACp equal (possibly up to order) with vertices in diagram *)
      Module[
         {
            i=1, t,
            whereToConj =
               DeleteCases[
                  Subsets[
                     Select[
                        Range[4, Length[fieldAssociation]],
                        (fieldAssociation[[#,2]] =!= SARAH`AntiField[fieldAssociation[[#,2]]])&
                     ]
                  ],
                  {}
               ]
         },
         (* We do a double sort because the order of fields in vertices or the order of
            vertices doesn't have to be same. *)
         While[Sort[Sort /@ verticesForFACp] =!= Sort[Sort /@ Drop[diagram, 3]],
            If[i > Length@whereToConj,
               Print["Error! Could not determine field association list"];
               Quit[1]
            ];
            t = {#, 2}& /@ (whereToConj[[i]]);
            fieldAssociation = MapAt[SARAH`AntiField, fieldAssociation, If[Length[t]===1, First@t, t]];
            verticesForFACp = verticesInFieldTypesForFACp /. f_[n_Integer] :> Field[n] /. fieldAssociation /. - e_ :> SARAH`AntiField[e];
            i++;
         ]
      ];

      (* set of unique indices used in names of vertices and indices *)
      indices = Table["Id"<>ToString@i, {i, Length@verticesForFACp}];

      (* create using declarations for SortCp'ed vertices *)
      verticesForFACp2 = (SortCp[Cp[Sequence@@#]])& /@ verticesForFACp;
      verticesForFACp2 = If[Length[#] == 2, List@@Last@#, List@@#]& /@ verticesForFACp2;
      cppVertices =
         "using vertex" <> ToString@#1 <> " = Vertex<" <>
            (StringJoin@Riffle[CXXDiagrams`CXXNameOfField /@ #2  ,", "] <> ">;\n")& @@@ Transpose[{indices, verticesForFACp2}];

      (* List of {integer, integer} -> Field[integer] *)
      externalEdges =
         Select[
            translation[[3]],
            (MatchQ[#, ({i_Integer, j_Integer} -> Field[_Integer]) /; (First@Sort[{i,j}] >= 1 && First@Sort[{i,j}] <= 3 && Last@Sort[{i,j}] > 3)])&
         ];

      verticesInFieldTypesForFACpCpSorted = MapIndexed[Permute[#1, FindPermutation[verticesForFACp[[First@#2]], verticesForFACp2[[First@#2]]]]&, verticesInFieldTypesForFACp];
        matchExternalFieldIndicesCode = StringJoin@@(With[{
          positions = Position[verticesInFieldTypesForFACpCpSorted /. -f_[i_] :> f[i], (Field[#] /. translation[[4]])[#]]
        },
If[Length@positions =!= 1, Quit[1]];
                  "const auto externalFieldIndicesIn" <> ToString[#] <>
                  " = vertex" <>
                  ToString@indices[[ positions[[1,1]] ]] <>
               "::template indices_of_field<" <>
                  ToString@Utils`MathIndexToCPP@positions[[1,2]] <>
                  ">(index" <>
                  ToString@indices[[ positions[[1,1]]]] <> ");\n"
                ]& /@ Range[3]
        );

      (* NOTE! Position has a weird feature that
         Position[{{x}}, x] = {{1, 1}} and
         Position[{{-x}}, x] = {{1, 1, 2}},
         hence the -field_->field rule *)
      mass =
         Map[
            (
               "const double mInternal" <> ToString[#-3] <> " {context.mass<" <>
               CXXNameOfField[
                  verticesForFACp2[[  Sequence@@First@Position[verticesInFieldTypesForFACpCpSorted /. -field_ :> field, _[#]]    ]]
               ] <> ">(" <>
               "vertex" <> ToString@indices[[First@First@Position[verticesInFieldTypesForFACpCpSorted/. -field_:>field, _[#]]]] <> "::template indices_of_field<" <>
               ToString@Utils`MathIndexToCPP[
                  Last@First@Position[verticesInFieldTypesForFACpCpSorted /. -field_:>field, _[#]]
               ] <>
               ">(index" <> ToString@indices[[First@First@Position[verticesInFieldTypesForFACpCpSorted /.-field_:>field, _[#]]]] <> "))};\n"
            )&,
            Range[4, Length@fieldAssociation ]
         ];

      mass = StringJoin@@mass;

      With[{
        positions = Position[verticesInFieldTypesForFACpCpSorted/. -f_[n_] :> Field[n] /. f_[n_] :> Field[n], Field[#]]
         },
        If[Length@positions =!= 2, Quit[1]];
            matchInternalFieldIndicesCode = matchInternalFieldIndicesCode <>
                "if (vertex" <> ToString@indices[[positions[[1,1]]]] <> "::template indices_of_field<" <>

                ToString@Utils`MathIndexToCPP@positions[[1,2]] <>

                ">(index"<> ToString@indices[[positions[[1,1]]]] <> ")" <>

                " != " <>
                "vertex" <> ToString@indices[[positions[[2,1]]]] <> "::template indices_of_field<" <>
                ToString@Utils`MathIndexToCPP@positions[[2,2]] <>

                ">(index"<> ToString@indices[[positions[[2,1]]]] <> "))\n" <> TextFormatting`IndentText["continue;\n"] <> "\n"
      ]& /@ Range[4, Length@fieldAssociation];

      (* body of nested loops over vertices indices *)
      Block[{ampCall =
                  "result += " <> ToString@symmetryFac <> " * " <> ToString@colorFac <> " * calculate_" <> translation[[1]] <> "(\n" <>
                  TextFormatting`IndentText[
                     (* external masses *)
                     FillMasses[decay] <> ",\n" <>
                     (* internal masses *)
                     StringJoin @@ Riffle[("mInternal" <> ToString@#)& /@ Range@CXXDiagrams`NumberOfPropagatorsInTopology[topology], ", "] <> ", " <> "\n" <>
                     (* couplings *)
                     TextFormatting`WrapLines[
                        StringJoin @@ Riffle[ToString /@ ConvertCouplingToCPP[#, fieldAssociation, verticesInFieldTypesForFACp, indices]& /@ translation[[-3]], ", "] <> ",\n"
                     ] <>
                     (* renormalization scale *)
                     "ren_scale" <>
                     (* if amplitude is UV divergent, take the finite part *)
                     If[!Last@translation === True, FinitePart=True; ",\nFinite", ""] <> ")"
                  ],
            fieldsInLoop = DeleteDuplicates[Map[
            verticesForFACp2[[  Sequence@@First@Position[verticesInFieldTypesForFACpCpSorted /. -field_ :> field, _[#]]    ]]&,
      Range[4, 3+Length@verticesForFACp]
      ] /. SARAH`bar|Susyno`LieGroups`conj -> Identity]

               },
         functionBody =
                  "// skip indices that don't match external indices\n" <>
                  matchExternalFieldIndicesCode <>
                     "\nif (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)\n" <>
                     TextFormatting`IndentText["continue;"] <>
                     "\n\n" <>

                  "// connect internal particles in vertices\n" <>
                  matchInternalFieldIndicesCode <>

                     StringJoin[(
                        "auto const vertex" <> # <> "Val = " <>
                        "vertex" <> # <> "::evaluate(index" <>
                           #  <> ", context);\n"
                     )& /@ indices] <>
                  "\nif (" <> StringJoin@Riffle[("!vertex" <> # <> "Val.isZero()")& /@ indices, " && "] <> ") {\n" <>
                  TextFormatting`IndentText[
                  "// internal masses\n" <>
                  mass <> "\n" <>
                  (* in some cases, we apply higher order corrections at the level of amplitude *)
                  If[
                     (* for H/A -> gamma gamma *)
                     (GetHiggsBoson[] === First@diagram || GetPseudoscalarHiggsBoson[] === First@diagram) && And @@ (TreeMasses`IsPhoton /@ Take[diagram, {2,3}]),
                     If[(* CP conserving Higgs sector *)
                        !SA`CPViolationHiggsSector &&
                        (* the quark loop amplitude *)
                        Length[fieldsInLoop] === 1 && ContainsAll[TreeMasses`GetSMQuarks[], fieldsInLoop],
                        "if (include_higher_order_corrections == SM_higher_order_corrections::enable &&\n" <>
                        TextFormatting`IndentText[
                           Module[{pos1, post2, res},
                              StringJoin@Riffle[
                              MapIndexed[
                              (pos1 = Position[#1, First@fieldsInLoop, 1];
                              pos2 = Position[#1, SARAH`bar[First@fieldsInLoop], 1];
                              If[MatchQ[pos1, {{_Integer}}] && MatchQ[pos2, {{_Integer}}],
                                 "vertexId" <> ToString@First@#2 <> "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP@First@First@pos1 <> ">(indexId" <> ToString@First@#2 <> ") == " <>
                                 "vertexId" <> ToString@First@#2 <> "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP@First@First@pos2 <> ">(indexId" <> ToString@First@#2 <> ")"
                              ])&, verticesForFACp]," &&\n"]
                           ] <>
                           ") {\n"
                        ] <>
                        TextFormatting`IndentText[
                           ampCall <> " * (1. + get_alphas(context)/Pi * delta_" <> Switch[First@diagram, GetHiggsBoson[], "h", GetPseudoscalarHiggsBoson[], "Ah"] <> "AA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale));\n"
                        ] <> "}\n" <>
                        "else {\n" <>
                        TextFormatting`IndentText[
                           ampCall <> ";\n"
                        ] <> "}\n",
                    If[!SA`CPViolationHiggsSector && Length[fieldsInLoop] === 1 && And@@Join[TreeMasses`IsScalar /@ fieldsInLoop, TreeMasses`ColorChargedQ /@ fieldsInLoop],
                     "\nif (include_higher_order_corrections == SM_higher_order_corrections::enable &&\n" <>
                     TextFormatting`IndentText[
                       Module[{pos1, post2, res},
                          StringJoin@Riffle[
                          MapIndexed[
                             (pos1 = Position[#1, First@fieldsInLoop, 1];
                          pos2 = Position[#1, Susyno`LieGroups`conj[First@fieldsInLoop], 1];
                          If[MatchQ[pos1, {{_Integer}}] && MatchQ[pos2, {{_Integer}}],
                              "vertexId" <> ToString@First@#2 <> "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP@First@First@pos1 <> ">(indexId" <> ToString@First@#2 <> ") == " <>
                              "vertexId" <> ToString@First@#2 <> "::template indices_of_field<" <> ToString@Utils`MathIndexToCPP@First@First@pos2 <> ">(indexId" <> ToString@First@#2 <> ")"
                          ])&, verticesForFACp]
                          , " &&\n"
                       ]
                       ] <>
                     ") {\n"] <>
                     TextFormatting`IndentText[
                        ampCall <> " * (1. + get_alphas(context)/Pi * " <>
                        ToString@N[Switch[TreeMasses`GetColorRepresentation@First@fieldsInLoop,
                           T, 4/3,
                           -T, 4/3,
                           O, 3,
                           _, Print["Error! Unknown color charge of scalar in 2-loop QCD corrections to H->gamma gamma"];Quit[1]
                         ], 16] <> " * delta_" <> Switch[First@diagram, GetHiggsBoson[], "h", GetPseudoscalarHiggsBoson[], "Ah"] <> "AA_2loopQCD_for_squark_loop(result.m_decay, mInternal1, ren_scale));\n"
                     ] <> "}\n" <>
                    "else {\n" <>
                    TextFormatting`IndentText[
                        ampCall <> ";\n"
                    ] <> "}\n", ampCall <> ";\n"
                    ]
                  ], ampCall <> ";\n"
                  ]] <> "}\n"
     ];

      {FinitePart,

         (* diagram information *)
         "\n// topology " <> FeynArtsTopologyName[topology] <>
            "\n// internal particles in the diagram: " <>  StringJoin[Riffle[ToString@Part[#, 2]& /@Drop[fieldAssociation, 3], ", "]] <> "\n" <>
            "{\n" <>
               TextFormatting`IndentText[

         (* usings for vertices *)
         cppVertices <>

         (* diagram symmetry factor *)
          "\nconstexpr double " <> ToString@symmetryFac <> " {" <>
             ToString @
               N[With[{topoName = FeynArtsTopologyName[topology]},

                (* weird FA factor *)
                If[MemberQ[{"T2", "T3", "T4", "T5", "T8", "T9", "T10"}, topoName], 2, 1] *
                If[topoName === "T9" || topoName === "T8",
                   If[(Field[5] /. fieldAssociation) === (SARAH`AntiField[Field[6] /. fieldAssociation]), 1, 1/2],
                   1
                ]
               ],16] <>
         "};\n" <>

         (* color factor *)
         "\nconstexpr " <>
         With[{cf = CXXDiagrams`ExtractColourFactor @
                CXXDiagrams`ColorFactorForDiagram[topology, diagram]},
            If[Head[cf] === Complex,
               (* complex colour factor *)
               "std::complex<double> " <> ToString@colorFac <> " " <>
                  ToString[N[#, 16]& /@ FSReIm @ cf],

               (* real color factor *)
               "double " <> ToString@colorFac <> " {" <>
                  ToString @ N[cf, 16] <> "}"
            ]
         ] <> ";\n" <>
         WrapCodeInLoop[indices, functionBody]
               ] <>
         "}\n"
      }
   ];

FillOneLoopDecayAmplitudeFormFactors[decay_FSParticleDecay, modelName_, structName_, paramsStruct_] :=
    Module[{oneLoopTopAndInsertion, body = "", FinitePart = False},

       (* list of elements like {topology, insertion (diagram)} *)
       oneLoopTopAndInsertion = Flatten[With[{topo = #[[1]], diags = #[[2]]}, {topo, #}& /@ diags]& /@ GetDecayTopologiesAndDiagramsAtLoopOrder[decay, 1], 1];

       With[{ret = WrapCodeInLoopOverInternalVertices[decay, Sequence @@ #]},
         FinitePart = (ret[[1]] || FinitePart);
         body = body <> ret[[2]];
       ]& /@ oneLoopTopAndInsertion;

       {FinitePart,
         "\n// ----------------- 1-loop contributions to the amplitude -----------------\n" <>
             body
       }
    ];

(* creates `calculate_amplitude` function
   that returns a total (sumed over internal insertions) 1-loop amplitude for a given external particles *)
CreateTotalAmplitudeSpecializationDef[decay_FSParticleDecay, modelName_] :=
   Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsNamespace = "fields",
            returnVar = "result", paramsStruct = "context", returnType = "",
            externalFieldsList, templatePars = "", args = "",
            body = ""},

           (* @todo StringPadRigh was introduced only in 10.1 *)
           (*WriteString["stdout", StringPadRight["   - Creating code for " <> ToString@initialParticle <> " -> " <> ToString@finalState, 64, "."]];*)

           (* decay amplitude type, e.g. Decay_amplitude_FSS *)
           returnType = GetDecayAmplitudeType[decay];

           (* template arguments *)
           externalFieldsList = Join[{initialParticle}, finalState];
           templatePars = "<" <> Utils`StringJoinWithSeparator[
              CXXDiagrams`CXXNameOfField[#]& /@ externalFieldsList, ", "] <> ">";

           (* function arguments *)
           args =
              "const context_base& " <> paramsStruct <> ",\n" <>
                 Utils`StringJoinWithSeparator[
                     MapIndexed[(CreateFieldIndices[#1, fieldsNamespace] <> " const& idx_" <> ToString[First[#2]])&, externalFieldsList],
                     ",\n"
                 ];

           (* function body *)

           body = "\n// amplitude type\n";
           body = body <> returnType <> " " <> returnVar <> ";\n";

           body = body <> "\n// external particles' masses\n";
           body = body <> FillDecayAmplitudeMasses[decay, modelName, returnVar, paramsStruct];

           (* @todo: it might be less verbose to by default initialize amplitude to 0 *)
           body = body <> "\n// set the initial value of an amplitude to 0\n";
           body = body <> ZeroDecayAmplitudeFormFactors[decay, returnVar];

           If[IsPossibleTreeLevelDecay[decay, True],
              body = body <> "// @todo correct prefactors\n" <> FillTreeLevelDecayAmplitudeFormFactors[decay, modelName, returnVar, paramsStruct] <> "\n";
             ];

           If[!IsPossibleTreeLevelDecay[decay, True] && IsPossibleOneLoopDecay[decay],
             With[{res = FillOneLoopDecayAmplitudeFormFactors[decay, modelName, returnVar, paramsStruct]},
                If[res[[1]],
                  body = body <> "\n// FormCalc's Finite variable\n";
                  body = body <>"constexpr double Finite {1.};\n"
                ];
                body = body <>"\nconst double ren_scale {result.m_decay};\n";
                body = body <> Last@res <> "\n";
             ]
             ];

           body = body <> "return " <> returnVar <> ";\n";

             "// " <> ToString@initialParticle <> " -> " <> ToString@finalState <> "\n" <>
                 "template<>\n" <>
                 returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
                 templatePars <>
                 "(\n" <> TextFormatting`IndentText[args <> ") const{\n"] <>
                 TextFormatting`IndentText[body] <>
                 "}\n"
   ];

GetHiggsBosonDecays[particleDecays_List] :=
    If[TreeMasses`GetHiggsBoson[] =!= Null, Select[particleDecays, (First[#] === TreeMasses`GetHiggsBoson[])&], {}];
GetPseudoscalarHiggsBosonDecays[particleDecays_List] :=
    If[TreeMasses`GetPseudoscalarHiggsBoson[] =!= Null, Select[particleDecays, (First[#] === TreeMasses`GetPseudoscalarHiggsBoson[])&], {}];

SelectDecayByFinalState[finalState_List, decays_List] :=
    Select[decays, (Sort[GetFinalState[#]] === Sort[finalState])&];

SelectDownQuarkDownQuarkFinalState[decays_List] :=
    Module[{downQuarkSymbol, result = {}},
           downQuarkSymbol = TreeMasses`GetDownQuark[1] /. field_[generation_] :> field;
           If[downQuarkSymbol =!= Null,
              result = SelectDecayByFinalState[{downQuarkSymbol, SARAH`AntiField[downQuarkSymbol]}, decays];
             ];
           result
          ];

SelectChargedLeptonChargedLeptonFinalState[decays_List] :=
    Module[{chargedLeptonSymbol, result = {}},
       chargedLeptonSymbol = TreeMasses`GetDownLepton[1] /. field_[generation_] :> field;
       If[chargedLeptonSymbol =!= Null,
          result = SelectDecayByFinalState[{chargedLeptonSymbol, SARAH`AntiField[chargedLeptonSymbol]}, decays];
       ];
       result
    ];

SelectGluonGluonFinalState[decays_List] :=
    Module[{gluonSymbol = TreeMasses`GetGluon[], result = {}},
           If[gluonSymbol =!= Null,
              result = SelectDecayByFinalState[{gluonSymbol, gluonSymbol}, decays];
             ];
           result
          ];

SelectHiggsHiggsFinalState[decays_List] :=
    Module[{higgsSymbol = TreeMasses`GetHiggsBoson[], result = {}},
           If[higgsSymbol =!= Null,
              result = SelectDecayByFinalState[{higgsSymbol, higgsSymbol}, decays];
             ];
           result
          ];

SelectPhotonPhotonFinalState[decays_List] :=
    Module[{photonSymbol = TreeMasses`GetPhoton[], result = {}},
           If[photonSymbol =!= Null,
              result = SelectDecayByFinalState[{photonSymbol, photonSymbol}, decays];
             ];
           result
          ];

SelectHiggsHiggsFinalState[decays_List] :=
    Module[{psSymbol = TreeMasses`GetPseudoscalarHiggsBoson[], result = {}},
           If[psSymbol =!= Null,
              result = SelectDecayByFinalState[{psSymbol, psSymbol}, decays];
             ];
           result
          ];

SelectUpQuarkUpQuarkFinalState[decays_List] :=
    Module[{upQuarkSymbol, result = {}},
           upQuarkSymbol = TreeMasses`GetUpQuark[1] /. field_[generation_] :> field;
           If[upQuarkSymbol =!= Null,
              result = SelectDecayByFinalState[{upQuarkSymbol, SARAH`AntiField[upQuarkSymbol]}, decays];
             ];
           result
          ];

SelectWWFinalState[decays_List] :=
    Module[{wBosonSymbol = TreeMasses`GetWBoson[], result = {}},
           If[wBosonSymbol =!= Null,
              result = SelectDecayByFinalState[{wBosonSymbol, SARAH`AntiField[wBosonSymbol]}, decays];
             ];
           result
          ];

SelectPhotonZFinalState[decays_List] :=
    Module[{zBosonSymbol = TreeMasses`GetZBoson[], photonSymbol = TreeMasses`GetPhoton[], result = {}},
           If[zBosonSymbol =!= Null && photonSymbol =!= Null,
              result = SelectDecayByFinalState[{zBosonSymbol, photonSymbol}, decays];
             ];
           result
          ];

SelectZZFinalState[decays_List] :=
    Module[{zBosonSymbol = TreeMasses`GetZBoson[], result = {}},
           If[zBosonSymbol =!= Null,
              result = SelectDecayByFinalState[{zBosonSymbol, zBosonSymbol}, decays];
             ];
           result
          ];

CreateHiggsToGluonGluonTotalAmplitudeFunction[hggDecay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[hggDecay], finalState = GetFinalState[hggDecay],
            fieldsList, returnType = "", args = "", templatePars = "", body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#]] <> "&")& /@ fieldsList, ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = returnType <> " result;\nreturn result;\n";
           "template<>\n" <> returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

CreateHiggsToGluonGluonTotalAmplitude[particleDecays_List, modelName_] :=
    Module[{higgsDecays, hggDecay, prototype = "", function = ""},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              hggDecay = SelectGluonGluonFinalState[Last[higgsDecays]];
              If[hggDecay =!= {},
                 hggDecay = First[hggDecay];
                 prototype = CreateTotalAmplitudeSpecializationDecl[hggDecay, modelName];
                 function  = CreateHiggsToGluonGluonTotalAmplitudeFunction[hggDecay, modelName]
                ];
             ];
           {prototype, function}
          ];

CreateHiggsToPhotonPhotonTotalAmplitudeFunction[hgamgamDecay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[hgamgamDecay], finalState = GetFinalState[hgamgamDecay],
            fieldsList, returnType = "", args = "", templatePars = "", body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#]] <> "&")& /@ fieldsList, ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = returnType <> " result;\nreturn result;\n";
           "template<>\n" <> returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

CreateHiggsToPhotonPhotonTotalAmplitude[particleDecays_List, modelName_] :=
    Module[{higgsDecays, hgamgamDecay, prototype = "", function = ""},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              hgamgamDecay = SelectPhotonPhotonFinalState[Last[higgsDecays]];
              If[hgamgamDecay =!= {},
                 hgamgamDecay = First[hgamgamDecay];
                 prototype = CreateTotalAmplitudeSpecializationDecl[hgamgamDecay, modelName];
                 function  = CreateHiggsToPhotonPhotonTotalAmplitudeFunction[hgamgamDecay, modelName]
                ];
             ];
           {prototype, function}
          ];

CreateHiggsToPhotonZTotalAmplitudeFunction[hgamgamDecay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[hgamgamDecay], finalState = GetFinalState[hgamgamDecay],
            fieldsList, returnType = "", args = "", templatePars = "", body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           returnType = GetDecayAmplitudeType[initialParticle, finalState];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#]] <> "&")& /@ fieldsList, ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = returnType <> " result;\nreturn result;\n";
           "template<>\n" <> returnType <> " CLASSNAME::" <> CreateTotalAmplitudeFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

CreateHiggsToPhotonZTotalAmplitude[particleDecays_List, modelName_] :=
    Module[{higgsDecays, hgamgamDecay, prototype = "", function = ""},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              hgamgamDecay = SelectPhotonPhotonFinalState[Last[higgsDecays]];
              If[hgamgamDecay =!= {},
                 hgamgamDecay = First[hgamgamDecay];
                 prototype = CreateTotalAmplitudeSpecializationDecl[hgamgamDecay, modelName];
                 function  = CreateHiggsToPhotonZTotalAmplitudeFunction[hgamgamDecay, modelName]
                ];
             ];
           {prototype, function}
          ];

CreateTotalAmplitudeSpecialization[decay_FSParticleDecay, modelName_] :=
    Module[{decl = "", def = "", vertices = {}},
           decl = CreateTotalAmplitudeSpecializationDecl[decay, modelName];
           def = CreateTotalAmplitudeSpecializationDef[decay, modelName];
           {decl, def}
          ];

CreateTotalAmplitudeSpecializations[particleDecays_List, modelName_] :=
    Module[{specializations, vertices = {}, listing = {},
            contextsToDistribute = {"SARAH`", "Susyno`LieGroups`", "FlexibleSUSY`", "CConversion`", "Himalaya`"}},
           Print[""];
           FSFancyLine[];
           Print["Creating a C++ code for decay amplitudes..."];
           If[FlexibleSUSY`FSEnableParallelism,
              ParallelEvaluate[(BeginPackage[#];EndPackage[];)& /@ contextsToDistribute, DistributedContexts->All];
              specializations =
                 AbsoluteTiming@ParallelMap[
                    (
                       (*If[!MemberQ[listing, GetInitialState[#]],
                          Print[""];
                          Print["Creating C++ code for ", GetInitialState[#], " decays..."];
                          AppendTo[listing, GetInitialState[#]];
                       ];*)
                       CreateTotalAmplitudeSpecialization[#, modelName]
                    )&,
                    Flatten[Last @@@ particleDecays, 1],
                    DistributedContexts -> All, Method -> "FinestGrained"
                 ],
              specializations =
                 AbsoluteTiming@Map[
                    (
                       If[!MemberQ[listing, GetInitialState[#]],
                          Print["Creating C++ code for ", GetInitialState[#], " decays..."];
                          AppendTo[listing, GetInitialState[#]];
                       ];
                       CreateTotalAmplitudeSpecialization[#, modelName]
                    )&,
                    Flatten[Last @@@ particleDecays, 1]
                 ]
           ];

           Print["The creation of C++ code for decays took ", Round[First@specializations, 0.1], "s"];
           specializations = Last@specializations;
           specializations = Select[specializations, (# =!= {} && # =!= {"", ""})&];
           Utils`StringJoinWithSeparator[#, "\n"]& /@ Transpose[specializations]
          ];

CreatePartialWidthSpecializationDecl[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, fieldsNamespace, args},
           fieldsList = Join[{initialParticle}, finalState];
           fieldsNamespace = If[modelName != "", modelName <> "_cxx_diagrams::fields", False];
           args = "const " <> modelName <> "_cxx_diagrams::context_base&, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[#, fieldsNamespace] <> "&")& /@ fieldsList, ", "];
           "template <>\n" <>
           "double " <> modelName <> "_decays::" <>
           CreateSpecializedPartialWidthCalculationName[initialParticle, finalState, fieldsNamespace] <>
           "(" <> args <> ") const;"
          ];

CreateIncludedPartialWidthSpecialization[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            declaration = "", includeStatement = ""},
           declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
           includeStatement = "#include \"decays/H_SM_decays/decay_" <>
                              SimplifiedName[initialParticle] <> "_to_" <>
                              StringJoin[SimplifiedName[# /. SARAH`bar|Susyno`LieGroups`conj -> Identity]& /@ finalState] <>
                              ".cpp\"";
           {declaration, includeStatement}
          ];

CreateHiggsToZZPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectZZFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToWWPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectWWFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

(*
CreateHiggsToGluonGluonPartialWidthFunction[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, args, templatePars, body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           fieldIndices = {"in_idx", "out1_idx", "out2_idx"};
           args = "const " <> modelName <> "_cxx_diagrams::context_base& context, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#1]] <> "& " <> #2)& @@@ Transpose[{fieldsList, fieldIndices}], ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body =
              "const auto amp = calculate_amplitude<H, G, G>(context, in_idx, out1_idx, out2_idx);\n" <>
              "const double mH = context.physical_mass<H>(in_idx);\n" <>
              "constexpr double ps {1./(8.*Pi)};\n" <>
              "constexpr double ps_symmetry {1./2.};\n" <>
              "constexpr double color_gen_sqr {8.0};\n" <>
              "const double flux = 0.5/mH;\n" <>
              "return  flux * color_gen_sqr * ps * ps_symmetry * amp.square();\n";
           "template <>\ndouble CLASSNAME::" <> CreateGenericGetPartialWidthFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateHiggsToPhotonPhotonPartialWidthFunction[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, args, templatePars, body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           fieldIndices = {"in_idx", "out1_idx", "out2_idx"};
           args = "const " <> modelName <> "_cxx_diagrams::context_base& context, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#1]] <> "& " <> #2)& @@@ Transpose[{fieldsList, fieldIndices}], ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body = 
              "const auto amp = calculate_amplitude<H, A, A>(context, in_idx, out1_idx, out2_idx);\n" <>
               "const double mH = context.physical_mass<H>(in_idx);\n" <>
               "constexpr double ps {1./(8.*Pi)};\n" <>
               "constexpr double ps_symmetry {1./2.};\n" <>
               "constexpr double color_gen_sqr {1.0};\n" <>
               "const double flux = 0.5/mH;\n" <>
               "return  flux * color_gen_sqr * ps * ps_symmetry * amp.square();\n";
           "template <>\ndouble CLASSNAME::" <> CreateGenericGetPartialWidthFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateHiggsToPhotonZPartialWidthFunction[decay_FSParticleDecay, modelName_] :=
    Module[{initialParticle = GetInitialState[decay], finalState = GetFinalState[decay],
            fieldsList, args, templatePars, body = ""},
           fieldsList = Join[{initialParticle}, finalState];
           fieldIndices = {"in_idx", "out1_idx", "out2_idx"};
           args = "const " <> modelName <> "_cxx_diagrams::context_base& context, " <>
                  Utils`StringJoinWithSeparator[("const " <> CreateFieldIndices[SimplifiedName[#1]] <> "& " <> #2)& @@@ Transpose[{fieldsList, fieldIndices}], ", "];
           templatePars = "<" <> Utils`StringJoinWithSeparator[SimplifiedName /@ fieldsList, ", "] <> ">";
           body =
              "const auto amp = calculate_amplitude<H, A, Z>(context, in_idx, out2_idx, out1_idx);\n" <>
              "const double mH = context.physical_mass<H>(in_idx);\n" <>
              "constexpr double ps {1./(8.*Pi) * beta(mIn, mOut1, mOut2)};\n" <>
              "constexpr double ps_symmetry {1.};\n" <>
              "constexpr double color_gen_sqr {1.0};\n" <>
              "const double flux = 0.5/mH;\n" <>
              "return  flux * color_gen_sqr * ps * ps_symmetry * amp.square();\n";
           "template <>\ndouble CLASSNAME::" <> CreateGenericGetPartialWidthFunctionName[] <>
           templatePars <> "(" <> args <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];
*)

CreateHiggsToGluonGluonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectGluonGluonFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreatePseudoscalarHiggsToGluonGluonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectGluonGluonFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToPhotonPhotonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectPhotonPhotonFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
              function = CreateHiggsToPhotonPhotonPartialWidthFunction[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToPhotonZPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, orderZPhotonFinalState, declaration = "", function = ""},
           decay = SelectPhotonZFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              orderZPhotonFinalState[finalState_] :=
                  Module[{zBosonSymbol = TreeMasses`GetZBoson[], photonSymbol = TreeMasses`GetPhoton[]},
                         If[finalState === {zBosonSymbol, photonSymbol},
                            finalState,
                            Reverse[finalState]
                           ]
                        ];
              decay = ReplacePart[decay, {{2}} -> orderZPhotonFinalState[GetFinalState[decay]]];
              declaration = CreatePartialWidthSpecializationDecl[decay, modelName];
              function = CreateHiggsToPhotonZPartialWidthFunction[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToHiggsHiggsPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectHiggsHiggsFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToPseudoscalarPseudoscalarPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectPseudoscalarPseudoscalarFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToUpQuarkUpQuarkPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectUpQuarkUpQuarkFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];
CreatePseudoscalarHiggsToUpQuarkUpQuarkPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectUpQuarkUpQuarkFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToDownQuarkDownQuarkPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectDownQuarkDownQuarkFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];
CreatePseudoscalarHiggsToDownQuarkDownQuarkPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
           decay = SelectDownQuarkDownQuarkFinalState[decaysList];
           If[decay =!= {},
              decay = First[decay];
              {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
             ];
           {declaration, function}
          ];

CreateHiggsToChargedLeptonChargedLeptonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
       decay = SelectChargedLeptonChargedLeptonFinalState[decaysList];
       If[decay =!= {},
          decay = First[decay];
          {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
       ];
       {declaration, function}
    ];
CreatePseudoscalarHiggsToChargedLeptonChargedLeptonPartialWidth[{higgsSymbol_, decaysList_}, modelName_] :=
    Module[{decay, declaration = "", function = ""},
       decay = SelectChargedLeptonChargedLeptonFinalState[decaysList];
       If[decay =!= {},
          decay = First[decay];
          {declaration, function} = CreateIncludedPartialWidthSpecialization[decay, modelName];
       ];
       {declaration, function}
    ];

CreateHiggsDecayPartialWidthSpecializations[particleDecays_, modelName_] :=
    Module[{higgsDecays, specializations = {}},
           higgsDecays = GetHiggsBosonDecays[particleDecays];
           If[higgsDecays =!= {},
              higgsDecays = First[higgsDecays];
              specializations = {CreateHiggsToZZPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToWWPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToGluonGluonPartialWidth[higgsDecays, modelName],
(*                                 CreateHiggsToPhotonPhotonPartialWidth[higgsDecays, modelName],*)
(*                                 CreateHiggsToPhotonZPartialWidth[higgsDecays, modelName],*)
                                 CreateHiggsToHiggsHiggsPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToUpQuarkUpQuarkPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToDownQuarkDownQuarkPartialWidth[higgsDecays, modelName],
                                 CreateHiggsToChargedLeptonChargedLeptonPartialWidth[higgsDecays, modelName]};
             ];
           specializations
          ];

CreatePseudoscalarHiggsDecayPartialWidthSpecializations[particleDecays_, modelName_] :=
    Module[{pseudoscalarHiggsDecays, specializations = {}},
           pseudoscalarHiggsDecays = GetPseudoscalarHiggsBosonDecays[particleDecays];
           If[pseudoscalarHiggsDecays =!= {},
              pseudoscalarHiggsDecays = First[pseudoscalarHiggsDecays];
              specializations =
                 {
                     CreatePseudoscalarHiggsToDownQuarkDownQuarkPartialWidth[pseudoscalarHiggsDecays, modelName],
                     CreatePseudoscalarHiggsToUpQuarkUpQuarkPartialWidth[pseudoscalarHiggsDecays, modelName],
                     CreatePseudoscalarHiggsToGluonGluonPartialWidth[pseudoscalarHiggsDecays, modelName],
                     CreatePseudoscalarHiggsToChargedLeptonChargedLeptonPartialWidth[pseudoscalarHiggsDecays, modelName]
                 }
              ];
           specializations
          ];

CreatePartialWidthSpecializations[particleDecays_List, modelName_] :=
    Module[{specializations},
           specializations = CreateHiggsDecayPartialWidthSpecializations[particleDecays, modelName];
           specializations = Join[specializations, CreatePseudoscalarHiggsDecayPartialWidthSpecializations[particleDecays, modelName]];
           specializations = Select[specializations, (# =!= {} && # =!= {"", ""})&];
           Utils`StringJoinWithSeparator[#, "\n"]& /@ Transpose[specializations]
          ];

CreateDecaysGetterFunctionName[particle_] :=
    "get_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysGetterFunction[particle_] :=
    Module[{dim, body = ""},
           dim = TreeMasses`GetDimension[particle];
           body = "return decay_table." <> CreateDecayTableEntryGetterName[particle] <>
                  "(" <> If[dim > 1, "i", ""] <> ");";
           "const Decays_list& " <> CreateDecaysGetterFunctionName[particle] <> "(" <>
           If[dim > 1, "int i", ""] <> ") const { " <> body <> " }"
          ];

CreateDecaysGetterFunctions[particles_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysGetterFunction /@ particles, "\n"];

CreateDecayTableEntryGetterName[particle_] :=
    "get_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecayTableEntryGetterPrototype[particle_] :=
    Module[{dim},
           dim = TreeMasses`GetDimension[particle];
           "Decays_list& " <> CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int", ""] <> ");"
          ];

CreateDecayTableEntryConstGetterPrototype[particle_] :=
    Module[{dim},
           dim = TreeMasses`GetDimension[particle];
           "const Decays_list& " <> CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int", ""] <> ") const;"
          ];

CreateDecayTableEntryGetterFunctionBody[particle_, rows_List] :=
    Module[{i, dimWithoutGoldstones, dimWithGoldstones, idxName = "gI1", errMsg = "", body = ""},
           dimWithoutGoldstones = TreeMasses`GetDimensionWithoutGoldstones[particle];
           dimWithGoldstones = TreeMasses`GetDimension[particle];
           If[dimWithoutGoldstones != Length[rows],
              Print["Error: number of rows (", Length[rows], ") does not match size of"];
              Print["    ", particle, " multiplet."];
              Quit[1];
             ];
           If[dimWithGoldstones == 1,
              body = "return decay_table[" <> ToString[rows[[1]]] <> "];\n";
              ,
              body = "switch (" <> idxName <> ") {\n";
              For[i = dimWithGoldstones - dimWithoutGoldstones, i < dimWithGoldstones, i++,
                  body = body <> "case " <> ToString[i] <> ": return decay_table[" <> ToString[rows[[i+1 - (dimWithGoldstones - dimWithoutGoldstones)]]] <> "]; break;\n";
                 ];
              body = body <> "}\n\n";
              errMsg = "std::ostringstream sstr;\n" <>
                       "sstr << \"invalid particle index \" << std::to_string(" <> idxName <> ") << '\\n';\n\n" <>
                       "throw OutOfBoundsError(sstr.str());\n";
              body = body <> errMsg;
             ];
           body
          ];

CreateDecayTableEntryGetterFunction[particle_, rows_List, scope_:"CLASSNAME"] :=
    Module[{dim, body},
           dim = TreeMasses`GetDimension[particle];
           body = CreateDecayTableEntryGetterFunctionBody[particle, rows];
           "Decays_list& " <> scope <> If[scope != "", "::", ""] <>
           CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int gI1", ""] <> ")\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateDecayTableEntryConstGetterFunction[particle_, rows_List, scope_:"CLASSNAME"] :=
    Module[{dim, body},
           dim = TreeMasses`GetDimension[particle];
           body = CreateDecayTableEntryGetterFunctionBody[particle, rows];
           "const Decays_list& " <> scope <> If[scope != "", "::", ""] <>
           CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int gI1", ""] <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateDecayTableEntryGetterFunction[particle_, row_Integer, scope_:"CLASSNAME"] :=
    CreateDecayTableEntryGetterFunction[particle, {row}, scope];

CreateDecayTableEntryConstGetterFunction[particle_, row_Integer, scope_:"CLASSNAME"] :=
    CreateDecayTableEntryConstGetterFunction[particle, {row}, scope];

CreateDecayTableGetterPrototypes[decayParticles_List] :=
    Utils`StringJoinWithSeparator[(CreateDecayTableEntryGetterPrototype[#] <> "\n" <>
                                   CreateDecayTableEntryConstGetterPrototype[#])& /@ decayParticles, "\n"];

CreateDecayTableGetterFunctions[decayParticles_List, scope_:"CLASSNAME"] :=
    Module[{dims, offsets, rowAssignments, defs = ""},
           dims = TreeMasses`GetDimensionWithoutGoldstones /@ decayParticles;
           offsets = If[Length[dims] == 1, {0}, Join[{0}, Accumulate[dims[[1;;-1]]]]];
           rowAssignments = MapIndexed[{decayParticles[[First[#2]]], Table[offsets[[First[#2]]] + i, {i, 0, #1 - 1}]}&, dims];
           defs = (CreateDecayTableEntryGetterFunction[#[[1]], #[[2]], scope] <> "\n\n" <>
                   CreateDecayTableEntryConstGetterFunction[#[[1]], #[[2]], scope])& /@ rowAssignments;
           Utils`StringJoinWithSeparator[defs, "\n"]
          ];

CreateDecayTableInitialization[decayParticles_List] :=
    Module[{dims, dimsWithoutGoldstones, starts, pdgCodes, initializerList = ""},
           dims = TreeMasses`GetDimension /@ decayParticles;
           dimsWithoutGoldstones = TreeMasses`GetDimensionWithoutGoldstones /@ decayParticles;
           starts = TreeMasses`GetDimensionStartSkippingGoldstones /@ decayParticles;
           pdgCodes = Parameters`GetPDGCodesForParticle /@ decayParticles;
           For[i = 1, i <= Length[decayParticles], i++,
               If[dims[[i]] != Length[pdgCodes[[i]]],
                  Print["Error: number of PDG codes does not match size of ", decayParticles[[i]], " multiplet."];
                  Quit[1];
                 ];
               If[dimsWithoutGoldstones[[i]] > 0,
                  initializerList = initializerList <> If[initializerList == "", "", ", "] <>
                                    Utils`StringJoinWithSeparator[("Decays_list(" <> ToString[#] <> ")")& /@ pdgCodes[[i, starts[[i]] ;;]], ", "];
                 ];
              ];
           ": decay_table({" <> initializerList <> "})"
          ];

End[];

EndPackage[];
