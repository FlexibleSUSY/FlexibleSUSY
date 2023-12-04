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

BeginPackage@"NPointFunctions`";
Begin@"`Private`";

AllTopologies[{2, 2}] = {
   treeS -> {1,0,1,0,0,1,0,1,0,1,0},
   treeT -> {1,0,0,1,1,0,0,1,0,1,0},
   treeU -> {1,0,0,1,0,1,1,0,0,1,0},
   treeAll -> {treeS, treeT, treeU},
   triangleT -> {1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0},
   inSelfT -> {1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,2,0,1,0,0,1,0},
   outSelfT -> {1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,2,0},
   penguinT -> {triangleT, inSelfT, outSelfT},
   boxS -> {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,1,0,1,0},
   boxT -> {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0},
   boxU -> {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,0},
   boxAll -> {boxS, boxT, boxU}
};

AdjaceTopology::usage = "Archivates the adjacency matrix of a given topology:
In[1]:= << FeynArts`;
        topologies = CreateTopologies[0, 1 -> 2];
        AdjaceTopology[topologies[[1]]]
Out[1]= {1, 1, 1, 0}";
AdjaceTopology[topology:_?IsTopology] :=
Module[{external, simple, graph, ordering, matrix},
   external = Tr@$externalFieldNumbers;
   simple = (#/. h_[i_, j_, _] :> h[i, j])&/@ topology;
   graph = List@@ (UndirectedEdge@@@ FeynArts`TopologySort@simple)/.
      FeynArts`Vertex[_][i_] :> i;
   ordering = Ordering@VertexList@graph;
   matrix = Normal[AdjacencyMatrix[graph]][[ordering, ordering]];
   Flatten@MapIndexed[Drop[#1, Max[external, #2-1]]&, matrix]
];
AdjaceTopology // secure;

DefineAllowedTopologies[] :=
Module[{allTopologies, single, combined},
   allTopologies = AllTopologies@$externalFieldNumbers;
   If[Head@allTopologies =!= List, Return[]];
   combined = Select[allTopologies, FreeQ[#, _Integer]&];
   single = Complement[allTopologies, combined];
   If[single =!= {}, DefineSingle/@ single];
   If[combined =!= {}, DefineCombined/@ combined];
];

DefineSingle[name_Symbol -> adjacencyVector:{__Integer}] :=
With[{function = name, vector = adjacencyVector},
   function[t:_?IsTopology] := AdjaceTopology@t === vector;
   function // secure;
];

DefineCombined[name_Symbol -> singleTopologies:{__Symbol}] :=
With[{function = name, list = singleTopologies},
   function[t:_?IsTopology] := Or@@ Through@list@t;
   function // secure;
];

GetExcludeTopologies[] :=
Once@Module[{all, name, default, loaded},
   default = {
      "Irreducible" -> (FreeQ[#, FeynArts`Internal]&),
      "Triangles"   -> (FreeQ[FeynArts`ToTree@#, FeynArts`Centre@Except@3]&)
   };
   loaded = If[MatchQ[#, {__}], #, {}]&[topologies@$loopNumber];
   loaded = Rule[SymbolName@First@#, With[{fun = Last@#}, fun@#&]] &/@ Utils`UnzipRules@loaded;

   all = Join[default, loaded];
   FeynArts`$ExcludeTopologies[name] = Function[Or@@Through[($expressionsToDerive/.all)@#]];
   name
];
GetExcludeTopologies // secure;

End[];
EndPackage[];
