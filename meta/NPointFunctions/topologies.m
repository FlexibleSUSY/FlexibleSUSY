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

(*          v-- Number of incoming particles.                                *)
(*             v-- Number of outgoing particles.                             *)
topologies[{2, 2}] = {
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
(* ^----^ Use these names in the 'settings.m' files for observables.         *)
(*          ^---------------------^ Use 'adjace' for a single topology.      *)
(*            ^--^ Use symbols from lhs. for a combination of topologies.    *)

adjace::usage = "
@brief Archivates the adjacency matrix of a given topology. Example:
       In[1]:= << FeynArts`;
               topologies = CreateTopologies[0, 1 -> 2];
               adjace[topologies[[1]]]
       Out[1]= {1, 1, 1, 0}";
(*              v-----------v Replace by _ in the notebook.                  *)
adjace[topology:type`topology] :=
Module[{external, simple, graph, ordering, matrix},
(* v------v Number of external particles.                                    *)
(*            v--------------------------v Replace by number in the notebook.*)
   external = Tr@`options`observable@Outer;
   simple = (#/. h_[i_, j_, _] :> h[i, j])&/@ topology;
   graph = List@@ (UndirectedEdge@@@ FeynArts`TopologySort@simple)/.
      FeynArts`Vertex[_][i_] :> i;
   ordering = Ordering@VertexList@graph;
   matrix = Normal[AdjacencyMatrix[graph]][[ordering, ordering]];
   Flatten@MapIndexed[Drop[#1, Max[external, #2-1]]&, matrix]
];
adjace // tools`secure;

define[topologies] :=
Module[{all, single, combined},
   all = topologies@`options`observable@Outer;
   If[Head@all =!= List, Return[]];
   combined = Select[all, FreeQ[#, _Integer]&];
   single = Complement[all, combined];
   If[single =!= {}, defineSingle/@ single];
   If[combined =!= {}, defineCombined/@ combined];
];

defineSingle[name_Symbol -> adjacencyVector:{__Integer}] :=
With[{function = name, vector = adjacencyVector},
   function[t:type`topology] := adjace@t === vector;
   function // tools`secure;
];
defineSingle // tools`secure;

defineCombined[name_Symbol -> singleTopologies:{__Symbol}] :=
With[{function = name, list = singleTopologies},
   function[t:type`topology] := Or@@ Through@list@t;
   function // tools`secure;
];
defineCombined // tools`secure;

getTopology[d:type`diagram] := First@d;
getTopology // tools`secure;

getExcludeTopologies::usage = "
@brief Registers a function, whose outcome (``True`` or everything else)
       determines whether the topology is kept or not.
@returns A name of generated function.";
getExcludeTopologies[] :=
Once@Module[{all, name, set, default},
   default = {
      "Irreducible" -> (FreeQ[#, FeynArts`Internal]&),
      "Triangles"   -> (FreeQ[FeynArts`ToTree@#, FeynArts`Centre@Except@3]&)
   };
   set = If[MatchQ[#, {__}], #, {}]&[topologies@`options`loops[]];
   set = Rule[SymbolName@First@#, With[{fun = Last@#}, fun@#&]]&/@ set;
   all = Join[default, set];
   FeynArts`$ExcludeTopologies[name] = Function[Or@@Through[
      (`options`processes[]/.all)@#]];
   name
];
getExcludeTopologies // tools`secure;

End[];
EndPackage[];
