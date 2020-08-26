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

BeginPackage["NPointFunctions`"];
Begin["NPointFunctions`internal`"];

getAdjacencyMatrix[topology:`type`topology] :=
Module[
   {
      propagatorPattern,needNewNumbers,adjacencies,adjacencyMatrix
   },

   propagatorPattern[i_,j_,f___] := _[_][_[_][i],_[_][j],f];

   needNewNumbers = And[
      Max@@(topology/.propagatorPattern[i_,j_,___]:>Sequence[i,j])>100,
      MatchQ[List@@topology,{propagatorPattern[_,_]..}]
   ];

   adjacencies = Tally[
      (List@@#)/.propagatorPattern[i_,j_,___]:>{{i,j},{j,i}}
   ] &@ If[needNewNumbers,FeynArts`TopologySort@#,#] &@ topology;

   adjacencyMatrix = Normal@SparseArray@Flatten[
      {#[[1,1]]->#[[2]],#[[1,2]]->#[[2]]} &/@ adjacencies
   ]
];

getAdjacencyMatrix // Utils`MakeUnknownInputDefinition;
getAdjacencyMatrix ~ SetAttributes ~ {Protected,Locked};

`topologyQ`pinguinT[topology:`type`topology] :=
Or[
   `topologyQ`trianglepinguinT@topology,
   `topologyQ`self1pinguinT@topology,
   `topologyQ`self3pinguinT@topology
];

`topologyQ`pinguinT // Utils`MakeUnknownInputDefinition;
`topologyQ`pinguinT ~ SetAttributes ~ {Protected,Locked};

`topologyQ`trianglepinguinT[topology:`type`topology] :=
getAdjacencyMatrix@topology === {
   {0,0,0,0,1,0,0,0},
   {0,0,0,0,0,1,0,0},
   {0,0,0,0,0,0,1,0},
   {0,0,0,0,0,1,0,0},
   {1,0,0,0,0,0,1,1},
   {0,1,0,1,0,0,0,1},
   {0,0,1,0,1,0,0,1},
   {0,0,0,0,1,1,1,0}
};

`topologyQ`trianglepinguinT // Utils`MakeUnknownInputDefinition;
`topologyQ`trianglepinguinT ~ SetAttributes ~ {Protected,Locked};

`topologyQ`self1pinguinT[topology:`type`topology] :=
getAdjacencyMatrix@topology === {
   {0,0,0,0,1,0,0,0},
   {0,0,0,0,0,1,0,0},
   {0,0,0,0,0,0,1,0},
   {0,0,0,0,0,1,0,0},
   {1,0,0,0,0,0,0,2},
   {0,1,0,1,0,0,1,0},
   {0,0,1,0,0,1,0,1},
   {0,0,0,0,2,0,1,0}
};

`topologyQ`self1pinguinT // Utils`MakeUnknownInputDefinition;
`topologyQ`self1pinguinT ~ SetAttributes ~ {Protected,Locked};

`topologyQ`self3pinguinT[topology:`type`topology] :=
getAdjacencyMatrix@topology === {
   {0,0,0,0,1,0,0,0},
   {0,0,0,0,0,1,0,0},
   {0,0,0,0,0,0,1,0},
   {0,0,0,0,0,1,0,0},
   {1,0,0,0,0,1,0,1},
   {0,1,0,1,1,0,0,0},
   {0,0,1,0,0,0,0,2},
   {0,0,0,0,1,0,2,0}};

`topologyQ`self3pinguinT // Utils`MakeUnknownInputDefinition;
`topologyQ`self3pinguinT ~ SetAttributes ~ {Protected,Locked};

`topologyQ`box[topology:`type`topology] :=
Or[
   `topologyQ`boxS@topology,
   `topologyQ`boxT@topology,
   `topologyQ`boxU@topology
];

`topologyQ`box // Utils`MakeUnknownInputDefinition;
`topologyQ`box ~ SetAttributes ~ {Protected,Locked};

`topologyQ`boxS[topology:`type`topology] :=
getAdjacencyMatrix@topology === {
   {0,0,0,0,1,0,0,0},
   {0,0,0,0,0,1,0,0},
   {0,0,0,0,0,0,1,0},
   {0,0,0,0,0,0,0,1},
   {1,0,0,0,0,1,1,0},
   {0,1,0,0,1,0,0,1},
   {0,0,1,0,1,0,0,1},
   {0,0,0,1,0,1,1,0}
};

`topologyQ`boxS // Utils`MakeUnknownInputDefinition;
`topologyQ`boxS ~ SetAttributes ~ {Protected,Locked};

`topologyQ`boxT[topology:`type`topology] :=
getAdjacencyMatrix@topology === {
   {0,0,0,0,1,0,0,0},
   {0,0,0,0,0,1,0,0},
   {0,0,0,0,0,0,1,0},
   {0,0,0,0,0,0,0,1},
   {1,0,0,0,0,1,0,1},
   {0,1,0,0,1,0,1,0},
   {0,0,1,0,0,1,0,1},
   {0,0,0,1,1,0,1,0}
};

`topologyQ`boxT // Utils`MakeUnknownInputDefinition;
`topologyQ`boxT ~ SetAttributes ~ {Protected,Locked};

`topologyQ`boxU[topology:`type`topology] :=
getAdjacencyMatrix@topology === {
   {0,0,0,0,1,0,0,0},
   {0,0,0,0,0,1,0,0},
   {0,0,0,0,0,0,1,0},
   {0,0,0,0,0,0,0,1},
   {1,0,0,0,0,0,1,1},
   {0,1,0,0,0,0,1,1},
   {0,0,1,0,1,1,0,0},
   {0,0,0,1,1,1,0,0}
};

`topologyQ`boxU // Utils`MakeUnknownInputDefinition;
`topologyQ`boxU ~ SetAttributes ~ {Protected,Locked};

topologyReplacements::usage = "
@brief List of topology replacement rules for a processes to keep.
@note R.h.s. should be pure functions of one argument.";
topologyReplacements =
{
   Irreducible -> (FreeQ[#,FeynArts`Internal]&), (*@todo something weird with this definition*)
   Triangles -> (FreeQ[FeynArts`ToTree@#,FeynArts`Centre@Except@3]&)
};
topologyReplacements ~ SetAttributes ~ {Protected, Locked};

getExcludedTopologies::usage =
"@brief Registers and returns a function, whose outcome - True or everything
        else - determines whether the topology is kept or discarded (see
        FeynArts manual).
@param keepProcesses Name(s) of processes to hold.
@returns _Symbol Generated name of topologies to hold.";
getExcludedTopologies[keepProcesses:{__Symbol}] :=
Module[
   {
      excludeTopologyName,
      rules = Join[topologyReplacements, `settings`topologyReplacements]
   },
   FeynArts`$ExcludeTopologies[excludeTopologyName] = Switch[Length@keepProcesses,
      1,  keepProcesses[[1]] /. rules,
      _, (Or @@ Through[(keepProcesses /. rules)@#])&
   ];
   excludeTopologyName];
getExcludedTopologies // Utils`MakeUnknownInputDefinition;
getExcludedTopologies ~ SetAttributes ~ {Protected, Locked};

End[];
EndPackage[];
