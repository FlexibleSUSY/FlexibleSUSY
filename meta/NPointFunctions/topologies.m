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

define[getAdjacencyMatrix, {topology:`type`topology} :>
   Module[{
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
   ]
];

define[`topologyQ`pinguinT, {topology:`type`topology} :>
   Or[
      `topologyQ`trianglepinguinT@topology,
      `topologyQ`self1pinguinT@topology,
      `topologyQ`self3pinguinT@topology
   ]
];

define[`topologyQ`trianglepinguinT, {topology:`type`topology} :>
   getAdjacencyMatrix@topology === {
      {0,0,0,0,1,0,0,0},
      {0,0,0,0,0,1,0,0},
      {0,0,0,0,0,0,1,0},
      {0,0,0,0,0,1,0,0},
      {1,0,0,0,0,0,1,1},
      {0,1,0,1,0,0,0,1},
      {0,0,1,0,1,0,0,1},
      {0,0,0,0,1,1,1,0}
   }
];

define[`topologyQ`self1pinguinT, {topology:`type`topology} :>
   getAdjacencyMatrix@topology === {
      {0,0,0,0,1,0,0,0},
      {0,0,0,0,0,1,0,0},
      {0,0,0,0,0,0,1,0},
      {0,0,0,0,0,1,0,0},
      {1,0,0,0,0,0,0,2},
      {0,1,0,1,0,0,1,0},
      {0,0,1,0,0,1,0,1},
      {0,0,0,0,2,0,1,0}
   }
];

define[`topologyQ`self3pinguinT, {topology:`type`topology} :>
   getAdjacencyMatrix@topology === {
      {0,0,0,0,1,0,0,0},
      {0,0,0,0,0,1,0,0},
      {0,0,0,0,0,0,1,0},
      {0,0,0,0,0,1,0,0},
      {1,0,0,0,0,1,0,1},
      {0,1,0,1,1,0,0,0},
      {0,0,1,0,0,0,0,2},
      {0,0,0,0,1,0,2,0}
   }
];

define[`topologyQ`box, {topology:`type`topology} :>
   Or[
      `topologyQ`boxS@topology,
      `topologyQ`boxT@topology,
      `topologyQ`boxU@topology
   ]
];

define[`topologyQ`boxS, {topology:`type`topology} :>
   getAdjacencyMatrix@topology === {
      {0,0,0,0,1,0,0,0},
      {0,0,0,0,0,1,0,0},
      {0,0,0,0,0,0,1,0},
      {0,0,0,0,0,0,0,1},
      {1,0,0,0,0,1,1,0},
      {0,1,0,0,1,0,0,1},
      {0,0,1,0,1,0,0,1},
      {0,0,0,1,0,1,1,0}
   }
];

define[`topologyQ`boxT, {topology:`type`topology} :>
   getAdjacencyMatrix@topology === {
      {0,0,0,0,1,0,0,0},
      {0,0,0,0,0,1,0,0},
      {0,0,0,0,0,0,1,0},
      {0,0,0,0,0,0,0,1},
      {1,0,0,0,0,1,0,1},
      {0,1,0,0,1,0,1,0},
      {0,0,1,0,0,1,0,1},
      {0,0,0,1,1,0,1,0}
   }
];

define[`topologyQ`boxU, {topology:`type`topology} :>
   getAdjacencyMatrix@topology === {
      {0,0,0,0,1,0,0,0},
      {0,0,0,0,0,1,0,0},
      {0,0,0,0,0,0,1,0},
      {0,0,0,0,0,0,0,1},
      {1,0,0,0,0,0,1,1},
      {0,1,0,0,0,0,1,1},
      {0,0,1,0,1,1,0,0},
      {0,0,0,1,1,1,0,0}
   }
];

Module[{
      topologyReplacements = {
         Irreducible -> (FreeQ[#,FeynArts`Internal]&), (*@todo something weird with this definition*)
         Triangles -> (FreeQ[FeynArts`ToTree@#,FeynArts`Centre@Except@3]&)
      }
   },

getExcludeTopologies::usage =
"@brief Registers and returns a function, whose outcome - True or everything
        else - determines whether the topology is kept or discarded (see
        FeynArts manual).
@param keepProcesses Name(s) of processes to hold.
@returns _Symbol Generated name of topologies to hold.
@todo Add a catcher for the case 2 or more.";
define[getExcludeTopologies, {keepProcesses:{__Symbol}} :>
   Module[{
         excludeTopologyName,
         rules = Join[topologyReplacements, `settings`topologyReplacements]
      },
      FeynArts`$ExcludeTopologies[excludeTopologyName] = Switch[Length@keepProcesses,
         1,  keepProcesses[[1]] /. rules,
         _, (Or @@ Through[(keepProcesses /. rules)@#])&
      ];
      excludeTopologyName]
];

];

End[];
EndPackage[];
