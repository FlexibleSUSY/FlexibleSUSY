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

adjace[topology:`type`topology] :=
Module[{propagatorPattern, needNewNumbers, adjacencies, matrix},

   propagatorPattern[i_,j_,f___] := _[_][_[_][i],_[_][j],f];

   needNewNumbers = And[
      Max@@(topology/.propagatorPattern[i_,j_,___]:>Sequence[i,j])>100,
      MatchQ[List@@topology,{propagatorPattern[_,_]..}]
   ];

   adjacencies = Tally[
      (List@@#)/.propagatorPattern[i_,j_,___]:>{{i,j},{j,i}}
   ] &@ If[needNewNumbers,FeynArts`TopologySort@#,#] &@ topology;

   matrix = Normal@SparseArray@Flatten[{#[[1,1]]->#[[2]],#[[1,2]]->#[[2]]} &/@ adjacencies];
   Flatten@Table[Drop[#[[i]], i-1], {i, Length@#}]&@matrix
];
adjace // Utils`MakeUnknownInputDefinition;
adjace ~ SetAttributes ~ {Protected, Locked};

define[`topologyQ`pinguinT, {t:`type`topology} :>
   Or[
      `topologyQ`trianglepinguinT@t,
      `topologyQ`self1pinguinT@t,
      `topologyQ`self3pinguinT@t
   ]
];

define[`topologyQ`trianglepinguinT, {t:`type`topology} :>
   adjace@t === {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0}];
define[`topologyQ`self1pinguinT, {t:`type`topology} :>
   adjace@t === {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,2,0,1,0,0,1,0}];
define[`topologyQ`self3pinguinT, {t:`type`topology} :>
   adjace@t === {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,0,2,0}];

define[`topologyQ`box, {t:`type`topology} :>
   Or[`topologyQ`boxS@t, `topologyQ`boxT@t, `topologyQ`boxU@t]];

define[`topologyQ`boxS, {t:`type`topology} :>
   adjace@t === {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,0,1,0,1,0}];
define[`topologyQ`boxT, {t:`type`topology} :>
   adjace@t === {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0}];
define[`topologyQ`boxU, {t:`type`topology} :>
   adjace@t === {0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,0,0}];

Module[{topologyReplacements},

topologyReplacements = {
   (*@todo something weird with this definition*)
   Irreducible -> (FreeQ[#, FeynArts`Internal]&),
   Triangles   -> (FreeQ[FeynArts`ToTree@#, FeynArts`Centre@Except@3]&)};

getExcludeTopologies::usage = "
@brief Registers a function, whose outcome (True or everything else) determines
       whether the topology is kept or discarded.
@param keep A list of processes to keep.
@returns A name of generated function.";
getExcludeTopologies[keep:{__Symbol}] :=
Module[{all, name},
   all = Join[topologyReplacements, `settings`topologyReplacements];
   FeynArts`$ExcludeTopologies[name] = Function[Or@@Through[(keep/.all)@#]];
   name];
getExcludeTopologies // Utils`MakeUnknownInputDefinition;
getExcludeTopologies ~ SetAttributes ~ {Protected, Locked};

];

End[];
EndPackage[];
