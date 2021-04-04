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

define::usage = "
@brief Defines a set of function with a given name in a safe way.
@param s A symbol, which represent a function name.
@param e A sequence of delayed rules. On lhs there is a list with pattern
       for a new function, on rhs there is function body.";
Module[{impl},
   impl[s:_Symbol, RuleDelayed[{p:___}, d:_]] := SetDelayed[s[p], d];
   impl ~ SetAttributes ~ {HoldAllComplete};
   define[s:_Symbol, e:RuleDelayed[{___},_]..] :=
   (  impl[s, ##] &@@@ Hold /@ {e};
      s // Utils`MakeUnknownInputDefinition;
      Protect@s;);];
define // Utils`MakeUnknownInputDefinition;
define ~ SetAttributes ~ {HoldAllComplete, Protected};

adjace[topology:`type`topology] :=
Module[{propagatorPattern, needNewNumbers, adjacencies, matrix, ext},
   ext = Count[topology, FeynArts`Incoming|FeynArts`Outgoing|FeynArts`External,
      Infinity, Heads -> True];
   propagatorPattern[i_,j_,f___] := _[_][_[_][i],_[_][j],f];
   needNewNumbers = And[
      Max@@(topology/.propagatorPattern[i_,j_,___]:>Sequence[i,j])>100,
      MatchQ[List@@topology,{propagatorPattern[_,_]..}]];
   adjacencies = Tally[
      (List@@#)/.propagatorPattern[i_,j_,___]:>{{i,j},{j,i}}] &@
         If[needNewNumbers,FeynArts`TopologySort@#,#] &@ topology;
   matrix = Normal@SparseArray@Flatten[{#[[1,1]]->#[[2]],#[[1,2]]->#[[2]]} &/@
      adjacencies];
   Flatten@Table[Drop[#[[i]], i-1+Max[ext-i+1, 0]], {i, Length@#}]&@matrix];
adjace // secure;

define[`topologyQ`penguinT, {t:`type`topology} :>
   Or[`topologyQ`trianglepenguinT@t,
      `topologyQ`self1penguinT@t,
      `topologyQ`self3penguinT@t]];

define[`topologyQ`trianglepenguinT, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,1,0}];
define[`topologyQ`self1penguinT, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,2,0,1,0,0,1,0}];
define[`topologyQ`self3penguinT, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,2,0}];

define[`topologyQ`penguinU, {t:`type`topology} :>
   Or[`topologyQ`trianglepenguinU@t,
      `topologyQ`self1penguinU@t,
      `topologyQ`self4penguinU@t]];

define[`topologyQ`trianglepenguinU, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1,0,1,0}];
define[`topologyQ`self1penguinU, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,2,0,1,0,0,1,0}];
define[`topologyQ`self4penguinU, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,2,0}];

define[`topologyQ`box, {t:`type`topology} :>
   Or[`topologyQ`boxS@t, `topologyQ`boxT@t, `topologyQ`boxU@t]];

define[`topologyQ`boxS, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,1,0,0,0,1,0,1,0}];
define[`topologyQ`boxT, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0}];
define[`topologyQ`boxU, {t:`type`topology} :>
   adjace@t === {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,0}];

getExcludeTopologies::usage = "
@brief Registers a function, whose outcome (``True`` or everything else)
       determines whether the topology is kept or not.
@param keep A list of processes to keep.
@returns A name of generated function.";
Module[{topologyReplacements},
   topologyReplacements = {
      "Irreducible" -> (FreeQ[#, FeynArts`Internal]&),
      "Triangles"   -> (FreeQ[FeynArts`ToTree@#, FeynArts`Centre@Except@3]&)};
   getExcludeTopologies[] :=
   Module[{all, name, set},
      set = If[#===Default, {}, #]&@`settings`topology;
      set = Rule[SymbolName@First@#, Last@#]&/@set;
      all = Join[topologyReplacements, set];
      FeynArts`$ExcludeTopologies[name] = Function[Or@@Through[
         ($Processes/.all)@#]];
      name];
   getExcludeTopologies // secure;];

End[];
EndPackage[];
