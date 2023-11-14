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

RemoveColors[expr_] := Delete[expr, Position[expr, _?IsColorIndex]];
RemoveColors // secure;

GenerateDiagrams[in_, out_] :=
Module[{topologies, diagrams},
   topologies = FeynArts`CreateTopologies[
      $loopNumber,
      Length@in -> Length@out,
      FeynArts`ExcludeTopologies -> GetExcludeTopologies[]
   ];
   Utils`AssertOrQuit[Length@topologies > 0, GenerateDiagrams::topologies];

   diagrams = FeynArts`InsertFields[
      topologies,
      ToExpression@in -> ToExpression@out
   ];
   Utils`AssertOrQuit[Length@diagrams > 0, GenerateDiagrams::diagrams];

   node[{Head@#}, Sequence@@#]&[diagrams] /.
      Rule[t_?IsTopology, rest_] :> node[t, rest] /.
      (h:_@Generic)[a__] :> Sequence@@(node@*h@@#&/@{a}) /.
      (h:_@Generic)[a_, rest__] :> Sequence[h@a, rest] /.
      (h:_@FeynArts`Classes)[a__] :> Sequence@@(node@*h/@{a})
];
GenerateDiagrams // secure;
GenerateDiagrams::topologies = "CreateTopologies generated zero topologies.";
GenerateDiagrams::diagrams = "InsertFields generated zero field insertions.";

GenerateColorlessAmplitudes[tree_?IsTree] :=
Module[{amps, generic, classes, i = 1, j = 1},
   amps = RemoveColors@FeynArts`CreateFeynAmp@ExtractDiagrams@tree;
   generic = Most/@List@@amps;
   classes = (Last/@List@@amps) /. (lhs_ -> _@rhs__) :> Sequence@@(Thread[lhs -> #]&/@{rhs});
   RemoveColors@tree /.
      node[e_?IsTreeHead, r__] :> node[Append[e, Head@amps], r] /.
      node[e_?IsGeneric, r__] :> node[Append[e, generic[[i++]]], r] /.
      node[e_?IsClasses] :> node@Append[e, classes[[j++]]]
];
GenerateColorlessAmplitudes // secure;

ExtractDiagrams[tree_?IsTree] :=
tree /.
   node[e_?IsTreeHead, rest__] :> First[e][rest] /.
   node[e_?IsTopology, rest__] :> Rule[e, FeynArts`Insertions[Generic][rest]]  /.
   node[e_?IsGeneric, rest__]  :> First[e] -> FeynArts`Insertions[FeynArts`Classes][rest] /.
   node[e_?IsClasses]          :> First[e];

ExtractAmplitudes[tree_?IsTree] :=
tree /.
   node[e_?IsTreeHead, rest__] :> Part[e, 2]@rest /.
   node[e_?IsTopology, rest__] :> rest /.
   node[e_?IsClasses] :> Last@e /.
   node[e_?IsGeneric, rest__] :> Append[Part[e, 2], wrap@rest];

wrap[data:{Rule[_, _]..}..] :=
   Module[{lhs, rhs},
      lhs = First/@First@{data};
      rhs = FeynArts`Insertions[FeynArts`Classes]@@(Last/@#&/@{data});
      lhs -> rhs];
wrap // secure;

GetFields[tree_?IsTree, Flatten] :=
   Flatten[GetFields@tree, 1];
GetFields[tree_?IsTree] :=
   tree /. node[e_?IsTreeHead, __] :> List@@(FeynArts`Process /. List@@First@e);
GetFields // secure;

ExportFeynArtsPaint[tree_?IsTree] :=
Module[{out = {}, directory, name},
   name = StringJoin[ToString /@ (
      FieldRules@Join[GetFields[tree, Flatten],
         $expressionsToDerive] /. e_@{_} :> e)
   ];
   directory = DirectoryName[FeynArts`$Model<>".mod"];
   FeynArts`Paint[ExtractDiagrams@tree,
      FeynArts`PaintLevel -> {Generic},
      FeynArts`ColumnsXRows -> 1,
      FeynArts`FieldNumbers -> True,
      FeynArts`SheetHeader -> None,
      FeynArts`Numbering -> FeynArts`Simple,
      DisplayFunction :> (AppendTo[out, #] &/@ Render[##, "JPG"] &)
   ];
   Put[out, FileNameJoin@{directory, name<>".m"}]
];
ExportFeynArtsPaint // secure;

RemoveNode::usage = "Removes node if both tQ[id] is True and fun[node, info] is True.";
RemoveNode[n:node[_?IsTopology, __], fun_, info__] :=
   n /. e:node[_?IsGeneric, __] :> RemoveNode[e, fun, info] /.
      node@_?IsTopology :> Sequence[];
RemoveNode[n:node[_?IsGeneric, __], fun_, info__] :=
   If[fun[#, info], # /. node[_?IsGeneric] :> Sequence[], ##&[]]&[n /. e:node[_?IsClasses] :> RemoveNode[e, fun, info]];
RemoveNode[n:node[_?IsClasses], fun_, info__] :=
   If[fun[n, info], n, ##&[]];
RemoveNode // secure;

TreeHead[tree_?IsTree] := tree[[1, 1]];
TreeHead // secure;

CombinatoricalFactors[tree_?IsTree] := CombinatoricalFactors /@ List@@ExtractAmplitudes@tree;
CombinatoricalFactors[_[_,_,_, generic_ -> _[_][classes__]]] :=
{classes}[[All, #[[1, 1]]]] /.
   {FeynArts`IndexDelta[___] -> 1, FeynArts`SumOver[__] -> 1} &@
      Position[generic, FeynArts`RelativeCF];
CombinatoricalFactors // secure;

ColorFactors[tree:_?IsTree] := FieldRules@Flatten[ColorFactors /@ List@@ExtractDiagrams@tree, 1];
ColorFactors[diagram:Rule[_[_][props__], _[_][_[__][rules__]->_,___]]] :=
Module[{propPatt, adjacencyMatrix, externalRules, genericDiagram},
   propPatt[i_, j_, f_] := _[_][_[_][i], _[_][j], f];
   adjacencyMatrix = Module[{adjs},
      adjs = Tally[{props} /. propPatt[i_, j_, _] :> {{i, j}, {j, i}}];
      Normal@SparseArray@Flatten[
         {#[[1,1]] -> #[[2]], #[[1,2]] -> #[[2]]} &/@ adjs]
   ];
   externalRules = Cases[{rules}, HoldPattern[_ -> _@__]];
   genericDiagram = Module[{fld},
      fld = Flatten[{props} /. propPatt[i_, j_, f_] :>
         {{j, i, -f},{i, j, f}}, 1];
      GatherBy[SortBy[fld, First], First] /. {_Integer, _Integer, f_} :>
         f
   ] /. Join[ {#} -> #&/@ First/@externalRules];
   Map[
      CXXDiagrams`ColourFactorForIndexedDiagramFromGraph[
         CXXDiagrams`IndexDiagramFromGraph[
            genericDiagram /. externalRules /. #,
            adjacencyMatrix
         ],
         adjacencyMatrix
      ]&,
      GetFieldInsertions[diagram, True],
      {2}
   ]
];
ColorFactors // secure;

End[];
EndPackage[];
