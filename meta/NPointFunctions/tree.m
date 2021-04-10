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

removeColors[set:`type`amplitudeSet] :=
   Delete[set, Position[set, FeynArts`Index[Global`Colour, _Integer]]];
removeColors // secure;

plant::usage = "
@brief Creates analytical form for amplitudes and forms a ``tree`` object
       from diagram-amplitude data.
@param d A set of diagrams.
@returns A ``tree`` object.";
plant[d:`type`diagramSet] :=
   Module[{diag, ampl, topologies, generic, gN, classes, cN, repl, tree, gray},
      gray = removeColors@FeynArts`CreateFeynAmp@d;
      {diag, ampl} = {List@@d, List@@gray} /. _FeynArts`Insertions -> List;
      topologies = First/@ diag;
      generic = First/@ #&/@ Last/@ diag;
      gN = Length/@ generic;
      classes = Last/@ #&/@ Last/@ diag;
      cN = Length/@ #&/@ classes;

      generic = MapThread[FeynArts`Insertions[Generic][#1, #2]&,
         {Flatten[generic, 1], Most/@List@@gray}];
      repl = Thread/@ Flatten[thread/@ Last/@ ampl, 1];
      classes = MapThread[FeynArts`Insertions[FeynArts`Classes][#1, #2]&,
         {Flatten[classes, 2], repl}];

      tree = node@{Head@d, Head@gray};
      Do[tree = grow[tree, topologies[[i]]];
         Do[tree = grow[tree, generic[[j]], {i}];
            Do[tree = grow[tree, classes[[k]], {i, j}];,
               {k, cN[[i, j]]}];
            classes = Drop[classes, cN[[i, j]]];,
            {j, gN[[i]]}];
         generic = Drop[generic, gN[[i]]];,
         {i, Length@diag}];
      tree];
plant // secure;

grow[expr_node, leaf_] :=
   Append[expr, node@leaf];
grow[expr_node, leaf_, pos:{__Integer}] :=
   Module[{place, elem},
      place = pos /. i_Integer :> i+1;
      elem = Extract[expr, place];
      ReplacePart[expr, place -> grow[elem, leaf]]];
grow::partw = "Invalid position `1`";
grow // secure;

diagrams[tree:`type`tree] :=
   tree /.
      node[e:`type`head, rest__] :> First[e]@rest /.
      node[e:`type`topology, rest__] :>
         Rule[e, FeynArts`Insertions[Generic][rest]]  /.
      node[e:`type`generic, rest__] :>
         First[e] -> FeynArts`Insertions[FeynArts`Classes]@rest /.
      node[e:`type`classes] :> First@e;
diagrams // secure;

amplitudes[tree:`type`tree] :=
   tree /.
      node[e:`type`head, rest__] :> Part[e, 2]@rest /.
      node[e:`type`topology, rest__] :> rest /.
      node[e:`type`classes] :> Last@e /.
      node[e:`type`generic, rest__] :> Append[Part[e, 2], wrap@rest];
amplitudes // secure;

fields[tree:`type`tree, Flatten] :=
   Flatten[fields@tree, 1];
fields[tree:`type`tree] :=
   tree /. node[e:`type`head, __] :> List@@(FeynArts`Process /. List@@First@e);
fields // secure;

picture[tree:`type`tree] :=
   Module[{out = {}, directory, name},
      name = StringJoin[ToString /@ (
         `rules`fields@Join[fields[tree, Flatten],
            `options`processes[]] /. e_@{_} :> e)];
      directory = DirectoryName[FeynArts`$Model<>".mod"];
      FeynArts`Paint[diagrams@tree,
         FeynArts`PaintLevel -> {Generic},
         FeynArts`ColumnsXRows -> 1,
         FeynArts`FieldNumbers -> True,
         FeynArts`SheetHeader -> None,
         FeynArts`Numbering -> FeynArts`Simple,
         DisplayFunction :> (AppendTo[out, #] &/@ Render[##, "JPG"] &)];
      Put[out, FileNameJoin@{directory, name<>".m"}]];
picture // secure;

thread[func_[lhs:{__}, rhs:{{__}..}]] :=
   Thread[func[Array[lhs&, Length@rhs], rhs]];
thread // secure;

wrap[data:{Rule[_, _]..}..] :=
   Module[{lhs, rhs},
      lhs = First/@First@{data};
      rhs = FeynArts`Insertions[FeynArts`Classes]@@(Last/@#&/@{data});
      lhs -> rhs];
wrap // secure;

erase::usage = "
@brief Erases nodes (or diagrams).
       Can be done *via*:

       1. Settings. Can be used to apply the following settings:

          ```settings`diagrams``
             It is used to remove diagrams *via* actions.
          ```settings`amplitudes``
             It is used to remove diagrams and amplitudes *via* actions.

       2. Topology
       3. Generic or classes level data

       Usually nodes are removed from the deepest to the highest level.

@param input A ``tree`` object or a set of diagrams (for ``settings``).
@param settings Data, which specifies replacements for each topology.
@param tQ A function to select a *topology* ``node[id, ___]`` if::

           In[1]:= tQ[id]
          Out[1]:= True

@param fun A function to erase *class* or *generic* ``node[id, ___]`` if::

           In[1]:= fun[id]
          Out[1]:= True

@param n A ``node``, which is checked by ``fun``.
@returns Nodes, cleaned by ``fun``.";
erase[input:`type`tree|`type`diagramSet, settings:{__Rule}|Default] :=
   Module[{res = input},
      Set[res, applyAction[res, #]]&/@ getActions@settings;
      res];
erase[input:`type`tree, tQ_, fun_] :=
   Module[{res = input},
      res = res /. e:node[t:`type`topology /; tQ@t, __] :> erase[e, fun];
      If[MatchQ[res, node@`type`head],
         Utils`AssertOrQuit[_, erase::empty, tQ, fun]];
      res];
erase[n:node[`type`topology, __], fun_] :=
   Module[{res = n},
      res = res /. e:node[`type`generic, __] :> erase[e, fun];
      res = res /. e:node@`type`topology :> Sequence[];
      If[fun@res, res, ##&[]]];
erase[n:node[`type`generic, __], fun_] :=
   Module[{res = n},
      res = res /. e:node@`type`classes :> erase[e, fun];
      res = res /. e:node@`type`generic :> Sequence[];
      If[fun@res, res, ##&[]]];
erase[n:node@`type`classes, fun_] :=
   If[fun@n, n, ##&[]];
erase::empty = "
After the usage of
tQ: \"`1`\"
fun: \"`2`\"
all topologies were erased.";
erase // secure

combinatoricalFactors::usage = "
@brief Finds combinatirical factors for a given input.
@param tree A ``tree`` object.
@returns ``List`` of combinatorical factors for a given ``tree``.";
combinatoricalFactors[tree:`type`tree] :=
   combinatoricalFactors /@ List@@amplitudes@tree;
combinatoricalFactors[_[_,_,_, generic_ -> _[_][classes__]]] :=
   {classes}[[All, #[[1, 1]]]] /.
      {FeynArts`IndexDelta[___] -> 1, FeynArts`SumOver[__] -> 1} &@
         Position[generic, FeynArts`RelativeCF];
combinatoricalFactors // secure;

colorFactors::usage = "
@brief Finds color factors for a given input.
@param tree A ``tree`` object.
@param diagram A diagram to work with.
@param props A ``Sequence`` of propagators. First numbers of vertices inside
       propagators are sorted by ``FeynArts``.
@param rules A ``Sequence`` of field replacement rules on generic level.
@returns ``List`` (for a given topology) of several ``List``
         (for generic fields) of colour factors: ``{{__}..}``.
@note External fields always come at first places in adjacency matrix.";
colorFactors[tree:`type`tree] :=
   `rules`fields@Flatten[colorFactors /@ List@@diagrams@tree, 1];
colorFactors[diagram:Rule[_[_][props__], _[_][_[__][rules__]->_,___]]] :=
Module[{propPatt, adjacencyMatrix, externalRules, genericDiagram},
   propPatt[i_, j_, f_] := _[_][_[_][i], _[_][j], f];
   adjacencyMatrix =
      Module[{adjs},
         adjs = Tally[{props} /. propPatt[i_, j_, _] :> {{i, j}, {j, i}}];
         Normal@SparseArray@Flatten[
            {#[[1,1]] -> #[[2]], #[[1,2]] -> #[[2]]} &/@ adjs]];
   externalRules = Cases[{rules}, HoldPattern[_ -> _@__]];
   genericDiagram =
      Module[{fld},
         fld = Flatten[{props} /. propPatt[i_, j_, f_] :>
            {{j, i, -f},{i, j, f}}, 1];
         GatherBy[SortBy[fld, First], First] /. {_Integer, _Integer, f_} :>
            f] /. Join[ {#} -> #&/@ First/@externalRules];
   Map[
      CXXDiagrams`ColourFactorForIndexedDiagramFromGraph[
         CXXDiagrams`IndexDiagramFromGraph[
            genericDiagram /. externalRules /. #,
            adjacencyMatrix],
         adjacencyMatrix]&,
      fieldInsertions[diagram, True],
      {2}]];
colorFactors // secure;

End[];
EndPackage[];
