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

lengthyQ::usage = "
@brief Checks, whether the ``Length`` of expression is zero or not.
       In the latter case returns input, otherwise writes an error message
       and stops the evaluation.
@param input Any expression.
@returns An input in case of non-zero ``Length`` of the input.";
lengthyQ[input:_] := With[{sym = Head@Unevaluated@input},
   If[Length@input =!= 0,
      input,
      sym::length = "The input has a zero lenght.";
      Utils`AssertOrQuit[_, sym::length];]];
lengthyQ // Utils`MakeUnknownInputDefinition;
lengthyQ ~ SetAttributes ~ {Protected, HoldFirst};

plant::usage = "
@brief Creates diagrams and amplitudes, forms a ``tree`` object.
@param in An expression, containing ``FeynArts`` fields (each as ``String``),
       required to be *incoming* fields for the amplitude.
@param out An expression, containing ``FeynArts`` fields (each as ``String``),
       required to be *outgoing* fields for the amplitude.
@param tree A ``tree`` object.
@returns A ``tree`` object.";
plant[in_, out_] :=
   Module[{topologies, diagrams},
      topologies = lengthyQ@FeynArts`CreateTopologies[
         `options`loops[],
         Length@in -> Length@out,
         FeynArts`ExcludeTopologies -> getExcludeTopologies[]];
      diagrams = lengthyQ@FeynArts`InsertFields[topologies,
         ToExpression@in -> ToExpression@out];

      node[{Head@#}, Sequence@@#]&[diagrams] /.
         Rule[t:`type`topology, rest_] :> node[t, rest] /.
         (h:_@Generic)@a__ :> Sequence@@(node@*h@@#&/@{a}) /.
         (h:_@Generic)[a_, rest__] :> Sequence[h@a, rest] /.
         (h:_@FeynArts`Classes)@a__ :> Sequence@@(node@*h/@{a})];
plant[tree:`type`tree] :=
   Module[{amps, generic, classes, i = 1, j = 1},
      amps = removeColors@FeynArts`CreateFeynAmp@diagrams@tree;
      generic = Most/@List@@amps;
      classes = (Last/@List@@amps) /. (lhs_ -> _@rhs__) :>
         Sequence@@(Thread[lhs -> #]&/@{rhs});
      tree /.
         node[e:`type`head, r__] :> node[Append[e, Head@amps], r] /.
         node[e:`type`generic, r__] :> node[Append[e, generic[[i++]]], r] /.
         node[e:`type`classes] :> node@Append[e, classes[[j++]]]];
plant // secure;

info[tree:`type`tree, str_String] :=
   (  Print@str;
      Print[" in total: ",
         Length@Cases[tree, `type`generic, Infinity], " Generic, ",
         Length@Cases[tree, `type`classes, Infinity], " Classes insertions"];
      tree);
info // secure;

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

wrap[data:{Rule[_, _]..}..] :=
   Module[{lhs, rhs},
      lhs = First/@First@{data};
      rhs = FeynArts`Insertions[FeynArts`Classes]@@(Last/@#&/@{data});
      lhs -> rhs];
wrap // secure;

cut::usage = "
@brief Removes nodes. Can be done *via*:

       1. Settings. Use to apply the following settings:

          ```settings`diagrams``
             Contains settings, specific for diagrams.
          ```settings`amplitudes``
             Contains settings, specific for diagrams and amplitudes.
       2. Topology
       3. Generic or classes level data
       Nodes are removed from the deepest to the highest level.
@param tree A ``tree`` object.
@param settings Data, which specifies replacements for topologies.
@param tQ A function to select a *topology* ``node[id, __]`` if::

           In[1]:= tQ[id]
          Out[1]:= True
@param fun A function to remove *class* or *generic* node if::

           In[1]:= fun[node, info]
          Out[1]:= True
@param info A ``Sequence`` of topology and topology list.
@param n A node to check.
@returns Nodes, cleaned by ``fun``.";
cut[tree:`type`tree, settings:{__Rule}|Default] :=
   Module[{res = tree},
      Set[res, applyAction[res, #]]&/@ getActions@settings;
      res];
cut[tree:`type`tree, tQ_, fun_] :=
   tree /.
      e:node[t:`type`topology /; tQ@t, __] :> cut[e, fun, t, head@tree];
cut[n:node[`type`topology, __], fun_, info__] :=
   n /. e:node[`type`generic, __] :> cut[e, fun, info] /.
      node@`type`topology :> Sequence[];
cut[n:node[`type`generic, __], fun_, info__] :=
   If[fun[#, info], # /. node@`type`generic :> Sequence[], ##&[]]&[
      n /. e:node@`type`classes :> cut[e, fun, info]];
cut[n:node@`type`classes, fun_, info__] := If[fun[n, info], n, ##&[]];
cut // secure;

LoopFields[node[id_, ___], info__] :=
   FeynArts`LoopFields[First@id, info];

TreeFields[node[id_, ___], info__] :=
   FeynArts`TreeFields[First@id, info];

head[tree:`type`tree] := tree[[1, 1]];
head // secure;

combinatoricalFactors::usage = "
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
