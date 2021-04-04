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

getActions::usage = "
@brief Converts ```settings`diagrams`` or ```settings`amplitudes`` to actions.
       The structure of settings is the following::

          {  <symbol> -> {  {  <if-in-$Processes action>
                               ...},
                            {  <if-not-in-$Processes action>
                               ...}}
             ..}

       <symbol>
         Any symbol, which defines the <action> to be done with diagrams,
         whether <symbol> is present or not in ``$Processes``.
       <action>
         Has the syntax of ``action`` for ``diagrams`` or ``amplitudes``.
@param settings An entry to define actions.
@returns A set of settings.";
getActions[settings:Default] := {};
getActions[settings:{Rule[_, {{___}, {___}}]..}] :=
Module[{positiveRules, negativeRules, discardProcesses, clean, parsed},
   parsed = Rule[SymbolName@First@#, Last@#]&/@settings;
   positiveRules = parsed /. Rule[s:_, {p:_, _}] :> Rule[s, p];
   negativeRules = parsed /. Rule[s:_, {_, n:_}] :> Rule[s, n];
   discardProcesses = Complement[parsed[[All, 1]], $Processes];
   clean = RuleDelayed[#, Sequence[]] &/@ $Processes;
   DeleteDuplicates@Join[
      DeleteDuplicates@Flatten[$Processes /. positiveRules, 1],
      DeleteDuplicates@Flatten[discardProcesses /. negativeRules, 1]] /.
         clean];
getActions // secure;

deleteClasses[a:`type`amplitudeSet,
   topoAmpList:`type`pickTopoAmp, classesToSave:`type`saveAmpClass] :=
Module[{i, numbersOfAmplitudes, numAmp, res = a},
   numbersOfAmplitudes = Cases[topoAmpList, Rule[True, {e:__}] :> e];
   Do[numAmp = Part[numbersOfAmplitudes,i];
      res[[numAmp,4,2]] = res[[numAmp,4,2]][[numAmp/.classesToSave]];,
      {i,Length@numbersOfAmplitudes}];
   res];
deleteClasses[d:`type`diagramSet,
   topoAmpList:`type`pickTopoAmp, classesToSave:`type`saveAmpClass] :=
Module[{i, currentClasses, res = d},
   Do[If[topoAmpList[[i,1]],
         currentClasses = topoAmpList[[i,2]]/.classesToSave;
         currentClasses = Array[Rule[{#,2}, res[[i,2,#,2]][[
            Part[currentClasses,#]]]]&,Length@currentClasses];
         res[[i, 2]] = ReplacePart[res[[i, 2]], currentClasses];];,
      {i,Length@topoAmpList}];
   res];
deleteClasses // secure;

getTruePositions::usage = "
@brief Converts a list with a boolean variables to the list of positions for
       all true entries.
@param list A list of booleans.
@returns A list of integers (or an empty one)."
getTruePositions[list:{(True|False)...}] := Flatten@Position[list, True];
getTruePositions // secure;

`action`amplitudes = Sequence[
   {diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet},
   text_String[topologyQ_, {function:UnsameQ, name_, value_}]];

getClassVariables[amp:`type`amplitude] := amp[[4, 1]];
getClassVariables // secure;

getClassInsertions[amp:`type`amplitude] := amp[[4, 2]];
getClassInsertions // secure;

getClassRules[amp:`type`amplitude] :=
   Thread[getClassVariables@amp -> Transpose[List@@getClassInsertions@amp]];
getClassRules // secure;

applyAction@`action`amplitudes :=
Module[{ daPairs, amplitudeNumbers, saveClassRules, viPairs, insertions, res},
   daPairs = getAmplitudeNumbers[diagrams, topologyQ];
   If[!Or@@daPairs[[All, 1]], Return@{diagrams, amplitudes}];
   amplitudeNumbers = Cases[daPairs, Rule[True, {e:__}] :> e];
   saveClassRules = Table[
      viPairs = getClassRules@amplitudes[[i]];
      insertions = Cases[viPairs, (name -> {e:__}) :> e, Infinity ];
      i -> getTruePositions[function[#, value] &/@ insertions] /. {} -> All,
      {i, amplitudeNumbers}];
   res =
   {  deleteClasses[diagrams, daPairs, saveClassRules],
      deleteClasses[amplitudes, daPairs, saveClassRules]};
   Print@text;
   printDiagramsInfo@res[[1]];
   res];

`action`diagrams = "
@brief Represents an action, which is supposed to be done with diagrams.
@param diagrams A set of diagrams.
@param s A string with comment.
@param t A topology or a set of topologies, where action should be applied.
@param c A criterion, which is passed to a FeynArts`DiagramSelect.";
`action`diagrams = Sequence[diagrams:`type`diagramSet,
   s_String[t:_Symbol|_Function, c_Function]];
`action`diagramsCompact = Sequence[d:`type`diagramSet,
   s_String[t:_@__Symbol, c_Function]];

secureSelect[input_, h_, c_] :=
Module[{selected},
   selected = FeynArts`DiagramSelect[h@input, c];
   If[0 === Length@selected, Return@Null];
   input[[1]] -> selected[[1, 2]]];
secureSelect // secure;

applyAction@`action`diagrams :=
Module[{d = diagrams, h = Head@diagrams},
   d = If[t@#[[1]], secureSelect[#, h, c], #] &/@ d;
   d = d /. Null :> Sequence[];
   Print@s;
   printDiagramsInfo@d;
   d];
applyAction@`action`diagramsCompact :=
   applyAction[d, s[Through[(Or@@t)@#]&, c]];

`action`sum = Sequence[d:`type`diagramSet,
   s_String[t_Symbol, {n_Integer, f:_}]];
`action`hold = Sequence[d:`type`diagramSet,
   s_String[t_Symbol, {Hold, e_Integer}]];
`action`append = Sequence[d:`type`diagramSet,
   s_String[t_Symbol, {Append, (f:`type`field)[n_Integer] :> e_Integer}]];

Module[{template},
   template[text:_, topologyQ:_, realization:_] :=
   If[topologyQ@getTopology@#,
      Print@text;
      getTopology@# -> List@@Map[realization&, getInsertions@#],
      (##&)[]]&;

   Module[{restrict},
      restrict[field_, number_] := {number -> Or[field, -field]};
      applyAction@`action`sum :=
      List@@Map[template[s, t, restrict[f@d, n]], d];];

   Module[{delete},
      delete[e_] := With[{pos = {{2*e}, {2*e-1}}}, Delete[#, pos]&];
      applyAction@`action`hold :=
      Sequence@@Map[template[s, t, delete@e], d];];

   Module[{append},
      append[field_, number:_, e_] :=
      With[{rhs = (First /@ getZeroMassRules[])[[2*e-1]]},
         Append[#, genericMass[field, number] :> rhs]&];
      applyAction@`action`append :=
      Sequence@@Map[template[s, t, append[f, n, e]], d];];];
applyAction // secure;

End[];
EndPackage[];
