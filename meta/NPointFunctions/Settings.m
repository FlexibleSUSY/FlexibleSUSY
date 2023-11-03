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

Utils`DynamicInclude@"PatternChecks.m";
Utils`DynamicInclude@"ModifyMasses.m";

BeginPackage@"NPointFunctions`";
Begin@"`Private`";

With[{dir = DirectoryName@$InputFileName},
LoadAllSettings[] :=
(
   BeginPackage@"NPointFunctions`";
   Begin@"`Private`";
   If[FileExistsQ@#, Get@#]&@FileNameJoin@
      {ParentDirectory@dir, "Observables", $observableName, "NPointFunctions.m"};
   End[];
   EndPackage[];
);
];

ApplyObservableSetting[tree_?IsTree, settingName:diagrams|amplitudes] :=
Module[{doPresent, doAbsent, absent, todos, res = tree},
   {doPresent, doAbsent} = If[Head@# === List, #, {}]&@ settingName[$loopNumber, #] &/@ {Present, Absent};
   {doPresent, doAbsent} = Utils`UnzipRules/@ {doPresent, doAbsent};
   {doPresent, doAbsent} = Rule[SymbolName@First@#, Last@#] &/@ #&/@ {doPresent, doAbsent};
   absent = Complement[First/@doAbsent, $expressionsToDerive];
   todos = DeleteDuplicates@Flatten[
      Join[$expressionsToDerive /. doPresent, absent /. doAbsent]
   ];
   todos = Select[todos, Head@# === Rule&];
   Set[res, ApplySetting[res, settingName, #]] &/@ todos;
   res
];

ApplySetting[tree_?IsTree, settingName:diagrams|amplitudes, tQ_ -> {str_String, fun_}] :=
Module[{res, isReplaced = False},
   res = tree /. e:node[t:_?IsTopology /; tQ@t, __] :>
      (
         isReplaced = True;
         RemoveNode[e, fun, t, TreeHead@tree]
      );
   If[isReplaced,
      PrintInfo[tQ, settingName, str,
         ToString@Length@Cases[res, _?IsGeneric, Infinity] <> "g, " <>
         ToString@Length@Cases[res, _?IsClasses, Infinity] <> "c"
      ]
   ];
   res
];

GetObservableSetting[order] :=
If[Head@order[] === List,
   order[],
   Reverse@Range@Tr@$externalFieldNumbers
];

ApplyObservableSetting[tree:_?IsTree, settingName:regularization|momenta|sum|mass] :=
Module[{res = {tree}, default, head},
   If[Head@settingName@$loopNumber === List,
      AppendTo[res, ApplySetting[tree, settingName, #]] &/@ Utils`UnzipRules@settingName@$loopNumber;
   ];
   default = Switch[settingName,
      regularization, $regularizationScheme,
      momenta, Automatic,
      sum, {},
      mass, {}
   ];
   head = Switch[settingName,
      regularization|momenta|sum, First,
      mass, Identity
   ];
   res = res /.
      node[_?IsGeneric, __] -> default /.
         node[_?IsTopology, rest__] :> rest /.
            node[_?IsTreeHead, rest__] :> {rest};
   DeleteDuplicates/@Transpose@res /.
      {default, rest__} :> head@{rest} /. {default} -> default
];

ApplyObservableSetting // secure;

ApplySetting[tree_?IsTree, settingName:regularization|momenta, tQ_ -> res_] :=
Module[{},
   tree /. node[t:_?IsTopology /; tQ@t, rest__] :>
      (
         (*PrintInfo[tQ, settingName, ""];*)
         node[t, rest] /. node[g_?IsGeneric, __] :> res
      )
];

ApplySetting[tree_?IsTree, settingName:sum, tQ_ -> {str_String, {int_Integer, fun_}}] :=
Module[{res},
   res = {int -> Or[fun@TreeHead@tree, -fun@TreeHead@tree]};
   tree /. node[t:_?IsTopology /; tQ@t, rest__] :>
      (
         PrintInfo[tQ, settingName, str];
         node[t, rest] /. node[g_?IsGeneric, __] :> res
      )
];

ApplySetting[tree_?IsTree, settingName:mass, tQ_ -> {str_String, mass_ :> _[i_Integer]}] :=
Module[{res},
   res = With[{rhs = MassRules[][[i, 1, 1]]}, Append[#, mass :> rhs]&];
   tree /. node[t:_?IsTopology/; tQ@t, rest__] :>
      (
         PrintInfo[tQ, settingName, str];
         node[t, rest] /. node[g_?IsGeneric, __] :> res
      )
];

ApplySetting[tree_?IsTree, settingName:mass, tQ_ -> {str_String, Hold :> _[i_Integer]}] :=
Module[{res},
   res = With[{pos = i}, ReplacePart[#, pos -> {}]&];
   tree /. node[t:_?IsTopology/; tQ@t, rest__] :>
      (
         PrintInfo[tQ, settingName, str];
         node[t, rest] /. node[g_?IsGeneric, __] :>  res
      )
];

ApplySetting // secure;

PrintInfo[tQ_Symbol, settingName_Symbol, str_String, finalStr_String:""] :=
Module[{set, top, inf, n},
   set = "Setting " <> SymbolName@settingName <> ": ";
   top = SymbolName@tQ;
   inf = If[str == "", "", ": "<>str];
   n = Tr[StringLength/@{set, top, inf, finalStr}];
   Print[set, top, inf, StringJoin@Array[" " &, 70-n]<>finalStr];
];

PrintInfo // secure;

LoopFields[node[id_, ___], info__] := FeynArts`LoopFields[First@id, info];
LoopFields // secure;

TreeFields[node[id_, ___], info__] := FeynArts`TreeFields[First@id, info];
TreeFields // secure;

Field[d_?IsTopologyListHead, i_Integer] := Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]];
Field // secure;

FieldPattern[d_?IsTopologyListHead, i_Integer] := Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]] /. _?IsGenerationIndex :> Blank[];

FieldPattern[d_?IsTopologyListHead, a:HoldPattern@Alternatives@__] := FieldPattern[d, #] &/@ a;
FieldPattern // secure;

InternalMass[f_?IsFeynArtsField, index:_Integer] := FeynArts`Mass[f@FeynArts`Index[Generic, index], _?IsMass];
InternalMass // secure;

ExternalMass[index:_Integer] := {index};
ExternalMass // secure;

End[];
EndPackage[];
