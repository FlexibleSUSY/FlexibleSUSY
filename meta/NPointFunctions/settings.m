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
Utils`DynamicInclude@"mass.m";

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

ApplyObservableSetting[tree_?IsTree, settingType:diagrams|amplitudes] :=
Module[{doPresent, doAbsent, absent, todos, res = tree},
   {doPresent, doAbsent} = If[Head@# === List, #, {}]&@ settingType[$loopNumber, #] &/@ {Present, Absent};
   {doPresent, doAbsent} = Utils`UnzipRules/@ {doPresent, doAbsent};
   {doPresent, doAbsent} = Rule[SymbolName@First@#, Last@#] &/@ #&/@ {doPresent, doAbsent};
   absent = Complement[First/@doAbsent, $expressionsToDerive];
   todos = DeleteDuplicates@Flatten[
      Join[$expressionsToDerive /. doPresent, absent /. doAbsent]
   ];
   todos = Select[todos, Head@# === Rule&];
   Set[res, applySetting[res, #, settingType]] &/@ todos;
   res
];

applySetting[tree_?IsTree, tQ_ -> {str_String, fun_}, diagrams|amplitudes] :=
   info[RemoveNode[tree, tQ, fun], str];

settings[order] :=
If[Head@order[] === List,
   order[],
   Reverse@Range@Tr@$externalFieldNumbers
];

settings[tree:_?IsTree, settings:regularization|momenta|sum|mass] :=
Module[{res = {tree}, default, head},
   If[Head@settings@$loopNumber === List,
      AppendTo[res, applySetting[tree, #]]&/@
         Utils`UnzipRules@settings@$loopNumber;
   ];
   default = Switch[settings,
      regularization, $regularizationScheme,
      momenta, Automatic,
      sum, {},
      mass, {}
   ];
   head = Switch[settings,
      regularization|momenta|sum, First,
      mass, Identity
   ];
   res = res /.
      node[_?IsGeneric, __] -> default /.
         node[_?IsTopology, rest__] :> rest /.
            node[type`head, rest__] :> {rest};
   DeleteDuplicates/@Transpose@res /.
      {default, rest__} :> head@{rest} /. {default} -> default
];

settings // secure;

(* -----v This is a generator function for 'applySetting'.                   *)
makeApply[pattern_, function:_Symbol] :=
(  Off@RuleDelayed::rhs;
   applySetting[tree:_?IsTree, pattern] :=
      Module[{once},
         once[arg_] := once@arg =
            If[# =!= {}, Print@@#]&@
               Cases[pattern, _String, Infinity, Heads -> True];

         tree /. node[t:_?IsTopology /; tQ@t, rest__] :>
            (once@_; node[t, rest] /. node[g_?IsGeneric, __] :>
               function[fun, g, t, head@tree])];
   On@RuleDelayed::rhs;
);

makeApply[tQ_ -> fun_, value];
value[val_, ___] := val;

makeApply[tQ_ -> {str_String, fun:{_Integer, _}}, restrict];
restrict[{int_, fun_}, __, head_] :=
   {int -> Or[fun[_, _, head], -fun[_, _, head]]};

applySetting[tree:_?IsTree, tQ_ -> {str_String, fun:{Append, _}}] :=
tree /. node[t:_?IsTopology/; tQ@t, rest__] :> (
   Print@str;
   node[t, rest] /. node[g_?IsGeneric, __] :>
   append[fun]
);

(* --v Append the result to the list of all mass rules [mass.m].             *)
(*                   v-----------v Head of expanded InternalMass.            *)
append[{Append, mass_FeynArts`Mass :> ExternalMass[i_Integer]}, ___] :=
With[{rhs = mass`rules[][[i, 1, 1]]}, Append[#, mass :> rhs]&];
append // secure;

applySetting[tree:_?IsTree, tQ_ -> {str_String, fun:{Hold, _}}] :=
tree /. node[t:_?IsTopology/; tQ@t, rest__] :> (
   Print@str;
   node[t, rest] /. node[g_?IsGeneric, __] :>
   hold[fun]
);

(*          v---------v Number of ext. particle.                             *)
(*          v---------v TODO(uukhas): replace.                               *)
hold[{Hold, {i_Integer}}, ___] :=
With[{pos = i}, ReplacePart[#, pos -> {}]&];
(*              ^---------^ First positions are for ext. particles [mass.m]. *)
hold // secure;

applySetting // secure;

(* Functions below are supposed to be used in [OBSERVABLE/settings.m].       *)

LoopFields[node[id_, ___], info__] :=
   FeynArts`LoopFields[First@id, info];

TreeFields[node[id_, ___], info__] :=
   FeynArts`TreeFields[First@id, info];

Field[d_?IsTopologyListHead, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]];

FieldPattern[d_?IsTopologyListHead, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]] /.
      type`generationIndex :> Blank[];
FieldPattern[d_?IsTopologyListHead, a:HoldPattern@Alternatives@__] :=
   FieldPattern[d, #] &/@ a;

InternalMass[f:type`field, index:_Integer] :=
   FeynArts`Mass[f@genericIndex@index, type`mass];

ExternalMass[index:_Integer] := {index};
(*                              ^-----^ TODO(uukhas): replace.               *)

LoopFields // secure;
TreeFields // secure;
Field // secure;
FieldPattern // secure;
InternalMass // secure;
ExternalMass // secure;
End[];
EndPackage[];
