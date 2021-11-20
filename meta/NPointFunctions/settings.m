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

Utils`DynamicInclude@"type.m";
Utils`DynamicInclude@"mass.m";

BeginPackage@"NPointFunctions`";
Begin@"`Private`";

settings::usage = "
@brief Loads default topologies from topologies.m file.
       Loads data from OBSERVABLE/settings.m file.
       Parses the following settings:

       diagrams[LOOPLEVEL, Plus] | diagrams[LOOPLEVEL, Minus]
       amplitudes[LOOPLEVEL, Plus] | amplitudes[LOOPLEVEL, Minus]
       regularization[LOOPLEVEL]
          Amplitudes for some topologies are calculated with bugs in CDR.
          One can override used scheme for desired topologies.
       momenta[LOOPLEVEL]
          Eliminate specific momenta in some topologies.
       sum[LOOPLEVEL]
          Skip summation over some indices of particles.
       massless[LOOPLEVEL]
          Do not put some masses to zero.
       order[]
          Overrides default fermion order.";
With[{dir = DirectoryName@$InputFileName},
   settings[] :=
      (
         define@topologies;
         BeginPackage@"NPointFunctions`";
         Begin@"`Private`";
         If[FileExistsQ@#, Get@#;]&@FileNameJoin@
            {dir, `options`observable[], "settings.m"};
         End[];
         EndPackage[];);];

settings[tree:type`tree, settings:diagrams|amplitudes] :=
(*       ^--^ This object is modified and returned.                          *)
(*                       ^------^ Are defined in OBSERVABLE/settings.m file. *)
Module[{doPresent, doAbsent, absent, todos, res = tree},
   {doPresent, doAbsent} = If[Head@# === List, #, {}]&@
      settings[`options`loops[], #]&/@ {Plus, Minus};
   {doPresent, doAbsent} = tools`unzipRule/@ {doPresent, doAbsent};
   {doPresent, doAbsent} = Rule[SymbolName@First@#, Last@#]&/@ #&/@
      {doPresent, doAbsent};
   absent = Complement[First/@doAbsent, `options`processes[]];
   todos = DeleteDuplicates@Flatten[
      Join[`options`processes[] /. doPresent, absent /. doAbsent]
   ];
   todos = Select[todos, Head@# === Rule&];
   Set[res, applySetting[res, #]]&/@ todos;
   res
];

applySetting[tree:type`tree, tQ_ -> {str_String, fun_}] :=
(* --------^ For diagrams and amplitudes.                                    *)
   info[cut[tree, tQ, fun], str];

settings[order] :=
If[Head@order[] === List,
   order[],
   Reverse@Range@Tr@`options`observable@Outer
];

settings[tree:type`tree, settings:regularization|momenta|sum|massless] :=
Module[{res = {tree}, default, head},
   If[Head@settings@`options`loops[] === List,
      AppendTo[res, applySetting[tree, #]]&/@
         tools`unzipRule@settings@`options`loops[];
   ];
   default = Switch[settings,
      regularization, `options`scheme[],
      momenta, Automatic,
      sum, {},
      massless, {}
   ];
   head = Switch[settings,
      regularization|momenta|sum, First,
      massless, Identity
   ];
   res = res /.
      node[type`generic, __] -> default /.
      node[type`topology, rest__] :> rest /.
      node[type`head, rest__] :> {rest};
   DeleteDuplicates/@Transpose@res /.
      {default, rest__} :> head@{rest} /. {default} -> default
];

settings // tools`secure;


makeApply[pattern_, function:_Symbol] :=
(* -----^ This is a generator function for 'applySetting'.                   *)
(  Off@RuleDelayed::rhs;
   applySetting[tree:type`tree, pattern] :=
      Module[{once},
         once[arg_] := once@arg =
            If[# =!= {}, Print@@#]&@
               Cases[pattern, _String, Infinity, Heads -> True];

         tree /. node[t:type`topology /; tQ@t, rest__] :>
            (once@_; node[t, rest] /. node[g:type`generic, __] :>
               function[fun, g, t, head@tree])];
   On@RuleDelayed::rhs;
);

makeApply[tQ_ -> fun_, value];
value[val_, ___] := val;

makeApply[tQ_ -> {str_String, fun:{_Integer, _}}, restrict];
restrict[{int_, fun_}, __, head_] :=
   {int -> Or[fun[_, _, head], -fun[_, _, head]]};

(*                           v- Selects on which topology RHS is applied.    *)
(*                                   v-- Will be printed during evaluation.  *)
(*                                             v---------v Use in settings.m.*)
applySetting[tree:type`tree, tQ_ -> {str_String, fun:{Append, _}}] :=
tree /. node[t:type`topology/; tQ@t, rest__] :> (
   Print@str;
   node[t, rest] /. node[g:type`generic, __] :>
   append[fun, g, t, head@tree]
);
(*                             ^--^ For topology nodes, allowed by tQ.       *)
(*                       ^------------^ For all generic nodes.               *)
(* ^----^ Apply this function.                                               *)

(* --v Acts on the list, returned by mass`rules[].                           *)
append[{Append, (f:type`field)[n_Integer] :> e_Integer}, ___] :=
With[{rhs = (First/@ Flatten@mass`rules[])[[2*e-1]]},
   Append[#, genericMass[f, n] :> rhs]&
];
(* ^ Because rhs should be evaluated, but it is used with :>.                 *)

makeApply[tQ_ -> {str_String, fun:{Hold, _}}, hold];
hold[{Hold, e_}, ___] :=
   With[{pos = {{2*e}, {2*e-1}}}, Delete[#, pos]&];

(* Functions below are supposed to be used in OBSERVABLE/settings.m          *)

LoopFields[node[id_, ___], info__] :=
   FeynArts`LoopFields[First@id, info];

TreeFields[node[id_, ___], info__] :=
   FeynArts`TreeFields[First@id, info];

FieldPattern[d:Head@type`diagramSet, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]] /.
      type`generationIndex :> Blank[];
FieldPattern[d:Head@type`diagramSet, a:HoldPattern@Alternatives@__] :=
   FieldPattern[d, #] &/@ a;

applySetting // tools`secure;
LoopFields // tools`secure;
TreeFields // tools`secure;
FieldPattern // tools`secure;
End[];
EndPackage[];
