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

settings::usage = "
@brief Converts `diagrams[_, _]` and `amplitudes[_, _]` from
       OBSERVABLE/settings.m files into operations, which remove diagrams
       with requested properties.
@returns A set of operations.";
settings[settings:diagrams|amplitudes] :=
Module[{doPresent, doAbsent, absent, todos},
   {doPresent, doAbsent} = If[Head@# === List, #, {}]&@
      settings[`options`loops[], #]&/@ {Plus, Minus};
   {doPresent, doAbsent} = tools`unzipRule/@ {doPresent, doAbsent};
   {doPresent, doAbsent} = Rule[SymbolName@First@#, Last@#]&/@ #&/@
      {doPresent, doAbsent};
   absent = Complement[First/@doAbsent, `options`processes[]];
   todos = DeleteDuplicates@Flatten[
      Join[`options`processes[] /. doPresent, absent /. doAbsent],
      1
   ];
   Select[todos, Head@# === List&]
];

settings[order] :=
If[Head@order[] === List,
   order[],
   Reverse@Range@Tr@`options`observable@Outer
];

settings[tree:type`tree, settings:regularization|momenta|sum|massless] :=
Module[{res = {tree}, default, head},
   If[Head@settings@`options`loops[] === List,
      AppendTo[res, applySetting[tree, #]]&/@
(*       v-----------------------v TODO: apply tools`unzip.                  *)
         settings@`options`loops[];
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

applySetting[tree:type`tree, {str_String, tQ_, fun_}] :=
   info[cut[tree, tQ, fun], str];

makeApply[pattern_, function:_Symbol] :=
   (  Off@RuleDelayed::rhs;
      applySetting[tree:type`tree, pattern] :=
         Module[{once},
            once[arg_] := once@arg =
               If[# =!= {}, Print@@#]&@
                  Cases[pattern, _String, Infinity, Heads -> True];
            tree /. node[t:type`topology /; tQ@t, rest__] :>
               (once@_; node[t, rest] /. node[g:type`generic, __] :>
                  function[fun, g, t, head@tree])];
      On@RuleDelayed::rhs;);

makeApply[tQ_ -> fun_, value];
makeApply[{str_String, tQ_, fun:{_Integer, _}}, restrict];
makeApply[{str_String, tQ_, fun:{Append, _}}, append];
makeApply[{str_String, tQ_, fun:{Hold, _}}, hold];
applySetting // tools`secure;

value[val_, ___] := val;
restrict[{int_, fun_}, __, head_] :=
   {int -> Or[fun[_, _, head], -fun[_, _, head]]};
append[{Append, (f:type`field)[n_Integer] :> e_Integer}, ___] :=
   With[{rhs = (First /@ getZeroMassRules[])[[2*e-1]]},
      Append[#, genericMass[f, n] :> rhs]&];
hold[{Hold, e_}, ___] :=
   With[{pos = {{2*e}, {2*e-1}}}, Delete[#, pos]&];

LoopFields[node[id_, ___], info__] :=
   FeynArts`LoopFields[First@id, info];
LoopFields // tools`secure;

TreeFields[node[id_, ___], info__] :=
   FeynArts`TreeFields[First@id, info];
TreeFields // tools`secure;

FieldPattern[d:Head@type`diagramSet, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]] /.
      type`generationIndex :> Blank[];
FieldPattern[d:Head@type`diagramSet, a:HoldPattern@Alternatives@__] :=
   FieldPattern[d, #] &/@ a;
FieldPattern // tools`secure;

End[];
EndPackage[];
