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

BeginPackage@"mass`";

rules::usage = "
@brief Sets or gets a set of rules to nullify masses of external particles.
@note Amplitudes are taken, because they do not have colour structures already.
@returns A list of list of rules to nullify masses of external particles.";

Begin@"`Private`";

externalMasses[set:type`fc`amplitudeSet] :=
   Flatten[List@@NPointFunctions`Private`process@set, 1][[All, 3]];
(*               ^-----------------------------^ TODO(uukhas)                *)
externalMasses[tree:type`tree] :=
   FeynArts`Mass[# /. -1 -> 1]&/@ NPointFunctions`Private`fields[tree, Flatten];
(*                                ^----------------------------^ TODO(uukhas)*)
externalMasses // tools`secure;

(*       v--v            v- FormCalc creates new masses - we also need them. *)
Module[{data},
   rules[tree:type`tree, fc:type`fc`amplitudeSet] :=
   data = Partition[
      RuleDelayed[#, 0]&/@ Riffle[externalMasses@tree, externalMasses@fc],
      2
   ];

   rules[] := (
      Utils`AssertOrQuit[Head@data =!= Symbol, rules::errNotSet];
      data);
];
rules::errNotSet = "Call mass`rules[...] first.";

rules // tools`secure;

End[];
Block[{$ContextPath}, EndPackage[]];
