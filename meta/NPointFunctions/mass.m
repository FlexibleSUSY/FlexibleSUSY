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

Utils`DynamicInclude@"tools.m";
Utils`DynamicInclude@"type.m";

BeginPackage@"mass`";

rules::usage = "
@brief Sets or gets a set of rules to nullify masses of external particles.
@note Amplitudes are taken, because they do not have colour structures already.
@returns A list of list of rules to nullify masses of external particles.";

modify::usage = "
@brief Sets the masses of external particles to zero according to
       $zeroExternalMomenta.";

Begin@"`Private`";

externalMasses[set:type`fc`amplitudeSet] :=
   Flatten[List@@NPointFunctions`Private`process@set, 1][[All, 3]];

externalMasses[tree:type`tree] :=
   FeynArts`Mass[# /. -1 -> 1]&/@ NPointFunctions`Private`fields[tree, Flatten];

externalMasses // secure;

(*                       v- FormCalc creates new masses - we also need them. *)
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
rules // secure;
rules::errNotSet = "Call mass`rules[...] first.";

modify[{generic_, chains_, subs_}, tree:type`tree, NPointFunctions`ExceptLoops] :=
Module[{names, loops, uniqueIntegrals, hideInt, showInt, massRules, new},
   tools`subWrite@"Applying subexpressions ... ";
   new = generic //. subs;
   tools`subWrite@"done\n";

   names = ToExpression/@ Names@RegularExpression@"LoopTools`[ABCD]\\d+i*";
   loops = Alternatives@@ ( #[__] &/@ names );
   uniqueIntegrals = DeleteDuplicates@Cases[new, loops, Infinity];
   hideInt = Rule[#, Unique@"loopIntegral"] &/@ uniqueIntegrals;
   showInt = hideInt /. Rule[x_, y_] -> Rule[y, x];

   massRules = Through[(Composition@@ #&/@
      NPointFunctions`Private`settings[tree, NPointFunctions`Private`mass])@rules[]];

   {
      MapThread[
         (#1//. Flatten@#2/. showInt)&,
         {new/. hideInt/. FormCalc`Pair[_, _]-> 0, massRules}
      ],
      NPointFunctions`Private`zeroMomenta@chains /. Flatten@rules[],
      {}
   }
];

modify[{generic_, chains_, subs_}, _, True] :=
Module[{zeroedRules, new},
   zeroedRules = Cases[subs, Rule[_, pair:FormCalc`Pair[_, _]] :> (pair->0)];
   {new, zeroedRules} = tools`zeroRules[subs, zeroedRules];
   {
      generic /. zeroedRules,
      NPointFunctions`Private`zeroMomenta@chains,
      new
   }
];

modify[e:{_, _, _}, __] := e;

modify // secure;

End[];
Block[{$ContextPath}, EndPackage[]];
