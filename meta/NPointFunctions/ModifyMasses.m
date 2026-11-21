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

BeginPackage@"NPointFunctions`";
Begin@"`Private`";

ExternalMasses[tree_?IsTree] := FeynArts`Mass[# /. -1 -> 1]&/@ GetFields[tree, Flatten];
ExternalMasses[set_?IsFormCalcSet] := Flatten[List@@GetProcess@set, 1][[All, 3]];
ExternalMasses // secure;


Module[{data},
MassRules[tree_?IsTree, fc_?IsFormCalcSet] :=
data = Partition[
   RuleDelayed[#, 0]&/@ Riffle[ExternalMasses@tree, ExternalMasses@fc],
   2
];
MassRules[] := (
   Utils`AssertOrQuit[Head@data =!= Symbol, MassRules::errNotSet];
   data
);
];
MassRules // secure;
MassRules::errNotSet = "Call MassRules[...] first.";

ModifyMasses[{generic_, chains_, subs_}, tree_?IsTree, NPointFunctions`ExceptLoops] :=
Module[{names, loops, uniqueIntegrals, hideInt, showInt, massRules, new},
   subWrite@"Applying subexpressions ... ";
   new = generic //. subs;
   subWrite@"done\n";

   names = ToExpression/@ Names@RegularExpression@"LoopTools`[ABCD]\\d+i*";
   loops = Alternatives@@ ( #[__] &/@ names );
   uniqueIntegrals = DeleteDuplicates@Cases[new, loops, Infinity];
   hideInt = Rule[#, Unique@"loopIntegral"] &/@ uniqueIntegrals;
   showInt = hideInt /. Rule[x_, y_] -> Rule[y, x];

   massRules = Through[(Composition@@ #&/@ ApplyObservableSetting[tree, mass])@MassRules[]];

   {
      MapThread[
         (#1//. Flatten@#2/. showInt)&,
         {new/. hideInt/. FormCalc`Pair[_, _]-> 0, massRules}
      ],
      ZeroMomentaInChains@chains /. Flatten@MassRules[],
      {}
   }
];

ModifyMasses[{generic_, chains_, subs_}, _, True] :=
Module[{zeroedRules, new},
   zeroedRules = Cases[subs, Rule[_, pair:FormCalc`Pair[_, _]] :> (pair->0)];
   {new, zeroedRules} = ZeroRules[subs, zeroedRules];
   {
      generic /. zeroedRules,
      ZeroMomentaInChains@chains,
      new
   }
];

ModifyMasses[e:{_, _, _}, __] := e;

ModifyMasses // secure;

ZeroRules[nonzeroRules:{Rule[_, _]...}, zeroRules:{Rule[_, 0]...}] :=
Module[{newNonzero, newZeroRules},
   newNonzero = Thread[
      Rule[nonzeroRules[[All,1]], nonzeroRules[[All,2]]/. zeroRules]];
   If[newNonzero === nonzeroRules, Return[{nonzeroRules, zeroRules}]];
   newZeroRules = Cases[newNonzero, HoldPattern[_-> 0]];
   newNonzero = Complement[newNonzero, newZeroRules];
   ZeroRules[newNonzero, Join[zeroRules, newZeroRules]]
];
ZeroRules // secure;

End[];
EndPackage[];
