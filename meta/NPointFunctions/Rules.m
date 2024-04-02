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

Module[{once},
FieldData[] := once@_;
once[arg_] := once@arg =
Module[{regex, lines, rules, names},
   regex = "(\\w+): ([SFVU])\\[(\\d+)\\]";
   lines = Utils`ReadLinesInFile@$particleNamesFile;
   rules = Rule[First@#, Sequence[Last@#, First@#]] &/@ Get@$particleNamespaceFile;
   names = StringCases[lines, RegularExpression@regex :> {"$1","FeynArts`$2","$3"}];
   Flatten[names, 1] /. rules
];
FieldData // secure;
];

Module[{once},
FieldRules[expression_] := expression //. FieldRules[];
FieldRules[] := once@_;
once[arg_] := once@arg =
Module[{bose = FeynArts`S|FeynArts`V, fermi = FeynArts`U|FeynArts`F, data},
   data = Map[ToExpression, {#1<>#2, #3, #4}&@@#&/@FieldData[], 2];
   Join[
      data /. {n_, t_, i_} :> Rule[t@i, n],
      data /. {n_, t_, i_} :> RuleDelayed[t[i, {ind__}], n@{ind}],
      data /.
      {
         {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n], Susyno`LieGroups`conj@n],
         {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n], SARAH`bar@n]
      },
      data /.
      {
         {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n@{ind__}], Susyno`LieGroups`conj@n@{ind}],
         {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n@{ind__}], SARAH`bar@n@{ind}]
      },
      {
         ind_?IsGenerationIndex :> Symbol["SARAH`gt" <> ToString@Last@ind],
         (* Color structures are handled separately *)
         (h_)[{_?IsColorIndex}] :> h(*Symbol["SARAH`ct" <> ToString@Last@ind]*),
         (h_)[{_?IsGluonIndex}] :> h(*Symbol["SARAH`ct" <> ToString@Last@ind]*)
      },
      {
         FeynArts`S -> GenericS,
         FeynArts`F -> GenericF,
         FeynArts`V -> GenericV,
         FeynArts`U -> GenericU
      },
      {
         Times[-1,field:_GenericS|_GenericV] :> Susyno`LieGroups`conj@field,
         Times[-1,field:_GenericF|_GenericU] :> SARAH`bar@field
      }
   ]
];
FieldRules // secure;
];

Module[{once},
SubexpressionMassRules[] := once@_;
once[arg_] := once@arg =
Module[{faMasses, sarahNames, massRules},
   faMasses = Symbol["Global`Mass" <> #[[2]]] &/@ FieldData[];
   sarahNames = Symbol[#[[1]] <> #[[2]]] &/@ FieldData[];
   massRules = MapThread[
      {
         #1[index_] :> SARAH`Mass@#2@{Symbol["SARAH`gt" <> StringTake[SymbolName@index, -1]]},
         #1[indices__] :> SARAH`Mass@#2@indices,
         #1 :> SARAH`Mass@#2
      }&,
      {faMasses, sarahNames}
   ];
   Append[Flatten@massRules, FeynArts`Mass[field_, _ : Null] :> SARAH`Mass@field]
];
SubexpressionMassRules // secure;
];

SubexpressionCouplingRules[] :=
Module[{PL, PR, MT, FV, g, md, v, i, p},
   v = Global`FourVector;
   p = FeynArts`Mom;
   i = FeynArts`KI1;
   l = NPointFunctions`LorentzIndex;
   {PL, PR} = FeynArts`NonCommutative@Global`ChiralityProjector@#|
      FeynArts`NonCommutative[Global`DiracMatrix@i@3,
         Global`ChiralityProjector@#] &/@ {-1, 1};
   Off@RuleDelayed::rhs;
   MT[i1_, i2_] := Global`MetricTensor[i@i1_Integer, i@i2_Integer];
   If[Head@FormCalc`$FormCalc === Real,
      FV[i1_, i2_, Repeated[_, {0, 1}]] := p@i1_Integer - p@i2_Integer;,
      FV[i1_, i2_] := v[p@i1_Integer - p@i2_Integer, i@3];
      FV[i1_, i2_, i3_] := v[p@i1_Integer - p@i2_Integer, i@i3_Integer];
   ];
   On@RuleDelayed::rhs;
   g[f_, i1:_Integer, i2:_Integer] := SARAH`g[l@f[[i1]], l@f[[i2]]];
   md[f_, i1:_Integer, i2:_Integer] := SARAH`Mom@f[[i1]] - SARAH`Mom@f[[i2]];
   md[f_, i1:_Integer, i2:_Integer, i3:_Integer] :=
      SARAH`Mom[f[[i1]], l@f[[i3]]] - SARAH`Mom[f[[i2]], l@f[[i3]]];

   With[{F = FeynArts`G[_][0][f__], S = SARAH`Cp@f},
      {
         F@1 :> S@1,
         F@PL :> S@SARAH`PL,
         F@PR :> S@SARAH`PR,
         F@MT[i1, i2] :> S@g[{f}, i1, i2],
         F@FV[i1, i2] :> S@md[{f}, i1, i2],
         F[FV[i2, i1, i3]*MT[i1, i2]+FV[i1, i3, i2]*MT[i1, i3]+FV[i3, i2, i1]*MT[i2, i3]] :>
         S[
           md[{f}, i2, i1, i3] * g[{f}, i1, i2],
           md[{f}, i1, i3, i2] * g[{f}, i1, i3],
           md[{f}, i3, i2, i1] * g[{f}, i2, i3]
         ]
      }
   ]
];
SubexpressionCouplingRules // secure;

With[{lt = Unique@"SARAH`lt"},
SubexpressionRules[expression_] :=
expression //. Join[
   FieldRules[],
   SubexpressionMassRules[],
   SubexpressionCouplingRules[],
   {
      FormCalc`Finite -> 1,
      FormCalc`Den[a_, b_] :> 1/(a-b),
      FormCalc`Pair[a_, b_] :> SARAH`sum[lt, 1, 4, SARAH`g[lt, lt]*Append[a, lt]*Append[b, lt]],
      f_?IsGenericField :> Head[f][GenericIndex@Last@Last@f],
      FormCalc`k[i_Integer, pairIndex___] :> SARAH`Mom[i, pairIndex]
   },
   {
      FormCalc`Spinor -> SARAH`DiracSpinor,
      FormCalc`Lor -> SARAH`Lorentz
   }
];
SubexpressionRules // secure;
];

AmplitudeRules[expression_] :=
SubexpressionRules[expression] //.
{
   FeynArts`SumOver[_, _, FeynArts`External] :> Sequence[],
   Times[e_, FeynArts`SumOver[i_, max_]] :> SARAH`sum[i, 1, max, e],
   Times[e_, FeynArts`SumOver[i, {min, max}]] :> SARAH`sum[i, min, max, e],
   SARAH`sum[i_, min_, max_, FeynArts`SumOver[_, n1_]] :> SARAH`sum[i, min, max, n1],
   SARAH`sum[i_, min_, max_, FeynArts`SumOver[_, {n1_, n2_}]] :> SARAH`sum[i, min, max, n2-n1],
   FeynArts`IndexSum -> Sum
};
AmplitudeRules // secure;

ExternalMomentaRules[tree_?IsTree, option:True|False|OperatorsOnly|ExceptLoops] :=
Module[{fsFields},
   Switch[option,
      True,
         {SARAH`Mom[_Integer,_] :> 0},
      OperatorsOnly|ExceptLoops,
         fsFields = FieldRules@GetFields[tree, Flatten];
         {SARAH`Mom[i_Integer, lorIndex_] :> SARAH`Mom[fsFields[[i]], lorIndex]},
      False,
         {}
   ]
];
ExternalMomentaRules // secure;

End[];
EndPackage[];
