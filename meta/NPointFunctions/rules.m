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
Begin@"`internal`";

$FieldNames =
Module[{regex, lines, rules, names},
   regex = "(\\w+): ([SFVU])\\[(\\d+)\\]";
   lines = Utils`ReadLinesInFile@$ParticleFile;
   rules = Rule[First@#, Sequence[Last@#, First@#]] &/@ Get@$ContextFile;
   names = StringCases[lines, RegularExpression@regex :>
      {"$1","FeynArts`$2","$3"}];
   Flatten[names, 1] /. rules];
$FieldNames // Protect;

$FieldRules =
Module[{bose = FeynArts`S|FeynArts`V, fermi = FeynArts`U|FeynArts`F},
   Join[
      # /. {n_, t_, i_} :> Rule[t@i, n],
      # /. {n_, t_, i_} :> RuleDelayed[t[i, {ind__}], n@{ind}],
      # /.
      {  {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n], Susyno`LieGroups`conj@n],
         {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n], SARAH`bar@n]},
      # /.
      {  {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n@{ind__}],
            Susyno`LieGroups`conj@n@{ind}],
         {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n@{ind__}],
            SARAH`bar@n@{ind}]},
      {  ind:`type`indexGeneration :> Symbol["SARAH`gt" <> ToString@Last@ind],
         ind:`type`indexCol :> Symbol["SARAH`ct" <> ToString@Last@ind],
         ind:`type`indexGlu :> Symbol["SARAH`ct" <> ToString@Last@ind]},
      {  FeynArts`S -> GenericS,
         FeynArts`F -> GenericF,
         FeynArts`V -> GenericV,
         FeynArts`U -> GenericU},
      {  Times[-1,field:_GenericS|_GenericV] :> Susyno`LieGroups`conj@field,
         Times[-1,field:_GenericF|_GenericU] :> SARAH`bar@field}]&@
         Map[ToExpression, {#[[1]]<>#[[2]], #[[3]], #[[4]]}&/@$FieldNames, 2]];
$FieldRules // Protect;

$MassRules =
Module[{faMasses, sarahNames, massRules},
   faMasses = Symbol["Global`Mass" <> #[[2]]] &/@ $FieldNames;
   sarahNames = Symbol[#[[1]] <> #[[2]]] &/@ $FieldNames;
   massRules = MapThread[
      {  #1[index_] :> SARAH`Mass@#2@{Symbol["SARAH`gt" <>
            StringTake[SymbolName@index, -1]]},
         #1[indices__] :> SARAH`Mass@#2@indices,
         #1 :> SARAH`Mass@#2} &,
      {faMasses, sarahNames}];
   Append[Flatten@massRules,
      FeynArts`Mass[field_, _ : Null] :> SARAH`Mass@field]];
$MassRules // Protect;

$CouplingRules =
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
   If[FormCalc`$FormCalc < 9.7,
      FV[i1_, i2_, Repeated[_, {0, 1}]] := p@i1_Integer - p@i2_Integer;,
      FV[i1_, i2_] := v[p@i1_Integer - p@i2_Integer, i@3];
      FV[i1_, i2_, i3_] := v[p@i1_Integer - p@i2_Integer, i@i3_Integer];];
   On@RuleDelayed::rhs;
   g[f_, i1:_Integer, i2:_Integer] := SARAH`g[l@f[[i1]], l@f[[i2]]];
   md[f_, i1:_Integer, i2:_Integer] := SARAH`Mom@f[[i1]] - SARAH`Mom@f[[i2]];
   md[f_, i1:_Integer, i2:_Integer, i3:_Integer] :=
      SARAH`Mom[f[[i1]], l@f[[i3]]] - SARAH`Mom[f[[i2]], l@f[[i3]]];

   With[{F = FeynArts`G[_][0][f__], S = SARAH`Cp@f},
      {  F@1 :> S@1,
         F@PL :> S@SARAH`PL,
         F@PR :> S@SARAH`PR,
         F@MT[i1, i2] :> S@g[{f}, i1, i2],
         F@FV[i1, i2] :> S@md[{f}, i1, i2],
         F[ FV[i2, i1, i3] * MT[i1, i2] +
            FV[i1, i3, i2] * MT[i1, i3] +
            FV[i3, i2, i1] * MT[i2, i3]] :>
         S[ md[{f}, i2, i1, i3] * g[{f}, i1, i2],
            md[{f}, i1, i3, i2] * g[{f}, i1, i3],
            md[{f}, i3, i2, i1] * g[{f}, i2, i3]]}]];
$CouplingRules // Protect;

$SubexpressionRules =
Join[$FieldRules, $MassRules, $CouplingRules,
   {  FormCalc`Finite -> 1,
      FormCalc`Den[a:_,b:_] :> 1/(a-b),
      FormCalc`Pair[a:_,b:_] :> SARAH`sum[#, 1, 4,
         SARAH`g[#, #]*Append[a, #]*Append[b, #]],
      f:`type`genericField :> Head[f][GenericIndex@Last@Last@f],
      FormCalc`k[i:_Integer, pairIndex:___] :> SARAH`Mom[i, pairIndex]} &@
         Unique@"SARAH`lt"];
$SubexpressionRules // Protect;

$AmplitudeRules = Join[
   $SubexpressionRules,
   {  FeynArts`SumOver[_,_,FeynArts`External] :> Sequence[],
      Times[e:_, FeynArts`SumOver[i:_Symbol, max:_Integer]] :>
         SARAH`sum[i, 1, max, e],
      Times[expr:_, FeynArts`SumOver[index:_Symbol,
         {min:_Integer, max:_Integer}]] :>
         SARAH`sum[index, min, max, expr],
      SARAH`sum[i:_Symbol, _Integer, max:_Integer,
         FeynArts`SumOver[_Symbol, max2:_Integer]] :>
         SARAH`sum[i, 1, max, max2],
      SARAH`sum[i:_Symbol, _Integer, max:_Integer,
         FeynArts`SumOver[_, {min2:_Integer, max2:_Integer}]] :>
         SARAH`sum[i, 1, max, max2-min2]},
   {FeynArts`IndexSum -> Sum}];
$AmplitudeRules // Protect;

End[];
EndPackage[];
