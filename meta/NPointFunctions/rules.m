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
   fieldData[] := once@_;
   once[arg_] := once@arg =
      Module[{regex, lines, rules, names},
         regex = "(\\w+): ([SFVU])\\[(\\d+)\\]";
         lines = Utils`ReadLinesInFile@`file`particles[];
         rules = Rule[First@#, Sequence[Last@#, First@#]] &/@
            Get@`file`contexts[];
         names = StringCases[lines, RegularExpression@regex :>
            {"$1","FeynArts`$2","$3"}];
         Flatten[names, 1] /. rules];];
fieldData // secure;

Module[{once},
   `rules`fields[] := once@_;
   once[arg_] := once@arg =
      Module[{bose = FeynArts`S|FeynArts`V, fermi = FeynArts`U|FeynArts`F,
            data, generation},
         data = Map[ToExpression, {#1<>#2, #3, #4}&@@#&/@fieldData[], 2];
         generation =
            Module[{filter, rules, indices},
               filter = (FeynArts`Indices -> e_) :> e;
               rules = {  FeynArts`Index -> Identity,
                          Global`Colour :> {},
                          Global`Gluon :> {}};
               indices = Cases[FeynArts`M$ClassesDescription,
                  filter,
                  Infinity];
               indices = DeleteDuplicates@Flatten[indices //. rules];
               FeynArts`Index[Alternatives@@indices, _Integer]];
         Join[
            data /. {n_, t_, i_} :>
               Rule[t@i, n],
            data /. {n_, t_, i_} :>
               RuleDelayed[t[i, {ind__}], n@{ind}],
            data /. {
               {n_, t:bose,  _} :>
                  RuleDelayed[Times[-1,f:n], Susyno`LieGroups`conj@n],
               {n_, t:fermi, _} :>
                  RuleDelayed[Times[-1,f:n], SARAH`bar@n]},
            data /.
            {  {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n@{ind__}],
                  Susyno`LieGroups`conj@n@{ind}],
               {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n@{ind__}],
                  SARAH`bar@n@{ind}]},
            {  ind:generation :>
                  Symbol["SARAH`gt" <> ToString@Last@ind],
               ind:`type`indexCol :>
                  Symbol["SARAH`ct" <> ToString@Last@ind],
               ind:`type`indexGlu :>
                  Symbol["SARAH`ct" <> ToString@Last@ind]},
            {  FeynArts`S -> GenericS,
               FeynArts`F -> GenericF,
               FeynArts`V -> GenericV,
               FeynArts`U -> GenericU},
            {  Times[-1,field:_GenericS|_GenericV] :>
                  Susyno`LieGroups`conj@field,
               Times[-1,field:_GenericF|_GenericU] :>
                  SARAH`bar@field}]];];
`rules`fields // secure;

Module[{once},
   `rules`mass[] := once@_;
   once[arg_] := once@arg =
      Module[{faMasses, sarahNames, massRules},
         faMasses = Symbol["Global`Mass" <> #[[2]]] &/@ fieldData[];
         sarahNames = Symbol[#[[1]] <> #[[2]]] &/@ fieldData[];
         massRules = MapThread[
            {  #1[index_] :> SARAH`Mass@#2@{Symbol["SARAH`gt" <>
                  StringTake[SymbolName@index, -1]]},
               #1[indices__] :> SARAH`Mass@#2@indices,
               #1 :> SARAH`Mass@#2} &,
            {faMasses, sarahNames}];
         Append[Flatten@massRules,
            FeynArts`Mass[field_, _ : Null] :> SARAH`Mass@field]];];
`rules`mass // secure;

`rules`couplings[] :=
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
   If[And@@MapThread[Less, {FormCalc`$FormCalc, {9, 7}}],
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
`rules`couplings // secure;

With[{lt = Unique@"SARAH`lt"},
   `rules`subexpressions[expression_] :=
      expression //. Join[
         `rules`fields[],
         `rules`mass[],
         `rules`couplings[],
         {  FormCalc`Finite -> 1,
            FormCalc`Den[a_, b_] :> 1/(a-b),
            FormCalc`Pair[a_, b_] :> SARAH`sum[lt, 1, 4,
               SARAH`g[lt, lt]*Append[a, lt]*Append[b, lt]],
            f:`type`genericField :> Head[f][GenericIndex@Last@Last@f],
            FormCalc`k[i_Integer, pairIndex___] :> SARAH`Mom[i, pairIndex]},
         {  FormCalc`Spinor -> SARAH`DiracSpinor,
            FormCalc`Lor -> SARAH`Lorentz}];];
`rules`subexpressions // secure;

`rules`amplitude[expression_] :=
   `rules`subexpressions[expression] //. {
      FeynArts`SumOver[_, _, FeynArts`External] :> Sequence[],
      Times[e_, FeynArts`SumOver[i_, max_]] :>
         SARAH`sum[i, 1, max, e],
      Times[e_, FeynArts`SumOver[i, {min, max}]] :>
         SARAH`sum[i, min, max, e],
      SARAH`sum[i_, min_, max_, FeynArts`SumOver[_, n1_]] :>
         SARAH`sum[i, min, max, n1],
      SARAH`sum[i_, min_, max_, FeynArts`SumOver[_, {n1_, n2_}]] :>
         SARAH`sum[i, min, max, n2-n1],
      FeynArts`IndexSum -> Sum};
`rules`amplitude // secure;

End[];
EndPackage[];
