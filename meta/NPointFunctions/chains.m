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

proceedChains[d:`type`diagramSet, a:`type`amplitudeSet, g:_] :=
Module[{abbr, subs, chains, generic},
   abbr = FormCalc`Abbr[] //. FormCalc`GenericList[];
   {chains, abbr} = {#, Complement[abbr, #]}&@ getChainRules@abbr;
   subs = FormCalc`Subexpr[] //. FormCalc`GenericList[] //. abbr;
   chains = simplifyChains@chains;
   chains = modifyChains[chains, d];
   {generic, chains} = makeChainsUnique@{g /. abbr, chains};
   chains = identifySpinors[chains, a];
   {generic, chains, subs}];
proceedChains // secure;

modifyChains::usage = "
@brief Transforms `settings`chains into rules and applies them onto expression.
       `settings`chains should have the following syntax:
          {  {  <symbol>..} ->  {  <rule>..}
             ..}
       <symbol> A symbol, which chooses (based on $ZeroMomenta) set of rules.
                Put different <symbol> in the same list, if you want to use the
                same rules for them.
       <rule> A rule with the following syntax:
                 <chain>[<something1>, <something2>...] -> 0
              <chain> A number of chain. Number is defined by `settings`order,
                      i.e. for {3, 1, 4, 2} the chain 1 is {3, 1} and chain 2
                      is {4, 2}.
              <something1> 6, -6, 7, -7, k@<integer>, l@<integer>. In the case
                           of last two, alternatives of first four is prepended.
              <something2> k@<integer>, l@<integer>, _, __, ___.
              <integer> Any integer number.
@param expression Any expression to modify.
@param set A set of diagrams.
@returns A modified expression.
@todo Currently a very specific set of rules is supported. In order to have
      chains on RHS of rules one needs to separate left and right parts and
      make different reveal functions.
@todo Write explanations about anticommutation rules in chains and other
      conventions.";
modifyChains[expression_, set:`type`diagramSet] :=
Module[{i = 0, rules, sp, L, reveal},
   If[`settings`chains === Default, Return@expression];
   Block[{k = FormCalc`k, l = FormCalc`Lor, ch = FormCalc`DiracChain},
      sp[mom_] := FormCalc`Spinor[k@mom, _, _];
      L[a_, e___ , b_] := L[ a,
         Switch[{e},
            {}, {}, {_Integer, ___}, {e}, {__}, {6|7|-6|-7, e}],
         b];
      L[a_, {e___}, b_] := ch[sp@a, e, sp@b];
      reveal@{a_, b_, c___} := Flatten@{i++; i[e:___] :> L[a, e, b], reveal@{c}};
      reveal@{} := Sequence[];
      chainRules = reveal@getFermionOrder@set;
      rules = $ZeroMomenta /. expandRules@`settings`chains /. chainRules;];
   Expand@expression //. rules];
modifyChains // secure;

simplifyChains::usage = "
@brief Simplifies some chains applying Dirac equation if $OnShell is True.
@param chain A chain to simplify.
@returns A simplified chain.";
If[$OnShell,
   simplifyChains[expr:_] :=
      (expr /. ch:FormCalc`DiracChain[__] :> simplifyChains@ch);
   simplifyChains[chain:_FormCalc`DiracChain] :=
   Module[{s = 6|7, a = -6|-7, ch = FormCalc`DiracChain, k = FormCalc`k,
         m, pair, sp, flip},
      m[FormCalc`Spinor[_, mass:_, type:_]] = type*mass;
      pair = FormCalc`Pair[k@#1,k@#2]&;
      sp[mom:_:_] = FormCalc`Spinor[k@mom, _, _];
      flip[7|-7] = 6;
      flip[6|-6] = 7;

      chain //. {
         ch[l:sp[j_],p:a,k[n_],k[i_],r:sp[n_]] :>
            pair[i,n]*ch[l,-p,r]-m[r]*ch[l,-p,k[i],r],
         ch[l:sp[n_],p:a,k[i_],k[n_],r:sp[j_]] :>
            pair[i,n]*ch[l,-p,r]-m[l]*ch[l,flip@p,k[i],r],
         ch[l:sp[],p:s,k[n_],r:sp[n_]] :> m[r]*ch[l,p,r],
         ch[l:sp[n_],p:s,k[n_],r:sp[]] :> m[l]*ch[l,flip@p,r]}];,
   simplifyChains[expr:_] := expr;];
simplifyChains // secure;

expandRules::usage = "
@brief Expands a set of compact rules into the full one.
@param rules A set of compact rules to expand.
@returns A set of rules.";
expandRules[rules:{Rule[{__Symbol}, _]..}] :=
   rules /. Rule[e:_, s:_] :> Sequence @@ (Rule[#, s] &/@ e);
expandRules // secure;

getChainRules::usage = "
@brief Finds a subset of rules inside a List, which represent Dirac chains. It
       is possible, because the naming convention for this abbreviation is fixed
       and it is given by encoded regular expression.
@param rules A list of rules.
@return A list of rules.";
getChainRules[rules:{Rule[_Symbol, _]...}] :=
Module[{regex},
   regex = RegularExpression@"[F][1-9][\\d]*";
   Cases[rules, e:Rule[_?(StringMatchQ[ToString@#, regex]&), _] :> e]];
getChainRules // secure;

makeChainsUnique::usage = "
@brief After manual simplification of dirac chains one can get duplicates. They
       have to be removed. Then chains acquire unique names.
@param list A list of expression to modify and a list of rules, containing
       chains.
@returns A list of expression and rules.";
makeChainsUnique[list:{expression_, rules:{Rule[_Symbol, _]...}}] :=
Module[{chains, chain, name, old, zero, rest, unique, erules},
   chains = Longest@HoldPattern@Times[FormCalc`DiracChain[__]..];
   chain = FormCalc`DiracChain[__];
   name = Rule[Symbol["NPointFunctions`internal`dc"<>ToString@#2[[1]]], #1]&;
   old = getChainRules@rules;
   zero = Cases[old, e:Rule[_, 0] :> e];
   rest = Complement[rules, old];
   old = Complement[old, zero];
   unique = MapIndexed[name, DeleteDuplicates@Cases[old, chains, Infinity]];
   If[unique == {},
      unique = MapIndexed[name, DeleteDuplicates@Cases[old, chain, Infinity]];];
   erules = (old /. (unique /. Rule[x_, y_] :> Rule[y, x]));
   {  expression /. zero /. erules,
      Join[unique, rest]} /. FormCalc`Mat -> NPointFunctions`internal`mat];
makeChainsUnique // secure;

mat::usage = "
@todo Sometimes this is needed, sometimes not. Why?";
mat[0] = 0;
mat[HoldPattern@Times[e__]] := Times@@mat/@{e};
mat[mass:_FeynArts`Mass] := mass;

identifySpinors::usage = "
@brief Inserts names of fermionic fields inside FormCalc`DicaChain structures.
@param rules List of rules to modify.
@param set A set of amplitudes.
@returns Modified rules with inserted fermion names.";
identifySpinors[rules:{Rule[_Symbol, _]...}, set:`type`amplitudeSet] :=
Module[{id, idf, ch = FormCalc`DiracChain, s = FormCalc`Spinor, k = FormCalc`k},
   id = MapIndexed[#2[[1]]->#1&, getField[set, All] //. getFieldRules[]];
   idf[ch[s[k[i1_], m1_, _], e___, s[k[i2_], m2_, _]]] :=
      ch[s[i1 /. id, k[i1], m1], e, s[i2 /. id, k[i2], m2]];
   rules /. ch:FormCalc`DiracChain[__] :> idf@ch];
identifySpinors // secure;

setZeroExternalMomentaInChains::usage = "
@brief Sets FormCalc`k[i] to zero inside fermionic chains.
@param expression Any expression.
@returns An expression with modified fermionic chains.";
setZeroExternalMomentaInChains[expression_] :=
   expression /. e:FormCalc`DiracChain[__] :> (e /. FormCalc`k[_] :> 0);
setZeroExternalMomentaInChains // secure;

End[];
EndPackage[];