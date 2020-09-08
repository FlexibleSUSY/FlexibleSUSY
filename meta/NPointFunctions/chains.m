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

BeginPackage["NPointFunctions`"];
Begin["NPointFunctions`internal`"];

simplifyChains::usage = "
@brief Simplifies some chains applying Dirac equation.
@param chain A chain to simplify.
@returns A simplified chain.";
simplifyChains[expr:_] := expr /. ch:FormCalc`DiracChain[__] :> simplifyChains@ch;
simplifyChains[chain:_FormCalc`DiracChain] :=
Module[{
      s = 6|7, a = -6|-7, ch = FormCalc`DiracChain, k = FormCalc`k,
      m, pair, sp, flip
   },
   m[FormCalc`Spinor[_, mass:_, type:_]] = type*mass;
   pair = FormCalc`Pair[k@#1,k@#2]&;
   sp[mom:_:_] = FormCalc`Spinor[k@mom, _, _];
   flip[7|-7] = 6;
   flip[6|-6] = 7;

   chain //. {
      ch[l:sp[j_],p:a,k[n_],k[i_],r:sp[n_]] :> pair[i,n]*ch[l,-p,r]-m[r]*ch[l,-p,k[i],r],
      ch[l:sp[n_],p:a,k[i_],k[n_],r:sp[j_]] :> pair[i,n]*ch[l,-p,r]-m[l]*ch[l,flip@p,k[i],r],
      ch[l:sp[],p:s,k[n_],r:sp[n_]] :> m[r]*ch[l,p,r],
      ch[l:sp[n_],p:s k[n_],r:sp[]] :> m[l]*ch[l,flip@p,r]
   }
];
simplifyChains // Utils`MakeUnknownInputDefinition;
simplifyChains ~ SetAttributes ~ {Protected, Locked};

modifyChains::usage = "
@brief Applies a set of process specific chain rules to some expression,
       depending on zeroMomentum option.
@param expression Any expression to modify.
@set A set of amplitudes or diagrams.
@zeroMomentum An value of ZeroExternalMomenta option.
@returns A modified expression.
@todo Currently a very specific set of rules is supported. In order to have
      chains on RHS of rules one needs to separate left and right parts and
      make different reveal functions.
@todo Write explanations about anticommutation rules in chains and other
      conventions.";
Module[{
      rules
   },
   modifyChains[
      expression_,
      set:`type`amplitudeSet|`type`diagramSet,
      zeroMomentum:_Symbol
   ] :=
   Expand@expression //.
   If[Head@rules === Symbol,
      Module[{
            i = 0, sp, L, reveal, ch = FormCalc`DiracChain
         },
         Block[{
            k = FormCalc`k,
            l = FormCalc`Lor
         },
            sp[mom_] := FormCalc`Spinor[k[mom], _, _];
            L[a_, e___ ,b_] := L[
               a,
               Switch[Length@{e},
                  0, {},
                  1, If[Head@Part[{e},1] === Integer, {e}, {6|7,e}],
                  _, If[Head@Part[{e},1] === Integer, {e}, {-6|-7,e}]
               ],
               b
            ];
            L[a_, {e___}, b_] := ch[sp@a, e, sp@b];
            reveal[{}] := Sequence[];
            reveal[{a_, b_, c___}] :=
               Flatten@{i++; i[e:___] :> L[a, e, b], reveal@{c}};
            chainRules = reveal@getFermionOrder@set;
            rules = zeroMomentum /. expandRules[`settings`chains] /. chainRules;
            If[Head@rules =!= List, rules = {}];
            rules
         ]
      ],
      rules
   ];
];
modifyChains // Utils`MakeUnknownInputDefinition;
modifyChains ~ SetAttributes ~ {Protected, Locked};

getChainRules::usage = "
@brief Finds a subset of rules inside a List, which represent Dirac chains. It
       is possible, because the naming convention for this abbreviation is fixed
       and it is given by encoded regular expression.
@param rules A list of rules.
@return A list of rules.";
getChainRules[rules:{Rule[_Symbol, _]...}] :=
Module[{
      regex = RegularExpression@"[F][1-9][\\d]*"
   },
   Cases[rules, e:Rule[_?(StringMatchQ[ToString@#, regex]&), _] :> e]
];
getChainRules // Utils`MakeUnknownInputDefinition;
getChainRules ~ SetAttributes ~ {Protected, Locked};

makeChainsUnique::usage = "
@brief After manual simplification of dirac chains one can get duplicates. They
       have to be removed. Then chains acquire unique names.
@param list A list of expression to modify and a list of rules, containing
            chains.
@returns A list of expression and rules.";
makeChainsUnique[list:{expression_, rules:{Rule[_Symbol, _]...}}] :=
Module[{
      chains = Longest@HoldPattern@Times[FormCalc`DiracChain[__]..],
      chain = FormCalc`DiracChain[__],
      name = Rule[Symbol["NPointFunctions`internal`dc"<>ToString@#1], #2]&,
      old, zero, rest, unique, erules
   },
   old = getChainRules@rules;
   zero = Cases[old, e:Rule[_, 0] :> e];
   rest = Complement[rules, old];
   old = Complement[old, zero];
   unique = foreach[name, DeleteDuplicates@Cases[old, chains, Infinity]];
   If[unique == {},
      unique = foreach[name, DeleteDuplicates@Cases[old, chain, Infinity]];
   ];

   erules = (old /. (unique /. Rule[x_, y_] :> Rule[y, x]));
   {
      expression /. zero /. erules,
      Join[unique, rest]
   } /. FormCalc`Mat -> NPointFunctions`internal`mat
];
makeChainsUnique // Utils`MakeUnknownInputDefinition;
makeChainsUnique ~ SetAttributes ~ {Protected, Locked};

mat::usage = "
@todo Sometimes this is needed, sometimes not. Why?";
mat[0] = 0;
mat[HoldPattern@Times[e__]] := Times@@mat/@{e};
mat[mass:_FeynArts`Mass] := mass;

identifySpinors::usage = "
@brief Inserts names of fermionic fields inside FormCalc`DicaChain structures.
@param rules List of rules to modify.
@param set A set of amplitudes or diagrams.
@returns Modified rules with inserted fermion names.
@note Dirac chains live only inside FormCalc`Abbr.
@note Should not be used with Automatic FormCalc`FermionOrder (later comment:
      why?).";
Module[{id, idf},
   identifySpinors[
      rules:{Rule[_Symbol, _]...},
      set:`type`amplitudeSet|`type`diagramSet
   ] :=
   (
      If[Head@id === Symbol,
         id = foreach[#1->#2&, getField[set, All] //. `rules`fieldNames];
         With[{
               ch = FormCalc`DiracChain, s = FormCalc`Spinor, k = FormCalc`k
            },
            idf[ch[s[k[i1_], m1_, _], e___, s[k[i2_], m2_, _]]] :=
               ch[s[i1 /. id, k[i1], m1], e, s[i2 /. id, k[i2], m2]];
         ];
      ];
      rules /. ch:FormCalc`DiracChain[__] :> idf@ch
   );

   identifySpinors // Utils`MakeUnknownInputDefinition;
   identifySpinors ~ SetAttributes ~ {Protected,Locked};
];

setZeroExternalMomentaInChains::usage = "
@brief Sets FormCalc`k[i] to zero inside fermionic chains.
@param expression Any expression.
@returns An expression with modified fermionic chains.";
setZeroExternalMomentaInChains[expression_] :=
   expression /. e:FormCalc`DiracChain[__] :> (e /. FormCalc`k[_] :> 0);
setZeroExternalMomentaInChains // Utils`MakeUnknownInputDefinition;
setZeroExternalMomentaInChains ~ SetAttributes ~ {Protected,Locked};

End[];
EndPackage[];
