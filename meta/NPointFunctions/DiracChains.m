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

ProceedChains[tree:_?IsTree, g:_] :=
Module[{abbr, subs, chains, generic},
   abbr = FormCalc`Abbr[] //. FormCalc`GenericList[];
   {chains, abbr} = {#, Complement[abbr, #]}&@ ExtractDiracChains@abbr;
   subs = FormCalc`Subexpr[] //. FormCalc`GenericList[] //. abbr;
   chains = SimplifyOnShellChains@chains;
   chains = ModifyDiracChains@chains;
   {generic, chains} = MakeChainsUnique@{g /. abbr, chains};
   chains = InsertFermionNames[tree, chains];
   If[chains =!= {},
      {generic, chains, subs},
      {generic, chains, Join[subs, abbr]}
   ]
];
ProceedChains // secure;

ModifyDiracChains[expression_] :=
Module[{i = 0, rules, sp, L, reveal},
   If[Not@MatchQ[chains@$loopNumber, {__Rule}],
      Return@expression
   ];

   Block[{k = FormCalc`k, l = FormCalc`Lor, ch = DiracChain},
      sp[mom_] := FormCalc`Spinor[k@mom, _, _];
      L[a_, e___ , b_] :=
         L[a,
            Switch[{e},
               {}, {},
               {_Integer, ___}, {e},
               {__}, {6|7|-6|-7, e}
            ],
         b
         ];
      L[a_, {e___}, b_] := ch[sp@a, e, sp@b];
      reveal@{a_, b_, c___} := Flatten@{i++; i[e___] :> L[a, e, b], reveal@{c}};
      reveal@{} := Sequence[];
      chainRules = reveal@GetObservableSetting@order;
      rules = $zeroExternalMomenta /. Utils`UnzipRules@chains@$loopNumber /. chainRules;
   ];
   Expand@expression //. rules
];
ModifyDiracChains // secure;

SimplifyOnShellChains[expr_] :=
If[$onShell,
   expr /. ch:DiracChain[__] :> SimplifyOnShellChains@ch,
   expr
];

SimplifyOnShellChains[chain_DiracChain] :=
If[$onShell,
   Module[{s = 6|7, a = -6|-7, ch = DiracChain, k = FormCalc`k, m, pair, sp, flip},
      m[FormCalc`Spinor[_, mass_, type_]] = type*mass;
      pair = FormCalc`Pair[k@#1, k@#2]&;
      sp[mom:_:_] = FormCalc`Spinor[k@mom, _, _];
      flip[7|-7] = 6;
      flip[6|-6] = 7;

      chain //.
      {
         ch[l:sp[j_],p:a,k[n_],k[i_],r:sp[n_]] :> pair[i,n]*ch[l,-p,r]-m[r]*ch[l,-p,k[i],r],
         ch[l:sp[n_],p:a,k[i_],k[n_],r:sp[j_]] :> pair[i,n]*ch[l,-p,r]-m[l]*ch[l,flip@p,k[i],r],
         ch[l:sp[],p:s,k[n_],r:sp[n_]] :>         m[r]*ch[l,p,r],
         ch[l:sp[n_],p:s,k[n_],r:sp[]] :>         m[l]*ch[l,flip@p,r]
      }
   ],
   chain
];
SimplifyOnShellChains // secure;

ExtractDiracChains[rules:{Rule[_Symbol, _]...}] :=
Module[{regex},
   regex = RegularExpression@"[F][1-9][\\d]*";
   Cases[rules, e:Rule[_?(StringMatchQ[ToString@#, regex]&), _] :> e]
];
ExtractDiracChains // secure;

MakeChainsUnique[list:{expression_, rules:{Rule[_Symbol, _]...}}] :=
Module[{chains, chain, name, old, zero, rest, unique, erules},
   chains = Longest@HoldPattern@Times[DiracChain[__]..];
   chain = DiracChain[__];
   name = Rule[Symbol["NPointFunctions`DiracChain"<>ToString@#2[[1]]], #1]&;
   old = ExtractDiracChains@rules;
   zero = Cases[old, e:Rule[_, 0] :> e];
   rest = Complement[rules, old];
   old = Complement[old, zero];
   unique = MapIndexed[name, DeleteDuplicates@Cases[old, chains, Infinity]];
   If[unique == {},
      unique = MapIndexed[name, DeleteDuplicates@Cases[old, chain, Infinity]];];
   erules = (old /. (unique /. Rule[x_, y_] :> Rule[y, x]));
   {
      expression /. zero /. erules,
      Join[unique, rest]
   } /. FormCalc`Mat -> Mat
];
MakeChainsUnique // secure;

Mat[0] = 0;
Mat[HoldPattern@Times[e__]] := Times@@Mat/@{e};
Mat[mass:_FeynArts`Mass] := mass;

InsertFermionNames[tree:_?IsTree, rules:{Rule[_Symbol, _]...}] :=
Module[{id, idf, ch = DiracChain, s = FormCalc`Spinor, k = FormCalc`k},
   id = FieldRules@MapIndexed[#2[[1]]->#1&, GetFields[tree, Flatten]];
   idf[ch[s[k[i1_], m1_, _], e___, s[k[i2_], m2_, _]]] :=
      ch[s[i1 /. id, k[i1], m1], e, s[i2 /. id, k[i2], m2]];
   rules /. ch:DiracChain[__] :> idf@ch
];
InsertFermionNames // secure;

ZeroMomentaInChains[expression_] :=
   expression /. e:DiracChain[__] :> (e /. FormCalc`k[_] :> 0);
ZeroMomentaInChains // secure;

End[];
EndPackage[];
