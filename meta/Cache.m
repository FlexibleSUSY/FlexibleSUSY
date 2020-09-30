(* ::Package:: *)

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

BeginPackage["Cache`", {"SARAH`", "Utils`"}];

GetVertex::usage = "
@brief Is used instead of SARAH`Vertex without additional options.
@param input A list of three or four elements without check, whether they are
       particles or not.
@returns An expression, which would be returned by SARAH`Vertex.";
GetVertex::errGetOnce = "
Date should be loaded once.";
GetVertex::errPutOnce = "
Data should be saved once.";
GetVertex::errNotDef = "
Variable `1` should be defined before initialization.";
GetVertex::errNoFile = "
File `1` does not exist, but required.";

Begin["`internal`"];

With[{AOQ = Utils`AssertOrQuit, NAME = "CachedVertices.m",
      M = SARAH`Mom, G = SARAH`g, B = SARAH`bar, PL = SARAH`PL, PR = SARAH`PR,
      L1 = SARAH`lt1, L2 = SARAH`lt2, L3 = SARAH`lt3, L4 = SARAH`lt4,
      GT = SARAH`getType, GG = SARAH`g[_, _], LP = SARAH`LorentzProduct,
      MD = HoldPattern[SARAH`Mom[_, _] - SARAH`Mom[_, _]]},

Module[{putOnce, defineRules,

      CheckInitialization, getOnce, dir, file3, file4, rules, impl,

      out, save = {}, vertexType, perm, indexRules,
      particlesC, orderC, vertexC, orderN, vertexN, idxCurrent,
      idxGeneric,

      DefineCanonical, CanonicalType, DefineSaved, IndexRules,
      OrderedIndices, ClearIndices, GetIndices, FillIndices, ZeroType,
      RuleType, CatchSU3, CatchZero, Lorentz, Swap, Fermion,
      NewIndices, PrepareToSave, RestoreIndices, Body, Momenta,


      fieldsN, rest, restoreRules, back
   },


   CheckInitialization[] := (
      AOQ[Head@getOnce === Symbol, GetVertex::errGetOnce];
      AOQ[Head@# =!= Symbol, GetVertex::errNotDef, #] &@ SA`subUnitaryCondition;
      AOQ[Head@# =!= Symbol, GetVertex::errNotDef, #] &@ SARAH`subDependences;
      AOQ[Head@# =!= Symbol, GetVertex::errNotDef, #] &@ SARAH`$sarahCurrentOutputMainDir;

      dir = FileNameJoin@{SARAH`$sarahCurrentOutputMainDir, "EWSB", "Vertices"};
      file3 = FileNameJoin@{dir, "VertexList3.m"};
      file4 = FileNameJoin@{dir, "VertexList4.m"};

      AOQ[FileExistsQ@file3, GetVertex::errNoFile, file3];
      AOQ[FileExistsQ@file4, GetVertex::errNoFile, file4]
   );

   GetVertex /: Get[GetVertex] :=
   If[CheckInitialization[],

      DefineCanonical /@ Get@file3;
      DefineCanonical /@ Get@file4;

      rules = ReleaseHold[(First /@ DownValues@impl) /. _[e : {__}] :> Rule[Sort@e, e]];

      If[FileExistsQ@#, DefineSaved /@ Get@#] &@ FileNameJoin@{dir, NAME};

      Unprotect@GetVertex;

      GetVertex[input : {Repeated[_, {2}]}] := SARAH`Vertex@input;

      GetVertex[input : {Repeated[_, {3, 4}]}] :=
      With[{particlesN = ClearIndices@input},
         {vertexType, out} = If[Head@# === List,
            #,

            vertexN = Simplify[SARAH`Vertex@particlesN] /. SA`subUnitaryCondition;
            PrepareToSave@particlesN;
            DefineSaved@Last@save
         ] &@ impl[particlesN];
         idxGeneric = OrderedIndices@First@out;
         idxCurrent = FillIndices[{idxGeneric, OrderedIndices@input}];
         restoreRules = RestoreIndices@{idxGeneric, idxCurrent};

         out = Switch[vertexType,
            Dynamic|Plus|Minus|PlusMinus, Lorentz@#,
            Identity|Ordering|Null|None, #
         ] & [out /. restoreRules];

         CatchZero@CatchSU3@out
      ];

      GetVertex /: Put[GetVertex] :=
      If[AOQ[Head@putOnce === Symbol, GetVertex::errPutOnce],

         If[save =!= {},
            If[FileExistsQ@#,
               Put[Sort@Join[Get@#, save], #];,
               Put[Sort@save, #]
            ] &@ FileNameJoin@{dir, NAME};
         ];
         putOnce = {};
      ];

      GetVertex // Utils`MakeUnknownInputDefinition;
      GetVertex ~ SetAttributes ~ {Protected, Locked};

      getOnce = {};
   ];
   Protect@GetVertex;

   PrepareToSave[particlesN_] := If[FreeQ[rules, Sort@particlesN],

      out = vertexN /. SARAH`subDependences // TrigReduce;
      AppendTo[save, {First@vertexN, ZeroType@out}];,

      particlesC = Sort[particlesN] /. rules;
      orderC = Ordering@particlesC;
      orderN = Ordering@particlesN;
      perm = orderC[[Ordering@orderN]];
      {vertexType, vertexC} = impl@particlesC;

      indexRules = MapThread[Rule, {IndexRules[vertexC, orderC], IndexRules[vertexN, orderN]}];

      AppendTo[save, {particlesC, perm, indexRules, RuleType[vertexC /. indexRules, vertexN]}];
   ];

   IndexRules[vertex_, order_] := GetIndices@Part[vertex, 1, order];

   OrderedIndices[expr : {__}] := GetIndices /@ expr;

   GetIndices[expr_] := Flatten@Cases[expr, e : {__} :> e, Infinity];

   FillIndices[idxs_] := MapThread[PadRight[#2, Length@#1, Drop[#1, Length@#2]] &, idxs];

   RestoreIndices[idxs_] := Flatten@MapThread[
      Thread@Rule[#1, PadRight[#2, Length@#1, Remove]] &, idxs
   ] /. Rule[_, Remove] :> (## &[]);

   ClearIndices[expr : {__}] := expr /. e_[{__}] :> e;

   NewIndices[ind_, {idxs___}] := ToExpression /@
      StringReplace[ToString/@{idxs}, DigitCharacter.. :> ToString@ind];

   DefineCanonical[v_] := With[{lhs = ClearIndices@Part[v, 1]},
      impl[lhs] = {CanonicalType@v, v /. SA`subUnitaryCondition};
   ];

   Fermion[{{e1_, PL}, {e2_, PR}}] :=
      {{-e2, PL}, {-e1, PR}};
   Fermion[{{e1_, LP[g1_, PL]}, {e2_, LP[g2_, PR]}}] :=
      {{-e2, LP[g2, PL]}, {-e1, LP[g1, PR]}};
   Fermion[{{e1_, LP[g1_, PL]}, {e2_, PR}}] :=
      {{-e2, PL}, {-e1, LP[g1, PR]}};

   Momenta[Plus, {{e_, d:MD}}] := {{e, d}};
   Momenta[Minus, {{e_, d:MD}}] := {{-e, -d}};

   Lorentz[{p_, {e_, d_M}}] :=
   Module[{ind, pos},
      back = Swap /@ restoreRules;
      ind = Flatten@Join[#[GT /@ p, SARAH`S | SARAH`V], #[Head /@ p, B]] &@ Position;
      ind = Complement[{1, 2, 3}, ind][[1]];
      {p, {e, d /. {i__} :> NewIndices[ind, {i} /. back]}}
   ];

   Lorentz[{p_, {e_, d:MD}}] :=
   Module[{newE, newD, asthere = Simplify[e*d]},
      If[FreeQ[restoreRules, _Integer],
         back = Swap /@ restoreRules;
         {newE, newD} = {Coefficient[asthere, #], #} &@ Cases[asthere, MD][[1]];

         If[-M[p[[1]] /. back, L3] + M[p[[2]] /. back, L3] === newD,
            {p, -{newE, newD}},
            {p, {newE, newD}}
         ],

         SARAH`Vertex[First@out /. restoreRules]
      ]

   ];

   Lorentz[{p_, {e_, d_ /; Count[d, GG, Infinity] === 3}}] :=
      Module[{newE = 1, newD = 1, asthere = Simplify[e*d], c12},
      Do[
         If[FreeQ[#, L1 | L2 | L3], newE = newE*#, newD = newD*#] &@ asthere[[i]],
         {i, Length@asthere}
      ];
      c12 = Coefficient[newD, G[L1, L2]];

      If[Head[DeleteCases[DeleteCases[c12, M[A_[{c___, L2, b___}], L3], 3], M[Susyno`LieGroups`conj[A_[{c___, L2, b___}]], L3], 3]] === Plus,
         {p, -{newE, newD}},
         {p, {newE, newD}}
      ]
   ];

   Lorentz[{p_, e : Repeated[{_, _}, {3}]}] := {p,
      Flatten@Cases[{e}, {_, G[L1, L2]*G[L3, L4]}],
      Flatten@Cases[{e}, {_, G[L1, L3]*G[L2, L4]}],
      Flatten@Cases[{e}, {_, G[L1, L4]*G[L2, L3]}]
   };

   CanonicalType[v_] := Switch[v,
      {_, {_, _M}}, Dynamic,
      {_, {_, MD}}, PlusMinus,
      {_, p_ /; Count[p, GG, Infinity] === 3}, Dynamic,
      {_, Repeated[{_, _}, {3}]}, Dynamic,
      _, None
   ];

   RuleType[{_, {_, _M}}, _] := Dynamic;
   RuleType[{_, {_, MD}}, _] :=
   If[Simplify@Last[vertexN/Body[particlesC, perm, indexRules]] === {-1, -1},
      Minus,
      Plus
   ];
   RuleType[{_, p_ /; Count[p, GG, Infinity] === 3}, _] := Dynamic;
   RuleType[{Repeated[_, {4}]}, _] := Dynamic;
   RuleType[vertexC_, vertexN_] :=
      If[TrigExpand@Expand@Rest[vertexC] === TrigExpand@Expand@Rest[vertexN], Identity, Ordering];

   CatchSU3[expr_] := expr /. { SARAH`fSU3[a_, a_, b_] -> 0, SARAH`fSU3[a_, b_, b_] -> 0};

   CatchZero[{p_}] := {p};
   CatchZero[{p_, {0, d_Integer}}] := {p, {0, d}};
   CatchZero[{p_, {0, d:MD}}] := {p, {0, d}};
   CatchZero[{p_, {0, _}}] := {p, {0, 1}};
   CatchZero[{p_, {0, _}, {0, _}}] := {p, {0, PL}, {0, PR}};
   CatchZero[{p_, Repeated[{}, {3}]}] := {p,
      {0, G[L1, L2]*G[L3, L4]},
      {0, G[L1, L3]*G[L2, L4]},
      {0, G[L1, L4]*G[L3, L2]}
   };
   CatchZero[{p_, v__}] := {p, v};

   ZeroType[{_}] := 0;
   ZeroType[{_, {0, _}}] := 1;
   ZeroType[{_, Repeated[{0, _}, {2}]}] := 2;
   ZeroType[{_, Repeated[{0, _}, {3}]}] := 3;

   DefineSaved[{particles_, 0}] :=
   With[{lhs = ClearIndices@particles},
      impl[lhs] = {Null, {particles}}
   ];

   DefineSaved[{particles_, 1}] :=
   With[{lhs = ClearIndices@particles},
      impl[lhs] = {Null, {particles, {0, 1}}}
   ];

   DefineSaved[{particles_, 2}] :=
   With[{lhs = ClearIndices@particles},
      impl[lhs] = {Null, {particles, {0, PL}, {0, PR}}}
   ];

   DefineSaved[{particles_, 3}] :=
   With[{lhs = ClearIndices@particles},
      impl[lhs] = {Null, {particles,
         {0, G[L1, L2]*G[L3, L4]},
         {0, G[L1, L3]*G[L2, L4]},
         {0, G[L1, L4]*G[L3, L2]}
      }}
   ];

   defineRules =
   {
      Ordering -> Fermion,
      Dynamic -> Identity,
      Plus -> (Momenta[Plus, #]&),
      Minus -> (Momenta[Minus, #]&)
   };

   DefineSaved[{particles_, ordering_, indexRules_,
      function:Identity|Ordering|Dynamic|Plus|Minus}] :=
   With[{
         lhs = Part[particles, ordering],
         rhs = (
            fieldsN = Part[First@# /. indexRules, ordering];
            rest = (function /. defineRules)[Rest@# /. indexRules];
            Prepend[rest, fieldsN]
         ) &@ Last@impl@particles
      },
      impl[lhs] = {function, rhs}
   ];

   Body[particles_, ordering_, indexRules_] := (
      fieldsN = Part[First@# /. indexRules, ordering];
      Prepend[Rest@# /. indexRules, fieldsN]
   ) &@ Last@impl@particles;

   Swap[x_ -> y_] := y -> x;

]; (* Module *)
]; (* With *)

End[];
EndPackage[];
