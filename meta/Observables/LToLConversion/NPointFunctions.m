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

BeginPackage@"LToLConversion`";
create::usage = "";
Begin@"`Private`";

create[obs:_[in_@_ -> _, _, con_, loopN_]] :=
Module[{npfVertices, npfHeader, npfDefinition, calculateDefinition, prototype},
   {npfVertices, npfHeader, npfDefinition} = generate@obs;

   prototype = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " " <> Observables`GetObservablePrototype@obs;

   {calculateDefinition, npfDefinition} = FixedPoint[
      StringReplace[#,
         {
            "@L@"   -> CConversion`ToValidCSymbolString@in,
            "@U@"   -> CConversion`ToValidCSymbolString@SARAH`UpQuark,
            "@D@"   -> CConversion`ToValidCSymbolString@SARAH`DownQuark,
            "@Ph@"  -> CConversion`ToValidCSymbolString@SARAH`Photon,
            "@N@"   -> CConversion`ToValidCSymbolString@loopN,
            "@con@" -> CConversion`ToValidCSymbolString@con,
            "@photon_penguin@" -> If[loopN === 1, "calculate_@L@_@L@_@Ph@_form_factors", "zero"],
            "@classU@" -> "conversion_@L@@U@_to_@L@@U@_@con@@N@loop",
            "@classD@" -> "conversion_@L@@D@_to_@L@@D@_@con@@N@loop"
         }
      ]&,
      {
         prototype <> " {
         return forge_conversion<
            fields::@L@, fields::@U@, fields::@D@, fields::@Ph@,
            @photon_penguin@,
            npointfunctions::@classU@,
            npointfunctions::@classD@
         >(in, out, n, model, parameters, qedqcd);\n}",
         npfDefinition
      }
   ];

   {
      npfVertices,
      {npfHeader, npfDefinition},
      {prototype <> ";", calculateDefinition}
   }
];

create[manyObservables_List] :=
Module[{unique},
   unique = DeleteDuplicates[manyObservables /. f_@_Integer -> f@_];
   {
      DeleteDuplicates[Join@@#[[All,1]]],
      {#[[1,2,1]], StringRiffle[#[[All,2,2]], "\n\n"]},
      {StringRiffle[#[[All,3,1]], "\n\n"], StringRiffle[#[[All,3,2]], "\n\n"]}
   }&[create/@unique]
];

parseSynonyms[_[__, con_, loopN_]] :=
Module[{parsed, result},
   parsed = SymbolName/@If[Head@# === List, #, {#}]&@con;
   result = Switch[{loopN, parsed},
      {0, {"All"}},       {Vectors, Scalars},
      {0, {"NoScalars"}}, {Vectors},
      {1, {"All"}},       {Vectors, Scalars, Boxes},
      {1, {"NoScalars"}}, {Vectors, Boxes},
      {1, {"Penguins"}},  {Vectors, Scalars},
      _, Symbol/@parsed
   ];
   Print["Contributions: ", SymbolName/@result];
   Print["Loop level: ", loopN];
   result
];

cleanLeftovers[npf_] := npf /. {
   SARAH`sum[__] -> 0,
   LoopTools`B0i[i_, _, mm__] :> LoopTools`B0i[i, 0, mm],
   LoopTools`C0i[i_, Repeated[_, {3}], mm__] :> LoopTools`C0i[i, Sequence@@Array[0&, 3], mm],
   LoopTools`D0i[i_, Repeated[_, {6}], mm__] :> LoopTools`D0i[i, Sequence@@Array[0&, 6], mm]
};

generate[obs:_[in_@_ -> _, __, loopN_]] :=
Module[{npfU, npfD, fields, keep, dim6, codeU, codeD},
   Utils`FSFancyLine[];
   keep = parseSynonyms@obs;
   {npfU, npfD} = NPointFunctions`NPointFunction[
      {in, #}, {in, #},
      NPointFunctions`OnShellFlag -> True,
      NPointFunctions`UseCache -> False,
      NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
      NPointFunctions`KeepProcesses -> keep,
      NPointFunctions`LoopLevel -> loopN,
      NPointFunctions`Observable -> obs] &/@ {SARAH`UpQuark, SARAH`DownQuark};
   {npfU, npfD} = cleanLeftovers/@{npfU, npfD};

   fields[SARAH`UpQuark] = Flatten@NPointFunctions`GetProcess@npfU;
   fields[SARAH`DownQuark] = Flatten@NPointFunctions`GetProcess@npfD;
   dim6[q_] := Module[{sp, dc, l = SARAH`Lorentz, R = 6, L = 7},
      sp[f_, n_] := SARAH`DiracSpinor[fields[f][[n]], 0, 0];
      dc[a_, b__, c_] := NPointFunctions`DiracChain[sp[q, a], b, sp[q, c]];
      {
         "S_LL" -> dc[3,L,1] dc[4,L,2],
         "S_LR" -> dc[3,L,1] dc[4,R,2],
         "S_RL" -> dc[3,R,1] dc[4,L,2],
         "S_RR" -> dc[3,R,1] dc[4,R,2],
         "V_LL" -> dc[3,R,l@1,1] dc[4,R,l@1,2],
         "V_LR" -> dc[3,R,l@1,1] dc[4,L,l@1,2],
         "V_RL" -> dc[3,L,l@1,1] dc[4,R,l@1,2],
         "V_RR" -> dc[3,L,l@1,1] dc[4,L,l@1,2],
         "T_LL" -> dc[3,-L,l@1,l@2,1] dc[4,-L,l@1,l@2,2],
         "T_RR" -> dc[3,-R,l@1,l@2,1] dc[4,-R,l@1,l@2,2]
      }
   ];
   npfU = WilsonCoeffs`InterfaceToMatching[npfU, dim6@SARAH`UpQuark];
   npfD = WilsonCoeffs`InterfaceToMatching[npfD, dim6@SARAH`DownQuark];

   codeU = NPointFunctions`CreateCXXFunctions[npfU, "@classU@", SARAH`Delta, dim6@SARAH`UpQuark][[2]];
   codeD = NPointFunctions`CreateCXXFunctions[npfD, "@classD@", SARAH`Delta, dim6@SARAH`DownQuark][[2]];
   Utils`FSFancyLine[];

   {
      DeleteDuplicates@Join[
         NPointFunctions`VerticesForNPointFunction@npfU,
         NPointFunctions`VerticesForNPointFunction@npfD
      ],
      NPointFunctions`CreateCXXHeaders[],
      codeU<>"\n\n"<>codeD
   }
];

End[];
Block[{$ContextPath}, EndPackage[]];
