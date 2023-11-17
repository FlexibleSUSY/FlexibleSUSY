Begin@"FlexibleSUSY`Private`";

FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`BrLTo3L, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      exportFields = {},
      npfDefinitions = "", obsPrototypes = "", obsDefinitions = "",
      ffvV = {}, npfV = {}, fermions = {}
   },

   If[observables =!= {} && FlexibleSUSY`FSFeynArtsAvailable && FlexibleSUSY`FSFormCalcAvailable,
   Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];

      exportFields = {#, #} &/@ observables[[All, 1]];
      fermions =     {SARAH`bar@#, #} &/@ observables[[All, 1]];
      ffvV = Flatten/@Tuples@{fermions, {TreeMasses`GetPhoton[]}};

      {npfV, npfDefinitions, obsPrototypes, obsDefinitions} = create@observables;
   ];

   WriteOut`ReplaceInFiles[files,
      {
         "@npointfunctions_headers@" -> NPointFunctions`CreateCXXHeaders[],
         "@npointfunctions_definitions@" -> npfDefinitions,
         "@calculate_prototypes@" -> obsPrototypes,
         "@calculate_definitions@" -> obsDefinitions,
         "@include_guard@" -> SymbolName@obs,
         "@namespace@" -> Observables`GetObservableNamespace@obs,
         "@filename@" -> Observables`GetObservableFileName@obs,
         "@get_MSUSY@" -> TextFormatting`IndentText@TextFormatting`WrapLines@AMM`AMMGetMSUSY[],
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];
   {
      "FFV fields" -> DeleteDuplicates@exportFields,
      "C++ vertices" -> DeleteDuplicates@Join[ffvV, npfV]
   }
];

create[manyObservables:{__FlexibleSUSYObservable`BrLTo3L}] := {
   DeleteDuplicates[Join@@#[[All,1]]],
   StringRiffle[#[[All, 2]], "\n\n"],
   StringRiffle[#[[All, 3]], "\n\n"],
   StringRiffle[#[[All, 4]], "\n\n"]
}&[create/@DeleteDuplicates[manyObservables /. Rule[_, {_, _, _}] -> Rule[_, {_, _, _}]]];

create[obs:FlexibleSUSYObservable`BrLTo3L[lep_, __, loopN_]] :=
Module[{npfVertices = {}, npfCode = "", prototype, definition, scalars, vectors, boxes},
   {{npfVertices, npfCode}, {scalars, vectors, boxes}} = generate@obs;

   prototype = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " " <> Observables`GetObservablePrototype@obs;

   definition = FixedPoint[
      StringReplace[#,
         {
            "@L@"   -> CConversion`ToValidCSymbolString@lep,
            "@Ph@"  -> CConversion`ToValidCSymbolString@SARAH`Photon,
            "@photon_penguin@" -> If[loopN === 1, "calculate_@L@_@L@_@Ph@_form_factors", "zero"],
            "@class_scalars@" -> scalars,
            "@class_vectors@" -> vectors,
            "@class_boxes@"   -> boxes
         }
      ]&,
      prototype <> " {
      return forge<
         fields::@L@, fields::@Ph@,
         @photon_penguin@,
         @class_scalars@,
         @class_vectors@,
         @class_boxes@
      >(nI, nO, nA, model, qedqcd);\n}"
   ];

   {npfVertices, npfCode, prototype <> ";", definition}
];

Module[{unevaluatedContributions = Tuples@{{0, 1}, {Scalars, Vectors, Boxes}}},
generate[obs:FlexibleSUSYObservable`BrLTo3L[lep_, _, con_, loopN_]] :=
Module[{
      parsed, keep, npfVertices = {}, npfCode = "",
      extraVertices, extraCode, forge, funName, npf, operatorRules
   },
   parsed = SymbolName/@If[Head@# === List, #, {#}]&@con;
   keep = Switch[{loopN, parsed},
      {0, {"All"}},       {Vectors, Scalars},
      {0, {"NoScalars"}}, {Vectors},
      {1, {"All"}},       {Vectors, Scalars, Boxes},
      {1, {"NoScalars"}}, {Vectors, Boxes},
      {1, {"Penguins"}},  {Vectors, Scalars},
      _, Symbol/@parsed
   ];

   If[Not@MemberQ[keep, #],
      forge@# = "zero",
      funName = ToLowerCase@SymbolName@# <>ToString@loopN<>"loop";
      forge@# = "npointfunctions::" <> funName;
      If[MemberQ[unevaluatedContributions, {loopN, #}],
         unevaluatedContributions = unevaluatedContributions /. {loopN, #} :> Sequence[];
         Print["Contribution: "<>SymbolName@#];
         Print["Loop level: ", loopN];
         npf = NPointFunctions`NPointFunction[{lep, lep}, {lep, lep},
            NPointFunctions`OnShellFlag -> True,
            NPointFunctions`UseCache -> True,
            NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
            NPointFunctions`KeepProcesses -> {#},
            NPointFunctions`LoopLevel -> loopN,
            NPointFunctions`Observable -> obs,
            NPointFunctions`Regularize -> FlexibleSUSY`FSRenormalizationScheme
         ];
         npf = npf /. {
            SARAH`sum[__] -> 0,
            LoopTools`B0i[i_, _, mm__] :> LoopTools`B0i[i, 0, mm],
            LoopTools`C0i[i_, Repeated[_, {3}], mm__] :> LoopTools`C0i[i, Sequence@@Array[0&, 3], mm],
            LoopTools`D0i[i_, Repeated[_, {6}], mm__] :> LoopTools`D0i[i, Sequence@@Array[0&, 6], mm]
         };

         {npf, operatorRules} = match@npf;
         extraVertices = NPointFunctions`VerticesForNPointFunction@npf;
         extraCode = NPointFunctions`CreateCXXFunctions[npf, funName, Identity, operatorRules][[2]];

         npfVertices = Join[npfVertices, extraVertices];
         npfCode = npfCode <> "\n" <> extraCode;
      ];
   ] &/@ {Scalars, Vectors, Boxes};

   {{npfVertices, npfCode}, Last/@DownValues@forge}
];
];

match[npf_] :=
Module[{fields, sp, dc, dim6},
   fields = Flatten@NPointFunctions`GetProcess@npf;
   sp[i_] := SARAH`DiracSpinor[fields[[i]], 0, 0];
   dc[a_, b__, c_] := NPointFunctions`DiracChain[sp@a, b, sp@c];
   dim6 = With[{l = SARAH`Lorentz, R = 6, L = 7},
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
   {WilsonCoeffs`InterfaceToMatching[npf, dim6], dim6}
];

End[];
