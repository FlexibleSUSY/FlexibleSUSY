Begin@"FlexibleSUSY`Private`";

WriteClass[obs:FlexibleSUSYObservable`LToLConversion, slha_, files_] :=
Module[{
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      fieldsFFV = {}, verticesFFV = {}, additionalVertices = {},
      prototypes = "", npfHeaders = "", definitions = "", npfDefinitions = "",
      newRules
   },

   newRules = {
      "@LToLConversion_fill@" -> "slha_io.fill(ltolconversion_settings);",
      "@LToLConversion_init@" -> "LToLConversion_settings ltolconversion_settings;",
      "@LToLConversion_class_name@" -> "ltolconversion_settings,",
      "@LToLConversion_class_declaration@" -> "class LToLConversion_settings;",
      "@LToLConversion_named_argument@" -> "const LToLConversion_settings& ltolconversion_settings,",
      (* Live inside librarylink. *)
      "@LToLConversion_private@" -> "LToLConversion_settings ltolconversion_settings{}; ///< LToLConversion settings",
      "@LToLConversion_setter@" -> "void set_ltolconversion_settings(const LToLConversion_settings& s) { ltolconversion_settings = s; }",
      "@LToLConversion_set_data@" -> "data.set_ltolconversion_settings(ltolconversion_settings);",
      "@LToLConversion_set_slha@" -> "slha_io.set_LToLConversion_settings(ltolconversion_settings);",
      "@LToLConversion_reset@" -> "ltolconversion_settings.reset();"
   };
   If[observables === {} || !FlexibleSUSY`FSFeynArtsAvailable || !FlexibleSUSY`FSFormCalcAvailable,
      newRules = newRules /. Rule[x_, _]:> Rule[x, ""];
   ];

   If[observables =!= {} && FlexibleSUSY`FSFeynArtsAvailable && FlexibleSUSY`FSFormCalcAvailable,
      Print["Creating ", SymbolName@obs, " class ..."];

      fieldsFFV = {#, #} &/@ observables[[All, 1]];
      verticesFFV = Flatten/@Tuples@{{SARAH`bar@#, #}&/@TreeMasses`GetSMQuarks[], {TreeMasses`GetPhoton[]}};

      {additionalVertices, npfDefinitions, prototypes, definitions} = create@observables;
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@npointfunctions_headers@" -> NPointFunctions`CreateCXXHeaders[],
         "@npointfunctions_definitions@" -> npfDefinitions,
         "@calculate_prototypes@" -> prototypes,
         "@calculate_definitions@" -> definitions,
         "@include_guard@" -> SymbolName@obs,
         "@namespace@" -> Observables`GetObservableNamespace@obs,
         "@filename@" -> Observables`GetObservableFileName@obs,
         "@get_MSUSY@" -> TextFormatting`IndentText@TextFormatting`WrapLines@AMM`AMMGetMSUSY[],
         Sequence@@GeneralReplacementRules[]
      }
   ];
   {
      "FFV fields" -> DeleteDuplicates@fieldsFFV,
      "C++ vertices" -> DeleteDuplicates@Join[verticesFFV, additionalVertices],
      "C++ replacements" -> newRules
   }
];

create[manyObservables:{__FlexibleSUSYObservable`LToLConversion}] := {
   DeleteDuplicates[Join@@#[[All,1]]],
   StringRiffle[#[[All, 2]], "\n\n"],
   StringRiffle[#[[All, 3]], "\n\n"],
   StringRiffle[#[[All, 4]], "\n\n"]
}&[create/@DeleteDuplicates[manyObservables /. Rule[_Integer, _Integer] -> Rule[_, _]]];

create[obs:FlexibleSUSYObservable`LToLConversion[in_, __, con_, loopN_]] :=
Module[{npfVertices, npfDefinition, calculateDefinition, prototype},
   {npfVertices, npfDefinition} = generate@obs;

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

   {npfVertices, npfDefinition, prototype <> ";", calculateDefinition}
];

generate[obs:FlexibleSUSYObservable`LToLConversion[in_, __, con_, loopN_]] :=
Module[{npfU, npfD, fields, keep, dim6, codeU, codeD, parsed},
   Utils`FSFancyLine[];
   parsed = SymbolName/@If[Head@# === List, #, {#}]&@con;
   keep = Switch[{loopN, parsed},
      {0, {"All"}},       {Vectors, Scalars},
      {0, {"NoScalars"}}, {Vectors},
      {1, {"All"}},       {Vectors, Scalars, Boxes},
      {1, {"NoScalars"}}, {Vectors, Boxes},
      {1, {"Penguins"}},  {Vectors, Scalars},
      _, Symbol/@parsed
   ];
   Print["Contributions: ", SymbolName/@keep];
   Print["Loop level: ", loopN];

   {npfU, npfD} = NPointFunctions`NPointFunction[
      {in, #}, {in, #},
      NPointFunctions`OnShellFlag -> True,
      NPointFunctions`UseCache -> False,
      NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
      NPointFunctions`KeepProcesses -> keep,
      NPointFunctions`LoopLevel -> loopN,
      NPointFunctions`Observable -> obs,
      NPointFunctions`Regularize -> FlexibleSUSY`FSRenormalizationScheme
   ] &/@ {SARAH`UpQuark, SARAH`DownQuark};
   {npfU, npfD} = {npfU, npfD} /. {
      SARAH`sum[__] -> 0,
      LoopTools`B0i[i_, _, mm__] :> LoopTools`B0i[i, 0, mm],
      LoopTools`C0i[i_, Repeated[_, {3}], mm__] :> LoopTools`C0i[i, Sequence@@Array[0&, 3], mm],
      LoopTools`D0i[i_, Repeated[_, {6}], mm__] :> LoopTools`D0i[i, Sequence@@Array[0&, 6], mm]
   };

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
      codeU<>"\n\n"<>codeD
   }
];

End[];
