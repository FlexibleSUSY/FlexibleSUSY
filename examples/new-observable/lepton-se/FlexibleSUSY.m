FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleLeptonSE, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototypes = "", definitions = "", npfHeaders = "", npfDefinitions = "",
      npfOutput, npfVertices
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];
      (* Task 1: combining prototypes and filling definitions. *)
      observables = DeleteDuplicates[observables /. f_@_Integer -> f@_];

      prototypes = TextFormatting`ReplaceCXXTokens[
         "@type@ @prototype@;",
         {
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ]&/@observables;

      npfOutput = Module[{field, contr, npf, basis, code, vertices},
         field = Head@First@#;
         contr = Last@#;
         npf = NPointFunctions`NPointFunction[
            {field}, (* Incoming particles *)
            {field}, (* Outgoing particles *)
            NPointFunctions`UseCache            -> False,
            NPointFunctions`OnShellFlag         -> True,
            NPointFunctions`ZeroExternalMomenta -> True,
            NPointFunctions`LoopLevel           -> 1,
            NPointFunctions`Regularize          -> FlexibleSUSY`FSRenormalizationScheme,
            NPointFunctions`Observable          -> obs[],
            NPointFunctions`KeepProcesses       -> If[Head@contr === List, contr, {contr}]
         ];

         basis =
         {
            "left_wilson"  -> NPointFunctions`DiracChain[SARAH`DiracSpinor[field[{SARAH`gt2}], 0, 0], 7, SARAH`DiracSpinor[field[{SARAH`gt1}], 0, 0]],
            "right_wilson" -> NPointFunctions`DiracChain[SARAH`DiracSpinor[field[{SARAH`gt2}], 0, 0], 6, SARAH`DiracSpinor[field[{SARAH`gt1}], 0, 0]]
         };

         npf = WilsonCoeffs`InterfaceToMatching[npf, basis];
         code = NPointFunctions`CreateCXXFunctions[npf, "se_"<>CConversion`ToValidCSymbolString[contr], Identity, basis][[2]];
         vertices = NPointFunctions`VerticesForNPointFunction@npf;
         {vertices, code}
      ]&/@observables;

      {npfVertices, npfDefinitions} = Transpose[npfOutput];

      definitions = TextFormatting`ReplaceCXXTokens["
         @type@ @prototype@ {
            const auto npf = npointfunctions::se_@npf_name@(model, {idx, idx}, {});
            return {0, npf[0], npf[1]};
         }",
         {
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#,
            "@npf_name@"  -> CConversion`ToValidCSymbolString[Last[#]]
         }
      ]&/@observables;

      prototypes     = StringRiffle[DeleteDuplicates@prototypes,     "\n\n"];
      definitions    = StringRiffle[DeleteDuplicates@definitions,    "\n\n"];
      npfHeaders     = NPointFunctions`CreateCXXHeaders[];
      npfDefinitions = StringRiffle[DeleteDuplicates@npfDefinitions, "\n\n"];
   ];

   (* Task 2: filling templates and moving them into models/Ma/observables/. *)
   WriteOut`ReplaceInFiles[
      files,
      {
         "@calculate_prototypes@"        -> prototypes,
         "@calculate_definitions@"       -> definitions,
         "@npointfunctions_headers@"     -> npfHeaders,
         "@npointfunctions_definitions@" -> npfDefinitions,
         "@include_guard@"               -> SymbolName@obs,
         "@namespace@"                   -> Observables`GetObservableNamespace@obs,
         "@filename@"                    -> Observables`GetObservableFileName@obs,
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];

   (* Task 3: returning something to the outside world. *)
   {
      "C++ vertices" -> Flatten[npfVertices, 1]
   }
];
