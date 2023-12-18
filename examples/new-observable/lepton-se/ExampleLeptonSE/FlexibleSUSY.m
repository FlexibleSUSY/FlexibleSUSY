FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleLeptonSE, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototypes = {}, definitions = {}, npfDefinitions = {}, cxxVertices = {}, npfHeaders = ""
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];
      prototypes = TextFormatting`ReplaceCXXTokens[
         "@type@ @prototype@;",
         {
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ] &/@ observables;

      (* Task 1: filling definitions. *)
      observables = DeleteDuplicates[observables /. f_@_Integer -> f@_];

      Module[{field, contr, npf, basis, name},
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
         name = "se_"<>CConversion`ToValidCSymbolString[contr];

         AppendTo[cxxVertices, NPointFunctions`VerticesForNPointFunction@npf];
         AppendTo[npfDefinitions, NPointFunctions`CreateCXXFunctions[npf, name, Identity, basis][[2]]];
         AppendTo[definitions,
            TextFormatting`ReplaceCXXTokens["
               @type@ @prototype@ {
                  const auto npf = npointfunctions::@name@(model, {gen, gen}, {});
                  return {0, npf[0], npf[1]};
               }",
               {
                  "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
                  "@prototype@" -> Observables`GetObservablePrototype@#,
                  "@name@"      -> name
               }
            ]
         ];
      ] &/@ observables;

      npfHeaders     = NPointFunctions`CreateCXXHeaders[];
   ];

   (* Task 2: filling templates and moving them into models/Ma/observables/. *)
   WriteOut`ReplaceInFiles[
      files,
      {
         "@npf_headers@"           -> npfHeaders,
         "@npf_definitions@"       -> StringRiffle[DeleteDuplicates[npfDefinitions], "\n\n"],
         "@calculate_prototypes@"  -> StringRiffle[DeleteDuplicates[prototypes],     "\n\n"],
         "@calculate_definitions@" -> StringRiffle[DeleteDuplicates[definitions],    "\n\n"],
         "@include_guard@"         -> SymbolName@obs,
         "@namespace@"             -> Observables`GetObservableNamespace@obs,
         "@filename@"              -> Observables`GetObservableFileName@obs,
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];

   (* Task 3: returning something to the outside world. *)
   {
      "C++ vertices" -> Flatten[cxxVertices, 1]
   }
];
