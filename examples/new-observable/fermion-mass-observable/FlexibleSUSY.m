FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleFermionMass, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototypes = "", definitions = "", npfHeaders = "", npfDefinitions = ""
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];
      (* Task 1: combining prototypes and filling definitions. *)
      prototypes = TextFormatting`ReplaceCXXTokens[
         "@type@ @prototype@;",
         {
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ]&/@observables;

      definitions = TextFormatting`ReplaceCXXTokens["
         @type@ @prototype@ {
            return forge<@type@, fields::@fermion@>(idx, model, qedqcd);
         }",
         {
            "@fermion@"   -> SymbolName@Head@#[[1]],
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ]&/@observables;

      prototypes  = StringRiffle[DeleteDuplicates@prototypes,  "\n\n"];
      definitions = StringRiffle[DeleteDuplicates@definitions, "\n\n"];
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
   {}
];
