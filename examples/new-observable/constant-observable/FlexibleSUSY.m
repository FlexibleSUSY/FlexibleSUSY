FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleConstantObservable, slha_, files_] :=
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
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@observables[[1]],
            "@prototype@" -> Observables`GetObservablePrototype@observables[[1]]
         }
      ];

      definitions = TextFormatting`ReplaceCXXTokens["
         @type@ @prototype@ {
            @type@ res {con};
            return res;
         }",
         {
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@observables[[1]],
            "@prototype@" -> Observables`GetObservablePrototype@observables[[1]]
         }
      ];
   ];

   (* Task 2: filling templates and moving them into models/Ma /observables/. *)
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
