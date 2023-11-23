FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleConstantObservable, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototype = "", definition = ""
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];

      prototype = TextFormatting`ReplaceCXXTokens["@type@ @prototype@;",
         {
            "@type@" -> CConversion`CreateCType@Observables`GetObservableType@observables[[1]],
            "@prototype@" -> Observables`GetObservablePrototype@observables[[1]]
         }
      ];

      definition = TextFormatting`ReplaceCXXTokens["
         @type@ @prototype@ {
            @type@ res {con};
            return res;
         }",
         {
            "@type@" -> CConversion`CreateCType@Observables`GetObservableType@observables[[1]],
            "@prototype@" -> Observables`GetObservablePrototype@observables[[1]]
         }
      ];
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@calculate_prototypes@"  -> prototype,
         "@calculate_definitions@" -> definition,
         "@include_guard@"         -> SymbolName@obs,
         "@namespace@"             -> Observables`GetObservableNamespace@obs,
         "@filename@"              -> Observables`GetObservableFileName@obs,
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];
   {}
];
