FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleFermionMass, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototypes = "", definitions = ""
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];

      prototypes = TextFormatting`ReplaceCXXTokens["@type@ @prototype@;",
         {
            "@type@" -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ]&/@observables;

      definitions = TextFormatting`ReplaceCXXTokens["
         @type@ @prototype@ {
            return forge<@type@, fields::@fermion@>(idx, model, qedqcd);
         }",
         {
            "@fermion@" -> SymbolName@Head@#[[1]],
            "@type@" -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ]&/@observables;
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@calculate_prototypes@"  -> StringRiffle[DeleteDuplicates@prototypes, "\n\n"],
         "@calculate_definitions@" -> StringRiffle[DeleteDuplicates@definitions, "\n\n"],
         "@include_guard@"         -> SymbolName@obs,
         "@namespace@"             -> Observables`GetObservableNamespace@obs,
         "@filename@"              -> Observables`GetObservableFileName@obs,
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];
   {}
];
