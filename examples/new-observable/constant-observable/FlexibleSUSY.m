FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleConstantObservable, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototypes = "", definitions = ""
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];

      prototypes =
         CConversion`CreateCType@Observables`GetObservableType@# <> " " <>
         Observables`GetObservablePrototype@# <> ";"&/@observables;

      definitions =
         CConversion`CreateCType@Observables`GetObservableType@# <> " " <>
         Observables`GetObservablePrototype@# <> " {\n" <>
         "   " <> CConversion`CreateCType@Observables`GetObservableType@# <> " res {"<> ToString@#[[1]]<> "};\n" <>
         "   return res;\n}"&/@observables;
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
