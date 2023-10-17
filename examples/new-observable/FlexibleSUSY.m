FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`MyNewObservable, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      npfHeaders = "", npfDefinitions = "", prototypes = "", definitions = ""
   },

   If[observables =!= {},
      Print["Creating ", SymbolName@obs, " class ..."];

      prototypes = CConversion`CreateCType@Observables`GetObservableType@# <> " " <>
         Observables`GetObservablePrototype@# <> ";"&/@observables;

      definitions = CConversion`CreateCType@Observables`GetObservableType@# <> " " <>
         Observables`GetObservablePrototype@# <> " {\n" <>
         "   " <> CConversion`CreateCType@Observables`GetObservableType@# <> " res;\n" <>
         "   return res;\n}"&/@observables;
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@npointfunctions_headers@" -> npfHeaders,
         "@npointfunctions_definitions@" -> npfDefinitions,
         "@calculate_prototypes@" -> DeleteDuplicates@StringRiffle[prototypes, "\n"],
         "@calculate_definitions@" -> DeleteDuplicates@StringRiffle[definitions, "\n"],
         "@include_guard@" -> SymbolName@obs,
         "@namespace@" -> Observables`GetObservableNamespace@obs,
         "@filename@" -> Observables`GetObservableFileName@obs,
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];
   {}
];
