FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`ExampleFermionMass, slha_, files_] :=
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
      ]&/@observables;

      (* Task 1: filling definitions. *)
      definitions = TextFormatting`ReplaceCXXTokens["
         @type@ @prototype@ {
            return forge<@type@, fields::@fermion@>(gen, model, qedqcd);
         }",
         {
            "@fermion@"   -> SymbolName@Head@#[[1]],
            "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
            "@prototype@" -> Observables`GetObservablePrototype@#
         }
      ]&/@observables;
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
   {}
];
