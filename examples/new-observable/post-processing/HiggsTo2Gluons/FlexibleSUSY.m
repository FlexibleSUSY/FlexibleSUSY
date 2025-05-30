FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`HiggsTo2Gluons, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables[slha], _obs],
      prototypes = {}, definitions = {}, npfDefinitions = {}, npfVertices = {}, npfHeaders = ""
   },

   If[observables =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName[obs] <> " class ..."];
      prototypes = TextFormatting`ReplaceCXXTokens[
         "@type@ @prototype@;",
         {
            "@type@"      -> CConversion`CreateCType[Observables`GetObservableType[#]],
            "@prototype@" -> Observables`GetObservablePrototype[#]
         }
      ] &/@ observables;

      (* Task 1: filling definitions. *)

      observables = DeleteDuplicates[observables /. f_[_Integer] -> f[_]];

      Module[{higgs, gluon, contr, npf, npfName},
         higgs  = #[[1, 1]];
         gluon  = #[[1, 2, 1]];
         contr  = #[[2]];
         npf = NPointFunctions`NPointFunction[
            {higgs},        (* Incoming particles *)
            {gluon, gluon}, (* Outgoing particles *)
            NPointFunctions`UseCache            -> False,
            NPointFunctions`OnShellFlag         -> True,
            NPointFunctions`ZeroExternalMomenta -> False,
            NPointFunctions`LoopLevel           -> 1,
            NPointFunctions`Regularize          -> FlexibleSUSY`FSRenormalizationScheme,
            NPointFunctions`Observable          -> obs[],
            NPointFunctions`KeepProcesses       -> If[Head[contr] === List, contr, {contr}]
         ];
         npf = NPointFunctions`ApplySubexpressions[npf];
         npf = npf /. {
            SARAH`sum[__, FormCalc`ec[2, l_] FormCalc`ec[3, l_] SARAH`g[__]] :>         "e2e3",
            SARAH`sum[__, FormCalc`ec[2, l_] SARAH`Mom[3, l_] SARAH`g[__]]   :>         "e2m3",
            SARAH`sum[__, FormCalc`ec[3, l_] SARAH`Mom[2, l_] SARAH`g[__]]   :>         "e3m2",
            FormCalc`Eps[FormCalc`ec[2], FormCalc`ec[3], SARAH`Mom[2], SARAH`Mom[3]] :> "eps",
            SARAH`sum[__, SARAH`Mom[2, l_] SARAH`Mom[3, l_] SARAH`g[__]] :> SARAH`Mass[higgs]^2/2
         };

         npf = WilsonCoeffs`InterfaceToMatching[npf, {"eps", "e2e3", "e2m3" "e3m2"}];

         npfName = TextFormatting`ReplaceCXXTokens[
            "@higgs@to2@gluon@@contr@",
            {
               "@higgs@"  -> SymbolName@higgs,
               "@gluon@"  -> SymbolName@gluon,
               "@contr@"  -> CConversion`ToValidCSymbolString@contr
            }
         ];
         AppendTo[npfDefinitions,
            NPointFunctions`NPFDefinitions[npf, npfName, SARAH`Delta, {"eps", "e2e3", "e2m3_e3m2"}]
         ];
         AppendTo[npfVertices, NPointFunctions`VerticesForNPointFunction@npf];

         AppendTo[definitions,
            TextFormatting`ReplaceCXXTokens["
               @type@ @prototype@ {
                  const auto npf = npointfunctions::@npf_name@(model, {}, {});
                  // Calculations go here
                  @type@ res {};
                  return res;
               }",
               {
                  "@type@"      -> CConversion`CreateCType@Observables`GetObservableType@#,
                  "@prototype@" -> Observables`GetObservablePrototype@#,
                  "@npf_name@"  -> npfName
               }
            ]
         ];
      ] &/@ observables;

      npfHeaders = NPointFunctions`CreateCXXHeaders[];
   ];

   (* Task 2: filling templates and moving them into models/Ma/observables/. *)
   WriteOut`ReplaceInFiles[
      files,
      {
         "@npf_headers@"           -> npfHeaders,
         "@npf_definitions@"       -> StringRiffle[DeleteDuplicates[npfDefinitions], "\n\n"],
         "@calculate_prototypes@"  -> StringRiffle[DeleteDuplicates[prototypes],     "\n\n"],
         "@calculate_definitions@" -> StringRiffle[DeleteDuplicates[definitions],    "\n\n"],
         "@include_guard@"         -> SymbolName[obs],
         "@namespace@"             -> Observables`GetObservableNamespace[obs],
         "@filename@"              -> Observables`GetObservableFileName[obs],
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];

   (* Task 3: returning something to the outside world. *)
   {
      "C++ vertices" -> Flatten[npfVertices, 1]
   }
];
