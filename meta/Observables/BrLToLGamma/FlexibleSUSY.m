Begin@"FlexibleSUSY`Private`";

WriteClass[obs:FlexibleSUSYObservable`BrLToLGamma, slha_, files_] :=
Module[{
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      LToLGammaFields, properStates, wrongFields,
      prototypes = "", definitions = ""
   },

   LToLGammaFields = DeleteDuplicates[observables /.
      {
         Rule[i_[_], {o_[_], v_}] :> Rule[i, {o, v}],
         obs -> Sequence
      }
   ];

   properStates = Cases[
      LToLGammaFields,
      Rule[a_?IsLepton, {b_?IsLepton, c_ /; c === GetPhoton[]}] -> (a -> {b, c})
   ];
   wrongFields = Complement[LToLGammaFields, properStates];
   If[wrongFields =!= {},
      Utils`FSFancyWarning[
         "BrLToLGamma works only for leptons and a photon.",
         " Removing requested process(es): ",
         Riffle[ToString/@wrongFields, ", "]
      ];
      LToLGammaFields = properStates;
   ];

   If[LToLGammaFields =!= {},
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];

      {prototypes, definitions} =
         StringJoin@@@(Riffle[#, "\n\n"] &/@ Transpose[createBrLToLGamma/@ LToLGammaFields]);
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@calculate_prototypes@" -> prototypes,
         "@calculate_definitions@" -> definitions,
         "@include_guard@" -> SymbolName@obs,
         "@namespace@" -> Observables`GetObservableNamespace@obs,
         "@filename@" -> Observables`GetObservableFileName@obs,
         "@get_MSUSY@" -> TextFormatting`IndentText@TextFormatting`WrapLines@AMM`AMMGetMLCP[],
         Sequence@@GeneralReplacementRules[]
      }
   ];

   {
      "FFV fields" -> DeleteDuplicates@LToLGammaFields
   }
];

createBrLToLGamma[inFermion_ -> {outFermion_, spectator_}] :=
    Module[{prototype, definition,
            numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
            numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
            numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator],
            (* don't remove potential SM contributions to l -> l gamma *)
            discardSMcontributions = CConversion`CreateCBoolValue[False]},

        prototype =
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <>
            CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "(\n" <>
            If[TreeMasses`GetDimension[inFermion] =!= 1,
               "int generationIndex1, ",
               ""
            ] <>
            If[TreeMasses`GetDimension[outFermion] =!= 1,
               "int generationIndex2, \n",
               ""
            ] <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates_running& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input);";

      definition =
            (* calculate observable using form factors *)
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> " (\n" <>
               IndentText[
                  If[TreeMasses`GetDimension[inFermion] =!= 1, "int generationIndex1, ", ""] <>
                  If[TreeMasses`GetDimension[outFermion] =!= 1, "int generationIndex2, ", ""] <> "\n" <>
                  "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates_running& model_, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input) {\n\n"
               ] <>
               (* write routine for mu to e gamma *)
               IndentText[
                  FlexibleSUSY`FSModelName<>"_mass_eigenstates_running model(model_);\n"<>
                  "try {\n"<>
                  IndentText[
                     "run_to_MSUSY(model);\n"<>
                     "model.get_physical().clear();\n"<>
                     "model.solve_ewsb();\n"
                  ]<>
                  "} catch (const Error& e) {\n"<>
                  IndentText[
                     "ERROR(\""<>FlexibleSUSY`FSModelName<>
                        "_l_to_lgamma:\" << e.what_detailed());\n"<>
                     "return std::numeric_limits<double>::quiet_NaN();\n"
                  ]<>
                  "}\n"<>
                  "context_base context {model};\n" <>
                  "std::array<int, " <> ToString @ numberOfIndices1 <>
                     "> indices1 = {" <>
                     (* TODO(wkotlarski): Specify indices correctly *)
                       If[TreeMasses`GetDimension[inFermion] =!= 1,
                          " generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringRiffle[Table[" 0", {numberOfIndices1}], ","] <> " ",
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[outFermion] =!= 1,
                          " generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringRiffle[Table[" 0", {numberOfIndices2}], ","] <> " ",
                             ""]
                         ] <> "};\n\n" <>
                    "const auto form_factors = calculate_" <> CXXNameOfField[inFermion] <> "_"
                   <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors "<>
                   "(" <> If[TreeMasses`GetDimension[inFermion] =!= 1,
                            "generationIndex1, ",
                            ""] <>
                         If[TreeMasses`GetDimension[outFermion] =!= 1,
                            "generationIndex2, ",
                            ""] <>
                  "model, " <> discardSMcontributions <> ");\n" <>
                  (* Dominik suggest that the phase space prefactor should use pole masses  so we get them from the input file *)
                  "double leptonInMassOS;\n" <>
                  "switch (generationIndex1) {\n" <>
                  IndentText[
                     "case 0: leptonInMassOS = qedqcd.displayMass(softsusy::mElectron); break;\n" <>
                     "case 1: leptonInMassOS = qedqcd.displayMass(softsusy::mMuon);     break;\n" <>
                     "case 2: leptonInMassOS = qedqcd.displayMass(softsusy::mTau);      break;\n" <>
                     "default: throw std::invalid_argument(\"Unrecognized lepton\");\n"
                  ] <>
                  "}\n" <>
                  "\n" <>
                  "// eq. 51 of arXiv:hep-ph/9510309 (note that we include 'e' in the definition of form_factor)\n" <>
                  "const double partial_width = pow(leptonInMassOS,5)/(16.0*Pi) * (std::norm(form_factors[2]) + std::norm(form_factors[3]));\n" <>

                  "const double total_width = lepton_total_decay_width<" <>
                     CXXNameOfField[inFermion] <> ", " <> CXXNameOfField[outFermion] <>
                     ">(indices1, indices2, model, qedqcd);\n" <>
                  "\nreturn partial_width/total_width;\n"
               ] <> "}";

        {prototype, definition}
    ];

End[];
