(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["Observables`", {"FlexibleSUSY`", "SARAH`", "BetaFunction`", "Parameters`", "TreeMasses`", "Utils`", "CConversion`", "TextFormatting`", "GM2Calc`"}];

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { AMM, AMMUncertainty, aMuonGM2Calc, aMuonGM2CalcUncertainty,
                  EDM, bsgamma};
End[];

DefineObservable::usage="Defines all C++ names for an observable";
GetObservablesHeaders::usage="Return the names of C++ headers with #include";
GetObservableNamespace::usage="Return the C++ context name of an observable";
GetObservableFileName::usage="Returns the C++ file name of a given observable";
GetObservablePrototype::usage="Returns the C++ prototype as string";
GetRequestedObservables::usage="";
CountNumberOfObservables::usage="";
CreateObservablesDefinitions::usage="";
CreateObservablesInitialization::usage="";
CreateSetAndDisplayObservablesFunctions::usage="";
CreateClearObservablesFunction::usage="";
CalculateObservables::usage="";
GetObservableName::usage="returns name of observable in Observables struct";
GetObservableType::usage="returns type of observable";
GetObservableDescription::usage="returns description of observable.";
IsObservable::usage = "Returns true if given symbol is an observable.";

Begin["`Private`"];

IsObservable[sym_] :=
    MemberQ[FlexibleSUSYObservable`FSObservables, sym] || \
    (Or @@ (MatchQ[sym, #[__]]& /@ FlexibleSUSYObservable`FSObservables));

GetRequestedObservables[blocks_] :=
    Module[{observables, dim, test},
           observables = DeleteDuplicates[Cases[blocks, a_?IsObservable :> a, {0, Infinity}]];
           observables
          ];

GetObservableName[FlexibleSUSYObservable`AMM[p_[idx_]]] := GetObservableName[FlexibleSUSYObservable`AMM[p]] <> "_" <> ToString[idx];
GetObservableName[FlexibleSUSYObservable`AMM[p_]] := "amm_" <> CConversion`ToValidCSymbolString[p];
GetObservableName[FlexibleSUSYObservable`AMMUncertainty[p_]] := "amm_uncertainty_" <> CConversion`ToValidCSymbolString[p];
GetObservableName[FlexibleSUSYObservable`AMMUncertainty[p_[idx_]]] := GetObservableName[FlexibleSUSYObservable`AMMUncertainty[p]] <> "_" <> ToString[idx];
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon_gm2calc";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "a_muon_gm2calc_uncertainty";
GetObservableName[FlexibleSUSYObservable`EDM[p_[idx_]]] := GetObservableName[FlexibleSUSYObservable`EDM[p]] <> "_" <> ToString[idx];
GetObservableName[FlexibleSUSYObservable`EDM[p_]]       := "edm_" <> CConversion`ToValidCSymbolString[p];

GetObservableName[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := "b_to_s_gamma";

GetObservableDescription[FlexibleSUSYObservable`AMM[p_[idx_]]] := "Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx+1] <> ") (calculated with FlexibleSUSY)";
GetObservableDescription[FlexibleSUSYObservable`AMM[p_]] := "Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> " (calculated with FlexibleSUSY)";
GetObservableDescription[FlexibleSUSYObservable`AMMUncertainty[p_]] := "uncertainty of Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> " (calculated with FlexibleSUSY)";
GetObservableDescription[FlexibleSUSYObservable`AMMUncertainty[p_[idx_]]] := "uncertainty of Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx+1] <> ") (calculated with FlexibleSUSY)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon = (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "uncertainty of (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[FlexibleSUSYObservable`EDM[p_[idx_]]] := "electric dipole moment of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx+1] <> ") [1/GeV]";
GetObservableDescription[FlexibleSUSYObservable`EDM[p_]]       := "electric dipole moment of " <> CConversion`ToValidCSymbolString[p] <> " [1/GeV]";

GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := "calculates the Wilson coefficients C7 and C8 for b -> s gamma";

GetObservableType[FlexibleSUSYObservable`AMM[_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`AMMUncertainty[_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`EDM[p_]] := CConversion`ScalarType[CConversion`realScalarCType];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := CConversion`ScalarType[CConversion`realScalarCType];

CountNumberOfObservables[observables_List] :=
    Module[{i, number = 0},
           For[i = 1, i <= Length[observables], i++,
               If[IsObservable[observables[[i]]],
                  number += BetaFunction`CountNumberOfParameters[GetObservableType[observables[[i]]]];,
                  Utils`FSFancyWarning["Ignoring invalid observable ", observables[[i]]];
               ];
           ];
           number
    ];

CreateObservablesDefinitions[observables_List] :=
    Module[{i, type, name, description, definitions = ""},
           For[i = 1, i <= Length[observables], i++,
               If[IsObservable[observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  description = GetObservableDescription[observables[[i]]];
                  type = CConversion`CreateCType[GetObservableType[observables[[i]]]];
                  definitions = definitions <> type <> " " <> name <> "; ///< " <> description <> "\n";,
                  Utils`FSFancyWarning["Ignoring invalid observable ", observables[[i]]];
               ];
           ];
           definitions
    ];

CreateObservablesInitialization[observables_List] :=
    Module[{i, name, type, init = ""},
           For[i = 1, i <= Length[observables], i++,
               If[IsObservable[observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  type = GetObservableType[observables[[i]]];
                  If[init == "",
                     init = ": " <> CConversion`CreateDefaultConstructor[name, type] <> "\n";,
                     init = init <> ", " <> CConversion`CreateDefaultConstructor[name, type] <> "\n";
                    ];,
                  Utils`FSFancyWarning["Ignoring invalid observable ", observables[[i]]];
               ];
           ];
           init
   ];

CreateSetAndDisplayObservablesFunctions[observables_List] :=
    Module[{numObservables, i, name, type, paramCount = 0, nAssignments, assignment,
            display = "", displayNames = "", set = ""},
           numObservables = CountNumberOfObservables[observables];
           If[numObservables != 0,
              display = "Eigen::ArrayXd vec(" <> FlexibleSUSY`FSModelName
                        <> "_observables::NUMBER_OF_OBSERVABLES);\n\n";
              displayNames = "std::vector<std::string> names("
                             <> FlexibleSUSY`FSModelName
                             <> "_observables::NUMBER_OF_OBSERVABLES);\n\n";
              set = "assert(vec.rows() == " <> FlexibleSUSY`FSModelName
                    <> "_observables::NUMBER_OF_OBSERVABLES);\n\n";
              For[i = 1, i <= Length[observables], i++,
                  If[IsObservable[observables[[i]]],
                     name = GetObservableName[observables[[i]]];
                     type = GetObservableType[observables[[i]]];
                     {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type, "vec"];
                     set = set <> assignment;
                     {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type, "vec"];
                     display = display <> assignment;
                     {assignment, nAssignments} = Parameters`CreateStdVectorNamesAssignment[name, paramCount, type];
                     displayNames = displayNames <> assignment;
                     paramCount += nAssignments;,
                     Utils`FSFancyWarning["Ignoring invalid observable ", observables[[i]]];
                  ];
               ];,
               display = "Eigen::ArrayXd vec(1);\n\nvec(0) = 0.;\n";
               set = "";
               displayNames = "std::vector<std::string> names(1);\n\n"
                              <> "names[0] = \"no observables defined\";\n";
             ];
           {display, displayNames, set}
          ];

CreateClearObservablesFunction[observables_List] :=
    Module[{i, name, type, result = ""},
           For[i = 1, i <= Length[observables], i++,
               If[IsObservable[observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  type = GetObservableType[observables[[i]]];
                  result = result <> CConversion`SetToDefault[name, type];,
                  Utils`FSFancyWarning["Ignoring invalid observable ", observables[[i]]];
               ];
           ];
           result
          ];

CalculateObservable[FlexibleSUSYObservable`AMM[p_], structName_String] :=
    structName <> ".AMM0(" <> CConversion`ToValidCSymbolString[p] <> ") = " <>
      FlexibleSUSY`FSModelName <> "_amm::calculate_amm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, settings);";

CalculateObservable[FlexibleSUSYObservable`AMM[p_[idx_]], structName_String] :=
    structName <> ".AMM1(" <> CConversion`ToValidCSymbolString[p] <> ", " <> ToString[idx] <> ") = " <>
      FlexibleSUSY`FSModelName <> "_amm::calculate_amm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, settings," <> ToString[idx] <> ");";

CalculateObservable[FlexibleSUSYObservable`AMMUncertainty[p_], structName_String] :=
    structName <> ".AMMUNCERTAINTY0(" <> CConversion`ToValidCSymbolString[p] <> ") = " <> FlexibleSUSY`FSModelName <> "_amm::calculate_amm_uncertainty<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, settings);";

CalculateObservable[FlexibleSUSYObservable`AMMUncertainty[p_[idx_]], structName_String] :=
    structName <> ".AMMUNCERTAINTY1(" <> CConversion`ToValidCSymbolString[p] <> ", " <> ToString[idx] <> ") = " <> FlexibleSUSY`FSModelName <> "_amm::calculate_amm_uncertainty<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, settings, " <> ToString[idx] <> ");";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    "#ifdef ENABLE_GM2CALC\n" <>
    structName <> ".AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);\n" <>
    "#endif";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty, structName_String] :=
    "#ifdef ENABLE_GM2CALC\n" <>
    structName <> ".AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);\n" <>
    "#endif";

CalculateObservable[FlexibleSUSYObservable`EDM[p_], structName_String] :=
    Module[{pStr = CConversion`ToValidCSymbolString[p]},
           structName <> ".EDM0(" <> pStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_edm::calculate_edm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd);"
          ];

CalculateObservable[FlexibleSUSYObservable`EDM[p_[idx_]], structName_String] :=
    Module[{pStr = CConversion`ToValidCSymbolString[p],
            idxStr = ToString[idx]},
           structName <> ".EDM1(" <> pStr <> ", " <> idxStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_edm::calculate_edm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, " <> idxStr <> ");"
          ];

(* TODO: move Wilson Coefficients to a different block *)
CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`bsgamma, structName_String] :=
    structName <> ".BSGAMMA = Re(" <> FlexibleSUSY`FSModelName <> "_b_to_s_gamma::calculate_b_to_s_gamma(MODEL, qedqcd)[0]);";

FillInterfaceData[{}] := "";

FillInterfaceData[obs_List] :=
    Module[{filled = ""},
           If[MemberQ[obs,FlexibleSUSYObservable`aMuonGM2Calc] ||
              MemberQ[obs,FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
              filled = filled <> GM2Calc`FillGM2CalcInterfaceData["gm2calc_data"];
             ];
           filled
          ];

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?IsObservable :> a, {0, Infinity}];
           FillInterfaceData[observables] <> "\n" <>
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

Options@DefineObservable = {
   InsertionFunction -> CConversion`ToValidCSymbolString,
   GetObservableType -> Unset,
   GetObservableName -> Unset,
   GetObservablePrototype -> Unset,
   GetObservableDescription -> Unset,
   CalculateObservable -> Unset,
   GetObservableFileName -> Unset,
   GetObservableNamespace -> Unset
};

DefineObservable[obs_@pattern___, OptionsPattern[]] :=
Module[{stringPattern, patternNames, uniqueNames, lhsRepl, rhsRepl, warn,
      obsStr = SymbolName@obs, extraCalc, nameAux, specialWords},
   warn := Utils`FSFancyWarning[#," for ", ToString@obs, " might not be specified."]&;
   extraCalc := StringReplace[#, "$(" ~~ Shortest[x__] ~~ ")" :> ToString@ToExpression[x]]&;

   stringPattern = ToString@FullForm@{pattern};
   patternNames = DeleteDuplicates@StringCases[stringPattern, "Pattern[" ~~ Shortest@x__ ~~ "," :> x];
   uniqueNames = Unique[#<>"$"]&/@patternNames;

   lhsRepl = MapThread[Rule, {patternNames, ToString/@uniqueNames}];
   rhsRepl = MapThread[
      RuleDelayed[#1, OptionValue[InsertionFunction]@#2]&,
      {patternNames, uniqueNames}
   ];

   nameAux = StringReplace[obsStr, nums___?DigitQ ~~ x_?UpperCaseQ :> "_" <> nums <> ToLowerCase[x]];
   Switch[OptionValue[GetObservableFileName],
      Unset,
         GetObservableFileName[obs | obsStr]  = StringDrop[nameAux, 1];,
      _String,
         GetObservableFileName[obs | obsStr]  = OptionValue@GetObservableFileName;
   ];
   Switch[OptionValue[GetObservableNamespace],
      Unset,
         GetObservableNamespace[obs | obsStr] = OptionValue[InsertionFunction][FlexibleSUSY`FSModelName <> nameAux];,
      _String,
         GetObservableNamespace[obs | obsStr] = OptionValue[InsertionFunction][FlexibleSUSY`FSModelName <> "_" <> OptionValue@GetObservableNamespace];
   ];

   specialWords = {
      "auto model" -> OptionValue[InsertionFunction]["const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model"],
      "auto qedqcd" -> OptionValue[InsertionFunction]["const softsusy::QedQcd& qedqcd"],
      "context" -> GetObservableNamespace[obs]
   };
   rhsRepl = Join[rhsRepl, specialWords];

   With[{args = Sequence@@ToExpression@StringReplace[stringPattern, lhsRepl],
         repl = rhsRepl,
         name = If[# === Unset, ToLowerCase@SymbolName@obs <> "_" <> StringRiffle[patternNames, "_"], #] &@ OptionValue@GetObservableName,
         description = OptionValue@GetObservableDescription,
         type = OptionValue@GetObservableType,
         calculate = OptionValue@CalculateObservable,
         prototype = OptionValue@GetObservablePrototype
      },
      AppendTo[FlexibleSUSYObservable`FSObservables, obs];

      GetObservableName@obs@args := extraCalc@StringReplace[name, repl];

      Switch[type,
         Unset,
            warn@"GetObservableType",
         {_Integer},
            GetObservableType@obs@args :=
               CConversion`ArrayType[CConversion`complexScalarCType, First@type],
         _,
            GetObservableType@obs@args := type
      ];

      Switch[description,
         Unset,
            GetObservableDescription@obs@args := StringReplace[extraCalc@StringReplace[name, repl], "_" -> " "];,
         _String,
            GetObservableDescription@obs@args := extraCalc@StringReplace[description, repl];
      ];

      If[prototype === Unset,
         warn@"GetObservablePrototype";,
         Module[{nameStr, argStr},
            nameStr = StringReplace[prototype, Longest["(" ~~ s___ ~~ ")"]  :> (argStr = s; "")];
            argStr = StringReplace[argStr, specialWords];
            GetObservablePrototype@obs@args := extraCalc@StringReplace[nameStr, repl] <> "(" <> argStr <> ")";
         ];
      ];

      Switch[calculate,
         Unset,
            Module[{nameStr, argStr, removeSub, newCalculate, niceStr = {}},
               nameStr = StringReplace[prototype, Longest["(" ~~ s___ ~~ ")"]  :> (argStr = s; "")];
               removeSub[pair_] := StringReplace[#, Longest[StringTake[pair, 1] ~~ ___ ~~ StringTake[pair, -1]] :> " "]&;
               (argStr = FixedPoint[removeSub@#, argStr]) &/@ {"()" , "[]", "{}", "<>"};
               argStr = StringReverse@StringReplace[argStr <> ",", {"&" -> " ", Whitespace ~~ "," :> ","}];
               StringReplace[argStr, Shortest["," ~~ x__ ~~ Whitespace] :> (AppendTo[niceStr, StringReverse[x]]; "")];
               newCalculate = nameStr <> "(" <> StringRiffle[Reverse@niceStr, ", "] <> ")";

               CalculateObservable[obs@args, structName:_String] :=
                  structName <> "." <> StringReplace[name, repl] <>
                  extraCalc@StringReplace[" = context::"<>newCalculate<>";", repl];
            ];,
         _String,
            CalculateObservable[obs@args, structName:_String] :=
               structName <> "." <> StringReplace[name, repl] <>
               extraCalc@StringReplace[" = context::"<>calculate<>";", repl];
      ];
   ];
];
Utils`MakeUnknownInputDefinition@DefineObservable;

GetObservablesHeaders[] =
Module[{files, obs},
   files = Utils`DynamicInclude@$observablesWildcard@"Observables.m";
   obs = StringSplit[files, $PathnameSeparator][[All, -2]];
   StringRiffle[
      "#include \"observables/"<>FlexibleSUSY`FSModelName<>"_"<>Observables`GetObservableFileName@#<>".hpp\""&/@obs,
      "\n"
   ]
];

End[];
EndPackage[];
