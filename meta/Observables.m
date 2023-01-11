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

BeginPackage["Observables`", {"FlexibleSUSY`", "SARAH`", "BetaFunction`", "Parameters`", "TreeMasses`", "Utils`", "CConversion`", "TextFormatting`"}];

(* GM2Calc interface parameters (gauge basis) *)
{ yukawaType, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6,
  lambda7, tanBeta, m122, zetau, zetad, zetal, deltau, deltad, deltal,
  piu, pid, pil };

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { AMM, AMMUncertainty, aMuonGM2Calc, aMuonGM2CalcUncertainty,
                  EDM, BrLToLGamma, bsgamma };

If[FlexibleSUSY`FSFeynArtsAvailable && FlexibleSUSY`FSFormCalcAvailable,
   AppendTo[FSObservables, FToFConversionInNucleus]
];

End[];

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
           test = Complement[
              Cases[observables, _FlexibleSUSYObservable`BrLToLGamma],
              Cases[observables, FlexibleSUSYObservable`BrLToLGamma[fin_?IsLepton -> {fout_?IsLepton, vout_ /; vout === GetPhoton[]}]]
                 ];
           If[test =!= {},
              Utils`FSFancyWarning[
                 "BrLToLGamma function works only for leptons and a photon.",
                 " Removing requested process(es): ", test
              ];
              observables = Complement[observables, test];
           ];
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
GetObservableName[FlexibleSUSYObservable`BrLToLGamma[pIn_[idxIn_] -> {pOut_[idxOut_], spectator_}]] := CConversion`ToValidCSymbolString[pIn] <> ToString[idxIn] <> "_to_" <> CConversion`ToValidCSymbolString[pOut] <> ToString[idxOut] <> "_" <> CConversion`ToValidCSymbolString[spectator];
GetObservableName[FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, spectator_}]] := CConversion`ToValidCSymbolString[pIn] <> "_to_" <> CConversion`ToValidCSymbolString[pOut] <> "_" <> CConversion`ToValidCSymbolString[spectator];
GetObservableName[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_[idxIn_] -> pOut_[idxOut_], nucleus_]] := CConversion`ToValidCSymbolString[pIn] <> "_to_" <> CConversion`ToValidCSymbolString[pOut] <> "_in_" <> ToString@nucleus;
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := "b_to_s_gamma";

GetObservableDescription[FlexibleSUSYObservable`AMM[p_[idx_]]] := "Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx+1] <> ") (calculated with FlexibleSUSY)";
GetObservableDescription[FlexibleSUSYObservable`AMM[p_]] := "Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> " (calculated with FlexibleSUSY)";
GetObservableDescription[FlexibleSUSYObservable`AMMUncertainty[p_]] := "uncertainty of Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> " (calculated with FlexibleSUSY)";
GetObservableDescription[FlexibleSUSYObservable`AMMUncertainty[p_[idx_]]] := "uncertainty of Delta(g-2)/2 of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx+1] <> ") (calculated with FlexibleSUSY)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon = (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "uncertainty of (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[FlexibleSUSYObservable`EDM[p_[idx_]]] := "electric dipole moment of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx+1] <> ") [1/GeV]";
GetObservableDescription[FlexibleSUSYObservable`EDM[p_]]       := "electric dipole moment of " <> CConversion`ToValidCSymbolString[p] <> " [1/GeV]";
GetObservableDescription[FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, _}]] :=
   "BR(" <> CConversion`ToValidCSymbolString[pIn] <> " -> " <>
   CConversion`ToValidCSymbolString[pOut] <> " " <>
   CConversion`ToValidCSymbolString[V] <> ")"  ;
GetObservableDescription[FlexibleSUSYObservable`BrLToLGamma[pIn_[idxIn_] -> {pOut_[idxOut_], V_}]] :=
   "BR(" <> CConversion`ToValidCSymbolString[pIn] <> ToString[idxIn] <> " -> " <>
      CConversion`ToValidCSymbolString[pOut] <> ToString[idxOut] <> " " <>
       CConversion`ToValidCSymbolString[V] <> ")"  ;
GetObservableDescription[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_ -> pOut_, nuc_]] :=
   "CR(" <> CConversion`ToValidCSymbolString[pIn] <> " -> " <>
      CConversion`ToValidCSymbolString[pOut] <> ", " <>
      ToString[nuc] <> ")/capture rate";
GetObservableDescription[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_[idxIn_] -> pOut_[idxOut_], nuc_]] :=
   "CR(" <> CConversion`ToValidCSymbolString[pIn] <> ToString[idxIn] <> " -> " <>
   CConversion`ToValidCSymbolString[pOut] <> ToString[idxOut] <> ", " <>
      ToString[nuc] <> ")/capture rate";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := "calculates the Wilson coefficients C7 and C8 for b -> s gamma";

GetObservableType[FlexibleSUSYObservable`AMM[_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`AMMUncertainty[_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`EDM[p_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, _}]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_[idxIn_] -> pOut_[idxOut_], _]] := CConversion`ScalarType[CConversion`realScalarCType];
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
      FlexibleSUSY`FSModelName <> "_amm::calculate_amm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd);";

CalculateObservable[FlexibleSUSYObservable`AMM[p_[idx_]], structName_String] :=
    structName <> ".AMM1(" <> CConversion`ToValidCSymbolString[p] <> ", " <> ToString[idx] <> ") = " <>
      FlexibleSUSY`FSModelName <> "_amm::calculate_amm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, " <> ToString[idx] <> ");";

CalculateObservable[FlexibleSUSYObservable`AMMUncertainty[p_], structName_String] :=
    structName <> ".AMMUNCERTAINTY0(" <> CConversion`ToValidCSymbolString[p] <> ") = " <> FlexibleSUSY`FSModelName <> "_amm::calculate_amm_uncertainty<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd);";

CalculateObservable[FlexibleSUSYObservable`AMMUncertainty[p_[idx_]], structName_String] :=
    structName <> ".AMMUNCERTAINTY1(" <> CConversion`ToValidCSymbolString[p] <> ", " <> ToString[idx] <> ") = " <> FlexibleSUSY`FSModelName <> "_amm::calculate_amm_uncertainty<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, " <> ToString[idx] <> ");";

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
           FlexibleSUSY`FSModelName <> "_edm::calculate_edm_" <> pStr <> "(MODEL);"
          ];

CalculateObservable[FlexibleSUSYObservable`EDM[p_[idx_]], structName_String] :=
    Module[{pStr = CConversion`ToValidCSymbolString[p],
            idxStr = ToString[idx]},
           structName <> ".EDM1(" <> pStr <> ", " <> idxStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_edm::calculate_edm<" <> CXXDiagrams`CXXNameOfField[p, prefixNamespace -> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">(MODEL, qedqcd, " <> idxStr <> ");"
          ];

CalculateObservable[FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, spectator_}], structName_String] :=
    Module[{pInStr = CConversion`ToValidCSymbolString[pIn], pOutStr = CConversion`ToValidCSymbolString[pOut],
    spec = CConversion`ToValidCSymbolString[spectator]},
           structName <> ".LToLGamma0(" <> pInStr <> ", " <> pOutStr <> ", " <> spec <> ") = " <>
           FlexibleSUSY`FSModelName <> "_l_to_lgamma::calculate_" <> pInStr <> "_to_" <> pOutStr <> "_" <> spec <> "(MODEL, qedqcd, physical_input);"
          ];

CalculateObservable[FlexibleSUSYObservable`BrLToLGamma[pIn_[idxIn_] -> {pOut_[idxOut_], spectator_}], structName_String] :=
    Module[{pInStr = CConversion`ToValidCSymbolString[pIn],
            pOutStr = CConversion`ToValidCSymbolString[pOut],
            idxInStr = ToString[idxIn],
            idxOutStr = ToString[idxOut],
            specStr = ToString[spectator]
    },
           structName <> ".LToLGamma1(" <> pInStr <> ", " <> idxInStr <> ", " <> pOutStr <> ", " <> idxOutStr <> ", " <> specStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_l_to_lgamma::calculate_" <> pInStr <> "_to_" <> pOutStr <> "_" <> specStr <> "(" <> idxInStr <> ", " <> idxOutStr <> ", MODEL, qedqcd, physical_input);"
          ];

CalculateObservable[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_ -> pOut_, nucleai_], structName_String] :=
    Module[{pInStr = CConversion`ToValidCSymbolString[pIn], pOutStr = CConversion`ToValidCSymbolString[pOut],
    nuc = CConversion`ToValidCSymbolString[nucleai]},
           structName <> ".FToFConversion0(" <> pInStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_f_to_f_conversion::calculate_" <> pInStr <> "_to_" <> pOutStr <> "_in_nucleus(MODEL);"
          ];

CalculateObservable[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_[idxIn_] -> pOut_[idxOut_], nucleus_], structName_String] :=
    Module[{pInStr = CConversion`ToValidCSymbolString[pIn],
            pOutStr = CConversion`ToValidCSymbolString[pOut],
            idxInStr = ToString[idxIn],
            idxOutStr = ToString[idxOut],
            nucleiStr = ToString[nucleus]
    },
           structName <> ".FToFConversion1(" <> pInStr <> ", " <> idxInStr <> ", " <> pOutStr <> ", " <> idxOutStr <> ", " <> nucleiStr <> ", " <> "qedqcd) = " <>
           FlexibleSUSY`FSModelName <> "_f_to_f_conversion::calculate_" <> pInStr <> "_to_" <> pOutStr <> "_in_nucleus(" <> idxInStr <> ", " <> idxOutStr <> ", "<> FlexibleSUSY`FSModelName <> "_f_to_f_conversion::Nucleus::" <> nucleiStr <> ", MODEL, qedqcd);"
          ];

(* TODO: move Wilson Coefficients to a different block *)
CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`bsgamma, structName_String] :=
    structName <> ".BSGAMMA = Re(" <> FlexibleSUSY`FSModelName <> "_b_to_s_gamma::calculate_b_to_s_gamma(MODEL, qedqcd)[0]);";

(* fill struct with model parameters to be passed to GM2Calc *)
FillGM2CalcInterfaceData[struct_String] :=
    Which[
        IsGM2CalcCompatibleMSSM[], FillGM2CalcMSSMNoFVInterfaceData[struct],
        IsGM2CalcCompatibleTHDM[], FillGM2CalcTHDMInterfaceData[struct, FlexibleSUSY`FSGM2CalcInput],
        True, Print["Error: This model is neither a MSSM-like nor a 2HDM-like model compatible with GM2Calc."]; Quit[1]
    ];

(* returns true, if model is an MSSM-like model compatible with GM2Calc *)
IsGM2CalcCompatibleMSSM[] :=
    Module[{w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
            mu, m1, m2, m3, mq2, mu2, md2, ml2, me2, tu, td, te, yu, yd, ye},
           w             = Parameters`GetParticleFromDescription["W-Boson"];
           pseudoscalar  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           smuon         = Parameters`GetParticleFromDescription["Smuon"];
           muonsneutrino = Parameters`GetParticleFromDescription["Muon Sneutrino"];
           chargino      = Parameters`GetParticleFromDescription["Charginos"];
           neutralino    = Parameters`GetParticleFromDescription["Neutralinos"];
           mu            = Parameters`GetParameterFromDescription["Mu-parameter"];
           m1            = Parameters`GetParameterFromDescription["Bino Mass parameter"];
           m2            = Parameters`GetParameterFromDescription["Wino Mass parameter"];
           m3            = Parameters`GetParameterFromDescription["Gluino Mass parameter"];
           mq2           = Parameters`GetParameterFromDescription["Softbreaking left Squark Mass"];
           mu2           = Parameters`GetParameterFromDescription["Softbreaking right Up-Squark Mass"];
           md2           = Parameters`GetParameterFromDescription["Softbreaking right Down-Squark Mass"];
           ml2           = Parameters`GetParameterFromDescription["Softbreaking left Slepton Mass"];
           me2           = Parameters`GetParameterFromDescription["Softbreaking right Slepton Mass"];
           tu            = Parameters`GetParameterFromDescription["Trilinear-Up-Coupling"];
           td            = Parameters`GetParameterFromDescription["Trilinear-Down-Coupling"];
           te            = Parameters`GetParameterFromDescription["Trilinear-Lepton-Coupling"];
           yu            = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd            = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye            = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           Not[MemberQ[{w, pseudoscalar, smuon, muonsneutrino,
                        chargino, neutralino, mu, m1, m2, m3, mq2, mu2,
                        md2, ml2, me2, tu, td, te, yu, yd, ye}, Null]]
    ];

(* fill struct with MSSM parameters to be passed to GM2Calc *)
FillGM2CalcMSSMNoFVInterfaceData[struct_String] :=
    Module[{mwStr, w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
            mu, m1, m2, m3, mq2, mu2, md2, ml2, me2, tu, td, te, yu, yd, ye},
           w             = Parameters`GetParticleFromDescription["W-Boson"];
           pseudoscalar  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           smuon         = Parameters`GetParticleFromDescription["Smuon"];
           muonsneutrino = Parameters`GetParticleFromDescription["Muon Sneutrino"];
           chargino      = Parameters`GetParticleFromDescription["Charginos"];
           neutralino    = Parameters`GetParticleFromDescription["Neutralinos"];
           mu            = Parameters`GetParameterFromDescription["Mu-parameter"];
           m1            = Parameters`GetParameterFromDescription["Bino Mass parameter"];
           m2            = Parameters`GetParameterFromDescription["Wino Mass parameter"];
           m3            = Parameters`GetParameterFromDescription["Gluino Mass parameter"];
           mq2           = Parameters`GetParameterFromDescription["Softbreaking left Squark Mass"];
           mu2           = Parameters`GetParameterFromDescription["Softbreaking right Up-Squark Mass"];
           md2           = Parameters`GetParameterFromDescription["Softbreaking right Down-Squark Mass"];
           ml2           = Parameters`GetParameterFromDescription["Softbreaking left Slepton Mass"];
           me2           = Parameters`GetParameterFromDescription["Softbreaking right Slepton Mass"];
           tu            = Parameters`GetParameterFromDescription["Trilinear-Up-Coupling"];
           td            = Parameters`GetParameterFromDescription["Trilinear-Down-Coupling"];
           te            = Parameters`GetParameterFromDescription["Trilinear-Lepton-Coupling"];
           yu            = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd            = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye            = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           mwStr         = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           (* compose string *)
           "#ifdef ENABLE_GM2CALC\n" <>
           "GM2Calc_MSSMNoFV_data " <> struct <> ";\n" <>
           struct <> ".scale = MODEL.get_scale();\n" <>
           struct <> ".alpha_em_MZ = ALPHA_EM_MZ;\n" <>
           struct <> ".alpha_em_0 = ALPHA_EM_0;\n" <>
           struct <> ".alpha_s_MZ = ALPHA_S_MZ;\n" <>
           struct <> ".MZ    = MZPole;\n" <>
           "if (!is_zero(" <> mwStr <> ")) {\n" <>
              TextFormatting`IndentText[struct <> ".MW = " <> mwStr <> ";"] <> "\n" <>
           "} else if (!is_zero(MWPole)) {\n" <>
              TextFormatting`IndentText[struct <> ".MW = MWPole;"] <> "\n}\n" <>
           struct <> ".mb_mb = MBMB;\n" <>
           struct <> ".MT    = MTPole;\n" <>
           struct <> ".MTau  = MTauPole;\n" <>
           struct <> ".MM    = MMPole;\n" <>
           struct <> ".MA0   = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[pseudoscalar][1]] <> ";\n" <>
           struct <> ".MSvm  = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[muonsneutrino]] <> ";\n" <>
           struct <> ".TB    = MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> "() / " <>
                              "MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> "();\n" <>
           struct <> ".Mu    = MODEL.get_" <> CConversion`RValueToCFormString[mu] <> "();\n" <>
           struct <> ".M1    = MODEL.get_" <> CConversion`RValueToCFormString[m1] <> "();\n" <>
           struct <> ".M2    = MODEL.get_" <> CConversion`RValueToCFormString[m2] <> "();\n" <>
           struct <> ".M3    = MODEL.get_" <> CConversion`RValueToCFormString[m3] <> "();\n" <>
           struct <> ".MSm   = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[smuon]] <> ";\n" <>
           struct <> ".MCha  = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[chargino]] <> ";\n" <>
           struct <> ".MChi  = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[neutralino]] <> ";\n" <>
           struct <> ".mq2   = MODEL.get_" <> CConversion`RValueToCFormString[mq2] <> "();\n" <>
           struct <> ".mu2   = MODEL.get_" <> CConversion`RValueToCFormString[mu2] <> "();\n" <>
           struct <> ".md2   = MODEL.get_" <> CConversion`RValueToCFormString[md2] <> "();\n" <>
           struct <> ".ml2   = MODEL.get_" <> CConversion`RValueToCFormString[ml2] <> "();\n" <>
           struct <> ".me2   = MODEL.get_" <> CConversion`RValueToCFormString[me2] <> "();\n" <>
           struct <> ".Au    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[tu] <> "(), MODEL.get_" <> CConversion`RValueToCFormString[yu] <> "());\n" <>
           struct <> ".Ad    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[td] <> "(), MODEL.get_" <> CConversion`RValueToCFormString[yd] <> "());\n" <>
           struct <> ".Ae    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[te] <> "(), MODEL.get_" <> CConversion`RValueToCFormString[ye] <> "());\n" <>
           "#endif\n\n"
          ];

(* returns true, if model is an THDM-like model compatible with GM2Calc *)
IsGM2CalcCompatibleTHDM[] :=
    Module[{w, h, a, hp, yu, yd, ye},
           w  = Parameters`GetParticleFromDescription["W-Boson"];
           h  = Parameters`GetParticleFromDescription["Higgs"];
           a  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           hp = Parameters`GetParticleFromDescription["Charged Higgs"];
           yu = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           Not[MemberQ[{w, h, a, hp, yu, yd, ye}, Null]]
    ];

(* fill struct with THDM parameters to be passed to GM2Calc *)
FillGM2CalcTHDMInterfaceData[struct_String, inputPars_List] :=
    Module[{w, mwStr, higgs, mhStr,
            inPars = Parameters`DecreaseIndexLiterals[inputPars],
            pars = First /@ Reverse /@ inputPars
           },
           w     = Parameters`GetParticleFromDescription["W-Boson"];
           higgs = Parameters`GetParticleFromDescription["Higgs"];
           mwStr = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           mhStr = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[higgs][0]];
           (* compose string *)
           "#ifdef ENABLE_GM2CALC\n" <>
           "GM2Calc_THDM_data " <> struct <> ";\n" <>
           "{\n" <>
              TextFormatting`IndentText[Parameters`CreateLocalConstRefs[pars]] <> "\n" <>
              TextFormatting`IndentText[
                 struct <> ".alpha_em_mz = ALPHA_EM_MZ;\n" <>
                 struct <> ".alpha_em_0 = ALPHA_EM_0;\n" <>
                 struct <> ".alpha_s_mz = ALPHA_S_MZ;\n" <>
                 "if (!is_zero(" <> mhStr <> ")) {\n" <>
                    TextFormatting`IndentText[struct <> ".mh = " <> mhStr <> ";"] <> "\n" <>
                 "} else if (!is_zero(MHPole)) {\n" <>
                    TextFormatting`IndentText[struct <> ".mh = MHPole;"] <> "\n}\n" <>
                 "if (!is_zero(" <> mwStr <> ")) {\n" <>
                    TextFormatting`IndentText[struct <> ".mw = " <> mwStr <> ";"] <> "\n" <>
                 "} else if (!is_zero(MWPole)) {\n" <>
                    TextFormatting`IndentText[struct <> ".mw = MWPole;"] <> "\n}\n" <>
                 struct <> ".mz = MZPole;\n" <>
                 struct <> ".mu(0) = MU2GeV;\n" <>
                 struct <> ".mu(1) = MCMC;\n" <>
                 struct <> ".mu(2) = MTPole;\n" <>
                 struct <> ".md(0) = MD2GeV;\n" <>
                 struct <> ".md(1) = MS2GeV;\n" <>
                 struct <> ".md(2) = MBMB;\n" <>
                 struct <> ".mv(0) = Mv1Pole;\n" <>
                 struct <> ".mv(1) = Mv2Pole;\n" <>
                 struct <> ".mv(2) = Mv3Pole;\n" <>
                 struct <> ".ml(0) = MEPole;\n" <>
                 struct <> ".ml(1) = MMPole;\n" <>
                 struct <> ".ml(2) = MTauPole;\n" <>
                 struct <> ".ckm = CKMInput;\n" <>
                 struct <> ".yukawa_type = " <> CConversion`RValueToCFormString[FlexibleSUSY`yukawaType /. inPars] <> ";\n" <>
                 struct <> ".lambda(0) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda1 /. inPars] <> ";\n" <>
                 struct <> ".lambda(1) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda2 /. inPars] <> ";\n" <>
                 struct <> ".lambda(2) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda3 /. inPars] <> ";\n" <>
                 struct <> ".lambda(3) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda4 /. inPars] <> ";\n" <>
                 struct <> ".lambda(4) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda5 /. inPars] <> ";\n" <>
                 struct <> ".lambda(5) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda6 /. inPars] <> ";\n" <>
                 struct <> ".lambda(6) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda7 /. inPars] <> ";\n" <>
                 struct <> ".tan_beta = " <> CConversion`RValueToCFormString[FlexibleSUSY`tanBeta /. inPars] <> ";\n" <>
                 struct <> ".m122 = " <> CConversion`RValueToCFormString[FlexibleSUSY`m122 /. inPars] <> ";\n" <>
                 struct <> ".zeta_u = " <> CConversion`RValueToCFormString[FlexibleSUSY`zetau /. inPars] <> ";\n" <>
                 struct <> ".zeta_d = " <> CConversion`RValueToCFormString[FlexibleSUSY`zetad /. inPars] <> ";\n" <>
                 struct <> ".zeta_l = " <> CConversion`RValueToCFormString[FlexibleSUSY`zetal /. inPars] <> ";\n" <>
                 struct <> ".delta_u = " <> CConversion`RValueToCFormString[FlexibleSUSY`deltau /. inPars] <> ";\n" <>
                 struct <> ".delta_d = " <> CConversion`RValueToCFormString[FlexibleSUSY`deltad /. inPars] <> ";\n" <>
                 struct <> ".delta_l = " <> CConversion`RValueToCFormString[FlexibleSUSY`deltal /. inPars] <> ";\n" <>
                 struct <> ".pi_u = " <> CConversion`RValueToCFormString[FlexibleSUSY`piu /. inPars] <> ";\n" <>
                 struct <> ".pi_d = " <> CConversion`RValueToCFormString[FlexibleSUSY`pid /. inPars] <> ";\n" <>
                 struct <> ".pi_l = " <> CConversion`RValueToCFormString[FlexibleSUSY`pil /. inPars] <> ";\n"
              ] <>
           "}\n" <>
           "#endif\n\n"
    ];

FillInterfaceData[{}] := "";

FillInterfaceData[obs_List] :=
    Module[{filled = ""},
           If[MemberQ[obs,FlexibleSUSYObservable`aMuonGM2Calc] ||
              MemberQ[obs,FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
              filled = filled <> FillGM2CalcInterfaceData["gm2calc_data"];
             ];
           filled
          ];

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?IsObservable :> a, {0, Infinity}];
           FillInterfaceData[observables] <> "\n" <>
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

End[];

EndPackage[];
