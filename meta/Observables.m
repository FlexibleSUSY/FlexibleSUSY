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

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { aMuon, aMuonUncertainty, aMuonGM2Calc, aMuonGM2CalcUncertainty,
                  CpHiggsPhotonPhoton, CpHiggsGluonGluon,
                  CpPseudoScalarPhotonPhoton, CpPseudoScalarGluonGluon,
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
           If[MemberQ[observables, FlexibleSUSYObservable`CpHiggsPhotonPhoton] ||
              MemberQ[observables, FlexibleSUSYObservable`CpHiggsGluonGluon],
              dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson]
              If[FreeQ[TreeMasses`GetParticles[], SARAH`HiggsBoson] ||
                 TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson] == 0,
                 Print["Warning: no physical Higgs boson found."];
                 Print["         Effective couplings for Higgs boson will not"];
                 Print["         be calculated."];
                 observables = DeleteCases[observables, a_ /; (a === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                                                               a === FlexibleSUSYObservable`CpHiggsGluonGluon)];
                ];
             ];
           If[MemberQ[observables, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] ||
              MemberQ[observables, FlexibleSUSYObservable`CpPseudoScalarGluonGluon],
              If[FreeQ[TreeMasses`GetParticles[], SARAH`PseudoScalar] ||
                 TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar] == 0,
                 Print["Warning: no physical pseudoscalar boson found."];
                 Print["         Effective couplings for pseudoscalar boson will not"];
                 Print["         be calculated."];
                 observables = DeleteCases[observables, a_ /; (a === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                                                               a === FlexibleSUSYObservable`CpPseudoScalarGluonGluon)];
                ];
             ];
test =       Complement[
              Cases[observables, _FlexibleSUSYObservable`BrLToLGamma],
              Cases[observables, FlexibleSUSYObservable`BrLToLGamma[fin_?IsLepton -> {fout_?IsLepton, vout_ /; vout === GetPhoton[]}]]
                 ];
           If[test =!= {},
              Print["Warning: BrLToLGamma function works only for leptons and a photon."];
              Print["         Removing requested process(es):"];
              Print["        " <> ToString@test];
              observables = Complement[observables, test];
           ];
           observables
          ];

GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuon] := "a_muon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonUncertainty] := "a_muon_uncertainty";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon_gm2calc";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "a_muon_gm2calc_uncertainty";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := "eff_cp_higgs_photon_photon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := "eff_cp_higgs_gluon_gluon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := "eff_cp_pseudoscalar_photon_photon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := "eff_cp_pseudoscalar_gluon_gluon";
GetObservableName[FlexibleSUSYObservable`EDM[p_[idx_]]] := GetObservableName[FlexibleSUSYObservable`EDM[p]] <> "_" <> ToString[idx];
GetObservableName[FlexibleSUSYObservable`EDM[p_]]       := "edm_" <> CConversion`ToValidCSymbolString[p];
GetObservableName[FlexibleSUSYObservable`BrLToLGamma[pIn_[_] -> {pOut_[_], spectator_}]] := CConversion`ToValidCSymbolString[pIn] <> "_to_" <> CConversion`ToValidCSymbolString[pOut] <> "_" <> CConversion`ToValidCSymbolString[spectator];
GetObservableName[FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, spectator_}]] := CConversion`ToValidCSymbolString[pIn] <> "_to_" <> CConversion`ToValidCSymbolString[pOut] <> "_" <> CConversion`ToValidCSymbolString[spectator];
GetObservableName[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_[idxIn_] -> pOut_[idxOut_], nucleus_]] := CConversion`ToValidCSymbolString[pIn] <> "_to_" <> CConversion`ToValidCSymbolString[pOut] <> "_in_" <> ToString@nucleus;
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := "b_to_s_gamma";

GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuon] := "a_muon = (g-2)/2 of the muon (calculated with FlexibleSUSY)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonUncertainty] := "uncertainty of a_muon = (g-2)/2 of the muon (calculated with FlexibleSUSY)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon = (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "uncertainty of (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := "effective H-Photon-Photon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := "effective H-Gluon-Gluon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := "effective A-Photon-Photon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := "effective A-Gluon-Gluon coupling";
GetObservableDescription[FlexibleSUSYObservable`EDM[p_[idx_]]] := "electric dipole moment of " <> CConversion`ToValidCSymbolString[p] <> "(" <> ToString[idx] <> ") [1/GeV]";
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

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuon] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonUncertainty] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`EDM[p_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, _}]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[FlexibleSUSYObservable`FToFConversionInNucleus[pIn_[idxIn_] -> pOut_[idxOut_], _]] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`bsgamma] := CConversion`ScalarType[CConversion`realScalarCType];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

CountNumberOfObservables[observables_List] :=
    Module[{i, number = 0},
           For[i = 1, i <= Length[observables], i++,
               If[IsObservable[observables[[i]]],
                  number += BetaFunction`CountNumberOfParameters[GetObservableType[observables[[i]]]];,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
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
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
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
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
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
                     Print["Warning: ignoring invalid observable ", observables[[i]]];
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
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuon, structName_String] :=
    structName <> ".AMU = " <> FlexibleSUSY`FSModelName <> "_a_muon::calculate_a_muon(MODEL, qedqcd);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonUncertainty, structName_String] :=
    structName <> ".AMUUNCERTAINTY = " <> FlexibleSUSY`FSModelName <> "_a_muon::calculate_a_muon_uncertainty(MODEL, qedqcd);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    "#ifdef ENABLE_GM2Calc\n" <>
    structName <> ".AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);\n" <>
    "#endif";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty, structName_String] :=
    "#ifdef ENABLE_GM2Calc\n" <>
    structName <> ".AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);\n" <>
    "#endif";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton, structName_String] :=
    Module[{i, type, dim, start, result = ""},
           type = GetObservableType[obs];
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim != 1,
              start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`HiggsBoson] - 1;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPHIGGSPHOTONPHOTON("
                           <> ToString[i-1] <> ") = effective_couplings.get_eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "("
                           <> ToString[start+i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim = TreeMasses`GetDimension[SARAH`HiggsBoson];
              If[dim == 1,
                 result = structName <> ".EFFCPHIGGSPHOTONPHOTON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "();",
                 start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`HiggsBoson] - 1;
                 result = structName <> ".EFFCPHIGGSPHOTONPHOTON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "("
                          <> ToString[start] <> ");"
                ];
             ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon, structName_String] :=
    Module[{i, type, dim, start, result = ""},
           type = GetObservableType[obs];
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim != 1,
              start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`HiggsBoson] - 1;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPHIGGSGLUONGLUON("
                           <> ToString[i-1] <> ") = effective_couplings.get_eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "("
                           <> ToString[start+i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim = TreeMasses`GetDimension[SARAH`HiggsBoson];
              If[dim == 1,
                 result = structName <> ".EFFCPHIGGSGLUONGLUON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "();",
                 start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`HiggsBoson] - 1;
                 result = structName <> ".EFFCPHIGGSGLUONGLUON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "("
                          <> ToString[start] <> ");"
                ];
             ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton, structName_String] :=
    Module[{i, type, dim, start, result = ""},
           type = GetObservableType[obs];
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim != 1,
              start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`PseudoScalar] - 1;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPPSEUDOSCALARPHOTONPHOTON("
                           <> ToString[i-1] <> ") = effective_couplings.get_eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "("
                           <> ToString[start+i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim = TreeMasses`GetDimension[SARAH`PseudoScalar];
              If[dim == 1,
                 result = structName <> ".EFFCPPSEUDOSCALARPHOTONPHOTON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "();",
                 start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`PseudoScalar] - 1;
                 result = structName <> ".EFFCPPSEUDOSCALARPHOTONPHOTON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "("
                          <> ToString[start] <> ");"
                ];
             ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon, structName_String] :=
    Module[{i, type, dim, start, result = ""},
           type = GetObservableType[obs];
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim != 1,
              start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`PseudoScalar] - 1;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPPSEUDOSCALARGLUONGLUON("
                           <> ToString[i-1] <> ") = effective_couplings.get_eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "("
                           <> ToString[start+i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim = TreeMasses`GetDimension[SARAH`PseudoScalar];
              If[dim == 1,
                 result = structName <> ".EFFCPPSEUDOSCALARGLUONGLUON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "();",
                 start = TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`PseudoScalar] - 1;
                 result = structName <> ".EFFCPPSEUDOSCALARGLUONGLUON = effective_couplings.get_eff_Cp"
                          <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                          <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "("
                          <> ToString[start] <> ");"
                ];
             ];
           result
          ];

CalculateObservable[FlexibleSUSYObservable`EDM[p_], structName_String] :=
    Module[{pStr = CConversion`ToValidCSymbolString[p]},
           structName <> ".EDM0(" <> pStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_edm::calculate_edm_" <> pStr <> "(MODEL);"
          ];

CalculateObservable[FlexibleSUSYObservable`EDM[p_[idx_]], structName_String] :=
    Module[{pStr = CConversion`ToValidCSymbolString[p],
            idxStr = ToString[idx]},
           structName <> ".EDM1(" <> pStr <> ", " <> idxStr <> ") = " <>
           FlexibleSUSY`FSModelName <> "_edm::calculate_edm_" <> pStr <> "(" <> idxStr <> ", MODEL);"
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

FillGM2CalcInterfaceData[struct_String] :=
    Module[{filling, mwStr,
            w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
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

           If[MemberQ[{w, pseudoscalar, smuon, muonsneutrino,
                       chargino, neutralino, mu, m1, m2, m3, mq2, mu2,
                       md2, ml2, me2, tu, td, te, yu, yd, ye}, Null],
              Print["Error: The GM2Calc addon cannot be used in this model, because it is not a MSSM-like model with sfermion flavour conservation. ",
                    "Please remove aMuonGM2Calc and aMuonGM2CalcUncertainty from the model file."];
              Quit[1];
           ];

           mwStr         = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           filling = \
           struct <> ".alpha_s_MZ = ALPHA_S_MZ;\n" <>
           struct <> ".MZ    = MZPole;\n" <>
           "if (!is_zero(" <> mwStr <> "))\n" <>
              TextFormatting`IndentText[struct <> ".MW = " <> mwStr <> ";"] <> "\n" <>
           "else if (!is_zero(MWPole))\n" <>
              TextFormatting`IndentText[struct <> ".MW = MWPole;"] <> "\n" <>
           struct <> ".mb_mb = MBMB;\n" <>
           struct <> ".MT    = MTPole;\n" <>
           struct <> ".MTau  = MTauPole;\n" <>
           struct <> ".MM    = MMPole;\n" <>
           struct <> ".MA0   = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[pseudoscalar][1]] <> ";\n" <>
           struct <> ".MSvm  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[muonsneutrino]] <> ";\n" <>
           struct <> ".MSm   = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[smuon]] <> ";\n" <>
           struct <> ".MCha  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[chargino]] <> ";\n" <>
           struct <> ".MChi  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[neutralino]] <> ";\n" <>
           struct <> ".scale = MODEL.get_scale();\n" <>
           struct <> ".TB    = MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> "() / " <>
                              "MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> "();\n" <>
           struct <> ".Mu    = MODEL.get_" <> CConversion`RValueToCFormString[mu] <> "();\n" <>
           struct <> ".M1    = MODEL.get_" <> CConversion`RValueToCFormString[m1] <> "();\n" <>
           struct <> ".M2    = MODEL.get_" <> CConversion`RValueToCFormString[m2] <> "();\n" <>
           struct <> ".M3    = MODEL.get_" <> CConversion`RValueToCFormString[m3] <> "();\n" <>
           struct <> ".mq2   = MODEL.get_" <> CConversion`RValueToCFormString[mq2] <> "();\n" <>
           struct <> ".mu2   = MODEL.get_" <> CConversion`RValueToCFormString[mu2] <> "();\n" <>
           struct <> ".md2   = MODEL.get_" <> CConversion`RValueToCFormString[md2] <> "();\n" <>
           struct <> ".ml2   = MODEL.get_" <> CConversion`RValueToCFormString[ml2] <> "();\n" <>
           struct <> ".me2   = MODEL.get_" <> CConversion`RValueToCFormString[me2] <> "();\n" <>
           struct <> ".Au    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[tu] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[yu] <> "());\n" <>
           struct <> ".Ad    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[td] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[yd] <> "());\n" <>
           struct <> ".Ae    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[te] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[ye] <> "());";
           "#ifdef ENABLE_GM2Calc\n" <>
           "GM2Calc_data " <> struct <> ";\n" <> filling <> "\n" <>
           "#endif\n\n"
          ];

FillEffectiveCouplingsInterfaceData[struct_String] :=
    Module[{result},
           result = FlexibleSUSY`FSModelName <> "_effective_couplings " <> struct <> "(model, qedqcd, physical_input);\n";
           result = result <> "effective_couplings.calculate_effective_couplings();\n"
          ];

FillInterfaceData[{}] := "";

FillInterfaceData[obs_List] :=
    Module[{filled = ""},
           If[MemberQ[obs,FlexibleSUSYObservable`aMuonGM2Calc] ||
              MemberQ[obs,FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
              filled = filled <> FillGM2CalcInterfaceData["gm2calc_data"];
             ];
           If[MemberQ[obs,FlexibleSUSYObservable`CpHiggsPhotonPhoton]         ||
              MemberQ[obs,FlexibleSUSYObservable`CpHiggsGluonGluon]           ||
              MemberQ[obs, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] ||
              MemberQ[obs, FlexibleSUSYObservable`CpPseudoScalarGluonGluon],
              filled = filled <> FillEffectiveCouplingsInterfaceData["effective_couplings"];
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
