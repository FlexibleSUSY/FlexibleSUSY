(* ::Package:: *)

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

BeginPackage["BrLToLGamma`", 
   {"SARAH`", "TextFormatting`", "TreeMasses`", "CXXDiagrams`", "CConversion`"}
];

CreateInterfaceFunctionForBrLToLGamma::usage = "";

Begin["`Private`"];

CreateInterfaceFunctionForBrLToLGamma[inFermion_ -> {outFermion_, spectator_}] :=
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
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input);";

      definition =
            (* calculate observable using form factors *)
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> " (\n" <>
               IndentText[
                  If[TreeMasses`GetDimension[inFermion] =!= 1, "int generationIndex1, ", ""] <>
                  If[TreeMasses`GetDimension[outFermion] =!= 1, "int generationIndex2, ", ""] <> "\n" <>
                  "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input) {\n\n"
               ] <>
               (* write routine for mu to e gamma *)
               IndentText[
                  "context_base context {model};\n" <>
                  "std::array<int, " <> ToString @ numberOfIndices1 <>
                     "> indices1 = {" <>
                     (* TODO: Specify indices correctly *)
                       If[TreeMasses`GetDimension[inFermion] =!= 1,
                          " generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ","] <> " ",
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
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ","] <> " ",
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
EndPackage[];
