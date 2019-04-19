
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

BeginPackage["BtoSGamma`",
   {"SARAH`", "TextFormatting`", "TreeMasses`", "CXXDiagrams`", "Utils`"}
];

CreateInterfaceBtoSGamma::usage = "";

Begin["`Private`"];


(** \brief Creates the code for the b_to_s_gamma.cpp.in template
 * \note We assume that in the quark multiplet the b quark is the 3rd
 * and the s quark is the 2nd element.
 **)
CreateInterfaceBtoSGamma[matchingOn_] :=
  Module[{definition="",
          inFermion = SARAH`DownQuark,
          outFermion = SARAH`DownQuark,
          dimensionInFermion = TreeMasses`GetDimension[inFermion],
          dimensionOutFermion = TreeMasses`GetDimension[outFermion]},
    Utils`AssertWithMessage[BooleanQ[matchingOn],
      "BtoSGamma`CreateInterfaceBtoSGamma[]: Error, argument must be either True or False."];

    If[matchingOn,
      definition = IndentText[If[dimensionInFermion =!= 1,
        "constexpr int b_quark_index = 2;\n",
        ""] <>
        If[dimensionOutFermion =!= 1,
        "constexpr int s_quark_index = 1;\n",
        ""] <>
        "const auto mb = context.mass<" <> CXXNameOfField[inFermion] <> ">(" <>
        If[dimensionInFermion =!= 1,
          "{2}",
          "{}"] <> ");\n" <>
        "const auto mW = context.mass<" <> CXXNameOfField[TreeMasses`GetWBoson[]] <> ">(" <>
        If[TreeMasses`GetDimension[TreeMasses`GetWBoson[]] =!= 1,
          "{0}",
          "{}"] <> ");\n" <>
        "const auto form_factors_VP = calculate_" <> CXXNameOfField[inFermion] <> "_" <>
        CXXNameOfField[inFermion] <> "_" <> CXXNameOfField[SARAH`Photon] <> "_form_factors(" <>
        If[dimensionInFermion =!= 1,
            "b_quark_index, s_quark_index, model);\n",
            "model);\n"] <>
        "const auto form_factors_VG = calculate_" <> CXXNameOfField[inFermion] <> "_" <>
        CXXNameOfField[inFermion] <> "_" <> CXXNameOfField[SARAH`Gluon] <> "_form_factors(" <>
        If[dimensionInFermion =!= 1,
            "b_quark_index, s_quark_index, model);\n\n",
            "model);\n\n"] <>
        "// choose normalization according to flavio\n" <>
        "const auto normC7 = -16*Sqr(Pi)*Sqr(mW)/(unit_charge(context)*Sqr(g2)*Vtb*Conj(Vts));\n" <>
        "const auto normC8 = -16*Sqr(Pi)*Sqr(mW)/(g3*Sqr(g2)*Vtb*Conj(Vts));\n" <>
        "const auto C7NP_bs = normC7 * form_factors_VP[3];\n" <>
        "const auto C7pNP_bs = normC7 * form_factors_VP[2];\n" <>
        "const auto C8NP_bs = normC8 * form_factors_VG[3];\n" <>
        "const auto C8pNP_bs = normC8 * form_factors_VG[2];\n\n" <>
        "write_wilsoncoeffs(C7NP_bs, C7pNP_bs, C8NP_bs, C8pNP_bs, model.get_scale());\n\n" <>
        "return {C7NP_bs, C7pNP_bs, C8NP_bs, C8pNP_bs};\n"
        ],
      definition = IndentText["return {0, 0, 0, 0};"]
    ];

    definition
  ];

End[];
EndPackage[];
