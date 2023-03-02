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

BeginPackage["Unitarity`", {"SARAH`", "Parameters`", "CConversion`", "TreeMasses`"}];

GetScatteringMatrix::usage = "";

Begin["`Private`"];

(* return a C++ lambda computing a given scattering eigenvalue in function of sqrt(s) *)
ExpressionToCPPLambda[expr__, istates_List, fstates_List] := Module[{params = Parameters`FindAllParametersClassified[expr], paramsCPP, mixingCPP},

   (* CPP definitions of parameters present in the expression *)
   paramsCPP =
      StringJoin[
         ("const auto " <> ToString@CConversion`ToValidCSymbol[#] <>
            " = model_.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "();\n")& /@ (Parameters`FSModelParameters /. params)
      ];
   (* definition of mixing matrices *)
   mixingCPP = ("auto " <> ToString[#] <> " = [&model_] (int i, int j) { return model_.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "(i,j); };\n")& /@ (Parameters`FSOutputParameters /. params);

   (* replace input parameters with their FS names *)
   newExpr = expr /. Thread[(Parameters`FSModelParameters /. params) -> CConversion`ToValidCSymbol /@ (Parameters`FSModelParameters /. params)];
   newExpr = newExpr /. Susyno`LieGroups`conj -> Conj;
   newExpr = ToString@CForm[newExpr];
   (* in expressions from SARAH masses are denoted as pmass(X)
      this converts them to proper C++ form *)
   newExpr = StringReplace[newExpr, "pmass(" ~~ f:Except[")"].. ~~ "(" ~~ i:Except[")"].. ~~ "))" :> "context.mass<" <> f <> ">({" <> ToString[i] <> "})"];
   newExpr = StringReplace[newExpr, "pmass(" ~~ f:Except[")"]..  ~~ ")" :> "context.mass<" <> f <> ">({})"];

If[expr === 0,
"[](double sqrtS, int in1, int in2, int out1, int out2) { return 0.; };\n\n",
"[&model" <> If[!FreeQ[expr, sChan], ",sChan", ""] <> If[!FreeQ[expr, tChan], ",tChan", ""] <> If[!FreeQ[expr, uChan], ",uChan", ""] <> If[!FreeQ[expr, qChan], ",qChan", ""] <> "](double sqrtS, int in1, int in2, int out1, int out2) {
      auto model_ = model;
      // couplings should be evaluated at the renormalization scale sqrt(s)
      // see comment around eq. 20 of 1805.07306
      model_.run_to(sqrtS);
      model_.solve_ewsb();
      const double s = Sqr(sqrtS);
      const " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::context_base context {model_};\n" <>
      "if (sqrtS < context.mass<" <> ToString[fstates[[1]] /. Susyno`LieGroups`conj->Identity] <> ">(" <> If[TreeMasses`GetDimension[fstates[[1]]]>1, "{out1}", "{}"] <> ") + context.mass<" <> ToString[fstates[[2]] /. Susyno`LieGroups`conj->Identity] <> ">(" <> If[TreeMasses`GetDimension[fstates[[2]]]>1, "{out2}", "{}"] <> ")) return 0.;\n" <>
      paramsCPP <>
      mixingCPP <>
      (* TODO: can this really be complex? *)
"     return " <> If[FreeQ[expr, Susyno`LieGroups`conj], newExpr, "std::real(" <> newExpr <> ")"] <> ";
   };\n\n"
]
];

GetScatteringMatrix[] := Module[{result},
   InitUnitarity[];
   (*RemoveParticlesFromScattering={Se , Sv, Sd, Su};*)

   a0 = Outer[GetScatteringDiagrams[#1 -> #2]&, scatteringPairs, scatteringPairs, 1];

   result = "";
   For[i=1, i<=Length[a0], i++,
      For[j=i, j<=Length[a0[[i]]], j++,
         result = result <> "// " <> ToString[scatteringPairs[[i]]] <> "->" <> ToString[scatteringPairs[[j]]] <> "\n" <>
                  "matrix[" <> ToString[i-1] <> "][" <> ToString[j-1] <> "] = " <> ExpressionToCPPLambda[a0[[i,j]], scatteringPairs[[i]], scatteringPairs[[j]]]
      ]
   ];
   {SparseArray[a0]["NonzeroPositions"], Dimensions[a0][[1]], result}
];

End[];
EndPackage[];

