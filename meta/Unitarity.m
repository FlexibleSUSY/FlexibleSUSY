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

replacePMASS[expr_] := Module[{temp},
   temp = ToString@CForm[expr];
   (* in expressions from SARAH masses are denoted as pmass(X)
      this converts them to proper C++ form *)
   temp = StringReplace[temp, "pmass(" ~~ f:Except[")"].. ~~ "(" ~~ i:DigitCharacter.. ~~ "))" :> "context.mass<" <> f <> ">({" <> ToString[ToExpression[i]-1] <> "})"];
   temp = StringReplace[temp, "pmass(" ~~ f:Except[")"].. ~~ "(" ~~ i:Except[")"].. ~~ "))" :> "context.mass<" <> f <> ">({" <> ToString[i] <> "-1})"];
   temp = StringReplace[temp, "pmass(" ~~ f:Except[")"]..  ~~ ")" :> "context.mass<" <> f <> ">({})"];
   temp
];

InfiniteS[a0Input_, generationSizes_] := Module[{params = Parameters`FindAllParametersClassified[a0Input], paramsCPP, mixingCPP, a0},
   (* CPP definitions of parameters present in the expression *)
   paramsCPP =
      StringJoin[
         ("const auto " <> ToString@CConversion`ToValidCSymbol[#] <>
            " = model.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "();\n")& /@ (Parameters`FSModelParameters /. params)
      ];
   (* definition of mixing matrices *)
   mixingCPP = ("auto " <> ToString[#] <> " = [&model] (int i, int j) { return model.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "(i-1,j-1); };\n")& /@ (Parameters`FSOutputParameters /. params);

   (* replace input parameters with their FS names *)
   a0 = a0Input /. Thread[(Parameters`FSModelParameters /. params) -> CConversion`ToValidCSymbol /@ (Parameters`FSModelParameters /. params)];
   resultInfinite = "";
   For[i=1, i<=Length[a0], i++,
      For[j=i, j<=Length[a0[[i]]], j++,
         If[a0[[i,j]] =!= 0,
            resultInfinite = resultInfinite <> "// " <> ToString[scatteringPairs[[i]]] <> "->" <> ToString[scatteringPairs[[j]]] <> "\n" <>
                     If[generationSizes[[i,j,1]] > 1, "for (int in1 = 1; in1 <= " <> ToString[generationSizes[[i,j,1]]] <> "; ++in1) {\n", ""] <>
                     If[generationSizes[[i,j,2]] > 1, "for (int in2 = 1; in2 <= " <> ToString[generationSizes[[i,j,2]]] <> "; ++in2) {\n", ""] <>
                     If[generationSizes[[i,j,3]] > 1, "for (int out1 = 1; out1 <= " <> ToString[generationSizes[[i,j,3]]] <> "; ++out1) {\n", ""] <>
                     If[generationSizes[[i,j,4]] > 1, "for (int out2 = 1; out2 <= " <> ToString[generationSizes[[i,j,4]]] <> "; ++out2) {\n", ""] <>
                     "matrix.coeffRef(" <> ToString[i-1] <> ", " <> ToString[j-1] <> ") = std::max(std::abs(" <> ToString@CForm[a0[[i,j]]] <> "), matrix.coeff(" <> ToString[i-1] <> ", " <> ToString[j-1] <> "));\n" <>
                     If[generationSizes[[i,j,1]] > 1, "}\n", ""] <>
                     If[generationSizes[[i,j,2]] > 1, "}\n", ""] <>
                     If[generationSizes[[i,j,3]] > 1, "}\n", ""] <>
                     If[generationSizes[[i,j,4]] > 1, "}\n", ""]
         ]
      ]
   ];
   paramsCPP <> mixingCPP <> resultInfinite
];

(* return a C++ lambda computing a given scattering eigenvalue in function of sqrt(s) *)
ExpressionToCPPLambda[expr__, istates_List, fstates_List] := Module[{params = Parameters`FindAllParametersClassified[expr], paramsCPP, mixingCPP},

   (* CPP definitions of parameters present in the expression *)
   paramsCPP =
      StringJoin[
         ("const auto " <> ToString@CConversion`ToValidCSymbol[#] <>
            " = model_.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "();\n")& /@ (Parameters`FSModelParameters /. params)
      ];
   (* definition of mixing matrices *)
   mixingCPP = ("auto " <> ToString[#] <> " = [&model_] (int i, int j) { return model_.get_" <> ToString@CConversion`ToValidCSymbol[#] <> "(i-1,j-1); };\n")& /@ (Parameters`FSOutputParameters /. params);

   (* replace input parameters with their FS names *)
   newExpr = expr /. Thread[(Parameters`FSModelParameters /. params) -> CConversion`ToValidCSymbol /@ (Parameters`FSModelParameters /. params)];
   newExpr = newExpr /. Susyno`LieGroups`conj -> Conj;
   newExpr = newExpr //. SARAH`sum[idx_, start_, stop_, exp_] :> FlexibleSUSY`SUM[idx, start, stop, exp];

"[&model" <> If[!FreeQ[expr, sChan], ",sChan", ""] <> If[!FreeQ[expr, tChan], ",tChan", ""] <> If[!FreeQ[expr, uChan], ",uChan", ""] <> If[!FreeQ[expr, qChan], ",qChan", ""] <> "](double sqrtS, int in1, int in2, int out1, int out2) {
      auto model_ = model;
      // couplings should be evaluated at the renormalization scale sqrt(s)
      // see comment around eq. 20 of 1805.07306
      model_.run_to(sqrtS);
      model_.solve_ewsb();
      const double s = Sqr(sqrtS);
      const " <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::context_base context {model_};\n" <>
      "if (sqrtS < context.mass<" <> ToString[fstates[[1]] /. Susyno`LieGroups`conj->Identity] <> ">(" <> If[TreeMasses`GetDimension[fstates[[1]]]>1, "{out1-1}", "{}"] <> ") + context.mass<" <> ToString[fstates[[2]] /. Susyno`LieGroups`conj->Identity] <> ">(" <> If[TreeMasses`GetDimension[fstates[[2]]]>1, "{out2-1}", "{}"] <> ")) return 0.;\n" <>
      paramsCPP <>
      mixingCPP <>
      (* TODO: can this really be complex? *)
"double res = 0.;\n" <>
With[{temp = Simplify@Coefficient[newExpr, qChan]},
If[temp =!= 0,
"if (qChan) {
   res += " <> replacePMASS@temp <> ";
}\n",
""
]] <>
With[{temp = Coefficient[newExpr, sChan]},
If[temp =!= 0,
"if (sChan) {
   res += " <> If[FreeQ[temp, Susyno`LieGroups`conj], "", "std::real("] <> replacePMASS@temp <> If[FreeQ[temp, Susyno`LieGroups`conj], "", ")"] <> ";
}\n",
""
]] <>
With[{temp = Coefficient[newExpr, tChan]},
If[temp =!= 0,
"if (tChan) {
   res += " <> If[FreeQ[temp, Susyno`LieGroups`conj], "", "std::real("] <> replacePMASS@temp <> If[FreeQ[temp, Susyno`LieGroups`conj], "", ")"] <> ";
}\n",
""
]] <>
With[{temp = Coefficient[newExpr, uChan]},
If[temp =!= 0,
"if (uChan) {
   res += " <> If[FreeQ[temp, Susyno`LieGroups`conj], "", "std::real("] <> replacePMASS@temp <> If[FreeQ[temp, Susyno`LieGroups`conj], "", ")"] <> ";
}\n",
""
]] <>
"     return " <> If[FreeQ[expr, Susyno`LieGroups`conj], "res", "std::real(res)"] <> ";
   };\n\n"
];

GetScatteringMatrix[] := Module[{result, generationSizes, a0, a0InfiniteS},
   InitUnitarity[];
   (*RemoveParticlesFromScattering={Se , Sv, Sd, Su};*)

   a0 = Outer[GetScatteringDiagrams[#1 -> #2]&, scatteringPairs, scatteringPairs, 1];
   a0InfiniteS = Map[Limit[# /. qChan -> 1 /. sChan|tChan|uChan -> 0, Susyno`LieGroups`s->Infinity]&, a0, {2}];
   generationSizes = Table[{i, j}, {i,1, Length[scatteringPairs]}, {j,1, Length[scatteringPairs]}];
   generationSizes = Apply[Join[TreeMasses`GetDimension /@ scatteringPairs[[#1]], TreeMasses`GetDimension /@ scatteringPairs[[#2]]]&, generationSizes, {2}];
   resultFinite = "";
   resultInfinite = "";
   For[i=1, i<=Length[a0], i++,
      For[j=i, j<=Length[a0[[i]]], j++,
         If[a0[[i,j]] =!= 0,
            resultFinite = resultFinite <> "// " <> ToString[scatteringPairs[[i]]] <> "->" <> ToString[scatteringPairs[[j]]] <> "\n" <>
                     "matrix[" <> ToString[i-1] <> "][" <> ToString[j-1] <> "] = " <> ExpressionToCPPLambda[a0[[i,j]], scatteringPairs[[i]], scatteringPairs[[j]]];
            resultInfinite = resultInfinite <> "// " <> ToString[scatteringPairs[[i]]] <> "->" <> ToString[scatteringPairs[[j]]] <> "\n" <>
                     "matrix[" <> ToString[i-1] <> "][" <> ToString[j-1] <> "] = " <> ToString@CForm[a0InfiniteS[[i,j]]] <> ";\n";
         ]
      ]
   ];
   {SparseArray[a0]["NonzeroPositions"], Dimensions[a0][[1]], generationSizes, InfiniteS[a0InfiniteS, generationSizes], resultFinite}
];

End[];
EndPackage[];

