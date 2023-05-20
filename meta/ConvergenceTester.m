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

BeginPackage["ConvergenceTester`", {"CConversion`", "TextFormatting`", "TreeMasses`", "Parameters`", "Utils`"}];

CreateCompareFunction::usage="";

Begin["`Private`"];

CountNumberOfParameters[FlexibleSUSY`M[particle_]] :=
    TreeMasses`GetDimension[particle];

CountNumberOfParameters[parameters_List] :=
    Plus @@ (CountNumberOfParameters /@ parameters);

CountNumberOfParameters[parameter_[__]] :=
    If[Parameters`IsRealParameter[parameter], 1, 2];

CountNumberOfParameters[parameter_] :=
    If[Parameters`IsRealParameter[parameter], 1, 2] * (Times @@ Parameters`GetParameterDimensions[parameter]);

(* Maps tensor indices to linear space *)
IndexMapping[{_}, {i1_String}] :=
    i1;

IndexMapping[d_List, i_List] :=
    IndexMapping[Take[d,Length[d]-1], Take[i,Length[i]-1]] <>
    " + " <> ToString[Times @@ Take[d,Length[d]-1]] <> "*" <> Last[i];

CalcDifference[FlexibleSUSY`M[particle_], offset_Integer, diff_String] :=
    Module[{result, body, dim, dimStart, esStr},
           dim = TreeMasses`GetDimension[particle];
           esStr = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[dim == 1,
              result = diff <> "[" <> ToString[offset] <> "] = " <>
                       "MaxRelDiff(OLD(" <> esStr <> "),NEW(" <> esStr <> "));\n";
              ,
              dimStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle] - 1;
              result = "for (int i = " <> ToString[dimStart] <>
                       "; i < " <> ToString[dim] <> "; ++i) {\n";
              body = diff <> "[i + " <> ToString[offset] <> "] = " <>
                     "MaxRelDiff(OLD1(" <> esStr <> ",i),NEW1(" <> esStr <> ",i));";
              result = result <> IndentText[body] <> "\n}\n";
             ];
           Return[result];
          ];

CalcDifference[parameter_[], offset_Integer, diff_String] :=
    Module[{result, parStr},
           parStr = ToValidCSymbolString[parameter];
           result = diff <> "[" <> ToString[offset] <> "] = " <>
                    "MaxRelDiff(OLD(" <> parStr <> "),NEW(" <> parStr <> "));\n";
           Return[result];
          ];

CalcDifference[parameter_[idx___Integer], offset_Integer, diff_String] :=
    Module[{result, parStr, dim, dimStr},
           dim = Length[{idx}];
           dimStr = ToString[dim];
           parStr = ToValidCSymbolString[parameter] <> "," <>
                    Utils`StringJoinWithSeparator[ToString[#-1]& /@ {idx},","];
           result = diff <> "[" <> ToString[offset] <> "] = " <>
                    "MaxRelDiff(OLD" <> dimStr <> "(" <> parStr <> "),NEW" <> dimStr <> "(" <> parStr <> "));\n";
           Return[result];
          ];

CalcDifference[parameter_ /; Parameters`GetParameterDimensions[parameter] == {1}, offset_Integer, diff_String] :=
    CalcDifference[parameter[], offset, diff];

CalcDifference[parameter_, offset_Integer, diff_String, {idx_List, pos_Integer, idxPool_List} /; pos > Length[idx]] :=
    Module[{body, dim, dimStr, parStr},
           dim = Length[idx];
           dimStr = ToString[dim];
           parStr = ToValidCSymbolString[parameter] <> "," <>
                    Utils`StringJoinWithSeparator[ToString /@ idxPool,","];
           body = diff <> "[" <> IndexMapping[idx,idxPool] <> " + " <> ToString[offset] <> "] = " <>
                  "MaxRelDiff(OLD" <> dimStr <> "(" <> parStr <> "),NEW" <> dimStr <> "(" <> parStr <> "));";
           Return[body];
          ];

CalcDifference[parameter_, offset_Integer, diff_String, {idx_List, pos_Integer, idxPool_List}] :=
    Module[{result, dim, dimStr, i},
           dim = Length[idx];
           dimStr = ToString[dim];
           i = idxPool[[pos]];
           result = "for (int " <> i <> " = 0; " <> i <> " < " <> ToString[idx[[pos]]] <> "; ++" <> i <> ") {\n" <>
                    IndentText[
                        CalcDifference[parameter, offset, diff, {idx, pos+1, idxPool}]
                    ] <>
                    "\n}";
           Return[result];
          ];

CalcDifference[parameter_, offset_Integer, diff_String] :=
    Module[{result, body, dim, parStr, idxPool},
           dim = Parameters`GetParameterDimensions[parameter];
           idxPool = Take[{"i", "j", "k", "l","m","n"}, Length[dim]];
           result = CalcDifference[parameter, offset, diff, {dim, 1, idxPool}] <> "\n";
           Return[result];
          ];

CreateCompareFunction[crit_ /; crit === Automatic] :=
    Module[{particles},
           If[SARAH`SupersymmetricModel,
              particles = TreeMasses`GetSusyParticles[];,
              particles = TreeMasses`GetParticles[];
             ];
           particles = Select[particles, (!TreeMasses`IsMassless[#] && !IsGhost[#] && !IsGoldstone[#])&];
           particles = FlexibleSUSY`M /@ particles;
           CreateCompareFunction[particles]
          ];

CreateCompareFunction[parameters_List] :=
    Module[{result, numberOfParameters, i, offset = 0},
           numberOfParameters = CountNumberOfParameters[parameters];
           If[numberOfParameters == 0,
              Print["Error: no parameters specified for the convergence test!"];
              Return["return 0.;"];
             ];
           ctype = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           result = "std::array<" <> ctype <> ", " <> ToString[numberOfParameters] <> "> diff{};\n\n";
           For[i = 1, i <= Length[parameters], i++,
               result = result <> CalcDifference[parameters[[i]], offset, "diff"];
               offset += CountNumberOfParameters[parameters[[i]]];
              ];
           If[offset != numberOfParameters,
              Print["Error: something is wrong with the counting of masses:"];
              Print["  numberOfParameters = ", numberOfParameters, ", offset = ", offset];
             ];
           result = result <>
                    "\nreturn *std::max_element(diff.cbegin(), diff.cend());\n";
           Return[result];
          ];

End[];

EndPackage[];
