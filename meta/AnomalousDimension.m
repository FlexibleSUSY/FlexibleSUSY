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

BeginPackage["AnomalousDimension`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Parameters`"}];

AnomalousDimension[];

CreateAnomDimFunctions::usage="";
CreateAnomDimPrototypes::usage="";
ConvertSarahAnomDim::usage="";

Begin["`Private`"];

GetName[AnomalousDimension[name_, type_, anom_List]] := name;

GetType[AnomalousDimension[name_, type_, anom_List]] := type;

GetAnomDim1Loop[AnomalousDimension[name_, type_, anom_List]] := anom[[1]];

GetAnomDim2Loop[AnomalousDimension[name_, type_, anom_List]] := anom[[2]];

GetAllAnomDims[AnomalousDimension[name_, type_, anom_List]] := anom;

CreateValidAnomDimName[names_List /; Length[names] == 2] :=
    Module[{name1, name2},
           name1 = GetHead[names[[1]]];
           name2 = GetHead[names[[2]]];
           Return[Symbol[ToString[name1] <> ToString[name2]]];
          ];

CreateValidAnomDimName[names_] :=
    Module[{},
           Print["Error: don't know how to create name for anomalous ",
                 "dimension from symbol: ", names];
           Quit[1];
          ];

CreateNameWithIndices[sym_Symbol, _CConversion`ScalarType] :=
    sym;

CreateNameWithIndices[sym_Symbol, _CConversion`VectorType] :=
    sym[Susyno`LieGroups`i1];

CreateNameWithIndices[sym_Symbol, _CConversion`MatrixType] :=
    sym[Susyno`LieGroups`i1, SARAH`i2];

CreateNameWithIndices[sym_, _] := sym;

GuessType[{field1_, field2_}] :=
    Module[{f1, f2, dim1, dim2, type},
           f1 = GetHead[field1];
           f2 = GetHead[field2];
           dim1 = TreeMasses`GetDimension[f1];
           dim2 = TreeMasses`GetDimension[f2];
           Which[dim1 > 1 && dim2 > 1,
                 Which[Parameters`AllModelParametersAreReal[],
                       type = CConversion`MatrixType[CConversion`realScalarCType, dim1, dim2];,
                       TreeMasses`IsRealScalar[f1] || TreeMasses`IsRealScalar[f2],
                       type = CConversion`MatrixType[CConversion`realScalarCType, dim1, dim2];,
                       TreeMasses`IsVector[f1] || TreeMasses`IsVector[f2],
                       type = CConversion`MatrixType[CConversion`realScalarCType, dim1, dim2];,
                       True,
                       type = CConversion`MatrixType[CConversion`complexScalarCType, dim1, dim2];
                      ];,
                 dim1 > 1,
                 type = CConversion`VectorType[CConversion`realScalarCType, dim1];,
                 dim2 > 1,
                 type = CConversion`VectorType[CConversion`realScalarCType, dim2];,
                 True,
                 type = CConversion`ScalarType[CConversion`realScalarCType];
                ];
           Return[type];
          ];

StripIndices[expr_, _CConversion`ScalarType] :=
    expr /. Parameters`StripIndicesRules[{Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4}, 2];

StripIndices[expr_, _CConversion`MatrixType] :=
    expr /. Parameters`StripIndicesRules[{Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4}, 2];

StripIndices[expr_, _CConversion`VectorType] :=
    expr /. Parameters`StripIndicesRules[{Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4}, 1];

(* Converts SARAH anomalous dimensions
 *
 * SARAH format:
 *   { {name1, name2}, one-loop, two-loop }
 *
 * Our format:
 *   AnomalousDimension[cName, type, {one-loop, two-loop}]
 *
 * @param gij list of SARAH-like formated anomalous dimensions
 *)
ConvertSarahAnomDim[gij_List] :=
    Module[{lst = {}, adim, i, name, type, expr, nameWithIndices},
           For[i = 1, i <= Length[gij], i++,
               adim = gij[[i]];
               (* adim[[1]] == {name1,name2}, adim[[2]] == 1-loop anom. dim *)
               name = CreateValidAnomDimName[adim[[1]]];
               type = GuessType[adim[[1]]];
               expr = Drop[adim, 1];
               (* append family indices to make ProtectTensorProducts work *)
               nameWithIndices = CreateNameWithIndices[name, type];
               (* protect tensor products *)
               expr = Simplify /@ ((CConversion`ProtectTensorProducts[#, nameWithIndices])& /@ expr);
               (* strip indices *)
               expr = StripIndices[expr /. Kronecker[Susyno`LieGroups`i1,SARAH`i2] :> CreateUnitMatrix[type], type];
               AppendTo[lst, AnomalousDimension[name, type, expr]];
              ];
           Return[lst];
          ];

(* create anomalous dimension prototypes *)
CreateAnomDimPrototypes[anomDim_AnomalousDimension] :=
    Module[{prototypes, name, type},
           name = ToValidCSymbolString[GetName[anomDim]];
           type = GetType[anomDim];
           prototypes = CreateGetterPrototype[name, CConversion`CreateCType[type]];
           Return[prototypes];
          ];

CreateAnomDimPrototypes[anomDim_List] :=
    Module[{prototypes = ""},
           (prototypes = prototypes <> CreateAnomDimPrototypes[#])& /@ anomDim;
           Return[prototypes];
          ];

CreateAnomDimFunction[anomDim_AnomalousDimension] :=
    Module[{def, body, type, name,
            exprOneLoop, exprTwoLoop, inputParsDecl},
           type = GetType[anomDim];
           name = ToValidCSymbolString[GetName[anomDim]];
           (* one-loop *)
           exprOneLoop = CConversion`oneOver16PiSqr * GetAnomDim1Loop[anomDim];
           body = "\nanomDim = " <> CastTo[RValueToCFormString[exprOneLoop],type] <> ";\n";
           (* two-loop *)
           If[Length[GetAllAnomDims[anomDim]] > 1,
              exprTwoLoop = CConversion`twoLoop * GetAnomDim2Loop[anomDim];
              If[exprTwoLoop =!= 0,
                 body = body <> "\nif (get_loops() > 1) {\n" <>
                        IndentText["anomDim += " <> CastTo[RValueToCFormString[exprTwoLoop],type]] <>
                        ";\n}\n";
                ];
             ];
           inputParsDecl = Parameters`CreateLocalConstRefsForInputParameters[{exprOneLoop, exprTwoLoop}];
           body = CreateDefaultDefinition["anomDim", type] <> ";\n" <>
                  inputParsDecl <>
                  body <>
                  "\nreturn anomDim;\n";
           def  = CConversion`CreateCType[type] <> " CLASSNAME::get_" <>
                  name <> "() const\n{\n" <> IndentText[body] <> "}\n\n";
           Return[def];
          ];

CreateAnomDimFunctions[anomDim_List] :=
    Module[{functions = ""},
           (functions = functions <> CreateAnomDimFunction[#])& /@ anomDim;
           Return[functions];
          ];

End[];

EndPackage[];
