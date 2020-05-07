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

BeginPackage["FSMathLink`", {"CConversion`", "Parameters`", "Utils`"}];

GetNumberOfInputParameterRules::usage = "";
GetNumberOfSpectrumEntries::usage = "";
PutInputParameters::usage = "";
SetInputParametersFromArguments::usage = "";
SetInputParameterDefaultArguments::usage = "";
SetInputParameterArguments::usage = "";
PutSpectrum::usage = "";
PutObservables::usage = "";

Begin["`Private`"];

GetNumberOfInputParameterRules[inputPars_List] :=
    Length[inputPars];

PutInputParameter[{par_, _}, linkName_String] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "MLPutRuleTo(" <> linkName <> ", " <>
           "INPUTPARAMETER(" <> parStr <> "), \"" <> parStr <> "\");\n"
          ];

PutInputParameters[inputPars_List, linkName_String] :=
    StringJoin[PutInputParameter[#, linkName]& /@ inputPars];

CreateComponent[CConversion`realScalarCType | CConversion`integerScalarCType, pars_String:"pars", count_String:"c"] :=
    pars <> "[" <> count <> "++];\n";
CreateComponent[CConversion`complexScalarCType, pars_String:"pars", count_String:"c"] :=
    CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
    "(" <> pars <> "[" <> count <> "], " <> pars <> "[" <> count <> "+1]); " <> count <> " += 2;\n";

SetInputParameterFromArguments[{par_, CConversion`ScalarType[st_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           "INPUTPARAMETER(" <> parStr <> ") = " <> CreateComponent[st]
          ];

SetInputParameterFromArguments[{par_, (CConversion`ArrayType | CConversion`VectorType)[st_,dim_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[("INPUTPARAMETER(" <> parStr <> "(" <> # <> ")) = " <>
                       CreateComponent[st])& /@ Table[ToString[i], {i, 0, dim-1}]]
          ];

SetInputParameterFromArguments[{par_, CConversion`MatrixType[st_,dim1_,dim2_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[Flatten[Outer[("INPUTPARAMETER(" <> parStr <> "(" <> #1 <> "," <> #2 <> ")) = " <>
                                     CreateComponent[st])&,
                                    Table[ToString[i], {i, 0, dim1-1}],
                                    Table[ToString[j], {j, 0, dim2-1}]], 1]]
          ];

SetInputParameterFromArguments[{par_, CConversion`TensorType[st_,dim1_,dim2_,dim3_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[Flatten[Outer[("INPUTPARAMETER(" <> parStr <> "(" <> #1 <> "," <> #2 <> "," <> #3 <> ")) = " <>
                                     CreateComponent[st])&,
                                    Table[ToString[i], {i, 0, dim1-1}],
                                    Table[ToString[j], {j, 0, dim2-1}],
                                    Table[ToString[k], {k, 0, dim3-1}]], 2]]
          ];

SetInputParameterFromArguments[{par_, CConversion`TensorType[st_,dim1_,dim2_,dim3_,dim4_]}] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par]},
           StringJoin[Flatten[Outer[("INPUTPARAMETER(" <> parStr <> "(" <> #1 <> "," <> #2 <> "," <> #3 <> "," <> #4 <> ")) = " <>
                                     CreateComponent[st])&,
                                    Table[ToString[i], {i, 0, dim1-1}],
                                    Table[ToString[j], {j, 0, dim2-1}],
                                    Table[ToString[k], {k, 0, dim3-1}],
                                    Table[ToString[l], {l, 0, dim4-1}]], 3]]
          ];

SetInputParametersFromArguments[inputPars_List] :=
    StringJoin[SetInputParameterFromArguments /@ inputPars];

SetInputParameterDefaultArgument[{par_, _[_,dims___]}] :=
    CConversion`ToValidCSymbol[par] -> Array[0&, {dims}];

SetInputParameterDefaultArguments[inputPars_List] :=
    Utils`StringJoinWithSeparator[ToString[SetInputParameterDefaultArgument[#]]& /@ inputPars, ",\n"];

ConcatIndices[idx__] := StringJoin[ToString /@ {idx}];

ParAndType[par_, CConversion`realScalarCType   ] := {HoldForm[OptionValue[par]], Real};
ParAndType[par_, CConversion`complexScalarCType] := Sequence[{HoldForm[Re[OptionValue[par]]], Real},
                                                             {HoldForm[Im[OptionValue[par]]], Real}];
ParAndType[par_, CConversion`integerScalarCType] := {HoldForm[OptionValue[par]], Integer};

ParAndType[par_, CConversion`realScalarCType, idx__   ] := {HoldForm[OptionValue[par][[idx]]], Real};
ParAndType[par_, CConversion`complexScalarCType, idx__] := Sequence[{HoldForm[Re[OptionValue[par][[idx]]]], Real},
                                                                    {HoldForm[Im[OptionValue[par][[idx]]]], Real}];
ParAndType[par_, CConversion`integerScalarCType, idx__] := {HoldForm[OptionValue[par][[idx]]], Integer};

SetInputParameterArgumentsAndType[{par_, CConversion`ScalarType[st_]}] :=
    {ParAndType[CConversion`ToValidCSymbol[par], st]};

SetInputParameterArgumentsAndType[{par_, (CConversion`ArrayType | CConversion`VectorType)[st_, dim_]}] :=
    Table[ParAndType[CConversion`ToValidCSymbol[par], st, i], {i, 1, dim}];

SetInputParameterArgumentsAndType[{par_, CConversion`MatrixType[st_, dim1_, dim2_]}] :=
    Flatten[Outer[ParAndType[CConversion`ToValidCSymbol[par], st, #1, #2] &,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}]], 1];

SetInputParameterArgumentsAndType[{par_, CConversion`TensorType[st_, dim1_, dim2_, dim3_]}] :=
    Flatten[Outer[ParAndType[CConversion`ToValidCSymbol[par], st, #1, #2, #] &,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}],
                  Table[k, {k, 1, dim3}]], 2];

SetInputParameterArgumentsAndType[{par_, CConversion`TensorType[st_, dim1_, dim2_, dim3_, dim4_]}] :=
    Flatten[Outer[ParAndType[CConversion`ToValidCSymbol[par], st, #1, #2, #3, #4] &,
                  Table[i, {i, 1, dim1}],
                  Table[j, {j, 1, dim2}],
                  Table[k, {k, 1, dim3}],
                  Table[l, {l, 1, dim4}]], 3];

SetInputParameterArgumentsAndTypes[inputPars_List] :=
    Join @@ SetInputParameterArgumentsAndType /@ inputPars;

SetInputParameterArguments[{}] := "";

SetInputParameterArguments[inputPars_List] :=
    ",\n" <> Utils`StringJoinWithSeparator[ToString[#[[1]]]& /@ SetInputParameterArgumentsAndTypes[inputPars], ",\n"];

GetNumberOfSpectrumEntries[pars_List] :=
    Length[pars];

(* returns all heads of a nested expression of the form f[g[h[x]]] -> {f,g,h} *)
GetHeads[h_[p___]] := Join[{h}, GetHeads[p]];
GetHeads[p___] := {};

RemoveNamespaces[strs_] :=
    StringReplace[strs, __ ~~ "`" .. -> ""];

HeadStr[par_] :=
    Module[{heads = RemoveNamespaces[ToString /@ GetHeads[par]]},
           If[heads === {}, "",
              ", {\"" <> Utils`StringJoinWithSeparator[heads, "\", \""] <> "\"}"
             ]
          ];

ToUTF8String[s_] :=
    StringJoin[("\\u" <> Utils`FSStringPadLeft[IntegerString[#,16], 4, "0"])& /@ ToCharacterCode[ToString[s]]];

ToValidWolframSymbolString[par_?CConversion`GreekQ] := ToUTF8String[par];
ToValidWolframSymbolString[par_] := ToString[par];

ToValidOutputParStr[FlexibleSUSY`Pole[par_]] := ToValidOutputParStr[par]; (* Pole[x] is not a valid parameter name *)
ToValidOutputParStr[p:FlexibleSUSY`M[_]] := CConversion`ToValidCSymbolString[p];
ToValidOutputParStr[FlexibleSUSY`SCALE] := "scale";
ToValidOutputParStr[par_] := CConversion`ToValidCSymbolString[par];

WithoutHeads[_[par_]] := WithoutHeads[par];
WithoutHeads[par_] := par;

GetMacroFor[FlexibleSUSY`Pole[par_]] := "PHYSICALPARAMETER";
GetMacroFor[par_] := "MODELPARAMETER";

PutParameter[par_, _, link_String] :=
    Module[{parStr = ToValidOutputParStr[par],
            parWithoutHeads = ToValidWolframSymbolString[WithoutHeads[par]]},
           "MLPutRuleTo(" <> link <> ", " <> GetMacroFor[par] <> "(" <> parStr <>
           "), \"" <> parWithoutHeads <> "\"" <> HeadStr[par] <> ");\n"
          ];

PutParameter[par_, link_String] :=
    PutParameter[par, Parameters`GetType[par /. FlexibleSUSY`Pole -> Identity], link];

PutSpectrum[pars_List, link_String] :=
    StringJoin[PutParameter[#,link]& /@ pars];

ObsToStr[obs_?NumericQ] := ToString[obs + 1]; (* +1 to convert to Mathematica's index convention *)
ObsToStr[obs_] := "\"" <> ToString[obs] <> "\"";

HeadToStr[sym_]    := "\"" <> ToString[sym] <> "\"";
HeadsToStr[{}]     := "";
HeadsToStr[l_List] := ", {" <> StringJoin[Riffle[HeadToStr /@ l, ", "]] <> "}";

PutObservable[obs_[sub_], type_, link_String, heads_:{}] :=
    PutObservable[sub, type, link, Join[heads, {obs}]];

PutObservable[obs_, type_, link_String, heads_:{}] :=
    "MLPutRuleTo(" <> link <> ", OBSERVABLE(" <>
    Observables`GetObservableName[Composition[Sequence @@ heads][obs]] <>
    "), " <> ObsToStr[obs] <> HeadsToStr[heads] <> ");\n";

PutObservables[obs_List, link_String] :=
    StringJoin[PutObservable[#, Observables`GetObservableType[#], link]& /@ obs];

End[];

EndPackage[];
