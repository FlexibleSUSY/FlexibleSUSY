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

BeginPackage["Traces`", {"SARAH`", "BetaFunction`", "CConversion`", "Parameters`"}];

CreateTraceRules::usage="takes a list of traces and returns a list of
rules to replace the traces by their C/C++ variables.";

ConvertSARAHTraces::usage="takes SARAH's `TraceAbbr' and returns a
two-component list, where the first entry is string of C/C++ variable
definitions that hold the trace values.  The second entry is a list of
rules to replace the traces by their C/C++ variables.";

CreateTraceDefs::usage="";
CreateTraceCalculation::usage="";
CreateSARAHTraceDefs::usage="";
CreateSARAHTraceCalculation::usage="";
CreateLocalCopiesOfTraces::usage="";
CreateLocalCopiesOfSARAHTraces::usage="";
FindSARAHTraces::usage="";

SARAHTrace;

Begin["`Private`"];

GetSARAHTraceName[Traces`SARAHTrace[name_, expr_]] := name;
GetSARAHTraceExpr[Traces`SARAHTrace[name_, expr_]] := expr;

ConvertToScalar[expr_] :=
    expr /. {
        SARAH`MatMul[a_, b_] :> SARAH`ScalarProd[a, b],
        SARAH`MatMul[a_, b__] :> SARAH`ScalarProd[a, SARAH`MatMul[b]]
            };

FindSARAHTraces[expr_, sarahTraces_List] :=
    Module[{traceSymbols},
           traceSymbols = DeleteDuplicates[Flatten[GetSARAHTraceName /@ sarahTraces]];
           Select[traceSymbols, (!FreeQ[expr,#])&]
          ];

ConvertSARAHTraces[abbrs_] :=
    Module[{i, j, name, expr, traces = {}},
           For[i = 1, i <= Length[abbrs], i++,
               For[j = 1, j <= Length[abbrs[[i]]], j++,
                   {name, expr} = abbrs[[i,j]];
                   (* replace MatMul by ScalarProd *)
                   expr = ConvertToScalar[expr];
                   AppendTo[traces, Traces`SARAHTrace[name, expr]];
                  ];
              ];
           Return[traces];
          ];

GetTraceType[trace_] :=
    CConversion`ScalarType @
    If[Parameters`IsRealExpression[trace] || Parameters`AllModelParametersAreReal[],
       CConversion`realScalarCType,
       CConversion`complexScalarCType
      ];

GetTraceCType[trace_] :=
    CConversion`CreateCType[GetTraceType[trace]];

FindMultipleTraces[list_List] :=
    Module[{traces = Flatten[Cases[list, trace[__], Infinity]]},
           First /@ Select[Tally[traces], (#[[2]] > 1)&]
          ];

FindAllTraces[list_List] :=
    DeleteDuplicates[Flatten[Cases[list, trace[__], Infinity]]];

(* returns all traces appearing at a given loop level *)
FindAllTracesAt[list_List, loopOrder_Integer] :=
    DeleteDuplicates[Flatten[
        Cases[BetaFunction`GetBeta[#,loopOrder]& /@ list, trace[__], {0,Infinity}]
    ]];

(* returns traces appearing at a given loop level, but not below *)
FindAllTracesOnlyAt[list_List, loopOrder_Integer] :=
    Module[{tracesAtThisLevel, tracesBelow, loopsBelow = Table[i,{i,1,loopOrder-1}]},
           tracesBelow = DeleteDuplicates @ Flatten[FindAllTracesAt[list,#]& /@ loopsBelow];
           tracesAtThisLevel = FindAllTracesAt[list, loopOrder];
           Complement[tracesAtThisLevel, tracesBelow]
          ];

CreateTraceRules[traces_List] :=
    (Rule[#, CConversion`ToValidCSymbol[#]])& /@ FindAllTraces[traces];

CreateLocalCopiesOfSARAHTraces[expr_, sarahTraces_List, structName_String] :=
    Module[{traces = FindSARAHTraces[expr, sarahTraces], traceExprs},
           tracesAndExprs = Select[sarahTraces, MemberQ[traces, GetSARAHTraceName[#]]&];
           StringJoin[("const " <> GetTraceCType[GetSARAHTraceExpr[#]] <>
                       " " <> CConversion`ToValidCSymbolString[GetSARAHTraceName[#]] <>
                       " = " <> structName <> "." <>
                       CConversion`ToValidCSymbolString[GetSARAHTraceName[#]] <>
                       ";\n")& /@ tracesAndExprs]
          ];

CreateLocalCopiesOfTraces[list_List, structName_String] :=
    StringJoin[("const " <> GetTraceCType[#] <>
               " " <> CConversion`ToValidCSymbolString[#] <>
               " = " <> structName <> "." <>
               CConversion`ToValidCSymbolString[#] <> ";\n")& /@ FindAllTraces[list]];

CreateTraceDefs[list_List] :=
    StringJoin[(GetTraceCType[#] <> " " <> CConversion`ToValidCSymbolString[#] <> "{};\n")& /@ FindAllTraces[list]];

CreateSARAHTraceDefs[list_List] :=
    StringJoin[(GetTraceCType[GetSARAHTraceExpr[#]] <> " " <>
                CConversion`ToValidCSymbolString[GetSARAHTraceName[#]] <> "{};\n")& /@ list];

CreateCastedTraceExprStr[expr_] :=
    CConversion`CastTo[RValueToCFormString[expr], GetTraceType[expr]];

CreateTraceCalculation[list_List, structName_String] :=
    Module[{traces, Def},
           traces = {FindAllTracesOnlyAt[list,1], FindAllTracesOnlyAt[list,2], FindAllTracesOnlyAt[list,3]};
           Def[expr_] := (structName <> "." <> CConversion`ToValidCSymbolString[expr] <>
                          " = " <> CreateCastedTraceExprStr[expr] <> ";\n");
           StringJoin /@ Map[Def, traces, {2}]
          ];

CreateSARAHTraceCalculation[list_List, structName_String] :=
    Parameters`CreateLocalConstRefsForInputParameters[list] <>
    "\n" <>
    StringJoin[(structName <> "." <>
               CConversion`ToValidCSymbolString[GetSARAHTraceName[#]] <>
               " = " <> CreateCastedTraceExprStr[GetSARAHTraceExpr[#]] <>
               ";\n")& /@ list];

End[];

EndPackage[];
