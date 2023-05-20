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

BeginPackage["Phases`", {"TextFormatting`", "CConversion`"}];

ConvertSarahPhases::usage="convert SARAH phases to our format";

CreatePhasesDefinition::usage="creates definitions of the phases";

CreatePhasesGetterSetters::usage="creates function definitions
for phase getters and setters";

CreatePhaseName::usage = "creates the name of a SARAH phase";

GetArg::usage = "returs the real arg of a given phase";

ClearPhases::usage="resets phases";

Begin["`Private`"];

GetArg[Exp[phase_]] := GetArg[phase];
GetArg[I phase_] := phase;
GetArg[phase_] := phase;

CreatePhaseName[phase_] :=
    CConversion`ToValidCSymbolString[GetArg[phase]];

ConvertSarahPhases[phases_List] :=
    DeleteDuplicates[(#[[2]])& /@ phases];

GetPhaseType[Exp[_]] := CConversion`ScalarType[CConversion`realScalarCType];
GetPhaseType[_]      := CConversion`ScalarType[CConversion`complexScalarCType];

CreatePhasesDefinition[name_String, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]] <> " " <>
    name <> "{1.,0.};\n";

CreatePhasesDefinition[name_String, type_] :=
    Parameters`CreateParameterDefinitionAndDefaultInitialize[{ name, type }];

CreatePhasesDefinition[phases_List] :=
    StringJoin[
       CreatePhasesDefinition[CreatePhaseName[#], GetPhaseType[#]]& /@ phases
    ];

CreatePhasesGetterSetters[phases_List] :=
    Module[{result = "", k, phaseName, type},
           For[k = 1, k <= Length[phases], k++,
               phaseName = CreatePhaseName[phases[[k]]];
               type = GetPhaseType[phases[[k]]];
               result = result <>
                        CConversion`CreateInlineSetter[phaseName, type] <>
                        CConversion`CreateInlineGetter[phaseName, phaseName, type];
              ];
           Return[result];
          ];

ClearPhase[p:Exp[phase_]] :=
    Phases`CreatePhaseName[p] <> " = 0.;\n";

ClearPhase[phase_] :=
    Phases`CreatePhaseName[phase] <> " = " <>
    CreateCType[CConversion`ScalarType[complexScalarCType]] <> "(1.,0.);\n";

ClearPhases[phases_List] :=
    StringJoin[ClearPhase /@ phases];

End[];

EndPackage[];
