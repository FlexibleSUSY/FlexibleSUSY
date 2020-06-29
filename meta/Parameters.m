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

BeginPackage["Parameters`", {"SARAH`", "CConversion`", "Utils`", "Phases`"}];

{ FSModelParameters, FSInputParameters, FSOutputParameters,
  FSPhysicalOutputParameters, FSPhases, FSDerivedParameters,
  FSExtraParameters };

ParameterDimensions::usage="option flag for setting the dimensions of an
input or extra parameter";
MassDimension::usage="option flag for setting the mass dimension of an
input or extra parameter";
InputParameter::usage="option flag for indicating parameter should be
treated as an input parameter";

FindSymbolDef::usage="";

CreateParameterDefinition::usage="";
CreateParameterDefinitionAndDefaultInitialize::usage="";
CreateSetAssignment::usage="";
CreateDisplayAssignment::usage="";
CreateParameterSARAHNames::usage="";
CreateParameterEnums::usage="";
CreateInputParameterEnum::usage="";
CreateInputParameterNames::usage="";
CreateExtraParameterEnum::usage="";
CreateExtraParameterNames::usage="";
CreateStdVectorNamesAssignment::usage="";

CreateExtraParameterArrayGetter::usage="";
CreateExtraParameterArraySetter::usage="";
CreateInputParameterArrayGetter::usage="";
CreateInputParameterArraySetter::usage="";
CreateModelParameterGetter::usage="";
CreateModelParameterSetter::usage="";
CreateDelegateModelParameterGetter::usage="";

CreateEnumName::usage="Creates enum symbol for given parameter";
DecomposeParameter::usage="decomposes parameter into its real components";

SetParameter::usage="set model parameter";
SetSMParameter::usage="sets a SM input parameter in the QedQcd class";
SetInputParameter::usage="set input parameter to value";
AddInputParameters::usage="add an input parameter";
AddExtraParameters::usage="add an extra parameter";
RemoveExtraParameters::usage="removes a parameter from the list of
known extra parameters";
SetPhases::usage="sets field phases";
GetPhases::usage="returns field phases";
SetPhase::usage="sets a phase to a value";

ApplyGUTNormalization::usage="Returns a list of repacement rules for
gauge couplings, which replace non-normalized gauge couplings
(e.g. gY) by normalized ones (e.g. g1).";

GetGUTNormalization::usage="Returns GUT normalization of the given
coupling";

CreateIndexReplacementRules::usage="";

AddRealParameter::usage="";
SetRealParameters::usage="";

SaveParameterLocally::usage="Save parameters in local variables";

GetType::usage="";
GetPhase::usage="";
HasPhase::usage="";
GetIntegerTypeFromDimension::usage="";
GetRealTypeFromDimension::usage="";
GetComplexTypeFromDimension::usage="";
GetParameterDimensions::usage="";
GetThirdGeneration::usage="returns parameter with third generation index";

GuessInputParameterType::usage="returns a guess for the type of the parameter";
GuessExtraParameterType::usage="returns a guess for the type of the parameter";

IsRealParameter::usage="";
IsComplexParameter::usage="";
IsRealExpression::usage="";
IsMatrix::usage="returns True if parameter is a matrix";
IsSymmetricMatrixParameter::usage="returns True if parameter is a matrix";
IsTensor::usage="returns True if parameter is a matrix";
IsParameter::usage="returns True if symbol is a model/input/output/phase parameter";
IsModelParameter::usage="returns True if parameter is a model parameter";
IsInputParameter::usage="returns True if parameter is an input parameter";
IsOutputParameter::usage="returns True if parameter is a defined output parameter";
IsIndex::usage="returns True if given symbol is an index";
IsPhase::usage="returns True if given symbol is a phase";
IsExtraParameter::usage="return True if parameter is an auxiliary parameter";
IsGaugeCoupling::usage="returns True if parameter is a gauge coupling.";
IsYukawaCoupling::usage="returns True if parameter is a Yukawa coupling.";
IsVEV::usage="returns True if parameter is a VEV.";

GetIndices::usage="returns list of indices from a given parameter";

AllModelParametersAreReal::usage="returns True if all model parameters
are real, False otherwise";

SetInputParameters::usage="";
SetModelParameters::usage="";
SetOutputParameters::usage="";
SetExtraParameters::usage="";

ApplyAuxiliaryParameterInfo::usage="Saves input and extra parameter properties"

CheckInputParameterDefinitions::usage="";

GetInputParameters::usage="";
GetInputParametersAndTypes::usage="";
GetModelParameters::usage="";
GetOutputParameters::usage="";
GetExtraParameters::usage="";
GetExtraParametersAndTypes::usage="";
GetModelParametersWithMassDimension::usage="Returns model parameters
with given mass dimension";
GetParametersWithMassDimension::usage="Returns model, input and extra
parameters with given mass dimension";
GetModelParameterMassDimension::usage="Returns mass dimension for given
model parameter";

FindMacro::usage="Returns preprocessor macro for parameter";
WrapPreprocessorMacroAround::usage="Applies preprocessor symbols
to parameters";

GetDependenceSPhenoSymbols::usage="Returns list of SARAH parameters
 for which a DependenceSPheno rule is defined";

GetDependenceSPhenoRules::usage="Returns list of replacement rules for
 SARAH parameters for which a DependenceSPheno rule is defined";

GetAllDependenceSPhenoRules::usage="Returns list of replacement rules
 for all DependenceSPheno rules"

GetOutputParameterDependencies::usage="Returns list of output
 parameters which appear in the given expression";

GetIntermediateOutputParameterDependencies::usage="Returns list of
 intermediate output parameters which appear in the given expression";

CreateLocalConstRefs::usage="creates local const references to model
parameters / input parameters.";

CreateLocalConstRefsForBetas::usage="";

CreateLocalConstRefsForInputParameters::usage="creates local const
references for all input parameters in the given expression.";

CreateLocalConstRefsForPhysicalParameters::usage="creates local const
references for all physical output parameters in the given
expression.";

FillInputParametersFromTuples::usage="";

DecreaseIndexLiterals::usage="@note
Definetly safe calls:
	f[ { Head1[Int1] , Head2[Int2], ... } ]
	f[ { Head1[Int1] , Head2[Int2], ... }, { HeadI1, HeadI2, ...} ]
Possible calls:
	1. f[ exprD ]
	2. f[ exprD, listH ]
where
	#1 exprD is expression which probably contains subexpressions (at any level) of the form HeadOfExpression[SomeInteger]
	#2 listH is List which has subexpressions (only at lowest level) of the form	HeadOfExpression
Output:
	1. Gets names ({HeadOfExpression1, HeadOfExpression2, ...}) of all InputParameters, ExtraParameters, allModelParameters, allOutputParameters and changes
	   in exprD any HeadOfExpression1[SomeInteger] to HeadOfExpression1[SomeInteger-1]
	2. Gets names ({HeadOfExpression1, HeadOfExpression2, ...}) in listH and changes
	   in exprD any HeadOfExpression1[SomeInteger] to HeadOfExpression1[SomeInteger-1]
";
IncreaseIndexLiterals::usage="";

DecreaseSumIndices::usage="";

ExpressionToString::usage="Converts an expression to a valid C++
string.";

GetEffectiveMu::usage="";
GetEffectiveMASqr::usage="";

GetParameterFromDescription::usage="Returns model parameter from a
given description string.";
GetParticleFromDescription::usage="Returns particle symbol from a
given description string.";
GetPDGCodesForParticle::usage="Returns the PDG codes for a particle."

NumberOfIndependentEntriesOfSymmetricMatrix::usage="Returns number of
independent parameters of a real symmetric nxn matrix";

ExpandExpressions::usage="";
AppendGenerationIndices::usage="";

StripIndices::usage="removes indices from model parameter";
CreateIndices::usage="creates indices in C++ expression for parameter";

StripIndicesRules::usage="removes given indices from a symbol";

StripSARAHIndicesRules::usage="removes SARAH-specific indices from a symbol";

ReplaceAllRespectingSARAHHeads::usage="applies given rules, respecting
SARAH heads SARAH`B, SARAH`L, SARAH`T and SARAH`Q";

FindAllParameters::usage = "returns list of all parameters contained
in the given expression";

FindAllParametersClassified::usage = "returns list of all parameters
 contained in the given expression, classified by their meaning.";

FindSLHABlock::usage = "returns SLHA input block name for given
 parameter";

Begin["`Private`"];

allInputParameters = {};
allExtraParameters = {};
allModelParameters = {};
allOutputParameters = {};
allPhases = {};

(* list storing mass dimensions for input and extra parameters,    *)
(* in the form {{0, {parameters ...}}, {1, {parameters ...}}, ...} *)
extraMassDimensions = {};

AddMassDimensionInfo[par_, dim_?IntegerQ] :=
    Module[{massDimensions, pos, known},
           massDimensions = First /@ extraMassDimensions;
           If[!MemberQ[massDimensions, dim],
              extraMassDimensions = Utils`ForceJoin[extraMassDimensions, {{dim, {par}}}];,
              pos = Position[massDimensions, dim];
              known = First[Extract[extraMassDimensions, pos]][[2]];
              extraMassDimensions = ReplacePart[extraMassDimensions,
                                                pos -> {dim, DeleteDuplicates[Join[known, {par}]]}];
             ];
          ];

AddMassDimensionInfo[par_, dim_] :=
    Print["Error: mass dimension for parameter ", par,
          " must be a number"];

GuessInputParameterType[Sign[par_]] :=
    CConversion`ScalarType[CConversion`integerScalarCType];
GuessInputParameterType[FlexibleSUSY`Phase[par_]] :=
    CConversion`ScalarType[CConversion`complexScalarCType];
GuessInputParameterType[par_] :=
    CConversion`ScalarType[CConversion`realScalarCType];

GuessExtraParameterType[Sign[par_]] :=
    CConversion`ScalarType[CConversion`integerScalarCType];
GuessExtraParameterType[FlexibleSUSY`Phase[par_]] :=
    CConversion`ScalarType[CConversion`complexScalarCType];
GuessExtraParameterType[par_] :=
    CConversion`ScalarType[CConversion`realScalarCType];

UpdateParameterInfo[currentPars_List, {par_, block_, type_}] :=
    Module[{parNames, pos, updatedPars},
           parNames = First /@ currentPars;
           If[!MemberQ[parNames, par],
              updatedPars = Utils`ForceJoin[currentPars, {{par, block, type}}];,
              pos = Position[parNames, par, 1];
              updatedPars = ReplacePart[currentPars, pos -> {par, block, type}];
             ];
           If[CConversion`GetElementType[type] === CConversion`realScalarCType ||
              CConversion`GetElementType[type] === CConversion`integerScalarCType,
              AddRealParameter[par];
             ];
           Return[updatedPars];
          ];

UpdateParameterInfo[currentPars_List, {par_, type_}] :=
    Module[{parNames, pos, updatedPars},
           parNames = First /@ currentPars;
           If[!MemberQ[parNames, par],
              updatedPars = Utils`ForceJoin[currentPars, {{par, type}}];,
              pos = Position[parNames, par, 1];
              updatedPars = ReplacePart[currentPars, pos -> {par, type}];
             ];
           If[CConversion`GetElementType[type] === CConversion`realScalarCType ||
              CConversion`GetElementType[type] === CConversion`integerScalarCType,
              AddRealParameter[par];
             ];
           Return[updatedPars];
          ];

SetStoredParameterSLHABlock[storedPars_List, par_, block_] :=
    Module[{pos, updated, updatedPars},
           pos = Position[storedPars, {par, __}];
           updated = Extract[storedPars, pos];
           If[MatchQ[block, {_, _}],
              updated = ({First[#], {ToString[block[[1]]], block[[2]]}, #[[3]]})& /@ updated,
              updated = ({First[#], ToString[block], #[[3]]})& /@ updated
             ];
           ReplacePart[storedPars, MapThread[Rule, {pos, updated}]]
          ];

SetStoredParameterType[storedPars_List, par_, type_] :=
    Module[{pos, updated},
           pos = Position[storedPars, {par, __}];
           updated = Extract[storedPars, pos] /. {par, block_, oldType_} :> {par, block, type};
           ReplacePart[storedPars, MapThread[Rule, {pos, updated}]]
          ];

SetStoredParameterDimensions[storedPars_List, par_, dims_] :=
    Module[{pos, updated},
           pos = Position[storedPars, {par, __}];
           updated = Extract[storedPars, pos];
           updated = ({First[#], #[[2]], If[CConversion`IsRealType[#[[3]]],
                                          If[CConversion`IsIntegerType[#[[3]]],
                                             GetIntegerTypeFromDimension[dims],
                                             GetRealTypeFromDimension[dims]
                                            ],
                                          GetComplexTypeFromDimension[dims]]})& /@ updated;
           ReplacePart[storedPars, MapThread[Rule, {pos, updated}]]
          ];

AddInputParameterInfo[{par_, block_, type_}] :=
    allInputParameters = UpdateParameterInfo[allInputParameters, {par, block, type}];

AddInputParameterInfo[{par_, type_}] :=
    allInputParameters = UpdateParameterInfo[allInputParameters, {par, {}, type}];

AddInputParameterInfo[par_] :=
    AddInputParameterInfo[{par, {}, GuessInputParameterType[par]}];

SetInputParameters[pars_List] :=
    (
     allInputParameters = {};
     AddInputParameterInfo /@ pars;
    )

AddInputParameters[pars_List] := AddInputParameterInfo /@ pars;

AddExtraParameterInfo[{par_, block_, type_}] :=
    allExtraParameters = UpdateParameterInfo[allExtraParameters, {par, block, type}];

AddExtraParameterInfo[{par_, type_}] :=
    allExtraParameters = UpdateParameterInfo[allExtraParameters, {par, {}, type}];

AddExtraParameterInfo[par_] :=
    AddExtraParameterInfo[{par, {}, GuessExtraParameterType[par]}];

SetExtraParameters[pars_List] :=
    (
     allExtraParameters = {};
     AddExtraParameterInfo /@ pars;
    )

AddExtraParameters[pars_List] := AddExtraParameterInfo /@ pars;

RemoveExtraParameters[par_] :=
    allExtraParameters = DeleteCases[allExtraParameters, {par, __}];

RemoveExtraParameters[pars_List] := RemoveExtraParameters /@ pars;

SetInputParameterType[par_?IsInputParameter, type_] :=
    allInputParameters = SetStoredParameterType[allInputParameters, par, type];

SetInputParameterType[par_, type_] :=
    Print["Error: ", par, " is not an input parameter!"];

SetInputParameterDimensions[par_?IsInputParameter, dims_] :=
    allInputParameters = SetStoredParameterDimensions[allInputParameters, par, dims];

SetInputParameterDimensions[par_, dims_] :=
    Print["Error: ", par, " is not an input parameter!"];

SetInputParameterSLHABlock[par_?IsInputParameter, block_] :=
    allInputParameters = SetStoredParameterSLHABlock[allInputParameters, par, block];

SetInputParameterSLHABlock[par_, block_] :=
    Print["Error: ", par, " is not an input parameter!"];

SetExtraParameterType[par_?IsExtraParameter, type_] :=
    allExtraParameters = SetStoredParameterType[allExtraParameters, par, type];

SetExtraParameterType[par_, type_] :=
    Print["Error: ", par, " is not a defined parameter!"];

SetExtraParameterDimensions[par_?IsExtraParameter, dims_] :=
    allExtraParameters = SetStoredParameterDimensions[allExtraParameters, par, dims];

SetExtraParameterDimensions[par_, dims_] :=
    Print["Error: ", par, " is not a defined parameter!"];

SetExtraParameterSLHABlock[par_?IsExtraParameter, block_] :=
    allExtraParameters = SetStoredParameterSLHABlock[allExtraParameters, par, block];

SetExtraParameterSLHABlock[par_, block_] :=
    Print["Error: ", par, " is not a defined parameter!"];

ProcessParameterInfo[{parameter_ /; (IsModelParameter[parameter] || IsOutputParameter[parameter]), {__}}] :=
    Block[{},
          Print["Warning: the properties of ", parameter, " are set"];
          Print["   in the SARAH model files and cannot be overridden."];
          Print["   Ignoring property settings for ", parameter];
         ];

ProcessParameterInfo[{parameter_?IsInputParameter, properties_List}] :=
    Module[{i, inputBlock, ignored = {}, validProperties = properties, property, setting},
           validProperties = Select[properties, (First[#] =!= InputParameter)&];
           For[i = 1, i <= Length[validProperties], i++,
               property = validProperties[[i, 1]];
               setting = validProperties[[i, 2]];
               Which[property === ParameterDimensions,
                     SetInputParameterDimensions[parameter, setting],
                     property === MassDimension,
                     AddMassDimensionInfo[parameter, setting],
                     property === SARAH`LesHouches,
                     SetInputParameterSLHABlock[parameter, setting],
                     True, Print["Warning: unrecognized property for parameter ", parameter, ": ", property]
                    ];
              ];
          ];

ProcessParameterInfo[{parameter_?IsExtraParameter, properties_List}] :=
    Module[{i, property, setting},
           For[i = 1, i <= Length[properties], i++,
               property = properties[[i, 1]];
               setting = properties[[i, 2]];
               Which[property === ParameterDimensions,
                     SetExtraParameterDimensions[parameter, setting],
                     property === MassDimension,
                     AddMassDimensionInfo[parameter, setting],
                     property === SARAH`LesHouches,
                     SetExtraParameterSLHABlock[parameter, setting],
                     property === InputParameter,
                     Print["Error: ", parameter, " is defined as an input parameter"];
                     Print["   but is being treated as an extra parameter."];
                     Quit[1];,
                     True, Print["Warning: unrecognized property for parameter ", parameter, ": ", property]
                    ];
              ];
          ];

ProcessParameterInfo[{parameter_, properties_List}] :=
    Block[{inputOptions, isInput = False},
          inputOptions = Select[properties, (First[#] === InputParameter)&];
          If[inputOptions =!= {},
             isInput = Last[Last[inputOptions]];
            ];
          If[isInput,
             AddInputParameterInfo[parameter],
             AddExtraParameterInfo[parameter]
            ];
          ProcessParameterInfo[{parameter, properties}];
         ];

ApplyAuxiliaryParameterInfo[properties_List] :=
    ProcessParameterInfo /@ properties;

CheckInputParameterDefinitions[] :=
    Module[{i, par, blockName, type},
           For[i = 1, i <= Length[allInputParameters], i++,
               par = allInputParameters[[i,1]];
               type = allInputParameters[[i,3]];
               (* with the exception of phases, complex input parameters are not currently supported *)
               If[!CConversion`IsRealType[type] && !MatchQ[par, FlexibleSUSY`Phase[_]],
                  Print["Error: ", par, " is defined to be complex,"];
                  Print["   but input parameters must be real."];
                  Print["   Please define ", par, " to be a real parameter."];
                  Quit[1];
                 ];
               If[MatchQ[allInputParameters[[i,2]], {_, _}],
                  blockName = allInputParameters[[i,2,1]],
                  blockName = allInputParameters[[i,2]]
                 ];
               blockName = ToString[blockName];
               If[blockName === "MINPAR" || blockName === "EXTPAR" || blockName === "IMEXTPAR",
                  If[!MatchQ[type, CConversion`ScalarType[_]],
                     Print["Error: ", par, " must be defined as a scalar"];
                     Print["   since it is defined in an SLHA1 input block."];
                     Quit[1];
                    ];
                  If[MatchQ[par, FlexibleSUSY`Phase[_]],
                     If[type =!= CConversion`ScalarType[CConversion`complexScalarCType],
                        Print["Error: ", par, " must be defined as a complex scalar"];
                        Print["   since it is defined in an SLHA1 input block."];
                        Quit[1];
                       ];,
                     If[type =!= CConversion`ScalarType[CConversion`realScalarCType] &&
                        type =!= CConversion`ScalarType[CConversion`integerScalarCType],
                        Print["Error: ", par, " must be defined as a real scalar"];
                        Print["   since it is defined in an SLHA1 input block."];
                        Quit[1];
                       ];
                    ];
                 ];
              ];
          ];

SetModelParameters[pars_List] := allModelParameters = DeleteDuplicates[pars];
SetOutputParameters[pars_List] := allOutputParameters = DeleteDuplicates[pars];
SetPhases[phases_List]        := allPhases = DeleteDuplicates[phases];

GetInputParameters[] := First /@ allInputParameters;
GetInputParametersAndTypes[] := allInputParameters;
GetModelParameters[] := allModelParameters;
GetOutputParameters[] := allOutputParameters;
GetPhases[] := allPhases;
GetExtraParameters[] := First /@ allExtraParameters;
GetExtraParametersAndTypes[] := allExtraParameters;

additionalRealParameters = {};

AddRealParameter[parameter_List] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, parameter]];

AddRealParameter[parameter_] :=
    additionalRealParameters = DeleteDuplicates[Join[additionalRealParameters, {parameter}]];

SetRealParameters[parameters_List] :=
    additionalRealParameters = parameters;

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<Parameters>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

FindSymbolDef[sym_, opt_:DependenceNum] :=
    Module[{symDef},
           symDef = Cases[SARAH`ParameterDefinitions,
                          {sym, {___, opt -> definition_, ___}} :> definition];
           If[Head[symDef] =!= List || symDef === {},
              Print["Error: Could not find definition of ",
                    sym, " in SARAH`ParameterDefinitions"];
              Return[0];
             ];
           If[Length[symDef] > 1,
              Print["Warning: ", sym, " defined multiple times"];
             ];
           symDef[[1]]
          ];

FindAllParametersFromList[expr_, parameters_List] :=
    Module[{symbols, compactExpr},
           compactExpr = RemoveProtectedHeads[expr];
           (* find all parameters with SARAH head *)
           symbols = DeleteDuplicates[Flatten[
               { Cases[compactExpr, SARAH`L[a_][__] /; MemberQ[parameters,SARAH`L[a]] :> SARAH`L[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`B[a_][__] /; MemberQ[parameters,SARAH`B[a]] :> SARAH`B[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`T[a_][__] /; MemberQ[parameters,SARAH`T[a]] :> SARAH`T[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`Q[a_][__] /; MemberQ[parameters,SARAH`Q[a]] :> SARAH`Q[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`L[a_]     /; MemberQ[parameters,SARAH`L[a]] :> SARAH`L[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`B[a_]     /; MemberQ[parameters,SARAH`B[a]] :> SARAH`B[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`T[a_]     /; MemberQ[parameters,SARAH`T[a]] :> SARAH`T[a], {0,Infinity}],
                 Cases[compactExpr, SARAH`Q[a_]     /; MemberQ[parameters,SARAH`Q[a]] :> SARAH`Q[a], {0,Infinity}]
               }]];
           (* remove parameters found from compactExpr *)
           compactExpr = compactExpr /. (RuleDelayed[#, Evaluate[CConversion`ToValidCSymbolString[#]]]& /@ symbols);
           (* find all parameters without SARAH head in compactExpr *)
           symbols = Join[symbols,
               { Cases[compactExpr, a_Symbol /; MemberQ[parameters,a], {0,Infinity}],
                 Cases[compactExpr, a_[__] /; MemberQ[parameters,a] :> a, {0,Infinity}],
                 Cases[compactExpr, FlexibleSUSY`M[a_]     /; MemberQ[parameters,FlexibleSUSY`M[a]], {0,Infinity}],
                 Cases[compactExpr, FlexibleSUSY`M[a_[__]] /; MemberQ[parameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
               }];
           DeleteDuplicates[Flatten[symbols]]
          ];

(* Returns all parameters within an expression *)
FindAllParameters[expr_, exceptions_:{}] :=
    Module[{allParameters, allOutPars},
           allOutPars = DeleteDuplicates[Flatten[
               Join[allOutputParameters,
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> FlexibleSUSY`M[a],
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> (FlexibleSUSY`M /@ {a})
                   ]]];
           allParameters = DeleteDuplicates[
               Join[allModelParameters, allOutPars,
                    GetInputParameters[], Phases`GetArg /@ allPhases,
                    GetDependenceSPhenoSymbols[], GetExtraParameters[]]];
           allParameters = DeleteCases[allParameters, p_ /; MemberQ[exceptions, p]];
           FindAllParametersFromList[expr, allParameters]
          ];

FindAllParametersClassified[expr_, exceptions_:{}] :=
    Module[{symbols = DeleteDuplicates[Flatten[FindAllParameters[expr, exceptions]]],
            inputPars, modelPars, outputPars, extraPars,
            poleMasses, phases, depNum, allOutPars},
           allOutPars = DeleteDuplicates[Flatten[
               Join[allOutputParameters,
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> FlexibleSUSY`M[a],
                    allOutputParameters /. FlexibleSUSY`M[{a__}] :> (FlexibleSUSY`M /@ {a})
                   ]]];
           poleMasses = {
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_]]     /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}],
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_[__]]] /; MemberQ[allOutputParameters,FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
                        };
           poleMasses   = DeleteDuplicates[Flatten[poleMasses]];
           inputPars    = DeleteDuplicates[Select[symbols, (MemberQ[GetInputParameters[],#])&]];
           modelPars    = DeleteDuplicates[Select[symbols, (MemberQ[allModelParameters,#])&]];
           outputPars   = DeleteDuplicates[Select[symbols, (MemberQ[allOutPars,#])&]];
           phases       = DeleteDuplicates[Select[symbols, (MemberQ[Phases`GetArg /@ allPhases,#])&]];
           depNum       = DeleteDuplicates[Select[symbols, (MemberQ[GetDependenceSPhenoSymbols[],#])&]];
           extraPars    = DeleteDuplicates[Select[symbols, (MemberQ[GetExtraParameters[],#])&]];
           {
               FSModelParameters -> modelPars,
               FSInputParameters -> inputPars,
               FSOutputParameters -> outputPars,
               FSPhysicalOutputParameters -> poleMasses,
               FSPhases -> phases,
               FSDerivedParameters -> depNum,
               FSExtraParameters -> extraPars
           }
          ];

ReplaceAllRespectingSARAHHeads[expr_, rules_] :=
    Module[{pars, parsWithoutHeads, removeHeadsRules,
            uniqueRules, uniqueExpr, uniqueSubs},
           pars = Parameters`FindAllParameters[First /@ rules];
           removeHeadsRules = { SARAH`L[p_][__] :> p, SARAH`L[p_] :> p,
                                SARAH`B[p_][__] :> p, SARAH`B[p_] :> p,
                                SARAH`T[p_][__] :> p, SARAH`T[p_] :> p,
                                SARAH`Q[p_][__] :> p, SARAH`Q[p_] :> p
                              };
           parsWithoutHeads = DeleteDuplicates[pars /. removeHeadsRules];
           uniqueRules = DeleteDuplicates @ Flatten[{
               Rule[SARAH`L[#], CConversion`ToValidCSymbol[SARAH`L[#]]],
               Rule[SARAH`B[#], CConversion`ToValidCSymbol[SARAH`B[#]]],
               Rule[SARAH`T[#], CConversion`ToValidCSymbol[SARAH`T[#]]],
               Rule[SARAH`Q[#], CConversion`ToValidCSymbol[SARAH`Q[#]]]
           }& /@ parsWithoutHeads];
           uniqueExpr = expr /. uniqueRules;
           uniqueSubs = rules /. uniqueRules;
           (uniqueExpr /. uniqueSubs) /. (Reverse /@ uniqueRules)
          ];

IsScalar[sym_] :=
    Length[SARAH`getDimParameters[sym]] === 1 || Length[SARAH`getDimParameters[sym]] == 0;

IsScalar[sym_?IsInputParameter] :=
    MatchQ[GetType[sym], CConversion`ScalarType[_]];

IsScalar[sym_?IsExtraParameter] :=
    MatchQ[GetType[sym], CConversion`ScalarType[_]];

IsMatrix[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsMatrix[sym];

IsMatrix[sym_] :=
    Length[SARAH`getDimParameters[sym]] === 2;

IsMatrix[sym_?IsInputParameter] :=
    MatchQ[GetType[sym], CConversion`MatrixType[__]];

IsMatrix[sym_?IsExtraParameter] :=
    MatchQ[GetType[sym], CConversion`MatrixType[__]];

IsSymmetricMatrixParameter[sym_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsSymmetricMatrixParameter[sym];

IsSymmetricMatrixParameter[sym_] :=
    IsMatrix[sym] && MemberQ[SARAH`ListSoftBreakingScalarMasses, sym];

IsTensor[sym_[Susyno`LieGroups`i1, SARAH`i2, SARAH`i3]] :=
    IsTensor[sym];

IsTensor[sym_?IsInputParameter] :=
    MatchQ[GetType[sym], CConversion`TensorType[__]];

IsTensor[sym_?IsExtraParameter] :=
    MatchQ[GetType[sym], CConversion`TensorType[__]];

IsTensor[sym_] :=
    Length[SARAH`getDimParameters[sym]] > 2;

AllModelParametersAreReal[] := MemberQ[SARAH`RealParameters, All];

sarahIndices = {
    SARAH`gt1, SARAH`gt2, SARAH`gt3, SARAH`gt4,
    Susyno`LieGroups`i1 , SARAH`i2 , SARAH`i3 , SARAH`i4
};

IsIndex[i_?NumberQ] := True;
IsIndex[i_ /; MemberQ[sarahIndices,i]] := True;
IsIndex[_] := False;
IsIndex[indices_List] := And @@ (IsIndex /@ indices);
IsIndex[indices__] := IsIndex[{indices}];

GetIndices[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] := {indices};
GetIndices[parameter_] := {};

IsPhase[parameter_] := MemberQ[Phases`GetArg /@ allPhases, Phases`GetArg[parameter]];

IsModelParameter[parameter_] := MemberQ[allModelParameters, parameter];
IsModelParameter[Re[parameter_]] := IsModelParameter[parameter];
IsModelParameter[Im[parameter_]] := IsModelParameter[parameter];
IsModelParameter[FlexibleSUSY`FSTemporary[parameter_]] := IsModelParameter[parameter];

IsModelParameter[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    IsModelParameter[parameter];

IsInputParameter[parameter_] := MemberQ[GetInputParameters[], parameter];

IsInputParameter[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    IsInputParameter[parameter];

IsOutputParameter[lst_List] := And @@ (IsOutputParameter /@ lst);
IsOutputParameter[sym_]     := MemberQ[GetOutputParameters[],sym];

IsExtraParameter[parameter_] := MemberQ[GetExtraParameters[], parameter];

IsExtraParameter[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    IsExtraParameter[parameter];

IsParameter[sym_] :=
    IsModelParameter[sym] ||
    IsInputParameter[sym] ||
    IsExtraParameter[sym] ||
    IsPhase[sym];

IsRealParameter[Re[sym_]] := True;
IsRealParameter[Im[sym_]] := True;
IsRealParameter[FlexibleSUSY`M[_]] := True;

IsRealParameter[sym_] :=
    (IsModelParameter[sym] && AllModelParametersAreReal[]) ||
    (IsInputParameter[sym] && CConversion`IsRealType[GetType[sym]]) ||
    (IsExtraParameter[sym] && CConversion`IsRealType[GetType[sym]]) ||
    MemberQ[Utils`ForceJoin[SARAH`realVar, additionalRealParameters, SARAH`RealParameters], sym];

IsRealParameter[sym_[indices__?IsIndex]] := IsRealParameter[sym];

IsComplexParameter[sym_] :=
    !IsRealParameter[sym];

IsRealExpression[parameter_?IsModelParameter] :=
    IsRealParameter[parameter];

IsRealExpression[parameter_?IsInputParameter] :=
    IsRealParameter[parameter];

IsRealExpression[parameter_?IsExtraParameter] :=
    IsRealParameter[parameter];

IsRealExpression[expr_?NumericQ] :=
    Element[expr, Reals];

IsRealExpression[_Complex] := False;

IsRealExpression[_Real] := True;

IsRealExpression[expr_[Susyno`LieGroups`i1, SARAH`i2]] :=
    IsRealExpression[expr];

IsRealExpression[HoldPattern[SARAH`Delta[_,_]]] := True;
IsRealExpression[HoldPattern[SARAH`ThetaStep[_,_]]] := True;
IsRealExpression[Cos[_]] := True;
IsRealExpression[Sin[_]] := True;
IsRealExpression[ArcCos[_]] := True;
IsRealExpression[ArcSin[_]] := True;
IsRealExpression[ArcTan[_]] := True;
IsRealExpression[Power[a_,b_]] := IsRealExpression[a] && IsRealExpression[b];
IsRealExpression[Susyno`LieGroups`conj[expr_]] :=
    IsRealExpression[expr];
IsRealExpression[SARAH`Conj[expr_]] := IsRealExpression[expr];
IsRealExpression[Conjugate[expr_]]  := IsRealExpression[expr];
IsRealExpression[Transpose[expr_]]  := IsRealExpression[expr];
IsRealExpression[SARAH`Tp[expr_]]   := IsRealExpression[expr];
IsRealExpression[SARAH`Adj[expr_]]  := IsRealExpression[expr];
IsRealExpression[SARAH`bar[expr_]]  := IsRealExpression[expr];

IsRealExpression[expr_Symbol] := IsRealParameter[expr];

IsRealExpression[Times[Conjugate[a_],a_]] := True;
IsRealExpression[Times[a_,Conjugate[a_]]] := True;
IsRealExpression[Times[Susyno`LieGroups`conj[a_],a_]] := True;
IsRealExpression[Times[a_,Susyno`LieGroups`conj[a_]]] := True;
IsRealExpression[Times[SARAH`Conj[a_],a_]]            := True;
IsRealExpression[Times[a_,SARAH`Conj[a_]]]            := True;

IsRealExpression[factors_Times] :=
    And @@ (IsRealExpression[#]& /@ (List @@ factors));

IsRealExpression[terms_Plus] :=
    And @@ (IsRealExpression[#]& /@ (List @@ terms));

IsHermitian[a_] :=
    Or[IsScalar[a] && IsRealParameter[a],
       IsSymmetricMatrixParameter[a] && IsRealParameter[a],
       MemberQ[SARAH`ListSoftBreakingScalarMasses, a]
      ];

(* helper function which calculates the adjoint of an expression *)
FSAdj[Susyno`LieGroups`conj[a_]] := SARAH`Tp[a];
FSAdj[SARAH`Tp[a_]] := a /; IsRealParameter[a];
FSAdj[SARAH`Tp[a_]] := Susyno`LieGroups`conj[a];
FSAdj[a_] := a /; IsHermitian[a];
FSAdj[a_] := SARAH`Adj[a];
FSAdj[a__] := Sequence @@ (FSAdj /@ Reverse[{a}]);

TraceEquality[SARAH`trace[a__], SARAH`trace[b__]] :=
    Or @@ (({a} === #)& /@ NestList[RotateLeft, {b}, Length[{b}] - 1]);

(* if all parameters in the trace are real, the trace is real *)
IsRealExpression[SARAH`trace[a__]] :=
    True /; And @@ (IsRealParameter /@ FindAllParameters[{a}]);

IsRealExpression[SARAH`trace[a__]] :=
    TraceEquality[SARAH`trace[a] /. SARAH`Adj -> FSAdj, SARAH`trace[FSAdj[a]]];

IsRealExpression[terms_SARAH`MatMul] :=
    And @@ (IsRealExpression[#]& /@ (List @@ terms));

IsRealExpression[sum[index_, start_, stop_, expr_]] :=
    IsRealExpression[expr];

IsRealExpression[otherwise_] := False;

HasPhase[particle_] :=
    MemberQ[First /@ SARAH`ParticlePhases, particle];

GetPhase[particle_ /; HasPhase[particle]] :=
    Cases[SARAH`ParticlePhases, {particle, phase_} :> phase][[1]];

GetPhase[_] := Null;

GetTypeFromDimension[sym_, {}] :=
    If[IsRealParameter[sym],
       CConversion`ScalarType[CConversion`realScalarCType],
       CConversion`ScalarType[CConversion`complexScalarCType]
      ];

GetTypeFromDimension[sym_, {0|1}] :=
    GetTypeFromDimension[sym, {}];

GetTypeFromDimension[sym_, {num_?NumberQ}] :=
    Module[{scalarType},
           scalarType = If[IsRealParameter[sym],
                           CConversion`realScalarCType,
                           CConversion`complexScalarCType
                          ];
           CConversion`VectorType[scalarType, num]
          ];

GetTypeFromDimension[sym_, {num1_?NumberQ, num2_?NumberQ}] :=
    If[IsRealParameter[sym],
       CConversion`MatrixType[CConversion`realScalarCType, num1, num2],
       CConversion`MatrixType[CConversion`complexScalarCType, num1, num2]
      ];

GetTypeFromDimension[sym_, {dims__} /; Length[{dims}] > 2 && (And @@ (NumberQ /@ {dims}))] :=
    If[IsRealParameter[sym],
       CConversion`TensorType[CConversion`realScalarCType, dims],
       CConversion`TensorType[CConversion`complexScalarCType, dims]
      ];

GetIntegerTypeFromDimension[{}] :=
    CConversion`ScalarType[CConversion`integerScalarCType];

GetIntegerTypeFromDimension[{0}] :=
    GetIntegerTypeFromDimension[{}];

GetIntegerTypeFromDimension[{1}] :=
    GetIntegerTypeFromDimension[{}];

GetIntegerTypeFromDimension[{num_?NumberQ}] :=
    CConversion`VectorType[CConversion`integerScalarCType, num];

GetIntegerTypeFromDimension[{num1_?NumberQ, num2_?NumberQ}] :=
    CConversion`MatrixType[CConversion`integerScalarCType, num1, num2];

GetIntegerTypeFromDimension[{dims__} /; Length[{dims}] > 2 && (And @@ (NumberQ /@ {dims}))] :=
    CConversion`TensorType[CConversion`integerScalarCType, dims];

GetRealTypeFromDimension[{}] :=
    CConversion`ScalarType[CConversion`realScalarCType];

GetRealTypeFromDimension[{0}] :=
    GetRealTypeFromDimension[{}];

GetRealTypeFromDimension[{1}] :=
    GetRealTypeFromDimension[{}];

GetRealTypeFromDimension[{num_?NumberQ}] :=
    CConversion`VectorType[CConversion`realScalarCType, num];

GetRealTypeFromDimension[{num1_?NumberQ, num2_?NumberQ}] :=
    CConversion`MatrixType[CConversion`realScalarCType, num1, num2];

GetRealTypeFromDimension[{dims__} /; Length[{dims}] > 2 && (And @@ (NumberQ /@ {dims}))] :=
    CConversion`TensorType[CConversion`realScalarCType, dims];

GetComplexTypeFromDimension[{}] :=
    CConversion`ScalarType[CConversion`complexScalarCType];

GetComplexTypeFromDimension[{0}] :=
    GetComplexTypeFromDimension[{}];

GetComplexTypeFromDimension[{1}] :=
    GetComplexTypeFromDimension[{}];

GetComplexTypeFromDimension[{num_?NumberQ}] :=
    CConversion`VectorType[CConversion`complexScalarCType, num];

GetComplexTypeFromDimension[{num1_?NumberQ, num2_?NumberQ}] :=
    CConversion`MatrixType[CConversion`complexScalarCType, num1, num2];

GetComplexTypeFromDimension[{dims__} /; Length[{dims}] > 2 && (And @@ (NumberQ /@ {dims}))] :=
    CConversion`TensorType[CConversion`complexScalarCType, dims];

GetType[sym_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    CConversion`GetScalarElementType[GetType[sym]];

GetType[FlexibleSUSY`SCALE] := GetRealTypeFromDimension[{}];

GetType[FlexibleSUSY`M[sym_]] :=
    GetRealTypeFromDimension[{SARAH`getGen[sym, FlexibleSUSY`FSEigenstates]}];

GetType[sym_?IsInputParameter] :=
    Cases[GetInputParametersAndTypes[], {sym, _, type_} :> type][[1]];

GetType[sym_?IsExtraParameter] :=
    Cases[GetExtraParametersAndTypes[], {sym, _, type_} :> type][[1]];

GetType[sym_] :=
    GetTypeFromDimension[sym, SARAH`getDimParameters[sym]];

GetParameterDimensions[sym_ /; (IsInputParameter[sym] || IsExtraParameter[sym])] :=
    Module[{type},
           type = GetType[sym];
           Switch[type,
                  CConversion`ScalarType[_], {1},
                  CConversion`VectorType[_, n_], {type[[2]]},
                  CConversion`ArrayType[_, n_], {type[[2]]},
                  CConversion`MatrixType[_, m_, n_], {type[[2]], type[[3]]},
                  CConversion`TensorType[_, indices__], List @@ Rest[type],
                  _, Print["Error: unknown parameter type: ", ToString[type]]; Quit[1];
                 ]
          ];

GetParameterDimensions[sym_] :=
    Module[{dim},
           dim = SARAH`getDimParameters[sym];
           Switch[dim,
                  {}, {1},
                  {0}, {1},
                  _, dim
                 ]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`ScalarType[_]}] := {};

CreateIndexReplacementRule[{parameter_, CConversion`VectorType[_,_] | CConversion`ArrayType[_,_]}] :=
    Module[{i},
           RuleDelayed @@ Rule[parameter[i_], parameter[i-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`MatrixType[_,_,_]}] :=
    Module[{i,j},
           RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`TensorType[_,_,_,_]}] :=
    Module[{i,j,k},
           RuleDelayed @@ Rule[parameter[i_,j_,k_], parameter[i-1,j-1,k-1]]
          ];

CreateIndexReplacementRule[{parameter_, CConversion`TensorType[_,_,_,_,_]}] :=
    Module[{i,j,k,l},
           RuleDelayed @@ Rule[parameter[i_,j_,k_,l_], parameter[i-1,j-1,k-1,l-1]]
          ];

CreateIndexReplacementRule[parameter_ /; (IsInputParameter[parameter] || IsExtraParameter[parameter])] :=
    CreateIndexReplacementRule[{parameter, GetType[parameter]}];

CreateIndexReplacementRule[parameter_] :=
    Module[{i,j,k,l, dim, rule},
           dim = SARAH`getDimParameters[parameter];
           rule = {};
           Switch[Length[dim],
                  1, rule = RuleDelayed @@ Rule[parameter[i_], parameter[i-1]];,
                  2, rule = RuleDelayed @@ Rule[parameter[i_,j_], parameter[i-1,j-1]];,
                  3, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_], parameter[i-1,j-1,k-1]];,
                  4, rule = RuleDelayed @@ Rule[parameter[i_,j_,k_,l_], parameter[i-1,j-1,k-1,l-1]];
                 ];
           rule
          ];

CreateIndexReplacementRules[pars_List] :=
    Flatten[CreateIndexReplacementRule /@ pars];

GetGUTNormalization[coupling_Symbol] :=
    Module[{pos, norm},
           pos = Position[SARAH`Gauge, coupling];
           If[pos =!= {},
              norm = SARAH`GUTren[pos[[1,1]]];
              If[NumericQ[norm],
                 Return[norm];
                ];
             ];
           Return[1];
          ];

ApplyGUTNormalization[] :=
    Module[{i, rules = {}, coupling},
           For[i = 1, i <= Length[SARAH`Gauge], i++,
               If[NumericQ[SARAH`GUTren[i]],
                  coupling = SARAH`Gauge[[i,4]];
                  AppendTo[rules, Rule[coupling, coupling SARAH`GUTren[i]]];
                 ];
              ];
           Return[rules];
          ];

CreateParameterDefinition[par_] :=
    CConversion`CreateCType[GetType[par]] <> " " <> CConversion`ToValidCSymbolString[par] <> ";\n";

CreateParameterDefinition[{par_, type_}] :=
    CConversion`CreateCType[type] <> " " <> CConversion`ToValidCSymbolString[par] <> ";\n";

CreateParameterDefinition[{par_, block_, type_}] :=
    CConversion`CreateCType[type] <> " " <> CConversion`ToValidCSymbolString[par] <> ";\n";

CreateParameterDefinitionAndDefaultInitialize[par_] :=
    CConversion`CreateCType[GetType[par]] <> " " <> CConversion`ToValidCSymbolString[par] <> "{};\n";

CreateParameterDefinitionAndDefaultInitialize[{par_, type_}] :=
    CConversion`CreateCType[type] <> " " <> CConversion`ToValidCSymbolString[par] <>
    CConversion`CreateDefaultAggregateInitialization[type] <> ";\n";

CreateParameterDefinitionAndDefaultInitialize[{par_, block_, type_}] :=
    CConversion`CreateCType[type] <> " " <> CConversion`ToValidCSymbolString[par] <>
    CConversion`CreateDefaultAggregateInitialization[type] <> ";\n";

CreateSetAssignment[name_, startIndex_, parameterType_, struct_:"pars"] :=
    Block[{},
          Print["Error: CreateSetAssignment: unknown parameter type: ", ToString[parameterType]];
          Quit[1];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType], struct_:"pars"] :=
    Module[{ass = ""},
           ass = name <> " = " <> struct <> "(" <> ToString[startIndex] <> ");\n";
           Return[{ass, 1}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`ScalarType[CConversion`complexScalarCType], struct_:"pars"] :=
    Module[{ass = "", type},
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           ass = name <> " = " <> type <> "(" <> struct <> "(" <> ToString[startIndex] <>
                 "), " <> struct <> "(" <> ToString[startIndex + 1] <> "));\n";
           Return[{ass, 2}];
          ];

CreateSetAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`realScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0},
           For[i = 0, i < rows, i++; count++,
               ass = ass <> name <> "(" <> ToString[i] <> ") = " <> struct <> "(" <>
                     ToString[startIndex + count] <> ");\n";
              ];
           If[rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, rows}];
          ];

CreateSetAssignment[name_, startIndex_, (CConversion`VectorType | CConversion`ArrayType)[CConversion`complexScalarCType, rows_], struct_:"pars"] :=
    Module[{ass = "", i, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < rows, i++; count+=2,
               ass = ass <> name <> "(" <> ToString[i] <> ") = " <>
                     type <> "(" <> struct <>
                     "(" <> ToString[startIndex + count    ] <> "), " <> struct <>
                     "(" <> ToString[startIndex + count + 1] <> "));\n";
              ];
           If[2 * rows != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0},
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count++,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j]
                         <> ") = " <> struct <> "(" <> ToString[startIndex + count] <> ");\n";
                  ];
              ];
           If[rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_], struct_:"pars"] :=
    Module[{ass = "", i, j, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < rows, i++,
               For[j = 0, j < cols, j++; count+=2,
                   ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <>
                         ") = " <> type <> "(" <> struct <>
                         "(" <> ToString[startIndex + count    ] <> "), " <> struct <>
                         "(" <> ToString[startIndex + count + 1] <> "));\n";
                  ];
              ];
           If[2 * rows * cols != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[2 rows * cols] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`TensorType[CConversion`realScalarCType, dim1_, dim2_, dim3_]] :=
    Module[{ass = "", i, j, k, count = 0},
           For[i = 0, i < dim1, i++,
               For[j = 0, j < dim2, j++,
                   For[k = 0, k < dim3, k++; count++,
                       ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <>
                             "," <> ToString[k] <>
                             ") = pars(" <> ToString[startIndex + count] <> ");\n";
                      ];
                  ];
              ];
           If[dim1 * dim2 * dim3 != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[dim1 dim2 dim3] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

CreateSetAssignment[name_, startIndex_, CConversion`TensorType[CConversion`complexScalarCType, dim1_, dim2_, dim3_]] :=
    Module[{ass = "", i, j, k, count = 0, type},
           type = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           For[i = 0, i < dim1, i++,
               For[j = 0, j < dim2, j++,
                   For[k = 0, k < dim3, k++; count+=2,
                       ass = ass <> name <> "(" <> ToString[i] <> "," <> ToString[j] <> "," <>
                             ToString[k] <> ") = " <> type <> "(" <>
                             "pars(" <> ToString[startIndex + count    ] <> "), " <>
                             "pars(" <> ToString[startIndex + count + 1] <> "));\n";
                      ];
                  ];
              ];
           If[2 * dim1 * dim2 * dim3 != count,
              Print["Error: CreateSetAssignment: something is wrong with the indices: "
                    <> ToString[2 dim1 dim2 dim3] <> " != " <> ToString[count]];];
           Return[{ass, count}];
          ];

(* assigns the parameter (str) to the container element *)
DisplayAssignParSet[str_String, startIndex_, currIdx_, struct_String] :=
    struct <> "(" <> ToString[startIndex + currIdx - 1] <> ") = " <> str <> ";\n";

(* creates parameter name and calls SetterFunc[] (= DisplayAssignParSet[] by default) *)
DisplayAssignPar[(h:(Re|Im))[(par_String)[idx__]], startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[ToString[h] <> "(" <> par <> "(" <>
         Utils`StringJoinWithSeparator[{idx}, ","] <>
         "))", startIndex, currIdx, struct];

DisplayAssignPar[(par_String)[idx__], startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[par <> "(" <>
         Utils`StringJoinWithSeparator[{idx}, ","] <>
         ")", startIndex, currIdx, struct];

DisplayAssignPar[(h:(Re|Im))[par_String], startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[ToString[h] <> "(" <> par <> ")", startIndex, currIdx, struct];

DisplayAssignPar[par_, startIndex_, currIdx_, struct_String, SetterFunc_:DisplayAssignParSet] :=
    SetterFunc[CConversion`RValueToCFormString[par], startIndex, currIdx, struct];

CreateDisplayAssignment[name_, startIndex_, type_, struct_:"pars"] :=
    { StringJoin @
        MapIndexed[DisplayAssignPar[#1,startIndex,First[#2],struct]&,
                   DecomposeParameter[name, type]],
      Length[DecomposeParameter[name, type]]
    };

(* assigns the string (str) to the std container element *)
DisplayAssignNameSet[str_String, startIndex_, currIdx_, struct_String] :=
    struct <> "[" <> ToString[startIndex + currIdx - 1] <> "] = \"" <> str <> "\";\n";

(* creates parameter name and calls DisplayAssignNameSet[] *)
DisplayAssignName[par_, startIndex_, currIdx_, struct_String] :=
    DisplayAssignPar[par, startIndex, currIdx, struct, DisplayAssignNameSet];

CreateStdVectorNamesAssignment[name_, startIndex_, type_, struct_:"names"] :=
    { StringJoin @
        MapIndexed[DisplayAssignName[#1,startIndex,First[#2],struct]&,
                   DecomposeParameter[name, type]],
      Length[DecomposeParameter[name, type]]
    };

CreateParameterSARAHNameStr[par_] :=
    "\"" <> CConversion`RValueToCFormString[par] <> "\"";

CreateParameterSARAHNames[name_, type_] :=
    Utils`StringJoinWithSeparator[CreateParameterSARAHNameStr /@ DecomposeParameter[name, type], ", "];

(* decomposes a parameter into its real pieces *)
DecomposeParameter[par_, parameterType_] :=
    Block[{},
          Print["Error: DecomposeParameter: unknown parameter type: ",
                ToString[parameterType]];
          Quit[1];
          ];

DecomposeParameter[name_, CConversion`ScalarType[CConversion`realScalarCType | CConversion`integerScalarCType]] :=
    { name };

DecomposeParameter[name_, CConversion`ScalarType[CConversion`complexScalarCType]] :=
    { Re[name], Im[name] };

DecomposeParameter[name_, (CConversion`VectorType|CConversion`ArrayType)[CConversion`realScalarCType, rows_]] :=
    Array[name, rows, 0];

DecomposeParameter[name_, (CConversion`VectorType|CConversion`ArrayType)[CConversion`complexScalarCType, rows_]] :=
    Flatten[{Re[#], Im[#]}& /@ Array[name, rows, 0]];

DecomposeParameter[name_, CConversion`MatrixType[CConversion`realScalarCType, rows_, cols_]] :=
    Flatten[Array[name, {rows, cols}, 0]];

DecomposeParameter[name_, CConversion`MatrixType[CConversion`complexScalarCType, rows_, cols_]] :=
    Flatten[{Re[#], Im[#]}& /@ Flatten[Array[name, {rows,cols}, 0]]];

DecomposeParameter[name_, CConversion`TensorType[CConversion`realScalarCType, dims__]] :=
    Flatten[Array[name, {dims}, 0]];

DecomposeParameter[name_, CConversion`TensorType[CConversion`complexScalarCType, dims__]] :=
    Flatten[{Re[#], Im[#]}& /@ Flatten[Array[name, {dims}, 0]]];

CreateEnumName[(h:(Re|Im))[par_]] :=
    CreateEnumName[h] <> CreateEnumName[par];

CreateEnumName[par_[idx__]] :=
    CreateEnumName[par] <>
    Utils`StringJoinWithSeparator[{idx}, "_", CConversion`ToValidCSymbolString];

CreateEnumName[par_] :=
    CConversion`ToValidCSymbolString[par];

CreateParameterEnums[name_, type_] :=
    Utils`StringJoinWithSeparator[CreateEnumName /@ DecomposeParameter[name, type], ", "];

CreateExtraParameterEnum[extraParameters_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParameterEnums[#, GetType[#]]& /@ extraParameters, ", "];
           If[Length[extraParameters] > 0, result = result <> ", ";];
           "enum Extra_parameters : int { " <> result <> "NUMBER_OF_EXTRA_PARAMETERS };\n"
          ];

CreateExtraParameterNames[extraParameters_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParameterSARAHNames[#,GetType[#]]& /@ extraParameters, ", "];
           "const std::array<std::string, NUMBER_OF_EXTRA_PARAMETERS> extra_parameter_names = {" <>
           result <> "};\n"
          ];

CreateInputParameterEnum[inputParameters_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParameterEnums[#[[1]],#[[3]]]& /@ inputParameters, ", "];
           If[Length[inputParameters] > 0, result = result <> ", ";];
           "enum Input_parameters : int { " <> result <> "NUMBER_OF_INPUT_PARAMETERS };\n"
          ];

CreateInputParameterNames[inputParameters_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParameterSARAHNames[#[[1]],#[[3]]]& /@ inputParameters, ", "];
           "const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> input_parameter_names = {" <>
           result <> "};\n"
          ];

SetInputParameter[parameter_, value_, wrapper_String, castToType_:None] :=
    Module[{parameterStr, valueStr},
           If[IsInputParameter[parameter],
              parameterStr = CConversion`ToValidCSymbolString[parameter];
              valueStr = CConversion`RValueToCFormString[value];
              If[wrapper == "",
                 parameterStr <> " = " <> CConversion`CastTo[valueStr,castToType] <> ";\n",
                 wrapper <> "(" <> parameterStr <> ") = " <> CConversion`CastTo[valueStr,castToType] <> ";\n"
                ]
              ,
              ""
             ]
          ];

SetPhase[phase_, value_, classPrefix_String] :=
    Module[{phaseStr, valueStr},
           If[IsPhase[phase],
              phaseStr = Phases`CreatePhaseName[phase];
              valueStr = CConversion`RValueToCFormString[value];
              classPrefix <> "set_" <> phaseStr <> "(" <> valueStr <> ");\n",
              ""
             ]
          ];

ConcatIndices[indices_List] :=
    Utils`StringJoinWithSeparator[ToString /@ indices,","];

CreateIndices[parameter_[indices__] /; And @@ (IsIndex /@ {indices})] :=
    "(" <> ConcatIndices[{indices}] <> ")";

CreateIndices[parameter_] := "";

AppendIfNotEmpty[str_String, sym_String] := If[str == "", "", str <> sym];

SetParameter[Re[parameter_], value_String, class_String, castToType_:None] :=
    Module[{parameterStr, indicesStr = "", newValue},
           If[GetIndices[parameter] =!= {},
              indicesStr = "(" <> ConcatIndices[GetIndices[parameter]] <> ")";
             ];
           parameterStr = CConversion`RValueToCFormString[StripIndices[parameter]];
           newValue = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
               "(" <> value <> ",Im(MODELPARAMETER(" <> parameterStr <> ")" <> indicesStr <> "))";
           SetParameter[parameter, newValue, class, castToType]
          ];

SetParameter[Im[parameter_], value_String, class_String, castToType_:None] :=
    Module[{parameterStr, indicesStr = "", newValue},
           If[GetIndices[parameter] =!= {},
              indicesStr = "(" <> ConcatIndices[GetIndices[parameter]] <> ")";
             ];
           parameterStr = CConversion`RValueToCFormString[StripIndices[parameter]];
           newValue = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <>
               "(Re(MODELPARAMETER(" <> parameterStr <> ")" <> indicesStr <> ")," <> value <> ")";
           SetParameter[parameter, newValue, class, castToType]
          ];

SetParameter[parameter_, value_String, class_String, castToType_:None] :=
    Module[{parameterStr, targetType = castToType},
           If[IsModelParameter[parameter] || IsExtraParameter[parameter],
              parameterStr = CConversion`ToValidCSymbolString[StripIndices[parameter]];
              (* if the parameter indices, we need to cast to the element type *)
              If[GetIndices[parameter] =!= {} && targetType =!= None,
                 targetType = CConversion`GetScalarElementType[targetType];
                ];
              class <> "set_" <> parameterStr <> "(" <>
              AppendIfNotEmpty[ConcatIndices[GetIndices[parameter]],","] <>
              CConversion`CastTo[value,targetType] <> ");\n"
              ,
              ""
             ]
          ];

SetParameter[parameter_, value_, class_String] :=
    SetParameter[parameter, CConversion`RValueToCFormString[value], class, GetType[parameter]]

SetParameter[Re[parameter_], value_, castToType_:CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> ".real(" <>
    CConversion`CastTo[CConversion`RValueToCFormString[value] <> CreateIndices[parameter],
                       castToType] <> ");\n";

SetParameter[Im[parameter_], value_, castToType_:CConversion`ScalarType[CConversion`realScalarCType]] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> ".imag(" <>
    CConversion`CastTo[CConversion`RValueToCFormString[value] <> CreateIndices[parameter],
                       castToType] <> ");\n";

SetParameter[parameter_, value_, castToType_:None] :=
    CConversion`ToValidCSymbolString[StripIndices[parameter]] <> CreateIndices[parameter] <> " = " <>
    CConversion`CastTo[CConversion`RValueToCFormString[value],castToType] <> ";\n";

SaveParameterLocally[parameters_List] :=
    StringJoin[SaveParameterLocally /@ parameters];

SaveParameterLocally[parameter_] :=
    Module[{ par, parStr, parStrSym },
           par = parameter /. { Re -> Identity, Im -> Identity, Abs -> Identity };
           parStr = CConversion`RValueToCFormString[par];
           parStrSym = CConversion`ToValidCSymbolString[par];
           "const auto save_" <> parStrSym <> "_raii = make_raii_save(" <> parStr <> ");\n"
          ];

RemoveProtectedHeads[expr_] :=
    expr /. { FlexibleSUSY`LowEnergyConstant[__] -> FlexibleSUSY`LowEnergyConstant[],
              FlexibleSUSY`Pole[__]  -> FlexibleSUSY`Pole[],
              FlexibleSUSY`BETA[__] -> FlexibleSUSY`BETA[],
              SARAH`Mass  -> FlexibleSUSY`M,
              SARAH`Mass2 -> FlexibleSUSY`M };

CreateRulesForProtectedHead[expr_, protectedHead_Symbol] :=
    Cases[expr, protectedHead[p__] :> Rule[protectedHead[p],Symbol["x$" <> ToString[Hash[p]]]], {0, Infinity}];

CreateRulesForProtectedHead[expr_, protectedHeads_List] :=
    Flatten @ Join[CreateRulesForProtectedHead[expr,#]& /@ protectedHeads];

FindMacro[par_] :=
    Which[IsModelParameter[par] , Global`MODELPARAMETER,
          IsOutputParameter[par], Global`MODELPARAMETER,
          IsPhase[par]          , Global`MODELPARAMETER,
          IsInputParameter[par] , Global`INPUTPARAMETER,
          IsExtraParameter[par] , Global`EXTRAPARAMETER,
          True                  , Identity
         ];

WrapPreprocessorMacroAround[expr_String, ___] := expr;

WrapPreprocessorMacroAround[expr_, protectedHeads_List:{FlexibleSUSY`Pole, SARAH`SM}] :=
    Module[{allPars, replacements, protectionRules, exprWithoutProtectedSymbols},
           allPars = Flatten[{FSModelParameters, FSInputParameters,
                              FSOutputParameters, FSPhysicalOutputParameters,
                              FSPhases, FSDerivedParameters, FSExtraParameters} /. FindAllParametersClassified[expr]];
           replacements = Join[
               RuleDelayed[#     , FindMacro[#][#]   ]& /@ allPars,
               RuleDelayed[#[i__], FindMacro[#][#][i]]& /@ allPars,
               {RuleDelayed[FlexibleSUSY`M[p_[i__]], FindMacro[FlexibleSUSY`M[p]][FlexibleSUSY`M[p]][i]]}
           ];
           protectionRules = CreateRulesForProtectedHead[expr, protectedHeads];
           exprWithoutProtectedSymbols = expr /. protectionRules;
           (* substitute back protected symbols *)
           exprWithoutProtectedSymbols /. replacements /. (Reverse /@ protectionRules)
          ];

DefineLocalConstCopy[parameter_, macro_String, prefix_String:""] :=
    "const auto " <> prefix <> CConversion`ToValidCSymbolString[parameter] <> " = " <>
    macro <> "(" <> CConversion`ToValidCSymbolString[parameter] <> ");\n";

PrivateCallLoopMassFunction[FlexibleSUSY`M[particle_Symbol]] :=
    "calculate_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]] <> "_pole();\n";

CalculateLocalPoleMasses[parameter_] :=
    "MODEL->" <> PrivateCallLoopMassFunction[parameter];

CreateLocalConstRefs[expr_] :=
    Module[{result = "", pars, inputSymbols, modelPars, outputPars,
            poleMasses, phases, depNum, extraPars},
           pars = FindAllParametersClassified[expr];
           inputSymbols = FSInputParameters /. pars;
           modelPars    = FSModelParameters /. pars;
           outputPars   = FSOutputParameters /. pars;
           phases       = FSPhases /. pars;
           depNum       = FSDerivedParameters /. pars;
           poleMasses   = FSPhysicalOutputParameters /. pars;
           extraPars    = FSExtraParameters /. pars;
           (result = result <> DefineLocalConstCopy[#,"INPUTPARAMETER"])& /@ inputSymbols;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ modelPars;
           (result = result <> DefineLocalConstCopy[#,"MODELPARAMETER"])& /@ outputPars;
           (result = result <> DefineLocalConstCopy[#,"PHASE"         ])& /@ phases;
           (result = result <> DefineLocalConstCopy[#,"DERIVEDPARAMETER"])& /@ depNum;
           (result = result <> DefineLocalConstCopy[#,"EXTRAPARAMETER"])& /@ extraPars;
           (result = result <> CalculateLocalPoleMasses[#])& /@ poleMasses;
           Return[result];
          ];

CreateLocalConstRefsForPhysicalParameters[expr_] :=
    Module[{result = "", pars, outputPars},
           pars = FindAllParametersClassified[expr];
           outputPars = FSOutputParameters /. pars;
           (result = result <> DefineLocalConstCopy[#,"PHYSICAL"])& /@ outputPars;
           Return[result];
          ];

DefineLocalConstCopyForBeta[{par_, -1}] :=
    DefineLocalConstCopy[par, "BETAPARAMETER", "beta_"];

DefineLocalConstCopyForBeta[{par_, loops_}] :=
    Module[{lstr = ToString[loops], parStr = CConversion`ToValidCSymbolString[par]},
           "const auto BETA1(" <> lstr <> "," <> parStr <>
           ") = BETAPARAMETER1(" <> lstr <> "," <> parStr <> ");\n"
          ];

CreateLocalConstRefsForBetas[expr_] :=
    Module[{result = "", exprFull, pars},
           (* add index to beta function for current loop level *)
           exprFull = expr /. FlexibleSUSY`BETA[p_] :> FlexibleSUSY`BETA[-1,p];
           pars = DeleteDuplicates @ \
                  Cases[exprFull, FlexibleSUSY`BETA[l_,p_] | FlexibleSUSY`BETA[l_,p_][___] :> {p,l}, {0, Infinity}];
           StringJoin[DefineLocalConstCopyForBeta /@ pars]
          ];

CreateLocalConstRefsForInputParameters[expr_, head_String:"INPUT"] :=
    Module[{result = "", pars, inputPars},
           pars = FindAllParametersClassified[expr];
           inputPars = FSInputParameters /. pars;
           (result = result <> DefineLocalConstCopy[#, head])& /@ inputPars;
           Return[result];
          ];

SetInputParameterTo[Re[par_], value_String] :=
    CConversion`ToValidCSymbolString[par] <> ".real(" <> value <> ");";

SetInputParameterTo[Im[par_], value_String] :=
    CConversion`ToValidCSymbolString[par] <> ".imag(" <> value <> ");";

SetInputParameterTo[par_, value_String] :=
    CConversion`ToValidCSymbolString[par] <> " = " <> value <> ";";

CreateCaseFromTuple[{key_?NumberQ, parameter_}] :=
    "case " <> ToString[key] <> ": input." <>
    SetInputParameterTo[parameter, "value"] <> " break;\n";

CreateCaseFromTuple[expr_] :=
    Block[{},
          Print["Error: not a valid {key,parameter} tuple: ", expr];
          ""
         ];

FillInputParametersFromTuples[minpar_List, blockName_String] :=
    Module[{result = ""},
           (result = result <> CreateCaseFromTuple[#])& /@ minpar;
           result = "switch (key) {\n" <> result <>
                    "default: WARNING(\"Unrecognized entry in block " <>
                    blockName <> ": \" << key); break;\n}\n";
           Return[result];
          ];

IncreaseIndex[ind_Integer, num_Integer] := ind + num;
IncreaseIndex[ind_, _]     := ind;
IncreaseIndices[a_[{ind__}], num_Integer] := a[IncreaseIndex[#,num]& /@ {ind}];
IncreaseIndices[a_[ind__], num_Integer] := a[Sequence @@ (IncreaseIndex[#,num]& /@ {ind})];
IncreaseIndices[a_, _]     := a;
IncreaseIndices[SARAH`Delta[a_, b_], num_Integer] :=
    CConversion`FSKroneckerDelta[IncreaseIndex[a,num], IncreaseIndex[b,num]];

IncreaseIndexLiterals[expr_] :=
    IncreaseIndexLiterals[expr, 1];

IncreaseIndexLiterals[expr_, num_Integer] :=
    IncreaseIndexLiterals[expr, num, Join[GetInputParameters[], GetExtraParameters[],
                                          allModelParameters, allOutputParameters]];

IncreaseIndexLiterals[expr_, num_Integer, heads_List] :=
    Module[{indexedSymbols, rules, allHeads},
           allHeads = Join[heads /. FlexibleSUSY`M -> Identity, {SARAH`Delta, SARAH`ThetaStep}];
           indexedSymbols = Extract[{expr}, Position[{expr}, s_[__] /; MemberQ[allHeads, s], Infinity]];
           rules = Rule[#, IncreaseIndices[#,num]] & /@ indexedSymbols;
           expr /. rules
          ];

DecreaseIndexLiterals[expr_] :=
    IncreaseIndexLiterals[expr, -1];

DecreaseIndexLiterals[expr_, heads_List] :=
    IncreaseIndexLiterals[expr, -1, heads];

DecreaseSumIndices[expr_] :=
    expr //. SARAH`sum[idx_, start_, stop_, exp_] :> FlexibleSUSY`SUM[idx, start - 1, stop - 1, exp];

ReplaceThetaStep[expr_] := expr /; FreeQ[expr,ThetaStep];

ReplaceThetaStep[expr_] :=
    Factor[
        Expand[expr] //. {
            Times[a__, ThetaStep[i1_, i2_], b__] :> FlexibleSUSY`IF[i1 < i2, a b, 0],
            Times[a__, ThetaStep[i1_, i2_]] :> FlexibleSUSY`IF[i1 < i2, a, 0],
            Times[ThetaStep[i1_, i2_], a__] :> FlexibleSUSY`IF[i1 < i2, a, 0]
        }
    ];

(* Converts an expression to a valid C++ string. *)
ExpressionToString[expr_] :=
    CConversion`RValueToCFormString[
        Simplify[Parameters`DecreaseIndexLiterals[Parameters`DecreaseSumIndices[ReplaceThetaStep[expr]]]]];

ExpressionToString[expr_, heads_] :=
    CConversion`RValueToCFormString[
        Simplify[Parameters`DecreaseIndexLiterals[Parameters`DecreaseSumIndices[ReplaceThetaStep[expr]], heads]]];

GetEffectiveMu[] :=
    Module[{},
           If[!ValueQ[FlexibleSUSY`EffectiveMu],
              Print["Error: effective Mu parameter not defined!"];
              Print["   Please set EffectiveMu to the expression of the"];
              Print["   effective Mu parameter."];
              Quit[1];
             ];
           FlexibleSUSY`EffectiveMu
          ];

GetEffectiveMASqr[] :=
    Module[{},
           If[!ValueQ[FlexibleSUSY`EffectiveMASqr],
              Print["Error: effective CP-odd Higgs mass not defined!"];
              Print["   Please set EffectiveMASqr to the expression of the"];
              Print["   effective CP-odd Higgs mass."];
              Quit[1];
             ];
           FlexibleSUSY`EffectiveMASqr
          ];

GetParameterFromDescription[description_String] :=
    Module[{parameter},
           parameter =Cases[SARAH`ParameterDefinitions,
                            {parameter_,
                             {___, SARAH`Description -> description, ___}} :>
                            parameter];
           If[Length[parameter] == 0,
              Print["Error: Parameter with description \"", description,
                    "\" not found."];
              Return[Null];
             ];
           If[Length[parameter] > 1,
              Print["Warning: Parameter with description \"", description,
                    "\" not unique."];
             ];
           parameter[[1]]
          ];

GetParticleFromDescription[description_String, eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{particle},
           particle =Cases[SARAH`ParticleDefinitions[eigenstates],
                            {particle_,
                             {___, SARAH`Description -> description, ___}} :>
                            particle];
           If[Length[particle] == 0,
              DebugPrint["Note: Particle with description \"", description,
                         "\" not found."];
              Return[Null];
             ];
           If[Length[particle] > 1,
              Print["Warning: Particle with description \"", description,
                    "\" not unique."];
             ];
           particle[[1]]
          ];

GetParticleFromDescription[multipletName_String, splitNames_List] :=
    Module[{result},
           result = GetParticleFromDescription[multipletName];
           If[result =!= Null, Return[{result}]];
           DeleteCases[GetParticleFromDescription /@ splitNames, Null]
          ];

GetPDGCodesForParticle[SARAH`bar[particle_]] :=
    -GetPDGCodesForParticle[particle];

GetPDGCodesForParticle[Susyno`LieGroups`conj[particle_]] :=
    -GetPDGCodesForParticle[particle];

GetPDGCodesForParticle[particle_] :=
    Module[{pdgList},
            pdgList = SARAH`getPDGList[particle];
            If[pdgList === None,
               pdgList = {};
              ];
           pdgList
          ];

NumberOfIndependentEntriesOfSymmetricMatrix[n_] := (n^2 + n) / 2;

AppendGenerationIndices[expr_List] :=
    AppendGenerationIndices /@ expr;

AppendGenerationIndices[expr_Symbol] :=
    Switch[SARAH`getDimParameters[expr],
           {}                          , expr,
           {0}                         , expr,
           {1}                         , expr,
           {idx_}                      , expr[SARAH`gt1],
           {idx1_, idx2_}              , expr[SARAH`gt1, SARAH`gt2],
           {idx1_, idx2_, idx3_}       , expr[SARAH`gt1, SARAH`gt2, SARAH`gt3],
           {idx1_, idx2_, idx3_, idx4_}, expr[SARAH`gt1, SARAH`gt2, SARAH`gt3, SARAH`gt4]
          ];

AppendGenerationIndices[expr_] := expr;

(*
 * Expands a list of expressions of the form
 *
 *   { 1 + A[SARAH`gt1] }
 *
 * to
 *
 *   { 1 + A[1], 1 + A[2], 1 + A[3] }
 *
 * where the indices SARAH`gt1 and SARAH`gt2 are assumed to run from 1
 * to their maximum value.
 *)
ExpandExpressions[eqs_List] :=
    Join[Flatten[ExpandExpressions /@ eqs]];

GetIdxRange[{idx_,par_}] :=
    Module[{dim = GetParameterDimensions[par]},
           Switch[dim,
                  {}  , {idx,1,1},
                  {0} , {idx,1,1},
                  {n_?NumberQ}, {idx,1,dim[[1]]},
                  {__}, {idx,1,1}
                 ]
          ];

ExpandExpressions[eq_] :=
    Module[{par, indexSymbols, indexRanges, indices = {SARAH`gt1, SARAH`gt2}},
           indexSymbols = DeleteDuplicates[
               Join @@ (Cases[eq, par_[idx_ /; !FreeQ[idx,#]] /;
                                  MemberQ[Join[GetModelParameters[], GetOutputParameters[]], par] :> {#,par},
                              {0,Infinity}]& /@ indices)];
           indexRanges  = DeleteDuplicates[GetIdxRange /@ indexSymbols];
           If[indexRanges === {},
              eq,
              Table[eq, Evaluate[Sequence @@ indexRanges]]
             ]
          ];

(* removes indices from model Parameter, taking SARAH's L, B, T, Q
   into account.  *)
StripIndices[par_[idx___] /; And @@ (NumberQ /@ {idx})] := par;

StripIndices[par_[idx___] /; MemberQ[Join[allModelParameters,allOutputParameters],par]] := par;

StripIndices[par_] := par;

StripIndicesRules[indices_List, numberOfIndices_] :=
    Flatten[{ RuleDelayed[a_[Sequence @@ #], a]& /@ Permutations[indices, {numberOfIndices}] }];

StripSARAHIndicesRules[numberOfIndices_] :=
    StripIndicesRules[sarahIndices, numberOfIndices];

ExtractParametersFromSARAHBetaLists[beta_List] :=
    StripIndices[First[#]]& /@ beta;

ExtractParametersFromSARAHBetaLists[_] := {};

GetModelParametersWithMassDimension[dim_?IntegerQ] :=
    Module[{dimPars},
           Switch[dim,
                  0, dimPars = Join[SARAH`BetaGauge, SARAH`BetaLijkl, SARAH`BetaYijk, SARAH`BetaQijkl];,
                  1, dimPars = Join[SARAH`BetaMuij, SARAH`BetaTijk, SARAH`BetaMi, SARAH`BetaDGi, SARAH`BetaVEV];,
                  2, dimPars = Join[SARAH`BetaLi, SARAH`BetaBij, SARAH`Betam2ij];,
                  3, dimPars = Join[SARAH`BetaLSi];,
                  _,
                  Print["Error: GetModelParametersWithMassDimension: ", dim,
                        " is not a valid dimension"];
                  Return[{}];
                 ];
           ExtractParametersFromSARAHBetaLists[dimPars]
          ];

GetParametersWithMassDimension[dim_?IntegerQ] :=
    Module[{dimPars},
           dimPars = GetModelParametersWithMassDimension[dim];
           dimPars = Join[dimPars, Flatten[Cases[extraMassDimensions, {dim, pars_} :> pars]]];
           DeleteDuplicates[dimPars]
          ];

GetModelParameterMassDimension[par_?IsModelParameter] :=
    Module[{i, parsList},
           For[i = 0, i <= 3, i++,
               parsList = GetModelParametersWithMassDimension[i];
               If[MemberQ[parsList, StripIndices[par]],
                  Return[i]
                 ];
              ];
           Print["Error: mass dimension for ", par, " not found!"];
           Quit[1];
          ];

GetModelParameterMassDimension[par_] :=
    Block[{},
          Print["Error: GetModelParameterMassDimension:", par, " is not a model parameter!"];
          Quit[1];
         ];

IsGaugeCoupling[par_] :=
    MemberQ[ExtractParametersFromSARAHBetaLists[SARAH`BetaGauge], par];

IsYukawaCoupling[par_] :=
    MemberQ[ExtractParametersFromSARAHBetaLists[SARAH`BetaYijk], par];

IsVEV[par_] :=
    MemberQ[ExtractParametersFromSARAHBetaLists[SARAH`BetaVEV], par];

AreLinearDependent[{eq1_, eq2_}, parameters_List] :=
    Module[{frac = Simplify[eq1/eq2 /. FlexibleSUSY`tadpole[_] -> 0],
            pars},
           (* ignore parameter heads Re[], Im[], Abs[], Phase[] *)
           pars = parameters /. { Re[p_] :> p, Im[p_] :> p,
                                  Abs[p_] :> p, FlexibleSUSY`Phase[p_] :> p };
           And @@ (FreeQ[frac,#]& /@ pars)
          ];

GetThirdGeneration[par_] :=
    Which[IsScalar[par], par,
          IsMatrix[par], par[2,2],
          True, Print["Warning: GetThirdGeneration[",par,"]: unknown type"]; par
         ];

GetSARAHParameters[] :=
    First /@ SARAH`SARAHparameters;

GetAllDependenceSPhenoSymbols[] :=
    DeleteDuplicates @ Flatten @
    Cases[SARAH`ParameterDefinitions,
          {parameter_, {___, SARAH`DependenceSPheno -> value:Except[None], ___}} :> parameter];

GetAllDependenceSPhenoRules[] :=
    Cases[SARAH`ParameterDefinitions,
          {parameter_, {___, SARAH`DependenceSPheno -> value:Except[None], ___}} :> RuleDelayed[parameter, value]];

GetDependenceSPhenoSymbols[] :=
    Module[{sarahPars = GetSARAHParameters[]},
           Select[GetAllDependenceSPhenoSymbols[], MemberQ[sarahPars,#]&]
          ];

GetDependenceSPhenoRules[] :=
    Module[{sarahPars = GetSARAHParameters[]},
           Select[GetAllDependenceSPhenoRules[], MemberQ[sarahPars,First[#]]&]
          ];

GetAllOutputParameterDependencies[expr_] :=
    Complement[Select[Join[GetSARAHParameters[],
                           GetDependenceSPhenoSymbols[]],
                      (!FreeQ[expr,#])&],
               GetModelParameters[]];

GetAllOutputParameterDependenciesReplaced[expr_] :=
    DeleteCases[GetAllOutputParameterDependencies[expr /. GetDependenceSPhenoRules[]], _?NumericQ];

GetOutputParameterDependencies[expr_] :=
    Select[GetOutputParameters[],
           (!FreeQ[GetAllOutputParameterDependenciesReplaced[expr],#])&];

GetExponent[a_^b_] := -I b;
GetExponent[a_]    := a;

GetIntermediateOutputParameterDependencies[expr_] :=
    Complement[
        GetAllOutputParameterDependenciesReplaced[expr],
        Join[GetOutputParameters[], GetInputParameters[], GetExponent /@ GetPhases[]]
    ];

CreateExtraParameterArrayGetter[{}] :=
    "return Eigen::ArrayXd();\n";

CreateExtraParameterArrayGetter[extraParameters_List] :=
    Module[{get = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[extraParameters], i++,
               par  = extraParameters[[i]];
               type = GetType[extraParameters[[i]]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type];
               get = get <> assignment;
               paramCount += nAssignments;
              ];
           get = "Eigen::ArrayXd pars(" <> ToString[paramCount] <> ");\n\n" <>
                 get <> "\n" <>
                 "return pars;";
           Return[get];
          ];

CreateExtraParameterArraySetter[extraParameters_List] :=
    Module[{set = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[extraParameters], i++,
               par  = extraParameters[[i]];
               type = GetType[extraParameters[[i]]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           Return[set];
          ];

CreateInputParameterArrayGetter[{}] :=
    "return Eigen::ArrayXd();\n";

CreateInputParameterArrayGetter[inputParameters_List] :=
    Module[{get = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type];
               get = get <> assignment;
               paramCount += nAssignments;
              ];
           get = "Eigen::ArrayXd pars(" <> ToString[paramCount] <> ");\n\n" <>
                 get <> "\n" <>
                 "return pars;";
           Return[get];
          ];

CreateInputParameterArraySetter[inputParameters_List] :=
    Module[{set = "", paramCount = 0, name = "", par,
            type, i, assignment = "", nAssignments = 0},
           For[i = 1, i <= Length[inputParameters], i++,
               par  = inputParameters[[i,1]];
               type = inputParameters[[i,2]];
               name = CConversion`ToValidCSymbolString[par];
               {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
               set = set <> assignment;
               paramCount += nAssignments;
              ];
           Return[set];
          ];

CreateModelParameterGetter[par_] :=
    Module[{name = CConversion`ToValidCSymbolString[par]},
           CConversion`CreateInlineGetters[name, name, GetType[par]]
          ];

CreateDelegateModelParameterGetter[par_, macro_String:"SUPER"] :=
    Module[{name = CConversion`ToValidCSymbolString[par]},
           CConversion`CreateInlineGetters[name, name, GetType[par], "", macro]
          ];

CreateModelParameterSetter[par_] :=
    Module[{name = CConversion`ToValidCSymbolString[par], type = GetType[par]},
           CConversion`CreateInlineSetters[name, type]
          ];

FindSLHABlock[blockList_List, par_] :=
    Module[{foundBlocks},
           foundBlocks = Cases[blockList, {par, block_, ___} :> block];
           If[foundBlocks === {},
              Print["Error: FindSLHABlock: no input block defined for parameter ", par];
              Quit[1];
             ];
           foundBlocks[[1]]
          ];

SetSMParameter[FlexibleSUSY`AlphaEMInvInput    , value_String, struct_String] := struct <> ".setAlphaEmInput(1./(" <> value <> "))";
SetSMParameter[FlexibleSUSY`GFermiInput        , value_String, struct_String] := struct <> ".setFermiConstant(" <> value <> ")";
SetSMParameter[FlexibleSUSY`AlphaSInput        , value_String, struct_String] := struct <> ".setAlphaSInput(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MZPoleInput        , value_String, struct_String] := struct <> ".setPoleMZ(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MBottomMbottomInput, value_String, struct_String] := struct <> ".setMbMb(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MTopPoleInput      , value_String, struct_String] := struct <> ".setPoleMt(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MTauPoleInput      , value_String, struct_String] := struct <> ".setPoleMtau(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MNeutrino3PoleInput, value_String, struct_String] := struct <> ".setNeutrinoPoleMass(3, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MWPoleInput        , value_String, struct_String] := struct <> ".setPoleMW(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MElectronPoleInput , value_String, struct_String] := struct <> ".setPoleMel(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MNeutrino1PoleInput, value_String, struct_String] := struct <> ".setNeutrinoPoleMass(1, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MMuonPoleInput     , value_String, struct_String] := struct <> ".setPoleMmuon(" <> value <> ")";
SetSMParameter[FlexibleSUSY`MNeutrino2PoleInput, value_String, struct_String] := struct <> ".setNeutrinoPoleMass(2, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MDown2GeVInput     , value_String, struct_String] := struct <> ".setMass(softsusy::mDown, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MUp2GeVInput       , value_String, struct_String] := struct <> ".setMass(softsusy::mUp, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MStrange2GeVInput  , value_String, struct_String] := struct <> ".setMass(softsusy::mStrange, " <> value <> ")";
SetSMParameter[FlexibleSUSY`MCharmMCharm       , value_String, struct_String] := struct <> ".setMass(softsusy::mCharm, " <> value <> ")";

End[];

EndPackage[];
