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

BeginPackage["Constraint`", {"CConversion`", "BetaFunction`", "Parameters`", "TextFormatting`", "TreeMasses`", "Utils`"}];

ApplyConstraints::usage="";
CalculateScale::usage="";
DefineInputParameters::usage="";
InitialGuessAtLowScaleGaugeCouplings::usage="";
IsFixed::usage="returns true if given parameter is fixed in given constraint";

FindFixedParametersFromConstraint::usage="Returns a list of all
parameters which are fixed by the given constraint";

CheckConstraint::usage="Checks a given constraint for syntax errors.";

SanityCheck::usage="Checks general features of the constraints.";

SetBetaFunctions::usage=""

RestrictScale::usage="";

CheckPerturbativityForParameters::usage="";

GetSMMatchingScale::usage="returns SM matching scale from low-energy data set";

InitialApplyConstraint::usage="apply constaints before calculating the mass spectrum";

Begin["`Private`"];

allBetaFunctions = {};

SetBetaFunctions[pars_List] := allBetaFunctions = pars;

ApplyConstraint[{parameter_, value_}, modelPrefix_String] :=
    Which[Parameters`IsModelParameter[parameter],
          Parameters`SetParameter[parameter, value, modelPrefix],
          Parameters`IsInputParameter[parameter],
          Parameters`SetInputParameter[parameter, value, "INPUTPARAMETER"],
          Parameters`IsPhase[parameter],
          Parameters`SetPhase[parameter, value, modelPrefix],
          Parameters`IsExtraParameter[parameter],
          Parameters`SetParameter[parameter, value, modelPrefix],
          True,
          Print["Error: ", parameter, " cannot be set in the constraint,",
                " because it is neither a model nor an input parameter!"];
          Quit[1];
          ""
         ];

ApplyConstraint[{parameter_ | parameter_[__] /; parameter === SARAH`UpYukawa,
                 value_ /; (!FreeQ[value, Global`upQuarksDRbar] || value === Automatic)},
                modelPrefix_String] :=
    "calculate_" <> CConversion`ToValidCSymbolString[parameter] <> "_DRbar();\n";

ApplyConstraint[{parameter_ | parameter_[__] /; parameter === SARAH`DownYukawa,
                 value_ /; (!FreeQ[value, Global`downQuarksDRbar] || value === Automatic)},
                modelPrefix_String] :=
    "calculate_" <> CConversion`ToValidCSymbolString[parameter] <> "_DRbar();\n";

ApplyConstraint[{parameter_ | parameter_[__] /; parameter === SARAH`ElectronYukawa,
                 value_ /; (!FreeQ[value, Global`downLeptonsDRbar] || value === Automatic)},
                modelPrefix_String] :=
    "calculate_" <> CConversion`ToValidCSymbolString[parameter] <> "_DRbar();\n";

ApplyConstraint[{parameter_,
                 value_ /; !FreeQ[value, Global`neutrinoDRbar]}, modelPrefix_String] :=
    Parameters`SetParameter[parameter, value, modelPrefix];

ApplyConstraint[{parameter_ /; !MemberQ[{SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa}, parameter], value_ /; value === Automatic}, modelPrefix_String] :=
    Block[{},
          Print["Error: cannot determine ", parameter, " automatically!"];
          Quit[1];
         ];

CreateStartPoint[parameters_List, name_String] :=
    Module[{dim, dimStr, startPoint = "", i},
           dim = Length[parameters];
           dimStr = ToString[dim];
           For[i = 1, i <= dim, i++,
               startPoint = startPoint <> If[i==1,"",", "] <>
                            CConversion`RValueToCFormString[parameters[[i]]];
              ];
           Parameters`CreateLocalConstRefs[parameters] <> "\n" <>
           "Eigen::VectorXd " <> name <> "(" <> dimStr <> ");\n" <>
           name <> " << " <> startPoint <> ";\n"
          ];

SetModelParametersFromVector[model_String, vector_String, parameters_List] :=
    Module[{result = "", i, gslElement},
           For[i = 1, i <= Length[parameters], i++,
               gslElement = vector <> "(" <> ToString[i-1] <> ")";
               result = result <> Parameters`SetParameter[parameters[[i]],gslElement,model];
              ];
           Return[result];
          ];

SetVectorFromExpressions[vector_String, expressions_List] :=
    Module[{result = "", i, gslElement},
           For[i = 1, i <= Length[expressions], i++,
               gslElement = vector <> "(" <> ToString[i-1] <> ") = " <>
                            CConversion`RValueToCFormString[expressions[[i]]] <> ";\n";
               result = result <> gslElement;
              ];
           Return[result];
          ];

CreateMinimizationFunctionWrapper[functionName_String, dim_Integer, parameters_List, function_, modelPrefix_String] :=
    Module[{type, stype},
           type  = CConversion`CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, 1]];
           stype = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
"auto " <> functionName <> " = [&](const "<> type <> "& x) {
" <> TextFormatting`IndentText[SetModelParametersFromVector[modelPrefix,"x",parameters]] <> "
   " <> modelPrefix <> "calculate_DRbar_masses();
" <> TextFormatting`IndentText[Parameters`CreateLocalConstRefs[function]] <> "
   return " <> CConversion`RValueToCFormString[function] <> ";
};
"
          ];

localFunctionWrapper = 0;

CreateSolverName[] := "solver_" <> ToString[localFunctionWrapper++];

ApplyConstraint[FlexibleSUSY`FSMinimize[parameters_List, function_], modelPrefix_String] :=
    Module[{callMinimizer, dim, dimStr, startPoint, functionWrapper, functionName},
           dim = Length[parameters];
           dimStr = ToString[dim];
           startPoint = CreateStartPoint[parameters, "start_point"];
           functionName = CreateSolverName[];
           functionWrapper = CreateMinimizationFunctionWrapper[functionName,dim,parameters,
                                                               Parameters`DecreaseIndexLiterals[function],
                                                               modelPrefix];
           callMinimizer = functionWrapper <> "\n" <> startPoint <>
                           "Minimizer<" <> dimStr <>
                           "> minimizer(" <> functionName <> ", 100, 1.0e-2);\n" <>
                           "const int status = minimizer.minimize(start_point);\n" <>
                           "VERBOSE_MSG(\"\\tminimizer status: \" << gsl_strerror(status));\n";
           "\n{" <> TextFormatting`IndentText[callMinimizer] <> "}\n"
          ];

CreateRootFinderFunctionWrapper[functionName_String, dim_Integer, parameters_List, function_List, modelPrefix_String] :=
    Module[{type},
           type = CConversion`CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, 1]];
"auto " <> functionName <> " = [&](const "<> type <> "& x) {
" <> TextFormatting`IndentText[SetModelParametersFromVector[modelPrefix,"x",parameters]] <> "
   " <> modelPrefix <> "calculate_DRbar_masses();
" <> TextFormatting`IndentText[Parameters`CreateLocalConstRefs[function]] <> "
   "<> type <> " f;
" <> TextFormatting`IndentText[SetVectorFromExpressions["f",function]] <> "
   return f;
};
"
          ];

ApplyConstraint[FlexibleSUSY`FSFindRoot[parameters_List, function_List], modelPrefix_String] :=
    Module[{callRootFinder, dim, dimStr, startPoint, functionWrapper, functionName},
           dim = Length[parameters];
           dimStr = ToString[dim];
           startPoint = CreateStartPoint[parameters, "start_point"];
           functionName = CreateSolverName[];
           functionWrapper = CreateRootFinderFunctionWrapper[functionName,dim,parameters,
                                                             Parameters`DecreaseIndexLiterals[function],
                                                             modelPrefix];
           callRootFinder = functionWrapper <> "\n" <> startPoint <>
                           "Root_finder<" <> dimStr <>
                           "> root_finder(" <> functionName <> ", 100, 1.0e-2);\n" <>
                           "const int status = root_finder.find_root(start_point);\n" <>
                           "VERBOSE_MSG(\"\\troot finder status: \" << gsl_strerror(status));\n";
           "\n{\n" <> TextFormatting`IndentText[callRootFinder] <> "}\n\n"
          ];

ApplyConstraint[FlexibleSUSY`FSSolveEWSBFor[___], modelPrefix_String] :=
    modelPrefix <> "solve_ewsb();\n";

ExpandRestrictParameter[FlexibleSUSY`FSRestrictParameter[p_, interval_, FlexibleSUSY`FSNoProblem, expr___]] :=
    ExpandRestrictParameter[FlexibleSUSY`FSRestrictParameter[p, interval, 0, expr]];

ExpandRestrictParameter[FlexibleSUSY`FSRestrictParameter[p_, interval_, problem_]] :=
    ExpandRestrictParameter[FlexibleSUSY`FSRestrictParameter[p, interval, problem, p]];

ExpandRestrictParameter[FlexibleSUSY`FSRestrictParameter[p_, {istart_, iend_}, problem_, expr_]] :=
    {p,
     FlexibleSUSY`WHICH[
         p < istart, expr + problem,
         p > iend  , expr + problem,
         True      , p
     ]
    };

ExpandRestrictParameter[FlexibleSUSY`FSInitialSetting[p_, expr_]] := {p, expr};

ApplyConstraint[patt:(FlexibleSUSY`FSRestrictParameter | FlexibleSUSY`FSInitialSetting)[__], modelPrefix_String] :=
    ApplyConstraint[ExpandRestrictParameter[patt], modelPrefix];

ApplyConstraint[Null, _] :=
    Block[{},
          Print["Error: Null is not a valid constraint setting!"];
          Print["   Maybe there is a trailing comma in the constraint list?"];
          Quit[1];
         ];

ApplyConstraint[p_, _] :=
    Block[{},
          Print["Error: This is not a valid constraint setting: ", p];
          Quit[1];
         ];

ApplyConstraints[settings_List, modelPrefix_String:"MODEL->"] :=
    Module[{result, noMacros, noTemp},
           noMacros = DeleteCases[
               settings,
               FlexibleSUSY`FSMinimize[__] | \
               FlexibleSUSY`FSFindRoot[__] | \
               FlexibleSUSY`FSSolveEWSBFor[__] | \
               {FlexibleSUSY`FSTemporary[_], _} | \
               FlexibleSUSY`FSRestrictParameter[__] | \
               FlexibleSUSY`FSInitialSetting[__]
           ];
           noTemp = DeleteCases[
               settings,
               {FlexibleSUSY`FSTemporary[_], _} | \
               FlexibleSUSY`FSRestrictParameter[__] | \
               FlexibleSUSY`FSInitialSetting[__]
           ];
           result = Parameters`CreateLocalConstRefs[(#[[2]])& /@ noMacros];
           result = result <> AddBetas[noMacros, modelPrefix];
           result = result <> "\n";
           (result = result <> ApplyConstraint[#, modelPrefix])& /@ noTemp;
           Return[result];
          ];

ContainsBetas[expr_] := !FreeQ[expr, FlexibleSUSY`BETA];

(* -1 = current beta-function loop order *)
MaxBetaLoopOrder[expr_] :=
    Sort @ Cases[expr /. FlexibleSUSY`BETA[p_] :> FlexibleSUSY`BETA[-1,p],
                 FlexibleSUSY`BETA[l_, _] | FlexibleSUSY`BETA[l_, _][___] :> l, {0, Infinity}];

CallCalcBeta[-1, modelPrefix_String] :=
    "const " <> FlexibleSUSY`FSModelName <> "_soft_parameters beta_functions(" <> modelPrefix <> "calc_beta());\n";

CallCalcBeta[l_?IntegerQ, modelPrefix_String] :=
    Module[{lstr = ToString[l]},
           "const " <> FlexibleSUSY`FSModelName <>
           "_soft_parameters beta_functions_" <> lstr <>
           "L(" <> modelPrefix <> "calc_beta(" <> lstr <> "));\n"
          ];

AddBetas[expr_?ContainsBetas, modelPrefix_String] :=
    StringJoin[CallCalcBeta[#, modelPrefix]& /@ MaxBetaLoopOrder[expr]] <>
    Parameters`CreateLocalConstRefsForBetas[expr];

AddBetas[_, modelPrefix_String] := "";

FindFixedParametersFromSetting[{parameter_, value_}] := Parameters`StripIndices[parameter];
FindFixedParametersFromSetting[FlexibleSUSY`FSMinimize[parameters_List, value_]] := parameters;
FindFixedParametersFromSetting[FlexibleSUSY`FSFindRoot[parameters_List, value_]] := parameters;
FindFixedParametersFromSetting[FlexibleSUSY`FSSolveEWSBFor[parameters_List]] := parameters;
FindFixedParametersFromSetting[_] := {};

FindFixedParametersFromConstraint[settings_List] :=
    DeleteDuplicates[Flatten[FindFixedParametersFromSetting /@ settings]];

CheckSetting[patt:(FlexibleSUSY`FSMinimize|FlexibleSUSY`FSFindRoot)[parameters_, value_],
             constraintName_String, __] :=
    Module[{unknownParameters},
           If[Head[parameters] =!= List,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   First parameter must be a list!"];
              Return[False];
             ];
           If[parameters === {},
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   List of parameters is empty!"];
              Return[False];
             ];
           If[MatchQ[patt, FlexibleSUSY`FSFindRoot[___]],
              If[!MatchQ[patt, FlexibleSUSY`FSFindRoot[_List, _List]],
                 Print["Error: Syntax error in constraint ", constraintName, ": ", InputForm[patt]];
                 Print["   Correct syntax: FSFindRoot[{a,b}, {f[a],f[b]}]"];
                 Return[False];
                 ,
                 If[Length[parameters] =!= Length[value],
                    Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
                    Print["   Argument lists must have the same length!"];
                    Return[False];
                   ];
                ];
             ];
           If[MatchQ[patt, FlexibleSUSY`FSMinimize[_List, _List]],
              Print["Error: Syntax error in constraint ", constraintName, ": ", InputForm[patt]];
              Print["   Correct syntax: FSMinimize[{a,b}, f[a,b]]"];
              Return[False];
             ];
           unknownParameters = Select[parameters, (!Parameters`IsModelParameter[#] &&
                                                   !Parameters`IsExtraParameter[#])&];
           If[unknownParameters =!= {},
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   Unknown parameters: ", unknownParameters];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:FlexibleSUSY`FSSolveEWSBFor[parameters_List],
             constraintName_String, __] :=
    Module[{unknownParameters = Select[parameters, (!Parameters`IsModelParameter[#] &&
                                                    !Parameters`IsExtraParameter[#])&]},
           If[unknownParameters =!= {},
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt],
                    "   Unknown parameters: ", unknownParameters];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:{parameter_[idx_Integer], value_}, constraintName_String, __] :=
    Module[{modelParameters, dim},
           modelParameters = Parameters`GetModelParameters[];
           If[!CheckSetting[{parameter, value}],
              Return[False];
             ];
           dim = SARAH`getDimParameters[parameter];
           If[Length[dim] =!= 1,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " has ", Length[dim],
                    " dimensions, but one index is accessed!"];
              Return[False];
             ];
           dim = dim[[1]];
           If[1 > idx || idx > dim,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["    ", parameter, " index out of range!",
                    " Allowed range: 1 ... ", dim];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:{parameter_[idx1_Integer, idx2_Integer], value_}, constraintName_String, __] :=
    Module[{modelParameters, dim},
           modelParameters = Parameters`GetModelParameters[];
           If[!CheckSetting[{parameter, value}],
              Return[False];
             ];
           dim = SARAH`getDimParameters[parameter];
           If[Length[dim] =!= 2,
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " has ", Length[dim],
                    " dimensions, but two indices are accessed!"];
              Return[False];
             ];
           If[1 > idx1 || idx1 > dim[[1]],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " first index out of range! ",
                    "Allowed range: 1 ... ", dim[[1]]];
              Return[False];
             ];
           If[1 > idx2 || idx2 > dim[[2]],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " second index out of range! ",
                    "Allowed range: 1 ... ", dim[[2]]];
              Return[False];
             ];
           True
          ];

CheckSetting[patt:{parameter_, value_}, constraintName_String, isInitial_] :=
    Module[{outputParameters, modelPars},
           outputParameters = Parameters`GetOutputParameters[];
           If[MemberQ[outputParameters, parameter],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", parameter, " is a output parameter!"];
              Return[False];
             ];
           If[isInitial,
              modelPars = Parameters`FSModelParameters /. Parameters`FindAllParametersClassified[value];
              If[Intersection[modelPars, Parameters`GetModelParameters[]] =!= {},
                 Print["Warning: In constraint ", constraintName, ": ", InputForm[patt]];
                 Print["   ", modelPars, " on the r.h.s. are model parameters, which may initially be zero!"];
                ];
             ];
           True
          ];

CheckSetting[patt:( FlexibleSUSY`FSRestrictParameter[p_,__] | FlexibleSUSY`FSInitialSetting[p_,__] ),
             constraintName_String, __] :=
    Module[{outputParameters},
           outputParameters = Parameters`GetOutputParameters[];
           If[MemberQ[outputParameters, p],
              Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
              Print["   ", p, " is a output parameter!"];
              Return[False];
             ];
           True
          ];

CheckSetting[patt_, constraintName_String, __] :=
    Module[{},
           Print["Error: In constraint ", constraintName, ": ", InputForm[patt]];
           Print["   This is not a valid constraint setting!"];
           If[patt === Null,
              Print["   Maybe there is a trailing comma in the constraint list?"];
             ];
           False
          ];

CheckConstraint[settings_List, constraintName_String, isInitial_:False] :=
    CheckSetting[#,constraintName,isInitial]& /@ settings;

SanityCheck[settings_List, constraintName_String:""] :=
    Module[{setParameters, y,
            yukawas = {SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa}},
           setParameters = #[[1]]& /@ settings;
           (* check for unset Yukawa couplings *)
           For[y = 1, y <= Length[yukawas], y++,
               If[(ValueQ /@ yukawas)[[y]] &&
                  FreeQ[setParameters, yukawas[[y]]],
                  Print["Warning: Yukawa coupling ", yukawas[[y]],
                        " not set",
                        If[constraintName != "", " in the " <> constraintName, ""],
                        "."];
                 ];
              ];
          ];

CalculateScale[Null, _] := "";

CalculateScale[False, _] :=
    "ERROR(\"scale condition is always false!\");\n";

CalculateScale[True, _] :=
    "WARNING(\"scale condition is always true!\");\n";

GetSMMatchingScale[FlexibleSUSY`LowEnergyConstant[FlexibleSUSY`MT], qedqcd_String] :=
    qedqcd <> ".displayPoleMt()";

GetSMMatchingScale[FlexibleSUSY`LowEnergyConstant[FlexibleSUSY`MZ], qedqcd_String] :=
    qedqcd <> ".displayPoleMZ()";

CalculateScale[s:FlexibleSUSY`LowEnergyConstant[_], scaleName_String] :=
    scaleName <> " = " <> GetSMMatchingScale[s, "qedqcd"] <> ";\n";

CalculateScale[expr_, scaleName_String] :=
    Module[{result},
           result = Parameters`CreateLocalConstRefs[expr];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr, scaleName];
           Return[result];
          ];

CreateBetasForParsIn[expr_] :=
    Module[{pars},
           pars = Parameters`FSModelParameters /. Parameters`FindAllParametersClassified[expr];
           FlexibleSUSY`BETA /@ pars
          ];

CalculateScale[expr_Equal, scaleName_String] :=
    Module[{scaleReset, result},
           result = "const double currentScale = model->get_scale();\n"
                    <> "const auto beta_functions(model->calc_beta());\n\n";
           result = result <> Parameters`CreateLocalConstRefs[expr];
           result = result <> Parameters`CreateLocalConstRefsForBetas[CreateBetasForParsIn[expr]];
           result = result <> "\n";
           result = result <> CalculateScaleFromExpr[expr, scaleName];
           scaleReset = "ERROR(\"Overflow error during calculation of scale: \"\n"
                        <> "      << strerror(errno) << '\\n'\n"
                        <> "      << \"   current scale = \" << currentScale << '\\n'\n"
                        <> "      << \"   new scale = \" << scale << '\\n'\n"
                        <> "      << \"   resetting scale to \" << get_initial_scale_guess());\n";
           scaleReset = "#ifdef ENABLE_VERBOSE\n" <> IndentText[scaleReset] <> "#endif\n";
           scaleReset = scaleReset <> IndentText["scale = get_initial_scale_guess();\n"]
                        <> IndentText["errno = 0;\n"];
           result = result <> "if (errno == ERANGE) {\n" <> scaleReset <> "}\n";
           Return[result];
          ];

GetIndexedParameter[BetaFunction`BetaFunction[name_, CConversion`ScalarType[_], _]] := name;

GetIndexedParameter[BetaFunction`BetaFunction[name_[__], CConversion`MatrixType[_,rows_,cols_], _]] :=
    Flatten[Table[name[i,j], {i,1,rows}, {j,1,cols}]];

GetListOfIndexedParameters[] :=
    Flatten[GetIndexedParameter[#]& /@ allBetaFunctions];

StripConditionalExpression[ConditionalExpression[expr_, ___]] :=
    expr;

StripConditionalExpression[expr_] :=
    expr;

(* Approximately find scale where equation
     expr1 == expr2
   is fulfilled, using the beta function for each parameter p
     p(MX) = p(mu) + BETA[p] log(MX/mu)
*)
CalculateScaleFromExprSymb[Equal[expr1_, expr2_]] :=
    Module[{result, F1, F2, betaFunctions, solution, scale, parameters, parSeq},
           parameters = GetListOfIndexedParameters[];
           parSeq = Sequence @@ parameters;
           F1[parSeq] := expr1;
           F2[parSeq] := expr2;
           betaFunctions = FlexibleSUSY`BETA /@ parameters;
           solution = Solve[Log[scale/Global`currentScale] *
                            (betaFunctions . D[F1[parSeq] - F2[parSeq],
                                               {parameters}])
                            == F2[parSeq] - F1[parSeq], scale];
           If[solution === {{}},
              Print["Error: no solution found for ", expr1 == expr2];
              result = Null;,
              result = FullSimplify[StripConditionalExpression[solution[[1,1,2]]]];
             ];
           Return[result];
          ];

CalculateScaleFromExprSymb[False] := Null;

CalculateScaleFromExprSymb[True] := Global`currentScale;

CalculateScaleFromExpr[Equal[expr1_, expr2_], scaleName_String] :=
    Module[{result, solution},
           solution = CalculateScaleFromExprSymb[Equal[expr1, expr2]];
           If[solution === Null,
              result = "ERROR(\"no solution found for the equation " <>
                        ToString[expr1] <> " == " <> ToString[expr2] <> "\");\n";
              ,
              result = scaleName <> " = " <> CConversion`RValueToCFormString[solution] <> ";\n";
             ];
           Return[result];
          ];

CalculateScaleFromExpr[expr_, scaleName_String] :=
    scaleName <> " = " <> CConversion`RValueToCFormString[Parameters`DecreaseIndexLiterals[expr, Parameters`GetOutputParameters[]]] <> ";\n";

DefineAndDefaultInitialize[{t:FlexibleSUSY`Phase[_], _}] :=
    CConversion`CreateCType[GuessExtraParameterType[t]] <> " " <>
    ToValidCSymbolString[t] <> "{1.,0.};\n";

DefineAndDefaultInitialize[{t:Sign[_], _}] :=
    CConversion`CreateCType[GuessExtraParameterType[t]] <> " " <>
    ToValidCSymbolString[t] <> "{1};\n";

DefineAndDefaultInitialize[p:{_,_}] :=
    Parameters`CreateParameterDefinitionAndDefaultInitialize[p];

DefineInputParameters[inputParameters_List] :=
    StringJoin[DefineAndDefaultInitialize /@ inputParameters];

InitialGuessAtLowScaleGaugeCouplings[] :=
    Module[{result = ""},
           If[ValueQ[SARAH`hyperchargeCoupling],
              result = result <> Parameters`SetParameter[SARAH`hyperchargeCoupling,
                                                         "Sqrt(4. * Pi * alpha_sm(0))", "MODEL->"];
             ];
           If[ValueQ[SARAH`leftCoupling],
              result = result <> Parameters`SetParameter[SARAH`leftCoupling,
                                                         "Sqrt(4. * Pi * alpha_sm(1))", "MODEL->"];
             ];
           If[ValueQ[SARAH`strongCoupling],
              result = result <> Parameters`SetParameter[SARAH`strongCoupling,
                                                         "Sqrt(4. * Pi * alpha_sm(2))", "MODEL->"];
             ];
           result
          ];

IsFixedIn[par_, {p_, _}] :=
    Parameters`StripIndices[par] === Parameters`StripIndices[p];

IsFixedIn[par_, FlexibleSUSY`FSMinimize[parameters_List, _]] :=
    MemberQ[Parameters`StripIndices /@ parameters, Parameters`StripIndices[par]];

IsFixedIn[par_, FlexibleSUSY`FSFindRoot[parameters_List, _]] :=
    MemberQ[Parameters`StripIndices /@ parameters, Parameters`StripIndices[par]];

IsFixedIn[par_, FlexibleSUSY`FSSolveEWSBFor[parameters___]] :=
    MemberQ[Parameters`StripIndices /@ Flatten[{parameters}], Parameters`StripIndices[par]];

IsFixedIn[par_, (FlexibleSUSY`FSRestrictParameter | FlexibleSUSY`FSInitialSetting)[p_, ___]] :=
    MemberQ[Parameters`StripIndices[p], Parameters`StripIndices[par]];

IsFixedIn[par_, p___] :=
    Block[{},
          Print["Error: This is not a valid constraint setting: ", p];
          Quit[1];
         ];

IsFixed[par_, constraint_List] :=
    Or @@ (IsFixedIn[par, #]& /@ constraint);

RestrictScale[{minimumScale_, maximumScale_}, scaleName_String:"scale"] :=
    Module[{result = "", value},
           If[NumericQ[minimumScale],
              value = CConversion`RValueToCFormString[minimumScale];
              result = result <>
                       "\
if (" <> scaleName <> " < " <> value <> ") {
#ifdef ENABLE_VERBOSE
" <> IndentText["WARNING(\"" <> scaleName <> " < " <> value <> "\");"] <> "
#endif
" <> IndentText[scaleName <> " = " <> value] <> ";
}
";
             ];
           If[NumericQ[maximumScale],
              value = CConversion`RValueToCFormString[maximumScale];
              result = result <>
                       "\
if (" <> scaleName <> " > " <> value <> ") {
#ifdef ENABLE_VERBOSE
" <> IndentText["WARNING(\"" <> scaleName <> " > " <> value <> "\");"] <> "
#endif
" <> IndentText[scaleName <> " = " <> value] <> ";
}
";
             ];
           Return[result];
          ];

PerturbativityCheckSnippet[par_, thresh_, model_String, problem_String] :=
    Module[{parStr, parEnum, threshStr},
           parStr = CConversion`RValueToCFormString[par];
           parEnum = FlexibleSUSY`FSModelName <> "_info::" <> Parameters`CreateEnumName[par];
           threshStr = CConversion`RValueToCFormString[thresh];
           "\
if (MaxAbsValue(" <> parStr <> ") > " <> threshStr <> ") {
   " <> problem <> " = true;
   " <> model <> "->get_problems().flag_non_perturbative_parameter(" <> parEnum <> ", MaxAbsValue(" <> parStr <> "), " <> model <> "->get_scale(), " <> threshStr <> ");
} else {
   " <> model <> "->get_problems().unflag_non_perturbative_parameter(" <> parEnum <> ");
}
"
          ];

CheckPerturbativityForParameter[par_, thresh_, model_String:"model", problem_String:"problem"] :=
    Module[{snippet, parStr, threshStr},
           Utils`StringJoinWithSeparator[
               PerturbativityCheckSnippet[#, thresh, model, problem]& /@ Parameters`DecomposeParameter[par, GetType[par]], "\n"]
          ];

CheckPerturbativityForParameters[pars_List, thresh_] :=
    Parameters`CreateLocalConstRefs[pars] <> "\n" <>
    StringJoin[CheckPerturbativityForParameter[#,thresh]& /@ pars];

SaveValue[par_[idx__]] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par],
            parStrIdx = CConversion`ToValidCSymbolString[par[idx]], oldVal},
           oldVal = "old_" <> parStrIdx;
           "const auto " <> oldVal <> " = MODELPARAMETER(" <> parStr <> ")" <>
           "(" <> Utils`StringJoinWithSeparator[ToString /@ {idx},","] <> ");\n" <>
           "const auto save_" <> parStrIdx <> " = make_raii_guard([this," <> oldVal <> "]{ " <>
           Parameters`SetParameter[par[idx], oldVal, "MODEL->"] <> " });\n"
          ];

SaveValue[par_] :=
    Module[{parStr = CConversion`ToValidCSymbolString[par], oldVal},
           oldVal = "old_" <> parStr;
           "const auto " <> oldVal <> " = MODELPARAMETER(" <> parStr <> ");\n" <>
           "const auto save_" <> parStr <> " = make_raii_guard([this," <> oldVal <> "]{ " <>
           Parameters`SetParameter[par, oldVal, "MODEL->"] <> " });\n"
          ];

SetTemporarily[settings_List] :=
    Module[{tempSettings = Cases[settings, {FlexibleSUSY`FSTemporary[p_], v_} :> {p,v}],
            set, savedVals},
           If[tempSettings === {}, Return[""];];
           set = ApplyConstraints[tempSettings];
           savedVals = Utils`StringJoinWithSeparator[SaveValue[#[[1]]]& /@ tempSettings, "\n"];
           "// temporary parameter re-definitons\n" <>
           savedVals <> "\n{\n" <> IndentText[set] <> "}"
          ];

SetRestrictions[settings_List] :=
    Module[{initSettings = Cases[settings, (FlexibleSUSY`FSRestrictParameter | FlexibleSUSY`FSInitialSetting)[__]]},
           If[initSettings === {},
              "",
              "// initial settings / parameter restrictions\n{\n" <> IndentText[
                  ApplyConstraints[ExpandRestrictParameter /@ initSettings]
              ] <> "\n}"
             ]
          ];

InitialApplyConstraint[settings_List] :=
    SetTemporarily[settings] <> "\n" <> SetRestrictions[settings];

End[];

EndPackage[];
