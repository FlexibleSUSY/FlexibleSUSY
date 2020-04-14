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

BeginPackage["EWSB`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "TreeMasses`", "WriteOut`", "Utils`"}];

FilterOutLinearDependentEqs::usage="returns linearly independent equations";

FilterOutIndependentEqs::usage = "returns equations that depend on the
given list of parameters.  I.e. equations, that do not depend on the
given list of parameters are omitted from the output.";

GetLinearlyIndependentEqs::usage="Removes linearly dependent EWSB equations
from a list of equations";

FindSolutionAndFreePhases::usage="Finds solution to the EWSB and free
phases / signs."

CreateEWSBEqPrototype::usage="creates C function prototype for a
given EWSB equation";

CreateEWSBEqFunction::usage="creates C function definition for a
given EWSB equation";

FillArrayWithEWSBEqs::usage="fills array of doubles with the values
of the EWSB equations";

FillInitialGuessArray::usage="fills a C array with initial values for the
EWSB eqs. solver";

CreateMemberTreeLevelEwsbSolver::usage="Converts tree-level EWSB solutions
to C form";

CreateTreeLevelEwsbSolver::usage="Converts tree-level EWSB solutions
to C form";

CreateEWSBRootFinders::usage="Creates comma separated list of GSL root
finders";

SetEWSBSolution::usage="sets the model parameters to the solution
provided by the solver";

SetTreeLevelSolution::usage="sets the model parameters to the
tree-level solution";

FillArrayWithParameters::usage="fill an array with parameters";

DivideTadpolesByVEV::usage="Divides an array of tadpoles by their
corresponding VEV";

CreateEwsbSolverWithTadpoles::usage="Solve EWSB eqs. including
tadpoles (one step, no iteration)";

GetEWSBParametersFromVector::usage="Create local copies of EWSB
 output parameters from an Eigen vector";

SetEWSBParametersFromLocalCopies::usage="Set model parameters from
local copies";

CreateEWSBParametersInitializationList::usage="Creates initialization
list with EWSB output parameters";

CreateEWSBParametersInitializationComma::usage="Creates initialization
list with EWSB output parameters";

CreateEWSBParametersInitialization::usage="Creates initialization
of EWSB output parameters";

GetValidEWSBInitialGuesses::usage="Remove invalid initial guess settings";

GetValidEWSBSubstitutions::usage="Remove invalid EWSB substitutions";

SetModelParametersFromEWSB::usage="Set model parameters from EWSB solution
using the given substitutions";

ApplyEWSBSubstitutions::usage="Set model parameters according to the
given list of substitutions";

SolveEWSBIgnoringFailures::usage="Solve EWSB conditions without
flagging a problem if no solution is found.";

Begin["`Private`"];

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<EWSB>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

AppearsInEquationOnlyAs[parameter_, equation_, function_] :=
    FreeQ[equation /. function[parameter] :> Unique[CConversion`ToValidCSymbolString[parameter]], parameter];

AppearsOnlySquaredInEquation[parameter_, equation_] :=
    AppearsInEquationOnlyAs[parameter, equation, Power[#,2]&];

AppearsOnlyAbsSquaredInEquation[parameter_, equation_] :=
    AppearsInEquationOnlyAs[parameter, equation, (Susyno`LieGroups`conj[#] #)&] ||
    AppearsInEquationOnlyAs[parameter, equation, (# Susyno`LieGroups`conj[#])&];

AppearsNotInEquation[parameter_, equation_] :=
    FreeQ[equation, parameter];

CheckInEquations[parameter_, statement_, equations_List] :=
    And @@ (statement[parameter,#]& /@ equations);

AreLinearDependent[{eq1_, eq2_}, parameters_List] :=
    Module[{frac = Simplify[eq1/eq2 /. FlexibleSUSY`tadpole[_] -> 0],
            pars},
           (* ignore parameter heads Re[], Im[], Abs[], Phase[] *)
           pars = parameters /. { Re[p_] :> p, Im[p_] :> p,
                                  Abs[p_] :> p, FlexibleSUSY`Phase[p_] :> p };
           And @@ (FreeQ[frac,#]& /@ pars)
          ];

FilterOutLinearDependentEqs[{}, _List] := {};

FilterOutLinearDependentEqs[{eq_}, _List] := {eq};

FilterOutLinearDependentEqs[{eq_, rest__}, parameters_List] :=
    If[Or @@ (AreLinearDependent[#,parameters]& /@ ({eq,#}& /@ {rest})),
       (* leave out eq and check rest *)
       FilterOutLinearDependentEqs[{rest}, parameters],
       (* keep eq and check rest *)
       {eq, Sequence @@ FilterOutLinearDependentEqs[{rest}, parameters]}
      ];

FilterOutIndependentEqs[eqs_List, pars_List] :=
    DeleteDuplicates @ Flatten @ Join[FilterOutIndependentEqs[eqs,#]& /@ pars];

FilterOutIndependentEqs[eqs_List, p_] :=
    Select[eqs, (!FreeQ[#,p])&];

GetLinearlyIndependentEqs[eqs_List, parameters_List, substitutions_List:{}] :=
    Module[{eqsToSolve, indepEqsToSolve, eqsToKeep},
           If[substitutions =!= {},
              eqsToSolve = Parameters`ReplaceAllRespectingSARAHHeads[eqs, substitutions];,
              eqsToSolve = eqs;
             ];
           indepEqsToSolve = FilterOutLinearDependentEqs[eqsToSolve, parameters];
           eqsToKeep = Position[eqsToSolve, p_ /; MemberQ[indepEqsToSolve, p]];
           Extract[eqs, eqsToKeep]
          ];

CreateEWSBEqPrototype[higgs_] :=
    Module[{result = "", i, ctype},
           ctype = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           For[i = 1, i <= TreeMasses`GetDimension[higgs], i++,
               result = result <> ctype <> " get_ewsb_eq_" <>
                        CConversion`ToValidCSymbolString[higgs] <>
                        "_" <> ToString[i] <> "() const;\n";
              ];
           Return[result];
          ];

(* Creates a get_ewsb_eq_<higgs>() function for each Higgs multiplet
   entry.

   This function assumes, that the ordering of the EWSB
   eqs. corresponds to the Higgs multiplet entries.  If there are less
   equations than Higgs multiplet entries, the remaining equations
   will be set to 0.  If there are more EWSB equations than Higgs
   bosons, the remaining equations will be ignored.
 *)
CreateEWSBEqFunction[higgs_, equation_List] :=
    Module[{result = "", body = "", dim, i, eq, dimEq, type, ctype},
           dim = TreeMasses`GetDimension[higgs];
           dimEq = Length[equation];
           type = CConversion`ScalarType[CConversion`realScalarCType];
           ctype = CConversion`CreateCType[type];
           If[dim < dimEq,
              Print["Warning: number of Higgs bosons (", dim,
                    ") < number of EWSB eqs. (",dimEq,")."
                    "The EWSB eqs. ", dim, "...", dimEq,
                    " will be ignored"];
             ];
           If[dim > dimEq,
              Print["Warning: number of physical Higgs bosons (", dim,
                    ") > number of EWSB eqs. (",dimEq,").",
                    "The EWSB eqs. for the fields ", higgs, "(n), n >= ",
                    dimEq, ", will be set to zero."];
             ];
           For[i = 1, i <= dim, i++,
               result = result <>
                        ctype <> " CLASSNAME::get_ewsb_eq_" <>
                        CConversion`ToValidCSymbolString[higgs] <>
                        "_" <> ToString[i] <> "() const\n{\n";
               If[i <= dimEq,
                  eq = equation[[i]];,
                  eq = 0;
                 ];
               body = Parameters`CreateLocalConstRefsForInputParameters[eq, "LOCALINPUT"] <>
                      "\n" <> ctype <> " result = " <>
                      CastTo[RValueToCFormString[eq], type] <> ";\n";
               body = body <> "\nreturn result;\n";
               body = IndentText[WrapLines[body]];
               result = result <> body <> "}\n\n"
              ];
           Return[result];
          ];

FindFreePhase[parameter_, freePhases_] :=
    Module[{phases},
           phases = Cases[freePhases, FlexibleSUSY`Phase[parameter] | Sign[parameter]];
           If[phases === {}, Null, phases[[1]]]
          ];

WrapPhase[phase_ /; phase === Null, str_String, macro_String:"LOCALINPUT"] :=
    str;

WrapPhase[phase_, str_String, macro_String:"LOCALINPUT"] :=
    macro <> "(" <> CConversion`ToValidCSymbolString[phase] <> ")*Abs(" <> str <> ")";

GetValueWithPhase[parameter_, gslIntputVector_String, index_Integer, freePhases_List] :=
    Module[{freePhase, gslInput},
           freePhase = FindFreePhase[parameter, freePhases];
           gslInput = gslIntputVector <> "(" <> ToString[index] <> ")";
           WrapPhase[freePhase, gslInput, "INPUT"]
          ];

FillArrayWithEWSBEqs[higgs_, gslOutputVector_String] :=
    Module[{i, result = "", par, dim},
           dim = TreeMasses`GetDimension[higgs];
           For[i = 1, i <= dim, i++,
               result = result <> gslOutputVector <> "[" <> ToString[i-1] <>
                        "] = " <> "get_ewsb_eq_" <>
                        CConversion`ToValidCSymbolString[higgs] <> "_" <>
                        ToString[i] <> "();\n";
              ];
           Return[result];
          ];

GetValidEWSBInitialGuesses[initialGuess_List] :=
    Module[{i},
           For[i = 1, i <= Length[initialGuess], i++,
               If[!MatchQ[initialGuess[[i]], {_,_}],
                  Print["Warning: ignoring invalid initial guess: ", initialGuess[[i]]];
                 ];
              ];
           Cases[initialGuess, {_,_}]
          ];

GetValidEWSBSubstitutions[substitutions_List] :=
    Module[{i},
           For[i = 1, i <= Length[substitutions], i++,
               If[!MatchQ[substitutions[[i]], {_,_}],
                  Print["Warning: ignoring invalid EWSB substitution: ", substitutions[[i]]];
                 ];
              ];
           Cases[substitutions, {_,_}]
          ];

InitialGuessFor[par_, initialGuesses_List:{}] :=
    Module[{guess},
           If[initialGuesses === {} || !MemberQ[#[[1]]& /@ initialGuesses, par],
              If[Parameters`IsRealParameter[par], guess = par, guess = Abs[par]];,
              guess = Cases[initialGuesses, {p_ /; p === par, val_} :> val];
              If[Length[guess] > 1,
                 Print["Warning: multiple initial guesses given for ", par];
                ];
              guess = First[guess];
             ];
           guess
          ];

FillInitialGuessArray[parametersFixedByEWSB_List, initialGuessValues_List:{}, arrayName_String:"x_init"] :=
    Module[{i, guesses, result = ""},
           guesses = InitialGuessFor[#, initialGuessValues]& /@ parametersFixedByEWSB;
           result = result <> Parameters`CreateLocalConstRefs[guesses];
           For[i = 1, i <= Length[guesses], i++,
               result = result <> arrayName <> "[" <> ToString[i-1] <> "] = " <>
                        CConversion`RValueToCFormString[guesses[[i]]] <>
                        ";\n";
              ];
           Return[result];
          ];

MakeParameterUnique[(Re|Im)[par_]] :=
    { MakeParameterUnique[par],
      Rule[Re[par], CConversion`ToValidCSymbol[Re[par]]],
      Rule[Im[par], CConversion`ToValidCSymbol[Im[par]]] };
MakeParameterUnique[SARAH`L[par_]] := Rule[SARAH`L[par], CConversion`ToValidCSymbol[SARAH`L[par]]];
MakeParameterUnique[SARAH`B[par_]] := Rule[SARAH`B[par], CConversion`ToValidCSymbol[SARAH`B[par]]];
MakeParameterUnique[SARAH`T[par_]] := Rule[SARAH`T[par], CConversion`ToValidCSymbol[SARAH`T[par]]];
MakeParameterUnique[SARAH`Q[par_]] := Rule[SARAH`Q[par], CConversion`ToValidCSymbol[SARAH`Q[par]]];
MakeParameterUnique[par_]          :=
    { MakeParameterUnique[SARAH`L[par]],
      MakeParameterUnique[SARAH`B[par]],
      MakeParameterUnique[SARAH`T[par]],
      MakeParameterUnique[SARAH`Q[par]] };

MakeParametersUnique[parameters_List] :=
    Flatten[MakeParameterUnique /@ parameters];

ComplexParameterReplacementRules[eqs_List, pars_List] :=
    Join[ComplexParameterReplacementRules[eqs,#]& /@ pars];

(* returns replacement rules which, if appied to eqs, lead to
   equations that are free of FlexibleSUSY`Phase[par] *)
ComplexParameterReplacementRules[eqs_List, par_] :=
    Module[{rules, replacedEqs},
           rules = {Rule[par,Abs[par] FlexibleSUSY`Phase[par]],
                    Rule[SARAH`Conj[par],Abs[par]/FlexibleSUSY`Phase[par]],
                    Rule[Susyno`LieGroups`conj[par],Abs[par]/FlexibleSUSY`Phase[par]]};
           replacedEqs = Simplify[eqs /. rules];
           If[FreeQ[replacedEqs, FlexibleSUSY`Phase[par]],
              rules,
              {}
             ]
          ];

SplitRealAndImagParts[eqs_List, pars_List] :=
    Module[{parsWithoutHeads, uniqueRules, uniqueEqs, uniquePars, result},
           parsWithoutHeads = Cases[pars, (Re | Im | SARAH`L | SARAH`B | SARAH`T | SARAH`Q)[p_] | p_ :> p];
           uniqueRules = DeleteDuplicates @ Flatten[{
               Rule[SARAH`L[#], CConversion`ToValidCSymbol[SARAH`L[#]]],
               Rule[SARAH`B[#], CConversion`ToValidCSymbol[SARAH`B[#]]],
               Rule[SARAH`T[#], CConversion`ToValidCSymbol[SARAH`T[#]]],
               Rule[SARAH`Q[#], CConversion`ToValidCSymbol[SARAH`Q[#]]]
           }& /@ parsWithoutHeads];
           uniqueEqs = eqs /. uniqueRules;
           uniquePars = pars /. uniqueRules;
           result = uniqueEqs /. DeleteDuplicates[Cases[uniquePars, Re[p_] | Im[p_] :> Rule[p,Re[p]+I Im[p]]]];
           result /. (Reverse /@ uniqueRules)
          ];

SimplifyEwsbEqs[equations_List, parametersFixedByEWSB_List] :=
    Module[{realParameters, complexParameters, simplificationRules,
            renamedEqs, splitEqs},
           DebugPrint["Splitting Re[] and Im[] within EWSB eqs. ..."];
           splitEqs = SplitRealAndImagParts[equations, parametersFixedByEWSB];
           realParameters = Select[parametersFixedByEWSB, Parameters`IsRealParameter[#]&];
           DebugPrint["real parameters: ", realParameters];
           complexParameters = Complement[parametersFixedByEWSB, realParameters];
           DebugPrint["complex parameters: ", complexParameters];
           (* make parameters unique *)
           uniqueParameters = MakeParametersUnique[parametersFixedByEWSB];
           DebugPrint["Making parameters unique via: ", uniqueParameters];
           realParameters = realParameters /. uniqueParameters;
           complexParameters = complexParameters /. uniqueParameters;
           renamedEqs = splitEqs /. uniqueParameters;
           DebugPrint["EWSB eqs. with unique parameters: ", renamedEqs];
           simplificationRules =
               Flatten[Join[{Rule[SARAH`Conj[#],#],
                             Rule[Susyno`LieGroups`conj[#],#]}& /@ realParameters,
                            ComplexParameterReplacementRules[renamedEqs, complexParameters]
                           ]
                      ];
           DebugPrint["Simplification rules: ", simplificationRules];
           (* substitute back *)
           uniqueParameters = Reverse /@ uniqueParameters;
           renamedEqs /. simplificationRules /. uniqueParameters
          ];

FindIndependentSubset[equations_List, {}] := {};

FindIndependentSubset[{}, parameters_List] := {};

FindIndependentSubset[equations_List, parameters_List] /; Length[parameters] > Length[equations] :=
    Module[{subsets = Subsets[parameters, {Length[equations]}], res},
           (* find independent subsets for all parameter combinations *)
           res = Select[FindIndependentSubset[equations, #]& /@ subsets, (# =!= {})&];
           (* select the largest one *)
           res = Sort[res, Length[First[#1]] > Length[First[#2]]&];
           If[res === {}, {}, First[res]]
          ];

FindIndependentSubset[{eq_}, {par_}] :=
    If[FreeQ[eq, par],
       {},
       {{ {eq}, {par} }}
      ];

FindIndependentSubset[equations_List, parameters_List] :=
    Module[{equationSubsets, numberOfEquations, parameterSubsets,
            numberOfParameters, e, p, result = {}, isFreeOf},
           numberOfEquations = Length[equations];
           numberOfParameters = Length[parameters];
           equationSubsets = Subsets[equations, {1, numberOfEquations - 1}];
           parameterSubsets = Subsets[parameters, {1, numberOfParameters - 1}];
           For[e = 1, e <= Length[equationSubsets], e++,
               For[p = 1, p <= Length[parameterSubsets], p++,
                   isFreeOf =
                   And @@ (FreeQ[equationSubsets[[e]], #] & /@
                           parameterSubsets[[p]]);
                   If[isFreeOf,
                      AppendTo[
                          result, {equationSubsets[[e]],
                                   Complement[parameters, parameterSubsets[[p]]]}]];
                  ];
              ];
           result = Select[result, (Length[#[[1]]] == Length[#[[2]]]) &];
           Return[result];
          ];

FindMinimumByteCount[{}] := {};

FindMinimumByteCount[lst_List] :=
    First[Sort[lst, ByteCount[#1] < ByteCount[#2]]];

IsNoSolution[expr_] :=
    Head[expr] === Solve || Flatten[expr] === {};

TimeConstrainedEliminate[eqs_List, par_] :=
    Module[{result},
           result = TimeConstrained[Eliminate[eqs, par],
                                    FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
           If[result === {},
              DebugPrint["unable to eliminate ", par, ": ", eqs];
             ];
           result
          ];

TimeConstrainedSolve[eq_, par_] :=
    Module[{result, independentEq = eq, Selector},
           If[Head[eq] === List,
              Selector = Function[e, If[Head[par] === List,
                                        And @@ (Function[p,FreeQ[e,p]] /@ par),
                                        FreeQ[e,par]
                                       ]];
              independentEq = Select[eq, (!Selector[#])&];,
              If[Head[par] =!= List,
                 If[FreeQ[eq, par],
                    Print["Error: parameter ", par, " does not appear in the equation"];
                    Return[{}];
                   ];
                ];
             ];
           result = TimeConstrained[Solve[independentEq, par],
                                    FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
           If[result === {} || result === {{}},
              Off[Reduce::nsmet];
              result = TimeConstrained[{ToRules[Reduce[independentEq, par]]},
                                        FlexibleSUSY`FSSolveEWSBTimeConstraint, {}];
              If[Length[result] == 1 &&
                 (Head[result[[1]]] === ToRules || Head[result[[1]]] === Reduce),
                 result = {};
                ];
              On[Reduce::nsmet];
             ];
           result
          ];

EliminateOneParameter[{}, _List] := {};

EliminateOneParameter[{eq_}, {p_}] :=
    Block[{},
          DebugPrint["eliminating ", p, ": ", eq];
          {TimeConstrainedSolve[eq, p]}
         ];

(* solves the two equations for `par' and returns the simpler solution *)
SolveRest[eq1_, eq2_, par_] :=
    Module[{rest = {}, solution},
           If[FreeQ[eq1, par],
              DebugPrint[par, " does not appear in equation: ", eq1];,
              DebugPrint["solving rest for ", par, ": ", eq1];
              solution = TimeConstrainedSolve[eq1, par];
              If[IsNoSolution[solution],
                 DebugPrint["Failed"];,
                 DebugPrint["Solution found: ", solution];
                 AppendTo[rest, solution];
                ];
             ];
           If[FreeQ[eq2, par],
              DebugPrint[par, " does not appear in equation: ", eq2];,
              DebugPrint["solving rest for ", par, ": ", eq2];
              solution = TimeConstrainedSolve[eq2, par];
              If[IsNoSolution[solution],
                 DebugPrint["Failed"];,
                 DebugPrint["Solution found: ", solution];
                 AppendTo[rest, solution];
                ];
             ];
           FindMinimumByteCount[rest]
          ];

EliminateOneParameter[{eq1_, eq2_}, {p1_, p2_}] :=
    Module[{reduction = {{}, {}}, solution, rest},
           DebugPrint["Trying to eliminate one of the parameters ",
                      {p1,p2}, " from the eqs.: ",
                      {eq1,eq2}];
           If[FreeQ[{eq1, eq2}, p1],
              Print["Error: EWSB output parameter ", p1, " does not appear in the EWSB eqs."];
              Return[{}];
             ];
           If[FreeQ[{eq1, eq2}, p2],
              Print["Error: EWSB output parameter ", p2, " does not appear in the EWSB eqs."];
              Return[{}];
             ];
           (* special case: no elimination needed *)
           If[!FreeQ[{eq1},p1] && FreeQ[{eq1},p2] &&
              !FreeQ[{eq2},p2] && FreeQ[{eq2},p1],
              DebugPrint["The two equations are independent of each other."];
              DebugPrint["Step 1: solving for ", p1, ": ", eq1];
              reduction[[1]] = TimeConstrainedSolve[{eq1}, p1];
              If[IsNoSolution[reduction[[1]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[1]]];
                ];
              DebugPrint["Step 2: solving for ", p2, ": ", eq2];
              reduction[[2]] = TimeConstrainedSolve[{eq2}, p2];
              If[IsNoSolution[reduction[[2]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[2]]];
                ];
              DebugPrint["Full solution: ", reduction];
              Return[reduction];
             ];
           If[!FreeQ[{eq1},p2] && FreeQ[{eq1},p1] &&
              !FreeQ[{eq2},p1] && FreeQ[{eq2},p2],
              DebugPrint["The two equations are independent of each other."];
              DebugPrint["Step 1: solving for ", p2, ": ", eq1];
              reduction[[1]] = TimeConstrainedSolve[{eq1}, p2];
              If[IsNoSolution[reduction[[1]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[1]]];
                ];
              DebugPrint["Step 2: solving for ", p1, ": ", eq2];
              reduction[[2]] = TimeConstrainedSolve[{eq2}, p1];
              If[IsNoSolution[reduction[[2]]],
                 DebugPrint["Failed"];
                 Return[{}];,
                 DebugPrint["Solution: ", reduction[[2]]];
                ];
              DebugPrint["Full solution: ", reduction];
              Return[reduction];
             ];
           DebugPrint["eliminate ", p1, " and solve for ", p2, ": ", {eq1, eq2}];
           reduction[[1]] =
           TimeConstrainedSolve[TimeConstrainedEliminate[{eq1, eq2}, p1], p2];
           DebugPrint["eliminate ", p2, " and solve for ", p1, ": ", {eq1, eq2}];
           reduction[[2]] =
           TimeConstrainedSolve[TimeConstrainedEliminate[{eq1, eq2}, p2], p1];
           If[IsNoSolution[reduction[[1]]] || IsNoSolution[reduction[[2]]],
              DebugPrint["Failed"];
              Return[{}];
             ];
           If[ByteCount[reduction[[1]]] <= ByteCount[reduction[[2]]],
              DebugPrint["continue with solution 1, because it is simpler:",
                         reduction[[1]]];
              rest = SolveRest[eq1, eq2, p1];
              If[IsNoSolution[rest],
                 DebugPrint["could not solve rest for solution 1"];
                 Return[{}];
                ];
              Return[{reduction[[1]], rest}];
              ,
              DebugPrint["continue with solution 2, because it is simpler:",
                         reduction[[2]]];
              rest = SolveRest[eq1, eq2, p2];
              If[IsNoSolution[rest],
                 DebugPrint["could not solve rest for solution 2"];
                 Return[{}];
                ];
              Return[{reduction[[2]], rest}];
             ];
          ];

EliminateOneParameter[equations_List, parameters_List] :=
    Module[{independentSubset, reducedEqs, reducedPars, reducedSolution,
            complementEq, complementPar, complementSolution,
            largestIndependentSubset, s},
           independentSubset = FindIndependentSubset[equations, parameters];
           If[independentSubset === {},
              DebugPrint["EWSB equations are not reducible"];
              Return[{}];
             ];
           (* search for the largest independent equation subset *)
           If[Length[independentSubset] == 1,
              largestIndependentSubset = independentSubset[[1]];,
              largestIndependentSubset = independentSubset[[1]];
              For[s = 1, s <= Length[independentSubset], s++,
                  If[Length[independentSubset[[s,2]]] > Length[largestIndependentSubset[[2]]],
                     largestIndependentSubset = independentSubset[[s]];
                    ];
                 ];
             ];
           reducedEqs = largestIndependentSubset[[1]];
           reducedPars = largestIndependentSubset[[2]];
           reducedSolution = EliminateOneParameter[reducedEqs, reducedPars];
           If[reducedSolution === {},
              DebugPrint["Could not solve reduced EWSB eqs. subset"];
              Return[{}];
             ];
           complementEq = Complement[equations, reducedEqs];
           complementPar = Complement[parameters, reducedPars];
           complementSolution = EliminateOneParameter[complementEq, complementPar];
           If[complementSolution === {},
              DebugPrint["Could not solve remaining EWSB eqs. subset"];
              Return[{}];
             ];
           DebugPrint["Solution = ", Join[reducedSolution, complementSolution]];
           Join[reducedSolution, complementSolution]
          ];

ToMathematicaSolutionFormat[{}] := {};

ToMathematicaSolutionFormat[sol_List] :=
    Tuples[Flatten /@ sol];

FindSolution[equations_List, parametersFixedByEWSB_List] :=
    Module[{simplifiedEqs, makeParsUnique, solution,
            uniquePars, uniqueEqs},
           DebugPrint["Simplifying the EWSB eqs. ..."];
           simplifiedEqs = SimplifyEwsbEqs[equations, parametersFixedByEWSB];
           simplifiedEqs = (# == 0)& /@ simplifiedEqs;
           DebugPrint["Simplified EWSB eqs.: ", simplifiedEqs];
           (* replace non-symbol parameters by unique symbols *)
           makeParsUnique = MakeParametersUnique[parametersFixedByEWSB];
           uniquePars = parametersFixedByEWSB /. makeParsUnique;
           uniqueEqs = simplifiedEqs /. makeParsUnique;
           DebugPrint["Eliminating the parameters ", uniquePars];
           solution = ToMathematicaSolutionFormat @ EliminateOneParameter[uniqueEqs, uniquePars];
           If[solution === {},
              DebugPrint["Trying Mathematica's Solve[] with time constraint of ",
                         FlexibleSUSY`FSSolveEWSBTimeConstraint, " seconds"];
              solution = TimeConstrainedSolve[uniqueEqs, uniquePars];
             ];
           (* substitute back unique parameters *)
           makeParsUnique = Reverse /@ makeParsUnique;
           solution /. makeParsUnique
          ];

StripSign[Times[int_?NumericQ,expr_]] := Abs[int] expr;

StripSign[expr_] := expr;

ReduceSolution[{}] := {{},{}};

ReduceSolution[{{}}] := {{},{}};

SignOrPhase[par_] :=
    If[Parameters`IsRealParameter[par],
       Sign[par],
       FlexibleSUSY`Phase[par]];

ReduceTwoSolutions[sol1_, sol2_] :=
    Module[{par, signOrPhase, reducedSolution},
           par = sol1[[1]];
           DebugPrint["Reducing solutions for ", par, "..."];
           signOrPhase = SignOrPhase[par];
           If[PossibleZeroQ[sol1[[2]] - sol2[[2]]],
              DebugPrint["The two solutions for ", par, " are identical"];
              Return[{sol1}];
             ];
           If[!PossibleZeroQ[sol1[[2]] + sol2[[2]]],
              Print["Warning: cannot reduce solution for ", par];
              Print["   because the two solutions are not related by a global sign."];
              Return[{}];
             ];
           DebugPrint["The two solutions for ", par,
                      " are related by a global sign/phase: ",
                      sol1, ", ", sol2];
           reducedSolution = sol1 /.
               Rule[p_, expr_] :>
               Rule[p, signOrPhase StripSign[expr]];
           DebugPrint["=> the reduced solution is: ",
                      {reducedSolution, signOrPhase}];
           {reducedSolution, signOrPhase}
          ];

ReduceSolution[{sol_}] := {{sol},{}};

ReduceSolution[{sol1_, sol2_}] :=
    Module[{reducedSolution = {}, freePhases = {}, s, red},
           DebugPrint["Reducing the two solutions: ", {sol1,sol2}];
           For[s = 1, s <= Length[sol1], s++,
               red = ReduceTwoSolutions[sol1[[s]], sol2[[s]]];
               Switch[Length[red],
                      1, AppendTo[reducedSolution, red[[1]]];,
                      2, AppendTo[reducedSolution, red[[1]]];
                         AppendTo[freePhases, red[[2]]];
                     ];
              ];
           If[Length[reducedSolution] != Length[sol1],
              Print["Warning: analytic reduction of EWSB solutions failed."];
              Return[{{},{}}];
             ];
           Return[{reducedSolution, freePhases}];
          ];

ReduceSolution[solution_List] :=
    Module[{},
           Print["Error: cannot reduce the solution ", solution];
           Print["   because there are more than two solutions"];
           {{},{}}
          ];

FindSolutionAndFreePhases[equations_List, parametersFixedByEWSB_List, outputFile_String:"", substitutions_List:{}] :=
    Module[{eqsToSolve, solution, reducedSolution, freePhases},
           If[substitutions =!= {},
              eqsToSolve = Parameters`ReplaceAllRespectingSARAHHeads[equations, substitutions];,
              eqsToSolve = equations;
             ];
           solution = FindSolution[eqsToSolve, parametersFixedByEWSB];
           {reducedSolution, freePhases} = ReduceSolution[solution];
           DebugPrint["The full, reduced solution to the EWSB eqs. is:",
                      reducedSolution];
           If[outputFile != "",
              If[reducedSolution === {},
                 DebugPrint["Writing full, non-reduced solution to file ", outputFile];
                 Put[solution, outputFile];,
                 DebugPrint["Writing reduced solution to file ", outputFile];
                 Put[reducedSolution, outputFile];
                ];
             ];
           Return[{Flatten[reducedSolution], freePhases}];
          ];

GetFreeParameter[Rule[a_,b_]] := b;
GetFreeParameter[RuleDelayed[a_,b_]] := b;
GetFreeParameter[p_] := Print["Error: not a rule: ", p];

GetFixedParameter[Rule[a_,b_]] := a;
GetFixedParameter[RuleDelayed[a_,b_]] := a;
GetFixedParameter[p_] := Print["Error: not a rule: ", p];

(* Replaces each symbol in `symbols' in `expr' by
   Transformator[symbol].  Heads to be protected can be given in
   `protectedHeads'.

   Example:

   In[]:= TransformSymbolsIn[a + A[a], a, Unique, {A}]
   Out[]= a$191 + A[a]
 *)
TransformSymbolsIn[expr_, symbols_List, Transformator_, protectedHeads_:{SARAH`L, SARAH`B, SARAH`T, SARAH`Q}] :=
    Module[{transformationRules, protectedSymbols, protectionRules},
           transformationRules =  Rule[#, Transformator[#]]& /@ symbols;
           protectedHeadsPattern = Alternatives @@ (Blank /@ protectedHeads);
           protectedSymbols = Complement[
               DeleteDuplicates @ Extract[expr, Position[expr, protectedHeadsPattern]],
               symbols
           ];
           protectionRules = Rule[#, Unique[CConversion`ToValidCSymbolString[#]]]& /@ protectedSymbols;
           transformationRules = transformationRules /. protectionRules;
           expr /. protectionRules /. transformationRules /. (Reverse /@ protectionRules)
          ];

(* returns list rhs of rules where the lhs symbols are replaced by unique symbols *)
RemoveFixedParameters[rules_List] :=
    Module[{fixedParameters, freeParameters},
           fixedParameters = GetFixedParameter /@ rules;
           freeParameters  = GetFreeParameter /@ rules;
           TransformSymbolsIn[freeParameters, fixedParameters, Unique[CConversion`ToValidCSymbolString[#]]&]
          ];

ReplaceFixedParametersBySymbolsInTarget[Rule[a_,b_], fixedParameters_List] :=
    Rule[a,
         TransformSymbolsIn[b, fixedParameters, CConversion`ToValidCSymbol[#]&] /.
         { Sign[p_]  :> Global`LOCALINPUT[CConversion`ToValidCSymbol[Sign[p]]],
           FlexibleSUSY`Phase[p_] :> Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Phase[p]]] }
        ];

ReplaceFixedParametersBySymbolsInTarget[rules_List, fixedParameters_List] :=
    ReplaceFixedParametersBySymbolsInTarget[#, fixedParameters]& /@ rules;

ReplaceFixedParametersBySymbolsInTarget[solution_List] :=
    ReplaceFixedParametersBySymbolsInTarget[solution, GetFixedParameter /@ solution];

CreateMemberTreeLevelEwsbSolver[solution_List, substitutions_List:{}] :=
    Module[{result = "", body = "", fixedPars,
            i, par, expr, parStr, oldParStr, reducedSolution,
            type},
           fixedPars = GetFixedParameter /@ solution;
           reducedSolution = solution;
           If[reducedSolution =!= {},
              (* create local const refs to input parameters appearing
                 in the solution *)
              reducedSolution = reducedSolution /. {
                  Sign[p_]               :> Global`LOCALINPUT[CConversion`ToValidCSymbol[Sign[p]]],
                  FlexibleSUSY`Phase[p_] :> Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Phase[p]]]
                                                   };
              result = Parameters`CreateLocalConstRefsForInputParameters[reducedSolution, "LOCALINPUT"] <> "\n";
              (* save old parameters *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = CConversion`CreateCType[CConversion`GetScalarElementType[Parameters`GetType[par]]];
                  parStr = CConversion`RValueToCFormString[par];
                  oldParStr = "old_" <> CConversion`ToValidCSymbolString[par];
                  result = result <>
                           "const " <> type <> " " <> oldParStr <> " = " <> parStr <> ";\n";
                 ];
              result = result <> "\n";
              (* write solution *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = CConversion`GetScalarElementType[Parameters`GetType[par]];
                  result = result <> Parameters`SetParameter[par, expr, type];
                 ];
              result = result <> "\n";
              (* check for errors *)
              result = result <> "const bool is_finite = ";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  parStr = CConversion`RValueToCFormString[par];
                  result = result <> "IsFinite(" <> parStr <> ")";
                  If[i != Length[reducedSolution],
                     result = result <> " && ";
                    ];
                 ];
              result = result <> ";\n\n";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  oldParStr = "old_" <> CConversion`ToValidCSymbolString[par];
                  body = body <> Parameters`SetParameter[par, oldParStr, None];
                 ];
              body = body <> "error = EWSB_solver::FAIL;\n";
              If[substitutions === {},
                 result = result <>
                          "if (!is_finite) {\n" <>
                          IndentText[body] <>
                          "}";,
                 result = result <>
                          "if (is_finite) {\n" <>
                          IndentText[WrapLines[SetModelParametersFromEWSB[fixedPars, substitutions]]] <>
                          "} else {\n" <>
                          IndentText[body] <>
                          "}";
                ];
              ,
              result = "error = solve_ewsb_tree_level();\n";
             ];
           Return[result];
          ];

CreateTreeLevelEwsbSolver[solution_List] :=
    Module[{result = "",
            i, par, expr, parStr, decls = "", reducedSolution,
            type},
           reducedSolution = solution;
           If[reducedSolution =!= {},
              (* create local const refs to input parameters appearing
                 in the solution *)
              reducedSolution = ReplaceFixedParametersBySymbolsInTarget[reducedSolution];
              result = Parameters`CreateLocalConstRefs[RemoveFixedParameters[reducedSolution]] <> "\n";
              (* save old parameters *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  type = CConversion`CreateCType[CConversion`GetScalarElementType[Parameters`GetType[par]]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> type <> " " <> parStr <> ";\n";
                 ];
              result = result <> "\n";
              (* write solution *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = CConversion`GetScalarElementType[Parameters`GetType[par]];
                  result = result <> CConversion`ToValidCSymbolString[par] <> " = " <>
                           CConversion`CastTo[CConversion`RValueToCFormString[expr], type] <> ";\n";
                 ];
              result = result <> "\n";
              ,
              result = "error = solve_iteratively_at(model, 0);\n";
             ];
           Return[result];
          ];

SetTreeLevelSolution[ewsbSolution_, substitutions_List:{}, struct_String:"model."] :=
    Module[{i, parametersFixedByEWSB, par, parStr, body = "", result = ""},
           If[ewsbSolution =!= {},
              parametersFixedByEWSB = #[[1]]& /@ ewsbSolution;
              result = result <> "const bool is_finite = ";
              For[i = 1, i <= Length[parametersFixedByEWSB], i++,
                  par    = parametersFixedByEWSB[[i]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> "IsFinite(" <> parStr <> ")";
                  If[i != Length[parametersFixedByEWSB],
                     result = result <> " && ";
                    ];
                 ];
              result = result <> ";\n\n";
              For[i = 1, i <= Length[parametersFixedByEWSB], i++,
                  par    = parametersFixedByEWSB[[i]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  body = body <> Parameters`SetParameter[par, parStr, struct, None];
                 ];
              result = result <>
                       "if (is_finite) {\n" <>
                       IndentText[body] <>
                       If[substitutions === {}, "",
                          IndentText[WrapLines[SetModelParametersFromEWSB[parametersFixedByEWSB, substitutions, struct]]]
                         ] <>
                       IndentText["model.get_problems().unflag_no_ewsb_tree_level();\n"] <>
                       "} else {\n" <>
                       IndentText["error = EWSB_solver::FAIL;\nmodel.get_problems().flag_no_ewsb_tree_level();\n"] <>
                       "}";
             ];
           result
          ];

CreateNewEWSBRootFinder[] :=
    "new Root_finder<number_of_ewsb_equations>(tadpole_stepper, number_of_iterations, precision, ";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`FPIRelative] :=
    "new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative<number_of_ewsb_equations> >(ewsb_stepper, number_of_iterations, fixed_point_iterator::Convergence_tester_relative<number_of_ewsb_equations>(precision))";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`FPIAbsolute] :=
    "new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_absolute<number_of_ewsb_equations> >(ewsb_stepper, number_of_iterations, fixed_point_iterator::Convergence_tester_absolute<number_of_ewsb_equations>(precision))";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`FPITadpole] :=
    "new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_tadpole<number_of_ewsb_equations> >(ewsb_stepper, number_of_iterations, fixed_point_iterator::Convergence_tester_tadpole<number_of_ewsb_equations>(precision, tadpole_stepper))";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLHybrid] :=
    CreateNewEWSBRootFinder[] <> "Root_finder<number_of_ewsb_equations>::GSLHybrid)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLHybridS] :=
    CreateNewEWSBRootFinder[] <> "Root_finder<number_of_ewsb_equations>::GSLHybridS)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLBroyden] :=
    CreateNewEWSBRootFinder[] <> "Root_finder<number_of_ewsb_equations>::GSLBroyden)";

CreateEWSBRootFinder[rootFinder_ /; rootFinder === FlexibleSUSY`GSLNewton] :=
    CreateNewEWSBRootFinder[] <> "Root_finder<number_of_ewsb_equations>::GSLNewton)";

CreateEWSBRootFinders[{}] :=
    Block[{},
          Print["Error: List of EWSB root finders must not be empty!"];
          Quit[1];
         ];

MakeUniquePtr[str_String, obj_String] :=
    "std::unique_ptr<" <> obj <> ">(" <> str <> ")";

CreateEWSBRootFinders[rootFinders_List] :=
    Utils`StringJoinWithSeparator[MakeUniquePtr[#,"EWSB_solver"]& /@ (CreateEWSBRootFinder /@ rootFinders), ",\n"];

ConvertToReal[par_] :=
    If[Parameters`IsRealParameter[par],
       CConversion`ToValidCSymbolString[par],
       "Abs(" <> CConversion`ToValidCSymbolString[par] <> ")"
      ];

FillArrayEntryWithParameter[arrayName_String, par_, idx_] :=
    arrayName <> "[" <> ToString[idx-1] <> "] = " <> ConvertToReal[par] <> ";\n";

FillArrayWithParameters[arrayName_String, parameters_List] :=
    Module[{result = "", i},
           For[i = 1, i <= Length[parameters], i++,
               result = result <> FillArrayEntryWithParameter[arrayName, parameters[[i]], i];
              ];
           result
          ];

DivideArrayEntryByParameter[arrayName_String, par_, idx_] :=
    arrayName <> "[" <> ToString[idx-1] <> "] /= " <> CConversion`ToValidCSymbolString[par] <> ";\n";

DivideTadpolesByVEV[arrayName_String, vevToTadpoleAssociation_List] :=
    Module[{result = "", i, vevs},
           vevs = #[[3]]& /@ vevToTadpoleAssociation;
           For[i = 1, i <= Length[vevs], i++,
               result = result <> DivideArrayEntryByParameter[arrayName, vevs[[i]], i];
              ];
           result
          ];

CreateEwsbSolverWithTadpoles[solution_List] :=
    Module[{result = "", i, par, expr, parStr, reducedSolution, rules, type},
           reducedSolution = solution /.
               FlexibleSUSY`tadpole[p_] :> CConversion`ReleaseHoldAt[HoldForm[FlexibleSUSY`tadpole[[p-1]]], {1,2}];
           If[reducedSolution =!= {},
              (* create local const refs to parameters appearing on RHS
                 in the solution *)
              reducedSolution = ReplaceFixedParametersBySymbolsInTarget[reducedSolution];
              result = Parameters`CreateLocalConstRefs[RemoveFixedParameters[reducedSolution]];
              (* define variables for new parameters *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = CConversion`CreateCType[CConversion`GetScalarElementType[Parameters`GetType[par]]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> type <> " " <> parStr <> ";\n";
                 ];
              result = result <> "\n";
              (* write solution *)
              For[i = 1, i <= Length[reducedSolution], i++,
                  par  = reducedSolution[[i,1]];
                  expr = reducedSolution[[i,2]];
                  type = CConversion`GetScalarElementType[Parameters`GetType[par]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> parStr <> " = " <>
                           CConversion`CastTo[CConversion`RValueToCFormString[expr],type] <> ";\n";
                 ];
              result = result <> "\n";
              (* check for errors *)
              result = result <> "const bool is_finite = ";
              For[i = 1, i <= Length[reducedSolution], i++,
                  par    = reducedSolution[[i,1]];
                  parStr = CConversion`ToValidCSymbolString[par];
                  result = result <> "IsFinite(" <> parStr <> ")";
                  If[i != Length[reducedSolution],
                     result = result <> " && ";
                    ];
                 ];
              result = result <> ";\n";
              ,
              result = "const bool is_finite = false;\n";
             ];
           Return[result];
          ];

GetEWSBParametersFromVector[parametersFixedByEWSB_List, freePhases_List,
                            vector_String] :=
    Module[{i, result = "", par, parStr, type},
           For[i = 1, i <= Length[parametersFixedByEWSB], i++,
               par = parametersFixedByEWSB[[i]];
               type = CConversion`CreateCType[CConversion`GetScalarElementType[Parameters`GetType[par]]];
               parStr = CConversion`ToValidCSymbolString[par];
               result = result <>
                        "const " <> type <> " " <> parStr <> " = " <>
                        GetValueWithPhase[par, vector, i-1, freePhases] <> ";\n";
              ];
           Return[result];
          ];

SetEWSBParametersFromLocalCopies[parameters_List, struct_String] :=
    Module[{result = ""},
           (result = result <> Parameters`SetParameter[#, CConversion`ToValidCSymbolString[#], struct])& /@ parameters;
           result
          ];

SetEWSBSolution[parametersFixedByEWSB_List, freePhases_List, func_String, class_String] :=
    GetEWSBParametersFromVector[parametersFixedByEWSB, freePhases, func] <>
    SetEWSBParametersFromLocalCopies[parametersFixedByEWSB, class];

CreateEWSBParametersInitializationComma[{}] := "";

CreateEWSBParametersInitializationComma[parameters_List] :=
    Utils`StringJoinWithSeparator[ConvertToReal /@ parameters, ", "];

CreateEWSBParametersInitializationList[parameters_List] :=
    "{" <> CreateEWSBParametersInitializationComma[parameters] <> "}";

SetEWSBParameter[par_, idx_, array_String] :=
    array <> "[" <> ToString[idx] <> "] = " <> ConvertToReal[par] <> ";\n";

CreateEWSBParametersInitialization[parameters_List, array_String] :=
    StringJoin[MapIndexed[SetEWSBParameter[#1,First[#2 - 1],array]&, parameters]];

(* @todo handle patterns here *)
SetModelParametersFromEWSB[parametersFixedByEWSB_List, substitutions_List, class_String] :=
    Module[{subs = substitutions, localPars, result = ""},
           subs = subs /. { RuleDelayed[Sign[p_] /; Parameters`IsInputParameter[Sign[p]],
                                        Global`LOCALINPUT[CConversion`ToValidCSymbol[Sign[p]]]],
                            RuleDelayed[FlexibleSUSY`Phase[p_] /; Parameters`IsInputParameter[FlexibleSUSY`Phase[p]],
                                        Global`LOCALINPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Phase[p]]]] };
           (result = result <> Parameters`SetParameter[#[[1]], #[[2]], class])& /@ subs;
           localPars = Parameters`FindAllParameters[#[[2]]& /@ subs, parametersFixedByEWSB];
           Parameters`CreateLocalConstRefs[localPars] <> result
          ];

(* @todo handle patterns here *)
ApplyEWSBSubstitutions[parametersFixedByEWSB_List, substitutions_List, class_String:"model."] :=
    Module[{pars, subs = substitutions, result = ""},
           subs = subs /. { RuleDelayed[Sign[p_] /; Parameters`IsInputParameter[Sign[p]],
                                        Global`INPUT[CConversion`ToValidCSymbol[Sign[p]]]],
                            RuleDelayed[FlexibleSUSY`Phase[p_] /; Parameters`IsInputParameter[FlexibleSUSY`Phase[p]],
                                        Global`INPUT[CConversion`ToValidCSymbol[FlexibleSUSY`Phase[p]]]] };
           (result = result <> Parameters`SetParameter[#[[1]], #[[2]], class])& /@ subs;
           pars = DeleteDuplicates[Parameters`FindAllParameters[#[[2]]& /@ subs]];
           pars = Select[pars, !MemberQ[parametersFixedByEWSB, #]&];
           Parameters`CreateLocalConstRefs[pars] <> result
          ];

SolveEWSBIgnoringFailures[loops_Integer] :=
    Module[{flagEWSB, unflagEWSB, warning, result,
            flagFunc = "flag_no_ewsb" <> If[loops == 0, "_tree_level", ""]},
           flagEWSB = "this->problems." <> flagFunc <> "();\n";
           unflagEWSB = "this->problems.un" <> flagFunc <> "();\n";
           body = "if (has_no_ewsb_flag) {\n" <> IndentText[flagEWSB]
                  <> "} else {\n" <> IndentText[unflagEWSB] <> "}\n";
           body = "[this, has_no_ewsb_flag] () {\n" <> IndentText[body] <> "}\n";
           result = "const bool has_no_ewsb_flag = problems.no_ewsb();\n";
           result = result <> "const auto save_ewsb_flag = make_raii_guard(\n"
                    <> IndentText[body] <> ");\n";
           result = result <> "problems.un" <> flagFunc <> "();\n";
           If[loops == 0,
              result = result <> "solve_ewsb_tree_level();\n";,
              result = result <> "const auto save_ewsb_loop_order = make_raii_save(ewsb_loop_order);\n"
                       <> "ewsb_loop_order = " <> ToString[loops] <> ";\n"
                       <> "solve_ewsb();\n";
             ];
           warning = "WARNING(\"solving EWSB at " <> ToString[loops]
                     <> "-loop order failed\");\n";
           warning = "if (problems.no_ewsb()) {\n" <> IndentText[warning] <> "}\n";
           warning = "#ifdef ENABLE_VERBOSE\n" <> IndentText[warning] <> "#endif";
           IndentText[result] <> warning
          ];

End[];

EndPackage[];
