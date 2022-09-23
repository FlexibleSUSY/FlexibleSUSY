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

BeginPackage["SemiAnalytic`", {"SARAH`", "CConversion`", "Constraint`", "EWSB`", "Parameters`",
                               "TextFormatting`", "Utils`", "WriteOut`"}];

SemiAnalyticSolution::usage="Head of a semi-analytic solution.
A semi-analytic solution to the RGEs has the structure
SemiAnalyticSolution[parameter, {basis}]";

CheckSemiAnalyticBoundaryConditions::usage="";
IsSemiAnalyticSetting::usage="";
IsBasisParameterSetting::usage="";
IsBoundaryConstraint::usage="Returns True if given constraint corresponds
to the scale at which the semi-analytic boundary values are input.";
IsSemiAnalyticConstraint::usage="Returns True if given constraint
corresponds to the scale at which the semi-analytic solutions are evaluated.";
SelectBoundaryConstraint::usage="";

SetSemiAnalyticParameters::usage="";
GetSemiAnalyticParameters::usage="";
GetBoundaryValueParameters::usage="";
GetBoundaryValueParameterSymbol::usage="Returns the symbol for a given boundary
value parameter.";

IsAllowedSemiAnalyticParameter::usage="";
IsSemiAnalyticParameter::usage="";

GetSemiAnalyticSolutions::usage="Constructs the semi-analytic
solutions implied by the given list of boundary conditions.";
ExpandSemiAnalyticSolutions::usage="Expands the given solutions
in terms of the coefficients.";
GetSemiAnalyticParameterSubstitutions::usage="";
CreateBoundaryValueParameters::usage="Creates new parameters
representing the boundary values.";
CreateCoefficientParameters::usage="Creates new parameters
representing the coefficients in the semi-analytic solutions.";

CreateSemiAnalyticSolutionsDefinitions::usage="";
CreateLocalBoundaryValuesDefinitions::usage="";
CreateSemiAnalyticCoefficientGetters::usage="";
ConstructTrialDatasets::usage="Returns a list of datasets of the form
{integer id, {pars}, {input values}}";
InitializeTrialInputValues::usage="";
CreateBasisEvaluators::usage="";
CreateLinearSystemSolvers::usage="";
CalculateCoefficients::usage="";
CreateSemiAnalyticCoefficientsCalculation::usage="";
CreateCoefficientsCalculations::usage="";

SetTreeLevelEWSBSolution::usage="";

ApplySemiAnalyticBoundaryConditions::usage="";
ReplacePreprocessorMacros::usage="";
GetSemiAnalyticEWSBSubstitutions::usage="";
EvaluateSemiAnalyticSolutions::usage="";
GetExtraBoundaryParametersToSave::usage="";
SaveBoundaryValueParameters::usage="";
SetBoundaryValueParametersFromModel::usage="";
GetModelCoefficients::usage="";
PrintModelCoefficients::usage="";

Begin["`Private`"];

allSemiAnalyticParameters = {};

GetName[SemiAnalyticSolution[name_, basis_List]] := name;

GetBasis[SemiAnalyticSolution[name_, basis_List]] := basis;

GetBoundaryValueParameters[solution_SemiAnalyticSolution] :=
    DeleteDuplicates[Flatten[Parameters`FindAllParameters[GetBasis[solution]]]];

GetBoundaryValueParameters[solutions_List] :=
    DeleteDuplicates[Flatten[(Parameters`FindAllParameters[GetBasis[#]])& /@ solutions]];

IsDimensionZero[par_] :=
    Module[{dimZeroPars},
           dimZeroPars = Parameters`GetParametersWithMassDimension[0];
           MemberQ[Parameters`StripIndices[#]& /@ dimZeroPars, Parameters`StripIndices[par]]
          ];

IsDimensionOne[par_] :=
    Module[{dimOnePars},
           dimOnePars = { SARAH`BetaTijk };
           If[SARAH`SupersymmetricModel,
              dimOnePars = Append[dimOnePars, SARAH`BetaMi];,
              dimOnePars = Append[dimOnePars, SARAH`BetaMuij];
             ];
           dimOnePars = (Parameters`StripIndices[#[[1]]])& /@ (Join @@ dimOnePars);
           MemberQ[dimOnePars, Parameters`StripIndices[par]]
          ];

IsDiracGauginoMass[par_] :=
    Module[{diracMasses = {}},
           If[SARAH`SupersymmetricModel,
              diracMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaDGi;
             ];
           MemberQ[diracMasses, Parameters`StripIndices[par]]
          ];

IsScalarMass[par_] :=
    Module[{scalarMasses},
           If[SARAH`SupersymmetricModel,
              scalarMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`Betam2ij;,
              scalarMasses = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
             ];
           MemberQ[scalarMasses, Parameters`StripIndices[par]]
          ];

IsSoftBilinear[par_] :=
    Module[{softBilinears = {}},
           If[SARAH`SupersymmetricModel,
              softBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
             ];
           MemberQ[softBilinears, Parameters`StripIndices[par]]
          ];

IsSoftLinear[par_] :=
    Module[{softLinears = {}},
           If[SARAH`SupersymmetricModel,
              softLinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaLSi;
             ];
           MemberQ[softLinears, Parameters`StripIndices[par]]
          ];

IsAllowedSemiAnalyticParameter[par_] :=
    Or[IsDimensionOne[par],
       IsDiracGauginoMass[par],
       IsScalarMass[par],
       IsSoftBilinear[par],
       IsSoftLinear[par]];

IsSemiAnalyticParameter[par_] := MemberQ[allSemiAnalyticParameters, Parameters`StripIndices[par]];

SetSemiAnalyticParameters[parameters_List] :=
    Module[{},
           allSemiAnalyticParameters = DeleteDuplicates[Select[parameters, (IsAllowedSemiAnalyticParameter[#])&]];
           allSemiAnalyticParameters = (Parameters`StripIndices[#])& /@ allSemiAnalyticParameters;
          ];

GetSemiAnalyticParameters[] := allSemiAnalyticParameters;

IsSemiAnalyticSetting[setting_] :=
    Intersection[Constraint`FindFixedParametersFromConstraint[{setting}],
                 allSemiAnalyticParameters] =!= {};

IsBasisParameterSetting[setting_, solutions_List] :=
    Module[{allBasisParameters},
           allBasisParameters = GetBoundaryValueParameters[solutions];
           allBasisParameters = DeleteCases[allBasisParameters,
                                            p_ /; (Parameters`IsModelParameter[p] && !IsAllowedSemiAnalyticParameter[p])];
           Intersection[Constraint`FindFixedParametersFromConstraint[{setting}],
                        allBasisParameters] =!= {}
          ];

RemoveUnusedSettings[constraints_List] := Select[constraints, IsSemiAnalyticSetting];

IsBoundaryConstraint[constraint_] :=
    Module[{sortedPars, fixedPars},
           sortedPars = Sort[allSemiAnalyticParameters];
           fixedPars = Sort[Intersection[sortedPars,
                                         Constraint`FindFixedParametersFromConstraint[constraint]]];
           fixedPars =!= {} && fixedPars === sortedPars
          ];

SelectBoundaryConstraint[constraints_List] :=
    Module[{validConstraints, result = {}},
           validConstraints = Select[constraints, IsBoundaryConstraint];
           If[Length[validConstraints] >= 1,
              result = validConstraints[[1]];
             ];
           RemoveUnusedSettings[result]
          ];

IsSemiAnalyticConstraint[settings_List] := MemberQ[settings, FlexibleSUSY`FSSolveEWSBFor[___]];

CheckSemiAnalyticBoundaryConditions[constraints_List] :=
    Module[{i, sortedPars, fixedPars, boundaryValueSettings, boundaryValuePars},
           sortedPars = Sort[allSemiAnalyticParameters];
           For[i = 1, i <= Length[constraints], i++,
               fixedPars = Sort[Intersection[sortedPars,
                                             Constraint`FindFixedParametersFromConstraint[constraints[[i]]]]];
               If[fixedPars =!= {} && fixedPars =!= sortedPars,
                  Print["Error: all semi-analytic parameters must be set at the same scale."];
                  Print["   The following parameters are not set: ", Complement[sortedPars, fixedPars]];
                  Return[False];
                 ];
              ];
           boundaryValueSettings = GetBoundaryValueSubstitutions[SelectBoundaryConstraint[constraints]];
           boundaryValuePars = Select[Parameters`FindAllParameters[#[[2]]& /@ boundaryValueSettings],
                                      !IsDimensionZero[#]&];
           And @@ (PolynomialQ[#[[2]], boundaryValuePars]& /@ boundaryValueSettings)
          ];

SelectParametersWithMassDimension[parameters_List, dim_?IntegerQ] :=
    Module[{allParameters},
           allParameters = Parameters`GetParametersWithMassDimension[dim];
           Select[parameters, MemberQ[allParameters, #]&]
          ];

GetParameterExpansionRules[par_] :=
    Module[{dims, numDims, i, indices, rules = {}},
           dims = Parameters`GetParameterDimensions[par];
           If[dims =!= {1},
              numDims = Length[dims];
              indices = Table["i" <> ToString[i], {i, 1, numDims}];
              (* do not expand when indices are explicitly given *)
              rules = { RuleDelayed @@ Rule[par[Sequence @@ (Pattern[Evaluate[Symbol[#]], Blank[]]& /@ indices)],
                                            par[Sequence @@ (Symbol[#] & /@ indices)]],
                        RuleDelayed @@ Rule[par, Table[par[Sequence @@ (Symbol[#] & /@ indices)],
                                                       Evaluate[Sequence @@ (Inner[Function[{idx, bound}, {Symbol[idx], bound}],
                                                                                   indices, dims, List])]]]
                       };
             ];
           rules
          ];

GetMacroExpansionRules[] :=
    { RuleDelayed[CConversion`UNITMATRIX[n_Integer], IdentityMatrix[n]],
      RuleDelayed[CConversion`UNITMATRIXCOMPLEX[n_Integer], IdentityMatrix[n]],
      RuleDelayed[CConversion`ZEROARRAY[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROARRAYCOMPLEX[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROVECTOR[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROVECTORCOMPLEX[n_Integer], Table[0, {n}]],
      RuleDelayed[CConversion`ZEROMATRIX[m_Integer, n_Integer], Table[0, {m}, {n}]],
      RuleDelayed[CConversion`ZEROMATRIXCOMPLEX[m_Integer, n_Integer], Table[0, {m}, {n}]],
      RuleDelayed[CConversion`ZEROTENSOR3[l_Integer, m_Integer, n_Integer], Table[0, {l}, {m}, {n}]],
      RuleDelayed[CConversion`ZEROTENSOR3COMPLEX[l_Integer, m_Integer, n_Integer], Table[0, {l}, {m}, {n}]],
      RuleDelayed[CConversion`ZEROTENSOR4[l_Integer, m_Integer, n_Integer, o_Integer], Table[0, {l}, {m}, {n}, {o}]],
      RuleDelayed[CConversion`ZEROTENSOR4COMPLEX[l_Integer, m_Integer, n_Integer, o_Integer], Table[0, {l}, {m}, {n}, {o}]],
      RuleDelayed[CConversion`ZEROTENSOR[dimensions__ /; And @@ (IntegerQ /@ {dimensions})],
                  Table[0, Evaluate[Sequence @@ Map[List, List[dimensions]]]]],
      RuleDelayed[CConversion`ZEROTENSORCOMPLEX[dimensions__ /; And @@ (IntegerQ /@ {dimensions})],
                  Table[0, Evaluate[Sequence @@ Map[List, List[dimensions]]]]]
    };

GetBoundaryValueFromSetting[par_, value_] :=
    Module[{dims, numDims, i, indices,
            parameters, parameterExpansions, lhs,
            rhs, result},
           dims = Parameters`GetParameterDimensions[par];
           parameters = Parameters`FindAllParameters[value];
           parameterExpansions = Flatten[GetParameterExpansionRules /@ parameters, 1];
           rhs = value /. Join[GetMacroExpansionRules[], parameterExpansions];
           lhs = par /. GetParameterExpansionRules[par];
           (* construct rules *)
           If[dims =!= {1},
              (* check dimensions are consistent *)
              If[Dimensions[rhs] =!= dims,
                 Print["Error: dimensions do not match in setting: " {par, value}];
                 Quit[1];
                ];
              numDims = Length[dims];
              indices = Table["i" <> ToString[i], {i, 1, numDims}];
              result = Flatten[Table[Rule[Part[lhs,Sequence @@ (Symbol[#]& /@ indices)],
                                          Part[rhs, Sequence @@ (Symbol[#]& /@ indices)]],
                                     Evaluate[Sequence @@ (Inner[Function[{idx, bound}, {Symbol[idx], bound}],
                                                                 indices, dims, List])]]];,
              (* @todo distinguish length 1 vectors from scalars? *)
              result = {Rule[lhs, rhs]};
             ];
           result
          ];

IsImplicitConstraint[FlexibleSUSY`FSMinimize[__]] := True;
IsImplicitConstraint[FlexibleSUSY`FSFindRoot[__]] := True;
IsImplicitConstraint[FlexibleSUSY`FSSolveEWSBFor[__]] := True;
IsImplicitConstraint[setting_] := False;

GetPlaceholders[FlexibleSUSY`FSMinimize[parameters_List, value_]] :=
    Symbol[CConversion`ToValidCSymbolString[#] <> "MinSol"]& /@ Select[parameters, IsSemiAnalyticParameter];
GetPlaceholders[FlexibleSUSY`FSFindRoot[parameters_List, value_]] :=
    Symbol[CConversion`ToValidCSymbolString[#] <> "RootSol"]& /@ Select[parameters, IsSemiAnalyticParameter];
GetPlaceholders[FlexibleSUSY`FSSolveEWSBFor[parameters_List]] :=
    Symbol[CConversion`ToValidCSymbolString[#] <> "EWSBSol"]& /@ Select[parameters, IsSemiAnalyticParameter];

(* Returns a list of rules of the form {parameter -> value, ...}
   determined from the provided constraint settings.

   The list of settings are assumed to apply at a single scale.
   The value in each rule is the final result after applying all
   of the settings in the list, in order.
*)
GetBoundaryValueSubstitutions[settings_List] :=
    Module[{i, noTemp, test, fixedPars, boundaryValues, rules, finalRules = {}},
           noTemp = DeleteCases[settings, {FlexibleSUSY`FSTemporary[_], _}];
           test = Function[{first, second}, first[[1]] === second[[1]]];
           For[i = 1, i <= Length[noTemp], i++,
               Switch[noTemp[[i]],
                      {_, _},
                      rules = GetBoundaryValueFromSetting[noTemp[[i,1]], noTemp[[i,2]]];
                      finalRules = Last[(finalRules = Utils`AppendOrReplaceInList[finalRules, #, test])& /@ rules];,
                      (* @todo intermediate boundary values can be added to a list of extra parameters *)
                      _?IsImplicitConstraint,
                      fixedPars = Select[Constraint`FindFixedParametersFromConstraint[{noTemp[[i]]}], IsSemiAnalyticParameter];
                      finalRules = Last[(finalRules = Utils`AppendOrReplaceInList[finalRules, Rule[#[[1]], #[[2]]], test])&
                                        /@ Utils`Zip[fixedPars, GetPlaceholders[noTemp[[i]]]]];
                     ];
              ];
           finalRules
          ];

SumOverElements[par_] :=
    Module[{i, dims, numDims, indices, result},
           dims = Parameters`GetParameterDimensions[par];
           If[dims === {1},
              result = par;,
              numDims = Length[dims];
              indices = Table["i" <> ToString[i], {i, 1, numDims}];
              result = Sum[par[Evaluate[Sequence @@ (Symbol[#]& /@ indices)]],
                           Evaluate[Sequence @@ (Inner[Function[{idx, bound}, {Symbol[idx], bound}],
                                                       indices, dims, List])]];
             ];
           result
          ];

(* @todo handle cases like TYu = Yu *)
GetCoefficientAndTerm[monomial_NumericQ, dimZeroPars_List] := {monomial, Null};

GetCoefficientAndTerm[monomial_Times, dimZeroPars_List] :=
    Module[{i, operands, coeff = {}, term = {}},
           operands = List @@ monomial;
           For[i = 1, i <= Length[operands], i++,
               If[NumericQ[operands[[i]]] || DependsAtMostOn[operands[[i]], dimZeroPars],
                  coeff = Append[coeff, operands[[i]]];,
                  term = Append[term, operands[[i]]];
                 ];
              ];
           Which[term === {} && coeff =!= {},
                 coeff = Times @@ coeff;
                 term = Null;,
                 term =!= {} && coeff === {},
                 coeff = 1;
                 term = Times @@ term;,
                 True,
                 coeff = Times @@ coeff;
                 term = Times @@ term;
                ];
           {coeff, term}
          ];

GetCoefficientAndTerm[monomial_, dimZeroPars_List] :=
    Module[{coeff, term},
           If[DependsAtMostOn[monomial, dimZeroPars],
              coeff = monomial;
              term = Null;,
              coeff = 1;
              term = monomial;
             ];
           {coeff, term}
          ];

GetSolutionBasis[pars_List, subs_List, addConjugateTerms_:False] :=
    Module[{i, dimZeroPars, dimensions, expr,
            monomials = {}, coeff, term, basis = {}},
           If[pars =!= {},
              dimZeroPars = Parameters`GetParametersWithMassDimension[0];
              dimensions = Parameters`GetParameterDimensions[#]& /@ pars;
              expr = Plus @@ (SumOverElements[#]& /@ pars);
              expr = Expand[expr /. subs];
              (* get independent basis elements from result *)
              If[Head[expr] === Plus,
                 monomials = List @@ expr;,
                 monomials = {expr};
                ];
              For[i = 1, i <= Length[monomials], i++,
                  {coeff, term} = GetCoefficientAndTerm[monomials[[i]], dimZeroPars];
                  If[term =!= Null,
                     basis = Append[basis, term];,
                     Print["Error: cannot interpret settings for parameters ", InputForm[pars]];
                     Print["   because a boundary value cannot be determined."];
                     Quit[1];
                    ];
                 ];
              If[addConjugateTerms,
                 basis = Join[basis, Conjugate /@ Select[basis, !Parameters`IsRealExpression[#]&]];
                ];
             ];
           DeleteDuplicates[basis]
          ];

GetComplexConjugateBasis[SemiAnalyticSolution[par_, basis_List]] :=
    SemiAnalyticSolution[par, (If[Parameters`IsRealExpression[#], #, Conjugate[#]])& /@ basis];

ListAllTermsOfForm[termTypes_List] :=
    Module[{i, sets = {}, possibleTerms, terms = {}},
           For[i = 1, i <= Length[termTypes], i++,
               sets = termTypes[[i]] /. {f_, n_Integer?Positive} :> Sequence @@ ConstantArray[f, n];
               sets = With[{set = #}, {#[[1]], #[[2]]}& /@ set]& /@ sets;
               possibleTerms = Tuples[sets];
               possibleTerms = Flatten[With[{tuple = #}, Tuples[#[[2]]& /@ tuple]]& /@ possibleTerms, 1];
               possibleTerms = DeleteDuplicates[Apply[Times, possibleTerms, 1]];
               terms = Join[terms, possibleTerms];
              ];
           DeleteDuplicates[terms]
          ];

GetLinearSystemSolutions[pars_List, subs_List, nonHomogeneousTerms_List:{}, addConjugateTerms_:False] :=
    Module[{basis, solns = {}},
           basis = GetSolutionBasis[pars, subs, addConjugateTerms];
           solns = SemiAnalyticSolution[#, basis]& /@ pars;
           If[nonHomogeneousTerms =!= {},
              solns = SemiAnalyticSolution[#[[1]],
                                           DeleteDuplicates[Join[#[[2]], ListAllTermsOfForm[nonHomogeneousTerms]]]]& /@ solns;
             ];
           solns
          ];

GetBoundaryConditionsFor[pars_List, subs_List] := Select[subs, MemberQ[pars, Parameters`StripIndices[#[[1]]]]&];

GetScalarMassesSolutions[scalarMasses_List, boundaryConditionSubs_List, dimOneSolns_List] :=
    Module[{extraTerms},
           If[SARAH`SupersymmetricModel,
              extraTerms = {{{dimOneSolns, 1}, {GetComplexConjugateBasis /@ dimOneSolns, 1}}};,
              extraTerms = {{{dimOneSolns, 2}}, {{GetComplexConjugateBasis /@ dimOneSolns, 2}},
                            {{dimOneSolns, 1}, {GetComplexConjugateBasis /@ dimOneSolns, 1}}};
             ];
           GetLinearSystemSolutions[scalarMasses, boundaryConditionSubs, extraTerms]
          ];

(* dimensionful SUSY parameters have trivial semi-analytic solutions *)
GetSUSYSemiAnalyticSolutions[susyBetas_List] :=
    Module[{pars},
           pars = Parameters`StripIndices[#[[1]]]& /@ susyBetas;
           pars = Flatten[(# /. GetParameterExpansionRules[#])& /@ pars];
           SemiAnalyticSolution[#, pars]& /@ pars
          ];

GetSoftBilinearsSolutions[softBilinears_List, boundaryConditionSubs_List, dimOneSolns_List] :=
    Module[{susyBilinears, extraTerms},
           susyBilinears = GetSUSYSemiAnalyticSolutions[SARAH`BetaMuij];
           extraTerms = {{{dimOneSolns, 1}, {susyBilinears, 1}}};
           GetLinearSystemSolutions[softBilinears, boundaryConditionSubs, extraTerms]
          ];

GetSoftLinearsSolutions[softLinears_List, boundaryConditionSubs_List, dimOneSolns_List,
                        diracSolns_List, scalarMassesSolns_List, softBilinearsSolns_List] :=
    Module[{susyBilinears, susyLinears, extraTerms},
           susyBilinears = GetSUSYSemiAnalyticSolutions[SARAH`BetaMuij];
           susyLinears = GetSUSYSemiAnalyticSolutions[SARAH`BetaLi];
           extraTerms = {{{susyLinears, 1}, {dimOneSolns, 1}},
                         {{susyBilinears, 2}, {dimOneSolns, 1}},
                         {{susyBilinears, 1}, {softBilinearsSolns, 1}},
                         {{GetComplexConjugateBasis /@ susyBilinears, 1},
                          {dimOneSolns, 1}, {GetComplexConjugateBasis /@ dimOneSolns, 1}},
                         {{GetComplexConjugateBasis /@ susyBilinears, 1}, {scalarMassesSolns, 1}},
                         {{GetComplexConjugateBasis /@ softBilinearsSolns, 1}, {dimOneSolns, 1}}
                        };
           If[diracSolns =!= {},
              extraTerms = Join[extraTerms,
                                {{{diracSolns, 1}, {scalarMassesSolns, 1}},
                                 {{GetComplexConjugateBasis /@ diracSolns, 2}, {dimOneSolns, 1}},
                                 {{diracSolns, 2}, {susyBilinears, 1}}
                                }
                               ];
             ];
           GetLinearSystemSolutions[softLinears, boundaryConditionSubs, extraTerms]
          ];

GetSemiAnalyticSolutions[settings_List] :=
    Module[{boundaryValSubs = {}, dimOnePars = {}, dimOneBCs = {}, dimOneSolns = {},
            diracPars = {}, diracBCs = {}, diracSolns = {}, scalarMasses = {},
            scalarMassesBCs = {}, scalarMassesSolns = {}, softBilinears = {},
            softBilinearsBCs = {}, softBilinearsSolns = {}, softLinears = {},
            softLinearBCs = {}, softLinearSolns = {}, result = {}},

           boundaryValSubs = GetBoundaryValueSubstitutions[settings];

           dimOnePars = SelectParametersWithMassDimension[allSemiAnalyticParameters, 1];
           If[SARAH`SupersymmetricModel,
              (* do Dirac gaugino masses separately *)
              diracPars = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaDGi;
              diracPars = Select[dimOnePars, MemberQ[diracPars, #]&];
              diracBCs = GetBoundaryConditionsFor[diracPars, boundaryValSubs];
              dimOnePars = Complement[dimOnePars, diracPars];
             ];
           dimOneBCs = GetBoundaryConditionsFor[dimOnePars, boundaryValSubs];
           dimOneSolns = GetLinearSystemSolutions[dimOnePars, dimOneBCs, !SARAH`SupersymmetricModel];
           result = Join[result, dimOneSolns];

           If[diracPars =!= {},
              diracSolns = GetLinearSystemSolutions[diracPars, diracBCs];
              result = Join[result, diracSolns];
             ];

           scalarMasses = SelectParametersWithMassDimension[allSemiAnalyticParameters, 2];
           If[SARAH`SupersymmetricModel,
              softBilinears = Parameters`StripIndices[#[[1]]]& /@ SARAH`BetaBij;
              softBilinears = Select[scalarMasses, MemberQ[softBilinears, #]&];
              softBilinearsBCs = GetBoundaryConditionsFor[softBilinears, boundaryValSubs];
              scalarMasses = Complement[scalarMasses, softBilinears];
             ];
           scalarMassesBCs = GetBoundaryConditionsFor[scalarMasses, boundaryValSubs];

           scalarMassesSolns
               = GetScalarMassesSolutions[scalarMasses, scalarMassesBCs, dimOneSolns];
           result = Join[result, scalarMassesSolns];

           If[SARAH`SupersymmetricModel,
              softBilinearsSolns
                  = GetSoftBilinearsSolutions[softBilinears, softBilinearsBCs, dimOneSolns];
              result = Join[result, softBilinearsSolns];
             ];

           If[SARAH`SupersymmetricModel,
              softLinears = SelectParametersWithMassDimension[allSemiAnalyticParameters, 3];
              softLinearsBCs = GetBoundaryConditionsFor[softLinears, boundaryValSubs];
              softLinearsSolns
                  = GetSoftLinearsSolutions[softLinears, softLinearsBCs, dimOneSolns,
                                            diracSolns, scalarMassesSolns, softBilinearsSolns];
              result = Join[result, softLinearsSolns];
             ];

           (* convert to C index convention in basis parameters *)
           (SemiAnalyticSolution[GetName[#], Parameters`DecreaseIndexLiterals[GetBasis[#]]])& /@ result
          ];

GetBoundaryValueParameterSymbol[parameter_?Parameters`IsInputParameter] := parameter;

GetBoundaryValueParameterSymbol[parameter_?Parameters`IsExtraParameter] := parameter;

GetBoundaryValueParameterSymbol[parameter_?Parameters`IsModelParameter] :=
    Symbol[CConversion`ToValidCSymbolString[parameter] <> "BV"];

GetBoundaryValueParameterSymbol[parameter_?Parameters`IsOutputParameter] :=
    Symbol[CConversion`ToValidCSymbolString[parameter] <> "BV"];

GetBoundaryValueParameterSymbol[parameter_?Parameters`IsPhase] :=
    Symbol[CConversion`ToValidCSymbolString[parameter] <> "BV"];

CreateBoundaryValueParameters[solutions_List] :=
    {GetBoundaryValueParameterSymbol[#], {}, Parameters`GetType[#]}& /@ (GetBoundaryValueParameters[solutions]);

CreateBoundaryValueParameterName[par_] :=
    CConversion`ToValidCSymbolString[GetBoundaryValueParameterSymbol[par]];

CreateCoefficients[SemiAnalyticSolution[par_, basis_]] :=
    Module[{i},
           Table[Symbol[CConversion`ToValidCSymbolString[par] <> "Coeff" <> ToString[i]], {i, 1, Length[basis]}]
          ];

CreateCoefficientParameters[solution_SemiAnalyticSolution] :=
    {#, {}, Parameters`GetType[GetName[solution]]}& /@ CreateCoefficients[solution];

CreateCoefficientParameters[solutions_List] :=
    Join[Sequence @@ (CreateCoefficientParameters /@ solutions)];

CreateCoefficientNames[solution_SemiAnalyticSolution] :=
    CConversion`ToValidCSymbolString /@ CreateCoefficients[solution];

CreateSemiAnalyticSolutionsDefinitions[solution_SemiAnalyticSolution] :=
    Module[{coeffs, defs = ""},
           coeffs = CreateCoefficientParameters[solution];
           (defs = defs <> Parameters`CreateParameterDefinitionAndDefaultInitialize[#])& /@ coeffs;
           Return[defs];
          ];

CreateSemiAnalyticSolutionsDefinitions[solutions_List] :=
    Module[{def = ""},
           (def = def <> CreateSemiAnalyticSolutionsDefinitions[#])& /@ solutions;
           Return[def];
          ];

CreateSemiAnalyticCoefficientGetters[solution_SemiAnalyticSolution] :=
    Module[{getters = "", coeffs},
           coeffs = CreateCoefficientParameters[solution];
           getters = (CConversion`CreateInlineGetters[CConversion`ToValidCSymbolString[#[[1]]],
                                                      CConversion`ToValidCSymbolString[#[[1]]],
                                                      #[[3]]])& /@ coeffs;
           StringJoin[getters]
          ];

CreateSemiAnalyticCoefficientGetters[solutions_List] :=
    Module[{getter = ""},
           (getter = getter <> CreateSemiAnalyticCoefficientGetters[#])& /@ solutions;
           Return[getter];
          ];

GetBoundaryValueFromStruct[parameter_, struct_String] :=
    Module[{parStr, type, name, body = "", function},
           type = Parameters`GetType[parameter];
           parStr = CreateBoundaryValueParameterName[parameter];
           name = CConversion`CreateGetterReturnType[type]
                  <> " get_" <> parStr <> "() const {";
           body = "return " <> struct <> "get_" <> parStr <> "();";
           name <> " " <> body <> " }\n"
          ];

SetBoundaryValueInStruct[parameter_, struct_String] :=
    Module[{parStr, type, name, body = "", function},
           type = Parameters`GetType[parameter];
           parStr = CreateBoundaryValueParameterName[parameter];
           name = "void set_" <> parStr <> "("
                  <> CConversion`CreateSetterInputType[type]
                  <> " " <> parStr <> "_) {";
           body = struct <> "set_" <> parStr <> "(" <> parStr <> "_);";
           name <> " " <> body <> " }\n"
          ];

GetCoefficientsFromStruct[solution_, struct_String] :=
    Module[{i, parameter, type, coeffStrs, name, body = "", functions = ""},
           parameter = GetName[solution];
           type = Parameters`GetType[parameter];
           coeffStrs = CreateCoefficientNames[solution];
           For[i = 1, i <= Length[coeffStrs], i++,
               name = CConversion`CreateGetterReturnType[type]
                      <> " get_" <> coeffStrs[[i]] <> "() const {";
               body = "return " <> struct <> "get_" <> coeffStrs[[i]] <> "();";
               functions = functions <> name <> " " <> body <> " }\n";
              ];
           Return[functions];
          ];

GetModelCoefficients[solutions_List, struct_String:"solutions."] :=
    Module[{result},
           result = GetCoefficientsFromStruct[#, struct]& /@ solutions;
           Utils`StringJoinWithSeparator[result, "\n"]
          ];

PrintModelCoefficients[solutions_List, streamName_String, struct_String:"solutions."] :=
    Module[{coeffs, result = ""},
           coeffs = CreateCoefficientParameters[solutions];
           result = streamName <> " << \"input_scale = \" << " <> struct <> "get_input_scale() << '\\n';\n";
           result = result <> streamName <> " << \"output_scale = \" << " <> struct
                    <> "get_output_scale() << '\\n';\n";
           result <> WriteOut`PrintExtraParameters[{#[[1]], #[[3]]}& /@ coeffs, streamName, Global`COEFFICIENT]
          ];

CreateLocalBoundaryValuesDefinitions[solutions_List] :=
    Module[{boundaryValues, defns},
           boundaryValues = {#, Parameters`GetType[#]}& /@ GetBoundaryValueParameters[solutions];
           defns = (Parameters`CreateParameterDefinitionAndDefaultInitialize[#])& /@ boundaryValues;
           StringJoin[defns] <> "\n"
          ];

ApplySettingLocally[{parameter_, value_}, modelPrefix_String] :=
    Which[Parameters`IsModelParameter[parameter],
          Parameters`SetParameter[parameter, value, modelPrefix],
          Parameters`IsInputParameter[parameter],
          Parameters`SetInputParameter[parameter, value, "INPUTPARAMETER"],
          Parameters`IsPhase[parameter],
          Parameters`SetPhase[parameter, value, modelPrefix],
          True,
          Print["Error: ", parameter, " is neither a model nor an input parameter!"];
          ""
         ];

ReplaceImplicitConstraints[settings_List] :=
    Module[{macroPosns, macros, fixedPars, values, replacements},
           macroPosns = Position[settings, FlexibleSUSY`FSMinimize[__] | \
                                           FlexibleSUSY`FSFindRoot[__] | \
                                           FlexibleSUSY`FSSolveEWSBFor[__]];
           macros = Extract[settings, macroPosns];
           fixedPars = Select[Constraint`FindFixedParametersFromConstraint[{#}],
                              IsSemiAnalyticParameter]& /@ macros;
           values = GetPlaceholders[#]& /@ macros;
           replacements = MapIndexed[(Utils`Zip[#1, values[[#2[[1]]]]])&, fixedPars];
           replacements = MapThread[Rule, {macroPosns, replacements}];
           ReplacePart[settings, Rule[#[[1]], Sequence @@ #[[2]]] & /@ replacements]
          ];

ApplySemiAnalyticBoundaryConditions[settings_List, solutions_List, modelPrefix_String:"model."] :=
    Module[{noMacros, boundaryValues, parameters, i,
            setBoundaryValues = "", result = ""},
           (* replace implicit constraints with placeholder values *)
           noMacros = ReplaceImplicitConstraints[settings];
           (* @todo handle temporary settings properly *)
           noMacros = DeleteCases[noMacros, {FlexibleSUSY`FSTemporary[_], _}];
           boundaryValues = GetBoundaryValueParameters[solutions];
           parameters = Select[Parameters`FindAllParameters[#[[2]]& /@ noMacros], !MemberQ[boundaryValues, #]&];
           (* in SUSY models, boundary values for the dimensionful SUSY parameters should also be set *)
           If[SARAH`SupersymmetricModel,
              noMacros = Join[noMacros, {#, #}& /@ (Select[boundaryValues,
                                                           (Parameters`IsModelParameter[#] && !IsAllowedSemiAnalyticParameter[#])&])];
             ];
           setBoundaryValues = ("const auto " <> CConversion`ToValidCSymbolString[#]
                                <> " = BOUNDARYVALUE(" <> CConversion`ToValidCSymbolString[#] <> ");\n")& /@ boundaryValues;
           (result = result <> ApplySettingLocally[#, modelPrefix])& /@ noMacros;
           Parameters`CreateLocalConstRefs[parameters] <> StringJoin[setBoundaryValues] <> result
          ];

ReplacePreprocessorMacros[expr_String, solutions_List] :=
    Module[{boundaryValues, coeffs, semiAnalyticPars},
           boundaryValues = CreateBoundaryValue /@ (GetBoundaryValueParameters[solutions]);
           coeffs = Flatten[CreateCoefficients /@ solutions];
           semiAnalyticPars = Join[boundaryValues, coeffs];
           macroRules = Flatten[{ Rule["EXTRAPARAMETER(" <> CConversion`ToValidCSymbolString[#] <> ")",
                                       "SEMIANALYTICPARAMETER(" <> CConversion`ToValidCSymbolString[#] <> ")"],
                                  Rule["MODELPARAMETER(" <> CConversion`ToValidCSymbolString[#] <> ")",
                                       "SEMIANALYTICPARAMETER(" <> CConversion`ToValidCSymbolString[#] <> ")"]
                                }& /@ semiAnalyticPars];
           StringReplace[expr, macroRules]
          ];

GetSubstitutionsWithIndices[{parameter_, replacement_}, basisPars_, coeffs_] :=
    Module[{type, numIndices, indices, coeffRules, indexedSubs, substitutions},
           type = Parameters`GetType[parameter];
           Which[MatchQ[type, CConversion`ScalarType[_]],
                 numIndices = 0,
                 MatchQ[type, CConversion`VectorType[_, _]],
                 numIndices = 1,
                 MatchQ[type, CConversion`ArrayType[_, _]],
                 numIndices = 1,
                 MatchQ[type, CConversion`MatrixType[_, _, _]],
                 numIndices = 2,
                 MatchQ[type, CConversion`TensorType[__]],
                 numIndices = Length[Rest[type]],
                 True,
                 Print["Error: unrecognized type: ", type];
                 Quit[1];
                ];
           substitutions = {{parameter, replacement}};
           If[numIndices =!= 0,
              indices = Table["i" <> ToString[i], {i, 1, numIndices}];
              coeffRules = RuleDelayed[#, #[Sequence @@ (Symbol[#] & /@ indices)]]& /@ coeffs;
              indexedSubs = {parameter[Sequence @@ (Pattern[Evaluate[Symbol[#]], Blank[]]& /@ indices)],
                             replacement /. coeffRules};
              substitutions = Append[substitutions, indexedSubs];
             ];
           substitutions
          ];

GetSemiAnalyticEWSBSubstitutions[solution_SemiAnalyticSolution] :=
    Module[{parameter, basisRules, basisPars, coeffs, replacement, result},
           parameter = GetName[solution];
           dim = Parameters`GetParameterDimensions[parameter];
           basisRules = Rule[#, GetBoundaryValueParameterSymbol[#]]& /@ (GetBoundaryValueParameters[solution]);
           basisPars = #[[2]]& /@ basisRules;
           coeffs = CreateCoefficients[solution];
           replacement = Dot[coeffs, GetBasis[solution]] /. basisRules;
           result = GetSubstitutionsWithIndices[{parameter, replacement}, basisPars, coeffs];
           (Rule @@ #)& /@ result
          ];

GetSemiAnalyticEWSBSubstitutions[solutions_List] :=
    Join[Sequence @@ (GetSemiAnalyticEWSBSubstitutions /@ solutions)]

GetDefaultSettings[parameters_List] :=
    Module[{generateIndices, defaultSettings},
           generateIndices[p_] :=
               Module[{dims, i},
                      dims = Parameters`GetParameterDimensions[p];
                      If[dims === {1},
                         {p},
                         ((p @@ #) & /@ Tuples[Table[i, {i, 0, #-1}] & /@ GetParameterDimensions[p]])
                        ]
                     ];
           {#, 0}& /@ (Flatten[generateIndices /@ parameters])
          ];

ReplaceIndicesFromMonomial[pars_List, monomial_] :=
    Module[{i, par, indexedParCases, parsWithIndices},
           parsWithIndices = Last[Reap[
               For[i = 1, i <= Length[pars], i++,
                   par = pars[[i]];
                   indexedParCases = Cases[monomial, p_[__] /; p === par, {0, Infinity}];
                   If[Parameters`GetParameterDimensions[par] =!= {1} &&
                      indexedParCases =!= {},
                      Sow[indexedParCases],
                      Sow[par]
                     ];
                  ];
               ]];
  Flatten[parsWithIndices]
  ];

GetParametersIncludingIndices[monomial_] :=
    Module[{pars},
           pars = Parameters`FindAllParameters[monomial];
           ReplaceIndicesFromMonomial[pars, monomial]
          ];

HaveEquivalentZeroSet[monomialOne_, monomialTwo_] :=
    Module[{monomialOnePars, monomialTwoPars},
           monomialOnePars = GetParametersIncludingIndices[monomialOne];
           monomialTwoPars = GetParametersIncludingIndices[monomialTwo];
           Sort[monomialOnePars] === Sort[monomialTwoPars]
          ];

RemoveConjugates[expr_] := expr //. (Conjugate | SARAH`Conj | Susyno`LieGroups`conj)[p_] :> p;

GetTrialInputsForSubset[monomials_List, defaultSettings_List] :=
    Module[{i, j, nEq, pars, equivalentRealTerms, realValues, complexValues, trialValues, evaluateSettings},
           pars = DeleteDuplicates[Flatten[GetParametersIncludingIndices /@ monomials]];
           incr = 1 / Length[pars];
           equivalentRealTerms = Gather[monomials, RemoveConjugates[#1] === RemoveConjugates[#2]&];
           trialValues = Reap[
                              For[i = 1, i <= Length[equivalentRealTerms], i++,
                                  nEq = Length[equivalentRealTerms[[i]]];
                                  realValues = MapIndexed[{#1, i + (First[#2] - 1) incr}&, pars];
                                  complexValues = Table[MapIndexed[{#1, i + (First[#2] - 1) incr +
                                                                        If[!Parameters`IsRealParameter[#],
                                                                           \[ImaginaryI] (j - 1 + (First[#2] - 1) incr),
                                                                           0
                                                                          ]}&, pars], {j, 2, nEq}];
                                  Sow[Join[{realValues}, complexValues]];
                                 ];
                             ];
           trialValues = Flatten[Last[trialValues], 2];
           evaluateSettings[values_] :=
               Module[{settings},
                      settings = Rule[{#[[1]], _}, {#[[1]], #[[2]]}] & /@ values;
                      defaultSettings /. settings
                     ];
           evaluateSettings /@ trialValues
          ];

(* @note assumes the solution is a polynomial in the boundary value parameters *)
GetRequiredBasisPoints[solution_SemiAnalyticSolution, defaultSettings_List] :=
    Module[{par, basis, trialValueRules = {}, independentSets, inputs},
           par = GetName[solution];
           basis = GetBasis[solution];
           independentSets = Gather[basis, HaveEquivalentZeroSet];
           inputs = Flatten[Join[GetTrialInputsForSubset[#, defaultSettings]& /@ independentSets], 1];
           CheckLinearlyIndependentInputs[solution, inputs];
           {par, inputs}
          ];

IsLinearlyDependentSet[vectors_List] := MatrixRank[vectors] < Length[vectors];

CheckLinearlyIndependentInputs[SemiAnalyticSolution[par_, basis_], inputs_List] :=
    Module[{inputsRules, basisValues},
           inputsRules = (Rule[#[[1]], #[[2]]]& /@ #)& /@ inputs;
           basisValues = (basis /. #)& /@ inputsRules;
           If[IsLinearlyDependentSet[basisValues],
              Print["Error: a linearly independent set of trial values could not be found!"];
              Print["   The semi-analytic coefficients for ", par,
                    " cannot be determined."];
              Quit[1];
             ];
          ];

AreEquivalentInputSets[inputSetOne_List, inputSetTwo_List] :=
    Module[{sortedSetOne, sortedSetTwo},
           sortedSetOne = Sort[inputSetOne];
           sortedSetTwo = Sort[inputSetTwo];
           sortedSetOne === sortedSetTwo
          ];

RequireSameInput[{parOne_, basisOne_List}, {parTwo_, basisTwo_List}] :=
    Module[{parOneType, parTwoType, sortedBasisOne, sortedBasisTwo},
           parOneType = CConversion`GetScalarElementType[Parameters`GetType[parOne]];
           parTwoType = CConversion`GetScalarElementType[Parameters`GetType[parTwo]];
           (parOneType === parTwoType) && AreEquivalentInputSets[basisOne, basisTwo]
          ];

ConstructTrialDatasets[solutions_List, trialValues_List:{}] :=
    Module[{boundaryValues, defaultSettings, requiredBases, datasets},
           boundaryValues = GetBoundaryValueParameters[solutions];
           defaultSettings = GetDefaultSettings[boundaryValues];
           requiredBases = GetRequiredBasisPoints[#, defaultSettings]& /@ solutions;
           datasets = Gather[requiredBases, RequireSameInput];
           MapIndexed[{First[#2], (#[[1]]) & /@ #1, First[(#[[2]]) & /@ #1]} &, datasets]
          ];

InitializeTrialInput[index_, basisValues_List, keys_List, struct_String:"trial_data"] :=
    Module[{result = ""},
           (result = result <> struct <> "[" <> ToString[index-1] <> "].boundary_values."
                    <> CConversion`ToValidCSymbolString[Parameters`StripIndices[#[[1]]]]
                    <> Parameters`CreateIndices[#[[1]]]
                    <> " = " <> CConversion`RValueToCFormString[#[[2]]]
                    <> ";\n")& /@ basisValues;
           (result = result <> struct <> "[" <> ToString[index-1] <> "].basis_sets.push_back("
                     <> ToString[#] <> ");\n")& /@ keys;
           result
          ];

CollectDatasetIndices[equivalentInputs_List] :=
    Module[{inputValues, indices},
           inputValues = First[equivalentInputs][[1]];
           indices = #[[2]] & /@ equivalentInputs;
           {inputValues, indices}
          ];

InitializeTrialInputValues[datasets_List] :=
    Module[{distinctInputs, numPoints, initialization = ""},
           distinctInputs = Flatten[With[{index = #[[1]], inputsList = #[[3]]}, {#, index} & /@
                                         inputsList] & /@ datasets, 1];
           distinctInputs = CollectDatasetIndices /@ Gather[distinctInputs, AreEquivalentInputSets[#1[[1]], #2[[1]]]&];
           initialization = MapIndexed[InitializeTrialInput[First[#2], #1[[1]], #1[[2]]]&, distinctInputs];
           {Length[distinctInputs], Utils`StringJoinWithSeparator[initialization, "\n"]}
          ];

CreateBasisEvaluator[name_String, parameters_List, basis_List] :=
    Module[{i, dim, boundaryValues, setBoundaryValues, returnType, body = "", result = ""},
           dim = Length[basis];
           boundaryValues = Parameters`FindAllParameters[basis];
           If[(And @@ (Parameters`IsRealParameter /@ parameters)) && (And @@ (Parameters`IsRealExpression /@ basis)),
              returnType = CConversion`MatrixType[CConversion`realScalarCType, 1, dim];,
              returnType = CConversion`MatrixType[CConversion`complexScalarCType, 1, dim];
             ];
           body = CConversion`CreateDefaultDefinition["result", returnType] <> ";\n";
           setBoundaryValues = ("const auto " <> CConversion`ToValidCSymbolString[#]
                                <> " = BOUNDARYVALUE(" <> CConversion`ToValidCSymbolString[#] <> ");\n")& /@ boundaryValues;
           body = body <> setBoundaryValues <> "\n";
           For[i = 1, i <= dim, i++,
               body = body <> "result(" <> ToString[i-1] <> ") = " <> CConversion`RValueToCFormString[basis[[i]]] <> ";\n";
              ];
           body = body <> "\nreturn result;\n";
           result = "auto " <> name <> " = [] (const Boundary_values& boundary_values) {\n";
           result <> IndentText[body] <> "};\n"
          ];

CreateBasisEvaluators[solutions_List] :=
    Module[{basisSets, bases, evaluators = ""},
           basisSets = Gather[{GetName[#], GetBasis[#]}& /@ solutions, RequireSameInput];
           bases = With[{equivalentBases = #}, {#[[1]]& /@ equivalentBases,
                                                Flatten[DeleteDuplicates[#[[2]]& /@ equivalentBases]]}]& /@ basisSets;
           evaluators = MapIndexed[CreateBasisEvaluator["basis_" <> ToString[First[#2]], #1[[1]], #1[[2]]]&, bases];
           Utils`StringJoinWithSeparator[evaluators, "\n"] <> "\n"
          ];

CreateLinearSystemSolver[dataset_, solutions_List] :=
    Module[{idx, inputs, pars, parTypes, matrixType, basisLengths, basisSize,
            solverName, evaluatorName, result = ""},
           idx = ToString[dataset[[1]]];
           pars = dataset[[2]];
           parTypes = CConversion`GetScalarElementType[Parameters`GetType[#]]& /@ dataset[[2]];
           If[MemberQ[parTypes, CConversion`ScalarType[CConversion`complexScalarCType]],
              matrixType = "Eigen::MatrixXcd",
              matrixType = "Eigen::MatrixXd"
             ];
           basisLengths = Length[GetBasis[#]]& /@ Select[solutions, MemberQ[pars, GetName[#]]&];
           If[Length[DeleteDuplicates[basisLengths]] =!= 1,
              Print["Error: invalid collection of parameters specified:", pars];
              Quit[1];
             ];
           basisSize = ToString[First[basisLengths]];
           solverName = "solver_" <> idx;
           evaluatorName = "basis_" <> idx;
           "auto " <> solverName <> " = create_solver<" <> matrixType
           <> "," <> basisSize <> ">(datasets[" <> idx <> "], " <> evaluatorName <> ");\n"
          ];

CreateLinearSystemSolvers[datasets_List, solutions_List] :=
    Module[{result = ""},
           (result = result <> CreateLinearSystemSolver[#, solutions])& /@ datasets;
           Return[result];
          ];

CalculateCoefficients[datasets_List] :=
    Module[{i, result = ""},
           For[i = 1, i <= Length[datasets], i++,
               (result = result <> "calculate_"
                         <> CConversion`ToValidCSymbolString[#]
                         <> "_coefficients(solver_" <> ToString[datasets[[i,1]]]
                         <> ", datasets[" <> ToString[datasets[[i,1]]]
                         <> "]);\n")& /@ datasets[[i,2]];
              ];
           Return[result];
          ];

ExpandSemiAnalyticSolution[solution_SemiAnalyticSolution] :=
    Module[{par, basis, coeffs},
           par = GetName[solution];
           basis = GetBasis[solution];
           coeffs = CreateCoefficients[solution];
           {par, Dot[coeffs, basis]}
          ];

ExpandSemiAnalyticSolutions[solutions_List] :=
    ExpandSemiAnalyticSolution /@ solutions;

GetSemiAnalyticParameterSubstitutions[solutions_List] :=
    Module[{boundaryValues, expanded},
           boundaryValueRules = Rule[#, GetBoundaryValueParameterSymbol[#]]& /@ (GetBoundaryValueParameters[solutions]);
           expanded = ExpandSemiAnalyticSolutions[solutions];
           {#[[1]], Parameters`ReplaceAllRespectingSARAHHeads[#[[2]], boundaryValueRules]}& /@ expanded
          ];

EvaluateSemiAnalyticSolution[solution_, class_String] :=
    Module[{parameter, basisRules, coeffs},
           parameter = GetName[solution];
           basisRules = Rule[#, GetBoundaryValueParameterSymbol[#]]& /@ (GetBoundaryValueParameters[solution]);
           coeffs = CreateCoefficients[solution];
           Parameters`SetParameter[parameter, Dot[coeffs, GetBasis[solution]] /. basisRules, class]
          ];

EvaluateSemiAnalyticSolutions[solutions_List, class_String:"model."] :=
    Module[{boundaryValues, result = ""},
           boundaryValues = GetBoundaryValueParameterSymbol /@ GetBoundaryValueParameters[solutions];
           (result = result <> EvaluateSemiAnalyticSolution[#, class])& /@ solutions;
           Return[Parameters`CreateLocalConstRefs[boundaryValues] <> "\n" <> result];
          ];

GetExtraBoundaryParametersToSave[solutions_List] :=
    Module[{basisPars},
           basisPars = GetBoundaryValueParameters[solutions];
           basisPars = Select[basisPars, (!(Parameters`IsInputParameter[#]
                                            || Parameters`IsExtraParameter[#]))&];
           {GetBoundaryValueParameterSymbol[#], {}, Parameters`GetType[#]}& /@ basisPars
          ];

SaveBoundaryValueParameter[parameter_, modelPrefix_String:"model->"] :=
    Module[{basisPar = GetBoundaryValueParameterSymbol[parameter],
            parStr = CConversion`ToValidCSymbolString[parameter],
            macro, value},
           Which[Parameters`IsModelParameter[parameter],
                 macro = "MODELPARAMETER",
                 Parameters`IsInputParameter[parameter],
                 macro = "INPUTPARAMETER",
                 Parameters`IsExtraParameter[parameter],
                 macro = "EXTRAPARAMETER",
                 Parameters`IsOutputParameter[parameter],
                 macro = "MODELPARAMETER",
                 Parameters`IsPhase[parameter],
                 macro = "PHASE",
                 True,
                 Print["Error: unrecognized parameter ", parameter];
                 Quit[1];
                ];
           value = macro <> "(" <> parStr <> ")";
           Parameters`SetParameter[basisPar, value, modelPrefix]
          ];

SaveBoundaryValueParameters[solutions_List] :=
    Module[{boundaryValueParameters, result = ""},
           boundaryValueParameters = GetBoundaryValueParameters[solutions];
           (result = result <> SaveBoundaryValueParameter[#])& /@ boundaryValueParameters;
           Return[result];
          ];

SetBoundaryValueParameterFromModel[parameter_] :=
    Module[{parStr = CConversion`ToValidCSymbolString[parameter], value},
           Which[Parameters`IsModelParameter[parameter],
                 macro = "MODELPARAMETER",
                 Parameters`IsInputParameter[parameter],
                 macro = "INPUTPARAMETER",
                 Parameters`IsExtraParameter[parameter],
                 macro = "EXTRAPARAMETER",
                 Parameters`IsOutputParameter[parameter],
                 macro = "MODELPARAMETER",
                 Parameters`IsPhase[parameter],
                 macro = "PHASE",
                 True,
                 Print["Error: unrecognized parameter ", parameter];
                 Quit[1];
                ];
           value = macro <> "(" <> parStr <> ")";
           Parameters`SetParameter[parStr, value]
          ];

SetBoundaryValueParametersFromModel[solutions_List, struct_String:"model."] :=
    Module[{boundaryValues, result = ""},
           boundaryValues = GetBoundaryValueParameters[solutions];
           (result = result <> SetBoundaryValueParameterFromModel[GetBoundaryValueParameterSymbol[#]])& /@ boundaryValues;
           Return[result];
          ];

SetTreeLevelEWSBSolution[ewsbSolution_, solutions_List, substitutions_List, struct_String:"model."] :=
    Module[{parametersFixedByEWSB, i, par, basisPar, parStr, body = "", result = ""},
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
              body = body <> "solutions->evaluate_solutions(model);\n";
              If[substitutions === {},
                 result = result <>
                          "if (is_finite) {\n" <>
                          IndentText[body] <>
                          "} else {\n" <>
                          IndentText["error = 1;\n"] <>
                          "}";,
                 result = result <>
                          "if (is_finite) {\n" <>
                          IndentText[body] <>
                          IndentText[WrapLines[EWSB`SetModelParametersFromEWSB[parametersFixedByEWSB, substitutions, struct]]] <>
                          "} else {\n" <>
                          IndentText["error = 1;\n"] <>
                          "}";
                ];
             ];
           result
          ];

GetAllIndexCombinations[bounds_List] := Tuples[Range[0, #-1]& /@ bounds];

CreateCoefficientsCalculation[solution_SemiAnalyticSolution] :=
    Module[{par = GetName[solution], basis = GetBasis[solution],
            parStr, parType, solverType, dims, i, indices = {}, numCols, coeffs, setCoeffs = "",
            matrixType, name, prototype, function = "", body = ""},
           parStr = CConversion`ToValidCSymbolString[par];
           parType = CConversion`CreateCType[CConversion`GetScalarElementType[Parameters`GetType[par]]];
           If[Parameters`IsRealParameter[par],
              solverType = "Eigen::JacobiSVD<Eigen::MatrixXd>",
              solverType = "Eigen::JacobiSVD<Eigen::MatrixXcd>"
             ];
           dims = Parameters`GetParameterDimensions[par];
           If[dims =!= {1},
              indices = Table["i" <> ToString[i], {i, 1, Length[dims]}];
             ];
           numCols = Times @@ dims;
           If[numCols === 1,
              matrixType = "Eigen::VectorX",
              matrixType = "Eigen::MatrixX"
             ];
           If[Parameters`IsRealParameter[par],
              matrixType = matrixType <> "d",
              matrixType = matrixType <> "cd"
             ];

           If[numCols === 1,
              body = "rhs(j) = data[j]->model.get_" <> parStr <> "();\n";,
              body = MapIndexed[("rhs(j, " <> ToString[First[#2-1]] <> ") = data[j]->model.get_"
                                 <> parStr <> "(" <> Utils`StringJoinWithSeparator[ToString /@ #1, ", "]
                                 <> ");\n")&, GetAllIndexCombinations[dims]];
              body = StringJoin[body];
             ];
           body = "for (std::size_t j = 0; j < n; ++j) {\n"
                  <> IndentText[body] <> "}\n";
           body = "const std::size_t n = data.size();\n"
                  <> matrixType <> " rhs(n" <> If[numCols =!= 1, ", "
                  <> ToString[numCols], ""] <> ");\n" <> body;
           body = body <> matrixType <> " solution = solver.solve(rhs);\n";

           coeffs = CConversion`ToValidCSymbolString[#]& /@ CreateCoefficients[solution];

           For[i = 1, i <= Length[coeffs], i++,
               If[numCols === 1,
                  setCoeffs = setCoeffs <> coeffs[[i]] <> " = solution(" <> ToString[i-1] <> ");\n",
                  setCoeffs = setCoeffs
                              <> StringJoin[MapIndexed[(coeffs[[i]] <> "("
                                                        <> Utils`StringJoinWithSeparator[ToString /@ #1, ", "]
                                                        <> ") = solution(" <> ToString[i-1] <> ", "
                                                        <> ToString[First[#2-1]] <> ");\n")&,
                                                       GetAllIndexCombinations[dims]]];
                 ];
              ];

           body = body <> setCoeffs;

           name = "calculate_" <> parStr <> "_coefficients";
           prototype = "void " <> name <> "(const " <> solverType <> "&, const Data_vector_t&);\n";

           function = "void " <> FlexibleSUSY`FSModelName <> "_semi_analytic_solutions::"
                      <> name <> "(\n"
                      <> IndentText["const " <> solverType <> "& solver, const "
                                    <> FlexibleSUSY`FSModelName
                                    <> "_semi_analytic_solutions::Data_vector_t& data)\n"] <> "{\n";
           function = function <> IndentText[body] <> "}\n";

           {prototype, function}
          ];

CreateCoefficientsCalculations[solutions_List] :=
    Module[{defs, prototypes = "", functions = ""},
           defs = CreateCoefficientsCalculation[#]& /@ solutions;
           prototypes = StringJoin[#[[1]]& /@ defs];
           functions = Utils`StringJoinWithSeparator[#[[2]]& /@ defs, "\n"];
           {prototypes, functions}
          ];

DependsAtMostOn[num_?NumericQ, pars_List] := True;
DependsAtMostOn[True, pars_List] := True;
DependsAtMostOn[False, pars_List] := True;
DependsAtMostOn[KroneckerDelta[i_,j_], pars_List] := True;
DependsAtMostOn[sym_Symbol, pars_List] := MemberQ[pars, sym];
DependsAtMostOn[sym_[indices__] /; And @@ (Parameters`IsIndex /@ {indices}), pars_List] := DependsAtMostOn[sym, pars];
DependsAtMostOn[expr_Times, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Plus, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Power, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[Factorial[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Conjugate[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SARAH`Conj[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Susyno`LieGroups`conj[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ConjugateTranspose[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Transpose[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SARAH`Tp[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SARAH`Adj[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Re[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Im[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Abs[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[AbsSqr[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[AbsSqrt[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[SignedAbsSqrt[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ZeroSqrt[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Arg[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sign[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Round[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Floor[expr_, mult_], pars_List] := DependsAtMostOn[expr, pars] && DependsAtMostOn[mult, pars];
DependsAtMostOn[Floor[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Ceiling[expr_, mult_], pars_List] := DependsAtMostOn[expr, pars] && DependsAtMostOn[mult, pars];
DependsAtMostOn[Ceiling[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sqrt[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cbrt[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[PolyLog[sub_, expr_], pars_List] := DependsAtMostOn[sub, pars] && DependsAtMostOn[expr, pars];
DependsAtMostOn[expr_Min, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Max, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_UnitStep, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_UnitBox, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_UnitTriangle, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Piecewise, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[HoldPattern[(Equal|Unequal|Greater|Less|GreaterEqual|LessEqual)[expr__]], pars_List] := DependsAtMostOn[List[expr], pars];
DependsAtMostOn[HoldPattern[(And|Or|Xor|Nand|Nor)[expr__]], pars_List] := DependsAtMostOn[List[expr], pars];
DependsAtMostOn[Not[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[expr_Which, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_If, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[expr_Switch, pars_List] := DependsAtMostOn[List @@ expr, pars];
DependsAtMostOn[Exp[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Log[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[FiniteLog[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sin[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cos[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Tan[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Csc[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sec[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cot[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSin[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCos[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcTan[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCsc[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSec[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCot[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sinh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Cosh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Tanh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Csch[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Sech[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Coth[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSinh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCosh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcTanh[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCsch[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcSech[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[ArcCoth[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[Total[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[IsClose[first_, second_], pars_List] := DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars];
DependsAtMostOn[IsClose[first_, second_, tol_], pars_List] :=
    DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars] && DependsAtMostOn[tol, pars];
DependsAtMostOn[IsCloseRel[first_, second_], pars_List] := DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars];
DependsAtMostOn[IsCloseRel[first_, second_, tol_], pars_List] :=
    DependsAtMostOn[first, pars] && DependsAtMostOn[second, pars] && DependsAtMostOn[tol, pars];
DependsAtMostOn[IsFinite[expr_], pars_List] := DependsAtMostOn[expr, pars];
DependsAtMostOn[expr_List, pars_List] := And @@ (DependsAtMostOn[#, pars]& /@ expr);
DependsAtMostOn[expr_, pars_List] := False;

End[];

EndPackage[];
