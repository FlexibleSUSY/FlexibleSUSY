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

BeginPackage["ThresholdCorrections`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Constraint`", "Vertices`", "LoopMasses`", "SelfEnergies`", "Utils`"}];

CalculateGaugeCouplings::usage="";
CalculateDeltaAlphaEm::usage="";
CalculateDeltaAlphaS::usage="";
CalculateThetaW::usage="";
GetParameter::usage="";
SetDRbarYukawaCouplingTop::usage="";
SetDRbarYukawaCouplingBottom::usage="";
SetDRbarYukawaCouplingElectron::usage="";
CalculateColorCoupling::usage="";
CalculateElectromagneticCoupling::usage="";
SetDRbarYukawaCouplings::usage="";
GetTwoLoopThresholdHeaders::usage="";
YukawaToMassPrefactor::usage="";

CalculateGaugeCouplings::MissingRelation = "Warning: Coupling `1` is not\
 releated to `2` via DependenceNum: `1` = `3`"

Begin["`Private`"];

DRbarConversion[SARAH`U[1]] := 0;

DRbarConversion[SARAH`SU[n_Integer]] := n/6;

DRbarConversion[group_] := Null;

CalculateColorCoupling[scheme_] :=
    CalculateCoupling[TreeMasses`FindColorGaugeGroup[], scheme];

CalculateElectromagneticCoupling[scheme_] :=
  Module[{conversion},
          conversion = Switch[scheme,
                              FlexibleSUSY`DRbar, 1/3,
                              FlexibleSUSY`MSbar, 0,
                              _, Message[CalculateCoupling::UnknownRenormalizationScheme, scheme]; 0
                             ];
         CalculateCoupling[{SARAH`electricCharge, SARAH`electricCharge, SARAH`U[1]}, scheme] + conversion
        ];

CalculateCoupling::UnknownRenormalizationScheme = "Unknown\
 renormalization scheme `1`.";

GetLorentzRepresentationFactor[particle_] :=
    Which[IsMajoranaFermion[particle], 2/3,
          IsDiracFermion[   particle], 4/3,
          IsRealScalar[     particle], 1/6,
          IsComplexScalar[  particle], 1/3,
          IsVector[         particle], 1/3,
          True                       , 1];

GetMultiplicityForUnbrokenGroups[particle_, group_] :=
    Times @@ Cases[SARAH`getIndizesWI[particle], {Except[SARAH`generation | group], i_} :> i];

(* Calculate threshold correction for a gauge coupling from SM
   (MS-bar) to a given renormalization scheme in the given model. *)
CalculateCoupling[{coupling_, name_, group_}, scheme_] :=
    Module[{susyParticles, prefactor, i, result = 0, particle,
            conversion},
           susyParticles = TreeMasses`GetSusyParticles[];
           For[i = 1, i <= Length[susyParticles], i++,
               particle = susyParticles[[i]];
               (* Eq.(11) and (13) of arXiv:1406.2319 *)
               prefactor =
                  GetLorentzRepresentationFactor[particle] \
                  GetMultiplicityForUnbrokenGroups[particle,name] \
                  If[coupling === SARAH`electricCharge,
                     TreeMasses`GetElectricCharge[particle]^2,
                     SA`Dynkin[particle, Position[SARAH`Gauge, name][[1, 1]]]
                    ];
               If[!NumericQ[prefactor], prefactor = 0];
               (* sum over generations *)
               If[GetDimension[particle] == 1,
                  result -= prefactor Global`FiniteLog[Abs[FlexibleSUSY`M[particle]/Global`currentScale]];,
                  result -= Sum[prefactor Global`FiniteLog[Abs[FlexibleSUSY`M[particle][i-1]/Global`currentScale]],
                                {i,TreeMasses`GetDimensionStartSkippingSMGoldstones[particle],GetDimension[particle]}];
                 ];
              ];
           conversion = Switch[scheme,
                               FlexibleSUSY`DRbar, DRbarConversion[group],
                               FlexibleSUSY`MSbar, 0,
                               _, Message[CalculateCoupling::UnknownRenormalizationScheme, scheme]; 0
                              ];
           result + conversion
          ];

CalculateDeltaAlphaEm[renormalizationScheme_] :=
    Module[{result, deltaSusy, deltaSM, prefactor, topQuark},
           topQuark = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           prefactor = Global`alphaEm / (2 Pi);
           deltaSM = -16/9 Global`FiniteLog[Abs[topQuark/Global`currentScale]];
           deltaSusy = CalculateElectromagneticCoupling[renormalizationScheme];
           result = Parameters`CreateLocalConstRefs[{deltaSusy, deltaSM}] <> "\n" <>
                    "const double delta_alpha_em_SM = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSM] <> ";\n\n" <>
                    "const double delta_alpha_em = " <>
                    CConversion`RValueToCFormString[prefactor * deltaSusy] <> ";\n\n" <>
                    "return delta_alpha_em + delta_alpha_em_SM;\n";
           Return[result];
          ];

CalculateDeltaAlpha2LSM[] :=
"if (model->get_thresholds() > 1 && model->get_threshold_corrections().alpha_s > 1) {\n" <>
IndentText["\
sm_fourloop_as::Parameters pars;
pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
pars.mt   = model->get_" <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[],True,True]] <> ";
pars.Q    = model->get_scale();

const auto das_1L = sm_fourloop_as::delta_alpha_s_1loop_as(pars);
const auto das_2L = sm_fourloop_as::delta_alpha_s_2loop_as_as(pars);

delta_alpha_s_2loop = das_2L - Sqr(das_1L);"
] <> "
}

";

CalculateDeltaAlpha3LSM[] :=
"if (model->get_thresholds() > 2 && model->get_threshold_corrections().alpha_s > 2) {\n" <>
IndentText["\
sm_fourloop_as::Parameters pars;
pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
pars.mt   = model->get_" <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[],True,True]] <> ";
pars.Q    = model->get_scale();

const auto das_1L = sm_fourloop_as::delta_alpha_s_1loop_as(pars);
const auto das_2L = sm_fourloop_as::delta_alpha_s_2loop_as_as(pars);
const auto das_3L = sm_fourloop_as::delta_alpha_s_3loop_as_as_as(pars);

delta_alpha_s_3loop = das_3L + Power3(das_1L) - 2. * das_1L * das_2L;"
] <> "
}

";

CalculateDeltaAlpha4LSM[] :=
"if (model->get_thresholds() > 3 && model->get_threshold_corrections().alpha_s > 3) {\n" <>
IndentText["\
sm_fourloop_as::Parameters pars;
pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
pars.mt   = model->get_" <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[],True,True]] <> ";
pars.Q    = model->get_scale();

const auto das_1L = sm_fourloop_as::delta_alpha_s_1loop_as(pars);
const auto das_2L = sm_fourloop_as::delta_alpha_s_2loop_as_as(pars);
const auto das_3L = sm_fourloop_as::delta_alpha_s_3loop_as_as_as(pars);
const auto das_4L = sm_fourloop_as::delta_alpha_s_4loop_as_as_as_as(pars);

delta_alpha_s_4loop = das_4L - 2. * das_1L * das_3L - Power2(das_2L) + 3. * Power2(das_1L) * das_2L - Power4(das_1L);"
] <> "
}

";

CalculateDeltaAlpha2LMSSM[] :=
"if (model->get_thresholds() > 1 && model->get_threshold_corrections().alpha_s > 1) {\n" <>
IndentText["\
" <> Parameters`CreateLocalConstRefs[Parameters`GetEffectiveMu[]] <> "

double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double msd_1, msd_2, theta_d;

model->" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <> ";
model->" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <> ";
model->" <> TreeMasses`CallGenerationHelperFunctionName[2, SARAH`BottomSquark, "msd_1", "msd_2", "theta_d"] <> ";

mssm_twoloop_as::Parameters pars;
pars.g3   = model->get_" <> CConversion`RValueToCFormString[SARAH`strongCoupling /. Parameters`ApplyGUTNormalization[]] <> "();
pars.yt   = model->get_" <> CConversion`RValueToCFormString[Parameters`GetThirdGeneration[SARAH`UpYukawa]] <> ";
pars.yb   = model->get_" <> CConversion`RValueToCFormString[Parameters`GetThirdGeneration[SARAH`DownYukawa]] <> ";
pars.mt   = model->get_" <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[],True,True]] <> ";
pars.mb   = model->get_" <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMBottomQuarkMultiplet[],True,True]] <> ";
pars.mg   = model->get_" <> CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]] <> "();
pars.mst1 = mst_1;
pars.mst2 = mst_2;
pars.msb1 = msb_1;
pars.msb2 = msb_2;
pars.msd1 = msd_1;
pars.msd2 = msd_2;
pars.xt   = Sin(2*theta_t) * (Sqr(mst_1) - Sqr(mst_2)) / (2. * pars.mt);
pars.xb   = Sin(2*theta_b) * (Sqr(msb_1) - Sqr(msb_2)) / (2. * pars.mb);
pars.mw   = model->get_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[SARAH`VectorW]] <> "();
pars.mz   = model->get_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[SARAH`VectorZ]] <> "();
pars.mh   = model->get_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[SARAH`HiggsBoson]] <>"(0);
pars.mH   = model->get_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[SARAH`HiggsBoson]] <> "(1);
pars.mC   = model->get_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[SARAH`ChargedHiggs]] <>
   "(" <> ToString[TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`ChargedHiggs]-1] <> ");
pars.mA   = model->get_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[SARAH`PseudoScalar]] <>
   "(" <> ToString[TreeMasses`GetDimensionStartSkippingGoldstones[SARAH`PseudoScalar]-1] <> ");
pars.mu   = " <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ";
pars.tb   = model->get_" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <>
   "() / model->get_" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> "();
pars.Q    = model->get_scale();

delta_alpha_s_2loop =
   - Sqr(delta_alpha_s_1loop)/4.
   - 2.*(
      + mssm_twoloop_as::delta_alpha_s_2loop_as_as(pars)
      + mssm_twoloop_as::delta_alpha_s_2loop_at_as(pars)
      + mssm_twoloop_as::delta_alpha_s_2loop_ab_as(pars)
     );"] <> "
}

";

CalculateDeltaAlphaS[renormalizationScheme_] :=
    Module[{result, deltaSusy, deltaSM, prefactor, topQuark},
           topQuark = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           prefactor = Global`alphaS / (2 Pi);
           deltaSM = - 2/3 Global`FiniteLog[Abs[topQuark/Global`currentScale]];
           deltaSusy = CalculateColorCoupling[renormalizationScheme];
           Parameters`CreateLocalConstRefs[{deltaSusy, deltaSM}] <> "\n" <>
           "const double delta_alpha_s_SM = " <>
           CConversion`RValueToCFormString[prefactor * deltaSM] <> ";\n\n" <>
           "const double delta_alpha_s = " <>
           CConversion`RValueToCFormString[prefactor * deltaSusy] <> ";\n\n" <>
           "const double delta_alpha_s_1loop = delta_alpha_s + delta_alpha_s_SM;\n" <>
           "double delta_alpha_s_2loop = 0.;\n" <>
           "double delta_alpha_s_3loop = 0.;\n" <>
           "double delta_alpha_s_4loop = 0.;\n\n" <>
           If[FlexibleSUSY`UseMSSMAlphaS2Loop === True, CalculateDeltaAlpha2LMSSM[], ""] <>
           If[FlexibleSUSY`UseSMAlphaS3Loop === True, CalculateDeltaAlpha2LSM[], ""] <>
           If[FlexibleSUSY`UseSMAlphaS3Loop === True, CalculateDeltaAlpha3LSM[], ""] <>
           If[FlexibleSUSY`UseSMAlphaS4Loop === True, CalculateDeltaAlpha4LSM[], ""] <>
           "return delta_alpha_s_1loop + delta_alpha_s_2loop + delta_alpha_s_3loop + delta_alpha_s_4loop;\n"
          ];

GetPrefactor[expr_Plus, _] := 1;

GetPrefactor[expr_Integer, _] := 1;

GetPrefactor[expr_Symbol, _] := 1;

GetPrefactor[expr_Times, yukawa_] :=
    Module[{factors, prefactors},
           factors = List @@ expr;
           prefactors = Select[factors, FreeQ[#, yukawa]&];
           Times @@ prefactors
          ];

ExtractSymbols[sym_?NumberQ] := {};

ExtractSymbols[sym_Symbol] := {sym};

ExtractSymbols[HoldPattern[SARAH`sum[_,_,_,expr_]]] := ExtractSymbols[expr];

ExtractSymbols[expr_Plus] := Flatten[ExtractSymbols /@ (List @@ expr)];

ExtractSymbols[expr_Times] := Flatten[ExtractSymbols /@ (List @@ expr)];

ExtractSymbols[sym_[_,_]] := {sym};

ToMatrixExpression[{}] := Null;

ToMatrixExpression[expr_ /; Head[expr] =!= List] := expr;

ToMatrixExpression[list_List] :=
    Module[{dim, symbol, matrix, i, k, diag, expression = Null,
            expandedList, permutations},
           dim = Length[list];
           symbol = ExtractSymbols[list[[1,1]]];
           (* expand sum[] *)
           expandedList = Expand[list //. SARAH`sum[idx_,start_,stop_,expr_] :> Sum[expr,{idx,start,stop}]];
           If[Length[symbol] == 1,
              symbol = symbol[[1]];
              matrix = Table[symbol[i,k], {i,1,dim}, {k,1,dim}];
              diag = DiagonalMatrix[Table[symbol[i,i], {i,1,dim}]];
              Which[matrix === expandedList, expression = symbol;,
                    Transpose[matrix] === expandedList, expression = SARAH`Tp[symbol];,
                    diag === expandedList, expression = FlexibleSUSY`Diag[symbol];
                   ];
              ,
              (* create all permutations of matrix products *)
              permutations = Permutations[symbol];
              For[i = 1, i <= Length[permutations], i++,
                  matrix = Table[#[i,k],{i,1,dim},{k,1,dim}]& /@ permutations[[i]];
                  matrix = Expand[Dot @@ matrix];
                  If[matrix === expandedList,
                     (* found a combination which yields our list *)
                     expression = SARAH`MatMul @@ permutations[[i]];
                     Break[];
                    ];
                 ];
             ];
           Return[expression];
          ];

(* Solve the equation #1 == #2 for #3 *)
InvertRelation[Transpose[sym_], expr_, sym_] :=
    {sym, SARAH`Tp[expr]};

InvertRelation[SARAH`Tp[sym_], expr_, sym_] :=
    {sym, SARAH`Tp[expr]};

InvertRelation[ConjugateTranspose[sym_], expr_, sym_] :=
    {sym, SARAH`Adj[expr]};

InvertRelation[SARAH`Adj[sym_], expr_, sym_] :=
    {sym, SARAH`Adj[expr]};

InvertRelation[FlexibleSUSY`Diag[sym_], expr_, sym_] :=
    {sym, FlexibleSUSY`Diag[expr]};

InvertRelation[sym_, expr_, sym_] :=
    {sym, expr};

InvertRelation[sym_[i1_,i2_], expr_, sym_] :=
    {sym[i1,i2], expr};

(* remove matrices from the left *)
InvertRelation[SARAH`MatMul[SARAH`Adj[U_],X___,sym_,V___], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[X,sym,V], SARAH`MatMul[U,expr], sym];

InvertRelation[SARAH`MatMul[U_,X___,sym_,V___], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[X,sym,V], SARAH`MatMul[SARAH`Adj[U],expr], sym];

(* remove matrices from the right *)
InvertRelation[SARAH`MatMul[sym_,V___,SARAH`Adj[U_]], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[sym,V], SARAH`MatMul[expr,U], sym];

InvertRelation[SARAH`MatMul[sym_,U__], expr_, sym_] :=
    InvertRelation[SARAH`MatMul[sym], SARAH`MatMul[expr,SARAH`Adj[U]], sym];

InvertRelation[SARAH`MatMul[sym_], expr_, sym_] :=
    InvertRelation[sym, expr, sym];

InvertRelation[sym_, expr_, other_] :=
    Block[{},
          Print["Error: InvertRelation: don't know how to solve equation: ",
                sym, " == ", expr, " for ", other];
          Quit[1];
         ];

InvertMassRelation[fermion_, yukawa_] :=
    Module[{massMatrix, polynom, prefactor, matrixExpression, dim},
           If[TreeMasses`IsUnmixed[fermion],
              massMatrix = TreeMasses`GetMassOfUnmixedParticle[fermion];
              massMatrix = TreeMasses`ReplaceDependencies[massMatrix];
              massMatrix = Vertices`StripGroupStructure[massMatrix, {SARAH`ct1, SARAH`ct2}];
              ,
              massMatrix = SARAH`MassMatrix[fermion];
             ];
           dim = Length[massMatrix];
           If[massMatrix === Table[0, {dim}, {dim}],
              Return[{yukawa,CConversion`ZEROMATRIX[dim,dim]}];
             ];
           polynom = Factor[massMatrix /. List -> Plus];
           prefactor = GetPrefactor[polynom, yukawa];
           matrixExpression = ToMatrixExpression[massMatrix / prefactor];
           If[matrixExpression === Null,
              Print["Error: could not convert list to matrix expression: ",
                    massMatrix / prefactor];
              Quit[1];
             ];
           InvertRelation[matrixExpression, fermion / prefactor, yukawa]
          ];

YukawaToMassPrefactor[fermion_, yukawa_] :=
    Module[{massMatrix, polynom},
           If[TreeMasses`IsUnmixed[fermion],
              massMatrix = TreeMasses`GetMassOfUnmixedParticle[fermion];
              massMatrix = TreeMasses`ReplaceDependencies[massMatrix];
              massMatrix = Vertices`StripGroupStructure[massMatrix, {SARAH`ct1, SARAH`ct2}];
              ,
              massMatrix = SARAH`MassMatrix[fermion];
             ];
           polynom = Factor[massMatrix /. List -> Plus];
           GetPrefactor[polynom, yukawa]
          ];

SetDRbarYukawaCouplingTop[settings_] :=
    SetDRbarYukawaCouplingFermion[TreeMasses`GetSMTopQuarkMultiplet[], SARAH`UpYukawa, Global`upQuarksDRbar, settings];

SetDRbarYukawaCouplingBottom[settings_] :=
    SetDRbarYukawaCouplingFermion[TreeMasses`GetSMBottomQuarkMultiplet[], SARAH`DownYukawa, Global`downQuarksDRbar, settings];

SetDRbarYukawaCouplingElectron[settings_] :=
    SetDRbarYukawaCouplingFermion[TreeMasses`GetSMTauLeptonMultiplet[], SARAH`ElectronYukawa, Global`downLeptonsDRbar, settings];

SetDRbarYukawaCouplingFermionMatrix[fermion_, yukawa_, mass_, setting_] :=
    Module[{y, f = setting},
           If[setting === Automatic,
              {y, f} = InvertMassRelation[fermion, yukawa];
              f = f /. fermion -> mass;
              ,
              y = yukawa;
             ];
           Parameters`CreateLocalConstRefs[f] <>
           Parameters`SetParameter[yukawa, f, "MODEL->"]
          ];

SetDRbarYukawaCouplings[] :=
    Module[{y, f, fermion, yukawa, mass, term = {0,0,0}, i},
           fermion = {TreeMasses`GetSMTopQuarkMultiplet[], TreeMasses`GetSMBottomQuarkMultiplet[], TreeMasses`GetSMTauLeptonMultiplet[]};
           yukawa  = {SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa};
           mass    = {Global`upQuarksDRbar, Global`downQuarksDRbar, Global`downLeptonsDRbar};
           For[i = 1, i <= 3, i++,
               {y, f} = InvertMassRelation[fermion[[i]], yukawa[[i]]];
               term[[i]] = f /. fermion[[i]] -> mass[[i]];
              ];
           Parameters`CreateLocalConstRefs[term] <>
           "model.set_" <> CConversion`ToValidCSymbolString[yukawa[[1]]] <> "(" <> CConversion`RValueToCFormString[term[[1]]] <> ");\n" <>
           "model.set_" <> CConversion`ToValidCSymbolString[yukawa[[2]]] <> "(" <> CConversion`RValueToCFormString[term[[2]]] <> ");\n" <>
           "model.set_" <> CConversion`ToValidCSymbolString[yukawa[[3]]] <> "(" <> CConversion`RValueToCFormString[term[[3]]] <> ");\n"
          ];

SetDRbarYukawaCouplingFermionElement[{y_, expr_}] :=
    Parameters`SetParameter[y, expr, "MODEL->"];

SetDRbarYukawaCouplingFermion[fermion_, yukawa_, mass_, settings_] :=
    Module[{f, p, i, result = ""},
           (* check for matrix setting *)
           f = Cases[settings, {yukawa, s_} :> s];
           If[Flatten[f] =!= {},
              For[i = 1, i <= Length[f], i++,
                  result = result <> SetDRbarYukawaCouplingFermionMatrix[fermion, yukawa, mass, f[[i]]];
                 ];
              Return[result];
             ];
           (* check for matrix element setting *)
           f = Cases[settings, p:{yukawa[__], s_} :> p];
           If[Flatten[f] =!= {},
              For[i = 1, i <= Length[f], i++,
                  result = result <> SetDRbarYukawaCouplingFermionElement[f[[i]]];
                 ];
              Return[Parameters`CreateLocalConstRefs[(#[[2]])& /@ f] <> result];
             ];
           ""
          ];

MultiplyBy[factor_ /; factor == 1] := "";

MultiplyBy[factor_] :=
    " * " <> CConversion`RValueToCFormString[factor];

GetParameter[par_[idx1_,idx2_], factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[par] <>
    "(" <> CConversion`RValueToCFormString[idx1] <> "," <>
    CConversion`RValueToCFormString[idx2] <> ")" <>
    MultiplyBy[factor];

GetParameter[FlexibleSUSY`M[par_], factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[FlexibleSUSY`M[par]] <> "()" <>
    MultiplyBy[factor];

GetParameter[par_[idx_], factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[par] <>
    "(" <> CConversion`RValueToCFormString[idx] <> ")" <>
    MultiplyBy[factor];

GetParameter[par_, factor_:1] :=
    "MODEL->get_" <> CConversion`RValueToCFormString[par] <> "()" <>
    MultiplyBy[factor];

CalculateThetaWFromFermiConstant[] :=
    Module[{},
    "\
" <> FlexibleSUSY`FSModelName <> "_weinberg_angle::Sm_parameters sm_pars;
sm_pars.fermi_constant = qedqcd.displayFermiConstant();
sm_pars.mw_pole = qedqcd.displayPoleMW();
sm_pars.mz_pole = qedqcd.displayPoleMZ();
sm_pars.mt_pole = qedqcd.displayPoleMt();
sm_pars.alpha_s = calculate_alpha_s_SM5_at(qedqcd, qedqcd.displayPoleMt());
sm_pars.higgs_index = higgs_idx;

const int number_of_iterations =
    std::max(20, static_cast<int>(std::abs(-log10(MODEL->get_precision()) * 10)));

" <> FlexibleSUSY`FSModelName <> "_weinberg_angle weinberg(MODEL, sm_pars);
weinberg.set_number_of_loops(MODEL->get_threshold_corrections().sin_theta_w);
weinberg.set_number_of_iterations(number_of_iterations);

try {
   const auto result = weinberg.calculate();
   THETAW = ArcSin(result.first);

   if (MODEL->get_thresholds() && MODEL->get_threshold_corrections().sin_theta_w > 0)
      qedqcd.setPoleMW(result.second);

   MODEL->get_problems().unflag_no_sinThetaW_convergence();
} catch (const Error& e) {
   VERBOSE_MSG(e.what_detailed());
   MODEL->get_problems().flag_no_sinThetaW_convergence();
}"
          ];

CalculateThetaWFromMW[expr_] :=
    Module[{subst, weinbergAngle, result},
           subst = { SARAH`Mass[SARAH`VectorW] -> FlexibleSUSY`MWDRbar,
                     SARAH`Mass[SARAH`VectorZ] -> FlexibleSUSY`MZDRbar,
                     SARAH`electricCharge      -> FlexibleSUSY`EDRbar };
           (* read weinberg angle from DependenceNum *)
           weinbergAngle = Parameters`FindSymbolDef[SARAH`Weinberg] /. subst;
           If[weinbergAngle === None || weinbergAngle === Null,
              weinbergAngle = expr /. subst;
              If[weinbergAngle === None || weinbergAngle === Null,
                 Print["Warning: No expression for the Weinberg angle defined, setting it to 0."];
                 weinbergAngle = 0;
                ];
             ];
           result = Parameters`CreateLocalConstRefs[{weinbergAngle}] <>
                    "THETAW = " <>
                    CConversion`RValueToCFormString[weinbergAngle] <> ";\n";
           Return[result];
          ];

CalculateThetaW[input_] :=
    Switch[input,
           FlexibleSUSY`FSFermiConstant, CalculateThetaWFromFermiConstant[],
           FlexibleSUSY`FSMassW, CalculateThetaWFromMW[FlexibleSUSY`FSWeakMixingAngleExpr],
           _,
           Print["Error: CalculateThetaW[", input, "]: unknown input ", input];
           Print["   Using default input: ", FlexibleSUSY`FSMassW];
           CalculateThetaWFromMW[FlexibleSUSY`FSWeakMixingAngleExpr]
          ];

WarnIfFreeQ[coupling_, expr_, sym_] :=
    If[FreeQ[expr, sym],
       Message[CalculateGaugeCouplings::MissingRelation, coupling, sym, expr];
      ];

CheckPerturbativityOf[coup_, name_String] :=
"
if (IsFinite(new_" <> name <> ")) {
   model->get_problems().unflag_non_perturbative_parameter(" <> FlexibleSUSY`FSModelName <> "_info::" <> CConversion`ToValidCSymbolString[coup] <> ");
} else {
   model->get_problems().flag_non_perturbative_parameter(
      " <> FlexibleSUSY`FSModelName <> "_info::" <> CConversion`ToValidCSymbolString[coup] <> ", new_" <> name <> ", get_scale());
   new_" <> name <> " = Electroweak_constants::" <> name <> ";
}";

CalculateGaugeCouplings[] :=
    Module[{subst, g1Def, g2Def, g3Def},
           subst = { SARAH`Mass[SARAH`VectorW] -> FlexibleSUSY`MWDRbar,
                     SARAH`Mass[SARAH`VectorZ] -> FlexibleSUSY`MZDRbar,
                     SARAH`electricCharge      -> FlexibleSUSY`EDRbar,
                     SARAH`Weinberg            -> FlexibleSUSY`ThetaWDRbar };
           g1Def = (Parameters`FindSymbolDef[SARAH`hyperchargeCoupling]
                    / Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]);
           g2Def = (Parameters`FindSymbolDef[SARAH`leftCoupling]
                    / Parameters`GetGUTNormalization[SARAH`leftCoupling]);
           g3Def = (Parameters`FindSymbolDef[SARAH`strongCoupling]
                    / Parameters`GetGUTNormalization[SARAH`strongCoupling]);
           WarnIfFreeQ[SARAH`hyperchargeCoupling, g1Def, SARAH`electricCharge];
           WarnIfFreeQ[SARAH`hyperchargeCoupling, g1Def, SARAH`Weinberg];
           WarnIfFreeQ[SARAH`leftCoupling       , g2Def, SARAH`electricCharge];
           WarnIfFreeQ[SARAH`leftCoupling       , g2Def, SARAH`Weinberg];
           WarnIfFreeQ[SARAH`strongCoupling     , g3Def, Parameters`GetParameterFromDescription["Alpha Strong"]];
           g1Def = g1Def /. subst;
           g2Def = g2Def /. subst;
           g3Def = g3Def /. subst;
           Parameters`CreateLocalConstRefs[{g1Def, g2Def, g3Def}] <>
           "new_g1 = " <> CConversion`RValueToCFormString[g1Def] <> ";\n" <>
           "new_g2 = " <> CConversion`RValueToCFormString[g2Def] <> ";\n" <>
           "new_g3 = " <> CConversion`RValueToCFormString[g3Def] <> ";\n" <>
           If[ValueQ[SARAH`hyperchargeCoupling], CheckPerturbativityOf[SARAH`hyperchargeCoupling, "g1"], ""] <> "\n" <>
           If[ValueQ[SARAH`leftCoupling]       , CheckPerturbativityOf[SARAH`leftCoupling, "g2"], ""]
          ];

GetTwoLoopThresholdHeaders[] :=
    Module[{result = ""},
           If[FlexibleSUSY`UseSMYukawa2Loop === True,
              result = result <> "#include \"sm_twoloop_mt.hpp\"\n";
             ];
           If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
              result = "#include \"mssm_twoloop_mb.hpp\"\n" <>
                       "#include \"mssm_twoloop_mt.hpp\"\n" <>
                       "#include \"mssm_twoloop_mtau.hpp\"\n";
             ];
           If[FlexibleSUSY`UseSMAlphaS3Loop === True ||
              FlexibleSUSY`UseSMAlphaS4Loop === True,
              result = result <> "#include \"sm_fourloop_as.hpp\"\n";
             ];
           If[FlexibleSUSY`UseMSSMAlphaS2Loop === True,
              result = result <> "#include \"mssm_twoloop_as.hpp\"\n";
             ];
           result
          ];

End[];

EndPackage[];
