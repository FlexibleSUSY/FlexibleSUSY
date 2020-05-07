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

BeginPackage["EffectiveCouplings`", {"SARAH`", "CConversion`", "Parameters`", "SelfEnergies`", "TreeMasses`", "TextFormatting`", "Utils`", "Vertices`", "Constraint`"}];

InitializeEffectiveCouplings::usage="";
InitializeMixingFromModelInput::usage="";
GetMixingMatrixFromModel::usage="";
GetNeededVerticesList::usage="";
CalculatePartialWidths::usage="";
CalculateQCDAmplitudeScalingFactors::usage="";
CalculateQCDScalingFactor::usage="";
CreateEffectiveCouplingsGetters::usage="";
CreateEffectiveCouplingsDefinitions::usage="";
CreateEffectiveCouplingsInit::usage="";
CreateEffectiveCouplingsCalculation::usage="";
CreateEffectiveCouplings::usage="";

(* @todo error handle when e.g. the strong coupling is not defined *)

Begin["`Private`"];

InitializeMixingFromModelInput[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{i, macro, symbol, result = ""},
           symbol = TreeMasses`GetMixingMatrixSymbol[massMatrix];
           If[SARAH`SupersymmetricModel,
              macro = "MODELPARAMETER";,
              macro = "PHYSICAL";
             ];
           If[symbol === Null,
              Return[""];
             ];
           If[Length[symbol] > 1,
              (result = result <> ", " <> CConversion`ToValidCSymbolString[#] <> "(" <> macro <> "("
                        <> CConversion`ToValidCSymbolString[#] <> "))")& /@ symbol;,
              result = ", " <> CConversion`ToValidCSymbolString[symbol]
                       <> "(" <> macro <> "(" <> CConversion`ToValidCSymbolString[symbol] <> "))";
             ];
           result
          ];

GetMixingMatrixFromModel[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{i, macro, symbol, result = ""},
           symbol = TreeMasses`GetMixingMatrixSymbol[massMatrix];
           If[symbol === Null,
              Return[""];
             ];
           If[SARAH`SupersymmetricModel,
              macro = "MODELPARAMETER";,
              macro = "PHYSICAL";
             ];
           If[Length[symbol] > 1,
              (result = result <> CConversion`ToValidCSymbolString[#] <> " = " <> macro <> "("
                        <> CConversion`ToValidCSymbolString[#] <> ");\n")& /@ symbol;,
              result = CConversion`ToValidCSymbolString[symbol] <> " = " <> macro <> "("
                       <> CConversion`ToValidCSymbolString[symbol] <> ");\n";
             ];
           result
          ];

CalculateQCDAmplitudeScalingFactors[] :=
    Module[{coeff, scalarQCD, fermionQCD, pseudoscalarQCD, body = "",
            scalarScalarLoopFactor = "", scalarFermionLoopFactor = "",
            pseudoscalarFermionLoopFactor = ""},
           scalarQCD = 1 + 2 SARAH`strongCoupling^2 / (3 Pi^2);
           body = "result = " <> CConversion`RValueToCFormString[scalarQCD] <> ";";
           scalarScalarLoopFactor = Parameters`CreateLocalConstRefs[{SARAH`strongCoupling}]
                                    <> "if (m_loop > m_decay) {\n"
                                    <> TextFormatting`IndentText[body] <> "\n}";
           scalarFermionLoopFactor = Parameters`CreateLocalConstRefs[{SARAH`strongCoupling}]
                                     <> "result = 1.0 + "
                                     <> CConversion`RValueToCFormString[SARAH`strongCoupling^2 / (4 Pi^2)]
                                     <> " * scalar_diphoton_fermion_loop(m_decay, m_loop);\n";
           pseudoscalarFermionLoopFactor = Parameters`CreateLocalConstRefs[{SARAH`strongCoupling}]
                                           <> "result = 1.0 + "
                                           <> CConversion`RValueToCFormString[SARAH`strongCoupling^2 / (4 Pi^2)]
                                           <> " * pseudoscalar_diphoton_fermion_loop(m_decay, m_loop);\n";
           scalarScalarLoopFactor = TextFormatting`IndentText[scalarScalarLoopFactor];
           scalarFermionLoopFactor = TextFormatting`IndentText[scalarFermionLoopFactor];
           pseudoscalarFermionLoopFactor = TextFormatting`IndentText[pseudoscalarFermionLoopFactor];
           {scalarScalarLoopFactor, scalarFermionLoopFactor, pseudoscalarFermionLoopFactor}
          ];

CalculateQCDScalingFactor[] :=
    Module[{nloQCD, nnloQCD, nnnloQCD, scalarFactor = "", pseudoscalarFactor = ""},
           (* NLO, NNLO and NNNLO contributions to scalar coupling *)
           nloQCD = (95 / 4 - 7 / 6 Symbol["Nf"]) SARAH`strongCoupling^2 / (4 Pi^2);
           nnloQCD = 149533 / 288 - 363 Zeta[2] / 8 - 495 Zeta[3] / 8 + 19 Symbol["l"] / 8;
           nnloQCD = nnloQCD + Symbol["Nf"] (-4157 / 72 + 11 Zeta[2] / 2 + 5 Zeta[3] / 4 + 2 Symbol["l"] / 3);
           nnloQCD = nnloQCD + Symbol["Nf"]^2 (127 / 108 - Zeta[2] / 6);
           nnloQCD = nnloQCD SARAH`strongCoupling^4 / (16 Pi^4);
           nnnloQCD = 467.683620788 + 122.440972222 Symbol["l"] + 10.9409722222 Symbol["l"]^2;
           nnnloQCD = nnnloQCD SARAH`strongCoupling^6 / (64 Pi^6);
           scalarFactor = scalarFactor <> "const double nlo_qcd = " <> CConversion`RValueToCFormString[nloQCD] <> ";\n";
           scalarFactor = scalarFactor <> "const double nnlo_qcd = " <> CConversion`RValueToCFormString[nnloQCD] <> ";\n";
           scalarFactor = scalarFactor <> "const double nnnlo_qcd = " <> CConversion`RValueToCFormString[nnnloQCD] <> ";\n";
           scalarFactor = Parameters`CreateLocalConstRefs[{nloQCD, nnloQCD, nnnloQCD}] <> "\n" <> scalarFactor;
           (* NLO, NNLO and NNNLO contributions to pseudoscalar coupling *)
           nloQCD = (97 / 4 - 7 / 6 Symbol["Nf"]) SARAH`strongCoupling^2 / (4 Pi^2);
           nnloQCD = (237311 / 864 - 529 Zeta[2] / 24 - 445 Zeta[3] / 8 + 5 Symbol["l"]);
           nnloQCD = nnloQCD SARAH`strongCoupling^4 / (16 Pi^4);
           nnnloQCD = 0;
           pseudoscalarFactor = pseudoscalarFactor <> "const double nlo_qcd = "
                                <> CConversion`RValueToCFormString[nloQCD] <> ";\n";
           pseudoscalarFactor = pseudoscalarFactor <> "const double nnlo_qcd = "
                                <> CConversion`RValueToCFormString[nnloQCD] <> ";\n";
           pseudoscalarFactor = pseudoscalarFactor <> "const double nnnlo_qcd = "
                                <> CConversion`RValueToCFormString[nnnloQCD] <> ";\n";
           pseudoscalarFactor = Parameters`CreateLocalConstRefs[{nloQCD, nnloQCD, nnnloQCD}] <> "\n" <> pseudoscalarFactor;
           {scalarFactor, pseudoscalarFactor}
          ];

GetAllowedCouplingsForModel[] :=
    Module[{dim,
            valid = {FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                     FlexibleSUSYObservable`CpHiggsGluonGluon,
                     FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                     FlexibleSUSYObservable`CpPseudoScalarGluonGluon}
           },
           If[FreeQ[TreeMasses`GetParticles[], SARAH`HiggsBoson],
              valid = DeleteCases[valid, a_ /; (a === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                                                a === FlexibleSUSYObservable`CpHiggsGluonGluon)];
             ];
           If[FreeQ[TreeMasses`GetParticles[], SARAH`PseudoScalar],
              valid = DeleteCases[valid, a_ /; (a === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                                                a === FlexibleSUSYObservable`CpPseudoScalarGluonGluon)];
             ];
           valid
          ];

GetExternalStates[couplingSymbol_] :=
    Module[{particle, vectorBoson},
           Which[couplingSymbol === FlexibleSUSYObservable`CpHiggsPhotonPhoton,
                 particle = SARAH`HiggsBoson;
                 vectorBoson = SARAH`VectorP;,
                 couplingSymbol === FlexibleSUSYObservable`CpHiggsGluonGluon,
                 particle = SARAH`HiggsBoson;
                 vectorBoson = SARAH`VectorG;,
                 couplingSymbol === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton,
                 particle = SARAH`PseudoScalar;
                 vectorBoson = SARAH`VectorP;,
                 couplingSymbol === FlexibleSUSYObservable`CpPseudoScalarGluonGluon,
                 particle = SARAH`PseudoScalar;
                 vectorBoson = SARAH`VectorG;,
                 True, particle = Null; vectorBoson = Null
                ];
           {particle, vectorBoson}
          ];

CalculatePartialWidths[couplings_List] :=
    Module[{couplingSymbol, dim, particle, vectorBoson, particlesStr,
            massStr, functionName, couplingName,
            body, prototypes = {}, functions = {}},
           For[i = 1, i <= Length[couplings], i++,
               couplingSymbol = couplings[[i,1]];
               {particle, vectorBoson} = GetExternalStates[couplingSymbol];
               particlesStr = CConversion`ToValidCSymbolString[particle]
                              <> CConversion`ToValidCSymbolString[vectorBoson]
                              <> CConversion`ToValidCSymbolString[vectorBoson];
               dim = TreeMasses`GetDimension[particle];
               functionName = "get_" <> particlesStr <> "_partial_width(";
               If[dim == 1,
                  functionName = functionName <> ") const";,
                  functionName = functionName <> "int gO1) const";
                 ];
               couplingName = "eff_Cp" <> particlesStr;
               massStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]];
               body = "const double mass = PHYSICAL(" <> massStr <> ")";
               If[dim != 1,
                  body = body <> "(gO1)";
                 ];
               body = body <> ";\n";
               If[vectorBoson === SARAH`VectorP,
                  body = body <> "return " <> CConversion`RValueToCFormString[1 / (64 Pi)]
                         <> " * Power(mass, 3.0) * AbsSqr(" <> couplingName
                         <> If[dim != 1, "(gO1)", ""] <> ");";,
                  body = body <> "return " <> CConversion`RValueToCFormString[1 / (8 Pi)]
                         <> " * Power(mass, 3.0) * AbsSqr(" <> couplingName
                         <> If[dim != 1, "(gO1)", ""] <> ");";
                 ];
               functions = Append[functions,
                                  "double " <> FlexibleSUSY`FSModelName <> "_effective_couplings::"
                                  <> functionName <> "\n{\n" <> TextFormatting`IndentText[body]
                                  <> "\n}\n"];
               prototypes = Append[prototypes,
                                   "double " <> functionName <> ";"];
              ];
           Utils`StringJoinWithSeparator[#, "\n"]& /@ {prototypes, functions}
          ];

NonZeroVertexQ[vertex_] := MemberQ[vertex[[2 ;;]][[All, 1]], Except[0]];

(* @todo extend to multiple non-Abelian groups *)
IsColorOrLorentzIndex[index_] := Vertices`SarahColorIndexQ[index] ||
                                 Vertices`SarahLorentzIndexQ[index];
StripColorAndLorentzIndices[p_Symbol] := p;
StripColorAndLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripColorAndLorentzIndices[p]];
StripColorAndLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripColorAndLorentzIndices[p]];
StripColorAndLorentzIndices[p_] :=
    Module[{remainingIndices},
           remainingIndices = Select[p[[1]], (!IsColorOrLorentzIndex[#])&];
           If[Length[remainingIndices] === 0,
              Head[p],
              Head[p][remainingIndices]
             ]
          ];
SetAttributes[StripColorAndLorentzIndices, {Listable}];

(* @todo this is very slow because each possible vertex must be calculated,
   this can be improved by either pre-calculating all vertices or saving
   previous results (e.g. define along the lines of f[p] := f[p] = ...)    *)
GetTwoBodyDecays[particle_] :=
    Module[{i, j, allParticles, combinations, vertex, fields, couplings,
            candidate, found = {}, decays = {}},
           allParticles = Select[TreeMasses`GetParticles[], !TreeMasses`IsGhost[#]&];
           combinations = Table[Sort[{SARAH`AntiField[particle],
                                      SARAH`AntiField[allParticles[[i]]],
                                      allParticles[[j]]}],
                                {i, 1, Length[allParticles]}, {j, 1, Length[allParticles]}];
           combinations = DeleteDuplicates[Flatten[combinations, 1]];
           For[i = 1, i <= Length[combinations], i++,
               vertex = SARAH`Vertex[combinations[[i]], UseDependences -> True];
               If[NonZeroVertexQ[vertex],
                  fields = First[vertex];
                  coupling = Rest[vertex];
                  If[Length[coupling] > 1,
                     coupling = Vertices`SortCp[SARAH`Cp @@ (StripColorAndLorentzIndices @ fields)][SARAH`PL];,
                     coupling = Vertices`SortCp[SARAH`Cp @@ (StripColorAndLorentzIndices @ fields)];
                    ];
                  candidate = Append[DeleteCases[fields /. head_[{__}] :> head, p_ /; p === SARAH`AntiField[particle], {0, Infinity}, 1], coupling];
                  If[FreeQ[found, C[candidate[[1]], candidate[[2]]]] &&
                     FreeQ[found, C[SARAH`AntiField[candidate[[1]]], SARAH`AntiField[candidate[[2]]]]] ||
                     SARAH`AntiField[particle] =!= particle,
                     If[((!TreeMasses`IsMassless[candidate[[1]]] || !TreeMasses`IsVector[candidate[[1]]]) &&
                        (!TreeMasses`IsMassless[candidate[[2]]] || !TreeMasses`IsVector[candidate[[2]]])) ||
                        !FreeQ[SARAH`AllowDecaysMasslessVectors, SARAH`RE[particle]],
                        decays = Append[decays, candidate];
                        found = Append[found, C[candidate[[1]], candidate[[2]]]];
                       ];
                    ];
                 ];
              ];
           decays
          ];

GetParticlesCouplingToVectorBoson[vector_] :=
    Module[{i, charge, allParticles, particles = {}},
           allParticles = Select[TreeMasses`GetParticles[], !TreeMasses`IsGhost[#]&];
           For[i = 1, i <= Length[allParticles], i++,
               If[vector === SARAH`VectorG,
                  (* @note could use defined functions in e.g. TreeMasses, plus check
                     for undefined group factors, but will do it this way for now
                     to ensure consistency with SARAH                                  *)
                  charge = SARAH`Vertex[{SARAH`AntiField[allParticles[[i]]],
                                         allParticles[[i]], SARAH`VectorG},
                                        UseDependences -> True][[2,1]];
                  If[charge =!= 0,
                     particles = Append[particles, allParticles[[i]]];
                    ];,
                  charge = TreeMasses`GetElectricCharge[allParticles[[i]]];
                  If[NumericQ[charge],
                     charge = {charge},
                     charge = Cases[SARAH`Vertex[{SARAH`AntiField[allParticles[[i]]],
                                                  allParticles[[i]], SARAH`VectorP},
                                                 UseDependences -> True][[2,1]], _?NumberQ];
                    ];
                  If[charge =!= {} && SARAH`AntiField[allParticles[[i]]] =!= allParticles[[i]] &&
                     charge =!= {0},
                     particles = Append[particles, allParticles[[i]]];
                    ];
                 ];
              ];
           particles
          ];

IsMasslessOrGoldstone[SARAH`bar[p_]] := IsMasslessOrGoldstone[p];
IsMasslessOrGoldstone[Susyno`LieGroups`conj[p_]] := IsMasslessOrGoldstone[p];
IsMasslessOrGoldstone[particle_] :=
    Module[{result},
           result = TreeMasses`IsMassless[particle] ||
                    (TreeMasses`IsGoldstone[particle] && TreeMasses`GetDimension[particle] == 1);
           result
          ];

InitializeEffectiveCouplings[] :=
    Module[{i, couplings, particle, vectorBoson,
            allParticles = {}, allVectorBosons = {},
            twoBodyDecays, vectorBosonInteractions,
            neededTwoBodyDecays, neededVectorBosonInteractions,
            neededCoups, result = {}},
           couplings = GetAllowedCouplingsForModel[];
           {allParticles, allVectorBosons} = DeleteDuplicates /@ {(#[[1]])& /@ (GetExternalStates[#]& /@ couplings),
                                                                  (#[[2]])& /@ (GetExternalStates[#]& /@ couplings)};
           twoBodyDecays = {#, GetTwoBodyDecays[#]}& /@ allParticles;
           vectorBosonInteractions = {#, GetParticlesCouplingToVectorBoson[#]}& /@ allVectorBosons;
           For[i = 1, i <= Length[couplings], i++,
                {particle, vectorBoson} = GetExternalStates[couplings[[i]]];
                neededTwoBodyDecays = First[Select[twoBodyDecays, (#[[1]] === particle)&]];
                neededVectorBosonInteractions = First[Select[vectorBosonInteractions, (#[[1]] === vectorBoson)&]];
                neededCoups = Select[neededTwoBodyDecays[[2]],
                                    (MemberQ[neededVectorBosonInteractions[[2]], #[[1]]] ||
                                     MemberQ[neededVectorBosonInteractions[[2]], #[[2]]])&];
                (* only keep vertices of the form pD -> p SARAH`AntiField[p] *)
                neededCoups = Cases[neededCoups, {p1_, p2_, _} /; p1 === SARAH`AntiField[p2]];
                (* filter out massless states and Goldstones *)
                neededCoups = Select[neededCoups, (!IsMasslessOrGoldstone[#[[1]]] && !IsMasslessOrGoldstone[#[[2]]])&];
                result = Append[result, {couplings[[i]], #[[3]]& /@ neededCoups}];
              ];
           result
          ];

GetNeededVerticesList[couplings_List] :=
    {Null[Null, Join[(#[[2]])& /@ couplings]]};

CreateEffectiveCouplingName[pIn_, pOut_] :=
    "eff_Cp" <> CConversion`ToValidCSymbolString[pIn] <> CConversion`ToValidCSymbolString[pOut] <> CConversion`ToValidCSymbolString[pOut];

CreateEffectiveCouplingsGetters[couplings_List] :=
    Module[{i, couplingSymbols, type, particle,
            vectorBoson, dim, couplingName, getters = ""},
           couplingSymbols = #[[1]]& /@ couplings;
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           For[i = 1, i <= Length[couplingSymbols], i++,
               {particle, vectorBoson} = GetExternalStates[couplingSymbols[[i]]];
               dim = TreeMasses`GetDimension[particle];
               couplingName = CreateEffectiveCouplingName[particle, vectorBoson];
               getters = getters <> type <> " get_" <> couplingName;
               If[dim == 1,
                  getters = getters <> "() const { return " <> couplingName <> "; }\n";,
                  getters = getters <> "(int gO1) const { return " <> couplingName <> "(gO1); }\n";
                 ];
              ];
           getters
          ];

CreateEffectiveCouplingsDefinitions[couplings_List] :=
    Module[{i, couplingSymbols, dim, type, particle, vectorBoson,
            couplingName, defs = ""},
           couplingSymbols = #[[1]]& /@ couplings;
           For[i = 1, i <= Length[couplingSymbols], i++,
               {particle, vectorBoson} = GetExternalStates[couplingSymbols[[i]]];
               couplingName = CreateEffectiveCouplingName[particle, vectorBoson];
               dim = TreeMasses`GetDimension[particle];
               If[dim == 1,
                  type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];,
                  type = CConversion`CreateCType[CConversion`ArrayType[CConversion`complexScalarCType, dim]];
                 ];
               defs = defs <> type <> " " <> couplingName <> ";\n";
              ];
           defs
          ];

CreateEffectiveCouplingsInit[couplings_List] :=
    Module[{i, couplingSymbols, particle,
            vectorBoson, couplingName, dim, type, init = ""},
           couplingSymbols = #[[1]]& /@ couplings;
           For[i = 1, i <= Length[couplingSymbols], i++,
               {particle, vectorBoson} = GetExternalStates[couplingSymbols[[i]]];
               couplingName = CreateEffectiveCouplingName[particle, vectorBoson];
               dim = TreeMasses`GetDimension[particle];
               If[dim == 1,
                  type = CConversion`ScalarType[CConversion`complexScalarCType];,
                  type = CConversion`ArrayType[CConversion`complexScalarCType, dim];
                 ];
               init = init <> ", " <> CConversion`CreateDefaultConstructor[couplingName, type];
              ];
           init
          ];

RunToDecayingParticleScale[scale_] :=
    Module[{body, result = ""},
           If[SARAH`SupersymmetricModel,
              body = "model.run_to(" <> scale <> ");\n";
              result = "if (rg_improve && scale > " <> scale <> ") {\n"
                       <> TextFormatting`IndentText[body] <> "}\n";
             ];
           result
          ];

RunStrongCouplingToScale[scale_] :=
    Module[{result = ""},
           If[ValueQ[SARAH`hyperchargeCoupling] && ValueQ[SARAH`leftCoupling] &&
              ValueQ[SARAH`strongCoupling],
              result = result <> "run_SM_strong_coupling_to(sm, " <> scale <> ");\n";
             ];
           result
          ];

CallEffectiveCouplingCalculation[couplingSymbol_] :=
    Module[{particle, vectorBoson, savedMass, dim, start, idx = "",
            body, couplingName, call = ""},
           {particle, vectorBoson} = GetExternalStates[couplingSymbol];
           savedMass = CConversion`RValueToCFormString[FlexibleSUSY`M[particle]];
           dim = TreeMasses`GetDimension[particle];
           start = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
           If[dim != 1 && start <= dim,
              idx = "gO1";
              savedMass = savedMass <> "(" <> idx <> ")";
             ];
           If[vectorBoson === SARAH`VectorP,
              body = RunStrongCouplingToScale["0.5 * " <> savedMass];,
              body = RunStrongCouplingToScale[savedMass];
             ];
           call = call <> body;
           couplingName = CreateEffectiveCouplingName[particle, vectorBoson];
           call = call <> "calculate_" <> couplingName;
           If[idx != "",
              call = call <> "(" <> idx <> ");";,
              call = call <> "();";
             ];
           call
          ];

CreateEffectiveCouplingsCalculation[couplings_List] :=
    Module[{i, couplingSymbols, particle, couplingsForParticles = {},
            pos, couplingList, mass,
            savedMass, dim, start, body, result = ""},
           couplingSymbols = #[[1]]& /@ couplings;
           For[i = 1, i <= Length[couplingSymbols], i++,
               particle = GetExternalStates[couplingSymbols[[i]]][[1]];
               If[FreeQ[couplingsForParticles, particle],
                  couplingsForParticles = Append[couplingsForParticles, {particle, {couplingSymbols[[i]]}}];,
                  pos = Position[couplingsForParticles, {particle, _List}][[1,1]];
                  couplingList = couplingsForParticles[[pos]] /. {p_, coups_} :> {p, Append[coups, couplingSymbols[[i]]]};
                  couplingsForParticles = ReplacePart[couplingsForParticles, pos -> couplingList];
                 ];
              ];
           For[i = 1, i <= Length[couplingsForParticles], i++,
               particle = couplingsForParticles[[i,1]];
               mass = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]];
               savedMass = "const auto " <> mass <> " = PHYSICAL(" <> mass <> ");\n";
               dim = TreeMasses`GetDimension[particle];
               start = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
               If[dim == 1 && !TreeMasses`IsGoldstone[particle],
                  body = RunToDecayingParticleScale[mass];
                  If[SARAH`SupersymmetricModel,
                     body = body <> "model.calculate_DRbar_masses();\n"
                                 <> "copy_mixing_matrices_from_model();\n";
                    ];
                  result = result <> savedMass <> body <> Utils`StringJoinWithSeparator[CallEffectiveCouplingCalculation[#]& /@ couplingsForParticles[[i,2]], "\n"] <> "\n\n";
                  ,
                  If[start <= dim,
                     body = RunToDecayingParticleScale[mass <> "(gO1)"];
                     If[SARAH`SupersymmetricModel,
                        body = body <> "model.calculate_DRbar_masses();\n"
                                    <> "copy_mixing_matrices_from_model();\n";
                       ];
                     result = result <> savedMass
                                     <> "for (int gO1 = " <> ToString[start-1] <> "; gO1 < " <> ToString[dim] <> "; ++gO1) {\n";
                     body = body <> Utils`StringJoinWithSeparator[CallEffectiveCouplingCalculation[#]& /@ couplingsForParticles[[i,2]], "\n"] <> "\n";
                     result = result <> TextFormatting`IndentText[body] <> "}\n\n";
                    ];
                 ];
              ];

           result = "const double scale = model.get_scale();\nconst Eigen::ArrayXd saved_parameters(model.get());\n\n"
                    <> "const double saved_mt = PHYSICAL("
                    <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[]]]
                    <> ");\nPHYSICAL("
                    <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[]]]
                    <> ") = qedqcd.displayPoleMt();\n\n"
                    <> result;
           result = result <> "PHYSICAL("
                           <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[]]]
                           <> ") = saved_mt;\n";
           result = result <> "model.set_scale(scale);\nmodel.set(saved_parameters);\n";

           result
          ];

CreateEffectiveCouplingPrototype[coupling_] :=
    Module[{couplingSymbol = coupling[[1]], particle, vectorBoson,
            dim, name, result = ""},
           {particle, vectorBoson} = GetExternalStates[couplingSymbol];
           If[particle =!= Null && vectorBoson =!= Null,
              dim = TreeMasses`GetDimension[particle];
              name = CreateEffectiveCouplingName[particle, vectorBoson];
              result = "void calculate_" <> name <> If[dim == 1, "();\n", "(int gO1);\n"];
             ];
           result
          ];

GetEffectiveVEV[] :=
    Module[{vev, parameters = {}, result = ""},
           If[SARAH`SupersymmetricModel,
              vev = Simplify[2 Sqrt[-SARAH`Vertex[{SARAH`VectorW, Susyno`LieGroups`conj[SARAH`VectorW]}][[2,1]]
                           / SARAH`leftCoupling^2] /. SARAH`sum[a_,b_,c_,d_] :> Sum[d,{a,b,c}]];
              vev = Parameters`DecreaseIndexLiterals[vev];
              parameters = Parameters`FindAllParameters[vev];
              result = "const auto vev = " <> CConversion`RValueToCFormString[vev] <> ";\n";,
              result = "const auto vev = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));\n";
           ];
           {result, parameters}
          ];

GetMultiplicity[vectorBoson_, internal_] := SARAH`ChargeFactor[vectorBoson, internal, internal];

GetParticleGenerationIndex[particle_, coupling_] :=
    Module[{dim, indexList, result = {}},
           dim = TreeMasses`GetDimension[particle];
           If[dim != 1,
              indexList = Flatten[Cases[Vertices`GetParticleList[coupling],
                                  p_[a_List] /; p === particle :> a, {0, Infinity}]];
              result = Select[indexList, StringMatchQ[ToString[#], "gt" ~~ __] &];
              If[Length[result] > 1,
                 result = {result[[1]]};
                ];
             ];
           result
          ];

(* @todo these are basically identical to those in SelfEnergies,
   it would be better to reuse the definitions there if possible  *)
GetParticleIndicesInCoupling[SARAH`Cp[a__]] := Flatten[Cases[{a}, List[__], Infinity]];

GetParticleIndicesInCoupling[SARAH`Cp[a__][_]] := GetParticleIndicesInCoupling[SARAH`Cp[a]];

CreateCouplingSymbol[coupling_] :=
    Module[{symbol, indices},
           indices = GetParticleIndicesInCoupling[coupling];
           symbol = ToValidCSymbol[coupling /. a_[List[__]] :> a];
           symbol[Sequence @@ indices]
          ];

CreateLocalConstRefsIgnoringMixings[expr_, mixings_List] :=
    Module[{symbols, poleMasses},
           symbols = Parameters`FindAllParameters[expr];
           poleMasses = {
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_]]     /; MemberQ[Parameters`GetOutputParameters[],FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}],
               Cases[expr, FlexibleSUSY`Pole[FlexibleSUSY`M[a_[__]]] /; MemberQ[Parameters`GetOutputParameters[],FlexibleSUSY`M[a]] :> FlexibleSUSY`M[a], {0,Infinity}]
                        };
           symbols = DeleteDuplicates[Flatten[symbols]];
           symbols = Complement[symbols, mixings];
           Parameters`CreateLocalConstRefs[symbols]
          ];

CreateNeededCouplingFunction[coupling_, expr_, mixings_List] :=
    Module[{symbol, prototype = "", definition = "",
            indices = {}, localExpr, body = "", functionName = "", i,
            type, typeStr},
           indices = GetParticleIndicesInCoupling[coupling];
           symbol = CreateCouplingSymbol[coupling];
           functionName = CConversion`ToValidCSymbolString[CConversion`GetHead[symbol]];
           functionName = functionName <> "(";
           For[i = 1, i <= Length[indices], i++,
               If[i > 1, functionName = functionName <> ", ";];
               functionName = functionName <> "int ";
               If[!IntegerQ[indices[[i]]] && !FreeQ[expr, indices[[i]]],
                  functionName = functionName <> CConversion`ToValidCSymbolString[indices[[i]]];
                 ];
              ];
           functionName = functionName <> ")";
           If[Parameters`IsRealExpression[expr],
              type = CConversion`ScalarType[CConversion`realScalarCType];,
              type = CConversion`ScalarType[CConversion`complexScalarCType];];
           typeStr = CConversion`CreateCType[type];
           prototype = typeStr <> " " <> functionName <> " const;\n";
           definition = typeStr <> " " <> FlexibleSUSY`FSModelName
                        <> "_effective_couplings::" <> functionName <> " const\n{\n";
           localExpr = expr /. (Rule[#[],#]& /@ Parameters`GetDependenceSPhenoSymbols[]);
           localExpr = localExpr /. Parameters`GetDependenceSPhenoRules[];
           body = CreateLocalConstRefsIgnoringMixings[localExpr, mixings] <> "\n" <>
                  "const " <> typeStr <> " result = " <>
                  Parameters`ExpressionToString[localExpr] <> ";\n\n" <>
                  "return result;\n";
           body = TextFormatting`IndentText[TextFormatting`WrapLines[body]];
           definition = definition <> body <> "}\n";
           {prototype, definition}
          ];

CreateNeededVertexExpressions[vertexRules_List, mixings_List] :=
    Module[{k, prototypes = "", defs = "", coupling, expr,
            p, d},
           For[k = 1, k <= Length[vertexRules], k++,
               coupling = Vertices`ToCp[vertexRules[[k,1]]];
               expr = vertexRules[[k,2]];
               {p, d} = CreateNeededCouplingFunction[coupling, expr, mixings];
               prototypes = prototypes <> p;
               defs = defs <> d <> "\n";
              ];
           {prototypes, defs}
          ];

HasColorCharge[particle_] :=
    Module[{dynkin},
           dynkin = SA`Dynkin[particle,Position[SARAH`Gauge, SARAH`strongCoupling][[1,1]]];
           If[!NumericQ[dynkin], dynkin = 0];
           dynkin == 1/2
          ];

CreateCouplingContribution[particle_, vectorBoson_, coupling_] :=
    Module[{i, internal, particleIndex, indices, dim, start, factor, qcdfactor,
            mass, massStr, couplingSymbol, couplingName,
            scaleFunction, body = "", result = "", parameters = {}},
           internal = DeleteCases[Vertices`GetParticleList[coupling] /. field_[{__}] :> field,
                                  p_ /; p === particle, 1];
           internal = First[internal /. {SARAH`bar[p_] :> p, Susyno`LieGroups`conj[p_] :> p}];
           dim = TreeMasses`GetDimension[internal];
           mass = FlexibleSUSY`M[internal];
           massStr = CConversion`ToValidCSymbolString[mass];
           If[dim != 1,
              massStr = massStr <> "(gI1)";
             ];
           parameters = Append[parameters, mass];
           Which[TreeMasses`IsScalar[internal],
                 factor = 1/2;
                 If[particle === SARAH`HiggsBoson,
                    scaleFunction = "AS0";
                    If[vectorBoson === SARAH`VectorP && HasColorCharge[internal],
                       qcdfactor = "scalar_scalar_qcd_factor(decay_mass, " <> massStr <> ")";,
                       qcdfactor = "";
                      ];,
                    Return[{"",{}}];
                   ];,
                 TreeMasses`IsVector[internal],
                 factor = -1/2;
                 If[particle === SARAH`HiggsBoson,
                    scaleFunction = "AS1";
                    qcdfactor = "";,
                    Return[{"",{}}];
                   ];,
                 TreeMasses`IsFermion[internal],
                 factor = 1;
                 If[particle === SARAH`HiggsBoson,
                    scaleFunction = "AS12";
                    If[vectorBoson === SARAH`VectorP && HasColorCharge[internal],
                       qcdfactor = "scalar_fermion_qcd_factor(decay_mass, " <> massStr <> ")";,
                       qcdfactor = "";
                      ];,
                    scaleFunction = "AP12";
                    If[vectorBoson === SARAH`VectorP && HasColorCharge[internal],
                       qcdfactor = "pseudoscalar_fermion_qcd_factor(decay_mass, " <> massStr <> ")";,
                       qcdfactor = "";
                      ];
                   ];
                ];
           If[vectorBoson === SARAH`VectorP,
              factor = factor * TreeMasses`GetElectricCharge[internal]^2 GetMultiplicity[vectorBoson, internal];
             ];
           indices = GetParticleIndicesInCoupling[coupling];
           particleIndex = GetParticleGenerationIndex[particle, coupling];
           indices = Replace[indices, p_ /; !MemberQ[particleIndex, p] -> SARAH`gI1, 1];
           indices = Replace[indices, p_ /; MemberQ[particleIndex, p] -> SARAH`gO1, 1];
           couplingSymbol = CConversion`ToValidCSymbol[coupling /. a_[List[__]] :> a];
           couplingSymbol = couplingSymbol[Sequence @@ indices];
           couplingName = CConversion`ToValidCSymbolString[CConversion`GetHead[couplingSymbol]];
           couplingName = couplingName <> "(";
           For[i = 1, i <= Length[indices], i++,
               If[i > 1, couplingName = couplingName <> ", ";];
               If[!IntegerQ[indices[[i]]],
                  couplingName = couplingName <> CConversion`ToValidCSymbolString[indices[[i]]];
                 ];
              ];
           couplingName = couplingName <> ")";
           body = "result += " <> If[factor != 1, CConversion`RValueToCFormString[factor] <> " * ", ""]
                    <> If[qcdfactor != "", qcdfactor <> " * ", ""] <> couplingName
                    <> " * vev * " <> scaleFunction <> "(decay_scale / Sqr(" <> massStr <> ")) / "
                    <> If[IsFermion[internal], massStr, "Sqr(" <> massStr <> ")"] <> ";";
           If[dim == 1,
              result = body <> "\n";,
              start = TreeMasses`GetDimensionStartSkippingGoldstones[internal];
              result = "for (int gI1 = " <> ToString[start - 1] <> "; gI1 < " <> ToString[dim] <> "; ++gI1) {\n";
              result = result <> TextFormatting`IndentText[body] <> "\n}\n";
             ];
           {result, parameters}
          ];

CreateEffectiveCouplingFunction[coupling_] :=
    Module[{i, couplingSymbol = coupling[[1]], neededCouplings = coupling[[2]],
            particle, vectorBoson, dim, type, name, savedMass, mass, mixingSymbol,
            mixingName, parameters = {}, poleMasses, currentLine, body = "", result = ""},
           {particle, vectorBoson} = GetExternalStates[couplingSymbol];
           If[particle =!= Null && vectorBoson =!= Null,
              name = CreateEffectiveCouplingName[particle, vectorBoson];
              dim = TreeMasses`GetDimension[particle];
              result = result <> "void " <> FlexibleSUSY`FSModelName
                       <> "_effective_couplings::calculate_" <> name <> "(";
              If[dim == 1,
                 result = result <> ")\n{\n";,
                 result = result <> "int gO1)\n{\n";
                ];

              mass = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]];
              savedMass = "const auto decay_mass = PHYSICAL(" <> mass <> ")";
              If[dim == 1,
                 savedMass = savedMass <> ";\n";,
                 savedMass = savedMass <> "(gO1);\n";
                ];
              body = body <> savedMass;
              body = body <> "const auto decay_scale = 0.25 * Sqr(decay_mass);\n";
              (* use physical mixing matrices for decaying particle *)
              mixingSymbol = TreeMasses`FindMixingMatrixSymbolFor[particle];
              If[mixingSymbol =!= Null,
                 mixingName = CConversion`ToValidCSymbolString[mixingSymbol];
                 body = body <> "const auto saved_" <> mixingName <> " = " <> mixingName <> ";\n";
                 body = body <> mixingName <> " = PHYSICAL(" <> mixingName <> ");\n\n";,
                 body = body <> "\n";
                ];
              {currentLine, parameters} = {#[[1]], Join[parameters, #[[2]]]}& @ (GetEffectiveVEV[]);
              body = body <> currentLine <> "\n";
              body = body <> CConversion`CreateDefaultDefinition["result", CConversion`ScalarType[CConversion`complexScalarCType]] <> ";\n";

              For[i = 1, i <= Length[neededCouplings], i++,
                  {currentLine, parameters} = {#[[1]], Join[parameters, #[[2]]]}& @ (CreateCouplingContribution[particle, vectorBoson, neededCouplings[[i]]]);
                  body = body <> currentLine;
                 ];

              Which[particle === SARAH`HiggsBoson && vectorBoson === SARAH`VectorG,
                    body = body <> "result *= 0.75;\n\n";
                    body = body <> "if (include_qcd_corrections) {\n"
                           <> TextFormatting`IndentText["result *= scalar_scaling_factor(decay_mass);"] <> "\n}\n";,
                    particle === SARAH`PseudoScalar && vectorBoson === SARAH`VectorP,
                    body = body <> "result *= 2.0;\n";,
                    particle === SARAH`PseudoScalar && vectorBoson === SARAH`VectorG,
                    body = body <> "result *= 1.5;\n\n";
                    body = body <> "if (include_qcd_corrections) {\n"
                           <> TextFormatting`IndentText["result *= pseudoscalar_scaling_factor(decay_mass);"] <> "\n}\n";
                   ];

              If[vectorBoson === SARAH`VectorG,
                 parameters = Append[parameters, SARAH`strongCoupling];
                 body = "const double alpha_s = " <> CConversion`RValueToCFormString[SARAH`strongCoupling^2 / (4 Pi)]
                        <> ";\n" <> body;
                ];

              (* use pole masses in loop functions *)
              If[SARAH`SupersymmetricModel,
                 body = Parameters`CreateLocalConstRefs[DeleteDuplicates[parameters]] <> body <> "\n";,
                 poleMasses = DeleteDuplicates[Select[parameters, Parameters`IsOutputParameter]];
                 body = Parameters`CreateLocalConstRefsForPhysicalParameters[poleMasses] <> body <> "\n";
                 body = Parameters`CreateLocalConstRefs[DeleteDuplicates[Complement[parameters, poleMasses]]] <> body <> "\n";
                ];

              If[vectorBoson === SARAH`VectorP,
                 body = body <> "result *= "
                        <> CConversion`RValueToCFormString[1 / (2^(3/4) Pi)]
                        <> " * physical_input.get(Physical_input::alpha_em_0) * Sqrt(qedqcd.displayFermiConstant());\n\n";,
                 body = body <> "result *= "
                        <> CConversion`RValueToCFormString[2^(1/4) / (3 Pi)]
                        <> " * alpha_s * Sqrt(qedqcd.displayFermiConstant());\n\n";
                ];

              (* restore saved mixing *)
              If[mixingSymbol =!= Null,
                 body = body <> mixingName <> " = saved_" <> mixingName <> ";\n";
                ];
              body = body <> name <> If[dim != 1, "(gO1) = ", " = "] <> "result;\n";

              result = result <> TextFormatting`IndentText[TextFormatting`WrapLines[body]] <> "\n}\n";
             ];
           result
          ];

CreateEffectiveCouplingsPrototypes[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingPrototype[#])& /@ couplings;
           result
          ];

CreateEffectiveCouplingsFunctions[couplings_List] :=
    Module[{result = ""},
           (result = result <> CreateEffectiveCouplingFunction[#] <> "\n")& /@ couplings;
           result
          ];

CreateEffectiveCouplings[couplings_List, massMatrices_List, vertexRules_List] :=
    Module[{mixings, relevantVertexRules, verticesPrototypes,
            verticesFunctions,prototypes = "", functions = ""},
           mixings = Cases[Flatten[TreeMasses`GetMixingMatrixSymbol[#]& /@ massMatrices], Except[Null]];
           relevantVertexRules = Cases[vertexRules, r:(Rule[a_,b_] /; !FreeQ[couplings,a]) :> r];
           {verticesPrototypes, verticesFunctions} =
               CreateNeededVertexExpressions[relevantVertexRules, mixings];
           prototypes = prototypes <> verticesPrototypes
                        <> CreateEffectiveCouplingsPrototypes[couplings];
           functions = functions <> verticesFunctions
                       <> CreateEffectiveCouplingsFunctions[couplings];
           {prototypes, functions}
          ];

End[];

EndPackage[];
