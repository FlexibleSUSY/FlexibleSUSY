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

BeginPackage["WeinbergAngle`", {"SARAH`", "CConversion`", "Parameters`", "SelfEnergies`",
                                "TextFormatting`", "ThresholdCorrections`", "TreeMasses`",
                                "Utils`", "Vertices`"}];

CheckMuonDecayRunning::usage="";
InitMuonDecay::usage="";

DefSMhyperCoupling::usage="";
DefSMleftCoupling::usage="";
GetBottomMass::usage="";
GetTopMass::usage="";
DefVZSelfEnergy::usage="";
DefVWSelfEnergy::usage="";
YukawaMatching::usage="";
DefVZVWSelfEnergies::usage="";
DeltaAlphaHatBSM::usage="";

DeltaRhoHat2LoopSM::usage="";
DeltaRHat2LoopSM::usage="";
RhoHatTree::usage="";
DeltaVBwave::usage="";
DeltaVBvertex::usage="";
DeltaVBbox::usage="";
CreateDeltaVBContributions::usage="";
GetNeutrinoIndex::usage="";
CreateDeltaVBCalculation::usage="";

Begin["`Private`"];

MuonDecayWorks = True;

CheckMuonDecayRunning[] := MuonDecayWorks;

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<WeinbergAngle>: ", Sequence @@ Utils`InputFormOfNonStrings /@ {msg}]];

CheckMuonDecayInputRequirements[] :=
    Module[{requiredSymbols, availPars, areDefined},
           requiredSymbols = {SARAH`VectorP, SARAH`VectorW, SARAH`VectorZ,
                              SARAH`hyperchargeCoupling, SARAH`leftCoupling, SARAH`strongCoupling,
                              SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa};
           availPars = Join[TreeMasses`GetParticles[],
                            Parameters`GetInputParameters[],
                            Parameters`GetModelParameters[],
                            Parameters`GetOutputParameters[]];
           areDefined = MemberQ[availPars, #]& /@ requiredSymbols;
           DebugPrint["Error: Unknown symbol: ", #]& /@
              Cases[Utils`Zip[areDefined, requiredSymbols], {False, p_} :> p];
           And @@ areDefined
          ];

(*prepare usage of muon decay method*)
InitMuonDecay[eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Module[{},
           MuonDecayWorks = CheckMuonDecayInputRequirements[];
           (*enable usage of routine SARAH`InsFields*)
           SA`CurrentStates = eigenstates;
           SARAH`InitVertexCalculation[eigenstates, False];
           SARAH`ReadVertexList[eigenstates, False, False, True];
           SARAH`MakeCouplingLists;
          ];

DefSMhyperCoupling[] :=
    Module[{result},
           result = "const auto gY = ";
           If[!MuonDecayWorks,
              Return[result <> "1.;"]];
           result = result <> ThresholdCorrections`GetParameter[SARAH`hyperchargeCoupling] <> " * ";
           result = result <> FlexibleSUSY`FSModelName <> "_info::normalization_";
           result = result <> CConversion`ToValidCSymbolString[SARAH`hyperchargeCoupling] <> ";";
           result
          ];

DefSMleftCoupling[] :=
    Module[{result},
           result = "const auto g2 = ";
           If[!MuonDecayWorks,
              Return[result <> "1.;"]];
           result = result <> ThresholdCorrections`GetParameter[SARAH`leftCoupling] <> " * ";
           result = result <> FlexibleSUSY`FSModelName <> "_info::normalization_";
           result = result <> CConversion`ToValidCSymbolString[SARAH`leftCoupling] <> ";";
           result
          ];

GetWPlusBoson[] :=
    Switch[TreeMasses`GetElectricCharge[SARAH`VectorW],
            1, SARAH`VectorW,
           -1, Susyno`LieGroups`conj[SARAH`VectorW],
            _, Print["Error: W Boson has charge of neither +1 nor -1"]; Null
          ];

GetBottomMass[] := ThresholdCorrections`GetParameter[TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]]];

GetTopMass[] := ThresholdCorrections`GetParameter[TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]]];

DefVZSelfEnergy[] :=
    Module[{result},
           result = "const auto pizzt   = ";
           If[!MuonDecayWorks,
              Return[result <> "0.;"]];
           result <> "Re(model->" <> SelfEnergies`CreateSelfEnergyFunctionName[SARAH`VectorZ, 1] <> "(p));"
          ];

DefVWSelfEnergy[] :=
    Module[{result},
           result = "const auto piwwt   = ";
           If[!MuonDecayWorks,
              Return[result <> "0.;"]];
           result <> "Re(model->" <> SelfEnergies`CreateSelfEnergyFunctionName[SARAH`VectorW, 1] <> "(p));"
          ];

YukawaMatching[] :=
    Module[{fermion, yukawa, result},
           fermion = {TreeMasses`GetSMTopQuarkMultiplet[],
                      TreeMasses`GetSMBottomQuarkMultiplet[],
                      TreeMasses`GetSMTauLeptonMultiplet[]};
           yukawa = {SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa};
           If[(Parameters`GetParameterDimensions /@ yukawa) != {{3, 3}, {3, 3}, {3, 3}},
              MuonDecayWorks = False;
              DebugPrint["Error: Not all SM Yukawas are 3x3 matrices"];];
           prefactor = Table[ThresholdCorrections`YukawaToMassPrefactor[fermion[[i]], yukawa[[i]]], {i, 3}];
           If[!(And @@ Not /@ NumericQ /@ prefactor),
              MuonDecayWorks = False;
              DebugPrint["Error: Prefactors of SM Yukawas cannot be determined"];];
           If[!MuonDecayWorks,
              Return[""]];
           result = Parameters`CreateLocalConstRefs[prefactor] <> "\n";
           result = result <> "sm.set_Yu(Re(model->get_";
           result = result <> CConversion`ToValidCSymbolString[yukawa[[1]]] <> "()*";
           result = result <> CConversion`RValueToCFormString[prefactor[[1]]/(Global`SMvev/Sqrt[2])] <> "));\n";
           result = result <> "sm.set_Yd(Re(model->get_";
           result = result <> CConversion`ToValidCSymbolString[yukawa[[2]]] <> "()*";
           result = result <> CConversion`RValueToCFormString[prefactor[[2]]/(Global`SMvev/Sqrt[2])] <> "));\n";
           result = result <> "sm.set_Ye(Re(model->get_";
           result = result <> CConversion`ToValidCSymbolString[yukawa[[3]]] <> "()*";
           result <> CConversion`RValueToCFormString[prefactor[[3]]/(Global`SMvev/Sqrt[2])] <> "));"
          ];

DefVZVWSelfEnergies[] :=
    Module[{result},
           result = "const double sigma_Z_MZ_Model = Re(model->";
           result = result <> SelfEnergies`CreateSelfEnergyFunctionName[SARAH`VectorZ, 1] <> "(mz));\n";
           result = result <> "const double sigma_W_MW_Model = Re(model->";
           result = result <> SelfEnergies`CreateSelfEnergyFunctionName[SARAH`VectorW, 1] <> "(mw));\n";
           result = result <> "const double sigma_W_0_Model  = Re(model->";
           result <> SelfEnergies`CreateSelfEnergyFunctionName[SARAH`VectorW, 1] <> "(0.));"
          ];

DeltaAlphaHatBSM[scheme_] :=
    Module[{dahatbsm, result},
           dahatbsm = ThresholdCorrections`CalculateElectromagneticCoupling[scheme];
           result = Parameters`CreateLocalConstRefs[dahatbsm] <> "\n";
           result = result <> "delta_alpha_hat_bsm += alpha_em/(2.*Pi)*(";
           result <> CConversion`RValueToCFormString[dahatbsm] <> ");"
          ];

FindMassZ2[masses_List] := FindMass2[masses, SARAH`VectorZ];

FindMass2[masses_List, particle_] :=
    Module[{massExpr},
           massExpr = Cases[masses, TreeMasses`FSMassMatrix[{mass_}, particle, ___] :> mass];
           If[Head[massExpr] =!= List || massExpr === {},
              MuonDecayWorks = False;
              DebugPrint["Error: Could not find mass of ", particle, " in masses list"];
              Return[1]];
           TreeMasses`ReplaceDependenciesReverse[massExpr[[1]]]
          ];

(*extracts squared tree-level mass of Z boson before mixing with additional Z`*)
UnmixedZMass2[] :=
    Module[{ZMassMatrix, extraGaugeCouplings, submatrixList, submatrix, mass2Eigenvalues},
           ZMassMatrix = SARAH`MassMatrix[SARAH`VectorZ];
           Assert[MatrixQ[ZMassMatrix]];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] &&
                                          FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixList = ZMassMatrix[[#, #]] & /@
                              Flatten[Table[{i, j}, {i, 1, Length[ZMassMatrix]},
                                            {j, i + 1, Length[ZMassMatrix]}], 1];
           submatrix = Cases[submatrixList, x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings)];
           If[Length[submatrix] != 1,
              MuonDecayWorks = False;
              DebugPrint["Error: Photon-Z mass matrix could not be identified"];
              Return[1]];
           mass2Eigenvalues = Eigenvalues[submatrix];
           If[Length[mass2Eigenvalues] != 2 || !MemberQ[mass2Eigenvalues, 0],
              MuonDecayWorks = False;
              DebugPrint["Error: Determination of UnmixedZMass2 failed"];
              Return[1]];
           Select[mass2Eigenvalues, # =!= 0 &][[1]] /. Parameters`ApplyGUTNormalization[]
          ];

(*extracts squared tree-level mass of W boson before mixing with additional W`*)
UnmixedWMass2[] :=
    Module[{WMassMatrix, extraGaugeCouplings, submatrixList, submatrix, mass2Eigenvalues},
           WMassMatrix = SARAH`MassMatrix[SARAH`VectorW];
           Assert[MatrixQ[WMassMatrix]];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] &&
                                          FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixList = WMassMatrix[[#, #]] & /@
                              Flatten[Table[{i, j}, {i, 1, Length[WMassMatrix]},
                                            {j, i + 1, Length[WMassMatrix]}], 1];
           submatrix = Cases[submatrixList, x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings)];
           If[Length[submatrix] != 1,
              MuonDecayWorks = False;
              DebugPrint["Error: W mass matrix could not be identified"];
              Return[0]];
           mass2Eigenvalues = Eigenvalues[submatrix];
           If[Length[DeleteDuplicates[mass2Eigenvalues]] != 1,
              MuonDecayWorks = False;
              DebugPrint["Error: Determination of UnmixedWMass2 failed"];
              Return[0]];
           mass2Eigenvalues[[1]] /. Parameters`ApplyGUTNormalization[]
          ];

(*checks whether the Gell-Mann-Nishijima relation is valid*)
GellMannNishijimaRelationHolds[] :=
    Module[{photonMassMatrix, extraGaugeCouplings, submatrixIndices,
            BW3pos, photonEigenSystem, photonVector},
           If[FreeQ[SARAH`Gauge, SARAH`hypercharge] || FreeQ[SARAH`Gauge, SARAH`left],
              DebugPrint["Error: hypercharge or left gauge group does not exist"];
              Return[False]];
           photonMassMatrix = SARAH`MassMatrix[SARAH`VectorP];
           Assert[MatrixQ[photonMassMatrix]];
           If[Length[photonMassMatrix] > 4,
              DebugPrint["Error: neutral vector boson mass matrix is too large to be diagonalized"];
              Return[False]];
           extraGaugeCouplings = Cases[SARAH`Gauge, x_ /; FreeQ[x, SARAH`hypercharge] &&
                                          FreeQ[x, SARAH`left] && FreeQ[x, SARAH`color] :> x[[4]]];
           submatrixIndices = Flatten[Table[{i, j}, {i, 1, Length[photonMassMatrix]},
                                            {j, i + 1, Length[photonMassMatrix]}], 1];
           BW3pos = Flatten[Extract[submatrixIndices,
                                    Position[photonMassMatrix[[#, #]] & /@ submatrixIndices,
                                             x_ /; And @@ (FreeQ[x, #] & /@ extraGaugeCouplings),
                                             {1}, Heads -> False]]];
           If[Length[BW3pos] != 2,
              DebugPrint["Error: Photon-Z mass matrix could not be identified"];
              Return[False]];
           photonEigenSystem = Eigensystem[photonMassMatrix];
           photonVector = Extract[photonEigenSystem[[2]], Position[photonEigenSystem[[1]], 0]];
           If[!MemberQ[Total[Abs[Part[#, Complement[Range[Length[photonMassMatrix]],
                                                    BW3pos]] & /@ photonVector], {2}], 0],
              DebugPrint["Error: SM-like photon could not be identified"];
              Return[False]];
           True
          ];

(*calculates rho_0 from SU(2)_L representations of the Higgs multiplets as in (16) from 0801.1345[hep-ph]*)
RhoZero[] :=
    Module[{hyperchargePos, leftPos, vevlist},
           If[!GellMannNishijimaRelationHolds[],
              MuonDecayWorks = False;
              DebugPrint["Error: the Gell-Mann-Nishijima relation does not hold"];
              Return[1]];
           hyperchargePos = Position[SARAH`Gauge, x_ /; !FreeQ[x, SARAH`hypercharge], {1}][[1, 1]];
           leftPos = Position[SARAH`Gauge, x_ /; !FreeQ[x, SARAH`left], {1}][[1, 1]];
           vevlist = SARAH`DEFINITION[SARAH`EWSB][SARAH`VEVs];
           (*extract isospin from SU(2)_left representation and its third component from
                Gell-Mann-Nishijima formula with given hypercharge and electric charge = 0*)
           vevlist = vevlist /. {fieldname_Symbol, vevinfo_List, comp1_List, __} :>
                        Flatten[{vevinfo Boole[ReleaseHold[TreeMasses`GetElectricCharge[comp1[[1]]]] == 0],
                                 (SA`DimensionGG[fieldname, leftPos] - 1) / 2,
                                 -SA`ChargeGG[fieldname, hyperchargePos]}];
           If[!FreeQ[vevlist, None],
              MuonDecayWorks = False;
              DebugPrint["Error: determination of electric charge did not work"];
              Return[1]];
           Simplify[Plus @@ ((#[[3]]^2 - #[[4]]^2 + #[[3]]) Abs[#[[1]] #[[2]] Sqrt[2]]^2 & /@ vevlist) /
                       Plus @@ (2 #[[4]]^2 Abs[#[[1]] #[[2]] Sqrt[2]]^2 & /@ vevlist),
                    Element[Alternatives @@ Cases[SARAH`DEFINITION[SARAH`EWSB][SARAH`VEVs][[All, 2, 1]],
                                                  x_ /; Parameters`IsRealParameter[x]], Reals]]
          ];

ExpressWeinbergAngleInTermsOfGaugeCouplings[] :=
    Simplify[ArcCos[Sqrt[UnmixedWMass2[] / UnmixedZMass2[] / RhoZero[]]]];

(*parameters and functions transferred directly to and defined on the C++ level*)
extPars = {SINTHETAW, RHOHATRATIO, GFERMI, MW, MZ, MT, ALPHAS, RHO2, DELTARHAT1LOOP, PIZZTMZ};
Do[Format[extPars[[i]], CForm] = Format[ToString[extPars[[i]]], OutputForm], {i, Length[extPars]}];

(*returns coefficients of 1 and gamma5 in Higgs-top-top vertices*)
HiggsTopVertices[higgsName_] :=
    Module[{indexRange, indexList, topQuark, higgsVertices},
           If[FreeQ[TreeMasses`GetParticles[], higgsName] ||
                 TreeMasses`GetDimensionWithoutGoldstones[higgsName] == 0,
              Return[{}]];
           indexRange = TreeMasses`GetParticleIndices[higgsName][[All, 2]];
           If[indexRange === {},
              indexRange = {1}];
           indexList = Flatten[Table @@ {Table[ToExpression["i" <> ToString[k]], {k, Length[indexRange]}],
                                         Sequence @@ Table[{ToExpression["i" <> ToString[k]], 1,
                                                            indexRange[[k]]}, {k, Length[indexRange]}]},
                               Length[indexRange] - 1];
           topQuark = Level[TreeMasses`GetUpQuark[{3}], {Boole[ListQ[TreeMasses`GetUpQuark[{3}]]]}][[1]];
           higgsVertices =
              Vertices`StripGroupStructure[SARAH`Vertex[{bar[topQuark], topQuark, higgsName[#]}] & /@
                                              indexList, SARAH`ctNr /@ Range[4]];
           higgsVertices = Cases[higgsVertices,
                                 {{__, higgsField_}, {coeffPL_, SARAH`PL}, {coeffPR_, SARAH`PR}} /;
                                    !TreeMasses`IsGoldstone[higgsField] :>
                                 {higgsField /. List -> Sequence,
                                  Simplify[(coeffPR + coeffPL)/2], Simplify[(coeffPR - coeffPL)/2]}];
           higgsVertices
          ];

(*generalizes Higgs dependent part of (C.5) and (C.6) in hep-ph/9606211 analogous to (C.9) and (C.10)*)
HiggsContributions2LoopSM[] :=
    Module[{higgsVEV = TreeMasses`GetSMVEVExpr[Undefined], higgsDep},
           If[higgsVEV === Undefined,
              MuonDecayWorks = False;
              DebugPrint["Error: SM like Higgs vev does not exist"];
              Return[0]];
           higgsDep = (Abs[#[[2]]]^2 - Abs[#[[3]]]^2) RHO2[FlexibleSUSY`M[#[[1]]]/MT] &;
           Simplify[3 (GFERMI MT higgsVEV / (8 Pi^2 Sqrt[2]))^2 *
                    (Plus @@ (higgsDep /@ Join[HiggsTopVertices[SARAH`HiggsBoson],
                                               HiggsTopVertices[SARAH`PseudoScalar]]))]
          ];

(*formula according to (C.6) from hep-ph/9606211*)
DeltaRhoHat2LoopSM[]:=
    Module[{gY, alphaDRbar, expr, result},
           If[!MuonDecayWorks,
              Return[""]];
           gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
           alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
           expr = (alphaDRbar ALPHAS / (4 Pi^2 SINTHETAW^2) *
                      (-2.145 MT^2/MW^2 + 1.262 Log[MT/MZ] - 2.24 - 0.85 MZ^2/MT^2) +
                   HiggsContributions2LoopSM[]) /
                  (1 + PIZZTMZ/MZ^2);
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> "deltaRhoHat2LoopSM = " <> Parameters`ExpressionToString[expr] <> ";";
           result
          ];

(*formula according to (C.5) from hep-ph/9606211*)
DeltaRHat2LoopSM[]:=
    Module[{gY, alphaDRbar, expr, result},
           If[!MuonDecayWorks,
              Return[""]];
           gY = SARAH`hyperchargeCoupling FlexibleSUSY`GUTNormalization[SARAH`hyperchargeCoupling];
           alphaDRbar = gY^2 SARAH`leftCoupling^2 / (4 Pi (gY^2 + SARAH`leftCoupling^2));
           expr = alphaDRbar ALPHAS / (4 Pi^2 SINTHETAW^2 (1 - SINTHETAW^2)) *
                     (2.145 MT^2/MZ^2 + 0.575 Log[MT/MZ] - 0.224 - 0.144 MZ^2/MT^2) -
                  HiggsContributions2LoopSM[] (1 - DELTARHAT1LOOP) RHOHATRATIO;
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> "deltaRHat2LoopSM = " <> Parameters`ExpressionToString[expr] <> ";";
           result
          ];

(*calculates tree-level value of rhohat parameter from umixed and mixed Z mass as well as RhoZero*)
RhoHatTree[]:=
    Module[{Zmass2unmixed, Zmass2mixed, expr, result},
           If[!MuonDecayWorks,
              Return[""]];
           Zmass2unmixed = UnmixedZMass2[];
           Zmass2mixed = FindMassZ2[TreeMasses`GetUnmixedParticleMasses[] /.
                                       Parameters`ApplyGUTNormalization[]];
           expr = Simplify[RhoZero[] Zmass2unmixed / Zmass2mixed /.
                              SARAH`Weinberg -> ExpressWeinbergAngleInTermsOfGaugeCouplings[],
                           SARAH`hyperchargeCoupling > 0 && SARAH`leftCoupling > 0];
           result = Parameters`CreateLocalConstRefs[expr] <> "\n";
           result = result <> "rhohat_tree = ";
           result = result <> CConversion`RValueToCFormString[expr] <> ";";
           result
          ];


(*functions for creation of wave-function renormalization, vertex and box corrections:*)

(*excludes diagrams in which an internal particle fulfills the condition excludeif*)
ExcludeDiagrams[diagrs_List, excludeif_:(False &)] :=
    Select[diagrs, !Or @@ (excludeif /@ (Cases[#, Rule[Internal[_], x_] :> x, Infinity])) &];

(*generates wave-function renormalization diagrams for given particle*)
GenerateDiagramsWave[particle_] :=
    Module[{couplings, insertrules, diagrs},
           couplings = {C[SARAH`External[1], SARAH`Internal[1], SARAH`AntiField[SARAH`Internal[2]]]};
           insertrules = {SARAH`External[1] -> particle, SARAH`Internal[1] -> SARAH`FieldToInsert[1],
                          SARAH`Internal[2] -> SARAH`FieldToInsert[2]};
           diagrs = SARAH`InsFields[{couplings /. insertrules, insertrules}];
           (*add indices for later summation*)
           diagrs = diagrs /. (Rule[SARAH`Internal[i_], x_] /; TreeMasses`GetDimension[x] > 1) :>
                                 Rule[SARAH`Internal[i], x[{ToExpression["SARAH`gI" <> ToString[i]]}]];
           diagrs = diagrs /. (Rule[SARAH`External[i_], x_] /; TreeMasses`GetDimension[x] > 1) :>
                                 Rule[SARAH`External[i], x[{ToExpression["SARAH`gO" <> ToString[i]]}]];
           diagrs = ({couplings /. #[[2]], #[[2]]}) & /@ diagrs;
           diagrs
          ];

(*calculates contribution from given wave-function renormalization diagram*)
WaveResult[diagr_List, includeGoldstones_] :=
    Module[{coupl, intparticles, intfermion, intscalar, result, intpartwithindex},
           coupl = (diagr[[1, 1]] /. C[a__] -> SARAH`Cp[a])[SARAH`PL];
           intparticles = ({SARAH`Internal[1], SARAH`Internal[2]} /. diagr[[2]]) /.
                             {SARAH`bar[p_] :> p, Susyno`LieGroups`conj[p_] :> p};
           If[Select[intparticles, TreeMasses`IsFermion] === {},
              Print["Warning: no internal fermion in wave function diagram"];
              Return[0]];
           intfermion = Select[intparticles, TreeMasses`IsFermion][[1]];
           If[Select[intparticles, TreeMasses`IsScalar] === {},
              Print["Warning: no internal scalar in wave function diagram"];
              Return[0]];
           intscalar = Select[intparticles, TreeMasses`IsScalar][[1]];
           result = -coupl Susyno`LieGroups`conj[coupl] *
                    SARAH`B1[0, SARAH`Mass2[intfermion], SARAH`Mass2[intscalar]];
           (*add sums over internal particles*)
           intpartwithindex = Reverse[Cases[intparticles, _[{_}]]];
           Do[result = FlexibleSUSY`SUM[
                          intpartwithindex[[i, 1, 1]],
                          If[includeGoldstones, 0,
                             TreeMasses`GetDimensionStartSkippingGoldstones[intpartwithindex[[i]]] - 1],
                          TreeMasses`GetDimension[intpartwithindex[[i]]] - 1,
                          result],
              {i, Length[intpartwithindex]}];
           result
          ];

(*combines generation of diagrams and calculation of their contributions*)
CompleteWaveResult[particle_, includeGoldstones_] :=
    Plus @@ (WaveResult[#, includeGoldstones] &) /@
       ExcludeDiagrams[GenerateDiagramsWave[particle],
                       If[includeGoldstones, TreeMasses`IsVector,
                          TreeMasses`IsVector[#] || TreeMasses`IsGoldstone[#] &]];

(*returns the complete wave-function renormalization part of deltaVB*)
DeltaVBwave[includeGoldstones_:False] :=
    Module[{neutrinofields, neutrinoresult, chargedleptonfields, chargedleptonresult},
           If[!MuonDecayWorks,
              Return[{}]];
           neutrinofields = TreeMasses`GetSMNeutralLeptons[];
           If[Length[neutrinofields] == 1,
              If[TreeMasses`GetDimension[neutrinofields[[1]]] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBwave does not work since there are not 3 neutrinos"];
                 Return[{}]];
              neutrinoresult =
                 {WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {SARAH`gO1}, neutrinofields[[1]]},
                                        CompleteWaveResult[neutrinofields[[1]], includeGoldstones]]},
              If[Length[neutrinofields] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBwave does not work since there are ",
                            "neither 1 nor 3 neutrino fields"];
                 Return[{}]];
              If[TreeMasses`GetDimension[neutrinofields[[1]]] != 1 ||
                    TreeMasses`GetDimension[neutrinofields[[2]]] != 1,
                 MuonDecayWorks = False;
                 DebugPrint["Error: definition of neutrino fields not supported"];
                 Return[{}]];
              neutrinoresult =
                 {WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {}, neutrinofields[[1]]},
                                        CompleteWaveResult[neutrinofields[[1]], includeGoldstones]],
                  WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {}, neutrinofields[[2]]},
                                        CompleteWaveResult[neutrinofields[[2]], includeGoldstones]]}];
           chargedleptonfields = TreeMasses`GetSMChargedLeptons[];
           If[Length[chargedleptonfields] == 1,
              If[TreeMasses`GetDimension[chargedleptonfields[[1]]] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBwave does not work since there are not 3 charged leptons"];
                 Return[{}]];
              chargedleptonresult =
                 {WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {SARAH`gO1}, chargedleptonfields[[1]]},
                                        CompleteWaveResult[chargedleptonfields[[1]], includeGoldstones]]},
              If[Length[chargedleptonfields] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBwave does not work since there are ",
                            "neither 1 nor 3 charged lepton fields"];
                 Return[{}]];
              If[TreeMasses`GetDimension[chargedleptonfields[[1]]] != 1 ||
                    TreeMasses`GetDimension[chargedleptonfields[[2]]] != 1,
                 MuonDecayWorks = False;
                 DebugPrint["Error: definition of charged lepton fields not supported"];
                 Return[{}]];
              chargedleptonresult =
                 {WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {}, chargedleptonfields[[1]]},
                                        CompleteWaveResult[chargedleptonfields[[1]], includeGoldstones]],
                  WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {}, chargedleptonfields[[2]]},
                                        CompleteWaveResult[chargedleptonfields[[2]], includeGoldstones]]}];
           Join[neutrinoresult, chargedleptonresult]
          ];

(*generates vertex diagrams for given external particles*)
GenerateDiagramsVertex[part1_, part2_, part3_] :=
    Module[{couplings, insertrules, diagrs},
           couplings = {C[SARAH`External[1], SARAH`AntiField[SARAH`Internal[2]], SARAH`Internal[3]],
                        C[SARAH`External[2], SARAH`Internal[1], SARAH`AntiField[SARAH`Internal[3]]],
                        C[SARAH`External[3], SARAH`AntiField[SARAH`Internal[1]], SARAH`Internal[2]]};
           insertrules =
              {SARAH`External[1] -> part1, SARAH`External[2] -> part2, SARAH`External[3] -> part3,
               SARAH`Internal[1] -> SARAH`FieldToInsert[1], SARAH`Internal[2] -> SARAH`FieldToInsert[2],
               SARAH`Internal[3] -> SARAH`FieldToInsert[3]};
           diagrs = SARAH`InsFields[{couplings /. insertrules, insertrules}];
           (*add indices for later summation*)
           diagrs = diagrs /. (Rule[SARAH`Internal[i_], x_] /; TreeMasses`GetDimension[x] > 1) :>
                                 Rule[SARAH`Internal[i], x[{ToExpression["SARAH`gI" <> ToString[i]]}]];
           diagrs = diagrs /. (Rule[SARAH`External[i_], x_] /; TreeMasses`GetDimension[x] > 1) :>
                                 Rule[SARAH`External[i], x[{ToExpression["SARAH`gO" <> ToString[i]]}]];
           diagrs = ({couplings /. #[[2]], #[[2]]}) & /@ diagrs;
           diagrs
          ];

(*True for Majorana fermions and outgoing Dirac fermions*)
IsOutgoingFermion[particle_] := TreeMasses`IsFermion[particle] &&
                                   (!FreeQ[particle, SARAH`bar] || SARAH`AntiField[particle] === particle);

(*True for Majorana Fermions and incoming Dirac fermions*)
IsIncomingFermion[particle_] := TreeMasses`IsFermion[particle] && FreeQ[particle, SARAH`bar];

(*calculates contribution from given vertex diagram with 1 internal fermion and 2 internal scalars*)
VertexResultFSS[diagr_List, includeGoldstones_] :=
    Module[{extparticles, extvectorindex, extoutindex, extinindex, inscalar, outscalar, couplSSV,
            scalarsoutin, factor, couplFFSout, couplFFSin, intparticles, intfermion, intscalars,
            result, intpartwithindex},
           extparticles = {SARAH`External[1], SARAH`External[2], SARAH`External[3]} /. diagr[[2]];
           extvectorindex = Position[extparticles, x_ /; TreeMasses`IsVector[x],
                                     {1}, Heads -> False][[1, 1]];
           extoutindex = Position[extparticles, x_ /; IsOutgoingFermion[x], {1}, Heads -> False][[1, 1]];
           extinindex = Complement[{1, 2, 3}, {extvectorindex, extoutindex}][[1]];
           (*scalars with incoming and outgoing momentum at SSV vertex*)
           inscalar = SARAH`AntiField[Select[List @@ diagr[[1, extinindex]], TreeMasses`IsScalar][[1]]];
           outscalar = SARAH`AntiField[Select[List @@ diagr[[1, extoutindex]], TreeMasses`IsScalar][[1]]];
           couplSSV = diagr[[1, extvectorindex]] /. C[a__] -> SARAH`Cp[a];
           scalarsoutin = List @@ Take[Vertices`SortCp[couplSSV], 2];
           (*momentum direction at SSV vertex has to be correct*)
           (*coupling convention assumes first scalar to be outgoing and second one to be incoming*)
           If[scalarsoutin === {outscalar, inscalar},
              factor = 1,
              If[scalarsoutin === {inscalar, outscalar},
                 factor = -1,
                 Print["Warning: scalar direction could not be determined"];
                 Return[0]]];
           couplSSV = factor couplSSV;
           couplFFSout = (diagr[[1, extoutindex]] /. C[a__] -> SARAH`Cp[a])[SARAH`PR];
           couplFFSin = (diagr[[1, extinindex]] /. C[a__] -> SARAH`Cp[a])[SARAH`PL];
           intparticles = ({SARAH`Internal[1], SARAH`Internal[2], SARAH`Internal[3]} /. diagr[[2]]) /.
                             {SARAH`bar[p_] :> p, Susyno`LieGroups`conj[p_] :> p};
           intfermion = Select[intparticles, TreeMasses`IsFermion][[1]];
           intscalars = Select[intparticles, TreeMasses`IsScalar];
           result = 1/2 couplSSV couplFFSout couplFFSin *
                    (1/2 + SARAH`B0[0, SARAH`Mass2[intscalars[[1]]], SARAH`Mass2[intscalars[[2]]]] +
                     SARAH`Mass2[intfermion] SARAH`C0[SARAH`Mass2[intfermion], SARAH`Mass2[intscalars[[1]]],
                                                      SARAH`Mass2[intscalars[[2]]]]);
           (*add sums over internal particles*)
           intpartwithindex = Reverse[Cases[intparticles, _[{_}]]];
           Do[result = FlexibleSUSY`SUM[
                          intpartwithindex[[i, 1, 1]],
                          If[includeGoldstones, 0,
                             TreeMasses`GetDimensionStartSkippingGoldstones[intpartwithindex[[i]]] - 1],
                          TreeMasses`GetDimension[intpartwithindex[[i]]] - 1,
                          result],
              {i, Length[intpartwithindex]}];
           result
          ];

(*calculates contribution from given vertex diagram with 2 internal fermions and 1 internal scalar*)
VertexResultFFS[diagr_List, includeGoldstones_] :=
    Module[{extparticles, extvectorindex, extoutindex, extinindex, fermiondirectok1, fermiondirectok2,
            needfermionflip, innaturalorder, orderedparticles, couplFFVPL, couplFFVPR,
            couplFFSout, couplFFSin, intparticles, intfermions, intscalar, result, intpartwithindex},
           extparticles = {SARAH`External[1], SARAH`External[2], SARAH`External[3]} /. diagr[[2]];
           extvectorindex = Position[extparticles, x_ /; TreeMasses`IsVector[x],
                                     {1}, Heads -> False][[1, 1]];
           extoutindex = Position[extparticles, x_ /; IsOutgoingFermion[x], {1}, Heads -> False][[1, 1]];
           extinindex = Complement[{1, 2, 3}, {extvectorindex, extoutindex}][[1]];
           (*is fermion flip necessary?*)
           fermiondirectok1 =
              Or @@ IsIncomingFermion /@
                 Complement[List @@ diagr[[1, extoutindex]], {SARAH`External[extoutindex]} /. diagr[[2]]];
           fermiondirectok2 =
              Or @@ IsOutgoingFermion /@
                 Complement[List @@ diagr[[1, extinindex]], {SARAH`External[extinindex]} /. diagr[[2]]];
           needfermionflip = !(fermiondirectok1 && fermiondirectok2);
           (*are fermions in natural order (outgoing before incoming)?*)
           innaturalorder = (IsOutgoingFermion[#[[1]]] && IsIncomingFermion[#[[2]]]) & @
                               Select[List @@ diagr[[1, extvectorindex]], TreeMasses`IsFermion];
           orderedparticles = If[innaturalorder,
                                 List @@ diagr[[1, extvectorindex]],
                                 Reverse[List @@ diagr[[1, extvectorindex]]]];
           (*use non-flipped or flipped FFV vertex appropriately*)
           couplFFVPL = (SARAH`Cp @@ If[needfermionflip,
                                        Reverse[orderedparticles],
                                        orderedparticles])[SARAH`PL];
           couplFFVPR = (SARAH`Cp @@ If[needfermionflip,
                                        Reverse[orderedparticles],
                                        orderedparticles])[SARAH`PR];
           couplFFSout = (diagr[[1, extoutindex]] /. C[a__] -> SARAH`Cp[a])[SARAH`PR];
           couplFFSin = (diagr[[1, extinindex]] /. C[a__] -> SARAH`Cp[a])[SARAH`PL];
           intparticles = ({SARAH`Internal[1], SARAH`Internal[2], SARAH`Internal[3]} /. diagr[[2]]) /.
                             {SARAH`bar[p_] :> p, Susyno`LieGroups`conj[p_] :> p};
           intfermions = Select[intparticles, TreeMasses`IsFermion];
           intscalar = Select[intparticles, TreeMasses`IsScalar][[1]];
           result = couplFFSout couplFFSin *
                    (-couplFFVPL FlexibleSUSY`M[intfermions[[1]]] FlexibleSUSY`M[intfermions[[2]]] *
                        SARAH`C0[SARAH`Mass2[intscalar], SARAH`Mass2[intfermions[[1]]],
                                 SARAH`Mass2[intfermions[[2]]]] +
                     1/2 couplFFVPR *
                        (-1/2 + SARAH`B0[0, SARAH`Mass2[intfermions[[1]]], SARAH`Mass2[intfermions[[2]]]] +
                         SARAH`Mass2[intscalar] SARAH`C0[SARAH`Mass2[intscalar],
                                                         SARAH`Mass2[intfermions[[1]]],
                                                         SARAH`Mass2[intfermions[[2]]]]));
           (*add sums over internal particles*)
           intpartwithindex = Reverse[Cases[intparticles, _[{_}]]];
           Do[result = FlexibleSUSY`SUM[
                          intpartwithindex[[i, 1, 1]],
                          If[includeGoldstones, 0,
                             TreeMasses`GetDimensionStartSkippingGoldstones[intpartwithindex[[i]]] - 1],
                          TreeMasses`GetDimension[intpartwithindex[[i]]] - 1,
                          result],
              {i, Length[intpartwithindex]}];
           result
          ];

(*calculates contribution from given vertex diagram*)
VertexResult[diagr_List, includeGoldstones_] :=
    Module[{intparticles, nFermions, nScalars},
           intparticles = {SARAH`Internal[1], SARAH`Internal[2], SARAH`Internal[3]} /. diagr[[2]];
           nFermions = Count[TreeMasses`IsFermion /@ intparticles, True];
           nScalars = Count[TreeMasses`IsScalar /@ intparticles, True];
           Switch[{nFermions, nScalars},
                  {1, 2}, VertexResultFSS[diagr, includeGoldstones],
                  {2, 1}, VertexResultFFS[diagr, includeGoldstones],
                  _, Print["Warning: vertex diagram type not supported"]; 0]
          ];

(*calculates tree-level vertex result for normalization of one-loop results*)
VertexTreeResult[part1_, part2_] :=
    Module[{part1withindex, part2withindex},
           If[TreeMasses`GetDimension[part1] > 1,
              part1withindex = part1[{SARAH`gO1}],
              part1withindex = part1];
           If[TreeMasses`GetDimension[part2] > 1,
              part2withindex = part2[{SARAH`gO2}],
              part2withindex = part2];
           If[TreeMasses`GetDimension[SARAH`VectorW] != 1,
              MuonDecayWorks = False;
              DebugPrint["Error: W boson is not defined uniquely or at all"];
              Return[1]];
           SARAH`Cp[part2withindex, part1withindex, GetWPlusBoson[]][SARAH`PL]
          ];

(*combines generation of diagrams and calculation of their contributions*)
CompleteVertexResult[part1_, part2_, includeGoldstones_] :=
    HoldForm[Evaluate[
       (Plus @@ (VertexResult[#, includeGoldstones] &) /@
          ExcludeDiagrams[GenerateDiagramsVertex[part1, part2, GetWPlusBoson[]],
                          If[includeGoldstones, TreeMasses`IsVector,
                             TreeMasses`IsVector[#] || TreeMasses`IsGoldstone[#] &]])]] /
    VertexTreeResult[part1, part2];

(*returns the complete vertex part of deltaVB*)
DeltaVBvertex[includeGoldstones_:False] :=
    Module[{neutrinofields, chargedleptonfields, result},
           If[!MuonDecayWorks,
              Return[{}]];
           neutrinofields = TreeMasses`GetSMNeutralLeptons[];
           chargedleptonfields = TreeMasses`GetSMChargedLeptons[];
           If[Length[neutrinofields] != Length[chargedleptonfields],
              MuonDecayWorks = False;
              DebugPrint["Error: DeltaVBvertex does not work since the numbers of ",
                         "neutrino and charged lepton fields are different"];
              Return[{}]];
           If[Length[neutrinofields] == 1,
              If[TreeMasses`GetDimension[neutrinofields[[1]]] != 3 ||
                    TreeMasses`GetDimension[chargedleptonfields[[1]]] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBvertex does not work since there are ",
                            "not 3 neutrinos and 3 charged leptons"];
                 Return[{}]];
              result = {WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, {SARAH`gO1, SARAH`gO2}},
                                              CompleteVertexResult[chargedleptonfields[[1]],
                                                                   SARAH`bar[neutrinofields[[1]]],
                                                                   includeGoldstones]]},
              If[Length[neutrinofields] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBvertex does not work since there are ",
                            "neither 1 nor 3 neutrino fields"];
                 Return[{}]];
              If[Or @@ ((TreeMasses`GetDimension[#] != 1) & /@
                           {neutrinofields[[1]], neutrinofields[[2]],
                            chargedleptonfields[[1]], chargedleptonfields[[2]]}),
                 MuonDecayWorks = False;
                 DebugPrint["Error: definition of neutrino or charged lepton fields not supported"];
                 Return[{}]];
              result = {WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, {}, chargedleptonfields[[1]],
                                               neutrinofields[[1]]},
                                              CompleteVertexResult[chargedleptonfields[[1]],
                                                                   SARAH`bar[neutrinofields[[1]]],
                                                                   includeGoldstones]],
                        WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, {}, chargedleptonfields[[2]],
                                               neutrinofields[[2]]},
                                              CompleteVertexResult[chargedleptonfields[[2]],
                                                                   SARAH`bar[neutrinofields[[2]]],
                                                                   includeGoldstones]]}];
           result
          ];

(*generates box diagrams for given external particles*)
GenerateDiagramsBox[part1_, part2_, part3_, part4_] :=
    Module[{couplings1, couplings2, couplings3, insertrules, diagrs1, diagrs2, diagrs3},
           couplings1 = {C[SARAH`External[1], SARAH`Internal[4], SARAH`AntiField[SARAH`Internal[1]]],
                         C[SARAH`External[2], SARAH`Internal[1], SARAH`AntiField[SARAH`Internal[2]]],
                         C[SARAH`External[3], SARAH`Internal[2], SARAH`AntiField[SARAH`Internal[3]]],
                         C[SARAH`External[4], SARAH`Internal[3], SARAH`AntiField[SARAH`Internal[4]]]};
           couplings2 = {C[SARAH`External[1], SARAH`Internal[4], SARAH`AntiField[SARAH`Internal[1]]],
                         C[SARAH`External[2], SARAH`Internal[1], SARAH`AntiField[SARAH`Internal[2]]],
                         C[SARAH`External[3], SARAH`Internal[3], SARAH`AntiField[SARAH`Internal[4]]],
                         C[SARAH`External[4], SARAH`Internal[2], SARAH`AntiField[SARAH`Internal[3]]]};
           couplings3 = {C[SARAH`External[1], SARAH`Internal[4], SARAH`AntiField[SARAH`Internal[1]]],
                         C[SARAH`External[2], SARAH`Internal[2], SARAH`AntiField[SARAH`Internal[3]]],
                         C[SARAH`External[3], SARAH`Internal[1], SARAH`AntiField[SARAH`Internal[2]]],
                         C[SARAH`External[4], SARAH`Internal[3], SARAH`AntiField[SARAH`Internal[4]]]};
           insertrules =
              {SARAH`External[1] -> part1, SARAH`External[2] -> part2,
               SARAH`External[3] -> part3, SARAH`External[4] -> part4,
               SARAH`Internal[1] -> SARAH`FieldToInsert[1], SARAH`Internal[2] -> SARAH`FieldToInsert[2],
               SARAH`Internal[3] -> SARAH`FieldToInsert[3], SARAH`Internal[4] -> SARAH`FieldToInsert[4]};
           diagrs1 = SARAH`InsFields[{couplings1 /. insertrules, insertrules}];
           diagrs2 = SARAH`InsFields[{couplings2 /. insertrules, insertrules}];
           diagrs3 = SARAH`InsFields[{couplings3 /. insertrules, insertrules}];
           (*add indices for later summation*)
           {diagrs1, diagrs2, diagrs3} =
              {diagrs1, diagrs2, diagrs3} /.
              (Rule[SARAH`Internal[i_], x_] /; TreeMasses`GetDimension[x] > 1) :>
                 Rule[SARAH`Internal[i], x[{ToExpression["SARAH`gI" <> ToString[i]]}]];
           {diagrs1, diagrs2, diagrs3} =
              {diagrs1, diagrs2, diagrs3} /.
              (Rule[SARAH`External[i_], x_] /; TreeMasses`GetDimension[x] > 1) :>
                 Rule[SARAH`External[i], x[{ToExpression["SARAH`gO" <> ToString[i]]}]];
           (*add topoNr to distinguish different topologies -> appropriate result can later be calculated*)
           diagrs1 = ({couplings1 /. #[[2]], Append[#[[2]], WeinbergAngle`topoNr -> 1]}) & /@ diagrs1;
           diagrs2 = ({couplings2 /. #[[2]], Append[#[[2]], WeinbergAngle`topoNr -> 2]}) & /@ diagrs2;
           diagrs3 = ({couplings3 /. #[[2]], Append[#[[2]], WeinbergAngle`topoNr -> 3]}) & /@ diagrs3;
           Join[diagrs1, diagrs2, diagrs3]
          ];

(*calculates contribution from given box diagram*)
BoxResult[diagr_List, includeGoldstones_] :=
    Module[{couplMu, couplMuNeutr, couplElNeutr, couplEl,
            intparticles, intfermions, toponr, result, intpartwithindex},
           couplMu = (diagr[[1, 1]] /. C[a__] -> SARAH`Cp[a])[SARAH`PL];
           couplMuNeutr = (diagr[[1, 2]] /. C[a__] -> SARAH`Cp[a])[SARAH`PR];
           couplElNeutr = (diagr[[1, 3]] /. C[a__] -> SARAH`Cp[a])[SARAH`PL];
           couplEl = (diagr[[1, 4]] /. C[a__] -> SARAH`Cp[a])[SARAH`PR];
           intparticles = ({SARAH`Internal[1], SARAH`Internal[2],
                            SARAH`Internal[3], SARAH`Internal[4]} /. diagr[[2]]) /.
                             {SARAH`bar[p_] :> p, Susyno`LieGroups`conj[p_] :> p};
           If[Length[Select[intparticles, TreeMasses`IsFermion]] != 2,
              Print["Warning: not 2 internal fermions in box diagram"];
              Return[0]];
           If[Length[Select[intparticles, TreeMasses`IsScalar]] != 2,
              Print["Warning: not 2 internal scalars in box diagram"];
              Return[0]];
           intfermions = Select[intparticles, TreeMasses`IsFermion];
           toponr = WeinbergAngle`topoNr /. diagr[[2]];
           result = couplMu couplMuNeutr couplElNeutr couplEl;
           If[toponr == 1,
              result = result * SARAH`D27[Sequence @@ SARAH`Mass2 /@ intparticles]];
           If[toponr == 2 && TreeMasses`IsFermion[intparticles[[1]]],
              result = result * (-1) * SARAH`D27[Sequence @@ SARAH`Mass2 /@ intparticles]];
           If[toponr == 2 && TreeMasses`IsScalar[intparticles[[1]]],
              result = result * 1/2 * FlexibleSUSY`M[intfermions[[1]]] FlexibleSUSY`M[intfermions[[2]]] *
                       SARAH`D0[Sequence @@ SARAH`Mass2 /@ intparticles]];
           If[toponr == 3 && TreeMasses`IsFermion[intparticles[[1]]],
              result = result * 1/2 * FlexibleSUSY`M[intfermions[[1]]] FlexibleSUSY`M[intfermions[[2]]] *
                       SARAH`D0[Sequence @@ SARAH`Mass2 /@ intparticles]];
           If[toponr == 3 && TreeMasses`IsScalar[intparticles[[1]]],
              result = result * (-1) * SARAH`D27[Sequence @@ SARAH`Mass2 /@ intparticles]];
           (*add sums over internal particles*)
           intpartwithindex = Reverse[Cases[intparticles, _[{_}]]];
           Do[result = FlexibleSUSY`SUM[
                          intpartwithindex[[i, 1, 1]],
                          If[includeGoldstones, 0,
                             TreeMasses`GetDimensionStartSkippingGoldstones[intpartwithindex[[i]]] - 1],
                          TreeMasses`GetDimension[intpartwithindex[[i]]] - 1,
                          result],
              {i, Length[intpartwithindex]}];
           result
          ];

(*combines generation of diagrams and calculation of their contributions*)
CompleteBoxResult[part1_, part2_, part3_, part4_, includeGoldstones_] :=
    Plus @@ (BoxResult[#, includeGoldstones] &) /@
       ExcludeDiagrams[GenerateDiagramsBox[part1, part2, part3, part4],
                       If[includeGoldstones, TreeMasses`IsVector,
                          TreeMasses`IsVector[#] || TreeMasses`IsGoldstone[#] &]];

(*returns the complete box part of deltaVB*)
DeltaVBbox[includeGoldstones_:False] :=
    Module[{neutrinofields, chargedleptonfields, result},
           If[!MuonDecayWorks,
              Return[{}]];
           neutrinofields = TreeMasses`GetSMNeutralLeptons[];
           chargedleptonfields = TreeMasses`GetSMChargedLeptons[];
           If[Length[neutrinofields] != Length[chargedleptonfields],
              MuonDecayWorks = False;
              DebugPrint["Error: DeltaVBbox does not work since the numbers of ",
                         "neutrino and charged lepton fields are different"];
              Return[{}]];
           If[Length[neutrinofields] == 1,
              If[TreeMasses`GetDimension[neutrinofields[[1]]] != 3 ||
                    TreeMasses`GetDimension[chargedleptonfields[[1]]] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBbox does not work since there are ",
                            "not 3 neutrinos and 3 charged leptons"];
                 Return[{}]];
              result = {WeinbergAngle`DeltaVB[{WeinbergAngle`fsbox, {SARAH`gO1, SARAH`gO2,
                                                                     SARAH`gO3, SARAH`gO4}},
                                              CompleteBoxResult[chargedleptonfields[[1]],
                                                                SARAH`bar[neutrinofields[[1]]],
                                                                neutrinofields[[1]],
                                                                SARAH`bar[chargedleptonfields[[1]]],
                                                                includeGoldstones]]},
              If[Length[neutrinofields] != 3,
                 MuonDecayWorks = False;
                 DebugPrint["Error: DeltaVBbox does not work since there are ",
                            "neither 1 nor 3 neutrino fields"];
                 Return[{}]];
              If[Or @@ ((TreeMasses`GetDimension[#] != 1) & /@
                           {neutrinofields[[1]], neutrinofields[[2]],
                            chargedleptonfields[[1]], chargedleptonfields[[2]]}),
                 MuonDecayWorks = False;
                 DebugPrint["Error: definition of neutrino or charged lepton fields not supported"];
                 Return[{}]];
              result = {WeinbergAngle`DeltaVB[{WeinbergAngle`fsbox, {}},
                                              CompleteBoxResult[chargedleptonfields[[2]],
                                                                SARAH`bar[neutrinofields[[2]]],
                                                                neutrinofields[[1]],
                                                                SARAH`bar[chargedleptonfields[[1]]],
                                                                includeGoldstones]]}];
           result
          ];

indextype = CConversion`CreateCType[CConversion`ScalarType[CConversion`integerScalarCType]];

AddIndices[{}] := "";

AddIndices[{ind_}] := indextype <> " " <> CConversion`ToValidCSymbolString[ind];

AddIndices[{ind1_, ind2_}] :=
    indextype <> " " <> CConversion`ToValidCSymbolString[ind1] <> ", " <>
    indextype <> " " <> CConversion`ToValidCSymbolString[ind2];

AddIndices[{ind1_, ind2_, ind3_, ind4_}] :=
    indextype <> " " <> CConversion`ToValidCSymbolString[ind1] <> ", " <>
    indextype <> " " <> CConversion`ToValidCSymbolString[ind2] <> ", " <>
    indextype <> " " <> CConversion`ToValidCSymbolString[ind3] <> ", " <>
    indextype <> " " <> CConversion`ToValidCSymbolString[ind4];

CreateContributionName[WeinbergAngle`DeltaVB[{type_, {___}}, _]] :=
    "delta_vb_" <> StringReplace[CConversion`ToValidCSymbolString[type],
                                 StartOfString ~~ "fs" ~~ rest_ :> rest];

CreateContributionName[WeinbergAngle`DeltaVB[{type_, {___}, spec_}, _]] :=
    "delta_vb_" <> StringReplace[CConversion`ToValidCSymbolString[type],
                                 StartOfString ~~ "fs" ~~ rest_ :> rest] <>
    "_" <> CConversion`ToValidCSymbolString[spec];

CreateContributionName[WeinbergAngle`DeltaVB[{type_, {___}, spec1_, spec2_}, _]] :=
    "delta_vb_" <> StringReplace[CConversion`ToValidCSymbolString[type],
                                 StartOfString ~~ "fs" ~~ rest_ :> rest] <>
    "_" <> CConversion`ToValidCSymbolString[spec1] <> "_" <> CConversion`ToValidCSymbolString[spec2];

CreateContributionPrototype[deltaVBcontri_WeinbergAngle`DeltaVB] :=
    CreateContributionName[deltaVBcontri] <> "(" <> AddIndices[deltaVBcontri[[1, 2]]] <> ") const";

(*creates C++ code for given part of deltaVB*)
(*based on CreateNPointFunction from SelfEnergies.m*)
CreateDeltaVBContribution[deltaVBcontri_WeinbergAngle`DeltaVB, vertexRules_List] :=
    Module[{expr, functionName, type, prototype, decl, body},
           expr = ReleaseHold[deltaVBcontri[[2]]];
           functionName = CreateContributionPrototype[deltaVBcontri];
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           prototype = type <> " " <> functionName <> ";\n";
           decl = "\n" <> type <> " CLASSNAME::" <> functionName <> "\n{\n";
           body = Parameters`CreateLocalConstRefs[expr] <> "\n";
           body = body <> "const " <> type <> " result = " <>
                  Parameters`ExpressionToString[expr /. vertexRules /. a_[List[i__]] :> a[i]] <> ";\n";
           body = body <> "\nreturn result;\n";
           body = TextFormatting`IndentText[TextFormatting`WrapLines[body]];
           decl = decl <> body <> "}\n";
           {prototype, decl}
          ];

PrintDeltaVBContributionName[WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {}, part_}, _]] :=
    "deltaVB wave-function contribution for field " <> CConversion`ToValidCSymbolString[part];

PrintDeltaVBContributionName[WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {idx_}, part_}, _]] :=
    "deltaVB wave-function contribution for field " <> CConversion`ToValidCSymbolString[part] <>
    "[" <> CConversion`ToValidCSymbolString[idx] <> "]";

PrintDeltaVBContributionName[WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, {}, part1_, part2_}, _]] :=
    "deltaVB vertex contribution for fields " <> CConversion`ToValidCSymbolString[part1] <>
    ", " <> CConversion`ToValidCSymbolString[part2];

PrintDeltaVBContributionName[WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, {__}}, _]] :=
    "deltaVB vertex contribution";

PrintDeltaVBContributionName[WeinbergAngle`DeltaVB[{WeinbergAngle`fsbox, {___}}, _]] :=
    "deltaVB box contribution";

(*creates C++ code for needed couplings and all parts of deltaVB*)
(*based on CreateNPointFunctions from SelfEnergies.m*)
CreateDeltaVBContributions[deltaVBcontris_List, vertexRules_List] :=
    Module[{relevantVertexRules, prototypes = "", defs = "", vertexFunctionNames = {}, p, d},
           Print["Converting vertex functions ..."];
           relevantVertexRules = Cases[vertexRules, r:(Rule[a_, b_] /; !FreeQ[deltaVBcontris, a])];
           {prototypes, defs, vertexFunctionNames} =
              SelfEnergies`CreateVertexExpressions[relevantVertexRules, False];
           Print["Generating C++ code for deltaVB ..."];
           For[k = 1, k <= Length[deltaVBcontris], k++,
               DebugPrint["   ", PrintDeltaVBContributionName[deltaVBcontris[[k]]]];
               {p, d} = CreateDeltaVBContribution[deltaVBcontris[[k]], vertexFunctionNames];
               prototypes = prototypes <> p;
               defs = defs <> d];
           {prototypes, defs}
          ];

GetNeutrinoIndex[] :=
    Module[{neutrinofield, chargedleptonfield, coupl, result = ""},
           If[!MuonDecayWorks,
              Return["return 0;"]];
           neutrinofield = TreeMasses`GetSMNeutralLeptons[];
           If[Length[neutrinofield] == 3,
              Return["return 0;"]];
           neutrinofield = neutrinofield[[1]];
           chargedleptonfield = TreeMasses`GetSMChargedLeptons[][[1]];
           coupl = SARAH`Cp[SARAH`bar[neutrinofield][{SARAH`gO2}], chargedleptonfield[{Global`FeIdx}],
                            GetWPlusBoson[]][SARAH`PL];
           (*follow vertex conventions:*)
           coupl = Vertices`SortCp[coupl];
           (*omit a possible minus sign:*)   
           If[MatchQ[coupl, Times[-1, _]], coupl = -coupl];
           coupl = SelfEnergies`CreateCouplingSymbol[coupl];
           For[k = 0, k <= 2, k++,
               result = result <> "const auto Cp" <> ToString[k] <> " = Abs(";
               result = result <> CConversion`RValueToCFormString[coupl /. SARAH`gO2 -> k] <> ");\n"];
           For[k = 0, k <= 2, k++,
               result = result <> "\nif (Cp" <> ToString[k] <> " >= std::max(";
               result = result <> Utils`StringJoinWithSeparator[
                                      ("Cp" <> ToString[#])& /@ Complement[{0,1,2},{k}], ","];
               result = result <> "))\n   return " <> ToString[k] <> ";"];
           result = result <> "\n\nthrow NonPerturbativeSinThetaW();";
           result
          ];

CreateContributionCall[deltaVBcontri_ /; MatchQ[deltaVBcontri,
       WeinbergAngle`DeltaVB[{_, {}, ___}, _]]] :=
    CreateContributionName[deltaVBcontri] <> "()";

CreateContributionCall[deltaVBcontri_ /; MatchQ[deltaVBcontri,
       WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {SARAH`gO1},
                              TreeMasses`GetSMChargedLeptons[][[1]]}, _]]] :=
    CreateContributionName[deltaVBcontri] <> "(0) + " <>
    CreateContributionName[deltaVBcontri] <> "(1)";

CreateContributionCall[deltaVBcontri_ /; MatchQ[deltaVBcontri,
       WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, {SARAH`gO1},
                              TreeMasses`GetSMNeutralLeptons[][[1]]}, _]]] :=
    CreateContributionName[deltaVBcontri] <> "(FveIdx) + " <>
    CreateContributionName[deltaVBcontri] <> "(FvmIdx)";

CreateContributionCall[deltaVBcontri_ /; MatchQ[deltaVBcontri,
       WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, {SARAH`gO1, SARAH`gO2}}, _]]] :=
    CreateContributionName[deltaVBcontri] <> "(0, FveIdx) + " <>
    CreateContributionName[deltaVBcontri] <> "(1, FvmIdx)";

CreateContributionCall[deltaVBcontri_ /; MatchQ[deltaVBcontri,
       WeinbergAngle`DeltaVB[{WeinbergAngle`fsbox, {SARAH`gO1, SARAH`gO2, SARAH`gO3, SARAH`gO4}}, _]]] :=
    CreateContributionName[deltaVBcontri] <> "(1, FvmIdx, FveIdx, 0)";

CreateContributionCall[0] := "0."; (*needed in case of an error*)

(*creates C++ code for calling the different functions contributing to deltaVB*)
CreateDeltaVBCalculation[deltaVBcontris_List] :=
    Module[{type, result = "", boxcontri, vertexcontris, wavecontris},
           If[!(TreeMasses`FindMixingMatrixSymbolFor[TreeMasses`GetSMNeutralLeptons[][[1]]] === Null),
              Print["Warning: neutrino mixing is not considered in muon decay"]];
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           boxcontri = Cases[deltaVBcontris, WeinbergAngle`DeltaVB[{WeinbergAngle`fsbox, __}, _]];
           If[boxcontri === {},
              boxcontri = 0,
              boxcontri = boxcontri[[1]]];
           vertexcontris = Cases[deltaVBcontris, WeinbergAngle`DeltaVB[{WeinbergAngle`fsvertex, __}, _]];
           If[vertexcontris === {},
              vertexcontris = {0}];
           wavecontris = Cases[deltaVBcontris, WeinbergAngle`DeltaVB[{WeinbergAngle`fswave, __}, _]];
           If[wavecontris === {},
              wavecontris = {0}];
           result = result <> "const " <> type <> " a1 = ";
           result = result <> CreateContributionCall[boxcontri] <> ";\n";
           result = result <> "const " <> type <> " deltaV =\n   ";
           For[k = 1, k <= Length[vertexcontris], k++,
               If[k > 1, result = result <> " + "];
               result = result <> CreateContributionCall[vertexcontris[[k]]]];
           result = result <> ";\n";
           result = result <> "const " <> type <> " deltaZ =\n   ";
           For[k = 1, k <= Length[wavecontris], k++,
               If[k > 1, result = result <> " + "];
               result = result <> CreateContributionCall[wavecontris[[k]]]];
           result <> ";"
          ];

End[];

EndPackage[];
