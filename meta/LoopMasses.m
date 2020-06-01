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

BeginPackage["LoopMasses`", {"SARAH`", "TextFormatting`",
                             "CConversion`", "TreeMasses`",
                             "SelfEnergies`", "TwoLoopQCD`", "ThreeLoopQCD`",
                             "Parameters`", "Utils`"}];

CreateOneLoopPoleMassFunctions::usage="";
CreateOneLoopPoleMassPrototypes::usage="";
CallAllPoleMassFunctions::usage="";
CreateRunningDRbarMassPrototypes::usage="";
CreateRunningDRbarMassFunctions::usage="";
CallCalculateDRbarMass::usage="";
CallPoleMassFunction::usage="";
CallThreadedPoleMassFunction::usage="";
CreateLoopMassFunctionName::usage="";

GetLoopCorrectedParticles::usage="Returns list of all particles that
get loop corrected masses.  These are all particles, except for
ghosts.";

CreateLSPFunctions::usage="";

Begin["`Private`"];

GetLoopCorrectedParticles[states_] :=
    Module[{particles},
           particles = GetParticles[states];
           Select[particles, (!IsGhost[#] && !IsGoldstone[#])&]
          ];

FillTadpoleMatrix[{}, _] := "";

FillTadpoleMatrix[tadpoles_List, matrixName_:"tadpoles"] :=
    Module[{result, dim, dimStr, particle, i, particleIndex, vev,
            tadpoleMatrixType},
           particle = tadpoles[[1,1]];
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           tadpoleMatrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType, dim, dim]];
           If[dim > 1,
              result = tadpoleMatrixType <> " " <> matrixName <> "(" <> dimStr <> "," <> dimStr <> ");\n";
              For[i = 1, i <= Length[tadpoles], i++,
                  particleIndex = ToString[tadpoles[[i,2]]];
                  vev           = ToValidCSymbolString[tadpoles[[i,3]]];
                  result = result <>
                           matrixName <> "(" <> particleIndex <> "," <> particleIndex <> ") = Re(" <>
                           CreateTadpoleFunctionName[particle] <> "(" <> particleIndex <> ")) / " <>
                           vev <> ";\n";
                 ];
              ,
              particleIndex = ToString[tadpoles[[1,2]]];
              vev           = ToValidCSymbolString[tadpoles[[1,3]]];
              result = "const double " <> matrixName <> " = " <>
                       CreateTadpoleFunctionName[particle] <> "(" <> particleIndex <> ");\n";
             ];
           Return[result];
          ];

FillMt2LStruct[] := "\
double mst_1, mst_2, theta_t;
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <> ";

mssm_twoloop_mt::Parameters pars;
pars.g3 = " <> CConversion`RValueToCFormString[SARAH`strongCoupling /. Parameters`ApplyGUTNormalization[]] <> ";
pars.mt = " <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[]]] <> ";
pars.mg = " <> CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]] <> ";
pars.mst1 = mst_1;
pars.mst2 = mst_2;
pars.msusy = " <> CConversion`RValueToCFormString[Sqrt[Sqrt[Abs[SARAH`SoftSquark[2,2]] Abs[SARAH`SoftDown[2,2]]]]] <> ";
pars.xt = Sin(2*theta_t) * (Sqr(mst_1) - Sqr(mst_2)) / (2. * pars.mt);
pars.Q = get_scale();";


AddMtPoleQCDCorrections[1, expr_] :=
    If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
"\
double qcd_1l = 0.;

{
" <> IndentText[FillMt2LStruct[]] <> "

   qcd_1l = - mssm_twoloop_mt::dMt_over_mt_1loop(pars);
}
"
       ,
"\
double qcd_1l = 0.;

{
   const double currentScale = get_scale();
   qcd_1l = " <> CConversion`RValueToCFormString[expr] <> ";
}
"
    ];

AddMtPoleQCDCorrections[2, expr_] :=
    If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
"\
double qcd_2l = 0.;

if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
" <> IndentText[FillMt2LStruct[]] <> "

   qcd_2l = - mssm_twoloop_mt::dMt_over_mt_2loop(pars);
}
"
       ,
 "\
double qcd_2l = 0.;

if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
   const double currentScale = get_scale();
   qcd_2l = " <> CConversion`RValueToCFormString[expr] <> ";
}
"
    ];

AddMtPoleQCDCorrections[3, expr_] := "\
double qcd_3l = 0.;

if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
   const double currentScale = get_scale();
   qcd_3l = " <> CConversion`RValueToCFormString[expr] <> ";
}
";

AddMtPoleQCDCorrections[4, expr_] := "\
double qcd_4l = 0.;

if (pole_mass_loop_order > 3 && TOP_POLE_QCD_CORRECTION > 2) {
   const double currentScale = get_scale();
   qcd_4l = " <> CConversion`RValueToCFormString[expr] <> ";
}
";

AddMbRun2LSQCDCorrections[] :=
    If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
"
if (get_thresholds() > 1 && threshold_corrections.mb > 1) {
   double mst_1, mst_2, theta_t;
   " <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <> ";
   double msb_1, msb_2, theta_b;
   " <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <> ";

   mssm_twoloop_mb::Parameters pars;
   pars.g3 = " <> CConversion`RValueToCFormString[SARAH`strongCoupling /. Parameters`ApplyGUTNormalization[]] <> ";
   pars.mt = " <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMTopQuarkMultiplet[]]] <> ";
   pars.mb = " <> CConversion`RValueToCFormString[TreeMasses`GetThirdGenerationMass[TreeMasses`GetSMBottomQuarkMultiplet[]]] <> ";
   pars.mg = " <> CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]] <> ";
   pars.mst1 = mst_1;
   pars.mst2 = mst_2;
   pars.msb1 = msb_1;
   pars.msb2 = msb_2;
   pars.msusy = " <> CConversion`RValueToCFormString[
       Power[
           AbsSqrt[SARAH`SoftSquark[0,0]] AbsSqrt[SARAH`SoftUp[0,0]] AbsSqrt[SARAH`SoftDown[0,0]]
           AbsSqrt[SARAH`SoftSquark[1,1]] AbsSqrt[SARAH`SoftUp[1,1]] AbsSqrt[SARAH`SoftDown[1,1]]
           , 1/6]
   ] <> ";
   pars.xt = Sin(2*theta_t) * (Sqr(mst_1) - Sqr(mst_2)) / (2. * pars.mt);
   pars.xb = Sin(2*theta_b) * (Sqr(msb_1) - Sqr(msb_2)) / (2. * pars.mb);
   pars.Q = get_scale();

   const double alpha_s = Sqr(pars.g3) / (4.*Pi);

   qcd_2l = mssm_twoloop_mb::delta_mb_2loop(pars)
            - alpha_s/(3.*Pi) * delta_mb_1loop;
}
"
       ,
       ""
      ];

(* bring logs into canonical form *)
SimpQCD[expr_] :=
    Collect[expr //. {
        Log[Global`currentScale^2/m_] :> -Log[m/Global`currentScale^2],
        Log[Global`currentScale^2/m_^2] :> -Log[m^2/Global`currentScale^2]
    }, {SARAH`strongCoupling}, Simplify];

Get[FileNameJoin[{"meta", "SM", "mf_3loop_qcd.m"}]];

Get3LQCDMtRelations[particle_, scale_] :=
    SimpQCD[
        k^3 { Coefficient[MfOvermf, k, 3], Coefficient[mfOverMf, k, 3] } //. {
            L -> Log[mf^2/scale^2],
            mf -> FlexibleSUSY`M[particle],
            k -> 1/(4Pi)^2,
            NL -> 5,
            NH -> 1,
            g3 -> SARAH`strongCoupling
        }
    ];

Get3LQCDMtOvermt[particle_, scale_] := Get3LQCDMtRelations[particle, scale][[1]];
Get3LQCDmtOverMt[particle_, scale_] := Get3LQCDMtRelations[particle, scale][[2]];

Get[FileNameJoin[{"meta", "SM", "mt_4loop_qcd.m"}]];

Get4LQCDMtRelations[particle_, scale_] :=
    SimpQCD[
        k^4 { Coefficient[MtOvermt, k, 4], Coefficient[mtOverMt, k, 4] } //. {
            L -> Log[mt^2/scale^2],
            mt -> FlexibleSUSY`M[particle],
            k -> 1/(4Pi)^2,
            g3 -> SARAH`strongCoupling
        }
    ];

Get4LQCDMtOvermt[particle_, scale_] := Get4LQCDMtRelations[particle, scale][[1]];
Get4LQCDmtOverMt[particle_, scale_] := Get4LQCDMtRelations[particle, scale][[2]];

Do1DimScalar[particle_, particleName_String, massName_String, massMatrixName_String,
             selfEnergyFunction_String, momentum_String, tadpole_String:""] :=
    "const double p = " <> momentum <> ";\n" <>
    "double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
    If[FlexibleSUSY`UseHiggs2LoopSM === True && particle === SARAH`HiggsBoson,
       "if (pole_mass_loop_order > 1)\n" <>
       IndentText["self_energy += self_energy_" <> particleName <> "_2loop(p);\n"], ""] <>
    If[FlexibleSUSY`UseHiggs3LoopSM === True && particle === SARAH`HiggsBoson,
       "if (pole_mass_loop_order > 2)\n" <>
       IndentText["self_energy += self_energy_" <> particleName <> "_3loop();\n"], ""] <>
    If[FlexibleSUSY`UseHiggs3LoopSplit === True && particle === SARAH`HiggsBoson,
       "if (pole_mass_loop_order > 2)\n" <>
       IndentText["self_energy += self_energy_" <> particleName <> "_3loop();\n"], ""] <>
    If[FlexibleSUSY`UseHiggs4LoopSM === True && particle === SARAH`HiggsBoson,
       "if (pole_mass_loop_order > 3)\n" <>
       IndentText["self_energy += self_energy_" <> particleName <> "_4loop();\n"], ""] <>
    "const double mass_sqr = " <> massMatrixName <> " - self_energy" <>
    If[tadpole == "", "", " + " <> tadpole] <> ";\n\n" <>
    "PHYSICAL(" <> massName <> ") = SignedAbsSqrt(mass_sqr);\n";

Do1DimFermion[particle_, massMatrixName_String, selfEnergyFunctionS_String,
              selfEnergyFunctionPL_String, selfEnergyFunctionPR_String, momentum_String, type_] :=
    "const double p = " <> momentum <> ";\n" <>
    "const " <> CreateCType[type] <> " self_energy_1  = " <> CastIfReal[selfEnergyFunctionS  <> "(p)",type] <> ";\n" <>
    "const " <> CreateCType[type] <> " self_energy_PL = " <> CastIfReal[selfEnergyFunctionPL <> "(p)",type] <> ";\n" <>
    "const " <> CreateCType[type] <> " self_energy_PR = " <> CastIfReal[selfEnergyFunctionPR <> "(p)",type] <> ";\n" <>
    "const auto M_loop = " <> massMatrixName <>
    " - self_energy_1 - " <> massMatrixName <> " * (self_energy_PL + self_energy_PR);\n" <>
    "PHYSICAL(" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> ") = " <>
    "calculate_singlet_mass(M_loop);\n";

Do1DimFermion[particle_ /; particle === TreeMasses`GetSMTopQuarkMultiplet[], massMatrixName_String,
              _String, _String, _String, momentum_String, type_] :=
    Module[{massName,
            topSelfEnergyFunctionS, topSelfEnergyFunctionPL, topSelfEnergyFunctionPR,
            qcdOneLoop, qcdTwoLoop, qcdThreeLoop, qcdFourLoop
           },
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           topSelfEnergyFunctionS  = SelfEnergies`CreateHeavySelfEnergyFunctionName[particle[1], 1];
           topSelfEnergyFunctionPL = SelfEnergies`CreateHeavySelfEnergyFunctionName[particle[PL], 1];
           topSelfEnergyFunctionPR = SelfEnergies`CreateHeavySelfEnergyFunctionName[particle[PR], 1];
           qcdOneLoop = N[SimpQCD[-TwoLoopQCD`GetDeltaMPoleOverMRunningQCDOneLoop[particle, Global`currentScale, FlexibleSUSY`FSRenormalizationScheme]]];
           qcdTwoLoop = N[SimpQCD[-TwoLoopQCD`GetDeltaMPoleOverMRunningQCDTwoLoop[particle, Global`currentScale, FlexibleSUSY`FSRenormalizationScheme]]];
           qcdThreeLoop = If[FlexibleSUSY`FSRenormalizationScheme === FlexibleSUSY`MSbar,
                             N @ SimpQCD @ Get3LQCDMtOvermt[particle, Global`currentScale],
                             0];
           qcdFourLoop  = If[FlexibleSUSY`FSRenormalizationScheme === FlexibleSUSY`MSbar,
                             N @ SimpQCD @ Get4LQCDMtOvermt[particle, Global`currentScale],
                             0];
           (* adding split-MSSM contributions if enabled *)
           If[FlexibleSUSY`UseHiggs3LoopSplit === True,
              qcdTwoLoop = N[SimpQCD[qcdTwoLoop - GetMtPoleOverMtMSbarSplitMSSM2L[particle, Global`currentScale]]];
             ];
           AddMtPoleQCDCorrections[1, qcdOneLoop] <> "\n" <>
           AddMtPoleQCDCorrections[2, qcdTwoLoop] <> "\n" <>
           AddMtPoleQCDCorrections[3, qcdThreeLoop] <> "\n" <>
           AddMtPoleQCDCorrections[4, qcdFourLoop] <> "
const double p = " <> momentum <> ";
const " <> CreateCType[type] <> " self_energy_1  = " <> CastIfReal[topSelfEnergyFunctionS  <> "(p)",type] <> ";
const " <> CreateCType[type] <> " self_energy_PL = " <> CastIfReal[topSelfEnergyFunctionPL <> "(p)",type] <> ";
const " <> CreateCType[type] <> " self_energy_PR = " <> CastIfReal[topSelfEnergyFunctionPR <> "(p)",type] <> ";
const auto M_loop = " <> massMatrixName <> "\
 - self_energy_1 - " <> massMatrixName <> " * (self_energy_PL + self_energy_PR)\
 - " <> massMatrixName <> " * (qcd_1l + qcd_2l + qcd_3l + qcd_4l);\n
PHYSICAL(" <> massName <> ") = calculate_singlet_mass(M_loop);\n"
          ];

Do1DimVector[particleName_String, massName_String, massMatrixName_String,
             selfEnergyFunction_String, momentum_String] :=
    "const double p = " <> momentum <> ";\n" <>
    "const double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
    "const double mass_sqr = " <> massMatrixName <> " - self_energy;\n\n" <>
    "if (mass_sqr < 0.) {\n" <>
    IndentText[TreeMasses`FlagPoleTachyon[particleName]] <> "\n}\n\n" <>
    "PHYSICAL(" <> massName <> ") = AbsSqrt(mass_sqr);\n";


(* ********** fast diagonalization routines ********** *)

CastIfReal[str_String, type_[CConversion`realScalarCType, ___]] :=
    "Re(" <> str <> ")";

CastIfReal[str_String, type_[___]] := str;

DoFastDiagonalization[particle_Symbol /; IsScalar[particle], tadpoles_List] :=
    Module[{result, dim, dimStr, massName, massNameReordered, particleName,
            mixingMatrix, selfEnergyFunction, reorderMasses,
            tadpoleMatrix, U, V, massMatrixStr, selfEnergyIsSymmetric,
            selfEnergyMatrixType, selfEnergyMatrixCType},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           massNameReordered = massName <> "_reordered";
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           massMatrixStr = "get_mass_matrix_" <> ToValidCSymbolString[particle];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle, 1];
           selfEnergyMatrixType = TreeMasses`GetMassMatrixType[particle];
           selfEnergyMatrixCType = CreateCType[selfEnergyMatrixType];
           tadpoleMatrix = FillTadpoleMatrix[tadpoles, "tadpoles"];
           reorderMasses = CreateCType[CConversion`ArrayType[realScalarCType, dim]] <> " " <>
                       massNameReordered <> "(" <> massName <> ");\n" <>
                       "reorder_vector(" <> massNameReordered <> ", " <>
                       "get_mass_matrix_" <> particleName <> "());\n";
           If[dim > 1,
              selfEnergyIsSymmetric = Length[Flatten[{mixingMatrix}]] === 1;
              result = reorderMasses <> "\n" <>
                       tadpoleMatrix <>
                       selfEnergyMatrixCType <> " self_energy;\n" <>
                       "for (int i1 = 0; i1 < " <> dimStr <>"; ++i1) {\n" <>
                       IndentText["for (int i2 = " <> If[selfEnergyIsSymmetric,"i1","0"] <>
                                  "; i2 < " <> dimStr <>"; ++i2) {\n" <>
                                  IndentText["const double p = AbsSqrt(" <> massNameReordered <> "(i1) * " <>
                                             massNameReordered <> "(i2));\n" <>
                                             "self_energy(i1,i2) = " <> CastIfReal[
                                             selfEnergyFunction <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n"] <>
                                  "}\n"
                                 ] <>
                       "}\n" <>
                       If[selfEnergyIsSymmetric, "Symmetrize(self_energy);\n", ""] <>
                       "const " <> selfEnergyMatrixCType <>" M_loop(" <> massMatrixStr <>
                       "() - self_energy" <>
                       If[tadpoleMatrix == "", "", " + tadpoles"] <> ");\n";
              If[Head[mixingMatrix] === List,
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          TreeMasses`CallSVDFunction[
                              particle, "M_loop", "PHYSICAL(" <> massName <> ")",
                              "PHYSICAL(" <> U <> ")", "PHYSICAL(" <> V <> ")"];
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 result = result <>
                          TreeMasses`CallDiagonalizeHermitianFunction[
                              particle, "M_loop", "PHYSICAL(" <> massName <> ")",
                              "PHYSICAL(" <> U <> ")"];
                ];
              result = result <>
                       "\n" <>
                       "PHYSICAL(" <> massName <> ") = SignedAbsSqrt(PHYSICAL(" <>
                       massName <> "));\n";
              ,
              result = "const " <> selfEnergyMatrixCType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       Do1DimScalar[particle, particleName, massName, "M_tree", selfEnergyFunction, massName,
                                    If[tadpoleMatrix == "", "", "tadpoles"]];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol /; IsFermion[particle], _] :=
    Module[{result, dim, dimStr, massName, mixingMatrix, U, V,
            massNameReordered, reorderMasses,
            selfEnergyFunctionS, selfEnergyFunctionPL, selfEnergyFunctionPR,
            massMatrixStr, selfEnergyMatrixType, selfEnergyMatrixCType, particleName},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           massNameReordered = massName <> "_reordered";
           massMatrixStr = "get_mass_matrix_" <> particleName;
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           selfEnergyMatrixType = TreeMasses`GetMassMatrixType[particle];
           selfEnergyMatrixCType = CreateCType[selfEnergyMatrixType];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1], 1];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL], 1];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR], 1];
           reorderMasses = CreateCType[CConversion`ArrayType[realScalarCType, dim]] <> " " <>
                       massNameReordered <> "(" <> massName <> ");\n" <>
                       "reorder_vector(" <> massNameReordered <> ", " <>
                       "get_mass_matrix_" <> particleName <> "());\n";
           If[dim > 1,
              result = reorderMasses <> "\n" <>
                       selfEnergyMatrixCType <> " self_energy_1;\n" <>
                       selfEnergyMatrixCType <> " self_energy_PL;\n" <>
                       selfEnergyMatrixCType <> " self_energy_PR;\n" <>
                       "for (int i1 = 0; i1 < " <> dimStr <>"; ++i1) {\n" <>
                       IndentText["for (int i2 = 0; i2 < " <> dimStr <>"; ++i2) {\n" <>
                                  IndentText["const double p = AbsSqrt(" <> massNameReordered <> "(i1) * " <>
                                             massNameReordered <> "(i2));\n" <>
                                             "self_energy_1(i1,i2)  = " <> CastIfReal[
                                             selfEnergyFunctionS <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n" <>
                                             "self_energy_PL(i1,i2) = " <> CastIfReal[
                                             selfEnergyFunctionPL <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n" <>
                                             "self_energy_PR(i1,i2) = " <> CastIfReal[
                                             selfEnergyFunctionPR <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n"
                                            ] <>
                                  "}\n"
                                 ] <>
                       "}\n" <>
                       "const " <> selfEnergyMatrixCType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       "const " <> selfEnergyMatrixCType <> " delta_M(- self_energy_PR * M_tree " <>
                       "- M_tree * self_energy_PL - self_energy_1);\n";
              If[IsMajoranaFermion[particle],
                 result = result <>
                          "const " <> selfEnergyMatrixCType <>
                          " M_loop(M_tree + 0.5 * (delta_M + delta_M.transpose()));\n";
                 ,
                 result = result <>
                          "const " <> selfEnergyMatrixCType <>
                          " M_loop(M_tree + delta_M);\n";
                ];
              If[Head[mixingMatrix] === List,
                 (* two mixing matrixs => SVD *)
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          TreeMasses`CallSVDFunction[
                              particle, "M_loop", "PHYSICAL(" <> massName <> ")",
                              "PHYSICAL(" <> U <> ")", "PHYSICAL(" <> V <> ")"];
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 result = result <>
                          TreeMasses`CallDiagonalizeSymmetricFunction[
                              particle, "M_loop", "PHYSICAL(" <> massName <> ")",
                              "PHYSICAL(" <> U <> ")"];
                ];
              ,
              (* for a dimension 1 fermion it plays not role if it's a
                 Majorana ferimion or not *)

              (* Note: For a 1-dimensional fermion multiplet SARAH
                 provides the self-energies in mass eigenstates, i.e.
                 the fermions at the external legs are multiplied by
                 their phase (= mixing matrix).  Therfore, M_tree must
                 be set to the (positive) tree-level mass.  M_tree
                 must not be set to the gauge eigenstate mass
                 parameter! *)

              result = "const " <> selfEnergyMatrixCType <> " M_tree(" <> massName <> ");\n" <>
                       Do1DimFermion[particle, "M_tree", selfEnergyFunctionS,
                                     selfEnergyFunctionPL, selfEnergyFunctionPR,
                                     massName, CConversion`GetScalarElementType[selfEnergyMatrixType]];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol /; IsVector[particle], _] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix,
            selfEnergyFunction, selfEnergyMatrixType, selfEnergyMatrixCType},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           mixingMatrix = ToValidCSymbolString[FindMixingMatrixSymbolFor[particle]];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle, 1];
           selfEnergyMatrixType = TreeMasses`GetMassMatrixType[particle];
           selfEnergyMatrixCType = CreateCType[selfEnergyMatrixType];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           If[dim > 1,
              result = "WARNING(\"diagonalization of " <> ToString[particle] <> " not implemented\");\n";
              ,
              result = "const " <> selfEnergyMatrixCType <> " M_tree(Sqr(" <> massName <> "));\n" <>
                       Do1DimVector[particleName, massName, "M_tree", selfEnergyFunction, massName];
             ];
           Return[result];
          ];

DoFastDiagonalization[particle_Symbol, _] :=
    "ERROR(\"fast diagonalization of " <> ToString[particle] <> " not implemented\");\n";

(* ********** medium diagonalization routines ********** *)

DoMediumDiagonalization[particle_Symbol /; IsScalar[particle], inputMomentum_, tadpole_List, calcEffPot_:True] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix, selfEnergyFunction,
            momentum = inputMomentum, U, V, Utemp, Vtemp, tadpoleMatrix, diagSnippet,
            massMatrixStr,
            selfEnergyMatrixType, selfEnergyMatrixCType, eigenArrayType,
            addHigherLoopHiggsContributions = "", calcHigherLoopHiggsContributions = "", n, l, k},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[inputMomentum == "", momentum = massName];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           mixingMatrixType = TreeMasses`GetMassMatrixType[particle];
           selfEnergyMatrixType = mixingMatrixType;
           mixingMatrixType = CreateCType[mixingMatrixType];
           selfEnergyMatrixCType = mixingMatrixType;
           eigenArrayType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]];
           (* create diagonalisation code snippet *)
           If[Head[mixingMatrix] === List,
              U = ToValidCSymbolString[mixingMatrix[[1]]];
              V = ToValidCSymbolString[mixingMatrix[[2]]];
              Utemp = "mix_" <> U;
              Vtemp = "mix_" <> V;
              diagSnippet = mixingMatrixType <> " " <> Utemp <> ", " <> Vtemp <> ";\n" <>
                            TreeMasses`CallSVDFunction[
                                particle, "M_loop", "eigen_values",
                                Utemp, Vtemp] <> "\n" <>
                            "PHYSICAL(" <> massName <> "(es)) = SignedAbsSqrt(eigen_values(es));\n" <>
                            "if (es == " <> ToString[GetDimensionStartSkippingGoldstones[particle]-1] <> ") {\n" <>
                            IndentText["PHYSICAL(" <> U <> ") = " <> Utemp <> ";\n" <>
                                       "PHYSICAL(" <> V <> ") = " <> Vtemp <> ";\n"] <>
                            "}\n";
              ,
              U = ToValidCSymbolString[mixingMatrix];
              Utemp = "mix_" <> U;
              diagSnippet = mixingMatrixType <> " " <> Utemp <> ";\n" <>
                            TreeMasses`CallDiagonalizeHermitianFunction[
                                particle, "M_loop", "eigen_values",
                                Utemp] <> "\n" <>
                            "PHYSICAL(" <> massName <> "(es)) = SignedAbsSqrt(eigen_values(es));\n";
              If[mixingMatrix =!= Null,
                 diagSnippet = diagSnippet <>
                               "if (es == " <> ToString[GetDimensionStartSkippingGoldstones[particle]-1] <> ")\n" <>
                               IndentText["PHYSICAL(" <> U <> ") = " <> Utemp <> ";\n"];
                ];
             ];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle, 1];
           tadpoleMatrix = FillTadpoleMatrix[tadpole, "tadpoles"];
           massMatrixStr = "get_mass_matrix_" <> ToValidCSymbolString[particle];
           (* fill self-energy and do diagonalisation *)
           If[dim > 1,
              If[(SARAH`UseHiggs2LoopMSSM === True ||
                  FlexibleSUSY`UseHiggs2LoopNMSSM === True ||
                  FlexibleSUSY`UseHiggs3LoopMSSM === True) &&
                 MemberQ[{SARAH`HiggsBoson, SARAH`PseudoScalar}, particle],
                 addHigherLoopHiggsContributions = "self_energy += self_energy_2l;\n";
                 If[calcEffPot,
                    calcHigherLoopHiggsContributions = CalcEffPot2L[particle];
                   ];
                ];
              If[FlexibleSUSY`UseHiggs3LoopMSSM === True && MemberQ[{SARAH`HiggsBoson}, particle],
                 addHigherLoopHiggsContributions = addHigherLoopHiggsContributions <> "self_energy += self_energy_3l;\n";
                 If[calcEffPot,
                    calcHigherLoopHiggsContributions = calcHigherLoopHiggsContributions <> CalcEffPot3L[particle];
                   ];
                ];
              result = tadpoleMatrix <>
                       "const " <> selfEnergyMatrixCType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       calcHigherLoopHiggsContributions <> "\n" <>
                       "for (int es = 0; es < " <> dimStr <> "; ++es) {\n" <>
                       IndentText["const double p = Abs(" <> momentum <> "(es));\n" <>
                                  selfEnergyMatrixCType <> " self_energy = " <> CastIfReal[selfEnergyFunction <> "(p)", selfEnergyMatrixType] <> ";\n" <>
                                  addHigherLoopHiggsContributions <>
                                  "const " <> selfEnergyMatrixCType <> " M_loop(M_tree - self_energy" <>
                                  If[tadpoleMatrix == "", "", " + tadpoles"] <> ");\n" <>
                                  eigenArrayType <> " eigen_values;\n" <>
                                  diagSnippet
                                 ] <>
                       "}\n";
              ,
              result = tadpoleMatrix <>
                       "const " <> selfEnergyMatrixCType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       Do1DimScalar[particle, particleName, massName, "M_tree", selfEnergyFunction, momentum,
                                    If[tadpoleMatrix == "", "", "tadpoles"]];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol /; IsFermion[particle], inputMomentum_, __] :=
    Module[{result, dim, dimStr, massName, mixingMatrix, U, V,
            selfEnergyFunctionS, selfEnergyFunctionPL, selfEnergyFunctionPR,
            momentum = inputMomentum, massMatrixStr,
            selfEnergyMatrixType, selfEnergyMatrixCType,
            eigenArrayType, mixingMatrixType, particleName,
            topSelfEnergyFunctionS, topSelfEnergyFunctionPL, topSelfEnergyFunctionPR,
            topTwoLoop = False, thirdGenMass, qcdCorrections = "",
            qcdOneLoop, qcdTwoLoop, qcdThreeLoop, qcdFourLoop
           },
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[inputMomentum == "", momentum = massName];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           mixingMatrix = FindMixingMatrixSymbolFor[particle];
           massMatrixStr = "get_mass_matrix_" <> particleName;
           selfEnergyMatrixType = TreeMasses`GetMassMatrixType[particle];
           selfEnergyMatrixCType = CreateCType[selfEnergyMatrixType];
           eigenArrayType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]];
           topTwoLoop = particle === TreeMasses`GetSMTopQuarkMultiplet[];
           If[topTwoLoop,
              topSelfEnergyFunctionS  = SelfEnergies`CreateHeavySelfEnergyFunctionName[particle[1], 1];
              topSelfEnergyFunctionPL = SelfEnergies`CreateHeavySelfEnergyFunctionName[particle[PL], 1];
              topSelfEnergyFunctionPR = SelfEnergies`CreateHeavySelfEnergyFunctionName[particle[PR], 1];
              thirdGenMass = TreeMasses`GetThirdGenerationMass[particle];
              qcdOneLoop = N[SimpQCD[-TwoLoopQCD`GetDeltaMPoleOverMRunningQCDOneLoop[particle, Global`currentScale, FlexibleSUSY`FSRenormalizationScheme]]];
              qcdTwoLoop = N[SimpQCD[-TwoLoopQCD`GetDeltaMPoleOverMRunningQCDTwoLoop[particle, Global`currentScale, FlexibleSUSY`FSRenormalizationScheme]]];
              (* adding split-MSSM contributions if enabled *)
              If[FlexibleSUSY`UseHiggs3LoopSplit === True,
                 qcdTwoLoop = N[SimpQCD[qcdTwoLoop - GetMtPoleOverMtMSbarSplitMSSM2L[particle, Global`currentScale]]];
                ];
              qcdThreeLoop = If[FlexibleSUSY`FSRenormalizationScheme === FlexibleSUSY`MSbar,
                                (* contribution to self-energy => neg. sign *)
                                N @ SimpQCD[-Get3LQCDMtOvermt[particle, Global`currentScale]],
                                0];
              qcdFourLoop  = If[FlexibleSUSY`FSRenormalizationScheme === FlexibleSUSY`MSbar,
                                (* contribution to self-energy => neg. sign *)
                                N @ SimpQCD[-Get4LQCDMtOvermt[particle, Global`currentScale]],
                                0];
              qcdCorrections = AddMtPoleQCDCorrections[1, qcdOneLoop /. FlexibleSUSY`M[particle] -> thirdGenMass] <> "\n" <>
                               AddMtPoleQCDCorrections[2, qcdTwoLoop /. FlexibleSUSY`M[particle] -> thirdGenMass] <> "\n" <>
                               AddMtPoleQCDCorrections[3, qcdThreeLoop /. FlexibleSUSY`M[particle] -> thirdGenMass] <> "\n" <>
                               AddMtPoleQCDCorrections[4, qcdFourLoop /. FlexibleSUSY`M[particle] -> thirdGenMass] <> "\n";
             ];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1], 1];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL], 1];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR], 1];
           If[dim > 1,
              result = qcdCorrections <>
                       "const " <> selfEnergyMatrixCType <> " M_tree(" <> massMatrixStr <> "());\n" <>
                       "for (int es = 0; es < " <> dimStr <> "; ++es) {\n" <>
                       IndentText["const double p = Abs(" <> momentum <> "(es));\n" <>
                                  If[topTwoLoop,
                                     selfEnergyMatrixCType <> " self_energy_1;\n" <>
                                     selfEnergyMatrixCType <> " self_energy_PL;\n" <>
                                     selfEnergyMatrixCType <> " self_energy_PR;\n" <>
                                     "for (int i1 = 0; i1 < " <> dimStr <>"; ++i1) {\n" <>
                                     IndentText["for (int i2 = 0; i2 < " <> dimStr <>"; ++i2) {\n" <>
                                                IndentText[
                                                    "if (i1 == 2 && i2 == 2) {\n" <>
                                                    IndentText["self_energy_1(i1,i2)  = " <> CastIfReal[topSelfEnergyFunctionS <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n" <>
                                                               "self_energy_PL(i1,i2) = " <> CastIfReal[topSelfEnergyFunctionPL <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n" <>
                                                               "self_energy_PR(i1,i2) = " <> CastIfReal[topSelfEnergyFunctionPR <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n"] <>
                                                    "} else {\n" <>
                                                    IndentText["self_energy_1(i1,i2)  = " <> CastIfReal[selfEnergyFunctionS <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n" <>
                                                               "self_energy_PL(i1,i2) = " <> CastIfReal[selfEnergyFunctionPL <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n" <>
                                                               "self_energy_PR(i1,i2) = " <> CastIfReal[selfEnergyFunctionPR <> "(p,i1,i2)", selfEnergyMatrixType] <> ";\n"] <>
                                                    "}\n"
                                                ] <>
                                                "}\n"
                                     ] <>
                                     "}\n"
                                     ,
                                     "const " <> selfEnergyMatrixCType <> " self_energy_1  = " <> CastIfReal[selfEnergyFunctionS  <> "(p)", selfEnergyMatrixType] <> ";\n" <>
                                     "const " <> selfEnergyMatrixCType <> " self_energy_PL = " <> CastIfReal[selfEnergyFunctionPL <> "(p)", selfEnergyMatrixType] <> ";\n" <>
                                     "const " <> selfEnergyMatrixCType <> " self_energy_PR = " <> CastIfReal[selfEnergyFunctionPR <> "(p)", selfEnergyMatrixType] <> ";\n"
                                    ] <>
                                  If[topTwoLoop,
                                     selfEnergyMatrixCType <> " delta_M(- self_energy_PR * M_tree " <>
                                     "- M_tree * self_energy_PL - self_energy_1);\n" <>
                                     "delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l + qcd_4l);\n"
                                     ,
                                     "const " <> selfEnergyMatrixCType <> " delta_M(- self_energy_PR * M_tree " <>
                                     "- M_tree * self_energy_PL - self_energy_1);\n"
                                    ]
                                 ];
              If[IsMajoranaFermion[particle],
                 result = result <>
                          IndentText["const " <> selfEnergyMatrixCType <> " M_loop(M_tree + 0.5 * (delta_M + delta_M.transpose()));\n"];
                 ,
                 result = result <>
                          IndentText["const " <> selfEnergyMatrixCType <> " M_loop(M_tree + delta_M);\n"];
                ];
              result = result <>
                       IndentText[eigenArrayType <> " eigen_values;\n"];
              If[Head[mixingMatrix] === List,
                 (* two mixing matrixs => SVD *)
                 U = ToValidCSymbolString[mixingMatrix[[1]]];
                 V = ToValidCSymbolString[mixingMatrix[[2]]];
                 result = result <>
                          IndentText["decltype(" <> U <> ") mix_" <> U <> ";\n" <>
                                     "decltype(" <> V <> ") mix_" <> V <> ";\n"];
                 result = result <>
                          TreeMasses`CallSVDFunction[
                              particle, "M_loop", "eigen_values",
                              "mix_" <> U, "mix_" <> V];
                 result = result <>
                          IndentText["if (es == 0) {\n" <>
                                     IndentText["PHYSICAL(" <> U <> ") = mix_" <> U <> ";\n" <>
                                                "PHYSICAL(" <> V <> ") = mix_" <> V <> ";\n"] <>
                                     "}\n"
                                    ];
                 ,
                 U = ToValidCSymbolString[mixingMatrix];
                 If[mixingMatrix =!= Null,
                    result = result <>
                             IndentText["decltype(" <> U <> ") mix_" <> U <> ";\n" <>
                                        TreeMasses`CallDiagonalizeSymmetricFunction[
                                            particle, "M_loop", "eigen_values",
                                            "mix_" <> U] <>
                                        "if (es == 0)\n" <>
                                        IndentText["PHYSICAL(" <> U <> ") = mix_" <> U <> ";\n"]
                                       ];
                    ,
                    mixingMatrixType = CreateCType[CConversion`MatrixType[CConversion`complexScalarCType, dim, dim]];
                    result = result <>
                             IndentText[mixingMatrixType <> " mix_" <> U <> ";\n" <>
                                        TreeMasses`CallDiagonalizeSymmetricFunction[
                                            particle, "M_loop", "eigen_values",
                                            "mix_" <> U]];
                   ];
                ];
              result = result <>
                       IndentText["PHYSICAL(" <> massName <>
                                  "(es)) = Abs(eigen_values(es));\n"];
              result = result <> "}\n";
              ,
              (* for a dimension 1 fermion it plays not role if it's a
                 Majorana fermion or not *)

              (* Note: For a 1-dimensional fermion multiplet SARAH
                 provides the self-energies in mass eigenstates, i.e.
                 the fermions at the external legs are multiplied by
                 their phase (= mixing matrix).  Therfore, M_tree must
                 be set to the (positive) tree-level mass.  M_tree
                 must not be set to the gauge eigenstate mass
                 parameter! *)

              result = "const " <> selfEnergyMatrixCType <> " M_tree(" <> massName <> ");\n" <>
                       Do1DimFermion[particle, "M_tree", selfEnergyFunctionS,
                                     selfEnergyFunctionPL, selfEnergyFunctionPR,
                                     momentum, CConversion`GetScalarElementType[selfEnergyMatrixType]];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol /; IsVector[particle], inputMomentum_, __] :=
    Module[{result, dim, dimStr, massName, particleName, mixingMatrix, selfEnergyFunction,
            momentum = inputMomentum, selfEnergyMatrixType, selfEnergyMatrixCType},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleName = ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[inputMomentum == "", momentum = massName];
           mixingMatrix = ToValidCSymbolString[FindMixingMatrixSymbolFor[particle]];
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle, 1];
           selfEnergyMatrixType = TreeMasses`GetMassMatrixType[particle];
           selfEnergyMatrixCType = CreateCType[selfEnergyMatrixType];
           If[IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0,
              If[dim == 1,
                 Return["PHYSICAL(" <> massName <> ") = 0.;\n"];,
                 Return["PHYSICAL(" <> massName <> ").setConstant(0.);\n"];
                ];
             ];
           If[dim > 1,
              result = "WARNING(\"diagonalization of " <> ToString[particle] <> " not implemented\");\n";
              ,
              result = "const " <> selfEnergyMatrixCType <> " M_tree(Sqr(" <> massName <> "));\n" <>
                       Do1DimVector[particleName, massName, "M_tree", selfEnergyFunction, momentum];
             ];
           Return[result];
          ];

DoMediumDiagonalization[particle_Symbol, __] := (
    Print["Error: medium diagonalization of ", particle, " not implemented"];
    "ERROR(\"medium diagonalization of " <> ToString[particle] <> " not implemented\");\n"
);


(* ********** high precision diagonalization routines ********** *)

CalcEffPot2L[particle_] :=
    Module[{dim = GetDimension[particle], dimStr, selfEnergyMatrixCType},
           dimStr = ToString[dim];
           selfEnergyMatrixCType = CreateCType[TreeMasses`GetMassMatrixType[particle]];
           "
// two-loop Higgs self-energy contributions
" <> selfEnergyMatrixCType <> " self_energy_2l(" <> selfEnergyMatrixCType <> "::Zero());

if (pole_mass_loop_order > 1) {
" <> IndentText["\
self_energy_2l = self_energy_" <> CConversion`ToValidCSymbolString[particle] <> "_2loop();
for (int i = 0; i < " <> dimStr <> "; i++) {
   for (int k = 0; k < " <> dimStr <> "; k++) {
      if (!std::isfinite(self_energy_2l(i,k))) {
         self_energy_2l(i,k) = 0.;
         problems.flag_bad_mass(" <> FlexibleSUSY`FSModelName <> "_info::" <> CConversion`ToValidCSymbolString[particle] <> ");
      }
   }
}
"] <> "}
"
          ];

CalcEffPot3L[particle_] :=
    Module[{selfEnergyMatrixCType},
           selfEnergyMatrixCType = CreateCType[TreeMasses`GetMassMatrixType[particle]];
           "
// three-loop Higgs self-energy contributions
" <> selfEnergyMatrixCType <> " self_energy_3l(" <> selfEnergyMatrixCType <> "::Zero());

try {
   if (pole_mass_loop_order > 2)
   " <> IndentText["self_energy_3l = self_energy_" <> CConversion`ToValidCSymbolString[particle] <> "_3loop();"] <> "
} catch (const flexiblesusy::Error& e) {
   WARNING(\"3-loop Higgs mass calculation failed: \" << e.what_detailed());
   problems.flag_bad_mass(" <> FlexibleSUSY`FSModelName <> "_info::" <> CConversion`ToValidCSymbolString[particle] <> ");
}
"
          ];


DoSlowDiagonalization[particle_Symbol, tadpole_] :=
    Module[{result, dim, dimStr, massName, inputMomenta, outputMomenta,
            body, particleStr, effPot = ""},
           dim = GetDimension[particle];
           dimStr = ToString[dim];
           particleStr = CConversion`ToValidCSymbolString[particle];
           massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           inputMomenta = "old_" <> massName;
           outputMomenta = "new_" <> massName;
           If[dim > 1 &&
              (SARAH`UseHiggs2LoopMSSM === True ||
               FlexibleSUSY`UseHiggs2LoopNMSSM === True) &&
              MemberQ[{SARAH`HiggsBoson, SARAH`PseudoScalar}, particle],
              effPot = effPot <> CalcEffPot2L[particle];
             ];
           If[dim > 1 &&
              FlexibleSUSY`UseHiggs3LoopMSSM === True &&
              MemberQ[{SARAH`HiggsBoson}, particle],
              effPot = effPot <> CalcEffPot3L[particle];
             ];
           body = DoMediumDiagonalization[particle, inputMomenta, tadpole, effPot === ""] <> "\n" <>
                  outputMomenta <> " = PHYSICAL(" <> massName <> ");\n" <>
                  "diff = MaxRelDiff(" <> outputMomenta <> ", " <> inputMomenta <> ");\n" <>
                  inputMomenta <> " = " <> outputMomenta <> ";\n" <>
                  "iteration++;\n";
           result = "const auto number_of_mass_iterations = get_number_of_mass_iterations();\n" <>
                    "int iteration = 0;\n" <>
                    "double diff = 0.0;\n" <>
                    "decltype(" <> massName <> ") " <>
                    inputMomenta  <> "(" <> massName <> "), " <>
                    outputMomenta <> "(" <> massName <> ");\n" <>
                    effPot <> "\n" <>
                    "do {\n" <>
                    IndentText[body] <>
                    "\
} while (diff > precision
         && iteration < number_of_mass_iterations);

if (diff > precision)
   problems.flag_no_pole_mass_convergence(" <> FlexibleSUSY`FSModelName <> "_info::" <> particleStr <> ");
else
   problems.unflag_no_pole_mass_convergence(" <> FlexibleSUSY`FSModelName <> "_info::" <> particleStr <> ");
";
           Return[result];
          ];


DoDiagonalization[particle_Symbol, FlexibleSUSY`LowPrecision, tadpole_] :=
    "// diagonalization with low precision\n" <> DoFastDiagonalization[particle, tadpole];

DoDiagonalization[particle_Symbol, FlexibleSUSY`MediumPrecision, tadpole_] :=
    "// diagonalization with medium precision\n" <> DoMediumDiagonalization[particle, "", tadpole, True];

DoDiagonalization[particle_Symbol, FlexibleSUSY`HighPrecision, tadpole_] :=
    "// diagonalization with high precision\n" <> DoSlowDiagonalization[particle, tadpole];

CreateLoopMassFunctionName[particle_Symbol] :=
    "calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <> "_pole";

CallPoleMassFunction[particle_Symbol, obj_:""] :=
    obj <> CreateLoopMassFunctionName[particle] <> "();\n";

CreateLoopMassPrototype[particle_Symbol] :=
    "void " <> CreateLoopMassFunctionName[particle] <> "();\n";

CreateLoopMassFunction[particle_Symbol, precision_Symbol, tadpole_] :=
    Module[{result, body = ""},
           If[!IsFermion[particle] &&
              !(IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0),
              body = "if (!force_output && problems.is_running_tachyon(" <> FlexibleSUSY`FSModelName <> "_info::" <> ToValidCSymbolString[particle] <> "))\n" <>
                     IndentText["return;"] <> "\n\n";
             ];
           body = body <> DoDiagonalization[particle, precision, tadpole];
           result = "void CLASSNAME::" <> CreateLoopMassFunctionName[particle] <>
                    "()\n{\n" <> IndentText[body] <> "}\n\n";
           Return[result];
          ];

(* return pole mass of a singlet as a function of p *)
Create1DimPoleMassPrototype[particle_Symbol] :=
    If[GetDimension[particle] > 1,
       Print["Warning: cannot generate extra pole mass"
             " calculation function for ", particle, ", because"
             " it has more than 1 generation"];
       "",
       "double " <> CreateLoopMassFunctionName[particle] <> "(double);\n"
      ];

(* return pole mass of a singlet as a function of p *)
Create1DimPoleMassFunction[particle_Symbol] :=
    Module[{result, body = "", particleName, massName, mTree},
           If[GetDimension[particle] > 1,
              Print["Warning: cannot generate extra pole mass"
                    " calculation function for ", particle, ", because"
                    " it has more than 1 generation"];
              Return[""];
             ];
           If[!(IsUnmixed[particle] && GetMassOfUnmixedParticle[particle] === 0),
              body = "if (!force_output && problems.is_running_tachyon(" <> FlexibleSUSY`FSModelName <> "_info::" <> ToValidCSymbolString[particle] <> "))\n" <>
                     IndentText["return 0.;"] <> "\n\n";
             ];
           If[!IsMassless[particle],
               particleName = ToValidCSymbolString[particle];
               massName = ToValidCSymbolString[FlexibleSUSY`M[particle]];
               selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle, 1];
               (* vector bosons are always unmixed -> make sure the
                  right mass eigenvalue is used.  The mass matrix might
                  contain mixing matrix elements, which can result in
                  the wrong mass eigenstate. *)
               If[IsVector[particle],
                  mTree = "Sqr(" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]] <> ")",
                  mTree = "get_mass_matrix_" <> particleName <> "()";
                 ];
               body = body <>
                      "const double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
                      "const double mass_sqr = " <> mTree <> " - self_energy;\n\n" <>
                      "if (mass_sqr < 0.) {\n" <>
                      IndentText[TreeMasses`FlagPoleTachyon[particleName]] <> "\n}\n\n" <>
                      "return AbsSqrt(mass_sqr);\n";
              ,
              body = "return 0.;\n";
             ];
           result = "double CLASSNAME::" <> CreateLoopMassFunctionName[particle] <>
                    "(double p)\n{\n" <> IndentText[body] <> "}\n\n";
           Return[result];
          ];

CreateOneLoopPoleMassFunctions[precision_List, oneLoopTadpoles_List, vevs_List] :=
    Module[{result = "", particle, prec, i = 1, f, d, tadpole, fieldsAndVevs = {}},
           (* create list that associates fields at vevs *)
           For[f = 1, f <= Length[oneLoopTadpoles], f++,
               particle = GetField[oneLoopTadpoles[[f]]];
               For[d = 1, d <= GetDimension[particle], d++; i++,
                   AppendTo[fieldsAndVevs, {particle, d, vevs[[i]]}];
                  ];
              ];
           For[i = 1, i <= Length[precision], i++,
               particle = precision[[i,1]];
               prec     = precision[[i,2]];
               tadpole  = Cases[fieldsAndVevs, {particle[_], __}];
               result   = result <> CreateLoopMassFunction[particle, prec, tadpole];
              ];
           If[ValueQ[SARAH`VectorW],
              result = result <> Create1DimPoleMassFunction[SARAH`VectorW];
             ];
           If[ValueQ[SARAH`VectorZ],
              result = result <> Create1DimPoleMassFunction[SARAH`VectorZ];
             ];
           Return[result];
          ];

CreateOneLoopPoleMassPrototypes[states_:FlexibleSUSY`FSEigenstates] :=
    Module[{particles, result = ""},
           particles = GetLoopCorrectedParticles[states];
           (result = result <> CreateLoopMassPrototype[#])& /@ particles;
           If[ValueQ[SARAH`VectorW],
              result = result <> Create1DimPoleMassPrototype[SARAH`VectorW];
             ];
           If[ValueQ[SARAH`VectorZ],
              result = result <> Create1DimPoleMassPrototype[SARAH`VectorZ];
             ];
           Return[result];
          ];

CallThreadedPoleMassFunction[particle_Symbol, ptr_:"this", pool_:"tp"] :=
    pool <> ".run_task([" <> ptr <> "] () { " <>
    If[ptr === "this", "", ptr <> "->"] <>
    CreateLoopMassFunctionName[particle] <> "(); });\n";

CallAllPoleMassFunctions[states_, enablePoleMassThreads_] :=
    Module[{particles, susyParticles, smParticles, callSusy,
            callSM, result},
           particles = GetLoopCorrectedParticles[states];
           smParticles = Select[particles, TreeMasses`IsSMParticle];
           susyParticles = Complement[particles, smParticles];
           If[enablePoleMassThreads =!= True,
              callSusy = StringJoin[CallPoleMassFunction /@ susyParticles];
              callSM   = StringJoin[CallPoleMassFunction /@ smParticles];
              result = "if (calculate_bsm_pole_masses) {\n" <>
                       IndentText[callSusy] <>
                       "}\n\n" <>
                       "if (calculate_sm_pole_masses) {\n" <>
                       IndentText[callSM] <>
                       "}\n";
              ,
              callSusy = StringJoin[CallThreadedPoleMassFunction /@ susyParticles];
              callSM   = StringJoin[CallThreadedPoleMassFunction /@ smParticles];
              result = "Thread_pool tp(std::min(std::thread::hardware_concurrency(), " <> ToString[Length[susyParticles] + Length[smParticles]] <> "u));\n\n" <>
                       "if (calculate_bsm_pole_masses) {\n" <>
                       IndentText[callSusy] <>
                       "}\n\n" <>
                       "if (calculate_sm_pole_masses) {\n" <>
                       IndentText[callSM] <>
                       "}\n";
             ];
           result
          ];

GetRunningOneLoopDRbarParticles[] :=
    Module[{downQuarks, upQuarks, downLeptons, upLeptons},
           upLeptons   = TreeMasses`GetSMNeutralLeptons[];
           downLeptons = TreeMasses`GetSMChargedLeptons[];
           upQuarks    = TreeMasses`GetSMUpQuarks[];
           downQuarks  = TreeMasses`GetSMDownQuarks[];
           Flatten[{upLeptons, downLeptons, upQuarks, downQuarks,
                    SARAH`VectorP, SARAH`VectorZ, SARAH`VectorW}]
          ];

(* returns conversion factor from MS-bar scheme to renormalizationScheme *)
GetConversionFactorMSbarTo[particle_,
                           renormalizationScheme_ /; renormalizationScheme === FlexibleSUSY`DRbar,
                           {alphaS_, gWeak_, gPrime_}
                          ] :=
    Which[(* down-type quarks *)
          TreeMasses`IsSMDownQuark[particle],
          (1 - alphaS / (3 Pi)
           + 3 gWeak^2 / (128 Pi^2)
           + 13 gPrime^2 / (1152 Pi^2)),
          (* down-type leptons *)
          TreeMasses`IsSMChargedLepton[particle],
          (1 - 3 (gPrime^2 - gWeak^2) / (128 Pi^2)),
          (* otherwise *)
          True, 1
         ];

GetConversionFactorMSbarTo[_,_,_] := 1;

CallCalculateDRbarMass[splitName_String, multipletName_String, index_Integer, resultMatrix_String, mass_String] :=
    Module[{cIdxStr = ToString[index-1], particle, optionalIndex = ""},
           particle = Parameters`GetParticleFromDescription[splitName];
           If[particle === Null,
              optionalIndex = ", " <> cIdxStr;
              particle = Parameters`GetParticleFromDescription[multipletName];
              ];
           resultMatrix <> "(" <> cIdxStr <> "," <> cIdxStr <> ") = " <>
           "MODEL->calculate_M" <> CConversion`ToValidCSymbolString[particle] <>
           "_DRbar(" <> mass <> optionalIndex <> ");"
          ];

CreateRunningDRbarMassPrototype[particle_ /; IsFermion[particle]] :=
    "double calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <>
    "_DRbar(double" <> If[TreeMasses`GetDimension[particle] > 1, ", int", ""] <> ") const;\n";

CreateRunningDRbarMassPrototype[particle_] :=
    "double calculate_" <> ToValidCSymbolString[FlexibleSUSY`M[particle]] <>
    "_DRbar(double) const;\n";

CreateRunningDRbarMassPrototypes[] :=
    Module[{result = "", particles},
           particles = GetRunningOneLoopDRbarParticles[];
           (result = result <> CreateRunningDRbarMassPrototype[#])& /@ particles;
           Return[result];
          ];

CreateRunningDRbarMassFunction[particle_ /; particle === TreeMasses`GetSMBottomQuarkMultiplet[], renormalizationScheme_] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, alphaS, drbarConversion, gPrime,
            dimParticle, treeLevelMass},
           dimParticle = TreeMasses`GetDimension[particle];
           treeLevelMass = TreeMasses`GetThirdGenerationMass[particle] /. a_[i_?IntegerQ] :> a[Global`idx];
           selfEnergyFunctionS  = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[1], 1];
           selfEnergyFunctionPL = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PL], 1];
           selfEnergyFunctionPR = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PR], 1];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double) const\n{\n";,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
                ];
              body = "return 0.0;\n";
              ,
              alphaS = SARAH`strongCoupling^2/(4 Pi);
              gPrime = SARAH`hyperchargeCoupling /. Parameters`ApplyGUTNormalization[];
              (* convert MSbar to DRbar mass hep-ph/0207126 *)
              drbarConversion = GetConversionFactorMSbarTo[particle, renormalizationScheme, {alphaS, SARAH`leftCoupling, gPrime}];
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar, int idx) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n";
                ];
              body = body <>
              "const double m_tree = " <> RValueToCFormString[treeLevelMass] <> ";\n" <>
              "const double drbar_conversion = " <> RValueToCFormString[drbarConversion] <> ";\n" <>
              "const double m_sm_drbar = m_sm_msbar * drbar_conversion;\n" <>
              "const double delta_mb_1loop = - self_energy_1/m_tree - self_energy_PL - self_energy_PR;\n" <>
              "double qcd_2l = 0.;\n" <>
              AddMbRun2LSQCDCorrections[] <> "\n" <>
              "const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mb_1loop + qcd_2l);\n\n" <>
              "return m_susy_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_ /; TreeMasses`IsSMChargedLepton[particle], renormalizationScheme_] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, drbarConversion, gPrime,
            dimParticle},
           dimParticle = TreeMasses`GetDimension[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[1], 1];
           selfEnergyFunctionPL = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PL], 1];
           selfEnergyFunctionPR = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PR], 1];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double) const\n{\n";,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
                ];
              body = "return 0.0;\n";
              ,
              (* convert MSbar to DRbar mass *)
              gPrime = SARAH`hyperchargeCoupling /. Parameters`ApplyGUTNormalization[];
              drbarConversion = GetConversionFactorMSbarTo[particle, renormalizationScheme, {SARAH`strongCoupling^2/(4 Pi), SARAH`leftCoupling, gPrime}];
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_sm_msbar, int idx) const\n{\n";
                 body = "const double p = m_sm_msbar;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n";
                ];
              body = body <>
              "const double drbar_conversion = " <> RValueToCFormString[drbarConversion] <> ";\n" <>
              "const double m_sm_drbar = m_sm_msbar * drbar_conversion;\n" <>
              "const double delta_mf_1loop = - self_energy_1/m_sm_drbar - self_energy_PL - self_energy_PR;\n\n" <>
              "const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mf_1loop);\n\n" <>
              "return m_susy_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateMSSM1LoopSQCDContributions[result_:"qcd_1l"] := "\
{
" <> IndentText[FillMt2LStruct[]] <> "

   " <> result <> " = - mssm_twoloop_mt::dMt_over_mt_1loop(pars);
}";

CreateMSSM2LoopSQCDContributions[result_:"qcd_2l"] :=
    FillMt2LStruct[] <> "

const double q_2l = mssm_twoloop_mt::dMt_over_mt_2loop(pars);

" <> result <> " = -q_2l + qcd_1l * qcd_1l;";

(* 2L QCD contribution from gluino in the split-MSSM,
   Eq. (4.7) of arxiv:1312.5220
 *)
GetMtPoleOverMtMSbarSplitMSSM2L[particle_, scale_] :=
    Module[{mg = FlexibleSUSY`M[SARAH`Gluino]},
           SARAH`strongCoupling^4/(4 Pi)^4 (
               89/9
               + 4 Log[mg^2/scale^2] (
                   13/3
                   + Log[mg^2/scale^2]
                   - 2 Log[FlexibleSUSY`M[particle]^2/scale^2]
               )
           )
          ];

CreateRunningDRbarMassFunction[particle_ /; particle === TreeMasses`GetSMTopQuarkMultiplet[], _] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, qcdOneLoop, qcdTwoLoop, qcdThreeLoop = 0, qcdFourLoop = 0,
            dimParticle, treeLevelMass},
           dimParticle = TreeMasses`GetDimension[particle];
           treeLevelMass = TreeMasses`GetThirdGenerationMass[particle] /. a_[i_?IntegerQ] :> a[Global`idx];
           selfEnergyFunctionS  = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[1], 1];
           selfEnergyFunctionPL = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PL], 1];
           selfEnergyFunctionPR = SelfEnergies`CreateHeavyRotatedSelfEnergyFunctionName[particle[PR], 1];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double) const\n{\n";,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
                ];
              body = "return 0.0;\n";
              ,
              qcdOneLoop = N[SimpQCD[-TwoLoopQCD`GetDeltaMPoleOverMRunningQCDOneLoop[particle, Global`currentScale, FlexibleSUSY`FSRenormalizationScheme]]];
              qcdTwoLoop = N[SimpQCD[TwoLoopQCD`GetDeltaMPoleOverMRunningQCDTwoLoop[particle, Global`currentScale, FlexibleSUSY`FSRenormalizationScheme]]];
              (* adding split-MSSM contributions if enabled *)
              If[FlexibleSUSY`UseHiggs3LoopSplit === True,
                 qcdTwoLoop = N[SimpQCD[qcdTwoLoop - GetMtPoleOverMtMSbarSplitMSSM2L[particle, Global`currentScale]]];
                ];
              If[FlexibleSUSY`UseYukawa3LoopQCD === True &&
                 FlexibleSUSY`FSRenormalizationScheme =!= FlexibleSUSY`MSbar,
                 Print["Warning: UseYukawa3LoopQCD == True, but the renormalization scheme is not MSbar!"];
                 Print["  The 3-loop QCD corrections to the top Yukawa coupling will be disabled."];
                ];
              If[FlexibleSUSY`UseYukawa4LoopQCD === True &&
                 FlexibleSUSY`FSRenormalizationScheme =!= FlexibleSUSY`MSbar,
                 Print["Warning: UseYukawa4LoopQCD == True, but the renormalization scheme is not MSbar!"];
                 Print["  The 4-loop QCD corrections to the top Yukawa coupling will be disabled."];
                ];
              If[FlexibleSUSY`UseYukawa3LoopQCD === True ||
                 FlexibleSUSY`UseYukawa3LoopQCD === Automatic,
                 If[FlexibleSUSY`FSRenormalizationScheme === FlexibleSUSY`MSbar,
                    qcdThreeLoop = N[Get3LQCDmtOverMt[particle, Global`currentScale]];
                   ];
                ];
              If[FlexibleSUSY`UseYukawa4LoopQCD === True ||
                 FlexibleSUSY`UseYukawa4LoopQCD === Automatic,
                 If[FlexibleSUSY`FSRenormalizationScheme === FlexibleSUSY`MSbar,
                    qcdFourLoop = N[Get4LQCDmtOverMt[particle, Global`currentScale]];
                   ];
                ];
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole, int idx) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n\n";
                ];
              body = body <>
              "const double currentScale = get_scale();\n" <>
              "double qcd_1l = 0., qcd_2l = 0., qcd_3l = 0., qcd_4l = 0.;\n" <>
              "double atas_S_2l = 0., atas_LR_2l = 0., atat_S_2l = 0., atat_LR_2l = 0.;\n\n" <>
                  If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
                     CreateMSSM1LoopSQCDContributions[],
                     "qcd_1l = " <> CConversion`RValueToCFormString[qcdOneLoop /. FlexibleSUSY`M[particle] -> treeLevelMass] <> ";"
                    ] <>
              "\n\n" <>
              "if (get_thresholds() > 1 && threshold_corrections.mt > 1) {\n" <>
              IndentText[
                  If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
                     CreateMSSM2LoopSQCDContributions[],
                     "const double q_2l = " <> CConversion`RValueToCFormString[qcdTwoLoop /. FlexibleSUSY`M[particle] -> treeLevelMass] <> ";\n\n" <>
                     "qcd_2l = -q_2l + qcd_1l * qcd_1l;"
                  ] <>
                  If[FlexibleSUSY`UseSMYukawa2Loop === True,
                     "

const double gs = " <> CConversion`RValueToCFormString[SARAH`strongCoupling] <> ";
const double yt = " <> CConversion`RValueToCFormString[Parameters`GetThirdGeneration[SARAH`UpYukawa]] <> ";
const double t = Sqr(" <> CConversion`RValueToCFormString[treeLevelMass] <> ");
const double h = Sqr(" <> CConversion`RValueToCFormString[TreeMasses`GetMass[SARAH`HiggsBoson]] <> ");
const double s = Sqr(p);
const double qq = Sqr(get_scale());

atas_S_2l  = sm_twoloop_mt::delta_mt_2loop_as_at_S_flexiblesusy(gs, yt, t, h, s, qq);
atas_LR_2l = sm_twoloop_mt::delta_mt_2loop_as_at_LR_flexiblesusy(gs, yt, t, h, s, qq);

atat_S_2l  = sm_twoloop_mt::delta_mt_2loop_at_at_S_flexiblesusy(yt, t, h, s, qq);
atat_LR_2l = sm_twoloop_mt::delta_mt_2loop_at_at_LR_flexiblesusy(yt, t, h, s, qq);"
                  , ""]
              ] <> "\n" <>
              "}\n\n" <>
              If[qcdThreeLoop =!= 0,
              "if (get_thresholds() > 2 && threshold_corrections.mt > 2) {\n" <>
                  IndentText["qcd_3l = " <> CConversion`RValueToCFormString[qcdThreeLoop /. FlexibleSUSY`M[particle] -> treeLevelMass] <> ";"] <> "\n" <>
              "}\n\n", ""] <>
              If[qcdFourLoop =!= 0,
              "if (get_thresholds() > 3 && threshold_corrections.mt > 3) {\n" <>
                  IndentText["qcd_4l = " <> CConversion`RValueToCFormString[qcdFourLoop /. FlexibleSUSY`M[particle] -> treeLevelMass] <> ";"] <> "\n" <>
              "}\n\n", ""] <>
              "const double m_susy_drbar = m_pole + self_energy_1 + atas_S_2l + atat_S_2l " <>
              "+ m_pole * (self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l + qcd_4l + atas_LR_2l + atat_LR_2l);\n\n" <>
              "return m_susy_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_ /; IsFermion[particle], _] :=
    Module[{result, body, selfEnergyFunctionS, selfEnergyFunctionPL,
            selfEnergyFunctionPR, name, dimParticle},
           dimParticle = TreeMasses`GetDimension[particle];
           selfEnergyFunctionS  = SelfEnergies`CreateSelfEnergyFunctionName[particle[1], 1];
           selfEnergyFunctionPL = SelfEnergies`CreateSelfEnergyFunctionName[particle[PL], 1];
           selfEnergyFunctionPR = SelfEnergies`CreateSelfEnergyFunctionName[particle[PR], 1];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double) const\n{\n";,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double, int) const\n{\n";
                ];
              body = "return 0.0;\n";
              ,
              If[dimParticle == 1,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p));\n";
                 ,
                 result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole, int idx) const\n{\n";
                 body = "const double p = m_pole;\n" <>
                 "const double self_energy_1  = Re(" <> selfEnergyFunctionS  <> "(p, idx, idx));\n" <>
                 "const double self_energy_PL = Re(" <> selfEnergyFunctionPL <> "(p, idx, idx));\n" <>
                 "const double self_energy_PR = Re(" <> selfEnergyFunctionPR <> "(p, idx, idx));\n";
                ];
              body = body <> "\n" <>
                     "const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL + self_energy_PR)";
              body = body <> ";\n\n";
              body = body <> "return m_drbar;\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunction[particle_, _] :=
    Module[{result, body, selfEnergyFunction, name, particleName},
           selfEnergyFunction = SelfEnergies`CreateSelfEnergyFunctionName[particle, 1];
           particleName = ToValidCSymbolString[particle];
           name = ToValidCSymbolString[FlexibleSUSY`M[particle]];
           If[IsMassless[particle],
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double) const\n{\n";
              body = "return 0.0;\n";
              ,
              result = "double CLASSNAME::calculate_" <> name <> "_DRbar(double m_pole) const\n{\n";
              body = "const double p = m_pole;\n" <>
              "const double self_energy = Re(" <> selfEnergyFunction <> "(p));\n" <>
              "const double mass_sqr = Sqr(m_pole) + self_energy;\n\n" <>
              "if (mass_sqr < 0.) {\n" <>
              IndentText[TreeMasses`FlagPoleTachyon[particleName] <>
                         "return m_pole;"] <> "\n}\n\n" <>
              "return AbsSqrt(mass_sqr);\n";
             ];
           Return[result <> IndentText[body] <> "}\n\n"];
          ];

CreateRunningDRbarMassFunctions[renormalizationScheme_:FlexibleSUSY`DRbar] :=
    Module[{result = "", particles},
           particles = GetRunningOneLoopDRbarParticles[];
           (result = result <> CreateRunningDRbarMassFunction[#,renormalizationScheme])& /@ particles;
           Return[result];
          ];

GetLightestMassEigenstate[FlexibleSUSY`M[mass_]] :=
    GetLightestMassEigenstate[mass];

GetLightestMassEigenstate[mass_] :=
    If[GetDimension[mass] == 1,
       FlexibleSUSY`M[mass],
       FlexibleSUSY`M[mass][GetDimensionStartSkippingGoldstones[mass] - 1]];

CreateLSPFunctions[{}] := {"", ""};

CreateLSPFunctions[masses_List] :=
    Module[{prototype, function, mass, info, particleType, m,
            comment},
           info = FlexibleSUSY`FSModelName <> "_info";
           particleType = info <> "::Particles";
           body = "double lsp_mass = std::numeric_limits<double>::max();
double tmp_mass;
particle_type = " <> info <> "::NUMBER_OF_PARTICLES;

";
           For[m = 1, m <= Length[masses], m++,
               mass = masses[[m]];
               body = body <> "\
tmp_mass = Abs(PHYSICAL(" <>
CConversion`RValueToCFormString[GetLightestMassEigenstate[mass]] <>
"));
if (tmp_mass < lsp_mass) {
" <> IndentText["\
lsp_mass = tmp_mass;
particle_type = " <> info <> "::" <>
CConversion`ToValidCSymbolString[mass /. FlexibleSUSY`M -> Identity] <>
";"] <>
"
}

";
              ];
           body = body <> "return lsp_mass;\n";
           prototype = "double get_lsp(" <> particleType <> "&) const;\n";
           comment = "\
/**
 * @brief finds the LSP and returns it's mass
 *
 * This function finds the lightest supersymmetric particle (LSP) and
 * returns it's mass.  The corresponding particle type is retured in
 * the reference parameter.  The list of potential LSPs is set in the
 * model file varible PotentialLSPParticles.  For this model it is set
 * to:
 * " <> ToString[masses] <> "
 *
 * @param particle_type particle type
 * @return mass of LSP
 */
";
           function = comment <>
                      "double CLASSNAME::get_lsp(" <> particleType <>
                      "& particle_type) const\n{\n" <> IndentText[body] <> "}\n";
           Return[{prototype, function}];
          ];

CreateLSPFunctions[_] := {"", ""};

End[];

EndPackage[];
