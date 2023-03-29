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

BeginPackage["FlexibleEFTHiggsMatching`", {"CConversion`", "TreeMasses`", "LoopMasses`", "Constraint`", "ThresholdCorrections`", "Parameters`", "Utils`"}];

CalculateMHiggsPoleOneMomentumIteration::usage = "";
Create3LoopMatching::usage = "";
CallMatch2LoopTopMass::usage = "";
CreateMt2Loop::usage = "";
CalculateRunningUpQuarkMasses::usage = "";
CalculateRunningDownQuarkMasses::usage = "";
CalculateRunningDownLeptonMasses::usage = "";
CalculateMUpQuarkPole1L::usage = "";
CalculateMDownQuarkPole1L::usage = "";
CalculateMDownLeptonPole1L::usage = "";
FillSMFermionPoleMasses::usage = "";
GetFixedBSMParameters::usage="Returns a list of the BSM parameters fixed by matching SM -> BSM.";
SetBSMParameters::usage = "";
SetGaugeLessLimit::usage = "";

Begin["`Private`"];

CalculateMHiggsPoleOneMomentumIteration[particle_] :=
    If[GetDimension[particle] == 1,
       "Mh2_pole = mh2_tree - self_energy - tadpole(0);"
       ,
"const auto M_loop = (mh2_tree - self_energy - " <> CreateCType[TreeMasses`GetMassMatrixType[particle]] <> "(tadpole.asDiagonal())).eval();

" <> CreateCType[CConversion`ArrayType[CConversion`realScalarCType, GetDimension[particle]]] <> " M_pole;
fs_diagonalize_hermitian(M_loop, M_pole);

Mh2_pole = M_pole(idx);"
      ];

Create3LoopMatching[] :=
Module[{g1str = CConversion`ToValidCSymbolString[SARAH`hyperchargeCoupling],
g2str = CConversion`ToValidCSymbolString[SARAH`leftCoupling],
g3str = ToString[TreeMasses`GetStrongCoupling[]]},
"
   using namespace flexiblesusy::sm_twoloophiggs;
      
   static const double gauge_less = 1e-10; 
   // approximates the gaugeless limit for gauge couplings

   auto model = model_input;
   auto model_gl = model;
   model_gl.get_problems().clear();
   model_gl.set_" <> g1str <> "(gauge_less / " <> ToString[FlexibleSUSY`FSModelName] <> "_info::normalization_g1);
   model_gl.set_" <> g2str <> "(gauge_less / " <> ToString[FlexibleSUSY`FSModelName] <> "_info::normalization_g2);
   auto model_no_g3 = model_gl; 
   model_no_g3.set_" <> g3str <> "(gauge_less / " <> ToString[FlexibleSUSY`FSModelName] <> "_info::normalization_g3);
   model_gl.calculate_DRbar_masses();
   model_no_g3.calculate_DRbar_masses();
   model.calculate_DRbar_masses();

   auto sm_0l = sm;
   auto sm_1l = sm;
   auto sm_0l_gl = sm;
   auto sm_1l_gl = sm;
   auto sm_0l_no_g3 = sm;
   auto sm_1l_gl_g3less = sm;
   auto sm_2l = sm;

   match_high_to_low_scale_sm_0l(sm_0l, model, idx);
   match_high_to_low_scale_sm_1l(sm_1l, model, idx);
   match_high_to_low_scale_sm_0l(sm_0l_gl, model_gl, idx);
   match_high_to_low_scale_sm_1l(sm_1l_gl, model_gl, idx);
   match_high_to_low_scale_sm_1l(sm_1l_gl_g3less, model_no_g3, idx);
   match_high_to_low_scale_sm_2l(sm_2l, model, idx);

   sm = sm_2l; 

   // calculation of 3-loop threshold corrections below
      
   const double lambda_2l = sm_2l.get_Lambdax();
   double delta_lambda_3l = 0.;

   const double v2 = Sqr(sm_0l_gl.get_v());
   const double yt = sm_0l_gl.get_Yu(2, 2);
   const double yt2 = Sqr(yt);
   const double gs = sm_0l_gl.get_g3();
   const double gs2 = Sqr(gs);
   const double gs4 = Sqr(gs2);
   const double MS = Sqrt(model_gl.get_MSt(0) * model_gl.get_MSt(1));
   const double Q2 = Sqr(sm_0l_gl.get_scale());
   const double mt = model_gl.get_MFt();
   const double mt2 = Sqr(mt);
   const double logmt = Log(mt2 / Q2);
   const double logmt2 = Sqr(Log(mt2 / Q2));
   const double logmt3 = logmt * logmt2;
   const double k = oneOver16PiSqr;
   const double k2 = twoLoop;
   const double k3 = threeLoop;

   const double delta_yt_1l = sm_1l_gl.get_Yu(2, 2) - sm_1l_gl_g3less.get_Yu(2, 2);
   const double delta_yt_2l = sm_2l.get_Yu(2, 2) - sm_1l.get_Yu(2, 2);
   const double sqr_delta_yt_1l = Sqr(delta_yt_1l);
   const double delta_g3_1l = sm_1l.get_g3() - sm_0l_gl.get_g3();

   const double S1_deriv_yt = delta_mh_1loop_at_sm_deriv_yt(0, sm_0l_gl.get_scale(),
                                                     sm_0l_gl.get_MFu(2),
                                                     sm_0l_gl.get_Yu(2, 2)) * delta_yt_2l;
  
   const double S1_deriv_yt2 =
      -0.5 * (24 * k * mt2 * (7. + 6. * logmt)) * sqr_delta_yt_1l;
   const double S2_deriv_yt =
      (64. * k2 * gs2 * mt2 * yt * (1. + 8. * logmt + 6. * logmt2)) *
      delta_yt_1l;
   const double S2_deriv_gs =
      (64. * k2 * gs * mt2 * yt2 * logmt * (1. + 3. * logmt)) * delta_g3_1l;

   const double mh2_conversion =
      S1_deriv_yt + S1_deriv_yt2 + S2_deriv_yt + S2_deriv_gs;           // sum in eq.(4.26c) JHEP07(2020)197
   const double mh2_sm_shift = -Re(sm_0l_gl.self_energy_hh_3loop());

   // calculate 2L self-energy using tree-level parameters
   Eigen::Matrix<double, 2, 2> self_energy_3l(
      Eigen::Matrix<double, 2, 2>::Zero());

   try {
      self_energy_3l = Re(model_gl.self_energy_hh_3loop());

      const auto Mh2_loop = (calculate_mh2_0l(model_gl) - self_energy_3l).eval();

      Eigen::Array<double, 2, 1> Mh2_pole;
      fs_diagonalize_hermitian(Mh2_loop, Mh2_pole);
      const double mh2_bsm_shift = Mh2_pole(idx) - Sqr(model_gl.get_Mhh(idx));

      delta_lambda_3l = (mh2_bsm_shift - mh2_sm_shift - mh2_conversion) / v2;   // eq. (4.28d) JHEP07(2020)197
   } catch (const flexiblesusy::Error &e){}

   const double lambda_3l = lambda_2l + delta_lambda_3l;
   sm.set_Lambdax(lambda_3l);
   sm.calculate_DRbar_masses();    
   "
];


CallMatch2LoopTopMass[] :=
Module[{},
"
const double delta_vev = sm_1l.get_v() - sm_0l.get_v();
const double yt_0l = sm_0l.get_Yu(2,2);
const double mt_2l = calculate_MFt_MSbar_2loop(sm_0l, model);
sm.set_Yu(2, 2, (Sqrt(2) * mt_2l / v - yt_0l * delta_vev / v));
"
];

CreateMt2Loop[] :=
Module[{g3str =  ToString[TreeMasses`GetStrongCoupling[]], 
mglustr = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]],
md2str = CConversion`RValueToCFormString[SARAH`SoftDown],
mq2str = CConversion`RValueToCFormString[SARAH`SoftSquark],
mtstr = CConversion`RValueToCFormString[TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]]] },
"double calculate_MFt_MSbar_2loop(
   const standard_model::Standard_model &sm_0l,
   const " <> ToString[FlexibleSUSY`FSModelName] <> "_mass_eigenstates &model_0l) {
   auto model = model_0l;
   auto sm = sm_0l;

   sm.set_pole_mass_loop_order(2);
   sm.set_top_QCD_order(1);
   model.set_pole_mass_loop_order(2);
   model.set_top_QCD_order(1);

   double mst_1, mst_2, theta_t;
   model."<> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>";

   const double Q2 = Sqr(model.get_scale());
   const double mt = sm_0l.get_MFu(2);   // is equal to the top quark mass in the MSSM
   const double mt2 = Sqr(mt);
   const double gs = model.get_" <> g3str <> "();
   const double gs2 = Sqr(gs);
   const double alpha_s = gs2 / (4. * Pi);

   const double k = oneOver16PiSqr;
   const double logmt = Log(mt2 / Q2);
   const double yt_SM = sm.get_Yu(2, 2);

   mssm_twoloop_mt::Parameters pars;
   pars.g3 = gs;
   pars.mt = mt;
   pars.mg = model.get_" <> mglustr <> "();
   pars.mst1 = mst_1;
   pars.mst2 = mst_2;
   pars.msusy = Sqrt(Sqrt(Abs(model.get_" <> md2str <> "(2, 2)))) *
                Sqrt(Sqrt(Abs(model.get_" <> mq2str <> "(2, 2))));
   pars.xt = Sin(2 * theta_t) * (Sqr(mst_1) - Sqr(mst_2)) / (2. * pars.mt);
   pars.Q = model.get_scale();

   const double delta_alpha_s = calculate_delta_alpha_s(alpha_s, model_0l);

   // top mass contribution at O(as) in MSSM and SM
   const double sqcd_1l_over_mt = mssm_twoloop_mt::dMt_over_mt_1loop(pars);
   const double qcd_1l_SM_over_mt =
      0.008443431970194815 * (4. - 3. * Log(mt2 / Q2)) * gs2;

   // top mass contribution at O(as) SM subscript
   const double qcd1l_mt = -4. / 3. * k * gs2 * mt * (2. + 3. * logmt);
   const double qcd1l_gs = 8. / 3. * k * gs2 * mt * (4. - 3. * logmt);

   const double delta_mt_over_mt = sqcd_1l_over_mt - qcd_1l_SM_over_mt;
   const double delta_gs_over_gs = -0.5 * delta_alpha_s;

   const double mfcon = delta_mt_over_mt * qcd1l_mt + delta_gs_over_gs * qcd1l_gs;

   model."<> CreateLoopMassFunctionName[TreeMasses`GetUpQuark[3,True]] <>"();
   sm.calculate_MFu_pole();

   const auto Mf_bsm = model.get_physical()."<> mtstr <>";
   const auto Mf_sm = sm.get_physical().MFu(2);
   const auto mf_sm_0l = sm.get_MFu(2);

   const double mf_sm = Mf_bsm - Mf_sm + mf_sm_0l - mfcon;   // eq. (4.18b)  JHEP07(2020)197

   return Abs(mf_sm);
}
      "
      ];



CalculateRunningUpQuarkMasses[] :=
    Module[{result = "", i, iStr, mq, mqFun},
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetUpQuark[i + 1, True]];
               If[Length[TreeMasses`GetSMUpQuarks[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "mf(" <> iStr <> "," <> iStr <> ") = model.get_M" <> mqFun <> ";\n";
              ];
           result
          ];

CalculateRunningDownQuarkMasses[] :=
    Module[{result = "", i, iStr, mq, mqFun},
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetDownQuark[i + 1, True]];
               If[Length[TreeMasses`GetSMDownQuarks[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "mf(" <> iStr <> "," <> iStr <> ") = model.get_M" <> mqFun <> ";\n";
              ];
           result
          ];

CalculateRunningDownLeptonMasses[] :=
    Module[{result = "", i, iStr, mq, mqFun},
           For[i = 0, i < 3, i++,
               iStr = ToString[i];
               mq = mqFun = CConversion`RValueToCFormString[TreeMasses`GetDownLepton[i + 1, True]];
               If[Length[TreeMasses`GetSMChargedLeptons[]] == 3, mqFun = mq <> "()"];
               result = result <>
                        "mf(" <> iStr <> "," <> iStr <> ") = model.get_M" <> mqFun <> ";\n";
              ];
           result
          ];

CalculateMFermionPole1L[name_String, GetList_, GetEntry_] :=
    Module[{result = "", i, iStr, mq, mqFun},
           If[Length[GetList[]] == 3,
              For[i = 0, i < 3, i++,
                  mq = CConversion`RValueToCFormString[GetEntry[i + 1, True]];
                  iStr = ToString[i];
                  result = result <>
"\
if (i == " <> iStr <> ") {
   const double p = model_0l.get_M" <> mq <> "();
   const auto self_energy_1  = Re(model_0l.self_energy_" <> mq <> "_1loop_1(p));
   const auto self_energy_PL = Re(model_0l.self_energy_" <> mq <> "_1loop_PL(p));
   const auto self_energy_PR = Re(model_0l.self_energy_" <> mq <> "_1loop_PR(p));
   const auto M_tree = model_0l.get_mass_matrix_" <> mq <> "();
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PR + self_energy_PL);

   m_pole = calculate_singlet_mass(M_loop);
}
"
                 ];
              ,
              mq = CConversion`RValueToCFormString[GetParticleFromDescription[name]];
              result =
"\
const double p = model_0l.get_M" <> mq <> "(i);
const auto self_energy_1  = Re(model_0l.self_energy_" <> mq <> "_1loop_1(p));
const auto self_energy_PL = Re(model_0l.self_energy_" <> mq <> "_1loop_PL(p));
const auto self_energy_PR = Re(model_0l.self_energy_" <> mq <> "_1loop_PR(p));
const auto M_tree = model_0l.get_mass_matrix_" <> mq <> "();
const auto M_loop = (M_tree - self_energy_PR * M_tree
                     - M_tree * self_energy_PL - self_energy_1).eval();

" <> CConversion`CreateCType[CConversion`ArrayType[CConversion`realScalarCType,3]] <> " M_pole;
fs_svd(M_loop, M_pole);

m_pole = M_pole(i);"
             ];
           result
          ];

CalculateMUpQuarkPole1L[]    := CalculateMFermionPole1L["Up-Quarks"  ,
                                                        TreeMasses`GetSMUpQuarks,
                                                        TreeMasses`GetUpQuark];
CalculateMDownQuarkPole1L[]  := CalculateMFermionPole1L["Down-Quarks",
                                                        TreeMasses`GetSMDownQuarks,
                                                        TreeMasses`GetDownQuark];
CalculateMDownLeptonPole1L[] := CalculateMFermionPole1L["Leptons",
                                                        TreeMasses`GetSMChargedLeptons,
                                                        TreeMasses`GetDownLepton];

FillSMFermionPoleMasses[] :=
    Module[{result = "", i, mq},
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetUpQuark[i + 1, True]];
               result = result <>
                        "this->model.get_physical().M" <> mq <> " = " <>
                        "eft.get_physical().MFu(" <> ToString[i] <> ");\n";
              ];
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetDownQuark[i + 1, True]];
               result = result <>
                        "this->model.get_physical().M" <> mq <> " = " <>
                        "eft.get_physical().MFd(" <> ToString[i] <> ");\n";
              ];
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetDownLepton[i + 1, True]];
               result = result <>
                        "this->model.get_physical().M" <> mq <> " = " <>
                        "eft.get_physical().MFe(" <> ToString[i] <> ");\n";
              ];
           result
          ];

(* find loop order at which parameter needs to be determined given the
   expression which determines the Higgs mass

   par = parameter

   exprs = list of expressions which determine Mh.
   exprs[[1]] contains the tree-level expression.
   exprs[[2]] contains the 1-loop expression.
   exprs[[3]] contains the 2-loop expression.

   towerLoopOrder = requested loop order of the tower
 *)
FindMHiggsLoopOrderFor[par_, exprs_List, towerLoopOrder_] :=
    Module[{i},
           For[i = 0, i < Length[exprs], i++,
               If[!FreeQ[exprs[[i+1]], par],
                  Return[Max[{0, towerLoopOrder - i}]];
                 ];
              ];
           0
          ];

GetFixedBSMParameters[susyScaleMatching_List] :=
    Intersection[Join[{SARAH`hyperchargeCoupling, SARAH`leftCoupling, SARAH`strongCoupling,
                       SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa},
                      First /@ susyScaleMatching],
                      Parameters`GetModelParameters[]
                     ];

SetBSMParameterAtLoopOrder[par_, lo_, struct_String] :=
    Module[{parName = CConversion`ToValidCSymbolString[par]},
           struct <> "set_" <> parName <> "(model_" <> ToString[lo] <> "l.get_" <> parName <> "());\n"
          ];

SetBSMParameters[susyScaleMatching_List, higgsMassMatrix_, struct_String:""] :=
    Module[{pars, loopOrder},
           (* BSM parameters fixed by matching SM -> BSM *)
           pars = GetFixedBSMParameters[susyScaleMatching];
           (* find matching loop orders for parameters for FlexibleEFTHiggs-1L *)
           loopOrder = FindMHiggsLoopOrderFor[#, {higgsMassMatrix, 0}, 1]& /@ pars;
           StringJoin[SetBSMParameterAtLoopOrder[#[[1]], #[[2]], struct]& /@ Utils`Zip[pars, loopOrder]]
          ];

GetSingletVEVInTermsOf[muEff_] :=
  Module[{vs = GetParameterFromDescription["Singlet-VEV"]},
    First[vs /. Solve[muEff == FlexibleSUSY`EffectiveMu, vs]]
  ]; 

(* @TODO: generalize gaugeless limit for any model 
   Set all dimensionless parameter to zero except (yt, yb, ytau, g3) in the MSSM.  
   In other models, need to check which 2-loop contributions to Mh are available.
    *)   
     
SetGaugeLessLimit[struct_String] :=
    Module[{result, 
            g1str = CConversion`ToValidCSymbolString[SARAH`hyperchargeCoupling],
            g2str = CConversion`ToValidCSymbolString[SARAH`leftCoupling],
            lambda = Parameters`GetParameterFromDescription["Singlet-Higgs-Interaction"],
            kappa = Parameters`GetParameterFromDescription["Singlet Self-Interaction"],
            vSexpr = GetSingletVEVInTermsOf[FlexibleSUSY`MuInput]},
    result ="\
      "<> struct <>".set_"<> g1str <>"(gauge_less / " <> ToString[FlexibleSUSY`FSModelName] <> "_info::normalization_"<> g1str <>");
      "<> struct <>".set_"<> g2str <>"(gauge_less / " <> ToString[FlexibleSUSY`FSModelName] <> "_info::normalization_"<> g2str <>");";
      If[lambda =!= Null, 
         result = result <> "\n"<> struct <>".set_"<> CConversion`ToValidCSymbolString[lambda] <>"(gauge_less);";
         result = result <> "\n"<> struct <>".set_TLambdax(Re(INPUTPARAMETER(ALambdaInput)*gauge_less));";
         result = result <> "\n"<> Parameters`CreateLocalConstRefs[vSexpr];
         result = result <> "\n"<> struct <>".set_vS("<>CConversion`RValueToCFormString[vSexpr]  <>");";

      ];
      
      If[kappa =!= Null, 
         result= result <> "\n"<> struct <>".set_"<> CConversion`ToValidCSymbolString[kappa] <>"(gauge_less);";
      ];
    result
    ]



End[];

EndPackage[];
