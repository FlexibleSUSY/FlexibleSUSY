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

CalculateMHiggs2LoopShift::usage = "Calculates 2-loop shift to the Higgs pole mass in the BSM model.";
CalculateMHiggs3LoopShift::usage = "Calculates 3-loop shift to the Higgs pole mass in the BSM model.";
Create2LoopMatching::usage = "Creates function body to calculate SM parameters from BSM parameters at 2-loop level.";
Create3LoopMatching::usage = "Creates function body to calculate SM parameters from BSM parameters at 3-loop level.";
CallMatch2LoopTopMass::usage = "Sets SM top Yukawa coupling to 2-loop value, determined from BSM model";
CreateSMMtop2LoopFunction::usage = "Creates a function that calculates the running MS-bar top quark mass in the SM from the BSM parameters";
CalculateRunningUpQuarkMasses::usage = "";
CalculateRunningDownQuarkMasses::usage = "";
CalculateRunningDownLeptonMasses::usage = "";
CalculateMUpQuarkPole1L::usage = "";
CalculateMDownQuarkPole1L::usage = "";
CalculateMDownLeptonPole1L::usage = "";
FillSMFermionPoleMasses::usage = "Set SM fermion pole masses from BSM fermion pole masses";
GetFixedBSMParameters::usage="Returns a list of the BSM parameters fixed by matching SM -> BSM.";
SetBSMParameters::usage = "";
SetLimit::usage = "applies a given limit to given model";

Begin["`Private`"];

CalculateMHiggs2LoopShift[inputModel_String, higgsBoson_, higgsIndex_String] :=
   Module[{modelNameStr = ToString[FlexibleSUSY`FSModelName],
	   higgsBosonStr = CConversion`RValueToCFormString[higgsBoson],
	   higgsMassStr = CConversion`RValueToCFormString[FlexibleSUSY`M[higgsBoson]]},
"const auto mh2_tree = calculate_mh2_0l(model);
const double mh2_0l = Sqr(model.get_" <> higgsMassStr <> "(" <> higgsIndex <> "));
const double p = 0.;

// calculate 2-loop self-energy using tree-level parameters in the
// gauge-less limit (g1 = g2 = 0)
const auto self_energy_2l = Re(model.self_energy_" <> higgsBosonStr <> "_2loop());

// calculate 2-loop contribution from momentum iteration using
// tree-level parameters in the gauge-less limit (g1 = g2 = 0) and
// Yukawa-less limit
decltype(model.self_energy_" <> higgsBosonStr <> "_1loop_deriv_p2(p)) self_energy_2l_mom;
self_energy_2l_mom.setZero();

if (mh2_1l_bsm_shift != 0.) {
   self_energy_2l_mom += Re(model.self_energy_" <> higgsBosonStr <> "_1loop_deriv_p2(p) * mh2_1l_bsm_shift);
}

const auto self_energy = (self_energy_2l + self_energy_2l_mom).eval();
const auto tadpole = calculate_tadpole_over_vevs(model, 2);
const auto mass_matrix = compose_loop_corrections(mh2_tree, self_energy, tadpole);

const double Mh2_pole = calculate_eigenvalue(mass_matrix, idx);

const double mh2_shift_2l = Mh2_pole - mh2_0l;

return mh2_shift_2l;"
];


CalculateMHiggs3LoopShift[inputModel_String, outputModel_String, higgsBoson_, higgsIndex_String] :=
   Module[{modelNameStr = ToString[FlexibleSUSY`FSModelName],
	   higgsMassStr = CConversion`RValueToCFormString[FlexibleSUSY`M[higgsBoson]]},
"double mh2_shift_3l = std::numeric_limits<double>::quiet_NaN();

try {
   // 3-loop self-energy, calculated using tree-level parameters
   const auto self_energy_3l = Re(model." <> SelfEnergies`CreateSelfEnergyFunctionName[higgsBoson, 3] <> "());
   const auto mass_matrix = (calculate_mh2_0l(model) - self_energy_3l).eval();
   const double Mh2_pole = calculate_eigenvalue(mass_matrix, idx);
   const double mh2_tree = Sqr(model.get_" <> higgsMassStr <> "(idx));

   // calculate 3-loop Higgs mass loop correction in the gauge-less limit
   mh2_shift_3l = Mh2_pole - mh2_tree;
} catch (const flexiblesusy::Error& e) {
   VERBOSE_MSG(\"Error: Calculation of 3-loop Higgs pole mass in the gauge-less limit in the " <> modelNameStr <> " at the matching scale failed: \" << e.what());
   " <> outputModel <> ".get_problems().flag_bad_mass(standard_model_info::hh);
}

return mh2_shift_3l;"
   ];


Create2LoopMatching[inputModel_String, outputModel_String, higgsBoson_, higgsIndex_String] :=
   Module[{modelNameStr = ToString[FlexibleSUSY`FSModelName],
	   higgsMassStr = CConversion`RValueToCFormString[FlexibleSUSY`M[higgsBoson]]},
"// calculate running parameters of BSM model
const auto model = [] (const auto& model_input) {
   auto model = model_input;
   model.calculate_DRbar_masses();
   return model;
}(model_input);
const auto model_gl = make_gaugeless_g1_g2(model_input);
const auto model_gl_no_g3 = make_gaugeless_g3(model_gl);
const auto model_gl_yl = make_yukawaless(model_gl);

const auto sm_0l = match_high_to_low_scale_sm_0l_copy(sm, model, idx);
const auto sm_1l = match_high_to_low_scale_sm_1l_copy(sm, model, idx);
const auto sm_0l_gl = match_high_to_low_scale_sm_0l_copy(sm, model_gl, idx);
const auto sm_1l_gl = match_high_to_low_scale_sm_1l_copy(sm, model_gl, idx);
const auto sm_0l_gl_no_g3 = match_high_to_low_scale_sm_0l_copy(sm, model_gl_no_g3, idx);
const auto sm_1l_gl_no_g3 = match_high_to_low_scale_sm_1l_copy(sm, model_gl_no_g3, idx);
const auto sm_0l_gl_yl = make_yukawaless(sm_0l_gl);

double delta_yt   = 0;
double delta_yb   = 0;
double delta_ytau = 0;
double delta_v_gl = 0;
double p          = 0;

if (model.HIGGS_2LOOP_CORRECTION_AT_AS) {
   delta_yt += sm_1l_gl.get_Yu(2,2) - sm_1l_gl_no_g3.get_Yu(2,2);
}

if (model.HIGGS_2LOOP_CORRECTION_AB_AS) {
   delta_yb += sm_1l_gl.get_Yd(2,2) - sm_1l_gl_no_g3.get_Yd(2,2);
}

if (model.HIGGS_2LOOP_CORRECTION_AT_AT) {
   delta_yt   += sm_1l_gl_no_g3.get_Yu(2,2) - sm_0l_gl.get_Yu(2,2);
   delta_yb   += sm_1l_gl_no_g3.get_Yd(2,2) - sm_0l_gl.get_Yd(2,2);
   delta_ytau += sm_1l.get_Ye(2,2) - sm_0l_gl.get_Ye(2,2);
   delta_v_gl += sm_1l.get_v() - sm_0l.get_v();
}

const double Q = sm_0l_gl.get_scale();
const double mt = sm_0l_gl.get_MFu(2);
const double mb = sm_0l_gl.get_MFd(2);
const double mtau = sm_0l_gl.get_MFe(2);
const double yt = sm_0l_gl.get_Yu(2,2);
const double yb = sm_0l_gl.get_Yd(2,2);
const double ytau = sm_0l_gl.get_Ye(2,2);
const double v = sm_0l.get_v();
const double v2 = Sqr(v);
const double lambda_1l = sm_1l.get_Lambdax();

// 2-loop parameter conversion terms from the sum in Eq.(4.26b)
// [arxiv:2003.04639] for P = {yt, yb, ytau, g3, v}. The terms in
// the sum that involve a derivative w.r.t. p^2 are not included
// here.
const double mh2_2l_parameter_conversion =
   + delta_yt * sm_twoloophiggs::delta_mh_1loop_at_sm_deriv_yt(p, Q, mt, yt)
   + delta_yb * sm_twoloophiggs::delta_mh_1loop_ab_sm_deriv_yb(p, Q, mb, yb)
   + delta_ytau * sm_twoloophiggs::delta_mh_1loop_atau_sm_deriv_ytau(p, Q, mtau, ytau)
   + delta_v_gl * (+ sm_twoloophiggs::delta_mh_1loop_at_sm_deriv_v(p, Q, mt, yt)
                   + sm_twoloophiggs::delta_mh_1loop_ab_sm_deriv_v(p, Q, mb, yb)
                   + sm_twoloophiggs::delta_mh_1loop_atau_sm_deriv_v(p, Q, mtau, ytau));

// Eq.(4.21) [arxiv:2003.04639] for n=1
const double mh2_1l_bsm_shift = calculate_Mh2_1l_shift(model_gl_yl, idx);
// Eq.(4.26a) [arxiv:2003.04639]
const double mh2_1l_sm_shift  = calculate_Mh2_1l_shift(sm_0l_gl_yl);
// Eq.(4.21) [arxiv:2003.04639] for n=2, including derivative term w.r.t. p^2 (momentum iteration)
const double mh2_2l_bsm_shift = calculate_Mh2_2l_shift(model_gl_yl, idx, model.HIGGS_2LOOP_CORRECTION_AT_AT ? mh2_1l_bsm_shift : 0.);
// Eq.(4.26b) [arxiv:2003.04639], including derivative term w.r.t. p^2 (momentum iteration)
const double mh2_2l_sm_shift  = calculate_Mh2_2l_shift(sm_0l_gl_yl, model.HIGGS_2LOOP_CORRECTION_AT_AT ? mh2_1l_bsm_shift : 0.);

// Eq.(4.28c) [arxiv:2003.04639]
const double delta_lambda_2l =
   (mh2_2l_bsm_shift - mh2_2l_sm_shift - mh2_2l_parameter_conversion)/v2
   + 2*delta_v_gl/(v*v2)*(mh2_1l_sm_shift - mh2_1l_bsm_shift);

// Eq.(4.28a) [arxiv:2003.04639] up to (including) 2-loop terms
const double lambda_2l = lambda_1l + delta_lambda_2l;

sm = sm_1l;
sm.set_Lambdax(lambda_2l);

if (yt_loop_order > 1) {
" <> TextFormatting`IndentText[If[FlexibleSUSY`UseHiggs3LoopMSSM === True || FlexibleSUSY`UseHiggs3LoopNMSSM === True, FlexibleEFTHiggsMatching`CallMatch2LoopTopMass[], ""]] <> "
}

sm.calculate_DRbar_masses();"
];

AttachBrackets[FlexibleSUSY`M[sym_[idx_]]] := FlexibleSUSY`M[sym][idx];
AttachBrackets[FlexibleSUSY`M[sym_]] := FlexibleSUSY`M[sym][];
AttachBrackets[sym_] := sym;

Create3LoopMatching[inputModel_String, outputModel_String, higgsBoson_, higgsIndex_String] :=
   Module[{modelNameStr = ToString[FlexibleSUSY`FSModelName],
	   higgsMassStr = CConversion`RValueToCFormString[FlexibleSUSY`M[higgsBoson]],
           mtstr = CConversion`RValueToCFormString[AttachBrackets[TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]]]]},
"// calculate running masses of the input model
const auto model = [] (const auto& model_input) {
   auto model = model_input;
   model.calculate_DRbar_masses();
   return model;
}(" <> inputModel <> ");
const auto model_gl = make_gaugeless_g1_g2(model);
const auto model_gl_no_g3 = make_gaugeless_g3(model_gl);

const auto sm_0l = match_high_to_low_scale_sm_0l_copy(" <> outputModel <> ", model, " <> higgsIndex <> ");
const auto sm_1l = match_high_to_low_scale_sm_1l_copy(" <> outputModel <> ", model, " <> higgsIndex <> ");
const auto sm_0l_gl = match_high_to_low_scale_sm_0l_copy(" <> outputModel <> ", model_gl, " <> higgsIndex <> ");
const auto sm_1l_gl = match_high_to_low_scale_sm_1l_copy(" <> outputModel <> ", model_gl, " <> higgsIndex <> ");
const auto sm_1l_gl_no_g3 = match_high_to_low_scale_sm_1l_copy(" <> outputModel <> ", model_gl_no_g3, " <> higgsIndex <> ");
const auto sm_2l = match_high_to_low_scale_sm_2l_copy(" <> outputModel <> ", model, " <> higgsIndex <> ", 2);

" <> outputModel <> " = sm_2l;

// calculation of 3-loop threshold corrections below

const double lambda_2l = sm_2l.get_Lambdax();

const double v2 = Sqr(sm_0l_gl.get_v());
const double yt = sm_0l_gl.get_Yu(2, 2);
const double yt2 = Sqr(yt);
const double g3 = sm_0l_gl.get_g3();
const double g32 = Sqr(g3);
const double g34 = Sqr(g32);
const double Q = sm_0l_gl.get_scale();
const double Q2 = Sqr(Q);
const double mt = model_gl.get_" <> mtstr <> ";
const double mt2 = Sqr(mt);
const double logmt = Log(mt2 / Q2);
const double logmt2 = Sqr(logmt);
const double logmt3 = logmt * logmt2;
const double k = oneOver16PiSqr;
const double k2 = twoLoop;
const double k3 = threeLoop;

const double delta_yt_1l = sm_1l_gl.get_Yu(2, 2) - sm_1l_gl_no_g3.get_Yu(2, 2);
const double delta_yt_2l = sm_2l.get_Yu(2, 2) - sm_1l.get_Yu(2, 2);
const double delta_g3_1l = sm_1l.get_g3() - g3;

// 1st derivative of 1-loop SM contribution to Mh w.r.t. yt, times Delta yt(2l)
const double S1_deriv_yt =
   delta_yt_2l * sm_twoloophiggs::delta_mh_1loop_at_sm_deriv_yt(0, Q, sm_0l_gl.get_MFu(2), yt);
// 2nd derivative of 1-loop SM contribution to Mh w.r.t. yt, times [Delta yt(1l)]^2
const double S1_deriv_yt2 =
   -0.5*(24*k*mt2*(7 + 6*logmt)) * Sqr(delta_yt_1l);
// 1st derivative of 2-loop SM contribution to Mh w.r.t. yt, times Delta yt(1l)
const double S2_deriv_yt =
   64*k2*g32*mt2*yt*(1 + 8*logmt + 6*logmt2) * delta_yt_1l;
// 1st derivative of 2-loop SM contribution to Mh w.r.t. g3, times Delta g3(1l)
const double S2_deriv_g3 =
   64*k2*g3*mt2*yt2*logmt*(1 + 3*logmt) * delta_g3_1l;

// 3-loop parameter conversion terms from thesum on the r.h.s.
// of Eq.(4.26c) [arxiv:2003.04639]
const double mh2_3l_parameter_conversion = S1_deriv_yt + S1_deriv_yt2 + S2_deriv_yt + S2_deriv_g3;

// Eq.(4.26c) [arxiv:2003.04639] w/o the sum
const double mh2_3l_sm_shift = -Re(sm_0l_gl.self_energy_hh_3loop());

// Eq.(4.21) [arxiv:2003.04639] for n=3
const double mh2_3l_bsm_shift = calculate_Mh2_3l_shift(sm, model_gl, idx);

// Eq.(4.28d) [arxiv:2003.04639], Delta lambda^(3l)
const double delta_lambda_3l = std::isfinite(mh2_3l_bsm_shift) ? (mh2_3l_bsm_shift - mh2_3l_sm_shift - mh2_3l_parameter_conversion)/v2 : 0.0;

// Eq.(4.28a) [arxiv:2003.04639] up to (including) 3-loop terms
const double lambda_3l = lambda_2l + delta_lambda_3l;

" <> outputModel <> ".set_Lambdax(lambda_3l);
" <> outputModel <> ".calculate_DRbar_masses();"
];

CallMatch2LoopTopMass[smStr_:"sm."] :=
smStr <> "set_Yu(2, 2, calculate_yt_sm_2l(sm_0l, sm_1l, model));"

CreateSMMtop2LoopFunction[] :=
Module[{g3str = ToString[TreeMasses`GetStrongCoupling[]], 
        mglustr = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]],
        md2str = CConversion`RValueToCFormString[SARAH`SoftDown],
        mq2str = CConversion`RValueToCFormString[SARAH`SoftSquark],
        mtstr = CConversion`RValueToCFormString[TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]]]},
"/**
 * Calculates SM top quark masses at tree-level and 2-loop level.
 * Returns the tuple (mt, Mt) with mt being the SM tree-level top
 * quark mass and Mt being the SM top quark pole mass calculated
 * at 2-loop level (including 2-loop QCD contributions).
 *
 * @param sm Standard Model parameters
 * @return tuple (mt, Mt)
 */
std::tuple<double, double> calculate_mt_sm_2l(
   const standard_model::Standard_model& sm)
{
   auto sm_2l = sm;
   auto loop_corrections = sm_2l.get_loop_corrections();
   loop_corrections.top_qcd = 1; // enable 2-loop QCD contributions

   sm_2l.set_pole_mass_loop_order(2);
   sm_2l.set_loop_corrections(loop_corrections);

   sm_2l.calculate_MFu_pole();

   const auto Mt = sm_2l.get_physical().MFu(2);
   const auto mt = sm_2l.get_MFu(2);

   return std::make_tuple(mt, Mt);
}

/**
 * Calculates stop masses and stop mixing angle.
 *
 * @param model BSM model parameters
 * @return tuple (mst_1, mst_2, theta_t)
 */
std::tuple<double, double, double> calculate_mstop(
   const " <> ToString[FlexibleSUSY`FSModelName] <> "_mass_eigenstates& model)
{
   double mst_1, mst_2, theta_t;
   model." <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <> ";

   return std::make_tuple(mst_1, mst_2, theta_t);
}

/**
 * Calculates top quark pole in the BSM model at 2-loop level, including
 * 2-loop QCD contributions.
 *
 * @param model BSM model parameters
 * @return 2-loop top quark mass in the BSM
 */
double calculate_Mt_bsm_2l(
   const " <> ToString[FlexibleSUSY`FSModelName] <> "_mass_eigenstates& model)
{
   auto model_2l = model;
   auto loop_corrections = model_2l.get_loop_corrections();
   loop_corrections.top_qcd = 1; // enable 2-loop QCD contributions

   model_2l.set_pole_mass_loop_order(2);
   model_2l.set_loop_corrections(loop_corrections);

   model_2l." <> CreateLoopMassFunctionName[Parameters`StripIndices[TreeMasses`GetUpQuark[3,True]]] <> "();

   const auto Mt_bsm = model_2l.get_physical()." <> mtstr <> ";

   return Mt_bsm;
}

/**
 * Calculates running MS-bar SM top quark mass, given the BSM model parameters.
 *
 * @note The given SM and BSM model parameters should be in the
 * gauge-less limit (g1 = g2 = 0).
 *
 * @param sm Standard Model parameters
 * @param model BSM model parameters
 * @param idx Higgs index (from the mass ordered Higgs multiplet in the high-scale model)
 */
double calculate_MFt_MSbar_sm_2l(
   const standard_model::Standard_model& sm,
   const " <> ToString[FlexibleSUSY`FSModelName] <> "_mass_eigenstates& model)
{
   const double Q = model.get_scale();
   const double Q2 = Sqr(Q);
   const double mt = sm.get_MFu(2); // is equal to the top quark mass in the MSSM
   const double mt2 = Sqr(mt);
   const double g3 = model.get_" <> g3str <> "();
   const double g32 = Sqr(g3);
   const double alpha_s = g32*oneOver4Pi;

   const double k = oneOver16PiSqr;
   const double logmt = Log(mt2 / Q2);
   const double yt_SM = sm.get_Yu(2, 2);

   const auto [mst_1, mst_2, theta_t] = calculate_mstop(model);

   const mssm_twoloop_mt::Parameters pars{
      .g3 = g3,
      .mt = mt,
      .mg = model.get_" <> mglustr <> "(),
      .mst1 = mst_1,
      .mst2 = mst_2,
      .msusy = Sqrt(Sqrt(Abs(model.get_" <> md2str <> "(2, 2)))) *
               Sqrt(Sqrt(Abs(model.get_" <> mq2str <> "(2, 2)))),
      .xt = Sin(2*theta_t)*(Sqr(mst_1) - Sqr(mst_2))/(2*mt),
      .Q = Q
   };

   // top mass contribution at O(as) in MSSM and SM
   const double sqcd_1l_over_mt = mssm_twoloop_mt::dMt_over_mt_1loop(pars);
   const double qcd_1l_SM_over_mt = k * g32 * 4./3 * (4 - 3*Log(mt2 / Q2));

   // top mass contribution at O(as) in the SM
   const double qcd1l_mt = -4./3 * k * g32 * mt * (2 + 3*logmt);
   const double qcd1l_g3 =  8./3 * k * g32 * mt * (4 - 3*logmt);

   const double delta_alpha_s = calculate_delta_alpha_s(alpha_s, model);
   const double delta_mt_over_mt = sqcd_1l_over_mt - qcd_1l_SM_over_mt;
   const double delta_g3_over_g3 = -0.5 * delta_alpha_s;

   // 2-loop contribution from parameter conversion
   const double mt_con = delta_mt_over_mt * qcd1l_mt + delta_g3_over_g3 * qcd1l_g3;

   // calculate top quark pole masses in the BSM and in the SM
   const double Mt_bsm = calculate_Mt_bsm_2l(model);
   const auto [mt_sm_0l, Mt_sm] = calculate_mt_sm_2l(sm);

   // determine SM MS-bar top quark mass from top quark pole mass matching,
   // Eq.(4.18b) [arxiv:2003.04639]
   const double mt_sm = Mt_bsm - Mt_sm + mt_sm_0l - mt_con;

   return Abs(mt_sm);
}

/**
 * Calculates SM top quark Yukawa coupling at 2-loop level.
 *
 * @param sm_0l SM parameters (determined from BSM parameters at tree-level)
 * @param sm_1l SM parameters (determined from BSM parameters at 1-loop level)
 * @param model BSM parameters
 */
double calculate_yt_sm_2l(
   const standard_model::Standard_model& sm_0l,
   const standard_model::Standard_model& sm_1l,
   const " <> ToString[FlexibleSUSY`FSModelName] <> "_mass_eigenstates& model)
{
   const double delta_vev = sm_1l.get_v() - sm_0l.get_v();
   const double yt_0l = sm_0l.get_Yu(2,2);
   const double mt_2l = calculate_MFt_MSbar_sm_2l(sm_0l, model);
   const double v     = sm_0l.get_v();
   const double yt_2l = sqrt2 * mt_2l / v - yt_0l * delta_vev / v;

   return yt_2l;
}"
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
"if (i == " <> iStr <> ") {
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
"const double p = model_0l.get_M" <> mq <> "(i);
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

FillSMFermionPoleMasses[eftStr_:"eft.", bsmStr_:"this->model."] :=
    Module[{result = "", i, mq},
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetUpQuark[i + 1, True]];
               result = result <>
                        bsmStr <> "get_physical().M" <> mq <> " = " <>
                        eftStr <> "get_physical().MFu(" <> ToString[i] <> ");\n";
              ];
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetDownQuark[i + 1, True]];
               result = result <>
                        bsmStr <> "get_physical().M" <> mq <> " = " <>
                        eftStr <> "get_physical().MFd(" <> ToString[i] <> ");\n";
              ];
           For[i = 0, i < 3, i++,
               mq = CConversion`RValueToCFormString[TreeMasses`GetDownLepton[i + 1, True]];
               result = result <>
                        bsmStr <> "get_physical().M" <> mq <> " = " <>
                        eftStr <> "get_physical().MFe(" <> ToString[i] <> ");\n";
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

(* set `parameter' to `value' for given model `struct' *)
FEFTApplyParameterSetting[{parameter_, value_}, struct_String] :=
    Which[Parameters`IsModelParameter[parameter],
          Parameters`SetParameter[parameter, value, struct],
          Parameters`IsInputParameter[parameter],
          Parameters`SetInputParameter[parameter, value, "INPUTPARAMETER"],
          Parameters`IsPhase[parameter],
          Parameters`SetPhase[parameter, value, struct],
          Parameters`IsExtraParameter[parameter],
          Parameters`SetParameter[parameter, value, struct],
          True,
          Print["Error: ", parameter, " cannot be set for the MSSM-limit,",
                " because it is neither a model nor an input parameter!"];
          Quit[1];
          ""
         ];

(* @todo: Set all dimensionless parameters to zero except (yt, yb, ytau, g3) in the MSSM.
   @todo: Check which 2-loop contributions to Mh are available.
 *)
SetLimit[struct_String, limit_List] :=
   Switch[limit,
	  {}, "",
	  {{_,_}..},
	  Parameters`CreateLocalConstRefs[(#[[2]])& /@ limit] <> "\n" <>
	  StringJoin[FEFTApplyParameterSetting[#, struct]& /@ limit],
	  _,
	  Print["Error: The list is malformed (neither empty nor a list of 2-component lists): ", limit];
	  Quit[1];
	  ""
   ];

End[];

EndPackage[];
