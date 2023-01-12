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

BeginPackage["GM2Calc`", {"Parameters`", "CConversion`", "TextFormatting`"}];

(* GM2Calc interface parameters (gauge basis) *)
{ yukawaType, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6,
  lambda7, tanBeta, m122, zetau, zetad, zetal, deltau, deltad, deltal,
  piu, pid, pil };

FillGM2CalcInterfaceData::usage="Fills GM2Calc interface struct with model parameters.";

Begin["Private`"];

(* fill struct with model parameters to be passed to GM2Calc *)
FillGM2CalcInterfaceData[struct_String] :=
    Which[
        IsGM2CalcCompatibleMSSM[], FillGM2CalcMSSMNoFVInterfaceData[struct],
        IsGM2CalcCompatibleTHDM[], FillGM2CalcTHDMInterfaceData[struct, FlexibleSUSY`FSGM2CalcInput],
        True, Print["Error: This model is neither a MSSM-like nor a 2HDM-like model compatible with GM2Calc."]; Quit[1]
    ];

(* returns true, if model is an MSSM-like model compatible with GM2Calc *)
IsGM2CalcCompatibleMSSM[] :=
    Module[{w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
            mu, m1, m2, m3, mq2, mu2, md2, ml2, me2, tu, td, te, yu, yd, ye},
           w             = Parameters`GetParticleFromDescription["W-Boson"];
           pseudoscalar  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           smuon         = Parameters`GetParticleFromDescription["Smuon"];
           muonsneutrino = Parameters`GetParticleFromDescription["Muon Sneutrino"];
           chargino      = Parameters`GetParticleFromDescription["Charginos"];
           neutralino    = Parameters`GetParticleFromDescription["Neutralinos"];
           mu            = Parameters`GetParameterFromDescription["Mu-parameter"];
           m1            = Parameters`GetParameterFromDescription["Bino Mass parameter"];
           m2            = Parameters`GetParameterFromDescription["Wino Mass parameter"];
           m3            = Parameters`GetParameterFromDescription["Gluino Mass parameter"];
           mq2           = Parameters`GetParameterFromDescription["Softbreaking left Squark Mass"];
           mu2           = Parameters`GetParameterFromDescription["Softbreaking right Up-Squark Mass"];
           md2           = Parameters`GetParameterFromDescription["Softbreaking right Down-Squark Mass"];
           ml2           = Parameters`GetParameterFromDescription["Softbreaking left Slepton Mass"];
           me2           = Parameters`GetParameterFromDescription["Softbreaking right Slepton Mass"];
           tu            = Parameters`GetParameterFromDescription["Trilinear-Up-Coupling"];
           td            = Parameters`GetParameterFromDescription["Trilinear-Down-Coupling"];
           te            = Parameters`GetParameterFromDescription["Trilinear-Lepton-Coupling"];
           yu            = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd            = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye            = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           Not[MemberQ[{w, pseudoscalar, smuon, muonsneutrino,
                        chargino, neutralino, mu, m1, m2, m3, mq2, mu2,
                        md2, ml2, me2, tu, td, te, yu, yd, ye}, Null]]
    ];

(* fill struct with MSSM parameters to be passed to GM2Calc *)
FillGM2CalcMSSMNoFVInterfaceData[struct_String] :=
    Module[{mwStr, w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
            mu, m1, m2, m3, mq2, mu2, md2, ml2, me2, tu, td, te, yu, yd, ye},
           w             = Parameters`GetParticleFromDescription["W-Boson"];
           pseudoscalar  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           smuon         = Parameters`GetParticleFromDescription["Smuon"];
           muonsneutrino = Parameters`GetParticleFromDescription["Muon Sneutrino"];
           chargino      = Parameters`GetParticleFromDescription["Charginos"];
           neutralino    = Parameters`GetParticleFromDescription["Neutralinos"];
           mu            = Parameters`GetParameterFromDescription["Mu-parameter"];
           m1            = Parameters`GetParameterFromDescription["Bino Mass parameter"];
           m2            = Parameters`GetParameterFromDescription["Wino Mass parameter"];
           m3            = Parameters`GetParameterFromDescription["Gluino Mass parameter"];
           mq2           = Parameters`GetParameterFromDescription["Softbreaking left Squark Mass"];
           mu2           = Parameters`GetParameterFromDescription["Softbreaking right Up-Squark Mass"];
           md2           = Parameters`GetParameterFromDescription["Softbreaking right Down-Squark Mass"];
           ml2           = Parameters`GetParameterFromDescription["Softbreaking left Slepton Mass"];
           me2           = Parameters`GetParameterFromDescription["Softbreaking right Slepton Mass"];
           tu            = Parameters`GetParameterFromDescription["Trilinear-Up-Coupling"];
           td            = Parameters`GetParameterFromDescription["Trilinear-Down-Coupling"];
           te            = Parameters`GetParameterFromDescription["Trilinear-Lepton-Coupling"];
           yu            = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd            = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye            = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           mwStr         = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           (* compose string *)
           "#ifdef ENABLE_GM2CALC\n" <>
           "GM2Calc_MSSMNoFV_data " <> struct <> ";\n" <>
           struct <> ".scale = MODEL.get_scale();\n" <>
           struct <> ".alpha_em_MZ = ALPHA_EM_MZ;\n" <>
           struct <> ".alpha_em_0 = ALPHA_EM_0;\n" <>
           struct <> ".alpha_s_MZ = ALPHA_S_MZ;\n" <>
           struct <> ".MZ    = MZPole;\n" <>
           "if (!is_zero(" <> mwStr <> ")) {\n" <>
              TextFormatting`IndentText[struct <> ".MW = " <> mwStr <> ";"] <> "\n" <>
           "} else if (!is_zero(MWPole)) {\n" <>
              TextFormatting`IndentText[struct <> ".MW = MWPole;"] <> "\n}\n" <>
           struct <> ".mb_mb = MBMB;\n" <>
           struct <> ".MT    = MTPole;\n" <>
           struct <> ".MTau  = MTauPole;\n" <>
           struct <> ".MM    = MMPole;\n" <>
           struct <> ".MA0   = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[pseudoscalar][1]] <> ";\n" <>
           struct <> ".MSvm  = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[muonsneutrino]] <> ";\n" <>
           struct <> ".TB    = MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> "() / " <>
                              "MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> "();\n" <>
           struct <> ".Mu    = MODEL.get_" <> CConversion`RValueToCFormString[mu] <> "();\n" <>
           struct <> ".M1    = MODEL.get_" <> CConversion`RValueToCFormString[m1] <> "();\n" <>
           struct <> ".M2    = MODEL.get_" <> CConversion`RValueToCFormString[m2] <> "();\n" <>
           struct <> ".M3    = MODEL.get_" <> CConversion`RValueToCFormString[m3] <> "();\n" <>
           struct <> ".MSm   = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[smuon]] <> ";\n" <>
           struct <> ".MCha  = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[chargino]] <> ";\n" <>
           struct <> ".MChi  = MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[neutralino]] <> ";\n" <>
           struct <> ".mq2   = MODEL.get_" <> CConversion`RValueToCFormString[mq2] <> "();\n" <>
           struct <> ".mu2   = MODEL.get_" <> CConversion`RValueToCFormString[mu2] <> "();\n" <>
           struct <> ".md2   = MODEL.get_" <> CConversion`RValueToCFormString[md2] <> "();\n" <>
           struct <> ".ml2   = MODEL.get_" <> CConversion`RValueToCFormString[ml2] <> "();\n" <>
           struct <> ".me2   = MODEL.get_" <> CConversion`RValueToCFormString[me2] <> "();\n" <>
           struct <> ".Au    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[tu] <> "(), MODEL.get_" <> CConversion`RValueToCFormString[yu] <> "());\n" <>
           struct <> ".Ad    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[td] <> "(), MODEL.get_" <> CConversion`RValueToCFormString[yd] <> "());\n" <>
           struct <> ".Ae    = div_safe(MODEL.get_" <> CConversion`RValueToCFormString[te] <> "(), MODEL.get_" <> CConversion`RValueToCFormString[ye] <> "());\n" <>
           "#endif\n\n"
          ];

(* returns true, if model is an THDM-like model compatible with GM2Calc *)
IsGM2CalcCompatibleTHDM[] :=
    Module[{w, h, a, hp, yu, yd, ye},
           w  = Parameters`GetParticleFromDescription["W-Boson"];
           h  = Parameters`GetParticleFromDescription["Higgs"];
           a  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           hp = Parameters`GetParticleFromDescription["Charged Higgs"];
           yu = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           Not[MemberQ[{w, h, a, hp, yu, yd, ye}, Null]]
    ];

(* fill struct with THDM parameters to be passed to GM2Calc *)
FillGM2CalcTHDMInterfaceData[struct_String, inputPars_List] :=
    Module[{w, mwStr, higgs, mhStr,
            inPars = Parameters`DecreaseIndexLiterals[inputPars],
            pars = First /@ Reverse /@ inputPars
           },
           w     = Parameters`GetParticleFromDescription["W-Boson"];
           higgs = Parameters`GetParticleFromDescription["Higgs"];
           mwStr = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           mhStr = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[higgs][0]];
           (* compose string *)
           "#ifdef ENABLE_GM2CALC\n" <>
           "GM2Calc_THDM_data " <> struct <> ";\n" <>
           "{\n" <>
              TextFormatting`IndentText[Parameters`CreateLocalConstRefs[pars]] <> "\n" <>
              TextFormatting`IndentText[
                 struct <> ".alpha_em_mz = ALPHA_EM_MZ;\n" <>
                 struct <> ".alpha_em_0 = ALPHA_EM_0;\n" <>
                 struct <> ".alpha_s_mz = ALPHA_S_MZ;\n" <>
                 "if (!is_zero(" <> mhStr <> ")) {\n" <>
                    TextFormatting`IndentText[struct <> ".mh = " <> mhStr <> ";"] <> "\n" <>
                 "} else if (!is_zero(MHPole)) {\n" <>
                    TextFormatting`IndentText[struct <> ".mh = MHPole;"] <> "\n}\n" <>
                 "if (!is_zero(" <> mwStr <> ")) {\n" <>
                    TextFormatting`IndentText[struct <> ".mw = " <> mwStr <> ";"] <> "\n" <>
                 "} else if (!is_zero(MWPole)) {\n" <>
                    TextFormatting`IndentText[struct <> ".mw = MWPole;"] <> "\n}\n" <>
                 struct <> ".mz = MZPole;\n" <>
                 struct <> ".mu(0) = MU2GeV;\n" <>
                 struct <> ".mu(1) = MCMC;\n" <>
                 struct <> ".mu(2) = MTPole;\n" <>
                 struct <> ".md(0) = MD2GeV;\n" <>
                 struct <> ".md(1) = MS2GeV;\n" <>
                 struct <> ".md(2) = MBMB;\n" <>
                 struct <> ".mv(0) = Mv1Pole;\n" <>
                 struct <> ".mv(1) = Mv2Pole;\n" <>
                 struct <> ".mv(2) = Mv3Pole;\n" <>
                 struct <> ".ml(0) = MEPole;\n" <>
                 struct <> ".ml(1) = MMPole;\n" <>
                 struct <> ".ml(2) = MTauPole;\n" <>
                 struct <> ".ckm = CKMInput;\n" <>
                 struct <> ".yukawa_type = " <> CConversion`RValueToCFormString[FlexibleSUSY`yukawaType /. inPars] <> ";\n" <>
                 struct <> ".lambda(0) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda1 /. inPars] <> ";\n" <>
                 struct <> ".lambda(1) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda2 /. inPars] <> ";\n" <>
                 struct <> ".lambda(2) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda3 /. inPars] <> ";\n" <>
                 struct <> ".lambda(3) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda4 /. inPars] <> ";\n" <>
                 struct <> ".lambda(4) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda5 /. inPars] <> ";\n" <>
                 struct <> ".lambda(5) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda6 /. inPars] <> ";\n" <>
                 struct <> ".lambda(6) = " <> CConversion`RValueToCFormString[FlexibleSUSY`lambda7 /. inPars] <> ";\n" <>
                 struct <> ".tan_beta = " <> CConversion`RValueToCFormString[FlexibleSUSY`tanBeta /. inPars] <> ";\n" <>
                 struct <> ".m122 = " <> CConversion`RValueToCFormString[FlexibleSUSY`m122 /. inPars] <> ";\n" <>
                 struct <> ".zeta_u = " <> CConversion`RValueToCFormString[FlexibleSUSY`zetau /. inPars] <> ";\n" <>
                 struct <> ".zeta_d = " <> CConversion`RValueToCFormString[FlexibleSUSY`zetad /. inPars] <> ";\n" <>
                 struct <> ".zeta_l = " <> CConversion`RValueToCFormString[FlexibleSUSY`zetal /. inPars] <> ";\n" <>
                 struct <> ".delta_u = " <> CConversion`RValueToCFormString[FlexibleSUSY`deltau /. inPars] <> ";\n" <>
                 struct <> ".delta_d = " <> CConversion`RValueToCFormString[FlexibleSUSY`deltad /. inPars] <> ";\n" <>
                 struct <> ".delta_l = " <> CConversion`RValueToCFormString[FlexibleSUSY`deltal /. inPars] <> ";\n" <>
                 struct <> ".pi_u = " <> CConversion`RValueToCFormString[FlexibleSUSY`piu /. inPars] <> ";\n" <>
                 struct <> ".pi_d = " <> CConversion`RValueToCFormString[FlexibleSUSY`pid /. inPars] <> ";\n" <>
                 struct <> ".pi_l = " <> CConversion`RValueToCFormString[FlexibleSUSY`pil /. inPars] <> ";\n"
              ] <>
           "}\n" <>
           "#endif\n\n"
    ];

End[];

EndPackage[];
