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
CalculateRunningUpQuarkMasses::usage = "";
CalculateRunningDownQuarkMasses::usage = "";
CalculateRunningDownLeptonMasses::usage = "";
CalculateMUpQuarkPole1L::usage = "";
CalculateMDownQuarkPole1L::usage = "";
CalculateMDownLeptonPole1L::usage = "";
FillSMFermionPoleMasses::usage = "";
GetFixedBSMParameters::usage="Returns a list of the BSM parameters fixed by matching SM -> BSM.";
SetBSMParameters::usage = "";

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

End[];

EndPackage[];
