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

ParameterDefinitions = { 

{g1,        { Description -> "Hypercharge-Coupling"}},
{g2,        { Description -> "Left-Coupling"}},
{g3,        { Description -> "Strong-Coupling"}},    
{AlphaS,    {Description -> "Alpha Strong"}},	
{e,         { Description -> "electric charge"}}, 

{Gf,        { Description -> "Fermi's constant"}},
{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},

{Yu,        { Description -> "Up-Yukawa-Coupling",
			 DependenceNum ->  Sqrt[2]/v2* {{Mass[Fu,1],0,0},
             									{0, Mass[Fu,2],0},
             									{0, 0, Mass[Fu,3]}}}}, 
             									
{Yd,        { Description -> "Down-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/v1* {{Mass[Fd,1],0,0},
             									{0, Mass[Fd,2],0},
             									{0, 0, Mass[Fd,3]}}}},
             									
{Ye,        { Description -> "Lepton-Yukawa-Coupling",
			  DependenceNum ->  Sqrt[2]/v1* {{Mass[Fe,1],0,0},
             									{0, Mass[Fe,2],0},
             									{0, 0, Mass[Fe,3]}}}}, 
{Xu,        { Description -> "wrong Up-Yukawa-Coupling",
              LaTeX -> "X_u", OutputName -> Xu,
              LesHouches -> Xu }},
{Xd,        { Description -> "wrong Down-Yukawa-Coupling",
              LaTeX -> "X_d", OutputName -> Xd,
              LesHouches -> Xd }},
{Xe,        { Description -> "wrong Lepton-Yukawa-Coupling",
              LaTeX -> "X_e", OutputName -> Xe,
              LesHouches -> Xe }},
                                                                            
                                                                           
{Lambda1,    { LaTeX -> "\\lambda_1",
               OutputName -> Lam1,
               LesHouches -> {HMIX,31}}},
{Lambda2,    { LaTeX -> "\\lambda_2",
               OutputName -> Lam2,
               LesHouches -> {HMIX,32}}},
{Lambda3,    { LaTeX -> "\\lambda_3",
               OutputName -> Lam3,
               LesHouches -> {HMIX,33}}},
{Lambda4,    { LaTeX -> "\\lambda_4",
               OutputName -> Lam4,
               LesHouches -> {HMIX,34}}},
{Lambda5,    { LaTeX -> "\\lambda_5",
               OutputName -> Lam5,
               LesHouches -> {HMIX,35}}},

{Lambda6,    { LaTeX -> "\\lambda_6",
               OutputName -> Lam6,
               LesHouches -> {HMIX,36}}},

{Lambda7,    { LaTeX -> "\\lambda_7",
               OutputName -> Lam7,
               LesHouches -> {HMIX,37}}},


{M112,    {    LaTeX -> "m^2_1",
               OutputName -> M112,
               LesHouches -> {HMIX,20}}},


{M222,    {    LaTeX -> "m^2_2",
               OutputName -> M222,
               LesHouches -> {HMIX,21}}},

{M122,    {    LaTeX -> "m^2_{12}",
               OutputName -> M122,
               LesHouches -> {HMIX,22}}},


{v1,        { Description -> "Down-VEV", LaTeX -> "v_1"}}, 
{v2,        { Description -> "Up-VEV", LaTeX -> "v_2"}},       
{v,         { Description -> "EW-VEV", DependenceSPheno -> None }},
             
{\[Beta],   { Description -> "Pseudo Scalar mixing angle"  }},             
{TanBeta,   { Description -> "Tan Beta" }},              
{\[Alpha],  { Description -> "Scalar mixing angle" }},  

{ZH,        { Description->"Scalar-Mixing-Matrix"}},
{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix"}},
{ZP,        { Description->"Charged-Mixing-Matrix"}},  


{ThetaW,    { Description -> "Weinberg-Angle"}}, 

{ZZ, {Description ->   "Photon-Z Mixing Matrix"}},
{ZW, {Description -> "W Mixing Matrix" }},


{Vu,        {Description ->"Left-Up-Mixing-Matrix"}},
{Vd,        {Description ->"Left-Down-Mixing-Matrix"}},
{Uu,        {Description ->"Right-Up-Mixing-Matrix"}},
{Ud,        {Description ->"Right-Down-Mixing-Matrix"}}, 
{Ve,        {Description ->"Left-Lepton-Mixing-Matrix"}},
{Ue,        {Description ->"Right-Lepton-Mixing-Matrix"}},

(* new parameters for THDM-II w/ extra Higgsions and Gauginos *)

{\[Mu],     { Description -> "Mu-parameter",
              Real -> True,
              LaTeX -> "\\mu_{\\text{MSSM}}",
              OutputName -> Mu }},

{ZN,        { Description->"Neutralino Mixing-Matrix" }},
{UP,        { Description->"Chargino-plus Mixing-Matrix"}},
{UM,        { Description->"Chargino-minus Mixing-Matrix"}},

{MassB,     { Description -> "Bino Mass parameter" }},
{MassWB,    { Description -> "Wino Mass parameter" }},
{MassG,     { Description -> "Gluino Mass parameter" }},
{PhaseGlu,  { Description -> "Gluino-Phase" }},
{g2up,      { Description -> "Higgs-Bino-Up-Higgsino-Coupling",
              Real -> True,
              OutputName -> g2up,
              LaTeX -> "\\tilde{g}_{2\\text{u}}^{\\prime}",
              LesHouches -> {SplitMSSM,1} }},
{g2u,       { Description -> "Up-Higgs-Wino-Up-Higgsino-Coupling",
              Real -> True,
              OutputName -> g2u,
              LaTeX -> "\\tilde{g}_{2\\text{u}}",
              LesHouches -> {SplitMSSM,2} }},
{g1dp,      { Description -> "Down-Higgs-Bino-Down-Higgsino-Coupling",
              Real -> True,
              OutputName -> g1dp,
              LaTeX -> "\\tilde{g}_{1\\text{d}}^{\\prime}",
              LesHouches -> {SplitMSSM,3} }},
{g1d,       { Description -> "Down-Higgs-Wino-Down-Higgsino-Coupling",
              Real -> True,
              OutputName -> g1d,
              LaTeX -> "\\tilde{g}_{1\\text{d}}",
              LesHouches -> {SplitMSSM,4} }},
{d2up,      { Description -> "wrong Higgs-Bino-Up-Higgsino-Coupling",
              Real -> True,
              OutputName -> d2up,
              LaTeX -> "\\tilde{\\gamma}_{2\\text{u}}^{\\prime}",
              LesHouches -> {SplitMSSM,5} }},
{d2u,       { Description -> "wrong Up-Higgs-Wino-Up-Higgsino-Coupling",
              Real -> True,
              OutputName -> d2u,
              LaTeX -> "\\tilde{\\gamma}_{2\\text{u}}",
              LesHouches -> {SplitMSSM,6} }},
{d1dp,      { Description -> "wrong Down-Higgs-Bino-Down-Higgsino-Coupling",
              Real -> True,
              OutputName -> d1dp,
              LaTeX -> "\\tilde{\\gamma}_{1\\text{d}}^{\\prime}",
              LesHouches -> {SplitMSSM,7} }},
{d1d,       { Description -> "wrong Down-Higgs-Wino-Down-Higgsino-Coupling",
              Real -> True,
              OutputName -> d1d,
              LaTeX -> "\\tilde{\\gamma}_{1\\text{d}}",
              LesHouches -> {SplitMSSM,8} }}
};
