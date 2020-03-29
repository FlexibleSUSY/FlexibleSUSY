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
    {AlphaS,    { Description -> "Alpha Strong"}},
    {e,         { Description -> "electric charge"}},
    {Gf,        { Description -> "Fermi's constant"}},
    {aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},
    {Yu,        { Description -> "Up-Yukawa-Coupling",
                  DependenceNum ->  Sqrt[2]/v* {{Mass[Fu,1],0,0},
                                                {0, Mass[Fu,2],0},
                                                {0, 0, Mass[Fu,3]}}}},
    {Yd,        { Description -> "Down-Yukawa-Coupling",
                  DependenceNum ->  Sqrt[2]/v* {{Mass[Fd,1],0,0},
                                                {0, Mass[Fd,2],0},
                                                {0, 0, Mass[Fd,3]}}}},
    {Ye,        { Description -> "Lepton-Yukawa-Coupling",
                  DependenceNum ->  Sqrt[2]/v* {{Mass[Fe,1],0,0},
                                                {0, Mass[Fe,2],0},
                                                {0, 0, Mass[Fe,3]}}}},
    {mu2,       { Description -> "SM Mu Parameter",
                  LaTeX -> "m^2",
                  OutputName -> mu2 }},
    {\[Lambda], { Description -> "SM Higgs Selfcouplings",
                  DependenceNum -> Mass[hh]^2/(2 v^2)}},
    {v,         { Description -> "EW-VEV",
                  DependenceNum -> Sqrt[4*Mass[VWp]^2/(g2^2)],
                  DependenceSPheno -> None  }},
    {ThetaW,    { Description -> "Weinberg-Angle",
                  DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]  }},
    {ZZ,        { Description -> "Photon-Z Mixing Matrix"}},
    {ZW,        { Description -> "W Mixing Matrix",
                  Dependence ->   1/Sqrt[2] {{1, 1},
                                             {\[ImaginaryI],-\[ImaginaryI]}} }},
    {Vu,        { Description ->"Left-Up-Mixing-Matrix"}},
    {Vd,        { Description ->"Left-Down-Mixing-Matrix"}},
    {Uu,        { Description ->"Right-Up-Mixing-Matrix"}},
    {Ud,        { Description ->"Right-Down-Mixing-Matrix"}},
    {Ve,        { Description ->"Left-Lepton-Mixing-Matrix"}},
    {Ue,        { Description ->"Right-Lepton-Mixing-Matrix"}},

(* split MSSM parameters *)

    {ZN,        { Description->"Neutralino Mixing-Matrix" }},
    {UP,        { Description->"Chargino-plus Mixing-Matrix"}},
    {UM,        { Description->"Chargino-minus Mixing-Matrix"}},
    {MassB,     { Description -> "Bino Mass parameter",
                  Real -> True }},
    {MassWB,    { Description -> "Wino Mass parameter",
                  Real -> True }},
    {MassG,     { Description -> "Gluino Mass parameter",
                  Real -> True }},
    {PhaseGlu,  { Description -> "Gluino-Phase" }},
    {\[Mu],     { Description -> "Mu-parameter",
                  Real -> True,
                  LaTeX -> "\\mu_{\\text{MSSM}}",
                  OutputName -> Mu }},
    {gYu,       { Description -> "Higgs-Bino-Up-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> gYu,
                  LaTeX -> "\\tilde{g}_{1\\text{u}}",
                  LesHouches -> {SplitMSSM,1} }},
    {g2u,       { Description -> "Higgs-Wino-Up-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> g2u,
                  LaTeX -> "\\tilde{g}_{2\\text{u}}",
                  LesHouches -> {SplitMSSM,2} }},
    {gYd,       { Description -> "Higgs-Bino-Down-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> gYd,
                  LaTeX -> "\\tilde{g}_{1\\text{d}}",
                  LesHouches -> {SplitMSSM,3} }},
    {g2d,       { Description -> "Higgs-Wino-Down-Higgsino-Coupling",
                  Real -> True,
                  OutputName -> g2d,
                  LaTeX -> "\\tilde{g}_{2\\text{d}}",
                  LesHouches -> {SplitMSSM,4} }}
};
