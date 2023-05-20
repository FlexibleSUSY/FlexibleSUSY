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
    {g1,        { Description -> "Hypercharge-Coupling",
                  GUTnormalization -> 1}},
    {g2,        { Description -> "Left-Coupling"}},
    {g3,        { Description -> "Strong-Coupling"}},
    {AlphaS,    { Description -> "Alpha Strong"}},
    {e,         { Description -> "electric charge"}},
    {Gf,        { Description -> "Fermi's constant"}},
    {aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},
    {Yu,        { Description -> "Up-Yukawa-Coupling"}},
    {Yd,        { Description -> "Down-Yukawa-Coupling"}},
    {Ye,        { Description -> "Lepton-Yukawa-Coupling"}},
    {mH2,       { Description -> "SM Mu Parameter"}},
    {mS2,       { Description -> "Singlet mass parameter"}},
    {\[Lambda]1,{ Description -> "SM Higgs Selfcouplings"}},
    {\[Lambda]2,{ Description -> "Singlet 4-selfcoupling"}},
    {\[Lambda]3,{ Description -> "Singlet 2S-2H coupling"}},
    {\[Lambda]4,{ Description -> "Singlet S-2H coupling"}},
    {\[Lambda]5,{ Description -> "Singlet 3-selfcoupling"}},
    {vH,        { Description -> "EW-VEV",
                  DependenceSPheno -> None }},
    {vS,        { Description -> "Singlet-VEV",
                  DependenceSPheno -> None  }},
    {ThetaW,    { Description -> "Weinberg-Angle",
                  DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]  }},
    {ZZ,        { Description -> "Photon-Z Mixing Matrix"}},
    {ZH,        { Description -> "Higgs Mixing Matrix",
                  LesHouches -> "ZH" }},
    {ZW,        { Description -> "W Mixing Matrix",
                  Dependence ->   1/Sqrt[2] {{1, 1},
                                             {\[ImaginaryI],-\[ImaginaryI]}} }},
    {Vu,        { Description ->"Left-Up-Mixing-Matrix"}},
    {Vd,        { Description ->"Left-Down-Mixing-Matrix"}},
    {Uu,        { Description ->"Right-Up-Mixing-Matrix"}},
    {Ud,        { Description ->"Right-Down-Mixing-Matrix"}},
    {Ve,        { Description ->"Left-Lepton-Mixing-Matrix"}},
    {Ue,        { Description ->"Right-Lepton-Mixing-Matrix"}}
};
