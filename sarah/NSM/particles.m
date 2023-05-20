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

ParticleDefinitions[GaugeES] = {
    {H0,  {    PDG -> {0},
               Width -> 0,
               Mass -> Automatic,
               FeynArtsNr -> 1,
               LaTeX -> "H^0",
               OutputName -> "H0" }},

    {Hp,  {    PDG -> {0},
               Width -> 0,
               Mass -> Automatic,
               FeynArtsNr -> 2,
               LaTeX -> "H^+",
               OutputName -> "Hp" }},

    {s0,  {    PDG -> {0},
               Width -> 0,
               Mass -> Automatic,
               FeynArtsNr -> 3,
               LaTeX -> "S^0",
               OutputName -> "S0" }},

    {VB,   { Description -> "B-Boson"}},
    {VG,   { Description -> "Gluon"}},
    {VWB,  { Description -> "W-Bosons"}},
    {gB,   { Description -> "B-Boson Ghost"}},
    {gG,   { Description -> "Gluon Ghost" }},
    {gWB,  { Description -> "W-Boson Ghost"}}

};

ParticleDefinitions[EWSB] = {
    {hh   ,  {  Description -> "Higgs",
                PDG -> {25,35},
                PDG.IX -> {101000001, 101000002} }},

    {Ah   ,  {  Description -> "Pseudo-Scalar Higgs",
                PDG -> {0},
                PDG.IX ->{0},
                Mass -> {0},
                Width -> {0} }},

    {Hp,     { Description -> "Charged Higgs",
               PDG -> {0},
               PDG.IX ->{0},
               Width -> {0},
               Mass -> {0},
               LaTeX -> {"H^+","H^-"},
               OutputName -> {"Hp","Hm"}
             }},

    {VP,   { Description -> "Photon"}},
    {VZ,   { Description -> "Z-Boson",
             Goldstone -> Ah }},
    {VG,   { Description -> "Gluon" }},
    {VWp,  { Description -> "W+ - Boson",
             Goldstone -> Hp }},
    {gP,   { Description -> "Photon Ghost"}},
    {gWp,  { Description -> "Positive W+ - Boson Ghost"}},
    {gWpC, { Description -> "Negative W+ - Boson Ghost" }},
    {gZ,   { Description -> "Z-Boson Ghost" }},
    {gG,   { Description -> "Gluon Ghost" }},

    {Fd,   { Description -> "Down-Quarks"}},
    {Fu,   { Description -> "Up-Quarks"}},
    {Fe,   { Description -> "Leptons" }},
    {Fv,   { Description -> "Neutrinos" }}
};

WeylFermionAndIndermediate = {
    {H,      { PDG -> {0},
               Width -> 0,
               Mass -> Automatic,
               LaTeX -> "H",
               OutputName -> "" }},
    {s,      { PDG -> {0},
               Width -> 0,
               Mass -> Automatic,
               LaTeX -> "S",
               OutputName -> "" }},

    {dR,     {LaTeX -> "d_R" }},
    {eR,     {LaTeX -> "e_R" }},
    {lep,    {LaTeX -> "l" }},
    {uR,     {LaTeX -> "u_R" }},
    {q,      {LaTeX -> "q" }},
    {eL,     {LaTeX -> "e_L" }},
    {dL,     {LaTeX -> "d_L" }},
    {uL,     {LaTeX -> "u_L" }},
    {vL,     {LaTeX -> "\\nu_L" }},

    {DR,     {LaTeX -> "D_R" }},
    {ER,     {LaTeX -> "E_R" }},
    {UR,     {LaTeX -> "U_R" }},
    {EL,     {LaTeX -> "E_L" }},
    {DL,     {LaTeX -> "D_L" }},
    {UL,     {LaTeX -> "U_L" }}
};
