(* ::Package:: *)




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
    {H0,  { PDG -> {0},
            Width -> 0,
            Mass -> Automatic,
            FeynArtsNr -> 1,
            LaTeX -> "H^0",
            OutputName -> "H0" }},

    {Hp,  { PDG -> {0},
            Width -> 0,
            Mass -> Automatic,
            FeynArtsNr -> 2,
            LaTeX -> "H^+",
            OutputName -> "Hp" }},

    {VB,   { Description -> "B-Boson"}},
    {VG,   { Description -> "Gluon"}},
    {VWB,  { Description -> "W-Bosons"}},
    {gB,   { Description -> "B-Boson Ghost"}},
    {gG,   { Description -> "Gluon Ghost" }},
    {gWB,  { Description -> "W-Boson Ghost"}},
    {Fd1,  { Description -> "Dirac Left Down-Quark"}},
    {Fd2,  { Description -> "Dirac Right Down-Quark"}},
    {Fu1,  { Description -> "Dirac Left Up-Quark"}},
    {Fu2,  { Description -> "Dirac Right Up-Quark"}},
    {Fe1,  { Description -> "Dirac Left Electron"}},
    {Fe2,  { Description -> "Dirac Right Electron"}},
    {Fv,   { Description -> "Neutrinos" }},

    (* split MSSM particles *)

    {Bino, { Description -> "Bino"}},
    {FH0,  { Description -> "Neutral Higgsinos",
             OutputName -> "FH0"}},
    {FHC,  { Description -> "Charged Higgsinos",
             OutputName -> "FHC"}},
    {FW0,  { Description -> "Neutral Wino",
             PDG -> {0},
             PDG.IX -> {0},
             Width -> 0,
             Mass -> Automatic,
             FeynArtsNr -> 666,
             LaTeX -> "\\tilde{W}^0",
             OutputName -> "FW0" }},
    {FWp,  { Description -> "Positive Wino",
             PDG -> {0},
             PDG.IX -> {0},
             Width -> 0,
             Mass -> Automatic,
             FeynArtsNr -> 667,
             LaTeX -> "\\tilde{W}^+",
             OutputName -> "FWp" }},
    {SeL1,  { Description -> "Left Selectron",
             PDG -> {0},
             PDG.IX -> {0},
             Mass -> Automatic,
             FeynArtsNr -> 13,
             LaTeX -> "\\tilde{e}_L",
             OutputName -> "eL" }},
    {SeL2,  { Description -> "Left Smuon",
             PDG -> {0},
             PDG.IX -> {0},
             Mass -> Automatic,
             FeynArtsNr -> 14,
             LaTeX -> "\\tilde{\\mu}_L",
             OutputName -> "muL" }},
    {SeR1,  { Description -> "Right Selectron",
             PDG -> {0},
             PDG.IX -> {0},
             Mass -> Automatic,
             FeynArtsNr -> 15,
             LaTeX -> "\\tilde{e}_R",
             OutputName -> "eR" }},
    {SeR2,  { Description -> "Right Smuon",
             PDG -> {0},
             PDG.IX -> {0},
             Mass -> Automatic,
             FeynArtsNr -> 16,
             LaTeX -> "\\tilde{\\mu}_R",
             OutputName -> "muR" }},
    {SvL1,  { Description -> "Left electron-sneutrino",
             PDG -> {0},
             PDG.IX -> {0},
             Mass -> Automatic,
             FeynArtsNr -> 11,
             LaTeX -> "\\tilde{\\nu}_{eL}",
             OutputName -> "enL" }},
    {SvL2,  { Description -> "Left muon-sneutrino",
             PDG -> {0},
             PDG.IX -> {0},
             Mass -> Automatic,
             FeynArtsNr -> 12,
             LaTeX -> "\\tilde{\\nu}_{\\muL}",
             OutputName -> "munL" }}
};

ParticleDefinitions[EWSB] = {
    {hh,   {  Description -> "Higgs",
              PDG -> {25},
              PDG.IX -> {101000001} }},

    {Ah,   {  Description -> "Pseudo-Scalar Higgs",
              PDG -> {0},
              PDG.IX ->{0},
              Mass -> {0},
              Width -> {0} }},

    {Hp,   { Description -> "Charged Higgs",
             PDG -> {0},
             PDG.IX ->{0},
             Width -> {0},
             Mass -> {0},
             LaTeX -> {"H^+","H^-"},
             OutputName -> {"Hp","Hm"} }},

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


    {Fd,   { Description -> "Down-Quarks" }},
    {Fu,   { Description -> "Up-Quarks" }},
    {Fe,   { Description -> "Leptons" }},
    {Fv,   { Description -> "Neutrinos" }},

    (* split MSSM particles *)

    {Chi,  { Description -> "Neutralinos" }},
    {Cha,  { Description -> "Charginos" }},
    {Se,   { Description -> "Sleptons",
             PDG ->  {1000011, 1000013, 2000011, 2000013},
             PDG.IX -> {-200000601, -200000602, -200000604, -200000605},
             Mass -> LesHouches,
             FeynArtsNr -> 12,
             ElectricCharge -> -1,
             LaTeX -> "\\tilde{e}",
             LHPC -> {5, "green"},
             OutputName -> "se" }},
    {Sv,   { Description -> "Sneutrinos",
             PDG ->  {1000012, 1000014},
             PDG.IX -> {200000001, 200000002},
             Mass -> LesHouches,
             FeynArtsNr -> 11,
             ElectricCharge -> 0,
             LaTeX -> "\\tilde{\\nu}",
             LHPC -> {5, "turquoise"},
             OutputName -> "sv" }}
};

WeylFermionAndIndermediate = {
     {H,      { PDG -> {0},
                Width -> 0,
                Mass -> Automatic,
                LaTeX -> "H",
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
     {UL,     {LaTeX -> "U_L" }},

    (* split MSSM particles *)

     {FW,     { Description -> "Wino field",
                PDG -> {0},
                Width -> 0,
                Mass -> Automatic,
                LaTeX -> "\\tilde{W}",
                OutputName -> "" }},
     {FB,     { Description -> "Bino field",
                PDG -> {0},
                Width -> 0,
                Mass -> Automatic,
                LaTeX -> "\\tilde{B}",
                OutputName -> "" }},
     {Hd,     { Description -> "Down-Higgs Superfield",
                PDG -> {0},
                Width -> 0,
                Mass -> Automatic,
                LaTeX -> "H_d",
                OutputName -> "" }},
     {Hu,     { Description -> "Up-Higgs Superfield",
                PDG -> {0},
                Width -> 0,
                Mass -> Automatic,
                LaTeX -> "H_u",
                OutputName -> "" }},

     {FHd0,   { Description -> "Neutral Down-Higgsino"}},
     {FHu0,   { Description -> "Neutral Up-Higgsino" }},
     {FHdm,   { Description -> "Charged Down-Higgsino"}},
     {FHup,   { Description -> "Charged Up-Higgsino"}},
     {L0,     { Description -> "Neutralino Weyl-Spinor"}},
     {Lm,     { Description -> "Negative Chargino Weyl-Spinor"}},
     {Lp,     { Description -> "Positive Chargino Weyl-Spinor"}},
     {fW0,    { Description ->"Neutral Wino" }},
     {fWm,    { Description ->"Negative Wino"}},
     {fWp,    { Description ->"Positive Wino"}},
     {fB,     { Description ->"Bino Weyl-Spinor"}},
     {SEL1,   { Description -> "Left 1st generation Slepton doublet" }},
     {SEL2,   { Description -> "Left 2nd generation Slepton doublet" }},
     {SER1,   { Description -> "Right 1st generation Slepton doublet" }},
     {SER2,   { Description -> "Right 2nd generation Slepton doublet" }}
};
