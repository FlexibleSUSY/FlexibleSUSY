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
    {VB,  { Description -> "B-Boson"}},
    {VG,  { Description -> "Gluon"}},
    {VWB, { Description -> "W-Bosons"}},
    {gB,  { Description -> "B-Boson Ghost"}},
    {gG,  { Description -> "Gluon Ghost" }},
    {gWB, { Description -> "W-Boson Ghost"}},
    {Fd1, { Description -> "Dirac Left Down-Quark"}},
    {Fd2, { Description -> "Dirac Right Down-Quark"}},
    {Fu1, { Description -> "Dirac Left Up-Quark"}},
    {Fu2, { Description -> "Dirac Right Up-Quark"}},
    {Fe1, { Description -> "Dirac Left Electron"}},
    {Fe2, { Description -> "Dirac Right Electron"}},
    {Fv,  { Description -> "Neutrinos" }},

    (* particles for THDM-II w/ Higgsinos and Gauginos *)

    {Glu,  { Description -> "Gluino"}},
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
             OutputName -> "FWp" }}
};


ParticleDefinitions[EWSB] = {
    {hh,   { Description -> "Higgs"}},
    {Ah,   { Description -> "Pseudo-Scalar Higgs"}},
    {Hm,   { Description -> "Charged Higgs"}},
    {VP,   { Description -> "Photon"}},
    {VZ,   { Description -> "Z-Boson"}},
    {VG,   { Description -> "Gluon" }},
    {VWm,  { Description -> "W-Boson",
             Goldstone -> Hm[{1}] }},
    {gP,   { Description -> "Photon Ghost"}},
    {gWm,  { Description -> "Negative W-Boson Ghost"}},
    {gWmC, { Description -> "Positive W-Boson Ghost" }},
    {gZ,   { Description -> "Z-Boson Ghost" }},
    {gG,   { Description -> "Gluon Ghost" }},
    {Fd,   { Description -> "Down-Quarks"}},
    {Fu,   { Description -> "Up-Quarks"}},
    {Fe,   { Description -> "Leptons" }},
    {Fv,   { Description -> "Neutrinos" }},

    (* particles for THDM-II w/ Higgsinos *)

    {Glu,  { Description -> "Gluino" }},
    {Chi,  { Description -> "Neutralinos"}},
    {Cha,  { Description -> "Charginos"}}
};


WeylFermionAndIndermediate = {
   {H,      {PDG -> {0},
             Width -> 0,
             Mass -> Automatic,
             LaTeX -> "H",
             OutputName -> "" }},
   {H10,    {LaTeX -> "H_1^0"}},
   {H20,    {LaTeX -> "H_2^0"}},
   {H1p,    {LaTeX -> "H_1^+"}},
   {H2p,    {LaTeX -> "H_2^+"}},
   {sigma1, {LaTeX -> "\\sigma_1"}},
   {sigma2, {LaTeX -> "\\sigma_2"}},
   {phi1,   {LaTeX -> "\\phi_1"}},
   {phi2,   {LaTeX -> "\\phi_2"}},
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

     {FG,     { Description -> "Gluino field",
                PDG -> {0},
                Width -> 0,
                Mass -> Automatic,
                LaTeX -> "\\tilde{g}",
                OutputName -> "" }},
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
     {fG,     { Description ->"Gluino Weyl-Spinor"}},
     {fW0,    { Description ->"Neutral Wino" }},
     {fWm,    { Description ->"Negative Wino"}},
     {fWp,    { Description ->"Positive Wino"}},
     {fB,     { Description ->"Bino Weyl-Spinor"}}
};
