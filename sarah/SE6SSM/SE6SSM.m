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

Off[General::spell]

Model`Name = "SE6SSM";
Model`NameLaTeX = "SE6SSM";
Model`Authors = "P.Athron,D.Harries";
Model`Date = "2013-11-10";

(*-------------------------------------------*)
(*   Particle Content                        *)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2], HParity};
HpM = {-1,-1,1};
HpP = {1,1,-1};

Global[[2]] = {Z[2], EParity};
EpM = -1;
EpP = 1;

(* Vector Superfields *)

Gauge[[1]] = {B,   U[1], hypercharge, g1, False, HpM, EpP};
Gauge[[2]] = {WB, SU[2], left,        g2, True,  HpM, EpP};
Gauge[[3]] = {G,  SU[3], color,       g3, False, HpM, EpP};
Gauge[[4]] = {Bp,  U[1], Ncharge,    g1p, False, HpM, EpP};

(* Chiral Superfields *)

SuperFields[[1]] = {q,  3, {uL,  dL},     1/6, 2, 3, 1 / Sqrt[40], HpM, EpP};
SuperFields[[2]] = {l,  3, {vL,  eL},    -1/2, 2, 1, 2 / Sqrt[40], HpM, EpP};
SuperFields[[3]] = {Hd, 1, {Hd0, Hdm},   -1/2, 2, 1, -3 / Sqrt[40], HpP, EpP};
SuperFields[[4]] = {Hu, 1, {Hup, Hu0},    1/2, 2, 1, -2 / Sqrt[40], HpP, EpP};

SuperFields[[5]] = {d, 3, conj[dR],    1/3, 1, -3, 2 / Sqrt[40], HpM, EpP};
SuperFields[[6]] = {u, 3, conj[uR],   -2/3, 1, -3, 1 / Sqrt[40], HpM, EpP};
SuperFields[[7]] = {e, 3, conj[eR],      1, 1,  1, 1 / Sqrt[40], HpM, EpP};
SuperFields[[8]] = {s, 1, sR,     0, 1,  1, 5 / Sqrt[40], HpP, EpP};
SuperFields[[9]] = {sbar, 1, sbarR,    0, 1,  1, -5 / Sqrt[40], HpP, EpP};

SuperFields[[10]] = {H1I, 2, {H1I0, H1Im},  -1/2, 2, 1, -3 / Sqrt[40], HpM, EpM};
SuperFields[[11]] = {H2I, 2, {H2Ip, H2I0},   1/2, 2, 1, -2 / Sqrt[40], HpM, EpM};
SuperFields[[12]] = {SI, 3, SIR,    0, 1,  1, 5 / Sqrt[40], HpM, EpM};

SuperFields[[13]] = {Dx, 3, DxL,  -1/3, 1, 3, -2 / Sqrt[40], HpM, EpM};
SuperFields[[14]] = {Dxbar, 3, conj[DxbarR],  1/3, 1, -3, -3 / Sqrt[40], HpM, EpM};

SuperFields[[15]] = {L4, 1, {L40, L4m},  -1/2, 2,  1, 2 / Sqrt[40], HpP, EpM};
SuperFields[[16]] = {L4bar, 1, {L4barp, L4bar0}, 1/2, 2,  1, -2 / Sqrt[40], HpP, EpM};
SuperFields[[17]] = {phi, 1, phiR, 0, 1, 1, 0, HpP, EpP};

NoU1Mixing = True;
AddMixedSofts = False;

(*------------------------------------------------------*)
(* Z_2^H Exact Superpotential                           *)
(*------------------------------------------------------*)

SuperPotential = Yu u.q.Hu - Yd d.q.Hd - Ye e.l.Hd + \[Lambda] s.Hu.Hd + \[Lambda]12 s.H2I.H1I + \[Kappa] s.Dx.Dxbar + \[Mu]L L4bar.L4 - \[Sigma] phi.s.sbar + \[Kappa]Pr/3 phi.phi.phi + \[Mu]Phi/2 phi.phi + \[Xi]F phi + fu SI.H1I.Hu + fd SI.Hd.H2I + gD q.L4.Dxbar + hE e.H1I.L4 + \[Sigma]L phi.L4.L4bar;

(*-------------------------------------------*)
(* Integrate Out or Delete Particles         *)
(*-------------------------------------------*)

IntegrateOut = {};
DeleteParticles = {};


(*----------------------------------------------*)
(*   ROTATIONS                                  *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB, conj[fB]},
    Wino -> {fWB, conj[fWB]},
    Glu -> {fG, conj[fG]},
    (*  H0 -> {FHd0, conj[FHu0]},
    HC -> {FHdm, conj[FHup]}, *)
    H01 -> {FHd0, 0},
    H02 -> {0, conj[FHu0]},
    HC1 -> {FHdm, 0},
    HC2 -> {0, conj[FHup]},
    Fd1 -> {FdL, 0},
    Fd2 -> {0, FdR},
    Fu1 -> {FuL, 0},
    Fu2 -> {0, FuR},
    Fe1 -> {FeL, 0},
    Fe2 -> {0, FeR},
    Fv -> {FvL, 0},
    FS1 -> {FsR, 0},
    FS2 -> {0, conj[FsR]},
    FSbar1 -> {FsbarR, 0},
    FSbar2 -> {0, conj[FsbarR]},
    Fphi1 -> {FphiR, 0},
    Fphi2 -> {0, conj[FphiR]},
    FBp -> {fBp, conj[fBp]},
    H0I1 -> {FH1I0, 0},
    H0I2 -> {0, conj[FH2I0]},
    HCI1 -> {FH1Im, 0},
    HCI2 -> {0, conj[FH2Ip]},
    FSI1 -> {FSIR, 0},
    FSI2 -> {0, conj[FSIR]},
    FDx1 -> {FDxL, 0},
    FDx2 -> {0, FDxbarR},
    Hp01 -> {FL40, 0},
    Hp02 -> {0, conj[FL4bar0]}
};

(* ----- After EWSB ----- *)

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB, VWB[3], VBp}, {VP, VZ, VZp}, ZZ},
    {{VWB[1], VWB[2]}, {VWm, conj[VWm]}, ZW},
    {{fWB[1], fWB[2], fWB[3]}, {fWm, fWp, fW0}, ZfW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= {
    {SHd0,   {vd,   1 / Sqrt[2]}, {sigmad,    \[ImaginaryI] / Sqrt[2]}, {phid,    1 / Sqrt[2]}},
    {SHu0,   {vu,   1 / Sqrt[2]}, {sigmau,    \[ImaginaryI] / Sqrt[2]}, {phiu,    1 / Sqrt[2]}},
    {SsR,    {vs,   1 / Sqrt[2]}, {sigmaS,    \[ImaginaryI] / Sqrt[2]}, {phiS,    1 / Sqrt[2]}},
    {SsbarR, {vsb,  1 / Sqrt[2]}, {sigmaSbar, \[ImaginaryI] / Sqrt[2]}, {phiSbar, 1 / Sqrt[2]}},
    {SphiR,  {vphi, 1 / Sqrt[2]}, {sigmaPhi,  \[ImaginaryI] / Sqrt[2]}, {phiPhi,  1 / Sqrt[2]}},
    {SH1I0,  {0, 0},              {sigmaH1I0, \[ImaginaryI] / Sqrt[2]}, {phiH1I0, 1 / Sqrt[2]}},
    {SH2I0,  {0, 0},              {sigmaH2I0, \[ImaginaryI] / Sqrt[2]}, {phiH2I0, 1 / Sqrt[2]}},
    {SSIR,   {0, 0},              {sigmaSIR,  \[ImaginaryI] / Sqrt[2]}, {phiSIR,  1 / Sqrt[2]}}
};

(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector] = {
    {{SdL, SdR}, {Sd, ZD}},
    {{SuL, SuR}, {Su, ZU}},
    {{SeL, SeR}, {Se, ZE}},
    {{SvL}, {Sv, ZV}},
    {{SDxL, SDxbarR}, {SDX, ZDX}},
    {{phid, phiu, phiS, phiSbar, phiPhi}, {hh, ZH}},
    {{sigmad, sigmau, sigmaS, sigmaSbar, sigmaPhi}, {Ah, ZA}},
    {{SHdm, conj[SHup]}, {Hpm, ZP}},
    {{fB, fW0, FHd0, FHu0, FsR, FsbarR, FphiR, fBp}, {L0, ZN}},
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm, UM}, {Lp, UP}}},
    {{{FeL}, {conj[FeR]}}, {{FEL, ZEL}, {FER, ZER}}},
    {{{FdL}, {conj[FdR]}}, {{FDL, ZDL}, {FDR, ZDR}}},
    {{{FuL}, {conj[FuR]}}, {{FUL, ZUL}, {FUR, ZUR}}},
    {{{FDxL}, {conj[FDxbarR]}}, {{FDXL, ZDXL}, {FDXR, ZDXR}}},
    {{phiH1I0, phiH2I0, phiSIR}, {hhI, ZHI}},
    {{sigmaH1I0, sigmaH2I0, sigmaSIR}, {AhI, ZAI}},
    {{SH1Im, conj[SH2Ip]}, {SHIPM, ZPI}},
    {{{FH1Im}, {FH2Ip}}, {{LmI, UMI}, {LpI, UPI}}},
    {{FH1I0, FH2I0, FSIR}, {L0I, ZNI}},
    {{SL40, conj[SL4bar0]}, {SHp0, ZHp}},
    {{SL4m, conj[SL4barp]}, {SHpp, ZPp}},
    {{FL40, FL4bar0}, {L0p, ZNp}}
};

DEFINITION[EWSB][Phases] = {
    {fG, PhaseGlu},
    {FL4m, PhaseChaP}
};

DEFINITION[EWSB][DiracSpinors] = {
    Fd -> {FDL, conj[FDR]},
    Fe -> {FEL, conj[FER]},
    Fu -> {FUL, conj[FUR]},
    Fv -> {FvL, 0},
    Chi -> {L0, conj[L0]},
    Cha -> {Lm, conj[Lp]},
    Glu -> {fG, conj[fG]},
    ChiI -> {L0I, conj[L0I]},
    ChaI -> {LmI, conj[LpI]},
    ChiP -> {L0p, conj[L0p]},
    ChaP -> {FL4m, conj[FL4barp]},
    FDX -> {FDXL, conj[FDXR]}
};
