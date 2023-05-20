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

Model`Name = "MSSMNoFV";
Model`NameLaTeX ="MSSM without Flavor violation";
Model`Authors = "F.Staub, A.Voigt";
Model`Date = "2014-04-08";


(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2],RParity};
RpM = {-1,-1,1};
RpP = {1,1,-1};

(* Vector Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False,RpM};
Gauge[[2]]={WB, SU[2], left,        g2,True, RpM};
Gauge[[3]]={G,  SU[3], color,       g3,False,RpM};


(* Chiral Superfields *)

SuperFields[[1]] = {q,  3, {uL0,  dL0},    1/6, 2, 3, RpM};
SuperFields[[2]] = {l,  3, {vL0,  eL0},   -1/2, 2, 1, RpM};
SuperFields[[3]] = {Hd, 1, {Hd0, Hdm},  -1/2, 2, 1, RpP};
SuperFields[[4]] = {Hu, 1, {Hup, Hu0},   1/2, 2, 1, RpP};

SuperFields[[5]] = {d, 3, conj[dR0],  1/3, 1, -3, RpM};
SuperFields[[6]] = {u, 3, conj[uR0], -2/3, 1, -3, RpM};
SuperFields[[7]] = {e, 3, conj[eR0],    1, 1,  1, RpM};


(*------------------------------------------------------*)
(* Superpotential *)
(*------------------------------------------------------*)

SuperPotential = Yu u.q.Hu - Yd d.q.Hd - Ye e.l.Hd + \[Mu] Hu.Hd;


(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};


(* ----- After EWSB ----- *)


(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3]},{VP,VZ},ZZ},
    {{VWB[1],VWB[2]},{VWm,conj[VWm]},ZW},
    {{fWB[1],fWB[2],fWB[3]},{fWm,fWp,fW0},ZfW}
};


(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
    {SHd0, {vd, 1/Sqrt[2]}, {sigmad, \[ImaginaryI]/Sqrt[2]},{phid, 1/Sqrt[2]}},
    {SHu0, {vu, 1/Sqrt[2]}, {sigmau, \[ImaginaryI]/Sqrt[2]},{phiu, 1/Sqrt[2]}}
};


(* ----- Flavors ---- *)

DEFINITION[EWSB][Flavors] = {
    {FdR0,{FdR,FsR,FbR}},
    {FdL0,{FdL,FsL,FbL}},
    {FuL0,{FuL,FcL,FtL}},
    {FuR0,{FuR,FcR,FtR}},
    {FvL0,{FveL,FvmL,FvtL}},
    {SdR0,{SdR,SsR,SbR}},
    {SdL0,{SdL,SsL,SbL}},
    {SuL0,{SuL,ScL,StL}},
    {SuR0,{SuR,ScR,StR}},
    {FeL0,{FeL,FmL,FtauL}},
    {FeR0,{FeR,FmR,FtauR}},
    {SeR0,{SeR,SmR,StauR}},
    {SeL0,{SeL,SmL,StauL}},
    {SvL0,{SveL,SvmL,SvtL}}
};


(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector] = {
    {{SdL, SdR}, {Sd, ZD}},
    {{SuL, SuR}, {Su, ZU}},
    {{SeL, SeR}, {Se, ZE}},
    {{SmL, SmR}, {Sm, ZM}},
    {{StauL, StauR}, {Stau, ZTau}},
    {{SsL, SsR}, {Ss, ZS}},
    {{ScL, ScR}, {Sc, ZC}},
    {{SbL, SbR}, {Sb, ZB}},
    {{StL, StR}, {St, ZT}},
    {{phid, phiu}, {hh, ZH}},
    {{sigmad, sigmau}, {Ah, ZA}},
    {{SHdm,conj[SHup]},{Hpm,ZP}},
    {{fB, fW0, FHd0, FHu0}, {L0, ZN}},
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}}
};

DEFINITION[EWSB][Phases] = {
    {fG, PhaseGlu}
};


(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors] = {
    Fd -> {FdL, FdR},
    Fb -> {FbL, FbR},
    Fs -> {FsL, FsR},
    Fc -> {FcL, FcR},
    Ft -> {FtL, FtR},
    Fm -> {FmL, FmR},
    Ftau -> {FtauL, FtauR},
    Fe -> {FeL, FeR},
    Fu -> {FuL, FuR},
    Fve -> {FveL, 0},
    Fvm -> {FvmL, 0},
    Fvt -> {FvtL, 0},
    Chi -> {L0, conj[L0]},
    Cha -> {Lm, conj[Lp]},
    Glu -> {fG, conj[fG]}
};


DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB, conj[fB]},
    Wino -> {fWB, conj[fWB]},
    H0 -> {FHd0, conj[FHu0]},
    HC -> {FHdm, conj[FHup]},
    Fd1 -> {FdL0, 0},
    Fd2 -> {0, FdR0},
    Fu1 -> {FuL0, 0},
    Fu2 -> {0, FuR0},
    Fe1 -> {FeL0, 0},
    Fe2 -> {0, FeR0},
    Fv -> {FvL0,0},
    Glu -> {fG, conj[fG]}
};
