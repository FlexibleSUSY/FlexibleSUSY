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

Model`Name = "HGTHDMII";
Model`NameLaTeX ="Two Higgs Doublet Model II w/ Higgsinos and Gauginos";
Model`Authors = "F.Staub,A.Voigt";
Model`Date = "2014-11-06";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2],RParity};
RpP = 1;
RpM = -1;

(* Gauge Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False,RpP};
Gauge[[2]]={WB, SU[2], left,        g2,True, RpP};
Gauge[[3]]={G,  SU[3], color,       g3,False,RpP};

(* Matter fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3, RpP};
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1, RpP};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3, RpP};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3, RpP};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1, RpP};

FermionFields[[6]]  = {FG, 1, fG          ,   0  , 1, 8, RpM};
FermionFields[[7]]  = {FW, 1, {{fW0/Sqrt[2],fWp},
                               {fWm,-fW0/Sqrt[2]}}, 0 , 3, 1, RpM};
FermionFields[[8]]  = {FB, 1, fB          ,   0  , 1, 1, RpM};
FermionFields[[9]]  = {Hd, 1, {FHd0, FHdm},  -1/2, 2, 1, RpM};
FermionFields[[10]] = {Hu, 1, {FHup, FHu0},   1/2, 2, 1, RpM};

ScalarFields[[1]]  = {H1, 1, {H1p, H10},  1/2, 2,  1, RpP};
ScalarFields[[2]]  = {H2, 1, {H2p, H20},  1/2, 2,  1, RpP};

(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][Additional] = {
    {LagHC      , { AddHC->True }},
    {LagNoHC    , { AddHC->False}},
    {LagSplit   , { AddHC->True }}
};

LagNoHC = -(M112 conj[H1].H1 + M222 conj[H2].H2 + Lambda1 conj[H1].H1.conj[H1].H1 + \
		Lambda2 conj[H2].H2.conj[H2].H2 + Lambda3 conj[H2].H2.conj[H1].H1 + Lambda4 conj[H2].H1.conj[H1].H2 );

LagHC = -(-M122 conj[H1].H2
          + Lambda5/2 conj[H2].H1.conj[H2].H1
          + Lambda6 conj[H1].H1.conj[H1].H2
          + Lambda7 conj[H2].H2.conj[H1].H2
          + Yd conj[H1].d.q + Ye conj[H1].e.l + Yu H2.u.q
          + Xd conj[H2].d.q + Xe conj[H2].e.l + Xu H1.u.q);

LagSplit = - MassG/2 FG.FG - MassWB/2 FW.FW - MassB/2 FB.FB - \[Mu] Hu.Hd \
    - g2u conj[H2].FW.Hu - (g2up/Sqrt[2]) conj[H2].FB.Hu \
    + g1d H1.FW.Hd - (g1dp/Sqrt[2]) H1.FB.Hd;

DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB , conj[fB] },
    FW0  -> {fW0, 0},
    FWp  -> {fWp, conj[fWm]},
    Glu  -> {fG , conj[fG] },
    FH0  -> {FHd0, conj[FHu0]},
    FHC  -> {FHdm, conj[FHup]},
    Fd1  -> {dL, 0},
    Fd2  -> {0, dR},
    Fu1  -> {uL, 0},
    Fu2  -> {0, uR},
    Fe1  -> {eL, 0},
    Fe2  -> {0, eR},
    Fv   -> {vL,0}
};

(* ----- After EWSB ----- *)

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3]}, {VP,VZ}, ZZ},
    {{VWB[1],VWB[2]}, {VWm,conj[VWm]}, ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
    {H10, {v1, 1/Sqrt[2]}, {sigma1, \[ImaginaryI]/Sqrt[2]}, {phi1, 1/Sqrt[2]}},
    {H20, {v2, 1/Sqrt[2]}, {sigma2, \[ImaginaryI]/Sqrt[2]}, {phi2, 1/Sqrt[2]}}
};

(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector] = {
    {{phi1, phi2}, {hh, ZH}},
    {{sigma1, sigma2}, {Ah, ZA}},
    {{conj[H1p],conj[H2p]},{Hm,ZP}},
    {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
    {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
    {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}},
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
    Fd -> {DL, conj[DR]},
    Fe -> {EL, conj[ER]},
    Fu -> {UL, conj[UR]},
    Fv -> {vL, 0},
    Chi -> {L0, conj[L0]},
    Cha -> {Lm, conj[Lp]},
    Glu -> {fG, conj[fG]}
};

DEFINITION[EWSB][GaugeES]={
    Fd1 -> {FdL, 0},
    Fd2 -> {0, FdR},
    Fu1 -> {Fu1, 0},
    Fu2 -> {0, Fu2},
    Fe1 -> {Fe1, 0},
    Fe2 -> {0, Fe2}
};
