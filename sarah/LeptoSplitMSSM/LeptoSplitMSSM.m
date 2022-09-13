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

Off[General::spell]

Model`Name = "LeptoSplitMSSM";
Model`NameLaTeX ="LeptoSplitMSSM";
Model`Authors = "F.Esser,A.Voigt";
Model`Date = "2018-12-04";

(* Implementation of the split-MSSM with extra light 1st and 2nd
   generation sleptons.  The convention of arXiv:1407.4081 is used. *)

(*-------------------------------------------*)
(*   Parameters in matrix form               *)
(*-------------------------------------------*)

me2 = {{me11, me12}, {me21, me22}};
ml2 = {{ml11, ml12}, {ml21, ml22}};
glllla = {{glllla11, glllla12}, {glllla12, glllla22}};
gllhha = {{gllhha11, gllhha12}, {gllhha21, gllhha22}};
gllhhb = {{gllhhb11, gllhhb12}, {gllhhb21, gllhhb22}};
geeee = {{geeee11, geeee12}, {geeee12, geeee22}};
geehh = {{geehh11, geehh12}, {geehh21, geehh22}};
glleea = {{glleea11, glleea12}, {glleea21, glleea22}};
gleh =  {{gleh11, gleh12}, {gleh21, gleh22}};


(*-------------------------------------------*)
(*   Particle Content                        *)
(*-------------------------------------------*)

(* Global symmetries *)

Global[[1]] = {Z[2],RParity};
RpP = 1;
RpM = -1;

(* Vector Superfields *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False,RpP};
Gauge[[2]]={WB, SU[2], left,        g2,True, RpP};
Gauge[[3]]={G,  SU[3], color,       g3,False,RpP};

(* Matter fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3, RpP};
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1, RpP};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3, RpP};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3, RpP};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1, RpP};

FermionFields[[7]]  = {FW, 1, {{fW0/Sqrt[2],fWp},
                               {fWm,-fW0/Sqrt[2]}}, 0 , 3, 1, RpM};
FermionFields[[8]]  = {FB, 1, fB          ,   0  , 1, 1, RpM};
FermionFields[[9]]  = {Hd, 1, {FHd0, FHdm},  -1/2, 2, 1, RpM};
FermionFields[[10]] = {Hu, 1, {FHup, FHu0},   1/2, 2, 1, RpM};

ScalarFields[[1]]  = {H  , 1, {Hp, H0}  ,   1/2, 2,  1, RpP};
ScalarFields[[2]]  = {SEL1, 1, {SvL1,SeL1}, -1/2, 2,  1, RpM};
ScalarFields[[3]]  = {SEL2, 1, {SvL2,SeL2}, -1/2, 2,  1, RpM};
ScalarFields[[4]]  = {SER1, 1, conj[SeR1],    1, 1,  1, RpM};
ScalarFields[[5]]  = {SER2, 1, conj[SeR2],    1, 1, 1, RpM};

(*----------------------------------------------*)
(*   ROTATIONS                                  *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] = {
    {LagHC    , {AddHC->True }},
    {LagNoHC  , {AddHC->False}},
    {LagSplit , {AddHC->True }},
    {LagReSlep, {AddHC->False}},
    {LagSlep  , {AddHC->True }}
};

LagNoHC = mu2 conj[H].H - 1/2 \[Lambda] conj[H].H.conj[H].H;

LagHC = - Yd conj[H].d.q - Ye conj[H].e.l + Yu H.u.q;

(* arXiv:1407.4081, Eq. (4), w/o gluino *)
LagSplit = \
    - MassWB/2 FW.FW - MassB/2 FB.FB - \[Mu] Hu.Hd \
    - g2u conj[H].FW.Hu - (g1u/Sqrt[2]) conj[H].FB.Hu \
    + g2d H.FW.Hd - (g1d/Sqrt[2]) H.FB.Hd;

(* The following terms are in principle also allowed, but omitted in
   this model *)
LagSplit2 = \
    - g2u Hu.conj[FW].conj[H] \
    + g2d FW.conj[H].conj[Hd] \
    - (g1u/Sqrt[2]) Hu.conj[FB].conj[H] \
    - (g1d/Sqrt[2]) FB.conj[H].conj[Hd] ;

(* real part of slepton Lagrangian *)
LagReSlep = (
    - ml2[[1,1]] conj[SEL1].SEL1 - ml2[[1,2]] conj[SEL1].SEL2 - ml2[[2,1]] conj[SEL2].SEL1 - ml2[[2,2]] conj[SEL2].SEL2
    - me2[[1,1]] conj[SER1].SER1 - me2[[1,2]] conj[SER1].SER2 - me2[[2,1]] conj[SER2].SER1 - me2[[2,2]] conj[SER2].SER2
    (* using completeness relation of pauli matrices, assuming real Yukawa couplings *)
    + glllla[[1,1]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL1].SEL1.conj[SEL1].SEL1
    + 2 glllla[[1,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL1].SEL1.conj[SEL2].SEL2
    + glllla[[2,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[SEL2].SEL2.conj[SEL2].SEL2
    + 2 gllllb Delta[lef1,lef4] Delta[lef2,lef3]conj[SEL1].SEL1.conj[SEL2].SEL2
    + gllhha[[1,1]] Delta[lef1,lef2] Delta[lef3,lef4] conj[H].H.SEL1.conj[SEL1]
    + gllhha[[1,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[H].H.SEL1.conj[SEL2]
    + gllhha[[2,1]] Delta[lef1,lef2] Delta[lef3,lef4] conj[H].H.SEL2.conj[SEL1]
    + gllhha[[2,2]] Delta[lef1,lef2] Delta[lef3,lef4] conj[H].H.SEL2.conj[SEL2]
    + gllhhb[[1,1]] Delta[lef1,lef4] Delta[lef2,lef3] conj[H].H.SEL1.conj[SEL1]
    + gllhhb[[1,2]] Delta[lef1,lef4] Delta[lef2,lef3] conj[H].H.SEL1.conj[SEL2]
    + gllhhb[[2,1]] Delta[lef1,lef4] Delta[lef2,lef3] conj[H].H.SEL2.conj[SEL1]
    + gllhhb[[2,2]] Delta[lef1,lef4] Delta[lef2,lef3] conj[H].H.SEL2.conj[SEL2]
    + geeee[[1,1]] conj[SER1].SER1.SER1.conj[SER1]
    + 2 geeee[[1,2]] conj[SER1].SER1.SER2.conj[SER2]
    + geeee[[2,2]] conj[SER2].SER2.SER2.conj[SER2]
    + geehh[[1,1]] conj[H].H.conj[SER1].SER1
    + geehh[[1,2]] conj[H].H.conj[SER1].SER2
    + geehh[[2,1]] conj[H].H.conj[SER2].SER1
    + geehh[[2,2]] conj[H].H.conj[SER2].SER2
    + glleea[[1,1]] conj[SEL1].SEL1.SER1.conj[SER1]
    + glleea[[1,2]] conj[SEL1].SEL1.SER2.conj[SER2]
    + glleea[[2,1]] conj[SEL2].SEL2.SER1.conj[SER2]
    + glleea[[2,2]] conj[SEL2].SEL2.SER2.conj[SER2]

);

(* part of slepton Lagrangian that needs to be conjugated *)
LagSlep = (
    + glleeb conj[SEL1].SER1.conj[SER2].SEL2
    + gleh[[1,1]] conj[H].SER1.SEL1
    + gleh[[2,1]] conj[H].SER1.SEL2
    + gleh[[1,2]] conj[H].SER2.SEL1
    + gleh[[2,2]] conj[H].SER2.SEL2
    - gdsle1 Hd.SEL1.e - gdsle2 Hd.SEL2.e
    - gdlse1 Hd.l.SER1 - gdlse2 Hd.l.SER2
    - gslwl1 conj[SEL1].FW.l - gslwl2 conj[SEL2].FW.l
    - gslbl1 conj[SEL1].FB.l - gslbl2 conj[SEL2].FB.l
    + gsebe1 conj[SER1].FB.e + gsebe2 conj[SER2].FB.e
);

DEFINITION[GaugeES][DiracSpinors] = {
    Bino -> {fB , conj[fB] },
    FW0  -> {fW0, 0},
    FWp  -> {fWp, conj[fWm]},
    FH0  -> {FHd0, conj[FHu0]},
    FHC  -> {FHdm, conj[FHup]},
    Fd1  -> {dL, 0},
    Fd2  -> {0, dR},
    Fu1  -> {uL, 0},
    Fu2  -> {0, uR},
    Fe1  -> {eL, 0},
    Fe2  -> {0, eR},
    Fv   -> {vL, 0}
};

(* ----- After EWSB ----- *)

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
    {{VB,VWB[3]}, {VP,VZ}, ZZ},
    {{VWB[1],VWB[2]}, {VWp,conj[VWp]}, ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
    {H0, {v,1}, {Ah,I/Sqrt[2]}, {hh,1/Sqrt[2]}}
};

(* ---- Mixings ---- *)

DEFINITION[EWSB][MatterSector] = {
    {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
    {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
    {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}},
    {{fB, fW0, FHd0, FHu0}, {L0, ZN}},
    {{{fWm, FHdm}, {fWp, FHup}}, {{Lm,UM}, {Lp,UP}}},
    {{SvL1, SvL2 }, {Sv, ZV}},
    {{SeL1, SeL2, SeR1, SeR2}, {Se, ZE}}
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
    Cha -> {Lm, conj[Lp]}
};
