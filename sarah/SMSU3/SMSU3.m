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

Model`Name = "SMSU3";
Model`NameLaTeX ="Standard Model";
Model`Authors = "F.Staub, A.Voigt";
Model`Date = "2014-11-06";

(* 2013-01-24: changed normalization of lambda term to convention of hep-ph/0207271 *)
(* 2013-06-24: using new name conventions (without inital "S" and "F" for scalar and matter fields) *)
(* 2013-09-01: changing to new conventions for FermionFields/MatterFields *)
(* 2014-11-06: Changed sign in Lagrangian to fit standard conventions *)
(* 2016-05-03: Changed sign of Yu *)



(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Groups *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False};
Gauge[[2]]={WB, SU[2], left,        g2,True};
Gauge[[3]]={G,  SU[3], color,       g3,False};
Gauge[[4]]={Co2,SU[3], color2,      g4,False};


(* Matter Fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3, 1};  
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1, 1};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3, 1};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3, 1};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1, 1};
FermionFields[[6]] = {f1,1, fL      ,       1, 1,  3, 3};
FermionFields[[7]] = {f2,1, fR      ,      -1, 1, -3,-3};

ScalarFields[[1]] =  {H, 1, {Hp, H0},     1/2, 2,  1, 1};


        
(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput]= {
	{LagHC, {AddHC->True}},
	{LagNoHC,{AddHC->False}}
};


LagNoHC = mu2 conj[H].H - 1/2 \[Lambda] conj[H].H.conj[H].H;
LagHC =  -(Yd conj[H].d.q + Ye conj[H].e.l + Yu u.q.H) \
         - mF1/2 conj[f1].f1 - mF2/2 conj[f2].f2 - mF3/2 f1.f2;



			  		  

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
};     
        
        
          	

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= 
{    {H0, {v, 1/Sqrt[2]}, {Ah, \[ImaginaryI]/Sqrt[2]},{hh, 1/Sqrt[2]}}     };
 

DEFINITION[EWSB][MatterSector]=   
    {{{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
     {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
     {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}};


(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  DL, conj[DR]},
 Fe ->{  EL, conj[ER]},
 Fu ->{  UL, conj[UR]},
 Ff ->{  fL, conj[fR]},
 Fv ->{  vL, 0}};

DEFINITION[EWSB][GaugeES]={
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};

DEFINITION[EWSB][Phases] = {
    {fL, PhaseFL},
    {fR, PhaseFR}
}; 

