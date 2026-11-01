Off[General::spell]

Model`Name = "LQS1_up_basis";
Model`NameLaTeX = "Scalar Leptoquark S1 (up-type quark mass eigenstate basis)";
Model`Authors = "U. Khasianevich, J. Wünsche";
Model`Date = "2021-10-20";

SetOptions[MakeAll, IncludeCalcHep-> False, IncludeWHIZARD-> False,
   IncludeSPheno-> False, IncludeUFO-> False, IncludeTeX-> False];

(* Implementing the model based on eq. 18 of 2104.03691. *)
(* Note: we use conjugated leptoquark. *)

(*------------------*)
(* Particle Content *)
(*------------------*)

(* Gauge Groups *)

Gauge[[1]]={B,   U[1], hypercharge, g1, False};
Gauge[[2]]={WB, SU[2], left,        g2, True};
Gauge[[3]]={G,  SU[3], color,       g3, False};

(* SM Fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3};
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1};
FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1};

ScalarFields[[1]] =  {H, 1, {Hp, H0},     1/2, 2,  1};

(* Scalar Leptoquark S1 *)
ScalarFields[[2]] =  {sc, 1, S1c,        -1/3, 1,  3};

(*------------*)
(* DEFINITION *)
(*------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput] =
{
   {LagSMYuk, {AddHC-> True}},
   {LagSMHiggs, {AddHC-> False}},
   {LagLQNoHC, {AddHC-> False}},
   {LagLQ, {AddHC-> True}}
};

LagSMYuk = -(Yd conj[H].d.q + Ye conj[H].e.l + Yu u.q.H );
LagSMHiggs = -mu2 conj[H].H - 1/2 \[Lambda] conj[H].H.conj[H].H ;
LagLQ = -(lamQL q.l.conj[sc] + lamql u.e.sc );
LagLQNoHC = -(MS12 conj[sc].sc + gHS1 conj[H].H.conj[sc].sc
                             + lamS1/2 conj[sc].sc.conj[sc].sc );

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{
  {{VB, VWB[3]}, {VP, VZ}, ZZ},
  {{VWB[1], VWB[2]},{VWp, conj[VWp]}, ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] =
{
   {H0, {v, 1/Sqrt[2]}, {Ah, I/Sqrt[2]}, {hh, 1/Sqrt[2]}}
};


DEFINITION[EWSB][MatterSector] =
{
   {{{dL}, {conj[dR]}}, {{DL, Vd}, {DR, Ud}}},
   {{{uL}, {conj[uR]}}, {{UL, Vu}, {UR, Uu}}},
   {{{eL}, {conj[eR]}}, {{EL, Ve}, {ER, Ue}}}
};

(*---------------*)
(* Dirac-Spinors *)
(*---------------*)

DEFINITION[EWSB][DiracSpinors]=
{
   Fd-> {DL, conj[DR]},
   Fe-> {EL, conj[ER]},
   Fu-> {UL, conj[UR]},
   Fv-> {vL, 0}
};

DEFINITION[EWSB][GaugeES]=
{
   Fd1-> {FdL, 0},
   Fd2-> {0, FdR},
   Fu1-> {FuL, 0},
   Fu2-> {0, FuR},
   Fe1-> {FeL, 0},
   Fe2-> {0, FeR}
};
