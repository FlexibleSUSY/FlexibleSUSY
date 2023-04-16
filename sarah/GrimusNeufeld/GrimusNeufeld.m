Off[General::spell];

Model`Name = "GrimusNeufeld";
Model`NameLaTeX = "Inert-like two Higgs doublet model with 1 Heavy neutrino breaking the Z2 symmetry";

(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Superfields *)

Gauge[[1]] = {B, U[1], hypercharge, g1, False};
Gauge[[2]] = {WB, SU[2], left, g2, True};
Gauge[[3]] = {G, SU[3], color, g3, False};

(* Chiral Superfields *)

FermionFields[[1]] = {q, 3, {uL, dL}, 1 / 6, 2, 3};
FermionFields[[2]] = {l, 3, {vL, eL}, -1 / 2, 2, 1};
FermionFields[[3]] = {d, 3, conj[dR], 1 / 3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR], -2 / 3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR], 1, 1, 1};
FermionFields[[6]] = {n, 1, conj[nR], 0, 1, 1};

ScalarFields[[1]] = {H1, 1, {H1p, H10}, 1 / 2, 2, 1};
ScalarFields[[2]] = {H2, 1, {H2p, H20}, 1 / 2, 2, 1};

(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates = {GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][Additional] = {
  {LagHC, {AddHC -> True}},
  {LagYukawan, {AddHC -> True}},
  {LagYukawae, {AddHC -> True}},
  {LagNoHC, {AddHC -> False}}
};


LagNoHC = -(M112 conj[H1].H1 + M222 conj[H2].H2 + 1 / 2 Lambda1 conj[H1].H1.conj[H1].H1 +
    1 / 2 Lambda2 conj[H2].H2.conj[H2].H2 + Lambda3 conj[H2].H2.conj[H1].H1 + Lambda4 conj[H2].H1.conj[H1].H2 );

LagHC = -( Lambda5 / 2 conj[H2].H1.conj[H2].H1 + Yd1 conj[H1].d.q - Yu1 H1.u.q );

LagYukawan = -( - Yn1 H1.n.l - Yn2 H2.n.l + 1 / 2 Mm n.n );

LagYukawae = - ( Ye1 conj[H1].e.l  );

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = {
  {{VB, VWB[3]}, {VP, VZ}, ZZ},
  {{VWB[1], VWB[2]}, {VWm, conj[VWm]}, ZW}
};

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs] = {
  {H10, {v, 1 / Sqrt[2]}, {sigma1, \[ImaginaryI] / Sqrt[2]}, {phi1, 1 / Sqrt[2]}},
  {H20, {0, 0}, {sigma2, \[ImaginaryI] / Sqrt[2]}, {phi2, 1 / Sqrt[2]}}
};

DEFINITION[EWSB][MatterSector] = {
  {{phi1, phi2}, {hh, ZH}},
  {{sigma2, sigma1}, {Ah, ZA}},
  {{conj[H1p], conj[H2p]}, {Hm, ZHm}},
  {{{dL}, {conj[dR]}}, {{DL, Vd}, {DR, Ud}}},
  {{{uL}, {conj[uR]}}, {{UL, Vu}, {UR, Uu}}},
  {{{eL}, {conj[eR]}}, {{EL, Ve}, {ER, Ue}}},
  {{vL, conj[nR]}, {VL, Un}}
};

(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][DiracSpinors] = {
  Fd -> {DL, conj[DR]},
  Fe -> {EL, conj[ER]},
  Fu -> {UL, conj[UR]},
  Fv -> {VL, conj[VL]}};

DEFINITION[EWSB][GaugeES] = {
  Fd1 -> {dL, 0},
  Fd2 -> {0, dR},
  Fu1 -> {uL, 0},
  Fu2 -> {0, uR},
  Fe1 -> {eL, 0},
  Fe2 -> {0, eR},
  Fv1 -> {vL, 0},
  Fv2 -> {0, nR}
};

