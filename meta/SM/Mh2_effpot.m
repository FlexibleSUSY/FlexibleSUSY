(* Effective Higgs potential in the SM from (4.34)-(4.38), (5.3) of [1508.00912] *)
(* k = 1/(4 Pi)^2 *)

(* Use SUSYHD/HSSUSY normalization of \[Lambda] *)
V0L = mu2/2 h^2 + (1/2) \[Lambda]/4 h^4;

V1L = -3 mt^4 (l - 3/2);

V2L = g3^2 Nc CF mt^4 (6 l^2 - 16 l + 18);

V3L = g3^4 Nc CF mt^4 (
    CG (-22/3 l^3 + 185/3 l^2 + (24 Zeta[3] - 1111/6) l
        + 2609/12 + 44/45 Pi^4 - 232/3 Zeta[3] + 16/3 Log[2]^2 (Pi^2 - Log[2]^2) - 128 PolyLog[4,1/2])
    + CF (-24 l^3 + 63 l^2 - (48 Zeta[3] + 121/2) l + 85/12 - 88/45 Pi^4
          + 192 Zeta[3] - 32/3 Log[2]^2 (Pi^2 - Log[2]^2) + 256 PolyLog[4,1/2])
    + TF (48 l - 232/3 + 96 Zeta[3])
    + TF nq (8/3 l^3 - 52/3 l^2 + 142/3 l - 161/3 - 64/3 Zeta[3])
);

V4L = g3^6 mt^4 (a4 l^4 + a3 l^3 + a2 l^2 + a1 l + a0);

Veff[h_] := V0L + k V1L + k^2 V2L + k^3 V3L + k^4 V4L /.
    { Nc -> 3, CF -> 4/3, CG -> 3, TF -> 1/2, nq -> 6, l -> Log[mt^2/Q^2] } /.
    { mt -> h yt/Sqrt[2] };

(* tadpole / v *)
tad = -1/v D[Veff[h], h] /. { h -> v } /. { v -> Sqrt[2] mt/yt };

(* self energy *)
sel = -D[Veff[h],{h,2}] /. { h -> Sqrt[2] mt/yt };

(* Delta Mh^2 *)
DMh2 = tad - sel /. {
    a4 -> aL4, (* 1380 *)
    a3 -> aL3, (* -9144 *)
    a2 -> aL2, (* 27699.06 *)
    a1 -> aL1, (* -54056.36 *)
    a0 -> aL0  (* 59366.97 *)
} // Collect[#, {k,Log[__]}, Simplify]&;

(* coefficients from Eq.(4.39) *)
aL4 = 1380;

aL3 = -9144;

aL2 = 30584 - 2400 Zeta[3];

aL1 = 27680/3 Zeta[3] - 63200/9 Zeta[5] - 1547146/27 - 208/9 Pi^4 +
    640/3 Log[2]^2 (Log[2]^2 - Pi^2) + 5120 PolyLog[4,1/2];

aL0 = 13820381/270 + 1747112 Zeta[3]/45 + 1984 Zeta[5]/9 -
    40288/9 Zeta[3]^2 - 298894/1215 Pi^4 - 1780/243 Pi^6 +
    5888/135 Log[2]^5 - 5888/81 Pi^2 Log[2]^3 - 36064/405 Pi^4 Log[2] +
    78464/81 Log[2]^2 (Log[2]^2 - Pi^2) + 627712/27 PolyLog[4,1/2] -
    47104/9 PolyLog[5,1/2];
