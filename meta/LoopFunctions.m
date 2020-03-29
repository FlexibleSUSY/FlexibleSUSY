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

BeginPackage["LoopFunctions`"];
EndPackage[];

(* loop functions *)
{A0, B0, B1, B00, B11, B22, B22tilde, C0, C1, C2, D0, D27, F, G, H};

LFFull::usage = "Returns explicit form of loop functions with full
 momentum dependence.  The convention depends on the value of
 $BPMZSign .

 The following functions are implemented:

   A0[m,Q]                     [arxiv:hep-ph/9606211, Eq. (B.5)]
   B0[p,m1,m2,Q]               [arxiv:hep-ph/9606211, Eq. (B.7)]
   B1[p,m1,m2,Q]               [arxiv:hep-ph/9606211, Eq. (B.5)]
   B00[p,m1,m2,Q]              [arxiv:hep-ph/9606211, Eq. (B.10)]
   B11[p,m1,m2,Q]              [arxiv:hep-ph/0709.1075, Eq. (4.6)]
   B22[p,m1,m2,Q]              [arxiv:hep-ph/9606211, Eq. (B.10)]
   B22tilde[p,m1,m2,Q]         [arxiv:hep-ph/9606211, Eq. (B.14)]
   F[p,m1,m2,Q]                [arxiv:hep-ph/9606211, Eq. (B.11)]
   G[p,m1,m2,Q]                [arxiv:hep-ph/9606211, Eq. (B.12)]
   H[p,m1,m2,Q]                [arxiv:hep-ph/9606211, Eq. (B.13)]
   C0[0,0,m1,m2,m3,Q]          [arxiv:hep-ph/9606211, Eq. (C.19)]
   C1[0,0,m1,m2,m3,Q]
   C2[0,0,m1,m2,m3,Q]
   D0[0,0,0,m1,m2,m3,m4,Q]     [arxiv:hep-ph/9606211, Eq. (C.21)]
   D27[0,0,0,m1,m2,m3,m4,Q]    [arxiv:hep-ph/9606211, Eq. (C.22)]

B00[p,m1,m2,Q] is an alias for B22[p,m1,m2,Q].
";

LFZeroMomentum::usage = "Returns loop functions at zero momentum.  The
 convention depends on the value of $BPMZSign .";

LFDivergence::usage = "Returns divergent part of loop functions.  The
 convention depends on the value of $BPMZSign .";

LFScaleDependence::usage = "Returns the renormalization scale
 dependent part of loop functions.  The convention depends on the
 value of $BPMZSign .";

Delta::usage = "1/eps - \[Gamma]_E + Log[4 Pi]";

$BPMZSign::usage = "If set to 1, loop functions are returned in BPMZ
 convention (arxiv:hep-ph/9606211).  If set to -1, loop functions are
 returned in Denner convention (arxiv:hep-ph/0709.1075).  The default
 is -1.";

$BPMZSign = -1;

Begin["LoopFunctions`Private`"];

LFZeroMomentum[] := {
    A0[m_,mu_]                                               :> A0impl[m,mu],
    B0[0,m1_,m2_,mu_]             | B0[m1_,m2_,mu_]          :> B0zero[m1,m2,mu],
    B1[0,m1_,m2_,mu_]             | B1[m1_,m2_,mu_]          :> B1zero[m1,m2,mu],
    B00[0,m1_,m2_,mu_]            | B00[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    B11[0,m1_,m2_,mu_]            | B11[m1_,m2_,mu_]         :> B11zero[m1,m2,mu],
    B22[0,m1_,m2_,mu_]            | B22[m1_,m2_,mu_]         :> B22zero[m1,m2,mu],
    B22tilde[0,m1_,m2_,mu_]       | B22tilde[m1_,m2_,mu_]    :> B22tildezero[m1,m2,mu],
    C0[0,0,m1_,m2_,m3_,mu_]       | C0[m1_,m2_,m3_,mu_]      :> C0zero[m1,m2,m3],
    C1[0,0,m1_,m2_,m3_,mu_]       | C1[m1_,m2_,m3_,mu_]      :> C1zero[m1,m2,m3],
    C2[0,0,m1_,m2_,m3_,mu_]       | C2[m1_,m2_,m3_,mu_]      :> C2zero[m1,m2,m3],
    D0[0,0,0,m1_,m2_,m3_,m4_,mu_] | D0[m1_,m2_,m3_,m4_,mu_]  :> D0zero[m1,m2,m3,m4],
    D27[0,0,0,m1_,m2_,m3_,m4_,mu_]| D27[m1_,m2_,m3_,m4_,mu_] :> D27zero[m1,m2,m3,m4],
    F[0,m1_,m2_,mu_]              | F[m1_,m2_,mu_]           :> Fzero[m1,m2,mu],
    G[0,m1_,m2_,mu_]              | G[m1_,m2_,mu_]           :> Gzero[m1,m2,mu],
    H[0,m1_,m2_,mu_]              | H[m1_,m2_,mu_]           :> Hzero[m1,m2,mu]
};

LFFull[] := {
    A0[m_,mu_]                          :> A0impl[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> B0impl[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> B1impl[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> B22impl[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> B11impl[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> B22impl[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> B22tildeimpl[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> C0impl[p1,p2,m1,m2,m3,mu],
    C1[p1_,p2_,m1_,m2_,m3_,mu_]         :> C1impl[p1,p2,m1,m2,m3],
    C2[p1_,p2_,m1_,m2_,m3_,mu_]         :> C2impl[p1,p2,m1,m2,m3],
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> D0impl[p1,p2,p3,m1,m2,m3,m4,mu],
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> D27impl[p1,p2,p3,m1,m2,m3,m4,mu],
    F[p_,m1_,m2_,mu_]                   :> Fimpl[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> Gimpl[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> Himpl[p,m1,m2,mu]
};

LFDivergence[] := {
    A0[m_,mu_]                          :> DivA0[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> DivB0[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> DivB1[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> DivB22[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> DivB11[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> DivB22[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> DivB22tilde[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C1[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C2[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> 0,
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> 0,
    F[p_,m1_,m2_,mu_]                   :> DivF[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> DivG[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> DivH[p,m1,m2,mu]
};

LFScaleDependence[] := {
    A0[m_,mu_]                          :> LogA0[m,mu],
    B0[p_,m1_,m2_,mu_]                  :> LogB0[p,m1,m2,mu],
    B1[p_,m1_,m2_,mu_]                  :> LogB1[p,m1,m2,mu],
    B00[p_,m1_,m2_,mu_]                 :> LogB22[p,m1,m2,mu],
    B11[p_,m1_,m2_,mu_]                 :> LogB11[p,m1,m2,mu],
    B22[p_,m1_,m2_,mu_]                 :> LogB22[p,m1,m2,mu],
    B22tilde[p_,m1_,m2_,mu_]            :> LogB22tilde[p,m1,m2,mu],
    C0[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C1[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    C2[p1_,p2_,m1_,m2_,m3_,mu_]         :> 0,
    D0[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_] :> 0,
    D27[p1_,p2_,p3_,m1_,m2_,m3_,m4_,mu_]:> 0,
    F[p_,m1_,m2_,mu_]                   :> LogF[p,m1,m2,mu],
    G[p_,m1_,m2_,mu_]                   :> LogG[p,m1,m2,mu],
    H[p_,m1_,m2_,mu_]                   :> LogH[p,m1,m2,mu]
};

(********************* A0 *********************)

(* A0 [arxiv:hep-ph/9606211 Eq. (B.5)] *)
A0impl[m_, mu_] :=
    If[PossibleZeroQ[m],
       0,
       m^2 (Delta + 1 + Log[mu^2/m^2])
      ];

DivA0[m_, mu_] := m^2 Delta;

LogA0[m_, mu_] :=
    If[PossibleZeroQ[m],
       0,
       m^2 Log[mu^2/m^2]
      ];

(********************* B0 *********************)

(* B0 [arxiv:hep-ph/9606211 Eq. (B.6)] *)
B0impl[p_, m1_, m2_, mu_] :=
    Which[PossibleZeroQ[p],
          B0zero[m1,m2,mu],
          PossibleZeroQ[p - m2] && PossibleZeroQ[m1],
          Delta + Log[mu^2/m2^2] + 2,
          PossibleZeroQ[p - m1] && PossibleZeroQ[m2],
          Delta + Log[mu^2/m1^2] + 2,
          PossibleZeroQ[m1] && PossibleZeroQ[m2],
          Delta - Log[(-p^2 - I eps)/mu^2] + 2,
          True,
          B0analytic[p,m1,m2,mu]
      ];

(* B0 for p = 0 *)
B0zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1 - m2],
          Delta + Log[mu^2/m2^2],
          PossibleZeroQ[m1],
          Delta + 1 + Log[mu^2/m2^2],
          PossibleZeroQ[m2],
          Delta + 1 + Log[mu^2/m1^2],
          True,
          Delta + 1 + (m1^2 Log[mu^2/m1^2] - m2^2 Log[mu^2/m2^2])/(m1^2 - m2^2)
   ];

(* B0 general (symmetric) form [arxiv:hep-ph/9606211 Eq. (B.7)] *)
B0analytic[p_, m1_, m2_, mu_] :=
    Module[{eps, fB, s, xp, xm},
           s = p^2 - m2^2 + m1^2;
           xp = (s + Sqrt[s^2 - 4 p^2 (m1^2 - I eps)]) / (2 p^2);
           xm = (s - Sqrt[s^2 - 4 p^2 (m1^2 - I eps)]) / (2 p^2);
           fB[x_] := Log[1-x] - x Log[1 - 1/x] - 1;
           Normal @ Series[Delta - Log[p^2/mu^2] - fB[xp] - fB[xm],
                           {eps, 0, 0},
                           Assumptions :> p > 0 && m1 > 0 && m2 > 0 && mu > 0]
          ];

(* B0 with explicit integration [arxiv:hep-ph/9606211 Eq. (B.6)] *)
B0integral[p_, m1_, m2_, mu_] :=
    Module[{eps},
           Normal @ Series[
               Delta - Integrate[Log[((1-x) m1^2 + x m2^2 - x (1-x) p^2 - I eps)/mu^2], {x,0,1}],
               {eps, 0, 0}
           ]
          ];

DivB0[_, _, _, _] := Delta;

LogB0[p_, _, _, mu_] := Log[mu^2/p^2];

(********************* B1 *********************)

(* B1 [arxiv:hep-ph/9606211 Eq. (B.9)] *)
B1impl[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B1zero[m1,m2,mu],
       $BPMZSign 1/(2 p^2) (A0impl[m2,mu] - A0impl[m1,mu]
                            + (p^2 + m1^2 - m2^2) B0impl[p,m1,m2,mu])
      ];

(* B1 for p = 0 *)
B1zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1],
          $BPMZSign (1/4 + Delta/2 + Log[mu^2/m2^2]/2),
          PossibleZeroQ[m2],
          $BPMZSign (3/4 + Delta/2 + Log[mu^2/m1^2]/2),
          PossibleZeroQ[m1 - m2],
          $BPMZSign (Delta + Log[mu^2/m2^2])/2,
          True,
          $BPMZSign 1/2 (Delta + 1 + Log[mu^2/m2^2]
                         + (m1^2/(m1^2 - m2^2))^2 Log[m2^2/m1^2]
                         + 1/2 (m1^2 + m2^2)/(m1^2 - m2^2))
         ];

DivB1[_, _, _, _] := $BPMZSign Delta / 2;

LogB1[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       $BPMZSign Log[mu^2/m2^2]/2,
       $BPMZSign 1/(2 p^2) (LogA0[m2,mu] - LogA0[m1,mu]
                            + (p^2 + m1^2 - m2^2) LogB0[p,m1,m2,mu])
      ];

(********************* B11 *********************)

(* B11 *)
B11impl[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B11zero[m1,m2,mu],
       1/(6 p^2) (
           2 A0impl[m2,mu]
           - 2 m1^2 B0impl[p,m1,m2,mu]
           + $BPMZSign 4 (p^2 - m2^2 + m1^2) B1impl[p,m1,m2,mu]
           - m1^2 - m2^2 + p^2/3)
      ];

(* B11 for p = 0 *)
(* DirectedInfinity[(m1^2 - m2^2)*(-2*m1^2 + 2*m2^2)] *)
B11zero[m1_, m2_, mu_] := Undefined;

DivB11[p_, m1_, m2_, _] := Delta / 3;

LogB11[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       Undefined,
       1/(6 p^2) (
           2 LogA0[m2,mu]
           - 2 m1^2 LogB0[p,m1,m2,mu]
           + $BPMZSign 4 (p^2 - m2^2 + m1^2) LogB1[p,m1,m2,mu])
      ];

(********************* B22 (= B00) *********************)

(* B22 [arxiv:hep-ph/9606211 Eq. (B.10)],
   identical to B00[p,m1,m2,mu] *)
B22impl[p_, m1_, m2_, mu_] :=
    If[PossibleZeroQ[p],
       B22zero[m1,m2,mu],
       1/6 (1/2 (A0impl[m1,mu] + A0impl[m2,mu])
            + (m1^2 + m2^2 - p^2/2) B0impl[p,m1,m2,mu]
            + (m2^2 - m1^2)/(2 p^2) (A0impl[m2,mu] - A0impl[m1,mu]
                                     - (m2^2 - m1^2) B0impl[p,m1,m2,mu])
            + m1^2 + m2^2 - p^2/3)
      ];

(* B22 for p = 0 *)
B22zero[m1_, m2_, mu_] :=
    Which[PossibleZeroQ[m1] && PossibleZeroQ[m2],
          0,
          PossibleZeroQ[m1],
          (m2^2*(9/2 + 3*Delta + 3*Log[mu^2/m2^2]))/12,
          PossibleZeroQ[m2],
          (m1^2*(9/2 + 3*Delta + 3*Log[mu^2/m1^2]))/12,
          PossibleZeroQ[m1 - m2],
          (m2^2*(1 + Delta + Log[mu^2/m2^2]))/2,
          True,
          (A0impl[m2,mu] + m1^2 B0zero[m1, m2, mu] + (m1^2 + m2^2)/2)/4
      ];

DivB22[p_, m1_, m2_, _] := Delta (3*m1^2 + 3*m2^2 - p^2)/12;

LogB22[p_, m1_, m2_, mu_] :=
    Which[PossibleZeroQ[p] && PossibleZeroQ[m1 - m2],
          m2^2 Log[mu^2/m2^2]/2,
          PossibleZeroQ[p],
          (m1^2*(3*m1^2 + m2^2)*Log[mu^2/m1^2] -
           m2^2*(m1^2 + 3*m2^2)*Log[mu^2/m2^2])/(12*(m1^2 - m2^2)),
          True,
          1/6 (1/2 (LogA0[m1,mu] + LogA0[m2,mu])
               + (m1^2 + m2^2 - p^2/2) LogB0[p,m1,m2,mu]
               + (m2^2 - m1^2)/(2 p^2) (LogA0[m2,mu] - LogA0[m1,mu]
                                        - (m2^2 - m1^2) LogB0[p,m1,m2,mu]))
         ];

(********************* C0 *********************)

C0impl[p1_, p2_, m1_, m2_, m3_, mu_] :=
    If[PossibleZeroQ[p1] && PossibleZeroQ[p2],
       C0zero[m1,m2,m3],
       C0analytic[p1, p2, m1, m2, m3, mu]
      ];

(* C0 for p = 0 [arxiv:hep-ph/9606211 Eq. (C.19)] *)
C0zero[m1_, m2_, m3_] := Which[
    PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3],
    -1/(2*m3^2),
    PossibleZeroQ[m1 - m2],
    (-m2^2 + m3^2 - m3^2*Log[m3^2/m2^2])/(m2^2 - m3^2)^2,
    PossibleZeroQ[m1 - m3],
    (m2^2 - m3^2 - m2^2*Log[m2^2/m3^2])/(m2^2 - m3^2)^2,
    PossibleZeroQ[m2 - m3],
    (m1^2 - m3^2 + m1^2*Log[m3^2/m1^2])/(m1^2 - m3^2)^2,
    (* general case *)
    True,
    (+ m2^2 / (m1^2 - m2^2) Log[m2^2/m1^2]
     - m3^2 / (m1^2 - m3^2) Log[m3^2/m1^2]) / (m2^2 - m3^2)
];

(* C0 for complex momenta [arxiv:hep-ph/0709.1075, Eq. (4.26)] *)
(* Note: This implementation is divergent for real momenta.    *)
(* For real momenta a different implementation should be used. *)
C0analytic[p1_, p2_, m1_, m2_, m3_, mu_] :=
    Module[{p21 = (p2 - p1), p12 = (p1 - p2), pjk, pki, pij, mi, mj, mk,
            Dilogs, y0, xi, yi, alpha, alphai, eps, kappa, eta, result},
           (* Källén function, [arxiv:hep-ph/0709.1075, Eq. (4.28)] *)
           kappa[x_, y_, z_] := Sqrt[x^2 + y^2 + z^2 - 2 (x y + y z + z x)];
           (* [arxiv:hep-ph/0709.1075, Eq. (4.30)] *)
           eta[a_, b_] := Log[a b] - Log[a] - Log[b];

           pjk[i_?IntegerQ] :=
               Which[i == 0, p12, (* j = 1, k = 2 *)
                     i == 1, p2,  (* j = 2, k = 0 *)
                     i == 2, p1   (* j = 0, k = 1 *)];
           pki[i_?IntegerQ] :=
               Which[i == 0, p2,  (* j = 1, k = 2 *)
                     i == 1, p1,  (* j = 2, k = 0 *)
                     i == 2, p12  (* j = 0, k = 1 *)];
           pij[i_?IntegerQ] :=
               Which[i == 0, p1,  (* j = 1, k = 2 *)
                     i == 1, p12, (* j = 2, k = 0 *)
                     i == 2, p2   (* j = 0, k = 1 *)];
           mi[i_?IntegerQ] := Which[i == 0, m1,
                                    i == 1, m2,
                                    i == 2, m3];
           mj[i_?IntegerQ] := Which[i == 0, m2, (* j = 1, k = 2 *)
                                    i == 1, m3, (* j = 2, k = 0 *)
                                    i == 2, m1  (* j = 0, k = 1 *)];
           mk[i_?IntegerQ] := Which[i == 0, m3, (* j = 1, k = 2 *)
                                    i == 1, m1, (* j = 2, k = 0 *)
                                    i == 2, m2  (* j = 0, k = 1 *)];
           y0[i_] := (pjk[i]^2 (pjk[i]^2 - pki[i]^2 - pij[i]^2 + 2 mi[i]^2 - mj[i]^2 - mk[i]^2)
                      - (pki[i]^2 - pij[i]^2) (mj[i]^2 - mk[i]^2)
                      + alpha (pjk[i]^2 - mj[i]^2 + mk[i]^2)) / (2 alpha pjk[i]^2);
           xi[i_, s_] := (pjk[i]^2 - mj[i]^2 - mk[i]^2 + s alphai[i]) / (2 pjk[i]^2);
           yi[i_, s_] := y0[i] - xi[i,s];
           alpha = kappa[p1^2, p21^2, p2^2];
           alphai[i_] := kappa[pjk[i]^2, mj[i]^2, mk[i]^2] (1 + I eps pjk[i]^2);
           Dilogs[i_, s_] := (
               PolyLog[2,(y0[i] - 1)/yi[i,s]] - PolyLog[2,y0[i]/yi[i,s]]
               + eta[1 - xi[i,s], 1/yi[i,s]] Log[(y0[i] - 1)/yi[i,s]]
               - eta[-xi[i,s], 1/yi[i,s]] Log[y0[i]/yi[i,s]]);

           result = Sum[Dilogs[i,+1] + Dilogs[i,-1]
                        - (eta[-xi[i,+1],-xi[i,-1]]
                           - eta[yi[i,+1],yi[i,-1]]
                           - 2 Pi I UnitStep[-Re[pjk[i]^2]] UnitStep[-Im[yi[i,+1] yi[i,-1]]]
                          ) Log[(1 - y0[i])/(-y0[i])],
                        {i,0,2}] / alpha;

           Normal @ Series[result, {eps, 0, 0},
                           Assumptions :> Element[p1, Complexes] || Element[p2, Complexes]]
          ];

(********************* C1 *********************)

C1impl[p1_, p2_, m1_, m2_, m3_] :=
    If[PossibleZeroQ[p1] && PossibleZeroQ[p2],
       C1zero[m1,m2,m3],
       NotImplemented
      ];

(* C1 for p = 0 *)
C1zero[m1_, m2_, m3_] :=
    Module[{t1 = m2^2/m1^2, t2 = m3^2/m1^2},
           Which[
               PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3],
               -1/(6 m1^2),
               PossibleZeroQ[m1 - m2],
               -(m2^4 - 4*m2^2*m3^2 + 3*m3^4 - 2*m3^4*Log[m3^2/m2^2])/(4*(m2^2 - m3^2)^3),
               PossibleZeroQ[m2 - m3],
               (3*m1^4 - 4*m1^2*m3^2 + m3^4 + 2*m1^4*Log[m3^2/m1^2])/(4*(m1^2 - m3^2)^3),
               PossibleZeroQ[m1 - m3],
               (-m2^4 + m3^4 + 2*m2^2*m3^2*Log[m2^2/m3^2])/(2*(m2^2 - m3^2)^3),
               (* general case *)
               True,
               -(t1 / (2 (t1 - 1) (t1 - t2))
                 - t1 (t1 - 2 t2 + t1 t2) Log[t1] / (2 (t1 - 1)^2 (t1 - t2)^2)
                 + (t2^2 - 2 t1 t2^2 + t1^2 t2^2) Log[t2] / (2 (t1 - 1)^2 (t1 - t2)^2 (t2 - 1))
                )/m1^2
                ]
          ];

(********************* C2 *********************)

C2impl[p1_, p2_, m1_, m2_, m3_] :=
    If[PossibleZeroQ[p1] && PossibleZeroQ[p2],
       C2zero[m1,m2,m3],
       NotImplemented
      ];

(* C2 for p = 0 *)
C2zero[m1_, m2_, m3_] :=
    Module[{t1 = m2^2/m1^2, t2 = m3^2/m1^2},
           Which[
               PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3],
               -1/(6 m1^2),
               PossibleZeroQ[m1 - m2],
               (-m2^4 + m3^4 + 2*m2^2*m3^2*Log[m2^2/m3^2])/(2*(m2^2 - m3^2)^3),
               PossibleZeroQ[m2 - m3],
               (3*m1^4 - 4*m1^2*m3^2 + m3^4 + 2*m1^4*Log[m3^2/m1^2])/(4*(m1^2 - m3^2)^3),
               PossibleZeroQ[m1 - m3],
               (3*m2^4 - 4*m2^2*m3^2 + m3^4 - 2*m2^4*Log[m2^2/m3^2])/(4*(m2^2 - m3^2)^3),
               (* general case *)
               True,
               -(- t2 / (2 (t1 - t2) (t2 - 1))
                 + Log[t1] / (2 (t1 - 1) (t2 - 1)^2)
                 + (2 t1 t2 - 2 t1^2 t2 - t2^2 + t1^2 t2^2) Log[t1/t2] / (2 (t1 - 1) (t1 - t2)^2 (t2 - 1)^2)
                )/m1^2
                ]
          ];

(********************* D0 *********************)

(* D0 *)
D0impl[p1_, p2_, p3_, m1_, m2_, m3_, m4_, mu_] :=
    If[PossibleZeroQ[p1] && PossibleZeroQ[p2] && PossibleZeroQ[p3],
       D0zero[m1, m2, m3, m4],
       NotImplemented
      ];

(* D0 for p = 0 [arxiv:hep-ph/9606211 Eq. (C.21)] *)
D0zero[m1_, m2_, m3_, m4_] := Which[
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3] && 
    PossibleZeroQ[m1 - m4],
   1/(6 m1^4),
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m2 - m3],
   (8*(m2^2 - m4^2)^2*(m2^2 + m4^2) + (-23*m2^4*m4^2 - 10*m2^2*m4^4 + 
         m4^6)*Log[m2^2/m4^2] + (-7*m2^4*m4^2 - 26*m2^2*m4^4 + m4^6)*
       Log[m4^2/m2^2])/(16*m2^2*(m2^2 - m4^2)^4),
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m2 - m4],
   (8*(m2^2 - m3^2)^2*(m2^2 + m3^2) + (7*m2^4*m3^2 + 26*m2^2*m3^4 - 
         m3^6)*Log[m2^2/m3^2] + (23*m2^4*m3^2 + 10*m2^2*m3^4 - m3^6)*
       Log[m3^2/m2^2])/(16*m2^2*(m2^2 - m3^2)^4),
   PossibleZeroQ[m1 - m2] && PossibleZeroQ[m3 - m4],
   (-8*(m2^2 - m4^2)^2 + (m2^4 - 8*m2^2*m4^2 - 5*m4^4)*
       Log[m2^2/m4^2] - (3*m2^4 + 8*m2^2*m4^2 + m4^4)*
       Log[m4^2/m2^2])/(4*(m2^2 - m4^2)^4),
   PossibleZeroQ[m1 - m2],
   (-((m2 - m3)*(m2 + m3)*(m2 - m4)*(m3 - m4)*(m2 + m4)*(m3 + 
           m4)) - (m2^4 - m3^2*m4^2)*(m3^2*Log[m3^2/m2^2] + 
         m4^2*Log[m2^2/m4^2]) + 
      m3^2*m4^2*(-2*m2^2 + m3^2 + m4^2)*
       Log[m4^2/m3^2])/((m2^2 - m3^2)^2*(m2^2 - m4^2)^2*(m3^2 - m4^2)),
   True,
   (C0zero[m1, m3, m4] - C0zero[m2, m3, m4])/(m1^2 - m2^2)];

(********************* D27 *********************)

(* D27 *)
D27impl[p1_, p2_, p3_, m1_, m2_, m3_, m4_, mu_] :=
    If[PossibleZeroQ[p1] && PossibleZeroQ[p2] && PossibleZeroQ[p3],
       D27zero[m1, m2, m3, m4],
       NotImplemented
      ];

(* D27 for p = 0 [arxiv:hep-ph/9606211 Eq. (C.22)] *)
D27zero[m1_, m2_, m3_, m4_] :=
   (m1^2 C0zero[m1, m3, m4] - m2^2 C0zero[m2, m3, m4])/(4 (m1^2 - m2^2));

(********************* F *********************)

(* F [arxiv:hep-ph/9606211 Eq. (B.11)] *)
Fimpl[p_, m1_, m2_, mu_] :=
    A0impl[m1,mu] - 2 A0impl[m2,mu] - (2 p^2 + 2 m1^2 - m2^2) B0impl[p,m1,m2,mu];

Fzero[m1_, m2_, mu_] := Fimpl[0,m1,m2,mu];

DivF[p_, m1_, m2_, mu_] := Delta (-m1^2 - m2^2 - 2*p^2);

LogF[p_, m1_, m2_, mu_] :=
    LogA0[m1,mu] - 2 LogA0[m2,mu] - (2 p^2 + 2 m1^2 - m2^2) LogB0[p,m1,m2,mu];

(********************* G *********************)

(* G [arxiv:hep-ph/9606211 Eq. (B.12)] *)
Gimpl[p_, m1_, m2_, mu_] :=
    (p^2 - m1^2 - m2^2) B0impl[p,m1,m2,mu] - A0impl[m1,mu] - A0impl[m2,mu];

Gzero[m1_, m2_, mu_] := Gimpl[0,m1,m2,mu];

DivG[p_, m1_, m2_, mu_] := Delta (-2*m1^2 - 2*m2^2 + p^2);

LogG[p_, m1_, m2_, mu_] :=
    (p^2 - m1^2 - m2^2) LogB0[p,m1,m2,mu] - LogA0[m1,mu] - LogA0[m2,mu];

(********************* H *********************)

(* H [arxiv:hep-ph/9606211 Eq. (B.13)] *)
Himpl[p_, m1_, m2_, mu_] := 4 B22impl[p,m1,m2,mu] + Gimpl[p,m1,m2,mu];

Hzero[m1_, m2_, mu_] := Himpl[0,m1,m2,mu];

DivH[p_, m1_, m2_, mu_] := Delta (-m1^2 - m2^2 + (2*p^2)/3);

LogH[p_, m1_, m2_, mu_] := 4 LogB22[p,m1,m2,mu] + LogG[p,m1,m2,mu];

(********************* ~B22 *********************)

(* ~B22 [arxiv:hep-ph/9606211 Eq. (B.14)] *)
B22tildeimpl[p_, m1_, m2_, mu_] := B22impl[p,m1,m2,mu] - A0impl[m1,mu]/4 - A0impl[m2,mu]/4;

B22tildezero[m1_, m2_, mu_] := B22tildeimpl[0,m1,m2,mu];

DivB22tilde[p_, m1_, m2_, mu_] := Delta (-p^2/12);

LogB22tilde[p_, m1_, m2_, mu_] := LogB22[p,m1,m2,mu] - LogA0[m1,mu]/4 - LogA0[m2,mu]/4;

End[];
