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

(* Implementation of THDM 1-loop threshold corrections from arxiv:0901.2065 *)

BeginPackage["THDMThresholds1L`"];
EndPackage[];

(* parameters *)
{gY, g2, Yu, Yd, Ye, Tu, Td, Te, msu, msd, mse, msq, msl, M1, M2, Mu, Q};

(* loop functions *)
{B0, DB0, C0, D0, D2tilde, D4tilde, W};

(* flags *)
{ flagSferm, flagIno, flagZdd, flagZud, flagZuu, flagdg, flagMSDRg2, flagMSDRlam };

(* options *)
{ coefficients, flags, loopFunctions, loopOrder, sumHead };

GetTHDMThresholds1L::usage = "Returns list of 1L threshold corrections
 for lambda_1 ... lambda_7 using the convention of
 arxiv:hep-ph/9307201 and arxiv:1508.00576.  The expressions have been
 taken from arxiv:0901.2065.

 Important note: The Yukawa couplings are transposed compared to
 SARAH's MSSM and THDM model files.

 To prevent the replacement of the coefficients by their explicite
 form, run

    GetTHDMThresholds1L[coefficients -> {}]

 To prevent the replacement of the flags marking the different
 contributions, run

    GetTHDMThresholds1L[flags -> {}]

 To prevent expanding the sums, run for example

    GetTHDMThresholds1L[sumHead -> Summation]

 To select the loop order, run

    GetTHDMThresholds1L[loopOrder -> {1,1}]

 {0,0} means no tree-level and no 1-loop.
 {1,0} means only tree-level.
 {0,1} means only 1-loop.
 {1,1} means tree-level plus 1-loop.

 By default, the loop functions are kept unevaluated.";

GetTHDMThresholds1LCoefficients::usage = "Returns list of replacement
 rules for the coefficients defined in Eq. (122), (124) and Tables 4-9
 of arxiv:0901.2065.";

GetTHDMThresholds1LFlags::usage = "Returns list of replacement rules
 for the flags which mark the different contributions.";

GetTHDMThresholds1LLoopFunctions::usage = "Returns list of replacement
 rules for the loop functions defined in Eq. (B16)-(B17) of
 Phys. Rev. D 84 (2011) 034030.

 Note: The loop functions B0' and W in Eq. (130) arxiv:0901.2065v1
 carry the wrong global sign.
";

Begin["THDMThresholds1L`Private`"];

GetTHDMThresholds1LFlags[] := {
    flagSferm    ->  1, (* Enable/disable sfermion contribution *)
    flagIno      ->  1, (* Enable/disable gaugino + Higgsino contribution *)
    flagZdd      ->  1, (* Enable/disable field renormalization Zdd *)
    flagZud      ->  0, (* Enable/disable field renormalization Zud *)
    flagZuu      ->  1, (* Enable/disable field renormalization Zuu *)
    flagdg       ->  1, (* Enable/disable gauge coupling renormalization *)
    flagMSDRg2   ->  1, (* Enable/disable DREG <-> DRED conversion of g2 *)
    flagMSDRlam  ->  1  (* Enable/disable DREG <-> DRED conversion of lambda_i *)
};

lamBar = lamHat = lamTree = lamIno = lamSferm = Table[Undef, {i, 1, 7}];

GetTHDMThresholds1LCoefficients[] := {
    Nc -> 3,
    kappa -> 1/(16 Pi^2),

    (* Eq. (122) *)
    as[1] -> -3/4,
    as[2] -> -3/4,
    as[3] -> -3/4,
    as[4] -> 0,
    as[5] -> 0,
    as[6] -> 0,
    as[7] -> 0,
    asp[1] -> -1/2,
    asp[2] -> -1/2,
    asp[3] -> -1/2,
    asp[4] -> 1,
    asp[5] -> 0,
    asp[6] -> 0,
    asp[7] -> 0,
    aspp[1] -> -1/4,
    aspp[2] -> -1/4,
    aspp[3] -> -1/4,
    aspp[4] -> 0,
    aspp[5] -> 0,
    aspp[6] -> 0,
    aspp[7] -> 0,

    (* Table 4 *)
    a2[1] -> Abs[M2]^2/2,
    a2[2] -> Abs[M2]^2/2,
    a2[3] -> 3 Abs[Mu]^2 + 5 Abs[M2]^2/2,
    a2[4] -> -3 Abs[Mu]^2 - 2 Abs[M2]^2,
    a2[6] -> 3 Mu M2,
    a2[7] -> 3 Mu M2,

    a4[1] -> 5/2,
    a4[2] -> 5/2,
    a4[3] -> 1/2,
    a4[4] -> 2,
    a4[6] -> 0,
    a4[7] -> 0,

    a2p[1] -> (M1 Conjugate[M2] + Conjugate[M1] M2)/2,
    a2p[2] -> (M1 Conjugate[M2] + Conjugate[M1] M2)/2,
    a2p[3] -> 2 Abs[Mu]^2 + (M1 Conjugate[M2] + Conjugate[M1] M2)/2,
    a2p[4] -> 2 Abs[Mu]^2 - M1 Conjugate[M2] - Conjugate[M1] M2,
    a2p[6] -> Mu (M1 + M2),
    a2p[7] -> Mu (M1 + M2),

    a4p[1] -> 1,
    a4p[2] -> 1,
    a4p[3] -> 1,
    a4p[4] -> -2,
    a4p[6] -> 0,
    a4p[7] -> 0,

    a2pp[1] -> Abs[M1]^2/2,
    a2pp[2] -> Abs[M1]^2/2,
    a2pp[3] -> Abs[Mu]^2 + Abs[M1]^2/2,
    a2pp[4] -> -Abs[Mu]^2,
    a2pp[6] -> Mu M1,
    a2pp[7] -> Mu M1,

    a4pp[1] -> 1/2,
    a4pp[2] -> 1/2,
    a4pp[3] -> 1/2,
    a4pp[4] -> 0,
    a4pp[6] -> 0,
    a4pp[7] -> 0,

    (* Table 7 *)
    c1p[6] -> -gY^2/2,
    c2p[6] -> 1,
    c3p[6] -> (gY^2 - g2^2)/4,
    c4p[6] -> 1,
    c5p[6] -> -gY^2/2,
    c6p[6] -> 3,
    c7p[6] -> (-gY^2 - 3 g2^2)/4,
    c8p[6] -> 3,
    c9p[6] -> (3 g2^2 - gY^2)/4,
    c10p[6] -> 0,
    c11p[6] -> gY^2,
    c12p[6] -> 0,

    c1p[7] -> gY^2/2,
    c2p[7] -> 0,
    c3p[7] -> (-gY^2 + g2^2)/4,
    c4p[7] -> 0,
    c5p[7] -> gY^2/2,
    c6p[7] -> 0,
    c7p[7] -> (gY^2 + 3 g2^2)/4,
    c8p[7] -> 0,
    c9p[7] -> (gY^2 - 3 g2^2)/4,
    c10p[7] -> 3,
    c11p[7] -> -gY^2,
    c12p[7] -> 3,

    (* Table 8 and 9 *)
    b1[1] -> -gY^4/4,
    b2[1] -> gY^2,
    b3[1] -> -1,
    b4[1] -> (-g2^4 - gY^4)/8,
    b5[1] -> (g2^2 - gY^2)/2,
    b6[1] -> -1,
    b7[1] -> -gY^4/12,
    b8[1] -> gY^2,
    b9[1] -> -3,
    b10[1] -> 0,
    b11[1] -> -gY^4/3,
    b12[1] -> 0,
    b13[1] -> 0,
    b14[1] -> (-9 g2^4 - gY^4)/24,
    b15[1] -> (3 g2^2 + gY^2)/2,
    b16[1] -> 0,
    b17[1] -> -3,
    b18[1] -> 0,
    b19[1] -> 0,

    b1[2] -> -gY^4/4,
    b2[2] -> 0,
    b3[2] -> 0,
    b4[2] -> (-g2^4 - gY^4)/8,
    b5[2] -> 0,
    b6[2] -> 0,
    b7[2] -> -gY^4/12,
    b8[2] -> 0,
    b9[2] -> 0,
    b10[2] -> 0,
    b11[2] -> -gY^4/3,
    b12[2] -> 2 gY^2,
    b13[2] -> -3,
    b14[2] -> (-9 g2^4 - gY^4)/24,
    b15[2] -> 0,
    b16[2] -> (3 g2^2 - gY^2)/2,
    b17[2] -> 0,
    b18[2] -> -3,
    b19[2] -> 0,

    b1[3] -> gY^4/4,
    b2[3] -> -gY^2/2,
    b3[3] -> 0,
    b4[3] -> (g2^4 + gY^4)/8,
    b5[3] -> (gY^2 - g2^2)/4,
    b6[3] -> 0,
    b7[3] -> gY^4/12,
    b8[3] -> -gY^2/2,
    b9[3] -> 0,
    b10[3] -> 0,
    b11[3] -> gY^4/3,
    b12[3] -> -gY^2,
    b13[3] -> 0,
    b14[3] -> (9 g2^4 + gY^4)/24,
    b15[3] -> (-3 g2^2 - gY^2)/4,
    b16[3] -> (-3 g2^2 + gY^2)/4,
    b17[3] -> 0,
    b18[3] -> 0,
    b19[3] -> 0,

    b1[4] -> 0,
    b2[4] -> 0,
    b3[4] -> 0,
    b4[4] -> -g2^4/4,
    b5[4] -> g2^2/2,
    b6[4] -> 0,
    b7[4] -> 0,
    b8[4] -> 0,
    b9[4] -> 0,
    b10[4] -> -3,
    b11[4] -> 0,
    b12[4] -> 0,
    b13[4] -> 0,
    b14[4] -> -3 g2^4/4,
    b15[4] -> 3 g2^2/2,
    b16[4] -> 3 g2^2/2,
    b17[4] -> 0,
    b18[4] -> 0,
    b19[4] -> -3,

    c1[1] -> 0,
    c2[1] -> gY^2,
    c3[1] -> 0,
    c4[1] -> -2,
    c5[1] -> 0,
    c6[1] -> (g2^2 - gY^2)/2,
    c7[1] -> 0,
    c8[1] -> -2,
    c9[1] -> 0,
    c10[1] -> gY^2,
    c11[1] -> 0,
    c12[1] -> -6,
    c13[1] -> 0,
    c14[1] -> 0,
    c15[1] -> (3 g2^2 + gY^2)/2,
    c16[1] -> -6,
    c17[1] -> 0,
    c18[1] -> 0,
    c19[1] -> (gY^2 - 3 g2^2)/2,
    c20[1] -> 0,
    c21[1] -> 0,
    c22[1] -> 0,
    c23[1] -> 0,
    c24[1] -> -2 gY^2,
    c25[1] -> 0,
    c26[1] -> 0,
    c27[1] -> 0,

    c1[2] -> -gY^2,
    c2[2] -> 0,
    c3[2] -> 0,
    c4[2] -> 0,
    c5[2] -> (gY^2 - g2^2)/2,
    c6[2] -> 0,
    c7[2] -> 0,
    c8[2] -> 0,
    c9[2] -> -gY^2,
    c10[2] -> 0,
    c11[2] -> 0,
    c12[2] -> 0,
    c13[2] -> -(3 g2^2 + gY^2)/2,
    c14[2] -> 0,
    c15[2] -> 0,
    c16[2] -> 0,
    c17[2] -> 0,
    c18[2] -> 0,
    c19[2] -> 0,
    c20[2] -> 0,
    c21[2] -> (3 g2^2 - gY^2)/2,
    c22[2] -> -6,
    c23[2] -> 0,
    c24[2] -> 0,
    c25[2] -> 0,
    c26[2] -> 2 gY^2,
    c27[2] -> -6,

    c1[3] -> gY^2/2,
    c2[3] -> -gY^2/2,
    c3[3] -> -1,
    c4[3] -> 0,
    c5[3] -> (-gY^2 + g2^2)/4,
    c6[3] -> (gY^2 - g2^2)/4,
    c7[3] -> -1,
    c8[3] -> 0,
    c9[3] -> gY^2/2,
    c10[3] -> -gY^2/2,
    c11[3] -> -3,
    c12[3] -> 0,
    c13[3] -> (3 g2^2 + gY^2)/4,
    c14[3] -> -3,
    c15[3] -> (-3 g2^2 - gY^2)/4,
    c16[3] -> 0,
    c17[3] -> 0,
    c18[3] -> 0,
    c19[3] -> (3 g2^2 - gY^2)/4,
    c20[3] -> -3,
    c21[3] -> (-3 g2^2 + gY^2)/4,
    c22[3] -> 0,
    c23[3] -> 0,
    c24[3] -> gY^2,
    c25[3] -> -3,
    c26[3] -> -gY^2,
    c27[3] -> 0,

    c1[4] -> 0,
    c2[4] -> 0,
    c3[4] -> 0,
    c4[4] -> 0,
    c5[4] -> -g2^2/2,
    c6[4] -> g2^2/2,
    c7[4] -> 1,
    c8[4] -> 0,
    c9[4] -> 0,
    c10[4] -> 0,
    c11[4] -> 0,
    c12[4] -> 0,
    c13[4] -> -3 g2^2/2,
    c14[4] -> 3,
    c15[4] -> 3 g2^2/2,
    c16[4] -> 0,
    c17[4] -> -3,
    c18[4] -> -3,
    c19[4] -> -3 g2^2/2,
    c20[4] -> 3,
    c21[4] -> 3 g2^2/2,
    c22[4] -> 0,
    c23[4] -> -3,
    c24[4] -> 0,
    c25[4] -> 0,
    c26[4] -> 0,
    c27[4] -> 0,

    (* Table 5 *)
    d1[1][i_, j_, k_, l_] :> (-Te[k, i] Te[l, j] Conjugate[Te[k, j]] Conjugate[Te[l, i]]),
    d1[2][i_, j_, k_, l_] :> (-Abs[Mu]^4 Ye[k, j] Ye[l, i] Conjugate[Ye[k, i]] Conjugate[Ye[l, j]]),
    d1[3][i_, j_, k_, l_] :> (-Abs[Mu]^2 Te[l, i] Ye[k, j] (Conjugate[Te[l, j]] Conjugate[Ye[k, i]] +
                                                            Conjugate[Te[k, i]] Conjugate[Ye[l, j]])),
    d1[4][i_, j_, k_, l_] :> (Abs[Mu]^2 Te[l, i] Conjugate[Te[k, i]] Ye[k, j] Conjugate[Ye[l, j]]),
    d1[5][i_, j_, k_, l_] :> (-Mu^2 Te[k, j] Te[l, i] Conjugate[Ye[k, i]] Conjugate[Ye[l, j]]),
    d1[6][i_, j_, k_, l_] :> (Mu Te[k, i] Te[l, j] Conjugate[Te[l, i]] Conjugate[Ye[k, j]]),
    d1[7][i_, j_, k_, l_] :> (Mu Abs[Mu]^2 Te[l, j] Ye[k, i] Conjugate[Ye[k, j]] Conjugate[Ye[l, i]]),

    (* Table 6 *)
    d2[1][i_, j_, k_, l_] :> (-3 Td[k, i] Td[l, j] Conjugate[Td[k, j]] Conjugate[Td[l, i]]),
    d2[2][i_, j_, k_, l_] :> (-3 Abs[Mu]^4 Yd[k, j] Yd[l, i] Conjugate[Yd[k, i]] Conjugate[Yd[l, j]]),
    d2[3][i_, j_, k_, l_] :> (-3 Abs[Mu]^2 Td[l, i] Yd[k, j] (Conjugate[Td[l, j]] Conjugate[Yd[k, i]] +
                                                              Conjugate[Td[k, i]] Conjugate[Yd[l, j]])),
    d2[4][i_, j_, k_, l_] :> (3 Abs[Mu]^2 Td[l, i] Conjugate[Td[k, i]] Yd[k, j] Conjugate[Yd[l, j]]),
    d2[5][i_, j_, k_, l_] :> (-3 Mu^2 Td[k, j] Td[l, i] Conjugate[Yd[k, i]] Conjugate[Yd[l, j]]),
    d2[6][i_, j_, k_, l_] :> (3 Mu Td[k, i] Td[l, j] Conjugate[Td[l, i]] Conjugate[Yd[k, j]]),
    d2[7][i_, j_, k_, l_] :> (3 Mu Abs[Mu]^2 Td[l, j] Yd[k, i] Conjugate[Yd[k, j]] Conjugate[Yd[l, i]]),

    d3[1][i_, j_, k_, l_] :> (-3 Abs[Mu]^4 Yu[i, l] Yu[j, k] Conjugate[Yu[i, k]] Conjugate[Yu[j, l]]),
    d3[2][i_, j_, k_, l_] :> (-3 Tu[i, k] Tu[j, l] Conjugate[Tu[i, l]] Conjugate[Tu[j, k]]),
    d3[3][i_, j_, k_, l_] :> (-3 Abs[Mu]^2 Tu[j, l] Yu[i, k] (Conjugate[Tu[j, k]] Conjugate[Yu[i, l]] +
                                                              Conjugate[Tu[i, l]] Conjugate[Yu[j, k]])),
    d3[4][i_, j_, k_, l_] :> (3 Abs[Mu]^2 Tu[j, l] Conjugate[Tu[i, l]] Yu[i, k] Conjugate[Yu[j, k]]),
    d3[5][i_, j_, k_, l_] :> (-3 Mu^2 Tu[i, k] Tu[j, l] Conjugate[Yu[i, l]] Conjugate[Yu[j, k]]),
    d3[6][i_, j_, k_, l_] :> (3 Mu Abs[Mu]^2 Tu[i, l] Yu[j, k] Conjugate[Yu[i, k]] Conjugate[Yu[j, l]]),
    d3[7][i_, j_, k_, l_] :> (3 Mu Tu[i, l] Tu[j, k] Conjugate[Tu[j, l]] Conjugate[Yu[i, k]]),

    (* Eq. (124) *)
    d4[1][__] :> 0,
    d4[2][__] :> 0,
    d4[3][__] :> 0,
    d4[5][__] :> 0,
    d4[6][__] :> 0,
    d4[7][__] :> 0,
    d4[4][i_, j_, k_, l_] :> (-3 (Td[k, i] Conjugate[Tu[k, l]] -
                                  Abs[Mu]^2 Yd[k, i] Conjugate[Yu[k, l]]) (Tu[j, l] Conjugate[
                                      Td[j, i]] - Abs[Mu]^2 Yu[j, l] Conjugate[Yd[j, i]]))
};

(* loop functions, Eq. (130)-(131) *)
GetTHDMThresholds1LLoopFunctions[] := {
    B0[m1_, m2_, mu_] :> If[PossibleZeroQ[m1 - m2],
       Log[mu^2/m2^2],
       1 + (m1^2 Log[mu^2/m1^2] - m2^2 Log[mu^2/m2^2])/(m1^2 - m2^2)
       ],

    DB0[m1_, m2_] :> If[PossibleZeroQ[m1 - m2],
       -1/(6*m2^2),
       -(m1^4 - m2^4 + 2 m1^2 m2^2 Log[m2^2/m1^2])/(2 (m1^2 - m2^2)^3)
       ],

    C0[m1_, m2_, m3_] :> Which[
       PossibleZeroQ[m1 - m2] && PossibleZeroQ[m1 - m3],
       -1/(2*m3^2),
       PossibleZeroQ[m1 - m2],
       (-m2^2 + m3^2 + 2*m3^2*Log[m2/m3])/(m2^2 - m3^2)^2,
       PossibleZeroQ[m1 - m3],
       (m2^2 - m3^2 - 2*m2^2*Log[m2/m3])/(m2^2 - m3^2)^2,
       PossibleZeroQ[m2 - m3],
       (m1^2 - m3^2 + m1^2*Log[m3^2/m1^2])/(m1^2 - m3^2)^2,
       (* general case *)
       True,
       (m1^2 m2^2 Log[m2^2/m1^2] + m3^2 m2^2 Log[m3^2/m2^2] +
          m1^2 m3^2 Log[m1^2/m3^2])/((m1^2 - m2^2) (m1^2 - m3^2) (m2^2 -
            m3^2))
       ],

    D0[m1_, m2_, m3_, m4_] :> Which[
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
       (C0[m1, m3, m4] - C0[m2, m3, m4])/(m1^2 - m2^2)],

    D2tilde[m1_, m2_, m3_, m4_] :>
      C0[m2, m3, m4] + m1^2 D0[m1, m2, m3, m4],

    D4tilde[m1_, m2_, m3_, m4_, mu_] :> (
       B0[m3, m4, mu] + (m1^2 + m2^2) C0[m2, m3, m4]
        + m1^4 D0[m1, m2, m3, m4]
       ),

    W[m1_, m2_, mu_] :> If[PossibleZeroQ[m1 - m2],
       -2/3 + 2 Log[mu^2/m2^2],
       (+ 2 Log[mu^2/m1^2]
        + Log[m2^2/m1^2] (2 m2^6 - 6 m1^2 m2^4)/(m1^2 - m2^2)^3
        + (m1^4 - 6 m2^2 m1^2 + m2^4)/(m1^2 - m2^2)^2)
       ]
};

(* counter-terms *)

dgY := -dZB/2 flagdg;
dg2 := -dZW/2 flagdg;
dZdd := flagSferm flagZdd dZddSferm + flagIno flagZdd dZddIno;
dZud := flagSferm flagZud dZudSferm + flagIno flagZud dZudIno;
dZuu := flagSferm flagZuu dZuuSferm + flagIno flagZuu dZuuIno;
(* Eq. (117) *)
dZW := Module[{i},
    g2^2 kappa/6 (
    - 4 flagMSDRg2 (* MS-bar/DR-bar conversion term *)
    + flagIno (4 Log[Abs[Mu]^2/Q^2] + 8 Log[M2^2/Q^2])
    + flagSferm (
        Summation[
            Log[msl[i]^2/Q^2] + Nc Log[msq[i]^2/Q^2], {i, 1, 3}]))];
(* Eq. (117) *)
dZB :=  Module[{i},
    gY^2 kappa/3 (
    flagIno 2 Log[Abs[Mu]^2/Q^2]
    + flagSferm (
        Summation[
            Log[mse[i]^2/Q^2] + Log[msl[i]^2/Q^2]/2 +
            4 Nc/9 Log[msu[i]^2/Q^2] + Nc/9 Log[msd[i]^2/Q^2] +
            Nc/18 Log[msq[i]^2/Q^2], {i, 1, 3}]))];
(* Eq. (118) *)
dZddSferm := Module[{i,j},
    kappa/2 Summation[
    3 DB0[msd[i], msq[j]] Td[j, i] Conjugate[Td[j, i]]
     + 3 Abs[Mu]^2 DB0[msu[i], msq[j]] Yu[j, i] Conjugate[Yu[j, i]]
     + DB0[mse[i], msl[j]] Te[j, i] Conjugate[Te[j, i]]
    , {j, 1, 3}, {i, 1, 3}]];
(* Eq. (119) *)
dZddIno := -kappa/8 (gY^2 W[Abs[M1], Abs[Mu], Q] +
     3 g2^2 W[Abs[M2], Abs[Mu], Q]);
(* Eq. (118) *)
dZudSferm :=  Module[{i,j},
    -kappa Summation[
    3 Conjugate[Mu] DB0[msd[i], msq[j]] Yd[j, i] Conjugate[Td[j, i]]
     + 3 Conjugate[Mu] DB0[msu[i], msq[j]] Yu[j, i] Conjugate[Tu[j, i]]
     + Conjugate[Mu] DB0[mse[i], msl[j]] Ye[j, i] Conjugate[Te[j, i]]
    , {j, 1, 3}, {i, 1, 3}]];
(* Eq. (119) *)
dZudIno := -kappa Conjugate[
    Mu] (gY^2 Conjugate[M1] DB0[Abs[M1], Abs[Mu]] +
     3 g2^2 Conjugate[M2] DB0[Abs[M2], Abs[Mu]]);
(* Eq. (118) *)
dZuuSferm :=  Module[{i,j},
    kappa/2 Summation[
    3 DB0[msu[i], msq[j]] Tu[j, i] Conjugate[Tu[j, i]]
     + 3 Abs[Mu]^2 DB0[msd[i], msq[j]] Yd[j, i] Conjugate[Yd[j, i]]
     + Abs[Mu]^2 DB0[mse[i], msl[j]] Ye[j, i] Conjugate[Ye[j, i]]
    , {j, 1, 3}, {i, 1, 3}]];
(* Eq. (119) *)
dZuuIno := dZddIno;

(* tree-level couplings, Eq. (21) *)
lamTree[[1]] = lamTree[[2]] = (g2^2 + gY^2)/4;
lamTree[[3]] = -lamTree[[1]];
lamTree[[4]] = g2^2/2;
lamTree[[5]] = lamTree[[6]] = lamTree[[7]] = 0;

(* Higgsino + gaugino contribution, Eq. (120) *)
lamIno[[5]] = (
   3 g2^4 Mu^2 M2^2 D0[M2, M2, Abs[Mu], Abs[Mu]]
    + 2 g2^2 gY^2 Mu^2 M1 M2 D0[M1, M2, Abs[Mu], Abs[Mu]]
    + gY^4 Mu^2 M1^2 D0[M1, M1, Abs[Mu], Abs[Mu]]
   );

(* Eq. (121) *)
lamIno123467[i_] := (
   g2^4 (flagMSDRlam as[i] + a2[i] D2tilde[M2, M2, Abs[Mu], Abs[Mu]] + 
       a4[i] D4tilde[M2, M2, Abs[Mu], Abs[Mu], Q])
    + g2^2 gY^2 (flagMSDRlam asp[i] + a2p[i] D2tilde[M1, M2, Abs[Mu], Abs[Mu]] + 
       a4p[i] D4tilde[M1, M2, Abs[Mu], Abs[Mu], Q])
    + gY^4 (flagMSDRlam aspp[i] + a2pp[i] D2tilde[M1, M1, Abs[Mu], Abs[Mu]] + 
       a4pp[i] D4tilde[M1, M1, Abs[Mu], Abs[Mu], Q])
   );

lamIno[[1]] = lamIno123467[1];
lamIno[[2]] = lamIno123467[2];
lamIno[[3]] = lamIno123467[3];
lamIno[[4]] = lamIno123467[4];
lamIno[[6]] = lamIno123467[6];
lamIno[[7]] = lamIno123467[7];

(* Eq. (127) *)
Yee[i_, j_]    := Module[{l1}, Summation[Conjugate[Ye[l1, i]] Ye[l1, j], {l1, 1, 3}]];
Yeebar[i_, j_] := Module[{l2}, Summation[Ye[i, l2] Conjugate[Ye[j, l2]], {l2, 1, 3}]];
Yuu[i_, j_]    := Module[{l3}, Summation[Conjugate[Yu[l3, i]] Yu[l3, j], {l3, 1, 3}]];
Yuubar[i_, j_] := Module[{l4}, Summation[Yu[i, l4] Conjugate[Yu[j, l4]], {l4, 1, 3}]];
Ydd[i_, j_]    := Module[{l5}, Summation[Conjugate[Yd[l5, i]] Yd[l5, j], {l5, 1, 3}]];
Yddbar[i_, j_] := Module[{l6}, Summation[Yd[i, l6] Conjugate[Yd[j, l6]], {l6, 1, 3}]];
Yud[i_, j_]    := Module[{l7}, Summation[Conjugate[Yu[l7, i]] Yd[l7, j], {l7, 1, 3}]];
Ydu[i_, j_]    := Module[{l8}, Summation[Conjugate[Yd[l8, i]] Yu[l8, j], {l8, 1, 3}]];

(* sfermion contribution *)

(* Eq. (125) *)
lamSlep1234[l_] := Module[{i,j,k},
    Summation[
   (b1[l] KroneckerDelta[i, j]
       + b2[l] Yee[i, i] KroneckerDelta[i, j]
       + b3[l] Yee[i, j] Yee[j, i]
      ) B0[mse[i], mse[j], Q]
    + (b4[l] KroneckerDelta[i, j]
       + b5[l] Yeebar[i, i] KroneckerDelta[i, j]
       + b6[l] Yeebar[i, j] Yeebar[j, i]
      ) B0[msl[i], msl[j], Q]
    + Summation[
     (c1[l] Abs[Mu]^2 Ye[k, i] Conjugate[Ye[k, i]] KroneckerDelta[i, j]
         + c2[l] Te[k, i] Conjugate[Te[k, i]] KroneckerDelta[i, j]
         + c3[l] Abs[Mu]^2 Ye[k, i] Conjugate[Ye[k, j]] Yee[i, j]
         + c4[l] Te[k, i] Conjugate[Te[k, j]] Yeebar[i, j]
        ) C0[mse[i], mse[j], msl[k]]
      + (c5[l] Abs[Mu]^2 Ye[j, i] Conjugate[Ye[j, i]] KroneckerDelta[
           j, k]
         + c6[l] Te[j, i] Conjugate[Te[j, i]] KroneckerDelta[j, k]
         + c7[l] Abs[Mu]^2 Ye[j, i] Conjugate[Ye[k, i]] Yeebar[k, j]
         + c8[l] Te[j, i] Conjugate[Te[k, i]] Yeebar[k, j]
        ) C0[mse[i], msl[j], msl[k]]
     , {k, 1, 3}], {i, 1, 3}, {j, 1, 3}]];

(* Eq. (126) of lambda_{1,2,3,4}(n) *)
lamSquark1234[n_] := Module[{i,j,k,l},
    Summation[
   (b7[n] KroneckerDelta[i, j]
       + b8[n] Ydd[i, i] KroneckerDelta[i, j]
       + b9[n] Ydd[i, j] Ydd[j, i]
      ) B0[msd[i], msd[j], Q]
    + b10[n] Ydu[i, j] Yud[j, i] B0[msd[i], msu[j], Q]
    + (b11[n] KroneckerDelta[i, j]
       + b12[n] Yuu[i, i] KroneckerDelta[i, j]
       + b13[n] Yuu[i, j] Yuu[j, i]
      ) B0[msu[i], msu[j], Q]
    + (b14[n] KroneckerDelta[i, j]
       + b15[n] Yddbar[i, i] KroneckerDelta[i, j]
       + b16[n] Yuubar[i, i] KroneckerDelta[i, j]
       + b17[n] Yddbar[i, j] Yddbar[j, i]
       + b18[n] Yuubar[i, j] Yuubar[j, i]
       + b19[n] Yddbar[i, j] Yuubar[j, i]
      ) B0[msq[i], msq[j], Q]
    + Summation[(c9[n] Abs[Mu]^2 Yd[k, i] Conjugate[
           Yd[k, i]] KroneckerDelta[i, j]
         + c10[n] Td[k, i] Conjugate[Td[k, i]] KroneckerDelta[i, j]
         + c11[n] Abs[Mu]^2 Yd[k, i] Conjugate[Yd[k, j]] Ydd[i, j]
         + c12[n] Td[k, i] Conjugate[Td[k, j]] Ydd[i, j]
        ) C0[msd[i], msd[j], msq[k]]
      + (c13[n] Abs[Mu]^2 Yd[j, i] Conjugate[Yd[j, i]] KroneckerDelta[
           j, k]
         + c14[n] Abs[Mu]^2 Yd[j, i] Conjugate[Yd[k, i]] Yddbar[k, j]
         + c15[n] Td[j, i] Conjugate[Td[j, i]] KroneckerDelta[j, k]
         + c16[n] Td[j, i] Conjugate[Td[k, i]] Yddbar[k, j]
         + c17[n] Td[j, i] Conjugate[Td[k, i]] Yuubar[k, j]
        ) C0[msd[i], msq[j], msq[k]]
      + (-c18[n] Yd[j, i] Conjugate[Yu[j, k]] Ydu[i, k] Abs[Mu]^2
         - c18[n] Yu[j, k] Conjugate[Yd[j, i]] Yud[k, i] Abs[Mu]^2
         + c18[n] Tu[j, k] Conjugate[Td[j, i]] Yud[k, i]
         + c18[n] Td[j, i] Conjugate[Tu[j, k]] Ydu[i, k]
        ) C0[msd[i], msq[j], msu[k]]
      + (c19[n] Abs[Mu]^2 Yu[i, k] Conjugate[Yu[i, k]] KroneckerDelta[
           i, j]
         + c20[n] Abs[Mu]^2 Yu[j, k] Conjugate[Yu[i, k]] Yuubar[i, j]
         + c21[n] Tu[i, k] Conjugate[Tu[i, k]] KroneckerDelta[i, j]
         + c22[n] Tu[i, k] Conjugate[Tu[j, k]] Yuubar[j, i]
         + c23[n] Tu[i, k] Conjugate[Tu[j, k]] Yddbar[j, i]
        ) C0[msq[i], msq[j], msu[k]]
      + (c24[n] Abs[Mu]^2 Yu[i, j] Conjugate[Yu[i, j]] KroneckerDelta[
           j, k]
         + c25[n] Abs[Mu]^2 Yu[i, k] Conjugate[Yu[i, j]] Yuu[k, j]
         + c26[n] Tu[i, j] Conjugate[Tu[i, j]] KroneckerDelta[j, k]
         + c27[n] Tu[i, j] Conjugate[Tu[i, k]] Yuu[j, k]
        ) C0[msq[i], msu[j], msu[k]]
      + Summation[
       d1[n][i, j, k, l] D0[mse[i], mse[j], msl[k], msl[l]]
        + d2[n][i, j, k, l] D0[msd[i], msd[j], msq[k], msq[l]]
        + d3[n][i, j, k, l] D0[msq[i], msq[j], msu[k], msu[l]]
        + d4[n][i, j, k, l] D0[msd[i], msq[j], msq[k], msu[l]]
       , {l, 1, 3}]
     , {k, 1, 3}]
   , {i, 1, 3}, {j, 1, 3}]];

(* Eq. (128) for lambda_{6,7}(n) *)
lamSferm67[n_] := Module[{i,j,k,l},
    Summation[
   (c1p[n] Mu Te[k, i] Conjugate[Ye[k, i]] KroneckerDelta[i, j] + 
       c2p[n] Mu Te[k, i] Conjugate[Ye[k, j]] Yee[i, j]) C0[mse[i], 
      mse[j], msl[k]]
    + (c3p[n] Mu Te[j, i] Conjugate[Ye[j, i]] KroneckerDelta[j, k] + 
       c4p[n] Mu Te[j, i] Conjugate[Ye[k, i]] Yeebar[k, j]) C0[mse[i],
       msl[j], msl[k]]
    + (c5p[n] Mu Td[k, i] Conjugate[Yd[k, i]] KroneckerDelta[i, j] + 
       c6p[n] Mu Td[k, i] Conjugate[Yd[k, j]] Ydd[i, j]) C0[msd[i], 
      msd[j], msq[k]]
    + (c7p[n] Mu Td[j, i] Conjugate[Yd[j, i]] KroneckerDelta[j, k] + 
       c8p[n] Mu Td[j, i] Conjugate[Yd[k, i]] Yddbar[k, j]) C0[msd[i],
       msq[j], msq[k]]
    + (c9p[n] Mu Tu[i, k] Conjugate[Yu[i, k]] KroneckerDelta[i, j] + 
       c10p[n] Mu Tu[i, k] Conjugate[Yu[j, k]] Yuubar[j, i]) C0[
      msq[i], msq[j], msu[k]]
    + (c11p[n] Mu Tu[i, j] Conjugate[Yu[i, j]] KroneckerDelta[j, k] + 
       c12p[n] Mu Tu[i, j] Conjugate[Yu[i, k]] Yuu[j, k]) C0[msq[i], 
      msu[j], msu[k]]
    + Summation[
     +d1[n][i, j, k, l] D0[mse[i], mse[j], msl[k], msl[l]]
      + d2[n][i, j, k, l] D0[msd[i], msd[j], msq[k], msq[l]]
      + d3[n][i, j, k, l] D0[msq[i], msq[j], msq[k], msq[l]]
     , {l, 1, 3}], {j, 1, 3}, {k, 1, 3}, {i, 1, 3}]];

lamSquark = {Undef, Undef, Undef, Undef};
lamSquark[[1]] = lamSquark1234[1];
lamSquark[[2]] = lamSquark1234[2];
lamSquark[[3]] = lamSquark1234[3];
lamSquark[[4]] = lamSquark1234[4];

(* above Eq. (125) *)
lamSferm[[1]] = lamSlep1234[1] + lamSquark[[1]];
lamSferm[[2]] = lamSlep1234[2] + lamSquark[[2]];
lamSferm[[3]] = lamSlep1234[3] + lamSquark[[3]];
lamSferm[[4]] = lamSlep1234[4] + lamSquark[[4]];
(* Eq. (123) *)
lamSferm[[5]] = Module[{i,j,k,l},
    Summation[
   d1[5][i, j, k, l] D0[mse[i], mse[j], msl[k], msl[l]]
    + d2[5][i, j, k, l] D0[msd[i], msd[j], msq[k], msq[l]]
    + d3[5][i, j, k, l] D0[msq[i], msq[j], msu[k], msu[l]]
   , {i, 1, 3}, {j, 1, 3}, {k, 1, 3}, {l, 1, 3}]];
lamSferm[[6]] = lamSferm67[6];
lamSferm[[7]] = lamSferm67[7];

(* Eq. (116) *)
lamHat := loopOrder0 lamTree +
    loopOrder1 kappa (flagIno lamIno + flagSferm lamSferm);

(* Eq. (71) *)
lamBar[[1]] := lamHat[[1]] +
    loopOrder1 ((g2^2 + gY^2) Re[dZdd] + (g2^2 dg2 + gY^2 dgY)/2);

lamBar[[2]] := lamHat[[2]] +
   loopOrder1 ((g2^2 + gY^2) Re[dZuu] + (g2^2 dg2 + gY^2 dgY)/2);

lamBar[[3]] := lamHat[[3]] +
   loopOrder1 (-(g2^2 + gY^2)/2 (Re[dZdd] + Re[dZuu]) - (g2^2 dg2 + gY^2 dgY)/2);

lamBar[[4]] = lamHat[[4]] +
   loopOrder1 (g2^2 (Re[dZdd] + Re[dZuu]) + g2^2 dg2);

lamBar[[5]] = lamHat[[5]];

lamBar[[6]] = lamHat[[6]] +
    loopOrder1 (-(g2^2 + gY^2)/4 Conjugate[dZud]);

lamBar[[7]] = lamHat[[7]] +
    loopOrder1 ((g2^2 + gY^2)/4 Conjugate[dZud]);

Options[GetTHDMThresholds1L] = {
    coefficients -> GetTHDMThresholds1LCoefficients[],
    flags -> GetTHDMThresholds1LFlags[],
    loopFunctions -> {},
    loopOrder -> {1,1},
    sumHead -> Sum
};

(* relation to Haber/Wagner/Lee convention p. 6 *)
GetTHDMThresholds1L[OptionsPattern[]] := {
    lamBar[[1]],
    lamBar[[2]],
    lamBar[[3]] + lamBar[[4]],
    -lamBar[[4]],
    lamBar[[5]],
    lamBar[[6]],
    lamBar[[7]]
} /.
    loopOrder0 -> OptionValue[loopOrder][[1]] /.
    loopOrder1 -> OptionValue[loopOrder][[2]] /.
    OptionValue[flags] /.
    OptionValue[coefficients] /.
    Summation -> OptionValue[sumHead] //.
    OptionValue[loopFunctions];

End[];
