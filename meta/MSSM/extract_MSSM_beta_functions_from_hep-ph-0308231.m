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

(*
   This script extracts the 3-loop MSSM beta-functions from
   http://www.liv.ac.uk/~dij/betas/allgennb.log .

   Author: Alexander Voigt

   Run it like this:

   math -run "<< extract_MSSM_beta_functions_from_hep-ph-0308231.m; Quit[]"
*)

betaStr = Import["allgennb.log"];

(* fix syntax error *)
betaStr = StringReplace[betaStr, "beta_htau^(3)) =" -> "beta_htau^(3) ="];
    
(* wrap expressions around brackets *)
betaStr = StringReplace[betaStr, {"=" -> "= (", ";" -> ");"}];

(* convert symbols into valid Wolfram symbols *)
betaStr = StringReplace[betaStr,
   {"_" -> "", "^(1)" -> "1", "^(2)" -> "2", "^(3)" -> "3"}];

matrixPattern = (
    (WordBoundary ~~ "Ye" ~~ WordBoundary) | "Yec" |
    (WordBoundary ~~ "Yt" ~~ WordBoundary) | "Ytc" |
    (WordBoundary ~~ "Yb" ~~ WordBoundary) | "Ybc" |
    (WordBoundary ~~ "he" ~~ WordBoundary) | "hec" |
    (WordBoundary ~~ "ht" ~~ WordBoundary) | "htc" |
    (WordBoundary ~~ "hb" ~~ WordBoundary) | "hbc" |
    (WordBoundary ~~ "mq" ~~ WordBoundary) |
    (WordBoundary ~~ "mt" ~~ WordBoundary) |
    (WordBoundary ~~ "mb" ~~ WordBoundary) |
    (WordBoundary ~~ "ml" ~~ WordBoundary) |
    (WordBoundary ~~ "me" ~~ WordBoundary));

(* convert matrix products to non-commutative products *)
betaStr = StringReplace[betaStr,
   {"tr(" ~~ p:(Except[")"]..) ~~ ")" :>
    "trace[" <> StringReplace[p, "*" -> ","] <> "]",
    p:(matrixPattern ~~ ("*" ~~ matrixPattern)..) :>
    "MatMul[" <> StringReplace[p, "*" -> ","] <> "]"}];

(* convert adjoint symbols to SARAH convention *)
betaStr = StringReplace[betaStr,
   {"Yec" -> "Adj[Ye]",
    "Ytc" -> "Adj[Yt]",
    "Ybc" -> "Adj[Yb]",
    "hec" -> "Adj[he]",
    "htc" -> "Adj[ht]",
    "hbc" -> "Adj[hb]" }];

ReleaseHold[ToExpression[betaStr, InputForm, Hold] /. Times :> MatMul];

scalarPat = g1 | g2 | g3 | M1 | M2 | M3 | mh1 | mh2 | trace[__];

repl = {
    a1 -> g1^2,
    a2 -> g2^2,
    a3 -> g3^2,
    k -> 6 Zeta[3],
    n5 -> 0,
    n10 -> 0,
    MatMul[] -> 1,
    MatMul[a___, 0, c___] -> 0,
    MatMul[a_] :> a,
    MatMul[a___, f:scalarPat, c___] :> f MatMul[a,c],
    MatMul[a___, f:Power[scalarPat,_], c___] :> f MatMul[a,c],
    MatMul[a___, b_?NumberQ, c___] :> b MatMul[a,c],
    MatMul[a___, Plus[b1_,b2_], c___] :> MatMul[a,b1,c] + MatMul[a,b2,c],
    MatMul[a___, HoldPattern[Times[b__]], c___] :> MatMul[a,b,c]
};

betaMu1 = Mu * (gammaH11 + gammaH21);
betaMu2 = Mu * (gammaH12 + gammaH22);
betaMu3 = Mu * (gammaH13 + gammaH23);

(* calculates gamma_1, Eq. (2.4b) arxiv:hep-ph/0408128 *)
Gamma1[gamma_] :=
    Module[{gp = gamma //. repl},
           DYuk[yuk_, tri_] := {
               Derivative[1, 0][trace][yuk, p_] :> trace[tri, p],
               Derivative[0, 1][trace][p_, yuk] :> trace[p, tri],
               Derivative[1, 0][trace][Adj[yuk], p_] -> 0,
               Derivative[0, 1][trace][p_, Adj[yuk]] -> 0,
               Derivative[1, 0, 0, 0][trace][yuk, p__    ] :> trace[tri, p],
               Derivative[0, 1, 0, 0][trace][p_, yuk, q__] :> trace[p, tri, q],
               Derivative[0, 0, 1, 0][trace][p__, yuk, q_] :> trace[p, tri, q],
               Derivative[0, 0, 0, 1][trace][p__, yuk    ] :> trace[p, tri],
               Derivative[1, 0, 0, 0][trace][Adj[yuk], p__    ] -> 0,
               Derivative[0, 1, 0, 0][trace][p_, Adj[yuk], q__] -> 0,
               Derivative[0, 0, 1, 0][trace][p__, Adj[yuk], q_] -> 0,
               Derivative[0, 0, 0, 1][trace][p__, Adj[yuk]    ] -> 0,
               Derivative[1, 0, 0, 0, 0, 0][trace][yuk, p__        ] :> trace[tri, p],
               Derivative[0, 1, 0, 0, 0, 0][trace][p_, yuk, q__    ] :> trace[p, tri, q],
               Derivative[0, 0, 1, 0, 0, 0][trace][p_, q_, yuk, r__] :> trace[p, q, tri, r],
               Derivative[0, 0, 0, 1, 0, 0][trace][p__, yuk, q_, r_] :> trace[p, tri, q, r],
               Derivative[0, 0, 0, 0, 1, 0][trace][p__, yuk, q_    ] :> trace[p, tri, q],
               Derivative[0, 0, 0, 0, 0, 1][trace][p__, yuk        ] :> trace[p, tri],
               Derivative[1, 0, 0, 0, 0, 0][trace][Adj[yuk], p__        ] -> 0,
               Derivative[0, 1, 0, 0, 0, 0][trace][p_, Adj[yuk], q__    ] -> 0,
               Derivative[0, 0, 1, 0, 0, 0][trace][p_, q_, Adj[yuk], r__] -> 0,
               Derivative[0, 0, 0, 1, 0, 0][trace][p__, Adj[yuk], q_, r_] -> 0,
               Derivative[0, 0, 0, 0, 1, 0][trace][p__, Adj[yuk], q_    ] -> 0,
               Derivative[0, 0, 0, 0, 0, 1][trace][p__, Adj[yuk]        ] -> 0
           };
           (
               + M1 g1/2 D[gp, g1]
               + M2 g2/2 D[gp, g2]
               + M3 g3/2 D[gp, g3]
               - (D[gp, Yt] /. DYuk[Yt, ht])
               - (D[gp, Yb] /. DYuk[Yb, hb])
               - (D[gp, Ye] /. DYuk[Ye, he])
           )
          ];

gamma1H11 = Gamma1[gammaH11];
gamma1H21 = Gamma1[gammaH21];
gamma1H12 = Gamma1[gammaH12];
gamma1H22 = Gamma1[gammaH22];
gamma1H13 = Gamma1[gammaH13];
gamma1H23 = Gamma1[gammaH23];

On[Assert];
Assert[FreeQ[gamma1H13, Derivative[__][__][__]]];
Assert[FreeQ[gamma1H23, Derivative[__][__][__]]];

(* calculate beta-functions of BMu, Eq. (2.3), arxiv:hep-ph/0408128 *)
betaBMu1 = Simplify[BMu (gammaH11 + gammaH21) - 2 Mu (gamma1H11 + gamma1H21)];
betaBMu2 = Simplify[BMu (gammaH12 + gammaH22) - 2 Mu (gamma1H12 + gamma1H22)];
betaBMu3 = Simplify[BMu (gammaH13 + gammaH23) - 2 Mu (gamma1H13 + gamma1H23)];

{g1 bg11  , g1 bg12  , g1 bg13  } //. repl >> "beta_g1.m";
{g2 bg21  , g2 bg22  , g2 bg23  } //. repl >> "beta_g2.m";
{g3 bg31  , g3 bg32  , g3 bg33  } //. repl >> "beta_g3.m";
{betat1   , betat2   , betat3   } //. repl >> "beta_Yu.m";
{betab1   , betab2   , betab3   } //. repl >> "beta_Yd.m";
{betae1   , betae2   , betae3   } //. repl >> "beta_Ye.m";
{betamq1  , betamq2  , betamq3  } //. repl >> "beta_mq2.m";
{betamt1  , betamt2  , betamt3  } //. repl >> "beta_mu2.m";
{betamb1  , betamb2  , betamb3  } //. repl >> "beta_md2.m";
{betamL1  , betamL2  , betamL3  } //. repl >> "beta_ml2.m";
{betamtau1, betamtau2, betamtau3} //. repl >> "beta_me2.m";
{betamH21 , betamH22 , betamH23 } //. repl >> "beta_mHu2.m";
{betamH11 , betamH12 , betamH13 } //. repl >> "beta_mHd2.m";
{betaht1  , betaht2  , betaht3  } //. repl >> "beta_TYu.m";
{betahb1  , betahb2  , betahb3  } //. repl >> "beta_TYd.m";
{betahtau1, betahtau2, betahtau3} //. repl >> "beta_TYe.m";
{betaM11  , betaM12  , betaM13  } //. repl >> "beta_M1.m";
{betaM21  , betaM22  , betaM23  } //. repl >> "beta_M2.m";
{betaM31  , betaM32  , betaM33  } //. repl >> "beta_M3.m";
{betaMu1  , betaMu2  , betaMu3  } //. repl >> "beta_Mu.m";
{betaBMu1 , betaBMu2 , betaBMu3 } //. repl >> "beta_BMu.m";

{gammaQ1  , gammaQ2  , gammaQ3  } //. repl >> "gamma_SqL.m";
{gammat1  , gammat2  , gammat3  } //. repl >> "gamma_SuR.m";
{gammab1  , gammab2  , gammab3  } //. repl >> "gamma_SdR.m";
{gammaL1  , gammaL2  , gammaL3  } //. repl >> "gamma_SlL.m";
{gammatau1, gammatau2, gammatau3} //. repl >> "gamma_SeR.m";
{gammaH11 , gammaH12 , gammaH13 } //. repl >> "gamma_SHd.m";
{gammaH21 , gammaH22 , gammaH23 } //. repl >> "gamma_SHu.m";
