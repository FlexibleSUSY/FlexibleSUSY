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

Needs["TestSuite`", "TestSuite.m"];

smDir = FileNameJoin[{Directory[], "meta", "SM"}];

Get[FileNameJoin[{smDir, "mt_2loop_gaugeless.m"}]];
Get[FileNameJoin[{smDir, "mt_4loop_qcd.m"}]];

L = At/t+1

dQCD = t Normal[Series[MtOvermt^2, {k,0,4}]];

delta1QCD = Coefficient[dQCD, k, 1]

delta2QCD = Coefficient[dQCD, k, 2]

delta3QCD = Coefficient[dQCD, k, 3]

delta4QCD = Coefficient[dQCD, k, 4]

(* squared top pole mass *)
s = (
    t
    + k (delta1QCD + delta1Yukawa)
    + k^2 (delta2QCD + g3^2 delta2mixed + delta2Yukawa)
    + k^3 delta3QCD
    + k^4 delta4QCD
)

(* Define B and A loop integrals in order to evaluate scale dependance analytically *)
ln[x_] := Log[x] - 2q

del[x_,y_,z_] := x^2+y^2+z^2-2x y-2x z-2y z

tau[x_,y_,z_] := (x+y-z+Sqrt[del[x,y,z]])/(2 x)

r[x_,y_,z_] := (x+y-z-Sqrt[del[x,y,z]])/(2 x)

A[x_] := x (ln[x] - 1)

B[x_,y_,s_] := 2 - r[s,x,y] ln[x] - tau[s,y,x] ln[y] + Sqrt[del[s,x,y]]/s Log[tau[x,y,s]]

simpBeta = {
    g1 -> 0,
    g2 -> 0,
    g\[Tau] -> 0,
    gb -> 0,
    gt -> Sqrt[2t]/v,
    \[Lambda] -> h/v^2
}

AddLoops[b_List] :=
    Total @ MapIndexed[#1 k^First[#2]&, b]

(* Load SM beta functions.
   Convention: h = lambda v^2, v ~ 245 GeV *)
betag3 = AddLoops[Get[FileNameJoin[{smDir, "beta_g3.m"}]] /. simpBeta];
betagt = AddLoops[Get[FileNameJoin[{smDir, "beta_gt.m"}]] /. simpBeta];
betal  = AddLoops[Get[FileNameJoin[{smDir, "beta_lambda.m"}]] /. simpBeta];
betav  = AddLoops[Get[FileNameJoin[{smDir, "beta_v.m"}]] /. simpBeta];

betah = Dt[\[Lambda] v^2] /. {
    Dt[\[Lambda]] -> betal,
    Dt[v] -> betav
} /. simpBeta

betat = Dt[gt^2 v^2/2] /. {
    Dt[gt] -> betagt,
    Dt[v] -> betav
} /. simpBeta

(* total derivative *)
ds = Dt[s] /. {
    Dt[k] -> 0,
    Dt[g3] -> betag3,
    Dt[h] -> betah,
    Dt[t] -> betat,
    Dt[v] -> betav,
    Dt[q] -> 1
}

dBhtdh = 1/del[t,h,t] ((h-2t)(Bht-2)+Ah+h-2At-2t)

dBhtdt = h/t 1/del[t,h,t] ((2t-h)(Bht-1)-Ah+2At)-1/t

replDeriv = {
    Dt[Ah] -> -2h+betah (Ah/h+1),
    Dt[At] -> -2t+betat (At/t+1),
    Dt[Bht] -> 2+dBhtdh betah + dBhtdt betat,
    Dt[B00] -> 2-betat/t,
    Dt[Bt0] -> 2,
    Dt[B0h] -> 2,
    Dt[I0h0] -> 2(Ah-h),
    Dt[Ih00] -> 2(Ah-h),
    Dt[I0t0] -> 2(At-t),
    Dt[Ihtt] -> 2(Ah-h)+4(At-t),
    Dt[Ihhh] -> 6(Ah-h),
    Dt[Thht] -> -2 Ah/h,
    Dt[Th00] -> -2 Ah/h,
    Dt[Th0t] -> -2 Ah/h,
    Dt[Tht0] -> -2 Ah/h,
    Dt[Tth0] -> -2 At/t,
    Dt[S0h0] -> t+2(Ah-h),
    Dt[Tbar0ht] -> 2-2Bht,
    Dt[Tbar000] -> 2-2B00,
    Dt[Tbar0t0] -> 2-2Bt0,
    Dt[Tbar00h] -> 2-2B0h,
    Dt[Uthtt] -> 2+2Bht,
    Dt[Uhtht] -> 2+2Bht,
    Dt[Uthhh] -> 2+2Bht,
    Dt[U0000] -> 2+2B00,
    Dt[U000h] -> 2+2B00,
    Dt[U00h0] -> 2+2B00,
    Dt[U000t] -> 2+2B00,
    Dt[U0tht] -> 2+2Bt0,
    Dt[U0tt0] -> 2+2Bt0,
    Dt[Uht00] -> 2+2Bht,
    Dt[Uhtt0] -> 2+2Bht,
    Dt[Ut0h0] -> 2+2Bt0,
    Dt[Uth00] -> 2+2Bht,
    Dt[M0ttht] -> 0,
    Dt[Mhhtth] -> 0,
    Dt[Mhttht] -> 0,
    Dt[M00t00] -> 0,
    Dt[M0t0h0] -> 0,
    Dt[Mh0tt0] -> 0,
    Dt[Mtt00h] -> 0,
    B00 -> (1-At/t)
}

ds = Normal[Series[ds //. replDeriv, {k,0,4}]]

TestEquality[Coefficient[ds,k,1], 0];

TestEquality[Coefficient[ds,k,2], 0];

TestEquality[Coefficient[ds,k^3 g3^6] // Simplify, 0];

TestEquality[Coefficient[ds,k^4 g3^8] // FullSimplify, 0];

PrintTestSummary[];
