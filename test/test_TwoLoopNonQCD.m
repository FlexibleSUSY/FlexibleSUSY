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

(* squared top pole mass *)
s = t + k (g3^2 delta1QCD + delta1Yukawa) + k^2 (g3^4 delta2QCD + g3^2 delta2mixed + delta2Yukawa)

(* Define B and A loop integrals in order to evaluate scale dependance analytically *)
ln[x_] := Log[x] - 2q

del[x_,y_,z_] := x^2+y^2+z^2-2x y-2x z-2y z

tau[x_,y_,z_] := (x+y-z+Sqrt[del[x,y,z]])/(2 x)

r[x_,y_,z_] := (x+y-z-Sqrt[del[x,y,z]])/(2 x)

A[x_] := x (ln[x] - 1)

B[x_,y_,s_] := 2 - r[s,x,y] ln[x] - tau[s,y,x] ln[y] + Sqrt[del[s,x,y]]/s Log[tau[x,y,s]]

simp = {
    g1 -> 0,
    g2 -> 0,
    g\[Tau] -> 0,
    gb -> 0,
    gt -> Sqrt[2t]/v,
    \[Lambda] -> h/v^2
}

(* SM beta functions
   Convention: h = lambda v^2
 *)

betag3 = (
    - 7*g3^3*k
    - 2*g3^3*(13*g3^2 + gt^2)*k^2
    + (5*g3^3*(13*g3^4 - 16*g3^2*gt^2 + 6*gt^4)*k^3)/2
    - (g3^3*k^4*(4*g3^4*gt^2*(6709 - 2448*Zeta[3]) + 
                 36*g3^2*gt^4*(-427 + 96*Zeta[3]) + 2*g3^6*(-63559 + 89896*Zeta[3]) + 
                 27*(20*gt^4*\[Lambda] - 12*gt^2*\[Lambda]^2 + gt^6*(141 + 8*Zeta[3]))))/36
) /. simp

betayt = (
    k*gt*(-8*g3^2 + (9*gt^2)/2)
    - (gt*k^2*(432*g3^4 - 144*g3^2*gt^2 + 48*gt^4 + 24*gt^2*\[Lambda] - 
               6*\[Lambda]^2))/4
    + (gt*k^3*(-6193500*g3^6 + 3637640*g3^4*gt^2 + 
               586028*gt^6 + 990000*gt^4*\[Lambda] + 9375*gt^2*\[Lambda]^2 - 
               45000*\[Lambda]^3 - 10000*g3^2*(157*gt^4 - 8*gt^2*\[Lambda])))/10000
    + k^4 (-(1/81) g3^8 gt (1379027 + 6048 \[Pi]^4 - 2277312 Zeta[3] + 561600 Zeta[5]))
) /. simp

betalambda = (
    + 6*k*(-gt^4 + gt^2*\[Lambda] + \[Lambda]^2)
    + k^2*(30*gt^6 - (3*gt^4*\[Lambda])/2 - 36*gt^2*\[Lambda]^2 - 39*\[Lambda]^3 - 
           8*g3^2*(4*gt^4 - 5*gt^2*\[Lambda]))
    + (k^3*(-1945192*gt^8 - 893528*gt^6*\[Lambda] + 3536520*gt^4*\[Lambda]^2 + 
            873000*gt^2*\[Lambda]^3 + 3005675*\[Lambda]^4 - 
            8*g3^4*(50201*gt^4 - 178484*gt^2*\[Lambda]) + 
            4*g3^2*(500988*gt^6 - 662866*gt^4*\[Lambda] + 80385*gt^2*\[Lambda]^2)))/4000

    + k^4 1/2 16/27 g3^6 gt^4 (13097 + 108 \[Pi]^4 + 44568 Zeta[3] - 47400 Zeta[5])
) /. simp

betav = (
    -3*gt^2*k*v
    - (k^2*v*(80*g3^2*gt^2 - 3*(9*gt^4 - 2*\[Lambda]^2)))/4
) /. simp

betah = 2v^2 betalambda + 2h/v betav

betat = Sqrt[2t] v betayt + 2t/v betav

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

ds = Normal[Series[ds //. replDeriv, {k,0,2}]]

TestEquality[ds, 0];

PrintTestSummary[];
