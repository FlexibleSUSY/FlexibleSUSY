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

Get[FileNameJoin[{smDir, "mt_4loop_qcd_and_non-qcd.m"}]];

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

betag3 = (-k 7*g3^3 +k^2 g3^3*((11*g1^2)/10+(9*g2^2)/2-2*(13*g3^2+gb^2+gt^2))+ k^3 (g3^3*(-523*g1^4+g1^2*(-9*g2^2+616*g3^2-303*gt^2)+15*(109*g2^4+3*g2^2*(56*g3^2-31*gt^2)+20*(13*g3^4-16*g3^2*gt^2+6*gt^4))))/120+k^4 (-(1/36) g3^3 (4 g3^4 gt^2 (6709-2448 Zeta[3])+36 g3^2 gt^4 (-427+96 Zeta[3])+2 g3^6 (-63559+89896 Zeta[3])+27 (20 gt^4 \[Lambda]-12 gt^2 \[Lambda]^2+gt^6 (141+8 Zeta[3]))))) /. simp;

betayt = (k gt*((-17*g1^2)/20-(9*g2^2)/4-8*g3^2+(3*gb^2)/2+(9*gt^2)/2+g\[Tau]^2)-k^2(gt*(-2374*g1^4+5*g1^2*(108*g2^2-304*g3^2-3*(7*gb^2+393*gt^2+150*g\[Tau]^2))+75*(92*g2^4-3*g2^2*(48*g3^2+33*gb^2+75*gt^2+10*g\[Tau]^2)+4*(432*g3^4+gb^4+48*gt^4-16*g3^2*(gb^2+9*gt^2)+9*gt^2*g\[Tau]^2+9*g\[Tau]^4+gb^2*(11*gt^2-5*g\[Tau]^2)+24*gt^2*\[Lambda]-6*\[Lambda]^2))))/1200+k^3(gt*(321980*g1^6+3396580*g2^6-5*g1^4*(17768*g2^2+89276*g3^2+97688*gt^2+5445*\[Lambda])-5*g2^4*(84288*g3^2-67960*gt^2+21375*\[Lambda])-10*g1^2*(9486*g2^4+30192*g3^4-36148*g3^2*gt^2+60925*gt^4+g2^2*(32100*g3^2-69658*gt^2-2925*\[Lambda])+12700*gt^2*\[Lambda]-4500*\[Lambda]^2)+2*(-6193500*g3^6+3637640*g3^4*gt^2+586028*gt^6+990000*gt^4*\[Lambda]+9375*gt^2*\[Lambda]^2-45000*\[Lambda]^3-10000*g3^2*(157*gt^4-8*gt^2*\[Lambda]))+10*g2^2*(147308*g3^4+96740*g3^2*gt^2-1125*(177*gt^4+60*gt^2*\[Lambda]-20*\[Lambda]^2))))/20000+k^4 2*1154.09 g3^8 gt) /. simp;

betalambda = 1/2 (k (27*g1^4+90*g1^2*(g2^2-2*\[Lambda])+25*(9*g2^4-36*g2^2*\[Lambda]-16*(3*gb^4+3*gt^4+g\[Tau]^4-3*gb^2*\[Lambda]-3*gt^2*\[Lambda]-g\[Tau]^2*\[Lambda]-3*\[Lambda]^2)))/100+k^2(-3411*g1^6-15*g1^4*(559*g2^2-60*gb^2+228*gt^2+300*g\[Tau]^2-629*\[Lambda])-25*g1^2*(289*g2^4-6*g2^2*(36*gb^2+84*gt^2+44*g\[Tau]^2+39*\[Lambda])-4*(16*gb^4-32*gt^4-48*g\[Tau]^4+25*gb^2*\[Lambda]+85*gt^2*\[Lambda]+75*g\[Tau]^2*\[Lambda]+108*\[Lambda]^2))+125*(305*g2^6+12*g2^2*\[Lambda]*(15*gb^2+15*gt^2+5*g\[Tau]^2+36*\[Lambda])-g2^4*(36*gb^2+36*gt^2+12*g\[Tau]^2+73*\[Lambda])-8*(-60*gb^6-60*gt^6-20*g\[Tau]^6+3*gt^4*\[Lambda]+g\[Tau]^4*\[Lambda]+72*gt^2*\[Lambda]^2+24*g\[Tau]^2*\[Lambda]^2+78*\[Lambda]^3+3*gb^4*(4*gt^2+\[Lambda])+16*g3^2*(4*gb^4+4*gt^4-5*gb^2*\[Lambda]-5*gt^2*\[Lambda])+6*gb^2*(2*gt^4+7*gt^2*\[Lambda]+12*\[Lambda]^2))))/1000+k^3(-60320*g1^8-4563640*g2^8-40*g1^6*(1543*g2^2-663*g3^2-11117*gt^2-14084*\[Lambda])+20*g2^6*(15072*g3^2+125000*gt^2+865483*\[Lambda])+80*g2^4*(7942*gt^4+6*g3^2*(1372*gt^2-2381*\[Lambda])-79916*gt^2*\[Lambda]-98785*\[Lambda]^2)+2*g1^4*(130000*g2^4+318960*gt^4+10*g3^2*(2032*gt^2-8381*\[Lambda])-748599*gt^2*\[Lambda]-927660*\[Lambda]^2+10*g2^2*(2210*g3^2+21254*gt^2+61753*\[Lambda]))-10*g1^2*(151556*g2^6-135720*gt^6+42030*gt^4*\[Lambda]+63869*gt^2*\[Lambda]^2+38745*\[Lambda]^3-4*g2^4*(1507*g3^2+13041*gt^2+39819*\[Lambda])-4*g3^2*(17570*gt^4+8727*gt^2*\[Lambda])+g2^2*(-45544*g3^2*gt^2+281424*gt^4-11230*gt^2*\[Lambda]+316640*\[Lambda]^2))-5*(1945192*gt^8+893528*gt^6*\[Lambda]-3536520*gt^4*\[Lambda]^2-873000*gt^2*\[Lambda]^3-3005675*\[Lambda]^4+8*g3^4*(50201*gt^4-178484*gt^2*\[Lambda])-4*g3^2*(500988*gt^6-662866*gt^4*\[Lambda]+80385*gt^2*\[Lambda]^2))+2*g2^2*(g3^2*(266980*gt^4+151443*gt^2*\[Lambda])+5*(296552*gt^6-10940*gt^4*\[Lambda]-359539*gt^2*\[Lambda]^2-193726*\[Lambda]^3)))/10000+k^4 2*8308.17*g3^6*gt^4) /. simp;

betav = (k (v*(3*g1^2*(3+Xi)+5*(-4*(3*gb^2+3*gt^2+g\[Tau]^2)+3*g2^2*(3+Xi))))/20+k^2(v*(3*g1^4*(-431+12*Xi+12*Xi^2)+10*g1^2*(9*g2^2*(-3+4*Xi+4*Xi^2)-2*(3*g\[Tau]^2*(25+4*Xi)+gb^2*(25+36*Xi)+gt^2*(85+36*Xi)))+25*(-12*g2^2*(3*gb^2+3*gt^2+g\[Tau]^2)*(5+4*Xi)+g2^4*(271+108*Xi)-8*(80*g3^2*(gb^2+gt^2)-3*(9*gb^4-2*gb^2*gt^2+9*gt^4+3*g\[Tau]^4-2*\[Lambda]^2)))))/800) /. simp

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
