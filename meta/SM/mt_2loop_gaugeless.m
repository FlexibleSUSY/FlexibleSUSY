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

 2-loop relation between top quark pole and Standard Model MS-bar mass
 from [1604.01134] in the gaugeless limit (g1 = g2 = 0).  Terms of the
 order O(at*as, at^2, lambda^n*at) are included.  Note, that the
 relation is gauge-independent.

 The squared top pole mass s is relatated to the squared MS-bar t as:

 s = t + k (g3^2 delta1QCD + delta1Yukawa) + k^2 (g3^4 delta2QCD + g3^2 delta2mixed + delta2Yukawa)

 We have defined:

 k = 1/(4Pi)^2;
 h = squared MS-bar Higgs mass
 t = gt^2 v^2/2;

 Ah = A[h];
 At = A[t];
 B00 = B[0,0,t] = (1-At/t);
 Bht = B[h,t,t];
 Bt0 = B[t,0,t];
 B0h = B[0,h,t];

 q = Log[Q]

 A[x_] := x (ln[x] - 1)

 ln[x_] := Log[x] - 2q

 B[x,y,s] := Eq.(3.6) [hep-ph/0307101], s = squared momentum

 *)

delta1QCD = t(32/3-8(At/t+1))

delta1Yukawa = (t (Ah-2 At+Bht h-B00 t-4 Bht t))/v^2

delta2QCD = t(2309/9+(16\[Pi]^2)/9+(32\[Pi]^2 Log[2])/9-(16Zeta[3])/3-204(At/t+1)+60(At/t+1)^2)

delta2mixed = 1/(9 h t v^2) (-24 At^2 h (h-t)+24 At t (Ah (3 h+5 t)+h ((5+2 Bht) h+2 (-12+Bht) t))+2 t (-60 Ah t^2-2 h^2 (3 Ihtt+t (-6-51 Bht-8 \[Pi]^2+72 M0ttht t+12 Tbar0ht-24 Th0t))-6 h^3 (1+Bht-4 M0ttht t+Th0t)+h t (-54 Ah+84 Ihtt+t (333-288 Bht-61 \[Pi]^2+24 M00t00 t+192 M0ttht t+96 Tbar0ht-60 Th0t-96 Uthtt))))

delta2Yukawa = 1/(48 v^4) (-72 Bht h^3-12 h^2 I0h0+12 h^2 Ih00+12 h^3 \[Pi]^2-336 h^2 t+528 Bht h^2 t-60 Bht^2 h^2 t+288 h I0h0 t+120 h I0h0 t-144 h Ih00 t-288 h Ihhh t+96 h Ihtt t-96 h Ih00 t-52 h^2 \[Pi]^2 t+96 h^2 q t+96 B0h h^2 q t-336 h S0h0 t+258 h t^2-768 Bht h t^2-72 B00 Bht h t^2+336 Bht^2 h t^2+96 Bt0 h t^2-96 Bht Bt0 h t^2+72 I0h0 t^2+336 I0t0 t^2-960 Ihtt t^2+432 h^2 Mhhtth t^2+48 h^2 Mhttht t^2-96 h^2 Mh0tt0 t^2-48 h^2 Mtt00h t^2+64 h \[Pi]^2 t^2-288 h q t^2-96 B00 h q t^2+96 B0h h q t^2+288 Bt0 h q t^2+72 S0h0 t^2+1131 t^3-444 B00 t^3+96 B00^2 t^3-2112 Bht t^3+48 B00 Bht t^3-384 Bht^2 t^3+(1728 Ihtt t^3)/h-192 h M0t0h0 t^3-1152 h Mhhtth t^3+384 h Mhttht t^3-96 \[Pi]^2 t^3-144 q t^3+144 Bt0 q t^3+(1152 t^4)/h+(2304 Bht t^4)/h-768 Mhttht t^4+(36 Ah^2 (h^2-6 h t-5 t^2))/h+(96 At^2 (h^2-18 h t+42 t^2))/h-(24 At (6 (-1+Bht) h^3-4 h^2 (-7+8 Bht-Bt0+3 q) t+h (-27-10 B00+64 Bht+18 q) t^2+48 (3-2 Bht) t^3))/h-1/h 24 Ah (3 h^3+h^2 (-35+2 B00-4 Bht+Bt0-4 q) t+(-25-4 B00+4 Bht) h t^2+(9+4 B00+32 Bht) t^3+At (7 h^2-8 h t+39 t^2))+96 t^3 Tbar000+48 h^2 t Tbar00h+48 h t^2 Tbar00h+96 h t^2 Tbar0t0+96 t^3 Tbar0t0-12 t^3 Tbar000+48 h t^2 Tbar0t0-24 t^3 Tbar0t0-120 h^2 t Th00+168 h t^2 Th00-96 t^3 Th00-72 h^3 Th0t+384 h^2 t Th0t-384 h t^2 Th0t+408 h t^2 Thht-408 t^3 Thht+48 h^2 t Tht0+24 h t^2 Tht0+24 h^2 t Tth0+48 h t^2 Tth0-144 t^3 Tth0+72 h t^2 U000h+576 t^3 U000t-12 t^3 U0000-48 h^2 t U00h0+144 h t^2 U00h0-120 h t^2 Uht00+240 t^3 Uht00-72 h^2 t Uhtht+384 t^3 Uhtht-96 h t^2 Uhtt0-96 h^2 t Uth00+288 h t^2 Uth00-288 h^2 t Uthhh+432 h t^2 Uthhh-864 h t^2 Uthtt+2688 t^3 Uthtt-48 h^2 t Uth00+240 h t^2 Uth00-24 h^2 t Ut0h0-96 h t^2 U0tht+96 h t^2 U0tt0+1/(h-t) 4 (3 h^4 (7-12 \[Pi]^2)+3 h^3 (22+19 \[Pi]^2) t+h^2 (33-23 \[Pi]^2) t^2+2 h (13 \[Pi]^2-72 (1+q)^2) t^3-24 (\[Pi]^2-6 (1+q)^2) t^4+3 h (h-t) t (h+2 t) Log[h]^2-6 t^2 (-2 h+3 t) (-h+4 t) Log[t]+36 (h-t) t^3 Log[t]^2+3 h Log[-h+t] (2 t (-5 h^2+h t+t^2)-(-h+t) (-4 h^2+h t+2 t^2) Log[-h+t])+6 h Log[h] (h t (h+5 t)+(-h+t) (2 t^2 Log[t]+h (-2 h+t) Log[-h+t]))+18 h^4 Log[-1+t/h]+6 h t (-h+t) (-(h+2 t) PolyLog[2,t/h]-(h-2 t) PolyLog[2,t/(-h+t)])))
