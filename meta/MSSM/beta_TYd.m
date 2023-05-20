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

{(-7*g1^2*hb)/15 - 3*g2^2*hb - (16*g3^2*hb)/3 + (14*g1^2*M1*Yb)/15 + 
  6*g2^2*M2*Yb + (32*g3^2*M3*Yb)/3 + 5*MatMul[hb, Adj[Yb], Yb] + 
  MatMul[hb, Adj[Yt], Yt] + 4*MatMul[Yb, Adj[Yb], hb] + 
  2*MatMul[Yb, Adj[Yt], ht] + 6*Yb*trace[hb, Adj[Yb]] + 
  2*Yb*trace[he, Adj[Ye]] + 3*hb*trace[Adj[Yb], Yb] + hb*trace[Adj[Ye], Ye], 
 (287*g1^4*hb)/90 + g1^2*g2^2*hb + (15*g2^4*hb)/2 + (8*g1^2*g3^2*hb)/9 + 
  8*g2^2*g3^2*hb - (16*g3^4*hb)/9 - (574*g1^4*M1*Yb)/45 - 2*g1^2*g2^2*M1*Yb - 
  (16*g1^2*g3^2*M1*Yb)/9 - 2*g1^2*g2^2*M2*Yb - 30*g2^4*M2*Yb - 
  16*g2^2*g3^2*M2*Yb - (16*g1^2*g3^2*M3*Yb)/9 - 16*g2^2*g3^2*M3*Yb + 
  (64*g3^4*M3*Yb)/9 + (6*g1^2*MatMul[hb, Adj[Yb], Yb])/5 + 
  12*g2^2*MatMul[hb, Adj[Yb], Yb] + (4*g1^2*MatMul[hb, Adj[Yt], Yt])/5 + 
  (6*g1^2*MatMul[Yb, Adj[Yb], hb])/5 + 6*g2^2*MatMul[Yb, Adj[Yb], hb] - 
  (8*g1^2*M1*MatMul[Yb, Adj[Yb], Yb])/5 - 
  12*g2^2*M2*MatMul[Yb, Adj[Yb], Yb] + (8*g1^2*MatMul[Yb, Adj[Yt], ht])/5 - 
  (8*g1^2*M1*MatMul[Yb, Adj[Yt], Yt])/5 - 
  6*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb] - 
  2*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt] - 
  8*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb] - 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb] - 
  4*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb] - 
  4*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht] - (4*g1^2*Yb*trace[hb, Adj[Yb]])/5 + 
  32*g3^2*Yb*trace[hb, Adj[Yb]] - 18*MatMul[Yb, Adj[Yb], Yb]*
   trace[hb, Adj[Yb]] + (12*g1^2*Yb*trace[he, Adj[Ye]])/5 - 
  6*MatMul[Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]] - 
  6*MatMul[Yb, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  (2*g1^2*hb*trace[Adj[Yb], Yb])/5 + 16*g3^2*hb*trace[Adj[Yb], Yb] - 
  32*g3^2*M3*Yb*trace[Adj[Yb], Yb] - 15*MatMul[hb, Adj[Yb], Yb]*
   trace[Adj[Yb], Yb] - 12*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb] + 
  (6*g1^2*hb*trace[Adj[Ye], Ye])/5 - 5*MatMul[hb, Adj[Yb], Yb]*
   trace[Adj[Ye], Ye] - 4*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  M1*((4*g1^2*Yb*trace[Adj[Yb], Yb])/5 - (12*g1^2*Yb*trace[Adj[Ye], Ye])/5) - 
  3*MatMul[hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  6*MatMul[Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  6*Yb*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
  36*Yb*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
  9*hb*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  12*Yb*trace[Adj[Ye], he, Adj[Ye], Ye] - 
  3*hb*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  6*Yb*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
  3*hb*trace[Adj[Yt], Yt, Adj[Yb], Yb], (194651*g1^6*hb)/6750 + 
  (109*g1^4*g2^2*hb)/50 + (17*g1^2*g2^4*hb)/2 + (345*g2^6*hb)/2 + 
  (3892*g1^4*g3^2*hb)/225 - (8*g1^2*g2^2*g3^2*hb)/5 + 140*g2^4*g3^2*hb + 
  (212*g1^2*g3^4*hb)/9 + 68*g2^2*g3^4*hb + (5440*g3^6*hb)/27 + 
  (16*g1^2*g2^2*g3^2*M1*Yb)/5 + (16*g1^2*g2^2*g3^2*M2*Yb)/5 + 
  (16*g1^2*g2^2*g3^2*M3*Yb)/5 - (1393*g1^6*MatMul[hb, Zeta[3]])/125 - 
  (189*g1^4*g2^2*MatMul[hb, Zeta[3]])/25 - (81*g1^2*g2^4*MatMul[hb, Zeta[3]])/
   5 + 315*g2^6*MatMul[hb, Zeta[3]] - (616*g1^4*g3^2*MatMul[hb, Zeta[3]])/
   25 - 216*g2^4*g3^2*MatMul[hb, Zeta[3]] - 
  (176*g1^2*g3^4*MatMul[hb, Zeta[3]])/5 - 144*g2^2*g3^4*MatMul[hb, Zeta[3]] + 
  640*g3^6*MatMul[hb, Zeta[3]] + M1*((-194651*g1^6*Yb)/1125 + 
    (8358*g1^6*MatMul[Yb, Zeta[3]])/125) + 
  M2*((-109*g1^4*g2^2*Yb)/25 + (378*g1^4*g2^2*MatMul[Yb, Zeta[3]])/25) + 
  M1*((-218*g1^4*g2^2*Yb)/25 + (756*g1^4*g2^2*MatMul[Yb, Zeta[3]])/25) + 
  M1*(-17*g1^2*g2^4*Yb + (162*g1^2*g2^4*MatMul[Yb, Zeta[3]])/5) + 
  M2*(-34*g1^2*g2^4*Yb + (324*g1^2*g2^4*MatMul[Yb, Zeta[3]])/5) + 
  M2*(-1035*g2^6*Yb - 1890*g2^6*MatMul[Yb, Zeta[3]]) + 
  M3*((-7784*g1^4*g3^2*Yb)/225 + (1232*g1^4*g3^2*MatMul[Yb, Zeta[3]])/25) + 
  M1*((-15568*g1^4*g3^2*Yb)/225 + (2464*g1^4*g3^2*MatMul[Yb, Zeta[3]])/25) + 
  M3*(-280*g2^4*g3^2*Yb + 432*g2^4*g3^2*MatMul[Yb, Zeta[3]]) + 
  M2*(-560*g2^4*g3^2*Yb + 864*g2^4*g3^2*MatMul[Yb, Zeta[3]]) + 
  M1*((-424*g1^2*g3^4*Yb)/9 + (352*g1^2*g3^4*MatMul[Yb, Zeta[3]])/5) + 
  M3*((-848*g1^2*g3^4*Yb)/9 + (704*g1^2*g3^4*MatMul[Yb, Zeta[3]])/5) + 
  M2*(-136*g2^2*g3^4*Yb + 288*g2^2*g3^4*MatMul[Yb, Zeta[3]]) + 
  M3*(-272*g2^2*g3^4*Yb + 576*g2^2*g3^4*MatMul[Yb, Zeta[3]]) + 
  M3*((-10880*g3^6*Yb)/9 - 3840*g3^6*MatMul[Yb, Zeta[3]]) + 
  (7*g1^4*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]])/150 + 
  (93*g1^2*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]])/5 - 
  (99*g2^4*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]])/2 + 
  (448*g1^2*g3^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]])/15 + 
  192*g2^2*g3^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]] - 
  (1360*g3^4*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]])/3 + 
  (143*g1^4*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]])/150 + 
  9*g1^2*g2^2*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]] - 
  (63*g2^4*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]])/3 + 
  (28*g1^4*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]])/75 + 
  12*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]] - 
  72*g2^4*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]] + 
  (416*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]])/15 + 
  96*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]] - 
  (1088*g3^4*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]])/3 + 
  (143*g1^4*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]])/75 + 
  18*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]] - 
  63*g2^4*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]] + 
  (256*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]])/15 - 
  (544*g3^4*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]])/3 + 
  (6*g1^2*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb], Zeta[3]] - 
  6*g1^2*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  18*g2^2*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
  (6*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] + 
  (18*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] + 
  30*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  6*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  36*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  36*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] - 
  (8639*g1^4*MatMul[hb, Adj[Yb], Yb])/300 - 
  (213*g1^2*g2^2*MatMul[hb, Adj[Yb], Yb])/10 - 
  (393*g2^4*MatMul[hb, Adj[Yb], Yb])/4 - 
  (44*g1^2*g3^2*MatMul[hb, Adj[Yb], Yb])/3 - 
  180*g2^2*g3^2*MatMul[hb, Adj[Yb], Yb] + (40*g3^4*MatMul[hb, Adj[Yb], Yb])/
   3 - (3767*g1^4*MatMul[hb, Adj[Yt], Yt])/300 - 
  (59*g1^2*g2^2*MatMul[hb, Adj[Yt], Yt])/10 - 
  (45*g2^4*MatMul[hb, Adj[Yt], Yt])/4 - 
  (68*g1^2*g3^2*MatMul[hb, Adj[Yt], Yt])/5 - 
  4*g2^2*g3^2*MatMul[hb, Adj[Yt], Yt] + (8*g3^4*MatMul[hb, Adj[Yt], Yt])/3 - 
  (1792*g1^4*MatMul[Yb, Adj[Yb], hb])/75 - 
  (84*g1^2*g2^2*MatMul[Yb, Adj[Yb], hb])/5 - 
  66*g2^4*MatMul[Yb, Adj[Yb], hb] - (224*g1^2*g3^2*MatMul[Yb, Adj[Yb], hb])/
   15 - 96*g2^2*g3^2*MatMul[Yb, Adj[Yb], hb] + 
  (32*g3^4*MatMul[Yb, Adj[Yb], hb])/3 + 
  M1*((-14*g1^4*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/25 + 
    (5269*g1^4*MatMul[Yb, Adj[Yb], Yb])/75) + 
  M1*((-102*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
    (127*g1^2*g2^2*MatMul[Yb, Adj[Yb], Yb])/5) + 
  M2*((-102*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
    (127*g1^2*g2^2*MatMul[Yb, Adj[Yb], Yb])/5) + 
  M2*(162*g2^4*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]] + 
    219*g2^4*MatMul[Yb, Adj[Yb], Yb]) + 
  M1*((-192*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
    (296*g1^2*g3^2*MatMul[Yb, Adj[Yb], Yb])/15) + 
  M3*((-192*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
    (296*g1^2*g3^2*MatMul[Yb, Adj[Yb], Yb])/15) + 
  M2*(-192*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]] + 
    184*g2^2*g3^2*MatMul[Yb, Adj[Yb], Yb]) + 
  M3*(-192*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]] + 
    184*g2^2*g3^2*MatMul[Yb, Adj[Yb], Yb]) + 
  M3*(1088*g3^4*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]] - 
    32*g3^4*MatMul[Yb, Adj[Yb], Yb]) - (3767*g1^4*MatMul[Yb, Adj[Yt], ht])/
   150 - (59*g1^2*g2^2*MatMul[Yb, Adj[Yt], ht])/5 - 
  (45*g2^4*MatMul[Yb, Adj[Yt], ht])/2 - 
  (136*g1^2*g3^2*MatMul[Yb, Adj[Yt], ht])/5 - 
  8*g2^2*g3^2*MatMul[Yb, Adj[Yt], ht] + (16*g3^4*MatMul[Yb, Adj[Yt], ht])/3 + 
  8*g2^2*g3^2*M2*MatMul[Yb, Adj[Yt], Yt] + 
  8*g2^2*g3^2*M3*MatMul[Yb, Adj[Yt], Yt] + 
  M1*((-286*g1^4*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/75 + 
    (3767*g1^4*MatMul[Yb, Adj[Yt], Yt])/75) + 
  M1*(-18*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]] + 
    (59*g1^2*g2^2*MatMul[Yb, Adj[Yt], Yt])/5) + 
  M2*(-18*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]] + 
    (59*g1^2*g2^2*MatMul[Yb, Adj[Yt], Yt])/5) + 
  M2*(126*g2^4*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]] + 
    45*g2^4*MatMul[Yb, Adj[Yt], Yt]) + 
  M1*((-256*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/15 + 
    (136*g1^2*g3^2*MatMul[Yb, Adj[Yt], Yt])/5) + 
  M3*((-256*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/15 + 
    (136*g1^2*g3^2*MatMul[Yb, Adj[Yt], Yt])/5) + 
  M3*((1088*g3^4*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/3 - 
    (32*g3^4*MatMul[Yb, Adj[Yt], Yt])/3) - 
  (g1^2*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb])/5 + 
  15*g2^2*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  64*g3^2*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb] - 
  (58*g1^2*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb])/15 + 
  18*g2^2*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb])/3 + 
  (11*g1^2*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt])/3 - 
  3*g2^2*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt] + 
  (64*g3^2*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (4*g1^2*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb])/15 + 
  12*g2^2*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
  (256*g3^2*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb])/3 + 
  (3*g1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb])/5 + 
  3*g2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb] + 
  64*g3^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb] - 
  (4*g1^2*M1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb])/15 - 
  12*g2^2*M2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb] - 
  (256*g3^2*M3*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb])/3 - 
  (58*g1^2*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb])/15 + 
  18*g2^2*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb])/3 + 
  (22*g1^2*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt])/3 - 
  6*g2^2*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt])/3 - 
  (29*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb])/15 + 
  9*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb] + 
  (64*g3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb])/3 - 
  (128*g3^2*M3*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb], Zeta[3]])/5 + 
    (58*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb])/15) + 
  M2*(36*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb], Zeta[3]] - 
    18*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]) + 
  (22*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht])/3 - 
  6*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht] + 
  (128*g3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht])/3 - 
  (128*g3^2*M3*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  M1*(12*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
    (22*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt])/3) + 
  M2*(-36*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
    6*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt]) + 
  12*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb], Yb] + 
  4*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  12*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  4*MatMul[Yb, Adj[Yb], hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], hb] + 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], hb] - 
  4*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb, Adj[Yb], Yb] + 
  8*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb, Adj[Yt], Yt] + 
  12*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt, Adj[Yb], Yb] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb, Adj[Yb], Yb] + 
  8*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb, Adj[Yt], Yt] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb], hb] + 
  8*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], ht] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht, Adj[Yb], Yb] + 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], hb] - 
  (63*g1^4*Yb*trace[hb, Adj[Yb]])/2 - (3*g1^2*g2^2*Yb*trace[hb, Adj[Yb]])/5 - 
  (315*g2^4*Yb*trace[hb, Adj[Yb]])/2 - (568*g1^2*g3^2*Yb*trace[hb, Adj[Yb]])/
   15 - 264*g2^2*g3^2*Yb*trace[hb, Adj[Yb]] - 
  (640*g3^4*Yb*trace[hb, Adj[Yb]])/3 - 
  (77*g1^4*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yb]])/25 - 
  18*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yb]] - 
  189*g2^4*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yb]] + 
  (224*g1^2*g3^2*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yb]])/5 + 
  288*g2^2*g3^2*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yb]] - 
  32*g3^4*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yb]] - 
  (108*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]])/5 - 
  108*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]] + 
  288*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]] + 
  (102*g1^2*MatMul[Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]])/5 + 
  90*g2^2*MatMul[Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]] - 
  48*g3^2*MatMul[Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]] + 
  24*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]] - 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[hb, Adj[Yb]] - 
  (633*g1^4*Yb*trace[he, Adj[Ye]])/10 - (81*g1^2*g2^2*Yb*trace[he, Adj[Ye]])/
   5 - (105*g2^4*Yb*trace[he, Adj[Ye]])/2 + 
  (81*g1^4*MatMul[Yb, Zeta[3]]*trace[he, Adj[Ye]])/25 + 
  (162*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[he, Adj[Ye]])/5 - 
  63*g2^4*MatMul[Yb, Zeta[3]]*trace[he, Adj[Ye]] + 
  (84*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[he, Adj[Ye]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[he, Adj[Ye]] - 
  (46*g1^2*MatMul[Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]])/5 + 
  30*g2^2*MatMul[Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  48*g3^2*MatMul[Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  8*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]] - 
  (182*g1^4*Yb*trace[ht, Adj[Yt]])/15 - 90*g2^4*Yb*trace[ht, Adj[Yt]] - 
  (320*g3^4*Yb*trace[ht, Adj[Yt]])/3 - 
  (48*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]])/5 + 
  96*g3^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]] + 
  4*g1^2*MatMul[Yb, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  36*g2^2*MatMul[Yb, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  16*g3^2*MatMul[Yb, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  24*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[ht, Adj[Yt]] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  (63*g1^4*hb*trace[Adj[Yb], Yb])/4 - (3*g1^2*g2^2*hb*trace[Adj[Yb], Yb])/
   10 - (315*g2^4*hb*trace[Adj[Yb], Yb])/4 - 
  (284*g1^2*g3^2*hb*trace[Adj[Yb], Yb])/15 - 
  132*g2^2*g3^2*hb*trace[Adj[Yb], Yb] - (320*g3^4*hb*trace[Adj[Yb], Yb])/3 - 
  (77*g1^4*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb])/50 - 
  9*g1^2*g2^2*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (189*g2^4*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb])/2 + 
  (112*g1^2*g3^2*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  144*g2^2*g3^2*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb] - 
  16*g3^4*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (84*g1^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  108*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  240*g3^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (78*g1^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  54*g2^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  192*g3^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (86*g1^2*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 + 
  72*g2^2*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  40*g3^2*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  (67*g1^2*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb])/5 + 
  63*g2^2*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
  32*g3^2*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb] + 
  18*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  12*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  24*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  18*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
  108*MatMul[Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  36*MatMul[Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  45*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb]^2 - 
  36*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb]^2 + 
  M1*((568*g1^2*g3^2*Yb*trace[Adj[Yb], Yb])/15 - 
    (224*g1^2*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb])/5) + 
  M3*((568*g1^2*g3^2*Yb*trace[Adj[Yb], Yb])/15 - 
    (224*g1^2*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb])/5) + 
  M2*(264*g2^2*g3^2*Yb*trace[Adj[Yb], Yb] - 288*g2^2*g3^2*MatMul[Yb, Zeta[3]]*
     trace[Adj[Yb], Yb]) + M3*(264*g2^2*g3^2*Yb*trace[Adj[Yb], Yb] - 
    288*g2^2*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb]) - 
  (633*g1^4*hb*trace[Adj[Ye], Ye])/20 - (81*g1^2*g2^2*hb*trace[Adj[Ye], Ye])/
   10 - (105*g2^4*hb*trace[Adj[Ye], Ye])/4 + 
  (81*g1^4*MatMul[hb, Zeta[3]]*trace[Adj[Ye], Ye])/50 + 
  (81*g1^2*g2^2*MatMul[hb, Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  (63*g2^4*MatMul[hb, Zeta[3]]*trace[Adj[Ye], Ye])/2 + 
  (72*g1^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  36*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye] + 
  (54*g1^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], hb], Zeta[3]]*trace[Adj[Ye], Ye] - 
  (38*g1^2*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
  24*g2^2*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  40*g3^2*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  (31*g1^2*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye])/5 + 
  21*g2^2*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  32*g3^2*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  6*MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  4*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  8*MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb]*trace[Adj[Ye], Ye] - 
  36*MatMul[Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  30*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  24*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  5*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye]^2 - 
  4*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye]^2 + 
  M1*((3*g1^2*g2^2*Yb*trace[Adj[Yb], Yb])/5 + 
    18*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb] + 
    (81*g1^2*g2^2*Yb*trace[Adj[Ye], Ye])/5 - 
    (162*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M2*((3*g1^2*g2^2*Yb*trace[Adj[Yb], Yb])/5 + 
    18*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb] + 
    (81*g1^2*g2^2*Yb*trace[Adj[Ye], Ye])/5 - 
    (162*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M1*((108*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/
     5 - (102*g1^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 - 
    (84*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
    (46*g1^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5) + 
  M2*(108*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
    90*g2^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
    36*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye] - 
    30*g2^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye]) + 
  M3*(-288*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
    48*g3^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
    48*g3^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye]) - 
  (91*g1^4*hb*trace[Adj[Yt], Yt])/15 - 45*g2^4*hb*trace[Adj[Yt], Yt] - 
  (160*g3^4*hb*trace[Adj[Yt], Yt])/3 - 
  (24*g1^2*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  48*g3^2*MatMul[MatMul[hb, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
  (48*g1^2*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  96*g3^2*MatMul[MatMul[Yb, Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt] + 
  2*g1^2*MatMul[hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  8*g3^2*MatMul[hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  4*g1^2*MatMul[Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  36*g2^2*MatMul[Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  16*g3^2*MatMul[Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  36*g2^2*M2*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  24*MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[Adj[Yt], Yt] + 
  6*MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  24*MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  36*MatMul[Yb, Adj[Yt], Yt]*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  9*MatMul[hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt]^2 - 
  18*MatMul[Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt]^2 + 
  M1*(63*g1^4*Yb*trace[Adj[Yb], Yb] + (154*g1^4*MatMul[Yb, Zeta[3]]*
      trace[Adj[Yb], Yb])/25 + (633*g1^4*Yb*trace[Adj[Ye], Ye])/5 - 
    (162*g1^4*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye])/25 + 
    (364*g1^4*Yb*trace[Adj[Yt], Yt])/15) + 
  M2*(315*g2^4*Yb*trace[Adj[Yb], Yb] + 378*g2^4*MatMul[Yb, Zeta[3]]*
     trace[Adj[Yb], Yb] + 105*g2^4*Yb*trace[Adj[Ye], Ye] + 
    126*g2^4*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye] + 
    180*g2^4*Yb*trace[Adj[Yt], Yt]) + 
  M3*((1280*g3^4*Yb*trace[Adj[Yb], Yb])/3 + 64*g3^4*MatMul[Yb, Zeta[3]]*
     trace[Adj[Yb], Yb] + (640*g3^4*Yb*trace[Adj[Yt], Yt])/3) + 
  M1*((48*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/
     5 - 4*g1^2*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt]) + 
  M3*(-96*g3^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] + 
    16*g3^2*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt]) - 
  (24*g1^2*Yb*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 + 
  36*g2^2*Yb*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  48*g3^2*Yb*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  (84*g1^2*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 - 
  96*g3^2*MatMul[Yb, Zeta[3]]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  36*MatMul[Yb, Adj[Yb], Yb]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], Yt]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  36*Yb*trace[Adj[Yt], Yt]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*g1^2*Yb*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  36*g2^2*Yb*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  288*g3^2*Yb*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  (216*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], hb, Adj[Yb], Yb])/5 + 
  216*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
  576*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  216*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  216*Yb*trace[Adj[Yb], Yb]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  72*Yb*trace[Adj[Ye], Ye]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  3*g1^2*hb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g2^2*hb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*g3^2*hb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  (54*g1^2*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb])/5 + 
  54*g2^2*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  144*g3^2*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  90*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  108*Yb*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*Yb*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  54*hb*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*hb*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*g1^2*Yb*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  12*g2^2*Yb*trace[Adj[Ye], he, Adj[Ye], Ye] - 
  (216*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], he, Adj[Ye], Ye])/5 + 
  72*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  72*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  72*Yb*trace[Adj[Yb], Yb]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  24*Yb*trace[Adj[Ye], Ye]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  9*g1^2*hb*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  3*g2^2*hb*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (54*g1^2*MatMul[hb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  18*g2^2*MatMul[hb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  30*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  24*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  36*Yb*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*Yb*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*hb*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*hb*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (24*g1^2*Yb*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 + 
  36*g2^2*Yb*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  48*g3^2*Yb*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  (84*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 - 
  96*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  36*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  36*Yb*trace[Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  72*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
  (12*g1^2*hb*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  18*g2^2*hb*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*g3^2*hb*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (42*g1^2*MatMul[hb, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
  48*g3^2*MatMul[hb, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  30*MatMul[hb, Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  6*MatMul[hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*MatMul[Yb, Adj[Yb], hb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*Yb*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*hb*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  M2*(-18*g2^2*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    108*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    6*g2^2*Yb*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 36*g2^2*MatMul[Yb, Zeta[3]]*
     trace[Adj[Ye], Ye, Adj[Ye], Ye] - 36*g2^2*Yb*trace[Adj[Yt], Yt, Adj[Yb], 
      Yb]) + M1*(-6*g1^2*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    (108*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb])/5 - 
    18*g1^2*Yb*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    (108*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
    (24*g1^2*Yb*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
    (84*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5) + 
  M3*(-144*g3^2*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    288*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    48*g3^2*Yb*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 96*g3^2*MatMul[Yb, Zeta[3]]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 18*MatMul[hb, Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 36*MatMul[Yb, Adj[Yt], ht]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  18*Yb*trace[hb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  18*Yb*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
  108*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
  3*hb*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*MatMul[hb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  6*Yb*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
  36*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
  hb*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*MatMul[hb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*Yb*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*Yb*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yb], Yb] + 
  9*hb*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb]}
