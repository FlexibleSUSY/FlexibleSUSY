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

{(-13*g1^2*ht)/15 - 3*g2^2*ht - (16*g3^2*ht)/3 + (26*g1^2*M1*Yt)/15 + 
  6*g2^2*M2*Yt + (32*g3^2*M3*Yt)/3 + MatMul[ht, Adj[Yb], Yb] + 
  5*MatMul[ht, Adj[Yt], Yt] + 2*MatMul[Yt, Adj[Yb], hb] + 
  4*MatMul[Yt, Adj[Yt], ht] + 6*Yt*trace[ht, Adj[Yt]] + 
  3*ht*trace[Adj[Yt], Yt], (2743*g1^4*ht)/450 + g1^2*g2^2*ht + 
  (15*g2^4*ht)/2 + (136*g1^2*g3^2*ht)/45 + 8*g2^2*g3^2*ht - (16*g3^4*ht)/9 - 
  (5486*g1^4*M1*Yt)/225 - 2*g1^2*g2^2*M1*Yt - (272*g1^2*g3^2*M1*Yt)/45 - 
  2*g1^2*g2^2*M2*Yt - 30*g2^4*M2*Yt - 16*g2^2*g3^2*M2*Yt - 
  (272*g1^2*g3^2*M3*Yt)/45 - 16*g2^2*g3^2*M3*Yt + (64*g3^4*M3*Yt)/9 + 
  (2*g1^2*MatMul[ht, Adj[Yb], Yb])/5 + 12*g2^2*MatMul[ht, Adj[Yt], Yt] + 
  (4*g1^2*MatMul[Yt, Adj[Yb], hb])/5 - (4*g1^2*M1*MatMul[Yt, Adj[Yb], Yb])/
   5 + (6*g1^2*MatMul[Yt, Adj[Yt], ht])/5 + 6*g2^2*MatMul[Yt, Adj[Yt], ht] - 
  (4*g1^2*M1*MatMul[Yt, Adj[Yt], Yt])/5 - 
  12*g2^2*M2*MatMul[Yt, Adj[Yt], Yt] - 
  2*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt] - 
  6*MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt] - 
  4*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb] - 
  4*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt] - 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb] - 
  2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht] - 
  8*MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt] - 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht] - 6*MatMul[Yt, Adj[Yb], Yb]*
   trace[hb, Adj[Yb]] - 2*MatMul[Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  (8*g1^2*Yt*trace[ht, Adj[Yt]])/5 + 32*g3^2*Yt*trace[ht, Adj[Yt]] - 
  18*MatMul[Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  3*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  6*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
  MatMul[ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 2*MatMul[Yt, Adj[Yb], hb]*
   trace[Adj[Ye], Ye] + (4*g1^2*ht*trace[Adj[Yt], Yt])/5 + 
  16*g3^2*ht*trace[Adj[Yt], Yt] - (8*g1^2*M1*Yt*trace[Adj[Yt], Yt])/5 - 
  32*g3^2*M3*Yt*trace[Adj[Yt], Yt] - 15*MatMul[ht, Adj[Yt], Yt]*
   trace[Adj[Yt], Yt] - 12*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  6*Yt*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
  6*Yt*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
  36*Yt*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
  3*ht*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  9*ht*trace[Adj[Yt], Yt, Adj[Yt], Yt], (352097*g1^6*ht)/6750 + 
  (379*g1^4*g2^2*ht)/50 + (17*g1^2*g2^4*ht)/2 + (345*g2^6*ht)/2 + 
  (5308*g1^4*g3^2*ht)/225 - (8*g1^2*g2^2*g3^2*ht)/5 + 140*g2^4*g3^2*ht + 
  (436*g1^2*g3^4*ht)/45 + 68*g2^2*g3^4*ht + (5440*g3^6*ht)/27 + 
  (16*g1^2*g2^2*g3^2*M1*Yt)/5 + (16*g1^2*g2^2*g3^2*M2*Yt)/5 + 
  (16*g1^2*g2^2*g3^2*M3*Yt)/5 - (2587*g1^6*MatMul[ht, Zeta[3]])/125 - 
  (351*g1^4*g2^2*MatMul[ht, Zeta[3]])/25 - (81*g1^2*g2^4*MatMul[ht, Zeta[3]])/
   5 + 315*g2^6*MatMul[ht, Zeta[3]] - (1144*g1^4*g3^2*MatMul[ht, Zeta[3]])/
   25 - 216*g2^4*g3^2*MatMul[ht, Zeta[3]] - 
  (176*g1^2*g3^4*MatMul[ht, Zeta[3]])/5 - 144*g2^2*g3^4*MatMul[ht, Zeta[3]] + 
  640*g3^6*MatMul[ht, Zeta[3]] + M1*((-352097*g1^6*Yt)/1125 + 
    (15522*g1^6*MatMul[Yt, Zeta[3]])/125) + 
  M2*((-379*g1^4*g2^2*Yt)/25 + (702*g1^4*g2^2*MatMul[Yt, Zeta[3]])/25) + 
  M1*((-758*g1^4*g2^2*Yt)/25 + (1404*g1^4*g2^2*MatMul[Yt, Zeta[3]])/25) + 
  M1*(-17*g1^2*g2^4*Yt + (162*g1^2*g2^4*MatMul[Yt, Zeta[3]])/5) + 
  M2*(-34*g1^2*g2^4*Yt + (324*g1^2*g2^4*MatMul[Yt, Zeta[3]])/5) + 
  M2*(-1035*g2^6*Yt - 1890*g2^6*MatMul[Yt, Zeta[3]]) + 
  M3*((-10616*g1^4*g3^2*Yt)/225 + (2288*g1^4*g3^2*MatMul[Yt, Zeta[3]])/25) + 
  M1*((-21232*g1^4*g3^2*Yt)/225 + (4576*g1^4*g3^2*MatMul[Yt, Zeta[3]])/25) + 
  M3*(-280*g2^4*g3^2*Yt + 432*g2^4*g3^2*MatMul[Yt, Zeta[3]]) + 
  M2*(-560*g2^4*g3^2*Yt + 864*g2^4*g3^2*MatMul[Yt, Zeta[3]]) + 
  M1*((-872*g1^2*g3^4*Yt)/45 + (352*g1^2*g3^4*MatMul[Yt, Zeta[3]])/5) + 
  M3*((-1744*g1^2*g3^4*Yt)/45 + (704*g1^2*g3^4*MatMul[Yt, Zeta[3]])/5) + 
  M2*(-136*g2^2*g3^4*Yt + 288*g2^2*g3^4*MatMul[Yt, Zeta[3]]) + 
  M3*(-272*g2^2*g3^4*Yt + 576*g2^2*g3^4*MatMul[Yt, Zeta[3]]) + 
  M3*((-10880*g3^6*Yt)/9 - 3840*g3^6*MatMul[Yt, Zeta[3]]) + 
  (7*g1^4*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]])/30 + 
  (9*g1^2*g2^2*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]])/5 - 
  (63*g2^4*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]])/3 - 
  (169*g1^4*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]])/30 + 
  (201*g1^2*g2^2*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]])/5 - 
  (99*g2^4*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]])/2 - 
  (64*g1^2*g3^2*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]])/3 + 
  192*g2^2*g3^2*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]] - 
  (1360*g3^4*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]])/3 + 
  (7*g1^4*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]])/15 + 
  (18*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]])/5 - 
  63*g2^4*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]] + 
  (256*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]])/15 - 
  (544*g3^4*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]])/3 - 
  (104*g1^4*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]])/75 + 
  (168*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]])/5 - 
  72*g2^4*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]] + 
  (32*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]])/15 + 
  96*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]] - 
  (1088*g3^4*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]])/3 - 
  (6*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt], Zeta[3]] + 
  6*g1^2*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
  18*g2^2*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]] + 
  (18*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht], Zeta[3]] - 
  6*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] + 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] + 
  6*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  30*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]] + 
  36*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  36*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] - 
  (633*g1^4*MatMul[ht, Adj[Yb], Yb])/100 - 
  (41*g1^2*g2^2*MatMul[ht, Adj[Yb], Yb])/10 - 
  (45*g2^4*MatMul[ht, Adj[Yb], Yb])/4 - 
  (76*g1^2*g3^2*MatMul[ht, Adj[Yb], Yb])/15 - 
  4*g2^2*g3^2*MatMul[ht, Adj[Yb], Yb] + (8*g3^4*MatMul[ht, Adj[Yb], Yb])/3 - 
  (2671*g1^4*MatMul[ht, Adj[Yt], Yt])/60 - 
  (327*g1^2*g2^2*MatMul[ht, Adj[Yt], Yt])/10 - 
  (393*g2^4*MatMul[ht, Adj[Yt], Yt])/4 - 
  (44*g1^2*g3^2*MatMul[ht, Adj[Yt], Yt])/3 - 
  180*g2^2*g3^2*MatMul[ht, Adj[Yt], Yt] + (40*g3^4*MatMul[ht, Adj[Yt], Yt])/
   3 - (633*g1^4*MatMul[Yt, Adj[Yb], hb])/50 - 
  (41*g1^2*g2^2*MatMul[Yt, Adj[Yb], hb])/5 - 
  (45*g2^4*MatMul[Yt, Adj[Yb], hb])/2 - 
  (152*g1^2*g3^2*MatMul[Yt, Adj[Yb], hb])/15 - 
  8*g2^2*g3^2*MatMul[Yt, Adj[Yb], hb] + (16*g3^4*MatMul[Yt, Adj[Yb], hb])/3 + 
  8*g2^2*g3^2*M2*MatMul[Yt, Adj[Yb], Yb] + 
  8*g2^2*g3^2*M3*MatMul[Yt, Adj[Yb], Yb] + 
  M1*((-14*g1^4*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]])/15 + 
    (633*g1^4*MatMul[Yt, Adj[Yb], Yb])/25) + 
  M1*((-18*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]])/5 + 
    (41*g1^2*g2^2*MatMul[Yt, Adj[Yb], Yb])/5) + 
  M2*((-18*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]])/5 + 
    (41*g1^2*g2^2*MatMul[Yt, Adj[Yb], Yb])/5) + 
  M2*(126*g2^4*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]] + 
    45*g2^4*MatMul[Yt, Adj[Yb], Yb]) + 
  M1*((-256*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]])/15 + 
    (152*g1^2*g3^2*MatMul[Yt, Adj[Yb], Yb])/15) + 
  M3*((-256*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]])/15 + 
    (152*g1^2*g3^2*MatMul[Yt, Adj[Yb], Yb])/15) + 
  M3*((1088*g3^4*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]])/3 - 
    (32*g3^4*MatMul[Yt, Adj[Yb], Yb])/3) - 
  (3082*g1^4*MatMul[Yt, Adj[Yt], ht])/75 - 
  (126*g1^2*g2^2*MatMul[Yt, Adj[Yt], ht])/5 - 
  66*g2^4*MatMul[Yt, Adj[Yt], ht] - (416*g1^2*g3^2*MatMul[Yt, Adj[Yt], ht])/
   15 - 96*g2^2*g3^2*MatMul[Yt, Adj[Yt], ht] + 
  (32*g3^4*MatMul[Yt, Adj[Yt], ht])/3 + 
  M1*((234*g1^4*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]])/25 + 
    (8561*g1^4*MatMul[Yt, Adj[Yt], Yt])/75) + 
  M1*((-246*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]])/5 + 
    (193*g1^2*g2^2*MatMul[Yt, Adj[Yt], Yt])/5) + 
  M2*((-246*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]])/5 + 
    (193*g1^2*g2^2*MatMul[Yt, Adj[Yt], Yt])/5) + 
  M2*(162*g2^4*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]] + 
    219*g2^4*MatMul[Yt, Adj[Yt], Yt]) + 
  M1*((64*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]])/5 + 
    (424*g1^2*g3^2*MatMul[Yt, Adj[Yt], Yt])/15) + 
  M3*((64*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]])/5 + 
    (424*g1^2*g3^2*MatMul[Yt, Adj[Yt], Yt])/15) + 
  M2*(-192*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]] + 
    184*g2^2*g3^2*MatMul[Yt, Adj[Yt], Yt]) + 
  M3*(-192*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]] + 
    184*g2^2*g3^2*MatMul[Yt, Adj[Yt], Yt]) + 
  M3*(1088*g3^4*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]] - 
    32*g3^4*MatMul[Yt, Adj[Yt], Yt]) + 
  (7*g1^2*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb])/15 - 
  3*g2^2*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb] + 
  (64*g3^2*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (38*g1^2*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt])/15 + 
  18*g2^2*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt])/3 + 
  3*g1^2*MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt] + 
  15*g2^2*MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt] + 
  64*g3^2*MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt] + 
  (14*g1^2*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb])/15 - 
  6*g2^2*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb])/3 + 
  (38*g1^2*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt])/15 + 
  18*g2^2*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt])/3 + 
  (14*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb])/15 - 
  6*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb] + 
  (128*g3^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb])/3 - 
  (128*g3^2*M3*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  M1*((12*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 - 
    (14*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb])/15) + 
  M2*(-36*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
    6*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb]) + 
  (19*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht])/15 + 
  9*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht] + 
  (64*g3^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht])/3 - 
  (128*g3^2*M3*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt], Zeta[3]])/5 - 
    (38*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt])/15) + 
  M2*(36*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt], Zeta[3]] - 
    18*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt]) + 
  (20*g1^2*MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt])/3 + 
  12*g2^2*MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt] + 
  (256*g3^2*MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt])/3 + 
  7*g1^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht] + 
  3*g2^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht] + 
  64*g3^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht] - 
  (20*g1^2*M1*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt])/3 - 
  12*g2^2*M2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt] - 
  (256*g3^2*M3*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  12*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt], Yt] + 
  4*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] - 
  4*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt], Yt] + 
  12*MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt] + 
  12*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt, Adj[Yb], Yb] - 
  4*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt, Adj[Yt], Yt] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb, Adj[Yt], Yt] + 
  6*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt], ht] + 
  8*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht, Adj[Yb], Yb] - 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht, Adj[Yt], Yt] + 
  8*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], hb] - 
  2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt], ht] + 
  4*MatMul[Yt, Adj[Yt], ht, Adj[Yb], Yb, Adj[Yt], Yt] + 
  12*MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt, Adj[Yt], Yt] + 
  4*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], hb, Adj[Yt], Yt] + 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], ht] + 
  12*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt] + 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], ht] - 
  (182*g1^4*Yt*trace[hb, Adj[Yb]])/15 - 90*g2^4*Yt*trace[hb, Adj[Yb]] - 
  (320*g3^4*Yt*trace[hb, Adj[Yb]])/3 - 
  (48*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]])/5 + 
  96*g3^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]] + 
  (32*g1^2*MatMul[Yt, Adj[Yb], Yb]*trace[hb, Adj[Yb]])/5 + 
  36*g2^2*MatMul[Yt, Adj[Yb], Yb]*trace[hb, Adj[Yb]] - 
  16*g3^2*MatMul[Yt, Adj[Yb], Yb]*trace[hb, Adj[Yb]] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]] + 
  24*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt]*trace[hb, Adj[Yb]] - 
  (78*g1^4*Yt*trace[he, Adj[Ye]])/5 - 30*g2^4*Yt*trace[he, Adj[Ye]] + 
  (24*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]]*trace[he, Adj[Ye]])/5 - 
  (16*g1^2*MatMul[Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]])/5 + 
  12*g2^2*MatMul[Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  16*g3^2*MatMul[Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  8*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt]*trace[he, Adj[Ye]] - 
  (171*g1^4*Yt*trace[ht, Adj[Yt]])/2 - (57*g1^2*g2^2*Yt*trace[ht, Adj[Yt]])/
   5 - (315*g2^4*Yt*trace[ht, Adj[Yt]])/2 - 
  (248*g1^2*g3^2*Yt*trace[ht, Adj[Yt]])/3 - 
  264*g2^2*g3^2*Yt*trace[ht, Adj[Yt]] - (640*g3^4*Yt*trace[ht, Adj[Yt]])/3 - 
  (13*g1^4*MatMul[Yt, Zeta[3]]*trace[ht, Adj[Yt]])/5 + 
  (126*g1^2*g2^2*MatMul[Yt, Zeta[3]]*trace[ht, Adj[Yt]])/5 - 
  189*g2^4*MatMul[Yt, Zeta[3]]*trace[ht, Adj[Yt]] + 
  (416*g1^2*g3^2*MatMul[Yt, Zeta[3]]*trace[ht, Adj[Yt]])/5 + 
  288*g2^2*g3^2*MatMul[Yt, Zeta[3]]*trace[ht, Adj[Yt]] - 
  32*g3^4*MatMul[Yt, Zeta[3]]*trace[ht, Adj[Yt]] + 
  (36*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]])/5 - 
  108*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]] + 
  288*g3^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]] + 
  18*g1^2*MatMul[Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  90*g2^2*MatMul[Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  48*g3^2*MatMul[Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  24*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
  (91*g1^4*ht*trace[Adj[Yb], Yb])/15 - 45*g2^4*ht*trace[Adj[Yb], Yb] - 
  (160*g3^4*ht*trace[Adj[Yb], Yb])/3 - 
  (24*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  48*g3^2*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (48*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  96*g3^2*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (16*g1^2*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 + 
  18*g2^2*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  8*g3^2*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  (32*g1^2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb])/5 + 
  36*g2^2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
  16*g3^2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb] + 
  6*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  24*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht]*trace[Adj[Yb], Yb] - 
  36*MatMul[Yt, Adj[Yb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  9*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb]^2 - 
  (39*g1^4*ht*trace[Adj[Ye], Ye])/5 - 15*g2^4*ht*trace[Adj[Ye], Ye] + 
  (12*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
  (24*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  (8*g1^2*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
  6*g2^2*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  8*g3^2*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  (16*g1^2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Ye], Ye])/5 + 
  12*g2^2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  16*g3^2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  2*MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  8*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  8*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht]*trace[Adj[Ye], Ye] - 
  12*MatMul[Yt, Adj[Yb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yt, Adj[Yb], Yb]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  6*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  MatMul[ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye]^2 - 
  2*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Ye], Ye]^2 + 
  M1*((48*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/
     5 - (32*g1^2*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 - 
    (24*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
    (16*g1^2*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5) + 
  M2*(-36*g2^2*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
    12*g2^2*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Ye], Ye]) + 
  M3*(-96*g3^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
    16*g3^2*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
    16*g3^2*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Ye], Ye]) - 
  (171*g1^4*ht*trace[Adj[Yt], Yt])/4 - (57*g1^2*g2^2*ht*trace[Adj[Yt], Yt])/
   10 - (315*g2^4*ht*trace[Adj[Yt], Yt])/4 - 
  (124*g1^2*g3^2*ht*trace[Adj[Yt], Yt])/3 - 
  132*g2^2*g3^2*ht*trace[Adj[Yt], Yt] - (320*g3^4*ht*trace[Adj[Yt], Yt])/3 - 
  (13*g1^4*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt])/10 + 
  (63*g1^2*g2^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  (189*g2^4*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt])/2 + 
  (208*g1^2*g3^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  144*g2^2*g3^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt] - 
  16*g3^4*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt] + 
  12*g1^2*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
  108*g2^2*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] + 
  240*g3^2*MatMul[MatMul[ht, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
  (6*g1^2*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  54*g2^2*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt] + 
  192*g3^2*MatMul[MatMul[Yt, Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt] + 
  16*g1^2*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  72*g2^2*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  40*g3^2*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  11*g1^2*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  63*g2^2*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  32*g3^2*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  12*MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  6*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  24*MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
  108*MatMul[Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  45*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt]^2 - 
  36*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt]^2 + 
  M1*((364*g1^4*Yt*trace[Adj[Yb], Yb])/15 + (156*g1^4*Yt*trace[Adj[Ye], Ye])/
     5 + 171*g1^4*Yt*trace[Adj[Yt], Yt] + 
    (26*g1^4*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M1*((57*g1^2*g2^2*Yt*trace[Adj[Yt], Yt])/5 - 
    (126*g1^2*g2^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M2*((57*g1^2*g2^2*Yt*trace[Adj[Yt], Yt])/5 - 
    (126*g1^2*g2^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M2*(180*g2^4*Yt*trace[Adj[Yb], Yb] + 60*g2^4*Yt*trace[Adj[Ye], Ye] + 
    315*g2^4*Yt*trace[Adj[Yt], Yt] + 378*g2^4*MatMul[Yt, Zeta[3]]*
     trace[Adj[Yt], Yt]) + M1*((248*g1^2*g3^2*Yt*trace[Adj[Yt], Yt])/3 - 
    (416*g1^2*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M3*((248*g1^2*g3^2*Yt*trace[Adj[Yt], Yt])/3 - 
    (416*g1^2*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M2*(264*g2^2*g3^2*Yt*trace[Adj[Yt], Yt] - 288*g2^2*g3^2*MatMul[Yt, Zeta[3]]*
     trace[Adj[Yt], Yt]) + M3*(264*g2^2*g3^2*Yt*trace[Adj[Yt], Yt] - 
    288*g2^2*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M3*((640*g3^4*Yt*trace[Adj[Yb], Yb])/3 + (1280*g3^4*Yt*trace[Adj[Yt], Yt])/
     3 + 64*g3^4*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M1*((-36*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/
     5 - 18*g1^2*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt]) + 
  M2*(108*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
    90*g2^2*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt]) + 
  M3*(-288*g3^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] + 
    48*g3^2*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt]) + 
  (12*g1^2*Yt*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 + 
  36*g2^2*Yt*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  48*g3^2*Yt*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  (12*g1^2*MatMul[Yt, Zeta[3]]*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 - 
  96*g3^2*MatMul[Yt, Zeta[3]]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yt, Adj[Yb], Yb]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  36*MatMul[Yt, Adj[Yt], Yt]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  36*Yt*trace[Adj[Yb], Yb]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*Yt*trace[Adj[Ye], Ye]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  72*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  18*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  6*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  (12*g1^2*Yt*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 + 
  36*g2^2*Yt*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  48*g3^2*Yt*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  (12*g1^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 - 
  96*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yb], Yb]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  36*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  36*Yt*trace[Adj[Yb], Yb]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  12*Yt*trace[Adj[Ye], Ye]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  (228*g1^2*Yt*trace[Adj[Yt], ht, Adj[Yt], Yt])/5 + 
  36*g2^2*Yt*trace[Adj[Yt], ht, Adj[Yt], Yt] + 
  288*g3^2*Yt*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
  (72*g1^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yt], Yt])/5 + 
  216*g2^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
  576*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yt], Yt] + 
  216*MatMul[Yt, Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yt], Yt] + 
  216*Yt*trace[Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yt], Yt] + 
  (6*g1^2*ht*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  18*g2^2*ht*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*g3^2*ht*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (6*g1^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
  48*g3^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  6*MatMul[ht, Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  30*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yb], hb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*Yt*trace[hb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*Yt*trace[he, Adj[Ye]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*ht*trace[Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  6*ht*trace[Adj[Ye], Ye]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (57*g1^2*ht*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
  9*g2^2*ht*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  72*g3^2*ht*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 
  (18*g1^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
  54*g2^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 
  144*g3^2*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  90*MatMul[ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  72*MatMul[Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  108*Yt*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  54*ht*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  M1*((-12*g1^2*Yt*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
    (12*g1^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
    (114*g1^2*Yt*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
    (36*g1^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5) + 
  M2*(-36*g2^2*Yt*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
    18*g2^2*Yt*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 108*g2^2*MatMul[Yt, Zeta[3]]*
     trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  M3*(-48*g3^2*Yt*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    96*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
    144*g3^2*Yt*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
    288*g3^2*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  18*Yt*trace[Adj[Yb], hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*Yt*trace[Adj[Yb], Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  9*ht*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*Yt*trace[Adj[Yt], Yt, Adj[Yb], hb, Adj[Yb], Yb] + 
  18*Yt*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt] + 
  108*MatMul[Yt, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt] + 
  3*ht*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt] + 
  18*MatMul[ht, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]}
