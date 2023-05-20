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

{(-9*g1^2*he)/5 - 3*g2^2*he + (18*g1^2*M1*Ye)/5 + 6*g2^2*M2*Ye + 
  5*MatMul[he, Adj[Ye], Ye] + 4*MatMul[Ye, Adj[Ye], he] + 
  6*Ye*trace[hb, Adj[Yb]] + 2*Ye*trace[he, Adj[Ye]] + 
  3*he*trace[Adj[Yb], Yb] + he*trace[Adj[Ye], Ye], 
 (27*g1^4*he)/2 + (9*g1^2*g2^2*he)/5 + (15*g2^4*he)/2 - 54*g1^4*M1*Ye - 
  (18*g1^2*g2^2*M1*Ye)/5 - (18*g1^2*g2^2*M2*Ye)/5 - 30*g2^4*M2*Ye - 
  (6*g1^2*MatMul[he, Adj[Ye], Ye])/5 + 12*g2^2*MatMul[he, Adj[Ye], Ye] + 
  (6*g1^2*MatMul[Ye, Adj[Ye], he])/5 + 6*g2^2*MatMul[Ye, Adj[Ye], he] - 
  12*g2^2*M2*MatMul[Ye, Adj[Ye], Ye] - 
  6*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye] - 
  8*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye] - 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he] - (4*g1^2*Ye*trace[hb, Adj[Yb]])/5 + 
  32*g3^2*Ye*trace[hb, Adj[Yb]] - 18*MatMul[Ye, Adj[Ye], Ye]*
   trace[hb, Adj[Yb]] + (12*g1^2*Ye*trace[he, Adj[Ye]])/5 - 
  6*MatMul[Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]] - 
  (2*g1^2*he*trace[Adj[Yb], Yb])/5 + 16*g3^2*he*trace[Adj[Yb], Yb] - 
  32*g3^2*M3*Ye*trace[Adj[Yb], Yb] - 15*MatMul[he, Adj[Ye], Ye]*
   trace[Adj[Yb], Yb] - 12*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb] + 
  (6*g1^2*he*trace[Adj[Ye], Ye])/5 - 5*MatMul[he, Adj[Ye], Ye]*
   trace[Adj[Ye], Ye] - 4*MatMul[Ye, Adj[Ye], he]*trace[Adj[Ye], Ye] + 
  M1*((4*g1^2*Ye*trace[Adj[Yb], Yb])/5 - (12*g1^2*Ye*trace[Adj[Ye], Ye])/5) - 
  6*Ye*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
  36*Ye*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
  9*he*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  12*Ye*trace[Adj[Ye], he, Adj[Ye], Ye] - 
  3*he*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  6*Ye*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
  3*he*trace[Adj[Yt], Yt, Adj[Yb], Yb], (24993*g1^6*he)/250 + 
  (837*g1^4*g2^2*he)/50 + (9*g1^2*g2^4*he)/2 + (345*g2^6*he)/2 + 
  (396*g1^4*g3^2*he)/5 + 180*g2^4*g3^2*he - (5373*g1^6*MatMul[he, Zeta[3]])/
   125 - (729*g1^4*g2^2*MatMul[he, Zeta[3]])/25 - 
  (81*g1^2*g2^4*MatMul[he, Zeta[3]])/5 + 315*g2^6*MatMul[he, Zeta[3]] - 
  (2376*g1^4*g3^2*MatMul[he, Zeta[3]])/25 - 
  216*g2^4*g3^2*MatMul[he, Zeta[3]] + 
  M1*((-74979*g1^6*Ye)/125 + (32238*g1^6*MatMul[Ye, Zeta[3]])/125) + 
  M2*((-837*g1^4*g2^2*Ye)/25 + (1458*g1^4*g2^2*MatMul[Ye, Zeta[3]])/25) + 
  M1*((-1674*g1^4*g2^2*Ye)/25 + (2916*g1^4*g2^2*MatMul[Ye, Zeta[3]])/25) + 
  M1*(-9*g1^2*g2^4*Ye + (162*g1^2*g2^4*MatMul[Ye, Zeta[3]])/5) + 
  M2*(-18*g1^2*g2^4*Ye + (324*g1^2*g2^4*MatMul[Ye, Zeta[3]])/5) + 
  M2*(-1035*g2^6*Ye - 1890*g2^6*MatMul[Ye, Zeta[3]]) + 
  M3*((-792*g1^4*g3^2*Ye)/5 + (4752*g1^4*g3^2*MatMul[Ye, Zeta[3]])/25) + 
  M1*((-1584*g1^4*g3^2*Ye)/5 + (9504*g1^4*g3^2*MatMul[Ye, Zeta[3]])/25) + 
  M3*(-360*g2^4*g3^2*Ye + 432*g2^4*g3^2*MatMul[Ye, Zeta[3]]) + 
  M2*(-720*g2^4*g3^2*Ye + 864*g2^4*g3^2*MatMul[Ye, Zeta[3]]) - 
  (1539*g1^4*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]])/50 + 
  81*g1^2*g2^2*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]] - 
  (99*g2^4*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]])/2 - 
  (324*g1^4*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]])/25 + 
  (324*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]])/5 - 
  72*g2^4*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]] + 
  (54*g1^2*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] - 
  (54*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he], Zeta[3]] + 
  30*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  36*MatMul[MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  36*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], he], Zeta[3]] - 
  (8757*g1^4*MatMul[he, Adj[Ye], Ye])/100 - 
  (621*g1^2*g2^2*MatMul[he, Adj[Ye], Ye])/10 - 
  (393*g2^4*MatMul[he, Adj[Ye], Ye])/4 - (2124*g1^4*MatMul[Ye, Adj[Ye], he])/
   25 - (216*g1^2*g2^2*MatMul[Ye, Adj[Ye], he])/5 - 
  66*g2^4*MatMul[Ye, Adj[Ye], he] + 
  M1*((1458*g1^4*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]])/25 + 
    (5751*g1^4*MatMul[Ye, Adj[Ye], Ye])/25) + 
  M1*((-486*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]])/5 + 
    (351*g1^2*g2^2*MatMul[Ye, Adj[Ye], Ye])/5) + 
  M2*((-486*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]])/5 + 
    (351*g1^2*g2^2*MatMul[Ye, Adj[Ye], Ye])/5) + 
  M2*(162*g2^4*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]] + 
    219*g2^4*MatMul[Ye, Adj[Ye], Ye]) + 
  (63*g1^2*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  15*g2^2*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye] + 
  (108*g1^2*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye])/5 + 
  12*g2^2*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye] + 
  (99*g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he])/5 + 
  3*g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he] - 
  (108*g1^2*M1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye])/5 - 
  12*g2^2*M2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], he] - 
  (301*g1^4*Ye*trace[hb, Adj[Yb]])/6 - (3*g1^2*g2^2*Ye*trace[hb, Adj[Yb]])/
   5 - (315*g2^4*Ye*trace[hb, Adj[Yb]])/2 - 
  (568*g1^2*g3^2*Ye*trace[hb, Adj[Yb]])/15 - 
  264*g2^2*g3^2*Ye*trace[hb, Adj[Yb]] - (320*g3^4*Ye*trace[hb, Adj[Yb]])/3 - 
  (77*g1^4*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yb]])/25 - 
  18*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yb]] - 
  189*g2^4*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yb]] + 
  (224*g1^2*g3^2*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yb]])/5 + 
  288*g2^2*g3^2*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yb]] - 
  32*g3^4*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yb]] - 
  (36*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[Yb]])/5 - 
  108*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[Yb]] + 
  288*g3^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[Yb]] + 
  (294*g1^2*MatMul[Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]])/5 + 
  90*g2^2*MatMul[Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]] - 
  192*g3^2*MatMul[Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]] + 
  24*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]] - 
  (873*g1^4*Ye*trace[he, Adj[Ye]])/10 - (81*g1^2*g2^2*Ye*trace[he, Adj[Ye]])/
   5 - (105*g2^4*Ye*trace[he, Adj[Ye]])/2 + 
  (81*g1^4*MatMul[Ye, Zeta[3]]*trace[he, Adj[Ye]])/25 + 
  (162*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[he, Adj[Ye]])/5 - 
  63*g2^4*MatMul[Ye, Zeta[3]]*trace[he, Adj[Ye]] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[he, Adj[Ye]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[he, Adj[Ye]] + 
  (18*g1^2*MatMul[Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]])/5 + 
  30*g2^2*MatMul[Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]] + 
  8*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]] - 
  (234*g1^4*Ye*trace[ht, Adj[Yt]])/5 - 90*g2^4*Ye*trace[ht, Adj[Yt]] - 
  (301*g1^4*he*trace[Adj[Yb], Yb])/12 - (3*g1^2*g2^2*he*trace[Adj[Yb], Yb])/
   10 - (315*g2^4*he*trace[Adj[Yb], Yb])/4 - 
  (284*g1^2*g3^2*he*trace[Adj[Yb], Yb])/15 - 
  132*g2^2*g3^2*he*trace[Adj[Yb], Yb] - (160*g3^4*he*trace[Adj[Yb], Yb])/3 - 
  (77*g1^4*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb])/50 - 
  9*g1^2*g2^2*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (189*g2^4*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb])/2 + 
  (112*g1^2*g3^2*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  144*g2^2*g3^2*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb] - 
  16*g3^4*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb] + 
  (24*g1^2*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  108*g2^2*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
  240*g3^2*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (78*g1^2*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  54*g2^2*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]]*trace[Adj[Yb], Yb] + 
  192*g3^2*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (254*g1^2*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yb], Yb])/5 + 
  72*g2^2*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  160*g3^2*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  (187*g1^2*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb])/5 + 
  63*g2^2*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb] - 
  128*g3^2*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb] + 
  18*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  24*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  18*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he]*trace[Adj[Yb], Yb] - 
  108*MatMul[Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  36*MatMul[Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  45*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yb], Yb]^2 - 
  36*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb]^2 + 
  M1*((568*g1^2*g3^2*Ye*trace[Adj[Yb], Yb])/15 - 
    (224*g1^2*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb])/5) + 
  M3*((568*g1^2*g3^2*Ye*trace[Adj[Yb], Yb])/15 - 
    (224*g1^2*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb])/5) + 
  M2*(264*g2^2*g3^2*Ye*trace[Adj[Yb], Yb] - 288*g2^2*g3^2*MatMul[Ye, Zeta[3]]*
     trace[Adj[Yb], Yb]) + M3*(264*g2^2*g3^2*Ye*trace[Adj[Yb], Yb] - 
    288*g2^2*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M3*((640*g3^4*Ye*trace[Adj[Yb], Yb])/3 + 64*g3^4*MatMul[Ye, Zeta[3]]*
     trace[Adj[Yb], Yb]) + 
  M3*(-288*g3^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
    192*g3^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb]) - 
  (873*g1^4*he*trace[Adj[Ye], Ye])/20 - (81*g1^2*g2^2*he*trace[Adj[Ye], Ye])/
   10 - (105*g2^4*he*trace[Adj[Ye], Ye])/4 + 
  (81*g1^4*MatMul[he, Zeta[3]]*trace[Adj[Ye], Ye])/50 + 
  (81*g1^2*g2^2*MatMul[he, Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  (63*g2^4*MatMul[he, Zeta[3]]*trace[Adj[Ye], Ye])/2 + 
  (108*g1^2*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  36*g2^2*MatMul[MatMul[he, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Ye], Ye] + 
  (54*g1^2*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  18*g2^2*MatMul[MatMul[Ye, Adj[Ye], he], Zeta[3]]*trace[Adj[Ye], Ye] + 
  (18*g1^2*MatMul[he, Adj[Ye], Ye]*trace[Adj[Ye], Ye])/5 + 
  24*g2^2*MatMul[he, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  (9*g1^2*MatMul[Ye, Adj[Ye], he]*trace[Adj[Ye], Ye])/5 + 
  21*g2^2*MatMul[Ye, Adj[Ye], he]*trace[Adj[Ye], Ye] + 
  6*MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  8*MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he]*trace[Adj[Ye], Ye] - 
  36*MatMul[Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  30*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  24*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  5*MatMul[he, Adj[Ye], Ye]*trace[Adj[Ye], Ye]^2 - 
  4*MatMul[Ye, Adj[Ye], he]*trace[Adj[Ye], Ye]^2 + 
  M1*((3*g1^2*g2^2*Ye*trace[Adj[Yb], Yb])/5 + 
    18*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb] + 
    (81*g1^2*g2^2*Ye*trace[Adj[Ye], Ye])/5 - 
    (162*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M2*((3*g1^2*g2^2*Ye*trace[Adj[Yb], Yb])/5 + 
    18*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb] + 
    (81*g1^2*g2^2*Ye*trace[Adj[Ye], Ye])/5 - 
    (162*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M1*((36*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb])/
     5 - (294*g1^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb])/5 - 
    (108*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Ye], Ye])/
     5 - (18*g1^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye])/5) + 
  M2*(108*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] - 
    90*g2^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
    36*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Ye], Ye] - 
    30*g2^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye]) - 
  (117*g1^4*he*trace[Adj[Yt], Yt])/5 - 45*g2^4*he*trace[Adj[Yt], Yt] + 
  M1*((301*g1^4*Ye*trace[Adj[Yb], Yb])/3 + 
    (154*g1^4*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb])/25 + 
    (873*g1^4*Ye*trace[Adj[Ye], Ye])/5 - 
    (162*g1^4*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye])/25 + 
    (468*g1^4*Ye*trace[Adj[Yt], Yt])/5) + 
  M2*(315*g2^4*Ye*trace[Adj[Yb], Yb] + 378*g2^4*MatMul[Ye, Zeta[3]]*
     trace[Adj[Yb], Yb] + 105*g2^4*Ye*trace[Adj[Ye], Ye] + 
    126*g2^4*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye] + 
    180*g2^4*Ye*trace[Adj[Yt], Yt]) - 
  (24*g1^2*Ye*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 + 
  36*g2^2*Ye*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  48*g3^2*Ye*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  (84*g1^2*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 - 
  96*g3^2*MatMul[Ye, Zeta[3]]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  36*MatMul[Ye, Adj[Ye], Ye]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  36*Ye*trace[Adj[Yt], Yt]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*g1^2*Ye*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  36*g2^2*Ye*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  288*g3^2*Ye*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  (216*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], hb, Adj[Yb], Yb])/5 + 
  216*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
  576*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  216*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  216*Ye*trace[Adj[Yb], Yb]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  72*Ye*trace[Adj[Ye], Ye]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  3*g1^2*he*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g2^2*he*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*g3^2*he*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  (54*g1^2*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb])/5 + 
  54*g2^2*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  144*g3^2*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  90*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  108*Ye*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*Ye*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  54*he*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*he*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*g1^2*Ye*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  12*g2^2*Ye*trace[Adj[Ye], he, Adj[Ye], Ye] - 
  (216*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], he, Adj[Ye], Ye])/5 + 
  72*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  72*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  72*Ye*trace[Adj[Yb], Yb]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  24*Ye*trace[Adj[Ye], Ye]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  9*g1^2*he*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  3*g2^2*he*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (54*g1^2*MatMul[he, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  18*g2^2*MatMul[he, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  30*MatMul[he, Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  24*MatMul[Ye, Adj[Ye], he]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  36*Ye*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*Ye*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*he*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*he*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (24*g1^2*Ye*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 + 
  36*g2^2*Ye*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  48*g3^2*Ye*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  (84*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 - 
  96*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  36*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  36*Ye*trace[Adj[Yt], Yt]*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
  (12*g1^2*he*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  18*g2^2*he*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*g3^2*he*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (42*g1^2*MatMul[he, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
  48*g3^2*MatMul[he, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  30*MatMul[he, Adj[Ye], Ye]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*MatMul[Ye, Adj[Ye], he]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*Ye*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*he*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  M2*(-18*g2^2*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    108*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    6*g2^2*Ye*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 36*g2^2*MatMul[Ye, Zeta[3]]*
     trace[Adj[Ye], Ye, Adj[Ye], Ye] - 36*g2^2*Ye*trace[Adj[Yt], Yt, Adj[Yb], 
      Yb]) + M1*(-6*g1^2*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    (108*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb])/5 - 
    18*g1^2*Ye*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    (108*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
    (24*g1^2*Ye*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
    (84*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5) + 
  M3*(-144*g3^2*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    288*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    48*g3^2*Ye*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 96*g3^2*MatMul[Ye, Zeta[3]]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb]) + 
  18*Ye*trace[hb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  18*Ye*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
  108*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] + 
  3*he*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*MatMul[he, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  6*Ye*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
  36*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] + 
  he*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*MatMul[he, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*Ye*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*Ye*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yb], Yb] + 
  9*he*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb]}
