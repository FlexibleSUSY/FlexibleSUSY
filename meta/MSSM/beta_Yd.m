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

{(-7*g1^2*Yb)/15 - 3*g2^2*Yb - (16*g3^2*Yb)/3 + 3*MatMul[Yb, Adj[Yb], Yb] + 
  MatMul[Yb, Adj[Yt], Yt] + 3*Yb*trace[Adj[Yb], Yb] + Yb*trace[Adj[Ye], Ye], 
 (287*g1^4*Yb)/90 + g1^2*g2^2*Yb + (15*g2^4*Yb)/2 + (8*g1^2*g3^2*Yb)/9 + 
  8*g2^2*g3^2*Yb - (16*g3^4*Yb)/9 + (4*g1^2*MatMul[Yb, Adj[Yb], Yb])/5 + 
  6*g2^2*MatMul[Yb, Adj[Yb], Yb] + (4*g1^2*MatMul[Yb, Adj[Yt], Yt])/5 - 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt] - (2*g1^2*Yb*trace[Adj[Yb], Yb])/5 + 
  16*g3^2*Yb*trace[Adj[Yb], Yb] - 9*MatMul[Yb, Adj[Yb], Yb]*
   trace[Adj[Yb], Yb] + (6*g1^2*Yb*trace[Adj[Ye], Ye])/5 - 
  3*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  3*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  9*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  3*Yb*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  3*Yb*trace[Adj[Yt], Yt, Adj[Yb], Yb], (194651*g1^6*Yb)/6750 + 
  (109*g1^4*g2^2*Yb)/50 + (17*g1^2*g2^4*Yb)/2 + (345*g2^6*Yb)/2 + 
  (3892*g1^4*g3^2*Yb)/225 - (8*g1^2*g2^2*g3^2*Yb)/5 + 140*g2^4*g3^2*Yb + 
  (212*g1^2*g3^4*Yb)/9 + 68*g2^2*g3^4*Yb + (5440*g3^6*Yb)/27 - 
  (1393*g1^6*MatMul[Yb, Zeta[3]])/125 - (189*g1^4*g2^2*MatMul[Yb, Zeta[3]])/
   25 - (81*g1^2*g2^4*MatMul[Yb, Zeta[3]])/5 + 315*g2^6*MatMul[Yb, Zeta[3]] - 
  (616*g1^4*g3^2*MatMul[Yb, Zeta[3]])/25 - 
  216*g2^4*g3^2*MatMul[Yb, Zeta[3]] - (176*g1^2*g3^4*MatMul[Yb, Zeta[3]])/5 - 
  144*g2^2*g3^4*MatMul[Yb, Zeta[3]] + 640*g3^6*MatMul[Yb, Zeta[3]] + 
  (7*g1^4*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/50 + 
  (51*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/5 - 
  (81*g2^4*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/2 + 
  (96*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
  96*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]] - 
  272*g3^4*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]] + 
  (143*g1^4*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/150 + 
  9*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]] - 
  (63*g2^4*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]])/3 + 
  (18*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb], Zeta[3]] - 
  6*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  18*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  6*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
  (5269*g1^4*MatMul[Yb, Adj[Yb], Yb])/300 - 
  (127*g1^2*g2^2*MatMul[Yb, Adj[Yb], Yb])/10 - 
  (219*g2^4*MatMul[Yb, Adj[Yb], Yb])/4 - 
  (148*g1^2*g3^2*MatMul[Yb, Adj[Yb], Yb])/15 - 
  92*g2^2*g3^2*MatMul[Yb, Adj[Yb], Yb] + 8*g3^4*MatMul[Yb, Adj[Yb], Yb] - 
  (3767*g1^4*MatMul[Yb, Adj[Yt], Yt])/300 - 
  (59*g1^2*g2^2*MatMul[Yb, Adj[Yt], Yt])/10 - 
  (45*g2^4*MatMul[Yb, Adj[Yt], Yt])/4 - 
  (68*g1^2*g3^2*MatMul[Yb, Adj[Yt], Yt])/5 - 
  4*g2^2*g3^2*MatMul[Yb, Adj[Yt], Yt] + (8*g3^4*MatMul[Yb, Adj[Yt], Yt])/3 + 
  (2*g1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb])/15 + 
  6*g2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb])/3 - 
  (29*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb])/15 + 
  9*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  (64*g3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb])/3 + 
  (11*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt])/3 - 
  3*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt] + 
  (64*g3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb], Yb] + 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] - 
  (63*g1^4*Yb*trace[Adj[Yb], Yb])/4 - (3*g1^2*g2^2*Yb*trace[Adj[Yb], Yb])/
   10 - (315*g2^4*Yb*trace[Adj[Yb], Yb])/4 - 
  (284*g1^2*g3^2*Yb*trace[Adj[Yb], Yb])/15 - 
  132*g2^2*g3^2*Yb*trace[Adj[Yb], Yb] - (320*g3^4*Yb*trace[Adj[Yb], Yb])/3 - 
  (77*g1^4*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb])/50 - 
  9*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (189*g2^4*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb])/2 + 
  (112*g1^2*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  144*g2^2*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb] - 
  16*g3^4*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (54*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  54*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  144*g3^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (51*g1^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 + 
  45*g2^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  24*g3^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  27*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb]^2 - 
  (633*g1^4*Yb*trace[Adj[Ye], Ye])/20 - (81*g1^2*g2^2*Yb*trace[Adj[Ye], Ye])/
   10 - (105*g2^4*Yb*trace[Adj[Ye], Ye])/4 + 
  (81*g1^4*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye])/50 + 
  (81*g1^2*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  (63*g2^4*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye])/2 + 
  (42*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye] - 
  (23*g1^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
  15*g2^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  24*g3^2*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  18*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  3*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye]^2 - 
  (91*g1^4*Yb*trace[Adj[Yt], Yt])/15 - 45*g2^4*Yb*trace[Adj[Yt], Yt] - 
  (160*g3^4*Yb*trace[Adj[Yt], Yt])/3 - 
  (24*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  48*g3^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] + 
  2*g1^2*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  8*g3^2*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb]*trace[Adj[Yt], Yt] + 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  9*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt]^2 + 
  3*g1^2*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g2^2*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*g3^2*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  (54*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb])/5 + 
  54*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  144*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  54*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  54*Yb*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*Yb*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g1^2*Yb*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  3*g2^2*Yb*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (54*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  18*g2^2*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*Yb*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*Yb*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (12*g1^2*Yb*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  18*g2^2*Yb*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*g3^2*Yb*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (42*g1^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
  48*g3^2*MatMul[Yb, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*MatMul[Yb, Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  6*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*Yb*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*MatMul[Yb, Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  3*Yb*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*MatMul[Yb, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  Yb*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*MatMul[Yb, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  9*Yb*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb]}
