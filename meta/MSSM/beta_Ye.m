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

{(-9*g1^2*Ye)/5 - 3*g2^2*Ye + 3*MatMul[Ye, Adj[Ye], Ye] + 
  3*Ye*trace[Adj[Yb], Yb] + Ye*trace[Adj[Ye], Ye], 
 (27*g1^4*Ye)/2 + (9*g1^2*g2^2*Ye)/5 + (15*g2^4*Ye)/2 + 
  6*g2^2*MatMul[Ye, Adj[Ye], Ye] - 4*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] - 
  (2*g1^2*Ye*trace[Adj[Yb], Yb])/5 + 16*g3^2*Ye*trace[Adj[Yb], Yb] - 
  9*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  (6*g1^2*Ye*trace[Adj[Ye], Ye])/5 - 3*MatMul[Ye, Adj[Ye], Ye]*
   trace[Adj[Ye], Ye] - 9*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  3*Ye*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  3*Ye*trace[Adj[Yt], Yt, Adj[Yb], Yb], (24993*g1^6*Ye)/250 + 
  (837*g1^4*g2^2*Ye)/50 + (9*g1^2*g2^4*Ye)/2 + (345*g2^6*Ye)/2 + 
  (396*g1^4*g3^2*Ye)/5 + 180*g2^4*g3^2*Ye - (5373*g1^6*MatMul[Ye, Zeta[3]])/
   125 - (729*g1^4*g2^2*MatMul[Ye, Zeta[3]])/25 - 
  (81*g1^2*g2^4*MatMul[Ye, Zeta[3]])/5 + 315*g2^6*MatMul[Ye, Zeta[3]] - 
  (2376*g1^4*g3^2*MatMul[Ye, Zeta[3]])/25 - 
  216*g2^4*g3^2*MatMul[Ye, Zeta[3]] - 
  (729*g1^4*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]])/50 + 
  (243*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]])/5 - 
  (81*g2^4*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]])/2 + 
  18*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] - 
  (5751*g1^4*MatMul[Ye, Adj[Ye], Ye])/100 - 
  (351*g1^2*g2^2*MatMul[Ye, Adj[Ye], Ye])/10 - 
  (219*g2^4*MatMul[Ye, Adj[Ye], Ye])/4 + 
  (54*g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  6*g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] - 
  (301*g1^4*Ye*trace[Adj[Yb], Yb])/12 - (3*g1^2*g2^2*Ye*trace[Adj[Yb], Yb])/
   10 - (315*g2^4*Ye*trace[Adj[Yb], Yb])/4 - 
  (284*g1^2*g3^2*Ye*trace[Adj[Yb], Yb])/15 - 
  132*g2^2*g3^2*Ye*trace[Adj[Yb], Yb] - (160*g3^4*Ye*trace[Adj[Yb], Yb])/3 - 
  (77*g1^4*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb])/50 - 
  9*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (189*g2^4*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb])/2 + 
  (112*g1^2*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  144*g2^2*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb] - 
  16*g3^4*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb] - 
  (18*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  54*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
  144*g3^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (147*g1^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb])/5 + 
  45*g2^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  96*g3^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  27*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb]^2 - 
  (873*g1^4*Ye*trace[Adj[Ye], Ye])/20 - (81*g1^2*g2^2*Ye*trace[Adj[Ye], Ye])/
   10 - (105*g2^4*Ye*trace[Adj[Ye], Ye])/4 + 
  (81*g1^4*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye])/50 + 
  (81*g1^2*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  (63*g2^4*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye])/2 + 
  (54*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  18*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Ye], Ye] + 
  (9*g1^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye])/5 + 
  15*g2^2*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  4*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye] - 
  18*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  3*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye]^2 - 
  (117*g1^4*Ye*trace[Adj[Yt], Yt])/5 - 45*g2^4*Ye*trace[Adj[Yt], Yt] + 
  3*g1^2*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g2^2*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*g3^2*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  (54*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb])/5 + 
  54*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  144*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  54*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  54*Ye*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*Ye*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g1^2*Ye*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  3*g2^2*Ye*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (54*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye])/5 + 
  18*g2^2*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*Ye*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*Ye*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (12*g1^2*Ye*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  18*g2^2*Ye*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*g3^2*Ye*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (42*g1^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 - 
  48*g3^2*MatMul[Ye, Zeta[3]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*MatMul[Ye, Adj[Ye], Ye]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*Ye*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  3*Ye*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  18*MatMul[Ye, Zeta[3]]*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  Ye*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*MatMul[Ye, Zeta[3]]*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  9*Ye*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb]}
