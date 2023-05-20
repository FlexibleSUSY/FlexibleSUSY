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

{-g1^2/30 - (3*g2^2)/2 - (8*g3^2)/3 + MatMul[Adj[Yb], Yb] + 
  MatMul[Adj[Yt], Yt], (199*g1^4)/900 + (g1^2*g2^2)/10 + (15*g2^4)/4 + 
  (8*g1^2*g3^2)/45 + 8*g2^2*g3^2 - (8*g3^4)/9 + 
  (2*g1^2*MatMul[Adj[Yb], Yb])/5 + (4*g1^2*MatMul[Adj[Yt], Yt])/5 - 
  2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] - 2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] - 
  3*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] - MatMul[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye] - 3*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt], 
 (28457*g1^6)/13500 + (11*g1^4*g2^2)/100 + (25*g1^2*g2^4)/4 + (345*g2^6)/4 + 
  (194*g1^4*g3^2)/225 - (8*g1^2*g2^2*g3^2)/5 + 50*g2^4*g3^2 + 
  (608*g1^2*g3^4)/45 + 8*g2^2*g3^4 + (2720*g3^6)/27 - 
  (633*g1^4*MatMul[Adj[Yb], Yb])/100 - (41*g1^2*g2^2*MatMul[Adj[Yb], Yb])/
   10 - (45*g2^4*MatMul[Adj[Yb], Yb])/4 - (76*g1^2*g3^2*MatMul[Adj[Yb], Yb])/
   15 - 4*g2^2*g3^2*MatMul[Adj[Yb], Yb] + (8*g3^4*MatMul[Adj[Yb], Yb])/3 - 
  (3767*g1^4*MatMul[Adj[Yt], Yt])/300 - (59*g1^2*g2^2*MatMul[Adj[Yt], Yt])/
   10 - (45*g2^4*MatMul[Adj[Yt], Yt])/4 - (68*g1^2*g3^2*MatMul[Adj[Yt], Yt])/
   5 - 4*g2^2*g3^2*MatMul[Adj[Yt], Yt] + (8*g3^4*MatMul[Adj[Yt], Yt])/3 + 
  (7*g1^4*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/30 + 
  (9*g1^2*g2^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/5 - 
  (63*g2^4*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/3 + 
  (143*g1^4*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/150 + 
  9*g1^2*g2^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]] - 
  (63*g2^4*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/3 - 
  (6*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] - 
  6*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  18*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  6*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  6*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  (7*g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/15 - 
  3*g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] + 
  (64*g3^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (11*g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 - 
  3*g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] + 
  (64*g3^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  4*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  4*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] - 
  (7*g1^4*trace[Adj[Yb], Yb])/30 - (45*g2^4*trace[Adj[Yb], Yb])/2 - 
  (80*g3^4*trace[Adj[Yb], Yb])/3 + 
  (16*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 + 
  18*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  8*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  (24*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  48*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  6*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  9*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]^2 - (3*g1^4*trace[Adj[Ye], Ye])/
   10 - (15*g2^4*trace[Adj[Ye], Ye])/2 - 
  (8*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
  6*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  8*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  (12*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
  2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  6*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]^2 - (13*g1^4*trace[Adj[Yt], Yt])/
   30 - (45*g2^4*trace[Adj[Yt], Yt])/2 - (80*g3^4*trace[Adj[Yt], Yt])/3 + 
  2*g1^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  8*g3^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  (24*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  48*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] + 
  6*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  9*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt]^2 + 
  18*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  6*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  6*MatMul[Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  6*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 
  (199*g1^6*Zeta[3])/250 - (27*g1^4*g2^2*Zeta[3])/50 - 
  (81*g1^2*g2^4*Zeta[3])/10 + (315*g2^6*Zeta[3])/2 - 
  (44*g1^4*g3^2*Zeta[3])/25 - 108*g2^4*g3^2*Zeta[3] - 
  (88*g1^2*g3^4*Zeta[3])/5 - 72*g2^2*g3^4*Zeta[3] + 320*g3^6*Zeta[3]}
