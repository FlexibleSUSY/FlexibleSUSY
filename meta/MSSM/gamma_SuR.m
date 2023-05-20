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

{(-8*g1^2)/15 - (8*g3^2)/3 + 2*MatMul[Yt, Adj[Yt]], 
 (856*g1^4)/225 + (128*g1^2*g3^2)/45 - (8*g3^4)/9 - 
  (2*g1^2*MatMul[Yt, Adj[Yt]])/5 + 6*g2^2*MatMul[Yt, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] - 2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]] - 
  6*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt], 
 (106868*g1^6)/3375 + (36*g1^4*g2^2)/5 + (2144*g1^4*g3^2)/225 - 
  (172*g1^2*g3^4)/45 + 60*g2^2*g3^4 + (2720*g3^6)/27 - 
  (799*g1^4*MatMul[Yt, Adj[Yt]])/50 - (67*g1^2*g2^2*MatMul[Yt, Adj[Yt]])/5 - 
  (87*g2^4*MatMul[Yt, Adj[Yt]])/2 - (8*g1^2*g3^2*MatMul[Yt, Adj[Yt]])/15 - 
  88*g2^2*g3^2*MatMul[Yt, Adj[Yt]] + (16*g3^4*MatMul[Yt, Adj[Yt]])/3 - 
  (247*g1^4*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/75 + 
  (78*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]] - 
  (224*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/15 + 
  96*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]] - 
  (544*g3^4*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/3 + 
  (18*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] + 
  6*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] - 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  (19*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/15 + 
  9*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] + 
  (64*g3^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/3 - 
  (g1^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  9*g2^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]] + 
  (64*g3^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  6*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]] - 
  (56*g1^4*trace[Adj[Yb], Yb])/15 - (80*g3^4*trace[Adj[Yb], Yb])/3 + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yb], Yb] - 
  (24*g1^4*trace[Adj[Ye], Ye])/5 + 4*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*
   trace[Adj[Ye], Ye] - (104*g1^4*trace[Adj[Yt], Yt])/15 - 
  (80*g3^4*trace[Adj[Yt], Yt])/3 + 7*g1^2*MatMul[Yt, Adj[Yt]]*
   trace[Adj[Yt], Yt] + 27*g2^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  16*g3^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  (42*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  54*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] + 
  96*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] - 
  6*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  18*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt]^2 + 
  12*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 
  (1592*g1^6*Zeta[3])/125 - (216*g1^4*g2^2*Zeta[3])/25 - 
  (704*g1^4*g3^2*Zeta[3])/25 - (88*g1^2*g3^4*Zeta[3])/5 - 
  72*g2^2*g3^4*Zeta[3] + 320*g3^6*Zeta[3]}
