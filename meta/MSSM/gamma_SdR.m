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

{(-2*g1^2)/15 - (8*g3^2)/3 + 2*MatMul[Yb, Adj[Yb]], 
 (202*g1^4)/225 + (32*g1^2*g3^2)/45 - (8*g3^4)/9 + 
  (2*g1^2*MatMul[Yb, Adj[Yb]])/5 + 6*g2^2*MatMul[Yb, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]] - 2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] - 
  6*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 2*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Ye], Ye], (5629*g1^6)/675 + (9*g1^4*g2^2)/5 + 
  (728*g1^4*g3^2)/225 + (452*g1^2*g3^4)/45 + 60*g2^2*g3^4 + (2720*g3^6)/27 - 
  (337*g1^4*MatMul[Yb, Adj[Yb]])/30 - (43*g1^2*g2^2*MatMul[Yb, Adj[Yb]])/5 - 
  (87*g2^4*MatMul[Yb, Adj[Yb]])/2 - (24*g1^2*g3^2*MatMul[Yb, Adj[Yb]])/5 - 
  88*g2^2*g3^2*MatMul[Yb, Adj[Yb]] + (16*g3^4*MatMul[Yb, Adj[Yb]])/3 - 
  (7*g1^4*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/75 + 
  (42*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]] + 
  (32*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3 + 
  96*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]] - 
  (544*g3^4*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3 + 
  (6*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  (18*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] - 
  (g1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  9*g2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]] + 
  (64*g3^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 - 
  (29*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/15 + 
  9*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] + 
  (64*g3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/3 + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] - 
  (14*g1^4*trace[Adj[Yb], Yb])/15 - (80*g3^4*trace[Adj[Yb], Yb])/3 + 
  7*g1^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  27*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  16*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  6*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
  54*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  18*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]^2 - 
  (6*g1^4*trace[Adj[Ye], Ye])/5 - 3*g1^2*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Ye], Ye] + 9*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  16*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  6*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye]^2 - (26*g1^4*trace[Adj[Yt], Yt])/
   15 - (80*g3^4*trace[Adj[Yt], Yt])/3 + 12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   trace[Adj[Yt], Yt] + 36*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], 
    Yb] + 12*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[Yb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  (398*g1^6*Zeta[3])/125 - (54*g1^4*g2^2*Zeta[3])/25 - 
  (176*g1^4*g3^2*Zeta[3])/25 - (88*g1^2*g3^4*Zeta[3])/5 - 
  72*g2^2*g3^4*Zeta[3] + 320*g3^6*Zeta[3]}
