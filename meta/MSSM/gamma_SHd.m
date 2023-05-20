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

{(-3*g1^2)/10 - (3*g2^2)/2 + 3*trace[Adj[Yb], Yb] + trace[Adj[Ye], Ye], 
 (207*g1^4)/100 + (9*g1^2*g2^2)/10 + (15*g2^4)/4 - 
  (2*g1^2*trace[Adj[Yb], Yb])/5 + 16*g3^2*trace[Adj[Yb], Yb] + 
  (6*g1^2*trace[Adj[Ye], Ye])/5 - 9*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  3*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 3*trace[Adj[Yt], Yt, Adj[Yb], Yb], 
 (1839*g1^6)/100 + (27*g1^4*g2^2)/100 + (9*g1^2*g2^4)/4 + (345*g2^6)/4 + 
  (66*g1^4*g3^2)/5 + 90*g2^4*g3^2 - (175*g1^4*trace[Adj[Yb], Yb])/12 - 
  (3*g1^2*g2^2*trace[Adj[Yb], Yb])/10 - (225*g2^4*trace[Adj[Yb], Yb])/4 - 
  (284*g1^2*g3^2*trace[Adj[Yb], Yb])/15 - 132*g2^2*g3^2*trace[Adj[Yb], Yb] - 
  (160*g3^4*trace[Adj[Yb], Yb])/3 - (603*g1^4*trace[Adj[Ye], Ye])/20 - 
  (81*g1^2*g2^2*trace[Adj[Ye], Ye])/10 - (75*g2^4*trace[Adj[Ye], Ye])/4 - 
  (39*g1^4*trace[Adj[Yt], Yt])/10 - (45*g2^4*trace[Adj[Yt], Yt])/2 + 
  3*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  9*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 54*trace[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 18*trace[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 9*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  3*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 18*trace[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*trace[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (12*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  18*g2^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 18*trace[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  3*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  9*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] - (1791*g1^6*Zeta[3])/250 - 
  (243*g1^4*g2^2*Zeta[3])/50 - (81*g1^2*g2^4*Zeta[3])/10 + 
  (315*g2^6*Zeta[3])/2 - (396*g1^4*g3^2*Zeta[3])/25 - 108*g2^4*g3^2*Zeta[3] - 
  (77*g1^4*trace[Adj[Yb], Yb]*Zeta[3])/50 - 9*g1^2*g2^2*trace[Adj[Yb], Yb]*
   Zeta[3] - (189*g2^4*trace[Adj[Yb], Yb]*Zeta[3])/2 + 
  (112*g1^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
  144*g2^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3] - 16*g3^4*trace[Adj[Yb], Yb]*
   Zeta[3] + (81*g1^4*trace[Adj[Ye], Ye]*Zeta[3])/50 + 
  (81*g1^2*g2^2*trace[Adj[Ye], Ye]*Zeta[3])/5 - 
  (63*g2^4*trace[Adj[Ye], Ye]*Zeta[3])/2 + 
  (54*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 + 
  54*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
  144*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
  (54*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 + 
  18*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  (42*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
  48*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] + 
  18*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3]}
