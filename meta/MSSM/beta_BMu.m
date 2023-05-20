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

{6*Mu*trace[Adj[Yb], hb] + 2*Mu*trace[Adj[Ye], he] + 
  (6*Mu*(g1^2*M1 + 5*g2^2*M2 + 5*trace[Adj[Yt], ht]))/5 + 
  BMu*((-3*g1^2)/5 - 3*g2^2 + 3*trace[Adj[Yb], Yb] + trace[Adj[Ye], Ye] + 
    3*trace[Adj[Yt], Yt]), 
 (-2*Mu*(207*g1^4*M1 + 45*g1^2*g2^2*M1 + 45*g1^2*g2^2*M2 + 375*g2^4*M2 + 
     10*(g1^2 - 40*g3^2)*trace[Adj[Yb], hb] - 10*(g1^2*M1 - 40*g3^2*M3)*
      trace[Adj[Yb], Yb] - 30*g1^2*trace[Adj[Ye], he] + 
     30*g1^2*M1*trace[Adj[Ye], Ye] - 20*g1^2*trace[Adj[Yt], ht] - 
     400*g3^2*trace[Adj[Yt], ht] + 20*g1^2*M1*trace[Adj[Yt], Yt] + 
     400*g3^2*M3*trace[Adj[Yt], Yt] + 225*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
     225*trace[Adj[Yb], Yb, Adj[Yb], hb] + 
     75*trace[Adj[Ye], he, Adj[Ye], Ye] + 
     75*trace[Adj[Ye], Ye, Adj[Ye], he] + 
     150*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
     225*trace[Adj[Yt], ht, Adj[Yt], Yt] + 
     150*trace[Adj[Yt], Yt, Adj[Yb], hb] + 
     225*trace[Adj[Yt], Yt, Adj[Yt], ht]))/25 + 
  BMu*((207*g1^4)/50 + (9*g1^2*g2^2)/5 + (15*g2^4)/2 - 
    (2*g1^2*trace[Adj[Yb], Yb])/5 + 16*g3^2*trace[Adj[Yb], Yb] + 
    (6*g1^2*trace[Adj[Ye], Ye])/5 + (4*g1^2*trace[Adj[Yt], Yt])/5 + 
    16*g3^2*trace[Adj[Yt], Yt] - 9*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
    3*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 6*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
    9*trace[Adj[Yt], Yt, Adj[Yt], Yt]), 
 BMu*((-1001*g1^4*trace[Adj[Yb], Yb])/60 - (3*g1^2*g2^2*trace[Adj[Yb], Yb])/
     10 - (315*g2^4*trace[Adj[Yb], Yb])/4 - 
    (284*g1^2*g3^2*trace[Adj[Yb], Yb])/15 - 132*g2^2*g3^2*
     trace[Adj[Yb], Yb] - (160*g3^4*trace[Adj[Yb], Yb])/3 - 
    (657*g1^4*trace[Adj[Ye], Ye])/20 - (81*g1^2*g2^2*trace[Adj[Ye], Ye])/10 - 
    (105*g2^4*trace[Adj[Ye], Ye])/4 - (2357*g1^4*trace[Adj[Yt], Yt])/60 - 
    (57*g1^2*g2^2*trace[Adj[Yt], Yt])/10 - (315*g2^4*trace[Adj[Yt], Yt])/4 - 
    (124*g1^2*g3^2*trace[Adj[Yt], Yt])/3 - 132*g2^2*g3^2*trace[Adj[Yt], Yt] - 
    (160*g3^4*trace[Adj[Yt], Yt])/3 + 3*g1^2*trace[Adj[Yb], Yb, Adj[Yb], 
      Yb] + 9*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    72*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 54*trace[Adj[Yb], Yb]*
     trace[Adj[Yb], Yb, Adj[Yb], Yb] + 18*trace[Adj[Ye], Ye]*
     trace[Adj[Yb], Yb, Adj[Yb], Yb] + 9*g1^2*trace[Adj[Ye], Ye, Adj[Ye], 
      Ye] + 3*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    18*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    6*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
    (6*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
    36*g2^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    48*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 18*trace[Adj[Yb], Yb]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*trace[Adj[Ye], Ye]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb] + 18*trace[Adj[Yt], Yt]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    (57*g1^2*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
    9*g2^2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
    72*g3^2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 54*trace[Adj[Yt], Yt]*
     trace[Adj[Yt], Yt, Adj[Yt], Yt] + 3*trace[Adj[Yb], Yb, Adj[Yb], Yb, 
      Adj[Yb], Yb] + 9*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
    trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
    9*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] + 
    3*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt] - 
    (77*g1^4*trace[Adj[Yb], Yb]*Zeta[3])/50 - 9*g1^2*g2^2*trace[Adj[Yb], Yb]*
     Zeta[3] - (189*g2^4*trace[Adj[Yb], Yb]*Zeta[3])/2 + 
    (112*g1^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
    144*g2^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3] - 16*g3^4*trace[Adj[Yb], Yb]*
     Zeta[3] + (81*g1^4*trace[Adj[Ye], Ye]*Zeta[3])/50 + 
    (81*g1^2*g2^2*trace[Adj[Ye], Ye]*Zeta[3])/5 - 
    (63*g2^4*trace[Adj[Ye], Ye]*Zeta[3])/2 - 
    (13*g1^4*trace[Adj[Yt], Yt]*Zeta[3])/10 + 
    (63*g1^2*g2^2*trace[Adj[Yt], Yt]*Zeta[3])/5 - 
    (189*g2^4*trace[Adj[Yt], Yt]*Zeta[3])/2 + 
    (208*g1^2*g3^2*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
    144*g2^2*g3^2*trace[Adj[Yt], Yt]*Zeta[3] - 16*g3^4*trace[Adj[Yt], Yt]*
     Zeta[3] + (54*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 + 
    54*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
    144*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
    (54*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 + 
    18*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
    (48*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
    96*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] - 
    (18*g1^2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3])/5 + 
    54*g2^2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] - 
    144*g3^2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
    18*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
    18*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
    2*((1839*g1^6)/100 - (1791*g1^6*Zeta[3])/250) + 
    2*((27*g1^4*g2^2)/100 - (243*g1^4*g2^2*Zeta[3])/50) + 
    2*((9*g1^2*g2^4)/4 - (81*g1^2*g2^4*Zeta[3])/10) + 
    2*((345*g2^6)/4 + (315*g2^6*Zeta[3])/2) + 
    2*((66*g1^4*g3^2)/5 - (396*g1^4*g3^2*Zeta[3])/25) + 
    2*(90*g2^4*g3^2 - 108*g2^4*g3^2*Zeta[3])) - 
  2*Mu*((5517*g1^6*M1)/50 + (27*g1^4*g2^2*M1)/25 + (9*g1^2*g2^4*M1)/2 + 
    (264*g1^4*g3^2*M1)/5 + (27*g1^4*g2^2*M2)/50 + 9*g1^2*g2^4*M2 + 
    (1035*g2^6*M2)/2 + 360*g2^4*g3^2*M2 + (132*g1^4*g3^2*M3)/5 + 
    180*g2^4*g3^2*M3 + (657*g1^4*trace[Adj[Ye], he])/20 + 
    (81*g1^2*g2^2*trace[Adj[Ye], he])/10 + (105*g2^4*trace[Adj[Ye], he])/4 - 
    (657*g1^4*M1*trace[Adj[Ye], Ye])/10 - 
    (81*g1^2*g2^2*M1*trace[Adj[Ye], Ye])/10 - 
    (81*g1^2*g2^2*M2*trace[Adj[Ye], Ye])/10 - 
    (105*g2^4*M2*trace[Adj[Ye], Ye])/2 + (2357*g1^4*trace[Adj[Yt], ht])/60 + 
    (57*g1^2*g2^2*trace[Adj[Yt], ht])/10 + (315*g2^4*trace[Adj[Yt], ht])/4 + 
    (124*g1^2*g3^2*trace[Adj[Yt], ht])/3 + 132*g2^2*g3^2*trace[Adj[Yt], ht] + 
    (160*g3^4*trace[Adj[Yt], ht])/3 - (2357*g1^4*M1*trace[Adj[Yt], Yt])/30 - 
    (57*g1^2*g2^2*M1*trace[Adj[Yt], Yt])/10 - 
    (124*g1^2*g3^2*M1*trace[Adj[Yt], Yt])/3 - 
    (57*g1^2*g2^2*M2*trace[Adj[Yt], Yt])/10 - 
    (315*g2^4*M2*trace[Adj[Yt], Yt])/2 - 132*g2^2*g3^2*M2*
     trace[Adj[Yt], Yt] - (124*g1^2*g3^2*M3*trace[Adj[Yt], Yt])/3 - 
    132*g2^2*g3^2*M3*trace[Adj[Yt], Yt] - (320*g3^4*M3*trace[Adj[Yt], Yt])/
     3 - 3*g1^2*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
    9*g2^2*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
    72*g3^2*trace[Adj[Yb], hb, Adj[Yb], Yb] - 18*trace[Adj[Ye], Ye]*
     trace[Adj[Yb], hb, Adj[Yb], Yb] - 3*g1^2*trace[Adj[Yb], Yb, Adj[Yb], 
      hb] - 9*g2^2*trace[Adj[Yb], Yb, Adj[Yb], hb] - 
    72*g3^2*trace[Adj[Yb], Yb, Adj[Yb], hb] - 18*trace[Adj[Ye], Ye]*
     trace[Adj[Yb], Yb, Adj[Yb], hb] + 3*g1^2*M1*trace[Adj[Yb], Yb, Adj[Yb], 
      Yb] + 9*g2^2*M2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    72*g3^2*M3*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 18*trace[Adj[Ye], he]*
     trace[Adj[Yb], Yb, Adj[Yb], Yb] - 9*g1^2*trace[Adj[Ye], he, Adj[Ye], 
      Ye] - 3*g2^2*trace[Adj[Ye], he, Adj[Ye], Ye] - 
    6*trace[Adj[Ye], Ye]*trace[Adj[Ye], he, Adj[Ye], Ye] - 
    9*g1^2*trace[Adj[Ye], Ye, Adj[Ye], he] - 
    3*g2^2*trace[Adj[Ye], Ye, Adj[Ye], he] - 6*trace[Adj[Ye], Ye]*
     trace[Adj[Ye], Ye, Adj[Ye], he] + 9*g1^2*M1*trace[Adj[Ye], Ye, Adj[Ye], 
      Ye] + 3*g2^2*M2*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
    6*trace[Adj[Ye], he]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    (6*g1^2*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 - 
    36*g2^2*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
    48*g3^2*trace[Adj[Yt], ht, Adj[Yb], Yb] - 6*trace[Adj[Ye], Ye]*
     trace[Adj[Yt], ht, Adj[Yb], Yb] - 18*trace[Adj[Yt], Yt]*
     trace[Adj[Yt], ht, Adj[Yb], Yb] - 
    (57*g1^2*trace[Adj[Yt], ht, Adj[Yt], Yt])/5 - 
    9*g2^2*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
    72*g3^2*trace[Adj[Yt], ht, Adj[Yt], Yt] - 54*trace[Adj[Yt], Yt]*
     trace[Adj[Yt], ht, Adj[Yt], Yt] + 
    (6*g1^2*trace[Adj[Yt], Yt, Adj[Yb], hb])/5 - 
    36*g2^2*trace[Adj[Yt], Yt, Adj[Yb], hb] - 
    48*g3^2*trace[Adj[Yt], Yt, Adj[Yb], hb] - 6*trace[Adj[Ye], Ye]*
     trace[Adj[Yt], Yt, Adj[Yb], hb] - 18*trace[Adj[Yt], Yt]*
     trace[Adj[Yt], Yt, Adj[Yb], hb] - 
    (6*g1^2*M1*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
    36*g2^2*M2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    48*g3^2*M3*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 6*trace[Adj[Ye], he]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb] - 18*trace[Adj[Yt], ht]*
     trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
    (57*g1^2*trace[Adj[Yt], Yt, Adj[Yt], ht])/5 - 
    9*g2^2*trace[Adj[Yt], Yt, Adj[Yt], ht] - 
    72*g3^2*trace[Adj[Yt], Yt, Adj[Yt], ht] - 54*trace[Adj[Yt], Yt]*
     trace[Adj[Yt], Yt, Adj[Yt], ht] + 
    (57*g1^2*M1*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
    9*g2^2*M2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
    72*g3^2*M3*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 54*trace[Adj[Yt], ht]*
     trace[Adj[Yt], Yt, Adj[Yt], Yt] - 3*trace[Adj[Yb], hb, Adj[Yb], Yb, 
      Adj[Yb], Yb] - 9*trace[Adj[Yb], hb, Adj[Yt], Yt, Adj[Yb], Yb] - 
    3*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb] - 
    3*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], hb] - 
    9*trace[Adj[Yb], Yb, Adj[Yt], ht, Adj[Yb], Yb] - 
    9*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], hb] - 
    trace[Adj[Ye], he, Adj[Ye], Ye, Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye] - 
    trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], he] - 
    9*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[Yb], Yb] - 
    3*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[Yt], Yt] - 
    9*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yb], Yb] - 
    3*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt] - 
    9*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], hb] - 
    3*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], ht] - 
    (5373*g1^6*M1*Zeta[3])/125 - (486*g1^4*g2^2*M1*Zeta[3])/25 - 
    (81*g1^2*g2^4*M1*Zeta[3])/5 - (1584*g1^4*g3^2*M1*Zeta[3])/25 - 
    (243*g1^4*g2^2*M2*Zeta[3])/25 - (162*g1^2*g2^4*M2*Zeta[3])/5 + 
    945*g2^6*M2*Zeta[3] - 432*g2^4*g3^2*M2*Zeta[3] - 
    (792*g1^4*g3^2*M3*Zeta[3])/25 - 216*g2^4*g3^2*M3*Zeta[3] - 
    (81*g1^4*trace[Adj[Ye], he]*Zeta[3])/50 - 
    (81*g1^2*g2^2*trace[Adj[Ye], he]*Zeta[3])/5 + 
    (63*g2^4*trace[Adj[Ye], he]*Zeta[3])/2 + 
    (81*g1^4*M1*trace[Adj[Ye], Ye]*Zeta[3])/25 + 
    (81*g1^2*g2^2*M1*trace[Adj[Ye], Ye]*Zeta[3])/5 + 
    (81*g1^2*g2^2*M2*trace[Adj[Ye], Ye]*Zeta[3])/5 - 
    63*g2^4*M2*trace[Adj[Ye], Ye]*Zeta[3] + 
    (13*g1^4*trace[Adj[Yt], ht]*Zeta[3])/10 - 
    (63*g1^2*g2^2*trace[Adj[Yt], ht]*Zeta[3])/5 + 
    (189*g2^4*trace[Adj[Yt], ht]*Zeta[3])/2 - 
    (208*g1^2*g3^2*trace[Adj[Yt], ht]*Zeta[3])/5 - 
    144*g2^2*g3^2*trace[Adj[Yt], ht]*Zeta[3] + 16*g3^4*trace[Adj[Yt], ht]*
     Zeta[3] - (13*g1^4*M1*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
    (63*g1^2*g2^2*M1*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
    (208*g1^2*g3^2*M1*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
    (63*g1^2*g2^2*M2*trace[Adj[Yt], Yt]*Zeta[3])/5 - 
    189*g2^4*M2*trace[Adj[Yt], Yt]*Zeta[3] + 144*g2^2*g3^2*M2*
     trace[Adj[Yt], Yt]*Zeta[3] + (208*g1^2*g3^2*M3*trace[Adj[Yt], Yt]*
      Zeta[3])/5 + 144*g2^2*g3^2*M3*trace[Adj[Yt], Yt]*Zeta[3] - 
    32*g3^4*M3*trace[Adj[Yt], Yt]*Zeta[3] - 
    (54*g1^2*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3])/5 - 
    54*g2^2*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] + 
    144*g3^2*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] - 
    (54*g1^2*trace[Adj[Yb], Yb, Adj[Yb], hb]*Zeta[3])/5 - 
    54*g2^2*trace[Adj[Yb], Yb, Adj[Yb], hb]*Zeta[3] + 
    144*g3^2*trace[Adj[Yb], Yb, Adj[Yb], hb]*Zeta[3] + 
    (54*g1^2*M1*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 + 
    54*g2^2*M2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
    144*g3^2*M3*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
    (54*g1^2*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3])/5 - 
    18*g2^2*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3] + 
    (54*g1^2*trace[Adj[Ye], Ye, Adj[Ye], he]*Zeta[3])/5 - 
    18*g2^2*trace[Adj[Ye], Ye, Adj[Ye], he]*Zeta[3] - 
    (54*g1^2*M1*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 + 
    18*g2^2*M2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] - 
    (48*g1^2*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3])/5 + 
    96*g3^2*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3] + 
    (18*g1^2*trace[Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3])/5 - 
    54*g2^2*trace[Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] + 
    144*g3^2*trace[Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] - 
    (48*g1^2*trace[Adj[Yt], Yt, Adj[Yb], hb]*Zeta[3])/5 + 
    96*g3^2*trace[Adj[Yt], Yt, Adj[Yb], hb]*Zeta[3] + 
    (48*g1^2*M1*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
    96*g3^2*M3*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] + 
    (18*g1^2*trace[Adj[Yt], Yt, Adj[Yt], ht]*Zeta[3])/5 - 
    54*g2^2*trace[Adj[Yt], Yt, Adj[Yt], ht]*Zeta[3] + 
    144*g3^2*trace[Adj[Yt], Yt, Adj[Yt], ht]*Zeta[3] - 
    (18*g1^2*M1*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3])/5 + 
    54*g2^2*M2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] - 
    144*g3^2*M3*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] - 
    18*trace[Adj[Yb], hb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
    18*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] - 
    18*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], hb]*Zeta[3] - 
    6*trace[Adj[Ye], he, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] - 
    6*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[Ye], Ye]*Zeta[3] - 
    6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], he]*Zeta[3] - 
    18*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] - 
    18*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] - 
    18*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], ht]*Zeta[3] + 
    trace[Adj[Yb], hb]*((1001*g1^4)/60 + (3*g1^2*g2^2)/10 + (315*g2^4)/4 + 
      (284*g1^2*g3^2)/15 + 132*g2^2*g3^2 + (160*g3^4)/3 - 
      54*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 18*trace[Adj[Ye], Ye, Adj[Ye], 
        Ye] - 18*trace[Adj[Yt], Yt, Adj[Yb], Yb] + (77*g1^4*Zeta[3])/50 + 
      9*g1^2*g2^2*Zeta[3] + (189*g2^4*Zeta[3])/2 - (112*g1^2*g3^2*Zeta[3])/
       5 - 144*g2^2*g3^2*Zeta[3] + 16*g3^4*Zeta[3]) + 
    trace[Adj[Yb], Yb]*((-1001*g1^4*M1)/30 - (3*g1^2*g2^2*M1)/10 - 
      (284*g1^2*g3^2*M1)/15 - (3*g1^2*g2^2*M2)/10 - (315*g2^4*M2)/2 - 
      132*g2^2*g3^2*M2 - (284*g1^2*g3^2*M3)/15 - 132*g2^2*g3^2*M3 - 
      (320*g3^4*M3)/3 - 54*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
      54*trace[Adj[Yb], Yb, Adj[Yb], hb] - 18*trace[Adj[Ye], he, Adj[Ye], 
        Ye] - 18*trace[Adj[Ye], Ye, Adj[Ye], he] - 
      18*trace[Adj[Yt], ht, Adj[Yb], Yb] - 18*trace[Adj[Yt], Yt, Adj[Yb], 
        hb] - (77*g1^4*M1*Zeta[3])/25 - 9*g1^2*g2^2*M1*Zeta[3] + 
      (112*g1^2*g3^2*M1*Zeta[3])/5 - 9*g1^2*g2^2*M2*Zeta[3] - 
      189*g2^4*M2*Zeta[3] + 144*g2^2*g3^2*M2*Zeta[3] + 
      (112*g1^2*g3^2*M3*Zeta[3])/5 + 144*g2^2*g3^2*M3*Zeta[3] - 
      32*g3^4*M3*Zeta[3]))}
