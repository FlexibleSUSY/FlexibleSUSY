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

{(-6*g1^2*M1^2)/5 - 6*g2^2*M2^2 + 6*trace[ht, Adj[ht]] + 
  6*mh2*trace[Adj[Yt], Yt] + 6*trace[Yt, Adj[Yt], mt] + 
  6*trace[Adj[Yt], Yt, mq], (621*g1^4*M1^2)/25 + (18*g1^2*g2^2*M1^2)/5 + 
  (18*g1^2*g2^2*M1*M2)/5 + (18*g1^2*g2^2*M2^2)/5 + 33*g2^4*M2^2 + 
  (9*g1^4*mh1)/25 + 3*g2^4*mh1 + (9*g1^4*mh2)/25 + 3*g2^4*mh2 + 
  (6*g1^4*trace[mb])/25 + (18*g1^4*trace[me])/25 + (9*g1^4*trace[ml])/25 + 
  3*g2^4*trace[ml] + (3*g1^4*trace[mq])/25 + 9*g2^4*trace[mq] + 
  (24*g1^4*trace[mt])/25 + (8*g1^2*trace[ht, Adj[ht]])/5 + 
  32*g3^2*trace[ht, Adj[ht]] + M1*((-8*g1^2*trace[ht, Adj[Yt]])/5 - 
    (8*g1^2*trace[Adj[ht], Yt])/5) + 
  M3*(-32*g3^2*trace[ht, Adj[Yt]] - 32*g3^2*trace[Adj[ht], Yt]) + 
  (16*g1^2*M1^2*trace[Adj[Yt], Yt])/5 + 64*g3^2*M3^2*trace[Adj[Yt], Yt] + 
  (8*g1^2*mh2*trace[Adj[Yt], Yt])/5 + 32*g3^2*mh2*trace[Adj[Yt], Yt] + 
  (8*g1^2*trace[Yt, Adj[Yt], mt])/5 + 32*g3^2*trace[Yt, Adj[Yt], mt] + 
  (8*g1^2*trace[Adj[Yt], Yt, mq])/5 + 32*g3^2*trace[Adj[Yt], Yt, mq] - 
  6*trace[hb, Adj[ht], Yt, Adj[Yb]] - 6*trace[Adj[hb], hb, Adj[Yt], Yt] - 
  6*trace[Adj[ht], ht, Adj[Yb], Yb] - 36*trace[Adj[ht], ht, Adj[Yt], Yt] - 
  6*trace[Adj[Yt], ht, Adj[hb], Yb] - 36*trace[Adj[Yt], ht, Adj[ht], Yt] - 
  6*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  6*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  36*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 
  6*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] - 
  6*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] - 
  36*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] - 
  6*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] - 
  6*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] - 
  36*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq], (-621*g1^6*mh1)/125 - 
  (27*g1^4*g2^2*mh1)/25 - (9*g1^2*g2^4*mh1)/5 - 3*g2^6*mh1 - 
  (621*g1^6*mh2)/125 - (27*g1^4*g2^2*mh2)/25 - (9*g1^2*g2^4*mh2)/5 - 
  3*g2^6*mh2 - (414*g1^6*trace[mb])/125 - (18*g1^4*g2^2*trace[mb])/25 - 
  (1242*g1^6*trace[me])/125 - (54*g1^4*g2^2*trace[me])/25 - 
  (621*g1^6*trace[ml])/125 - (27*g1^4*g2^2*trace[ml])/25 - 
  (9*g1^2*g2^4*trace[ml])/5 - 3*g2^6*trace[ml] - (207*g1^6*trace[mq])/125 - 
  (9*g1^4*g2^2*trace[mq])/25 - (27*g1^2*g2^4*trace[mq])/5 - 
  9*g2^6*trace[mq] - (1656*g1^6*trace[mt])/125 - 
  (72*g1^4*g2^2*trace[mt])/25 - (126*g1^4*trace[hb, Adj[hb]])/25 - 
  54*g2^4*trace[hb, Adj[hb]] - (162*g1^4*trace[he, Adj[he]])/25 - 
  18*g2^4*trace[he, Adj[he]] - (10849*g1^4*trace[ht, Adj[ht]])/150 - 
  (57*g1^2*g2^2*trace[ht, Adj[ht]])/5 - (243*g2^4*trace[ht, Adj[ht]])/2 - 
  (248*g1^2*g3^2*trace[ht, Adj[ht]])/3 - 264*g2^2*g3^2*trace[ht, Adj[ht]] - 
  (320*g3^4*trace[ht, Adj[ht]])/3 - (126*g1^4*mh1*trace[Adj[Yb], Yb])/25 - 
  54*g2^4*mh1*trace[Adj[Yb], Yb] - (162*g1^4*mh1*trace[Adj[Ye], Ye])/25 - 
  18*g2^4*mh1*trace[Adj[Ye], Ye] - (24*g1^4*mh1*trace[Adj[Yt], Yt])/25 - 
  (10993*g1^4*mh2*trace[Adj[Yt], Yt])/150 - 
  (57*g1^2*g2^2*mh2*trace[Adj[Yt], Yt])/5 - (243*g2^4*mh2*trace[Adj[Yt], Yt])/
   2 - (248*g1^2*g3^2*mh2*trace[Adj[Yt], Yt])/3 - 
  264*g2^2*g3^2*mh2*trace[Adj[Yt], Yt] - (320*g3^4*mh2*trace[Adj[Yt], Yt])/
   3 - (16*g1^4*trace[mb]*trace[Adj[Yt], Yt])/25 - 
  32*g3^4*trace[mb]*trace[Adj[Yt], Yt] - 
  (48*g1^4*trace[me]*trace[Adj[Yt], Yt])/25 - 
  (24*g1^4*trace[ml]*trace[Adj[Yt], Yt])/25 - 
  (8*g1^4*trace[mq]*trace[Adj[Yt], Yt])/25 - 
  64*g3^4*trace[mq]*trace[Adj[Yt], Yt] - 
  (64*g1^4*trace[mt]*trace[Adj[Yt], Yt])/25 - 
  32*g3^4*trace[mt]*trace[Adj[Yt], Yt] - (126*g1^4*trace[Yb, Adj[Yb], mb])/
   25 - 54*g2^4*trace[Yb, Adj[Yb], mb] - (162*g1^4*trace[Ye, Adj[Ye], me])/
   25 - 18*g2^4*trace[Ye, Adj[Ye], me] - (10849*g1^4*trace[Yt, Adj[Yt], mt])/
   150 - (57*g1^2*g2^2*trace[Yt, Adj[Yt], mt])/5 - 
  (243*g2^4*trace[Yt, Adj[Yt], mt])/2 - 
  (248*g1^2*g3^2*trace[Yt, Adj[Yt], mt])/3 - 
  264*g2^2*g3^2*trace[Yt, Adj[Yt], mt] - (320*g3^4*trace[Yt, Adj[Yt], mt])/
   3 - (126*g1^4*trace[Adj[Yb], Yb, mq])/25 - 
  54*g2^4*trace[Adj[Yb], Yb, mq] - (162*g1^4*trace[Adj[Ye], Ye, ml])/25 - 
  18*g2^4*trace[Adj[Ye], Ye, ml] - (10849*g1^4*trace[Adj[Yt], Yt, mq])/150 - 
  (57*g1^2*g2^2*trace[Adj[Yt], Yt, mq])/5 - (243*g2^4*trace[Adj[Yt], Yt, mq])/
   2 - (248*g1^2*g3^2*trace[Adj[Yt], Yt, mq])/3 - 
  264*g2^2*g3^2*trace[Adj[Yt], Yt, mq] - (320*g3^4*trace[Adj[Yt], Yt, mq])/
   3 + (12*g1^2*trace[hb, Adj[ht], Yt, Adj[Yb]])/5 + 
  36*g2^2*trace[hb, Adj[ht], Yt, Adj[Yb]] + 
  48*g3^2*trace[hb, Adj[ht], Yt, Adj[Yb]] + 36*trace[Adj[Yb], Yb]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 12*trace[Adj[Ye], Ye]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 36*trace[Adj[hb], Yb]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 12*trace[Adj[he], Ye]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  (12*g1^2*trace[Adj[hb], hb, Adj[Yt], Yt])/5 + 
  36*g2^2*trace[Adj[hb], hb, Adj[Yt], Yt] + 
  48*g3^2*trace[Adj[hb], hb, Adj[Yt], Yt] + 36*trace[Adj[Yb], Yb]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 12*trace[Adj[Ye], Ye]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 
  (12*g1^2*trace[Adj[ht], ht, Adj[Yb], Yb])/5 + 
  36*g2^2*trace[Adj[ht], ht, Adj[Yb], Yb] + 
  48*g3^2*trace[Adj[ht], ht, Adj[Yb], Yb] + 36*trace[Adj[Yb], Yb]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 12*trace[Adj[Ye], Ye]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 
  (228*g1^2*trace[Adj[ht], ht, Adj[Yt], Yt])/5 + 
  36*g2^2*trace[Adj[ht], ht, Adj[Yt], Yt] + 
  288*g3^2*trace[Adj[ht], ht, Adj[Yt], Yt] + 216*trace[Adj[Yt], Yt]*
   trace[Adj[ht], ht, Adj[Yt], Yt] + 36*trace[hb, Adj[Yb]]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 12*trace[he, Adj[Ye]]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 
  (12*g1^2*trace[Adj[Yt], ht, Adj[hb], Yb])/5 + 
  36*g2^2*trace[Adj[Yt], ht, Adj[hb], Yb] + 
  48*g3^2*trace[Adj[Yt], ht, Adj[hb], Yb] + 36*trace[Adj[Yb], Yb]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 12*trace[Adj[Ye], Ye]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 
  (228*g1^2*trace[Adj[Yt], ht, Adj[ht], Yt])/5 + 
  36*g2^2*trace[Adj[Yt], ht, Adj[ht], Yt] + 
  288*g3^2*trace[Adj[Yt], ht, Adj[ht], Yt] + 216*trace[Adj[Yt], Yt]*
   trace[Adj[Yt], ht, Adj[ht], Yt] + 36*trace[Adj[hb], Yb]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 12*trace[Adj[he], Ye]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 216*trace[Adj[ht], Yt]*
   trace[Adj[Yt], ht, Adj[Yt], Yt] + 36*trace[hb, Adj[Yb]]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] + 12*trace[he, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] + 216*trace[ht, Adj[Yt]]*
   trace[Adj[Yt], Yt, Adj[ht], Yt] + 
  (12*g1^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  36*g2^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  48*g3^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (12*g1^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  36*g2^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  48*g3^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*trace[hb, Adj[hb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*trace[he, Adj[he]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  72*mh1*trace[Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*mh2*trace[Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*mh1*trace[Adj[Ye], Ye]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*mh2*trace[Adj[Ye], Ye]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*trace[Yb, Adj[Yb], mb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*trace[Ye, Adj[Ye], me]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*trace[Adj[Yb], Yb, mq]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*trace[Adj[Ye], Ye, ml]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  (228*g1^2*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
  36*g2^2*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  288*g3^2*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  108*trace[ht, Adj[ht]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  324*mh2*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  108*trace[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  108*trace[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  (12*g1^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb])/5 + 
  36*g2^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  48*g3^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  36*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  12*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  (12*g1^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt])/5 + 
  36*g2^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  48*g3^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  36*trace[Adj[Yb], Yb]*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  12*trace[Adj[Ye], Ye]*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  (228*g1^2*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt])/5 + 
  36*g2^2*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
  288*g3^2*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
  216*trace[Adj[Yt], Yt]*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
  (12*g1^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq])/5 + 
  36*g2^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  48*g3^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  36*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  12*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  (12*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq])/5 + 
  36*g2^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  48*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  36*trace[Adj[Yb], Yb]*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  12*trace[Adj[Ye], Ye]*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  (228*g1^2*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq])/5 + 
  36*g2^2*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  288*g3^2*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  216*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  18*trace[Adj[hb], hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*trace[Adj[ht], Yt, Adj[Yb], hb, Adj[Yb], Yb] + 
  18*trace[Adj[Yb], hb, Adj[ht], Yt, Adj[Yb], Yb] + 
  18*trace[Adj[Yb], hb, Adj[Yt], Yt, Adj[hb], Yb] + 
  18*trace[Adj[Yb], Yb, Adj[ht], ht, Adj[Yb], Yb] + 
  18*trace[Adj[Yb], Yb, Adj[Yt], ht, Adj[hb], Yb] + 
  36*mh1*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*mh2*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*trace[Adj[Yt], ht, Adj[Yb], Yb, Adj[hb], Yb] + 
  18*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[ht], Yt] + 
  18*trace[Adj[Yt], Yt, Adj[hb], hb, Adj[Yb], Yb] + 
  18*trace[Adj[Yt], Yt, Adj[ht], ht, Adj[Yt], Yt] + 
  18*trace[Adj[Yt], Yt, Adj[Yb], hb, Adj[hb], Yb] + 
  18*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[ht], Yt] + 
  18*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt] + 
  18*trace[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  18*trace[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb], mb] + 
  18*trace[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt], mt] + 
  18*trace[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
  18*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  18*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  18*trace[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  18*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt, mq] - 
  (13*g1^4*trace[ht, Adj[ht]]*Zeta[3])/5 + 
  (126*g1^2*g2^2*trace[ht, Adj[ht]]*Zeta[3])/5 - 
  189*g2^4*trace[ht, Adj[ht]]*Zeta[3] + 
  (416*g1^2*g3^2*trace[ht, Adj[ht]]*Zeta[3])/5 + 
  288*g2^2*g3^2*trace[ht, Adj[ht]]*Zeta[3] - 32*g3^4*trace[ht, Adj[ht]]*
   Zeta[3] - (13*g1^4*mh2*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
  (126*g1^2*g2^2*mh2*trace[Adj[Yt], Yt]*Zeta[3])/5 - 
  189*g2^4*mh2*trace[Adj[Yt], Yt]*Zeta[3] + 
  (416*g1^2*g3^2*mh2*trace[Adj[Yt], Yt]*Zeta[3])/5 + 
  288*g2^2*g3^2*mh2*trace[Adj[Yt], Yt]*Zeta[3] - 
  32*g3^4*mh2*trace[Adj[Yt], Yt]*Zeta[3] - 
  (13*g1^4*trace[Yt, Adj[Yt], mt]*Zeta[3])/5 + 
  (126*g1^2*g2^2*trace[Yt, Adj[Yt], mt]*Zeta[3])/5 - 
  189*g2^4*trace[Yt, Adj[Yt], mt]*Zeta[3] + 
  (416*g1^2*g3^2*trace[Yt, Adj[Yt], mt]*Zeta[3])/5 + 
  288*g2^2*g3^2*trace[Yt, Adj[Yt], mt]*Zeta[3] - 
  32*g3^4*trace[Yt, Adj[Yt], mt]*Zeta[3] - 
  (13*g1^4*trace[Adj[Yt], Yt, mq]*Zeta[3])/5 + 
  (126*g1^2*g2^2*trace[Adj[Yt], Yt, mq]*Zeta[3])/5 - 
  189*g2^4*trace[Adj[Yt], Yt, mq]*Zeta[3] + 
  (416*g1^2*g3^2*trace[Adj[Yt], Yt, mq]*Zeta[3])/5 + 
  288*g2^2*g3^2*trace[Adj[Yt], Yt, mq]*Zeta[3] - 
  32*g3^4*trace[Adj[Yt], Yt, mq]*Zeta[3] + 
  (12*g1^2*trace[hb, Adj[ht], Yt, Adj[Yb]]*Zeta[3])/5 - 
  96*g3^2*trace[hb, Adj[ht], Yt, Adj[Yb]]*Zeta[3] + 
  (12*g1^2*trace[Adj[hb], hb, Adj[Yt], Yt]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[hb], hb, Adj[Yt], Yt]*Zeta[3] + 
  (12*g1^2*trace[Adj[ht], ht, Adj[Yb], Yb]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[ht], ht, Adj[Yb], Yb]*Zeta[3] - 
  (72*g1^2*trace[Adj[ht], ht, Adj[Yt], Yt]*Zeta[3])/5 + 
  216*g2^2*trace[Adj[ht], ht, Adj[Yt], Yt]*Zeta[3] - 
  576*g3^2*trace[Adj[ht], ht, Adj[Yt], Yt]*Zeta[3] + 
  (12*g1^2*trace[Adj[Yt], ht, Adj[hb], Yb]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[Yt], ht, Adj[hb], Yb]*Zeta[3] - 
  (72*g1^2*trace[Adj[Yt], ht, Adj[ht], Yt]*Zeta[3])/5 + 
  216*g2^2*trace[Adj[Yt], ht, Adj[ht], Yt]*Zeta[3] - 
  576*g3^2*trace[Adj[Yt], ht, Adj[ht], Yt]*Zeta[3] + 
  (12*g1^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
  96*g3^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] + 
  (12*g1^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
  96*g3^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] - 
  (72*g1^2*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3])/5 + 
  216*g2^2*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] - 
  576*g3^2*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  (12*g1^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb]*Zeta[3])/5 - 
  96*g3^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb]*Zeta[3] + 
  (12*g1^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt]*Zeta[3])/5 - 
  96*g3^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt]*Zeta[3] - 
  (72*g1^2*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt]*Zeta[3])/5 + 
  216*g2^2*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt]*Zeta[3] - 
  576*g3^2*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt]*Zeta[3] + 
  (12*g1^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq]*Zeta[3] + 
  (12*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]*Zeta[3] - 
  (72*g1^2*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq]*Zeta[3])/5 + 
  216*g2^2*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq]*Zeta[3] - 
  576*g3^2*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq]*Zeta[3] + 
  108*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[ht], Yt]*Zeta[3] + 
  108*trace[Adj[Yt], Yt, Adj[ht], ht, Adj[Yt], Yt]*Zeta[3] + 
  108*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[ht], Yt]*Zeta[3] + 
  108*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3] + 
  108*trace[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], mt]*Zeta[3] + 
  108*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt, mq]*Zeta[3] + 
  M1^2*((55767*g1^6)/125 - (21492*g1^6*Zeta[3])/125) + 
  M1^2*((81*g1^4*g2^2)/25 - (1458*g1^4*g2^2*Zeta[3])/25) + 
  M1*M2*((54*g1^4*g2^2)/25 - (972*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((108*g1^4*g2^2)/25 - (486*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((171*g1^2*g2^4)/5 - (486*g1^2*g2^4*Zeta[3])/5) + 
  M1*M2*(18*g1^2*g2^4 - (324*g1^2*g2^4*Zeta[3])/5) + 
  M1^2*((72*g1^2*g2^4)/5 - (162*g1^2*g2^4*Zeta[3])/5) + 
  M2^2*(2157*g2^6 + 3780*g2^6*Zeta[3]) + 
  M1^2*((792*g1^4*g3^2)/5 - (4752*g1^4*g3^2*Zeta[3])/25) + 
  M1*M3*((528*g1^4*g3^2)/5 - (3168*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((1584*g1^4*g3^2)/25 - (1584*g1^4*g3^2*Zeta[3])/25) + 
  M2^2*(1080*g2^4*g3^2 - 1296*g2^4*g3^2*Zeta[3]) + 
  M2*M3*(720*g2^4*g3^2 - 864*g2^4*g3^2*Zeta[3]) + 
  M3^2*(432*g2^4*g3^2 - 432*g2^4*g3^2*Zeta[3]) + 
  M1*((42*g1^4*trace[hb, Adj[Yb]])/5 + (54*g1^4*trace[he, Adj[Ye]])/5 + 
    (2123*g1^4*trace[ht, Adj[Yt]])/15 + (42*g1^4*trace[Adj[hb], Yb])/5 + 
    (54*g1^4*trace[Adj[he], Ye])/5 + (2123*g1^4*trace[Adj[ht], Yt])/15 + 
    (26*g1^4*trace[ht, Adj[Yt]]*Zeta[3])/5 + 
    (26*g1^4*trace[Adj[ht], Yt]*Zeta[3])/5) + 
  M1*((57*g1^2*g2^2*trace[ht, Adj[Yt]])/5 + (57*g1^2*g2^2*trace[Adj[ht], Yt])/
     5 - (126*g1^2*g2^2*trace[ht, Adj[Yt]]*Zeta[3])/5 - 
    (126*g1^2*g2^2*trace[Adj[ht], Yt]*Zeta[3])/5) + 
  M2*((57*g1^2*g2^2*trace[ht, Adj[Yt]])/5 + (57*g1^2*g2^2*trace[Adj[ht], Yt])/
     5 - (126*g1^2*g2^2*trace[ht, Adj[Yt]]*Zeta[3])/5 - 
    (126*g1^2*g2^2*trace[Adj[ht], Yt]*Zeta[3])/5) + 
  M2*(90*g2^4*trace[hb, Adj[Yb]] + 30*g2^4*trace[he, Adj[Ye]] + 
    225*g2^4*trace[ht, Adj[Yt]] + 90*g2^4*trace[Adj[hb], Yb] + 
    30*g2^4*trace[Adj[he], Ye] + 225*g2^4*trace[Adj[ht], Yt] + 
    378*g2^4*trace[ht, Adj[Yt]]*Zeta[3] + 378*g2^4*trace[Adj[ht], Yt]*
     Zeta[3]) + M1*((248*g1^2*g3^2*trace[ht, Adj[Yt]])/3 + 
    (248*g1^2*g3^2*trace[Adj[ht], Yt])/3 - 
    (416*g1^2*g3^2*trace[ht, Adj[Yt]]*Zeta[3])/5 - 
    (416*g1^2*g3^2*trace[Adj[ht], Yt]*Zeta[3])/5) + 
  M3*((248*g1^2*g3^2*trace[ht, Adj[Yt]])/3 + 
    (248*g1^2*g3^2*trace[Adj[ht], Yt])/3 - 
    (416*g1^2*g3^2*trace[ht, Adj[Yt]]*Zeta[3])/5 - 
    (416*g1^2*g3^2*trace[Adj[ht], Yt]*Zeta[3])/5) + 
  M2*(264*g2^2*g3^2*trace[ht, Adj[Yt]] + 264*g2^2*g3^2*trace[Adj[ht], Yt] - 
    288*g2^2*g3^2*trace[ht, Adj[Yt]]*Zeta[3] - 
    288*g2^2*g3^2*trace[Adj[ht], Yt]*Zeta[3]) + 
  M3*(264*g2^2*g3^2*trace[ht, Adj[Yt]] + 264*g2^2*g3^2*trace[Adj[ht], Yt] - 
    288*g2^2*g3^2*trace[ht, Adj[Yt]]*Zeta[3] - 
    288*g2^2*g3^2*trace[Adj[ht], Yt]*Zeta[3]) + 
  M3*((640*g3^4*trace[ht, Adj[Yt]])/3 + (640*g3^4*trace[Adj[ht], Yt])/3 + 
    64*g3^4*trace[ht, Adj[Yt]]*Zeta[3] + 64*g3^4*trace[Adj[ht], Yt]*
     Zeta[3]) + M1^2*((-126*g1^4*trace[Adj[Yb], Yb])/5 - 
    (162*g1^4*trace[Adj[Ye], Ye])/5 - (2123*g1^4*trace[Adj[Yt], Yt])/5 - 
    (78*g1^4*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M1^2*((-114*g1^2*g2^2*trace[Adj[Yt], Yt])/5 + 
    (252*g1^2*g2^2*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M1*M2*((-114*g1^2*g2^2*trace[Adj[Yt], Yt])/5 + 
    (252*g1^2*g2^2*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M2^2*((-114*g1^2*g2^2*trace[Adj[Yt], Yt])/5 + 
    (252*g1^2*g2^2*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M2^2*(-270*g2^4*trace[Adj[Yb], Yb] - 90*g2^4*trace[Adj[Ye], Ye] - 
    675*g2^4*trace[Adj[Yt], Yt] - 1134*g2^4*trace[Adj[Yt], Yt]*Zeta[3]) + 
  M1^2*((-496*g1^2*g3^2*trace[Adj[Yt], Yt])/3 + 
    (832*g1^2*g3^2*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M1*M3*((-496*g1^2*g3^2*trace[Adj[Yt], Yt])/3 + 
    (832*g1^2*g3^2*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M3^2*((-496*g1^2*g3^2*trace[Adj[Yt], Yt])/3 + 
    (832*g1^2*g3^2*trace[Adj[Yt], Yt]*Zeta[3])/5) + 
  M2^2*(-528*g2^2*g3^2*trace[Adj[Yt], Yt] + 576*g2^2*g3^2*trace[Adj[Yt], Yt]*
     Zeta[3]) + M2*M3*(-528*g2^2*g3^2*trace[Adj[Yt], Yt] + 
    576*g2^2*g3^2*trace[Adj[Yt], Yt]*Zeta[3]) + 
  M3^2*(-528*g2^2*g3^2*trace[Adj[Yt], Yt] + 576*g2^2*g3^2*trace[Adj[Yt], Yt]*
     Zeta[3]) + M3^2*(-448*g3^4*trace[Adj[Yt], Yt] - 
    192*g3^4*trace[Adj[Yt], Yt]*Zeta[3]) + 
  M1*((-12*g1^2*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 - 
    (12*g1^2*trace[Adj[ht], Yt, Adj[Yb], Yb])/5 - 
    (12*g1^2*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 - 
    (228*g1^2*trace[Adj[Yt], ht, Adj[Yt], Yt])/5 - 
    (12*g1^2*trace[Adj[Yt], Yt, Adj[hb], Yb])/5 - 
    (228*g1^2*trace[Adj[Yt], Yt, Adj[ht], Yt])/5 - 
    (12*g1^2*trace[hb, Adj[Yt], Yt, Adj[Yb]]*Zeta[3])/5 - 
    (12*g1^2*trace[Adj[ht], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
    (12*g1^2*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3])/5 + 
    (72*g1^2*trace[Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3])/5 - 
    (12*g1^2*trace[Adj[Yt], Yt, Adj[hb], Yb]*Zeta[3])/5 + 
    (72*g1^2*trace[Adj[Yt], Yt, Adj[ht], Yt]*Zeta[3])/5) + 
  M2*(-36*g2^2*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
    36*g2^2*trace[Adj[ht], Yt, Adj[Yb], Yb] - 
    36*g2^2*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
    36*g2^2*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
    36*g2^2*trace[Adj[Yt], Yt, Adj[hb], Yb] - 
    36*g2^2*trace[Adj[Yt], Yt, Adj[ht], Yt] - 
    216*g2^2*trace[Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] - 
    216*g2^2*trace[Adj[Yt], Yt, Adj[ht], Yt]*Zeta[3]) + 
  M3*(-48*g3^2*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
    48*g3^2*trace[Adj[ht], Yt, Adj[Yb], Yb] - 
    48*g3^2*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
    288*g3^2*trace[Adj[Yt], ht, Adj[Yt], Yt] - 
    48*g3^2*trace[Adj[Yt], Yt, Adj[hb], Yb] - 
    288*g3^2*trace[Adj[Yt], Yt, Adj[ht], Yt] + 
    96*g3^2*trace[hb, Adj[Yt], Yt, Adj[Yb]]*Zeta[3] + 
    96*g3^2*trace[Adj[ht], Yt, Adj[Yb], Yb]*Zeta[3] + 
    96*g3^2*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3] + 
    576*g3^2*trace[Adj[Yt], ht, Adj[Yt], Yt]*Zeta[3] + 
    96*g3^2*trace[Adj[Yt], Yt, Adj[hb], Yb]*Zeta[3] + 
    576*g3^2*trace[Adj[Yt], Yt, Adj[ht], Yt]*Zeta[3]) + 
  M1^2*((24*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
    (228*g1^2*trace[Adj[Yt], Yt, Adj[Yt], Yt])/5 + 
    (24*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
    (72*g1^2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3])/5) + 
  M2^2*(72*g2^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    36*g2^2*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
    216*g2^2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3]) + 
  M3^2*(96*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    288*g3^2*trace[Adj[Yt], Yt, Adj[Yt], Yt] - 
    192*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] - 
    576*g3^2*trace[Adj[Yt], Yt, Adj[Yt], Yt]*Zeta[3])}
