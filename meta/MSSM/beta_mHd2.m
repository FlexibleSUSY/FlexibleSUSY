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

{(-6*g1^2*M1^2)/5 - 6*g2^2*M2^2 + 6*trace[hb, Adj[hb]] + 
  2*trace[he, Adj[he]] + 6*mh1*trace[Adj[Yb], Yb] + 
  2*mh1*trace[Adj[Ye], Ye] + 6*trace[Yb, Adj[Yb], mb] + 
  2*trace[Ye, Adj[Ye], me] + 6*trace[Adj[Yb], Yb, mq] + 
  2*trace[Adj[Ye], Ye, ml], (621*g1^4*M1^2)/25 + (18*g1^2*g2^2*M1^2)/5 + 
  (18*g1^2*g2^2*M1*M2)/5 + (18*g1^2*g2^2*M2^2)/5 + 33*g2^4*M2^2 + 
  (9*g1^4*mh1)/25 + 3*g2^4*mh1 + (9*g1^4*mh2)/25 + 3*g2^4*mh2 + 
  (6*g1^4*trace[mb])/25 + (18*g1^4*trace[me])/25 + (9*g1^4*trace[ml])/25 + 
  3*g2^4*trace[ml] + (3*g1^4*trace[mq])/25 + 9*g2^4*trace[mq] + 
  (24*g1^4*trace[mt])/25 - (4*g1^2*trace[hb, Adj[hb]])/5 + 
  32*g3^2*trace[hb, Adj[hb]] + (12*g1^2*trace[he, Adj[he]])/5 + 
  M3*(-32*g3^2*trace[hb, Adj[Yb]] - 32*g3^2*trace[Adj[hb], Yb]) + 
  M1*((4*g1^2*trace[hb, Adj[Yb]])/5 - (12*g1^2*trace[he, Adj[Ye]])/5 + 
    (4*g1^2*trace[Adj[hb], Yb])/5 - (12*g1^2*trace[Adj[he], Ye])/5) + 
  64*g3^2*M3^2*trace[Adj[Yb], Yb] - (4*g1^2*mh1*trace[Adj[Yb], Yb])/5 + 
  32*g3^2*mh1*trace[Adj[Yb], Yb] + (12*g1^2*mh1*trace[Adj[Ye], Ye])/5 + 
  M1^2*((-8*g1^2*trace[Adj[Yb], Yb])/5 + (24*g1^2*trace[Adj[Ye], Ye])/5) - 
  (4*g1^2*trace[Yb, Adj[Yb], mb])/5 + 32*g3^2*trace[Yb, Adj[Yb], mb] + 
  (12*g1^2*trace[Ye, Adj[Ye], me])/5 - (4*g1^2*trace[Adj[Yb], Yb, mq])/5 + 
  32*g3^2*trace[Adj[Yb], Yb, mq] + (12*g1^2*trace[Adj[Ye], Ye, ml])/5 - 
  6*trace[hb, Adj[ht], Yt, Adj[Yb]] - 36*trace[Adj[hb], hb, Adj[Yb], Yb] - 
  6*trace[Adj[hb], hb, Adj[Yt], Yt] - 12*trace[Adj[he], he, Adj[Ye], Ye] - 
  6*trace[Adj[ht], ht, Adj[Yb], Yb] - 36*trace[Adj[Yb], hb, Adj[hb], Yb] - 
  36*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] - 
  12*trace[Adj[Ye], he, Adj[he], Ye] - 
  12*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  6*trace[Adj[Yt], ht, Adj[hb], Yb] - 6*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  6*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  36*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] - 
  6*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] - 
  12*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] - 
  6*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] - 
  36*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] - 
  6*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] - 
  12*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] - 
  6*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq], (-621*g1^6*mh1)/125 - 
  (27*g1^4*g2^2*mh1)/25 - (9*g1^2*g2^4*mh1)/5 - 3*g2^6*mh1 - 
  (621*g1^6*mh2)/125 - (27*g1^4*g2^2*mh2)/25 - (9*g1^2*g2^4*mh2)/5 - 
  3*g2^6*mh2 - (414*g1^6*trace[mb])/125 - (18*g1^4*g2^2*trace[mb])/25 - 
  (1242*g1^6*trace[me])/125 - (54*g1^4*g2^2*trace[me])/25 - 
  (621*g1^6*trace[ml])/125 - (27*g1^4*g2^2*trace[ml])/25 - 
  (9*g1^2*g2^4*trace[ml])/5 - 3*g2^6*trace[ml] - (207*g1^6*trace[mq])/125 - 
  (9*g1^4*g2^2*trace[mq])/25 - (27*g1^2*g2^4*trace[mq])/5 - 
  9*g2^6*trace[mq] - (1656*g1^6*trace[mt])/125 - 
  (72*g1^4*g2^2*trace[mt])/25 - (4501*g1^4*trace[hb, Adj[hb]])/150 - 
  (3*g1^2*g2^2*trace[hb, Adj[hb]])/5 - (243*g2^4*trace[hb, Adj[hb]])/2 - 
  (568*g1^2*g3^2*trace[hb, Adj[hb]])/15 - 264*g2^2*g3^2*trace[hb, Adj[hb]] - 
  (320*g3^4*trace[hb, Adj[hb]])/3 - (3069*g1^4*trace[he, Adj[he]])/50 - 
  (81*g1^2*g2^2*trace[he, Adj[he]])/5 - (81*g2^4*trace[he, Adj[he]])/2 - 
  (234*g1^4*trace[ht, Adj[ht]])/25 - 54*g2^4*trace[ht, Adj[ht]] - 
  (4429*g1^4*mh1*trace[Adj[Yb], Yb])/150 - 
  (3*g1^2*g2^2*mh1*trace[Adj[Yb], Yb])/5 - (243*g2^4*mh1*trace[Adj[Yb], Yb])/
   2 - (568*g1^2*g3^2*mh1*trace[Adj[Yb], Yb])/15 - 
  264*g2^2*g3^2*mh1*trace[Adj[Yb], Yb] - (320*g3^4*mh1*trace[Adj[Yb], Yb])/
   3 + (12*g1^4*mh2*trace[Adj[Yb], Yb])/25 + 
  (8*g1^4*trace[mb]*trace[Adj[Yb], Yb])/25 - 
  32*g3^4*trace[mb]*trace[Adj[Yb], Yb] + 
  (24*g1^4*trace[me]*trace[Adj[Yb], Yb])/25 + 
  (12*g1^4*trace[ml]*trace[Adj[Yb], Yb])/25 + 
  (4*g1^4*trace[mq]*trace[Adj[Yb], Yb])/25 - 
  64*g3^4*trace[mq]*trace[Adj[Yb], Yb] + 
  (32*g1^4*trace[mt]*trace[Adj[Yb], Yb])/25 - 
  32*g3^4*trace[mt]*trace[Adj[Yb], Yb] - (3141*g1^4*mh1*trace[Adj[Ye], Ye])/
   50 - (81*g1^2*g2^2*mh1*trace[Adj[Ye], Ye])/5 - 
  (81*g2^4*mh1*trace[Adj[Ye], Ye])/2 - (36*g1^4*mh2*trace[Adj[Ye], Ye])/25 - 
  (24*g1^4*trace[mb]*trace[Adj[Ye], Ye])/25 - 
  (72*g1^4*trace[me]*trace[Adj[Ye], Ye])/25 - 
  (36*g1^4*trace[ml]*trace[Adj[Ye], Ye])/25 - 
  (12*g1^4*trace[mq]*trace[Adj[Ye], Ye])/25 - 
  (96*g1^4*trace[mt]*trace[Adj[Ye], Ye])/25 - 
  (234*g1^4*mh2*trace[Adj[Yt], Yt])/25 - 54*g2^4*mh2*trace[Adj[Yt], Yt] - 
  (4501*g1^4*trace[Yb, Adj[Yb], mb])/150 - 
  (3*g1^2*g2^2*trace[Yb, Adj[Yb], mb])/5 - (243*g2^4*trace[Yb, Adj[Yb], mb])/
   2 - (568*g1^2*g3^2*trace[Yb, Adj[Yb], mb])/15 - 
  264*g2^2*g3^2*trace[Yb, Adj[Yb], mb] - (320*g3^4*trace[Yb, Adj[Yb], mb])/
   3 - (3069*g1^4*trace[Ye, Adj[Ye], me])/50 - 
  (81*g1^2*g2^2*trace[Ye, Adj[Ye], me])/5 - (81*g2^4*trace[Ye, Adj[Ye], me])/
   2 - (234*g1^4*trace[Yt, Adj[Yt], mt])/25 - 
  54*g2^4*trace[Yt, Adj[Yt], mt] - (4501*g1^4*trace[Adj[Yb], Yb, mq])/150 - 
  (3*g1^2*g2^2*trace[Adj[Yb], Yb, mq])/5 - (243*g2^4*trace[Adj[Yb], Yb, mq])/
   2 - (568*g1^2*g3^2*trace[Adj[Yb], Yb, mq])/15 - 
  264*g2^2*g3^2*trace[Adj[Yb], Yb, mq] - (320*g3^4*trace[Adj[Yb], Yb, mq])/
   3 - (3069*g1^4*trace[Adj[Ye], Ye, ml])/50 - 
  (81*g1^2*g2^2*trace[Adj[Ye], Ye, ml])/5 - (81*g2^4*trace[Adj[Ye], Ye, ml])/
   2 - (234*g1^4*trace[Adj[Yt], Yt, mq])/25 - 
  54*g2^4*trace[Adj[Yt], Yt, mq] - (24*g1^2*trace[hb, Adj[ht], Yt, Adj[Yb]])/
   5 + 36*g2^2*trace[hb, Adj[ht], Yt, Adj[Yb]] + 
  48*g3^2*trace[hb, Adj[ht], Yt, Adj[Yb]] + 36*trace[Adj[Yt], Yt]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 36*trace[Adj[ht], Yt]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*g1^2*trace[Adj[hb], hb, Adj[Yb], Yb] + 
  36*g2^2*trace[Adj[hb], hb, Adj[Yb], Yb] + 
  288*g3^2*trace[Adj[hb], hb, Adj[Yb], Yb] + 216*trace[Adj[Yb], Yb]*
   trace[Adj[hb], hb, Adj[Yb], Yb] + 72*trace[Adj[Ye], Ye]*
   trace[Adj[hb], hb, Adj[Yb], Yb] - 
  (24*g1^2*trace[Adj[hb], hb, Adj[Yt], Yt])/5 + 
  36*g2^2*trace[Adj[hb], hb, Adj[Yt], Yt] + 
  48*g3^2*trace[Adj[hb], hb, Adj[Yt], Yt] + 36*trace[Adj[Yt], Yt]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 
  36*g1^2*trace[Adj[he], he, Adj[Ye], Ye] + 
  12*g2^2*trace[Adj[he], he, Adj[Ye], Ye] + 72*trace[Adj[Yb], Yb]*
   trace[Adj[he], he, Adj[Ye], Ye] + 24*trace[Adj[Ye], Ye]*
   trace[Adj[he], he, Adj[Ye], Ye] - 
  (24*g1^2*trace[Adj[ht], ht, Adj[Yb], Yb])/5 + 
  36*g2^2*trace[Adj[ht], ht, Adj[Yb], Yb] + 
  48*g3^2*trace[Adj[ht], ht, Adj[Yb], Yb] + 36*trace[Adj[Yt], Yt]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 36*trace[ht, Adj[Yt]]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 
  12*g1^2*trace[Adj[Yb], hb, Adj[hb], Yb] + 
  36*g2^2*trace[Adj[Yb], hb, Adj[hb], Yb] + 
  288*g3^2*trace[Adj[Yb], hb, Adj[hb], Yb] + 216*trace[Adj[Yb], Yb]*
   trace[Adj[Yb], hb, Adj[hb], Yb] + 72*trace[Adj[Ye], Ye]*
   trace[Adj[Yb], hb, Adj[hb], Yb] + 216*trace[Adj[hb], Yb]*
   trace[Adj[Yb], hb, Adj[Yb], Yb] + 72*trace[Adj[he], Ye]*
   trace[Adj[Yb], hb, Adj[Yb], Yb] + 216*trace[hb, Adj[Yb]]*
   trace[Adj[Yb], Yb, Adj[hb], Yb] + 72*trace[he, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[hb], Yb] + 
  12*g1^2*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*g2^2*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  288*g3^2*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  108*trace[hb, Adj[hb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*trace[he, Adj[he]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  324*mh1*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  108*mh1*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  108*trace[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*trace[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  108*trace[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*trace[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*g1^2*trace[Adj[Ye], he, Adj[he], Ye] + 
  12*g2^2*trace[Adj[Ye], he, Adj[he], Ye] + 72*trace[Adj[Yb], Yb]*
   trace[Adj[Ye], he, Adj[he], Ye] + 24*trace[Adj[Ye], Ye]*
   trace[Adj[Ye], he, Adj[he], Ye] + 72*trace[Adj[hb], Yb]*
   trace[Adj[Ye], he, Adj[Ye], Ye] + 24*trace[Adj[he], Ye]*
   trace[Adj[Ye], he, Adj[Ye], Ye] + 72*trace[hb, Adj[Yb]]*
   trace[Adj[Ye], Ye, Adj[he], Ye] + 24*trace[he, Adj[Ye]]*
   trace[Adj[Ye], Ye, Adj[he], Ye] + 
  36*g1^2*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*g2^2*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  36*trace[hb, Adj[hb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*trace[he, Adj[he]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  108*mh1*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  36*mh1*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  36*trace[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*trace[Ye, Adj[Ye], me]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  36*trace[Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*trace[Adj[Ye], Ye, ml]*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
  (24*g1^2*trace[Adj[Yt], ht, Adj[hb], Yb])/5 + 
  36*g2^2*trace[Adj[Yt], ht, Adj[hb], Yb] + 
  48*g3^2*trace[Adj[Yt], ht, Adj[hb], Yb] + 36*trace[Adj[Yt], Yt]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 36*trace[Adj[ht], Yt]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 36*trace[ht, Adj[Yt]]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] - 
  (24*g1^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  36*g2^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  48*g3^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
  (24*g1^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
  36*g2^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  48*g3^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*trace[ht, Adj[ht]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*mh1*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  72*mh2*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*trace[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*trace[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*g1^2*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  36*g2^2*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  288*g3^2*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  216*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  72*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] - 
  (24*g1^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb])/5 + 
  36*g2^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  48*g3^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  36*trace[Adj[Yt], Yt]*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  36*g1^2*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
  12*g2^2*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
  72*trace[Adj[Yb], Yb]*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
  24*trace[Adj[Ye], Ye]*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] - 
  (24*g1^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt])/5 + 
  36*g2^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  48*g3^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  36*trace[Adj[Yt], Yt]*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  12*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  36*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  288*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  216*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  72*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] - 
  (24*g1^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq])/5 + 
  36*g2^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  48*g3^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  36*trace[Adj[Yt], Yt]*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  36*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
  12*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
  72*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
  24*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] - 
  (24*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq])/5 + 
  36*g2^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  48*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  36*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  18*trace[hb, Adj[ht], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  18*trace[hb, Adj[Yt], Yt, Adj[ht], Yt, Adj[Yb]] + 
  18*trace[Adj[ht], ht, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*trace[Adj[ht], Yt, Adj[Yt], ht, Adj[Yb], Yb] + 
  18*trace[Adj[Yb], hb, Adj[Yb], Yb, Adj[hb], Yb] + 
  18*trace[Adj[Yb], Yb, Adj[hb], hb, Adj[Yb], Yb] + 
  18*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[hb], Yb] + 
  18*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb] + 
  6*trace[Adj[Ye], he, Adj[Ye], Ye, Adj[he], Ye] + 
  6*trace[Adj[Ye], Ye, Adj[he], he, Adj[Ye], Ye] + 
  6*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[he], Ye] + 
  6*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*trace[Adj[Yt], ht, Adj[ht], Yt, Adj[Yb], Yb] + 
  18*trace[Adj[Yt], ht, Adj[Yt], Yt, Adj[hb], Yb] + 
  18*trace[Adj[Yt], Yt, Adj[hb], hb, Adj[Yt], Yt] + 
  18*trace[Adj[Yt], Yt, Adj[ht], ht, Adj[Yb], Yb] + 
  18*trace[Adj[Yt], Yt, Adj[Yt], ht, Adj[hb], Yb] + 
  18*mh1*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] + 
  36*mh2*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb] + 
  18*trace[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  18*trace[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], mb] + 
  6*trace[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], me] + 
  18*trace[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt], mt] + 
  18*trace[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  18*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  18*trace[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  6*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
  18*trace[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  18*trace[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], Yb, mq] - 
  (77*g1^4*trace[hb, Adj[hb]]*Zeta[3])/25 - 18*g1^2*g2^2*trace[hb, Adj[hb]]*
   Zeta[3] - 189*g2^4*trace[hb, Adj[hb]]*Zeta[3] + 
  (224*g1^2*g3^2*trace[hb, Adj[hb]]*Zeta[3])/5 + 
  288*g2^2*g3^2*trace[hb, Adj[hb]]*Zeta[3] - 32*g3^4*trace[hb, Adj[hb]]*
   Zeta[3] + (81*g1^4*trace[he, Adj[he]]*Zeta[3])/25 + 
  (162*g1^2*g2^2*trace[he, Adj[he]]*Zeta[3])/5 - 
  63*g2^4*trace[he, Adj[he]]*Zeta[3] - 
  (77*g1^4*mh1*trace[Adj[Yb], Yb]*Zeta[3])/25 - 
  18*g1^2*g2^2*mh1*trace[Adj[Yb], Yb]*Zeta[3] - 
  189*g2^4*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 
  (224*g1^2*g3^2*mh1*trace[Adj[Yb], Yb]*Zeta[3])/5 + 
  288*g2^2*g3^2*mh1*trace[Adj[Yb], Yb]*Zeta[3] - 
  32*g3^4*mh1*trace[Adj[Yb], Yb]*Zeta[3] + 
  (81*g1^4*mh1*trace[Adj[Ye], Ye]*Zeta[3])/25 + 
  (162*g1^2*g2^2*mh1*trace[Adj[Ye], Ye]*Zeta[3])/5 - 
  63*g2^4*mh1*trace[Adj[Ye], Ye]*Zeta[3] - 
  (77*g1^4*trace[Yb, Adj[Yb], mb]*Zeta[3])/25 - 
  18*g1^2*g2^2*trace[Yb, Adj[Yb], mb]*Zeta[3] - 
  189*g2^4*trace[Yb, Adj[Yb], mb]*Zeta[3] + 
  (224*g1^2*g3^2*trace[Yb, Adj[Yb], mb]*Zeta[3])/5 + 
  288*g2^2*g3^2*trace[Yb, Adj[Yb], mb]*Zeta[3] - 
  32*g3^4*trace[Yb, Adj[Yb], mb]*Zeta[3] + 
  (81*g1^4*trace[Ye, Adj[Ye], me]*Zeta[3])/25 + 
  (162*g1^2*g2^2*trace[Ye, Adj[Ye], me]*Zeta[3])/5 - 
  63*g2^4*trace[Ye, Adj[Ye], me]*Zeta[3] - 
  (77*g1^4*trace[Adj[Yb], Yb, mq]*Zeta[3])/25 - 
  18*g1^2*g2^2*trace[Adj[Yb], Yb, mq]*Zeta[3] - 
  189*g2^4*trace[Adj[Yb], Yb, mq]*Zeta[3] + 
  (224*g1^2*g3^2*trace[Adj[Yb], Yb, mq]*Zeta[3])/5 + 
  288*g2^2*g3^2*trace[Adj[Yb], Yb, mq]*Zeta[3] - 
  32*g3^4*trace[Adj[Yb], Yb, mq]*Zeta[3] + 
  (81*g1^4*trace[Adj[Ye], Ye, ml]*Zeta[3])/25 + 
  (162*g1^2*g2^2*trace[Adj[Ye], Ye, ml]*Zeta[3])/5 - 
  63*g2^4*trace[Adj[Ye], Ye, ml]*Zeta[3] + 
  (84*g1^2*trace[hb, Adj[ht], Yt, Adj[Yb]]*Zeta[3])/5 - 
  96*g3^2*trace[hb, Adj[ht], Yt, Adj[Yb]]*Zeta[3] + 
  (216*g1^2*trace[Adj[hb], hb, Adj[Yb], Yb]*Zeta[3])/5 + 
  216*g2^2*trace[Adj[hb], hb, Adj[Yb], Yb]*Zeta[3] - 
  576*g3^2*trace[Adj[hb], hb, Adj[Yb], Yb]*Zeta[3] + 
  (84*g1^2*trace[Adj[hb], hb, Adj[Yt], Yt]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[hb], hb, Adj[Yt], Yt]*Zeta[3] - 
  (216*g1^2*trace[Adj[he], he, Adj[Ye], Ye]*Zeta[3])/5 + 
  72*g2^2*trace[Adj[he], he, Adj[Ye], Ye]*Zeta[3] + 
  (84*g1^2*trace[Adj[ht], ht, Adj[Yb], Yb]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[ht], ht, Adj[Yb], Yb]*Zeta[3] + 
  (216*g1^2*trace[Adj[Yb], hb, Adj[hb], Yb]*Zeta[3])/5 + 
  216*g2^2*trace[Adj[Yb], hb, Adj[hb], Yb]*Zeta[3] - 
  576*g3^2*trace[Adj[Yb], hb, Adj[hb], Yb]*Zeta[3] + 
  (216*g1^2*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 + 
  216*g2^2*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
  576*g3^2*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
  (216*g1^2*trace[Adj[Ye], he, Adj[he], Ye]*Zeta[3])/5 + 
  72*g2^2*trace[Adj[Ye], he, Adj[he], Ye]*Zeta[3] - 
  (216*g1^2*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 + 
  72*g2^2*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  (84*g1^2*trace[Adj[Yt], ht, Adj[hb], Yb]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[Yt], ht, Adj[hb], Yb]*Zeta[3] + 
  (84*g1^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
  96*g3^2*mh1*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] + 
  (84*g1^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
  96*g3^2*mh2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3] + 
  (216*g1^2*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb]*Zeta[3])/5 + 
  216*g2^2*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb]*Zeta[3] - 
  576*g3^2*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb]*Zeta[3] + 
  (84*g1^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb]*Zeta[3])/5 - 
  96*g3^2*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb]*Zeta[3] - 
  (216*g1^2*trace[Ye, Adj[Ye], Ye, Adj[Ye], me]*Zeta[3])/5 + 
  72*g2^2*trace[Ye, Adj[Ye], Ye, Adj[Ye], me]*Zeta[3] + 
  (84*g1^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt]*Zeta[3])/5 - 
  96*g3^2*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt]*Zeta[3] + 
  (216*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq]*Zeta[3])/5 + 
  216*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq]*Zeta[3] - 
  576*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq]*Zeta[3] + 
  (84*g1^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq]*Zeta[3] - 
  (216*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml]*Zeta[3])/5 + 
  72*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml]*Zeta[3] + 
  (84*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]*Zeta[3])/5 - 
  96*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq]*Zeta[3] + 
  108*trace[Adj[Yb], hb, Adj[Yb], Yb, Adj[hb], Yb]*Zeta[3] + 
  108*trace[Adj[Yb], Yb, Adj[hb], hb, Adj[Yb], Yb]*Zeta[3] + 
  108*trace[Adj[Yb], Yb, Adj[Yb], hb, Adj[hb], Yb]*Zeta[3] + 
  108*mh1*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
  36*trace[Adj[Ye], he, Adj[Ye], Ye, Adj[he], Ye]*Zeta[3] + 
  36*trace[Adj[Ye], Ye, Adj[he], he, Adj[Ye], Ye]*Zeta[3] + 
  36*trace[Adj[Ye], Ye, Adj[Ye], he, Adj[he], Ye]*Zeta[3] + 
  36*mh1*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3] + 
  108*trace[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], mb]*Zeta[3] + 
  36*trace[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], me]*Zeta[3] + 
  108*trace[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb, mq]*Zeta[3] + 
  36*trace[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye, ml]*Zeta[3] + 
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
  M1*((568*g1^2*g3^2*trace[hb, Adj[Yb]])/15 + 
    (568*g1^2*g3^2*trace[Adj[hb], Yb])/15 - 
    (224*g1^2*g3^2*trace[hb, Adj[Yb]]*Zeta[3])/5 - 
    (224*g1^2*g3^2*trace[Adj[hb], Yb]*Zeta[3])/5) + 
  M3*((568*g1^2*g3^2*trace[hb, Adj[Yb]])/15 + 
    (568*g1^2*g3^2*trace[Adj[hb], Yb])/15 - 
    (224*g1^2*g3^2*trace[hb, Adj[Yb]]*Zeta[3])/5 - 
    (224*g1^2*g3^2*trace[Adj[hb], Yb]*Zeta[3])/5) + 
  M2*(264*g2^2*g3^2*trace[hb, Adj[Yb]] + 264*g2^2*g3^2*trace[Adj[hb], Yb] - 
    288*g2^2*g3^2*trace[hb, Adj[Yb]]*Zeta[3] - 
    288*g2^2*g3^2*trace[Adj[hb], Yb]*Zeta[3]) + 
  M3*(264*g2^2*g3^2*trace[hb, Adj[Yb]] + 264*g2^2*g3^2*trace[Adj[hb], Yb] - 
    288*g2^2*g3^2*trace[hb, Adj[Yb]]*Zeta[3] - 
    288*g2^2*g3^2*trace[Adj[hb], Yb]*Zeta[3]) + 
  M3*((640*g3^4*trace[hb, Adj[Yb]])/3 + (640*g3^4*trace[Adj[hb], Yb])/3 + 
    64*g3^4*trace[hb, Adj[Yb]]*Zeta[3] + 64*g3^4*trace[Adj[hb], Yb]*
     Zeta[3]) + M1*((175*g1^4*trace[hb, Adj[Yb]])/3 + 
    (603*g1^4*trace[he, Adj[Ye]])/5 + (78*g1^4*trace[ht, Adj[Yt]])/5 + 
    (175*g1^4*trace[Adj[hb], Yb])/3 + (603*g1^4*trace[Adj[he], Ye])/5 + 
    (78*g1^4*trace[Adj[ht], Yt])/5 + (154*g1^4*trace[hb, Adj[Yb]]*Zeta[3])/
     25 - (162*g1^4*trace[he, Adj[Ye]]*Zeta[3])/25 + 
    (154*g1^4*trace[Adj[hb], Yb]*Zeta[3])/25 - 
    (162*g1^4*trace[Adj[he], Ye]*Zeta[3])/25) + 
  M1*((3*g1^2*g2^2*trace[hb, Adj[Yb]])/5 + (81*g1^2*g2^2*trace[he, Adj[Ye]])/
     5 + (3*g1^2*g2^2*trace[Adj[hb], Yb])/5 + 
    (81*g1^2*g2^2*trace[Adj[he], Ye])/5 + 18*g1^2*g2^2*trace[hb, Adj[Yb]]*
     Zeta[3] - (162*g1^2*g2^2*trace[he, Adj[Ye]]*Zeta[3])/5 + 
    18*g1^2*g2^2*trace[Adj[hb], Yb]*Zeta[3] - 
    (162*g1^2*g2^2*trace[Adj[he], Ye]*Zeta[3])/5) + 
  M2*((3*g1^2*g2^2*trace[hb, Adj[Yb]])/5 + (81*g1^2*g2^2*trace[he, Adj[Ye]])/
     5 + (3*g1^2*g2^2*trace[Adj[hb], Yb])/5 + 
    (81*g1^2*g2^2*trace[Adj[he], Ye])/5 + 18*g1^2*g2^2*trace[hb, Adj[Yb]]*
     Zeta[3] - (162*g1^2*g2^2*trace[he, Adj[Ye]]*Zeta[3])/5 + 
    18*g1^2*g2^2*trace[Adj[hb], Yb]*Zeta[3] - 
    (162*g1^2*g2^2*trace[Adj[he], Ye]*Zeta[3])/5) + 
  M2*(225*g2^4*trace[hb, Adj[Yb]] + 75*g2^4*trace[he, Adj[Ye]] + 
    90*g2^4*trace[ht, Adj[Yt]] + 225*g2^4*trace[Adj[hb], Yb] + 
    75*g2^4*trace[Adj[he], Ye] + 90*g2^4*trace[Adj[ht], Yt] + 
    378*g2^4*trace[hb, Adj[Yb]]*Zeta[3] + 126*g2^4*trace[he, Adj[Ye]]*
     Zeta[3] + 378*g2^4*trace[Adj[hb], Yb]*Zeta[3] + 
    126*g2^4*trace[Adj[he], Ye]*Zeta[3]) + 
  M1^2*((-1136*g1^2*g3^2*trace[Adj[Yb], Yb])/15 + 
    (448*g1^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3])/5) + 
  M1*M3*((-1136*g1^2*g3^2*trace[Adj[Yb], Yb])/15 + 
    (448*g1^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3])/5) + 
  M3^2*((-1136*g1^2*g3^2*trace[Adj[Yb], Yb])/15 + 
    (448*g1^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3])/5) + 
  M2^2*(-528*g2^2*g3^2*trace[Adj[Yb], Yb] + 576*g2^2*g3^2*trace[Adj[Yb], Yb]*
     Zeta[3]) + M2*M3*(-528*g2^2*g3^2*trace[Adj[Yb], Yb] + 
    576*g2^2*g3^2*trace[Adj[Yb], Yb]*Zeta[3]) + 
  M3^2*(-528*g2^2*g3^2*trace[Adj[Yb], Yb] + 576*g2^2*g3^2*trace[Adj[Yb], Yb]*
     Zeta[3]) + M3^2*(-448*g3^4*trace[Adj[Yb], Yb] - 
    192*g3^4*trace[Adj[Yb], Yb]*Zeta[3]) + 
  M1^2*(-175*g1^4*trace[Adj[Yb], Yb] - (1809*g1^4*trace[Adj[Ye], Ye])/5 - 
    (234*g1^4*trace[Adj[Yt], Yt])/5 - (462*g1^4*trace[Adj[Yb], Yb]*Zeta[3])/
     25 + (486*g1^4*trace[Adj[Ye], Ye]*Zeta[3])/25) + 
  M1^2*((-6*g1^2*g2^2*trace[Adj[Yb], Yb])/5 - 
    (162*g1^2*g2^2*trace[Adj[Ye], Ye])/5 - 36*g1^2*g2^2*trace[Adj[Yb], Yb]*
     Zeta[3] + (324*g1^2*g2^2*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  M1*M2*((-6*g1^2*g2^2*trace[Adj[Yb], Yb])/5 - 
    (162*g1^2*g2^2*trace[Adj[Ye], Ye])/5 - 36*g1^2*g2^2*trace[Adj[Yb], Yb]*
     Zeta[3] + (324*g1^2*g2^2*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  M2^2*((-6*g1^2*g2^2*trace[Adj[Yb], Yb])/5 - 
    (162*g1^2*g2^2*trace[Adj[Ye], Ye])/5 - 36*g1^2*g2^2*trace[Adj[Yb], Yb]*
     Zeta[3] + (324*g1^2*g2^2*trace[Adj[Ye], Ye]*Zeta[3])/5) + 
  M2^2*(-675*g2^4*trace[Adj[Yb], Yb] - 225*g2^4*trace[Adj[Ye], Ye] - 
    270*g2^4*trace[Adj[Yt], Yt] - 1134*g2^4*trace[Adj[Yb], Yb]*Zeta[3] - 
    378*g2^4*trace[Adj[Ye], Ye]*Zeta[3]) + 
  M2*(-36*g2^2*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
    36*g2^2*trace[Adj[ht], Yt, Adj[Yb], Yb] - 
    36*g2^2*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
    36*g2^2*trace[Adj[Yb], Yb, Adj[hb], Yb] - 
    12*g2^2*trace[Adj[Ye], he, Adj[Ye], Ye] - 
    12*g2^2*trace[Adj[Ye], Ye, Adj[he], Ye] - 
    36*g2^2*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
    36*g2^2*trace[Adj[Yt], Yt, Adj[hb], Yb] - 
    216*g2^2*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] - 
    216*g2^2*trace[Adj[Yb], Yb, Adj[hb], Yb]*Zeta[3] - 
    72*g2^2*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3] - 
    72*g2^2*trace[Adj[Ye], Ye, Adj[he], Ye]*Zeta[3]) + 
  M2^2*(36*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    12*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
    72*g2^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
    216*g2^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] + 
    72*g2^2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3]) + 
  M1*((24*g1^2*trace[hb, Adj[Yt], Yt, Adj[Yb]])/5 + 
    (24*g1^2*trace[Adj[ht], Yt, Adj[Yb], Yb])/5 - 
    12*g1^2*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
    12*g1^2*trace[Adj[Yb], Yb, Adj[hb], Yb] - 
    36*g1^2*trace[Adj[Ye], he, Adj[Ye], Ye] - 
    36*g1^2*trace[Adj[Ye], Ye, Adj[he], Ye] + 
    (24*g1^2*trace[Adj[Yt], ht, Adj[Yb], Yb])/5 + 
    (24*g1^2*trace[Adj[Yt], Yt, Adj[hb], Yb])/5 - 
    (84*g1^2*trace[hb, Adj[Yt], Yt, Adj[Yb]]*Zeta[3])/5 - 
    (84*g1^2*trace[Adj[ht], Yt, Adj[Yb], Yb]*Zeta[3])/5 - 
    (216*g1^2*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3])/5 - 
    (216*g1^2*trace[Adj[Yb], Yb, Adj[hb], Yb]*Zeta[3])/5 + 
    (216*g1^2*trace[Adj[Ye], he, Adj[Ye], Ye]*Zeta[3])/5 + 
    (216*g1^2*trace[Adj[Ye], Ye, Adj[he], Ye]*Zeta[3])/5 - 
    (84*g1^2*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3])/5 - 
    (84*g1^2*trace[Adj[Yt], Yt, Adj[hb], Yb]*Zeta[3])/5) + 
  M3*(-48*g3^2*trace[hb, Adj[Yt], Yt, Adj[Yb]] - 
    48*g3^2*trace[Adj[ht], Yt, Adj[Yb], Yb] - 
    288*g3^2*trace[Adj[Yb], hb, Adj[Yb], Yb] - 
    288*g3^2*trace[Adj[Yb], Yb, Adj[hb], Yb] - 
    48*g3^2*trace[Adj[Yt], ht, Adj[Yb], Yb] - 
    48*g3^2*trace[Adj[Yt], Yt, Adj[hb], Yb] + 
    96*g3^2*trace[hb, Adj[Yt], Yt, Adj[Yb]]*Zeta[3] + 
    96*g3^2*trace[Adj[ht], Yt, Adj[Yb], Yb]*Zeta[3] + 
    576*g3^2*trace[Adj[Yb], hb, Adj[Yb], Yb]*Zeta[3] + 
    576*g3^2*trace[Adj[Yb], Yb, Adj[hb], Yb]*Zeta[3] + 
    96*g3^2*trace[Adj[Yt], ht, Adj[Yb], Yb]*Zeta[3] + 
    96*g3^2*trace[Adj[Yt], Yt, Adj[hb], Yb]*Zeta[3]) + 
  M1^2*(12*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    36*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye] - 
    (48*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb])/5 + 
    (216*g1^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3])/5 - 
    (216*g1^2*trace[Adj[Ye], Ye, Adj[Ye], Ye]*Zeta[3])/5 + 
    (168*g1^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])/5) + 
  M3^2*(288*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
    96*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb] - 
    576*g3^2*trace[Adj[Yb], Yb, Adj[Yb], Yb]*Zeta[3] - 
    192*g3^2*trace[Adj[Yt], Yt, Adj[Yb], Yb]*Zeta[3])}
