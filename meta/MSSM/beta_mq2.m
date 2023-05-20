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

{(-2*g1^2*M1^2)/15 - 6*g2^2*M2^2 - (32*g3^2*M3^2)/3 + 2*MatMul[Adj[hb], hb] + 
  2*MatMul[Adj[ht], ht] + 2*mh1*MatMul[Adj[Yb], Yb] + 
  2*mh2*MatMul[Adj[Yt], Yt] + MatMul[mq, Adj[Yb], Yb] + 
  MatMul[mq, Adj[Yt], Yt] + 2*MatMul[Adj[Yb], mb, Yb] + 
  MatMul[Adj[Yb], Yb, mq] + 2*MatMul[Adj[Yt], mt, Yt] + 
  MatMul[Adj[Yt], Yt, mq], (199*g1^4*M1^2)/75 + (2*g1^2*g2^2*M1^2)/5 + 
  (32*g1^2*g3^2*M1^2)/45 + (2*g1^2*g2^2*M1*M2)/5 + (2*g1^2*g2^2*M2^2)/5 + 
  33*g2^4*M2^2 + 32*g2^2*g3^2*M2^2 + (32*g1^2*g3^2*M1*M3)/45 + 
  32*g2^2*g3^2*M2*M3 + (32*g1^2*g3^2*M3^2)/45 + 32*g2^2*g3^2*M3^2 - 
  (128*g3^4*M3^2)/3 + (g1^4*mh1)/25 + 3*g2^4*mh1 + (g1^4*mh2)/25 + 
  3*g2^4*mh2 + (4*g1^2*MatMul[Adj[hb], hb])/5 - 
  (4*g1^2*M1*MatMul[Adj[hb], Yb])/5 + (8*g1^2*MatMul[Adj[ht], ht])/5 - 
  (8*g1^2*M1*MatMul[Adj[ht], Yt])/5 - (4*g1^2*M1*MatMul[Adj[Yb], hb])/5 + 
  (8*g1^2*M1^2*MatMul[Adj[Yb], Yb])/5 + (4*g1^2*mh1*MatMul[Adj[Yb], Yb])/5 - 
  (8*g1^2*M1*MatMul[Adj[Yt], ht])/5 + (16*g1^2*M1^2*MatMul[Adj[Yt], Yt])/5 + 
  (8*g1^2*mh2*MatMul[Adj[Yt], Yt])/5 + (2*g1^2*MatMul[mq, Adj[Yb], Yb])/5 + 
  (4*g1^2*MatMul[mq, Adj[Yt], Yt])/5 + (4*g1^2*MatMul[Adj[Yb], mb, Yb])/5 + 
  (2*g1^2*MatMul[Adj[Yb], Yb, mq])/5 + (8*g1^2*MatMul[Adj[Yt], mt, Yt])/5 + 
  (4*g1^2*MatMul[Adj[Yt], Yt, mq])/5 - 4*MatMul[Adj[hb], hb, Adj[Yb], Yb] - 
  4*MatMul[Adj[hb], Yb, Adj[Yb], hb] - 4*MatMul[Adj[ht], ht, Adj[Yt], Yt] - 
  4*MatMul[Adj[ht], Yt, Adj[Yt], ht] - 4*MatMul[Adj[Yb], hb, Adj[hb], Yb] - 
  4*MatMul[Adj[Yb], Yb, Adj[hb], hb] - 
  8*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] - 
  4*MatMul[Adj[Yt], ht, Adj[ht], Yt] - 4*MatMul[Adj[Yt], Yt, Adj[ht], ht] - 
  8*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] - 
  2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb] - 
  2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt] - 
  4*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb] - 
  4*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb] - 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb] - 
  2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq] - 
  4*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt] - 
  4*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt] - 
  4*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt] - 
  2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq] + (2*g1^4*trace[mb])/75 + 
  (16*g3^4*trace[mb])/3 + (2*g1^4*trace[me])/25 + (g1^4*trace[ml])/25 + 
  3*g2^4*trace[ml] + (g1^4*trace[mq])/75 + 9*g2^4*trace[mq] + 
  (32*g3^4*trace[mq])/3 + (8*g1^4*trace[mt])/75 + (16*g3^4*trace[mt])/3 - 
  6*MatMul[Adj[Yb], Yb]*trace[hb, Adj[hb]] - 6*MatMul[Adj[hb], Yb]*
   trace[hb, Adj[Yb]] - 2*MatMul[Adj[Yb], Yb]*trace[he, Adj[he]] - 
  2*MatMul[Adj[hb], Yb]*trace[he, Adj[Ye]] - 6*MatMul[Adj[Yt], Yt]*
   trace[ht, Adj[ht]] - 6*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]] - 
  6*MatMul[Adj[Yb], hb]*trace[Adj[hb], Yb] - 2*MatMul[Adj[Yb], hb]*
   trace[Adj[he], Ye] - 6*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt] - 
  6*MatMul[Adj[hb], hb]*trace[Adj[Yb], Yb] - 12*mh1*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb] - 3*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  6*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb] - 
  3*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb] - 
  2*MatMul[Adj[hb], hb]*trace[Adj[Ye], Ye] - 4*mh1*MatMul[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye] - MatMul[mq, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Ye], Ye] - 
  MatMul[Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye] - 
  6*MatMul[Adj[ht], ht]*trace[Adj[Yt], Yt] - 12*mh2*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt] - 3*MatMul[mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  6*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] - 
  3*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] - 
  6*MatMul[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  2*MatMul[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
  6*MatMul[Adj[Yt], Yt]*trace[Yt, Adj[Yt], mt] - 
  6*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
  2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
  6*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt, mq], 
 (-32*g1^2*g2^2*g3^2*M1^2)/5 - (32*g1^2*g2^2*g3^2*M1*M2)/5 - 
  (32*g1^2*g2^2*g3^2*M2^2)/5 - (32*g1^2*g2^2*g3^2*M1*M3)/5 - 
  (32*g1^2*g2^2*g3^2*M2*M3)/5 - (32*g1^2*g2^2*g3^2*M3^2)/5 - 
  (199*g1^6*mh1)/375 - (3*g1^4*g2^2*mh1)/25 - (g1^2*g2^4*mh1)/5 - 
  3*g2^6*mh1 - (16*g1^4*g3^2*mh1)/75 - 16*g2^4*g3^2*mh1 - 
  (199*g1^6*mh2)/375 - (3*g1^4*g2^2*mh2)/25 - (g1^2*g2^4*mh2)/5 - 
  3*g2^6*mh2 - (16*g1^4*g3^2*mh2)/75 - 16*g2^4*g3^2*mh2 - 
  (633*g1^4*MatMul[Adj[hb], hb])/50 - (41*g1^2*g2^2*MatMul[Adj[hb], hb])/5 - 
  (45*g2^4*MatMul[Adj[hb], hb])/2 - (152*g1^2*g3^2*MatMul[Adj[hb], hb])/15 - 
  8*g2^2*g3^2*MatMul[Adj[hb], hb] + (16*g3^4*MatMul[Adj[hb], hb])/3 + 
  8*g2^2*g3^2*M2*MatMul[Adj[hb], Yb] + 8*g2^2*g3^2*M3*MatMul[Adj[hb], Yb] - 
  (3767*g1^4*MatMul[Adj[ht], ht])/150 - (59*g1^2*g2^2*MatMul[Adj[ht], ht])/
   5 - (45*g2^4*MatMul[Adj[ht], ht])/2 - (136*g1^2*g3^2*MatMul[Adj[ht], ht])/
   5 - 8*g2^2*g3^2*MatMul[Adj[ht], ht] + (16*g3^4*MatMul[Adj[ht], ht])/3 + 
  8*g2^2*g3^2*M2*MatMul[Adj[ht], Yt] + 8*g2^2*g3^2*M3*MatMul[Adj[ht], Yt] + 
  8*g2^2*g3^2*M2*MatMul[Adj[Yb], hb] + 8*g2^2*g3^2*M3*MatMul[Adj[Yb], hb] - 
  16*g2^2*g3^2*M2^2*MatMul[Adj[Yb], Yb] - 16*g2^2*g3^2*M2*M3*
   MatMul[Adj[Yb], Yb] - 16*g2^2*g3^2*M3^2*MatMul[Adj[Yb], Yb] - 
  (657*g1^4*mh1*MatMul[Adj[Yb], Yb])/50 - 
  (41*g1^2*g2^2*mh1*MatMul[Adj[Yb], Yb])/5 - 
  (45*g2^4*mh1*MatMul[Adj[Yb], Yb])/2 - 
  (152*g1^2*g3^2*mh1*MatMul[Adj[Yb], Yb])/15 - 
  8*g2^2*g3^2*mh1*MatMul[Adj[Yb], Yb] + (16*g3^4*mh1*MatMul[Adj[Yb], Yb])/3 - 
  (12*g1^4*mh2*MatMul[Adj[Yb], Yb])/25 + 8*g2^2*g3^2*M2*MatMul[Adj[Yt], ht] + 
  8*g2^2*g3^2*M3*MatMul[Adj[Yt], ht] - 16*g2^2*g3^2*M2^2*
   MatMul[Adj[Yt], Yt] - 16*g2^2*g3^2*M2*M3*MatMul[Adj[Yt], Yt] - 
  16*g2^2*g3^2*M3^2*MatMul[Adj[Yt], Yt] - (24*g1^4*mh1*MatMul[Adj[Yt], Yt])/
   25 - (3911*g1^4*mh2*MatMul[Adj[Yt], Yt])/150 - 
  (59*g1^2*g2^2*mh2*MatMul[Adj[Yt], Yt])/5 - 
  (45*g2^4*mh2*MatMul[Adj[Yt], Yt])/2 - 
  (136*g1^2*g3^2*mh2*MatMul[Adj[Yt], Yt])/5 - 
  8*g2^2*g3^2*mh2*MatMul[Adj[Yt], Yt] + (16*g3^4*mh2*MatMul[Adj[Yt], Yt])/3 + 
  (7*g1^4*MatMul[MatMul[Adj[hb], hb], Zeta[3]])/15 + 
  (18*g1^2*g2^2*MatMul[MatMul[Adj[hb], hb], Zeta[3]])/5 - 
  63*g2^4*MatMul[MatMul[Adj[hb], hb], Zeta[3]] + 
  (256*g1^2*g3^2*MatMul[MatMul[Adj[hb], hb], Zeta[3]])/15 - 
  (544*g3^4*MatMul[MatMul[Adj[hb], hb], Zeta[3]])/3 + 
  M1*((633*g1^4*MatMul[Adj[hb], Yb])/25 - 
    (14*g1^4*MatMul[MatMul[Adj[hb], Yb], Zeta[3]])/15) + 
  M1*((41*g1^2*g2^2*MatMul[Adj[hb], Yb])/5 - 
    (18*g1^2*g2^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]])/5) + 
  M2*((41*g1^2*g2^2*MatMul[Adj[hb], Yb])/5 - 
    (18*g1^2*g2^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]])/5) + 
  M2*(45*g2^4*MatMul[Adj[hb], Yb] + 126*g2^4*MatMul[MatMul[Adj[hb], Yb], 
      Zeta[3]]) + M1*((152*g1^2*g3^2*MatMul[Adj[hb], Yb])/15 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]])/15) + 
  M3*((152*g1^2*g3^2*MatMul[Adj[hb], Yb])/15 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]])/15) + 
  M3*((-32*g3^4*MatMul[Adj[hb], Yb])/3 + 
    (1088*g3^4*MatMul[MatMul[Adj[hb], Yb], Zeta[3]])/3) + 
  (143*g1^4*MatMul[MatMul[Adj[ht], ht], Zeta[3]])/75 + 
  18*g1^2*g2^2*MatMul[MatMul[Adj[ht], ht], Zeta[3]] - 
  63*g2^4*MatMul[MatMul[Adj[ht], ht], Zeta[3]] + 
  (256*g1^2*g3^2*MatMul[MatMul[Adj[ht], ht], Zeta[3]])/15 - 
  (544*g3^4*MatMul[MatMul[Adj[ht], ht], Zeta[3]])/3 + 
  M1*((3767*g1^4*MatMul[Adj[ht], Yt])/75 - 
    (286*g1^4*MatMul[MatMul[Adj[ht], Yt], Zeta[3]])/75) + 
  M1*((59*g1^2*g2^2*MatMul[Adj[ht], Yt])/5 - 
    18*g1^2*g2^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]]) + 
  M2*((59*g1^2*g2^2*MatMul[Adj[ht], Yt])/5 - 
    18*g1^2*g2^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]]) + 
  M2*(45*g2^4*MatMul[Adj[ht], Yt] + 126*g2^4*MatMul[MatMul[Adj[ht], Yt], 
      Zeta[3]]) + M1*((136*g1^2*g3^2*MatMul[Adj[ht], Yt])/5 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]])/15) + 
  M3*((136*g1^2*g3^2*MatMul[Adj[ht], Yt])/5 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]])/15) + 
  M3*((-32*g3^4*MatMul[Adj[ht], Yt])/3 + 
    (1088*g3^4*MatMul[MatMul[Adj[ht], Yt], Zeta[3]])/3) + 
  M1*((633*g1^4*MatMul[Adj[Yb], hb])/25 - 
    (14*g1^4*MatMul[MatMul[Adj[Yb], hb], Zeta[3]])/15) + 
  M1*((41*g1^2*g2^2*MatMul[Adj[Yb], hb])/5 - 
    (18*g1^2*g2^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]])/5) + 
  M2*((41*g1^2*g2^2*MatMul[Adj[Yb], hb])/5 - 
    (18*g1^2*g2^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]])/5) + 
  M2*(45*g2^4*MatMul[Adj[Yb], hb] + 126*g2^4*MatMul[MatMul[Adj[Yb], hb], 
      Zeta[3]]) + M1*((152*g1^2*g3^2*MatMul[Adj[Yb], hb])/15 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]])/15) + 
  M3*((152*g1^2*g3^2*MatMul[Adj[Yb], hb])/15 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]])/15) + 
  M3*((-32*g3^4*MatMul[Adj[Yb], hb])/3 + 
    (1088*g3^4*MatMul[MatMul[Adj[Yb], hb], Zeta[3]])/3) + 
  (7*g1^4*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/15 + 
  (18*g1^2*g2^2*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/5 - 
  63*g2^4*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]] + 
  (256*g1^2*g3^2*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/15 - 
  (544*g3^4*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/3 + 
  M1^2*((-1899*g1^4*MatMul[Adj[Yb], Yb])/25 + 
    (14*g1^4*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/5) + 
  M1^2*((-82*g1^2*g2^2*MatMul[Adj[Yb], Yb])/5 + 
    (36*g1^2*g2^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/5) + 
  M1*M2*((-82*g1^2*g2^2*MatMul[Adj[Yb], Yb])/5 + 
    (36*g1^2*g2^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/5) + 
  M2^2*((-82*g1^2*g2^2*MatMul[Adj[Yb], Yb])/5 + 
    (36*g1^2*g2^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/5) + 
  M2^2*(-135*g2^4*MatMul[Adj[Yb], Yb] - 378*g2^4*MatMul[MatMul[Adj[Yb], Yb], 
      Zeta[3]]) + M1^2*((-304*g1^2*g3^2*MatMul[Adj[Yb], Yb])/15 + 
    (512*g1^2*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/15) + 
  M1*M3*((-304*g1^2*g3^2*MatMul[Adj[Yb], Yb])/15 + 
    (512*g1^2*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/15) + 
  M3^2*((-304*g1^2*g3^2*MatMul[Adj[Yb], Yb])/15 + 
    (512*g1^2*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]])/15) + 
  M3^2*(32*g3^4*MatMul[Adj[Yb], Yb] - 1088*g3^4*MatMul[MatMul[Adj[Yb], Yb], 
      Zeta[3]]) + M1*((3767*g1^4*MatMul[Adj[Yt], ht])/75 - 
    (286*g1^4*MatMul[MatMul[Adj[Yt], ht], Zeta[3]])/75) + 
  M1*((59*g1^2*g2^2*MatMul[Adj[Yt], ht])/5 - 
    18*g1^2*g2^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]]) + 
  M2*((59*g1^2*g2^2*MatMul[Adj[Yt], ht])/5 - 
    18*g1^2*g2^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]]) + 
  M2*(45*g2^4*MatMul[Adj[Yt], ht] + 126*g2^4*MatMul[MatMul[Adj[Yt], ht], 
      Zeta[3]]) + M1*((136*g1^2*g3^2*MatMul[Adj[Yt], ht])/5 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]])/15) + 
  M3*((136*g1^2*g3^2*MatMul[Adj[Yt], ht])/5 - 
    (256*g1^2*g3^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]])/15) + 
  M3*((-32*g3^4*MatMul[Adj[Yt], ht])/3 + 
    (1088*g3^4*MatMul[MatMul[Adj[Yt], ht], Zeta[3]])/3) + 
  (143*g1^4*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/75 + 
  18*g1^2*g2^2*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]] - 
  63*g2^4*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]] + 
  (256*g1^2*g3^2*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/15 - 
  (544*g3^4*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/3 + 
  M1^2*((-3767*g1^4*MatMul[Adj[Yt], Yt])/25 + 
    (286*g1^4*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/25) + 
  M1^2*((-118*g1^2*g2^2*MatMul[Adj[Yt], Yt])/5 + 
    36*g1^2*g2^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]) + 
  M1*M2*((-118*g1^2*g2^2*MatMul[Adj[Yt], Yt])/5 + 
    36*g1^2*g2^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]) + 
  M2^2*((-118*g1^2*g2^2*MatMul[Adj[Yt], Yt])/5 + 
    36*g1^2*g2^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]) + 
  M2^2*(-135*g2^4*MatMul[Adj[Yt], Yt] - 378*g2^4*MatMul[MatMul[Adj[Yt], Yt], 
      Zeta[3]]) + M1^2*((-272*g1^2*g3^2*MatMul[Adj[Yt], Yt])/5 + 
    (512*g1^2*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/15) + 
  M1*M3*((-272*g1^2*g3^2*MatMul[Adj[Yt], Yt])/5 + 
    (512*g1^2*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/15) + 
  M3^2*((-272*g1^2*g3^2*MatMul[Adj[Yt], Yt])/5 + 
    (512*g1^2*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]])/15) + 
  M3^2*(32*g3^4*MatMul[Adj[Yt], Yt] - 1088*g3^4*MatMul[MatMul[Adj[Yt], Yt], 
      Zeta[3]]) + (7*g1^4*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]])/30 + 
  (9*g1^2*g2^2*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]])/5 - 
  (63*g2^4*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]])/3 + 
  (143*g1^4*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]])/150 + 
  9*g1^2*g2^2*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]] - 
  (63*g2^4*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]])/3 + 
  (7*g1^4*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]])/15 + 
  (18*g1^2*g2^2*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]])/5 - 
  63*g2^4*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]] + 
  (256*g1^2*g3^2*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]])/15 - 
  (544*g3^4*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]])/3 + 
  (7*g1^4*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]])/30 + 
  (9*g1^2*g2^2*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]])/5 - 
  (63*g2^4*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]])/3 + 
  (143*g1^4*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]])/75 + 
  18*g1^2*g2^2*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]] - 
  63*g2^4*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]] + 
  (256*g1^2*g3^2*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]])/15 - 
  (544*g3^4*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]])/3 + 
  (143*g1^4*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]])/150 + 
  9*g1^2*g2^2*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]] - 
  (63*g2^4*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]])/2 + 
  (128*g1^2*g3^2*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]])/15 - 
  (272*g3^4*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]])/3 - 
  (12*g1^2*MatMul[MatMul[Adj[hb], hb, Adj[Yb], Yb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[hb], hb, Adj[Yb], Yb], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Adj[hb], Yb, Adj[Yb], hb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[hb], Yb, Adj[Yb], hb], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[ht], ht, Adj[Yt], Yt], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[ht], ht, Adj[Yt], Yt], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[ht], Yt, Adj[Yt], ht], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[ht], Yt, Adj[Yt], ht], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Adj[Yb], hb, Adj[hb], Yb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Yb], hb, Adj[hb], Yb], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[hb], hb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[hb], hb], Zeta[3]] - 
  (24*g1^2*mh1*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
  72*g2^2*mh1*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[Yt], ht, Adj[ht], Yt], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[Yt], ht, Adj[ht], Yt], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[ht], ht], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[ht], ht], Zeta[3]] - 
  24*g1^2*mh2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  72*g2^2*mh2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
  (6*g1^2*MatMul[MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] - 
  6*g1^2*MatMul[MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  18*g2^2*MatMul[MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb], Zeta[3]] - 
  (12*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb], Zeta[3]] - 
  (6*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt], Zeta[3]] - 
  12*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt], Zeta[3]] + 
  36*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt], Zeta[3]] - 
  6*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq], Zeta[3]] + 
  18*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq], Zeta[3]] + 
  12*MatMul[MatMul[Adj[hb], hb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[hb], Yb, Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[hb], Yb, Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[ht], ht, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[ht], Yt, Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[ht], Yt, Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], hb, Adj[hb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], hb, Adj[Yb], Yb, Adj[hb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[hb], hb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[hb], Yb, Adj[Yb], hb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], hb, Adj[hb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[hb], hb], Zeta[3]] + 
  36*mh1*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], ht, Adj[ht], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], ht, Adj[Yt], Yt, Adj[ht], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[ht], ht, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[ht], Yt, Adj[Yt], ht], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], ht, Adj[ht], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[ht], ht], Zeta[3]] + 
  36*mh2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  6*MatMul[MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  6*MatMul[MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq, Adj[Yb], Yb], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], mb, Yb], Zeta[3]] + 
  6*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], Yb, mq], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq, Adj[Yt], Yt], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], mt, Yt], Zeta[3]] + 
  6*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], Yt, mq], Zeta[3]] - 
  (633*g1^4*MatMul[mq, Adj[Yb], Yb])/100 - 
  (41*g1^2*g2^2*MatMul[mq, Adj[Yb], Yb])/10 - 
  (45*g2^4*MatMul[mq, Adj[Yb], Yb])/4 - 
  (76*g1^2*g3^2*MatMul[mq, Adj[Yb], Yb])/15 - 
  4*g2^2*g3^2*MatMul[mq, Adj[Yb], Yb] + (8*g3^4*MatMul[mq, Adj[Yb], Yb])/3 - 
  (3767*g1^4*MatMul[mq, Adj[Yt], Yt])/300 - 
  (59*g1^2*g2^2*MatMul[mq, Adj[Yt], Yt])/10 - 
  (45*g2^4*MatMul[mq, Adj[Yt], Yt])/4 - 
  (68*g1^2*g3^2*MatMul[mq, Adj[Yt], Yt])/5 - 
  4*g2^2*g3^2*MatMul[mq, Adj[Yt], Yt] + (8*g3^4*MatMul[mq, Adj[Yt], Yt])/3 - 
  (633*g1^4*MatMul[Adj[Yb], mb, Yb])/50 - 
  (41*g1^2*g2^2*MatMul[Adj[Yb], mb, Yb])/5 - 
  (45*g2^4*MatMul[Adj[Yb], mb, Yb])/2 - 
  (152*g1^2*g3^2*MatMul[Adj[Yb], mb, Yb])/15 - 
  8*g2^2*g3^2*MatMul[Adj[Yb], mb, Yb] + (16*g3^4*MatMul[Adj[Yb], mb, Yb])/3 - 
  (633*g1^4*MatMul[Adj[Yb], Yb, mq])/100 - 
  (41*g1^2*g2^2*MatMul[Adj[Yb], Yb, mq])/10 - 
  (45*g2^4*MatMul[Adj[Yb], Yb, mq])/4 - 
  (76*g1^2*g3^2*MatMul[Adj[Yb], Yb, mq])/15 - 
  4*g2^2*g3^2*MatMul[Adj[Yb], Yb, mq] + (8*g3^4*MatMul[Adj[Yb], Yb, mq])/3 - 
  (3767*g1^4*MatMul[Adj[Yt], mt, Yt])/150 - 
  (59*g1^2*g2^2*MatMul[Adj[Yt], mt, Yt])/5 - 
  (45*g2^4*MatMul[Adj[Yt], mt, Yt])/2 - 
  (136*g1^2*g3^2*MatMul[Adj[Yt], mt, Yt])/5 - 
  8*g2^2*g3^2*MatMul[Adj[Yt], mt, Yt] + (16*g3^4*MatMul[Adj[Yt], mt, Yt])/3 - 
  (3767*g1^4*MatMul[Adj[Yt], Yt, mq])/300 - 
  (59*g1^2*g2^2*MatMul[Adj[Yt], Yt, mq])/10 - 
  (45*g2^4*MatMul[Adj[Yt], Yt, mq])/4 - 
  (68*g1^2*g3^2*MatMul[Adj[Yt], Yt, mq])/5 - 
  4*g2^2*g3^2*MatMul[Adj[Yt], Yt, mq] + (8*g3^4*MatMul[Adj[Yt], Yt, mq])/3 + 
  (14*g1^2*MatMul[Adj[hb], hb, Adj[Yb], Yb])/15 - 
  6*g2^2*MatMul[Adj[hb], hb, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[Adj[hb], hb, Adj[Yb], Yb])/3 + 
  (14*g1^2*MatMul[Adj[hb], Yb, Adj[Yb], hb])/15 - 
  6*g2^2*MatMul[Adj[hb], Yb, Adj[Yb], hb] + 
  (128*g3^2*MatMul[Adj[hb], Yb, Adj[Yb], hb])/3 - 
  (128*g3^2*M3*MatMul[Adj[hb], Yb, Adj[Yb], Yb])/3 + 
  M1*((12*g1^2*MatMul[MatMul[Adj[hb], Yb, Adj[Yb], Yb], Zeta[3]])/5 - 
    (14*g1^2*MatMul[Adj[hb], Yb, Adj[Yb], Yb])/15) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[hb], Yb, Adj[Yb], Yb], Zeta[3]] + 
    6*g2^2*MatMul[Adj[hb], Yb, Adj[Yb], Yb]) + 
  (22*g1^2*MatMul[Adj[ht], ht, Adj[Yt], Yt])/3 - 
  6*g2^2*MatMul[Adj[ht], ht, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[Adj[ht], ht, Adj[Yt], Yt])/3 + 
  (22*g1^2*MatMul[Adj[ht], Yt, Adj[Yt], ht])/3 - 
  6*g2^2*MatMul[Adj[ht], Yt, Adj[Yt], ht] + 
  (128*g3^2*MatMul[Adj[ht], Yt, Adj[Yt], ht])/3 - 
  (128*g3^2*M3*MatMul[Adj[ht], Yt, Adj[Yt], Yt])/3 + 
  M1*(12*g1^2*MatMul[MatMul[Adj[ht], Yt, Adj[Yt], Yt], Zeta[3]] - 
    (22*g1^2*MatMul[Adj[ht], Yt, Adj[Yt], Yt])/3) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[ht], Yt, Adj[Yt], Yt], Zeta[3]] + 
    6*g2^2*MatMul[Adj[ht], Yt, Adj[Yt], Yt]) + 
  (14*g1^2*MatMul[Adj[Yb], hb, Adj[hb], Yb])/15 - 
  6*g2^2*MatMul[Adj[Yb], hb, Adj[hb], Yb] + 
  (128*g3^2*MatMul[Adj[Yb], hb, Adj[hb], Yb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], hb, Adj[Yb], Yb])/3 + 
  M1*((12*g1^2*MatMul[MatMul[Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]])/5 - 
    (14*g1^2*MatMul[Adj[Yb], hb, Adj[Yb], Yb])/15) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Yb], hb, Adj[Yb], Yb], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Yb], hb, Adj[Yb], Yb]) + 
  (14*g1^2*MatMul[Adj[Yb], Yb, Adj[hb], hb])/15 - 
  6*g2^2*MatMul[Adj[Yb], Yb, Adj[hb], hb] + 
  (128*g3^2*MatMul[Adj[Yb], Yb, Adj[hb], hb])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yb], Yb, Adj[hb], Yb])/3 + 
  M1*((12*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[hb], Yb], Zeta[3]])/5 - 
    (14*g1^2*MatMul[Adj[Yb], Yb, Adj[hb], Yb])/15) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[hb], Yb], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Yb], Yb, Adj[hb], Yb]) - 
  (128*g3^2*M3*MatMul[Adj[Yb], Yb, Adj[Yb], hb])/3 + 
  M1*((12*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]])/5 - 
    (14*g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], hb])/15) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], hb], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], hb]) + 
  (256*g3^2*M3^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (28*g1^2*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/15 - 
  12*g2^2*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb] + 
  (256*g3^2*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  M1^2*((-24*g1^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]])/5 + 
    (28*g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb])/15) + 
  M2^2*(72*g2^2*MatMul[MatMul[Adj[Yb], Yb, Adj[Yb], Yb], Zeta[3]] - 
    12*g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]) + 
  (22*g1^2*MatMul[Adj[Yt], ht, Adj[ht], Yt])/3 - 
  6*g2^2*MatMul[Adj[Yt], ht, Adj[ht], Yt] + 
  (128*g3^2*MatMul[Adj[Yt], ht, Adj[ht], Yt])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], ht, Adj[Yt], Yt])/3 + 
  M1*(12*g1^2*MatMul[MatMul[Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] - 
    (22*g1^2*MatMul[Adj[Yt], ht, Adj[Yt], Yt])/3) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Yt], ht, Adj[Yt], Yt], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Yt], ht, Adj[Yt], Yt]) + 
  (22*g1^2*MatMul[Adj[Yt], Yt, Adj[ht], ht])/3 - 
  6*g2^2*MatMul[Adj[Yt], Yt, Adj[ht], ht] + 
  (128*g3^2*MatMul[Adj[Yt], Yt, Adj[ht], ht])/3 - 
  (128*g3^2*M3*MatMul[Adj[Yt], Yt, Adj[ht], Yt])/3 + 
  M1*(12*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[ht], Yt], Zeta[3]] - 
    (22*g1^2*MatMul[Adj[Yt], Yt, Adj[ht], Yt])/3) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[ht], Yt], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Yt], Yt, Adj[ht], Yt]) - 
  (128*g3^2*M3*MatMul[Adj[Yt], Yt, Adj[Yt], ht])/3 + 
  M1*(12*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] - 
    (22*g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], ht])/3) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], ht], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], ht]) + 
  (256*g3^2*M3^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (44*g1^2*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 - 
  12*g2^2*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt] + 
  (256*g3^2*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  M1^2*(-24*g1^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] + 
    (44*g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt])/3) + 
  M2^2*(72*g2^2*MatMul[MatMul[Adj[Yt], Yt, Adj[Yt], Yt], Zeta[3]] - 
    12*g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]) + 
  (7*g1^2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb])/15 - 
  3*g2^2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb] + 
  (64*g3^2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb])/3 + 
  (11*g1^2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt])/3 - 
  3*g2^2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt] + 
  (64*g3^2*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt])/3 + 
  (14*g1^2*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb])/15 - 
  6*g2^2*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb])/3 + 
  (14*g1^2*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb])/15 - 
  6*g2^2*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb] + 
  (128*g3^2*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb])/3 + 
  (14*g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb])/15 - 
  6*g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb] + 
  (128*g3^2*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb])/3 + 
  (7*g1^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq])/15 - 
  3*g2^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  (64*g3^2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq])/3 + 
  (22*g1^2*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt])/3 - 
  6*g2^2*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt])/3 + 
  (22*g1^2*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt])/3 - 
  6*g2^2*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt] + 
  (128*g3^2*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt])/3 + 
  (22*g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt])/3 - 
  6*g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt] + 
  (128*g3^2*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt])/3 + 
  (11*g1^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq])/3 - 
  3*g2^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  (64*g3^2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq])/3 + 
  8*MatMul[Adj[hb], hb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[hb], Yb, Adj[Yt], ht, Adj[Yb], Yb] + 
  8*MatMul[Adj[hb], Yb, Adj[Yt], Yt, Adj[Yb], hb] + 
  8*MatMul[Adj[ht], ht, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[ht], Yt, Adj[Yb], hb, Adj[Yt], Yt] + 
  8*MatMul[Adj[ht], Yt, Adj[Yb], Yb, Adj[Yt], ht] + 
  8*MatMul[Adj[Yb], hb, Adj[ht], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], hb, Adj[Yt], Yt, Adj[hb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[ht], ht, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[ht], Yt, Adj[Yb], hb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], ht, Adj[hb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[hb], hb] + 
  16*mh1*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*mh2*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yt], ht, Adj[hb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], ht, Adj[Yb], Yb, Adj[ht], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[hb], hb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[hb], Yb, Adj[Yt], ht] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], hb, Adj[ht], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[ht], ht] + 
  8*mh1*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  16*mh2*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  4*MatMul[mq, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  4*MatMul[mq, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yb], mb, Yb, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, mq, Adj[Yt], Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], mt, Yt, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, mq, Adj[Yb], Yb] + 
  8*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], mb, Yb] + 
  4*MatMul[Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  8*MatMul[Adj[Yt], mt, Yt, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, mq, Adj[Yb], Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], mb, Yb, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, mq, Adj[Yt], Yt] + 
  8*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], mt, Yt] + 
  4*MatMul[Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], Yt, mq] - 
  (398*g1^6*trace[mb])/1125 - (2*g1^4*g2^2*trace[mb])/25 - 
  (32*g1^4*g3^2*trace[mb])/225 - (16*g1^2*g3^4*trace[mb])/45 - 
  16*g2^2*g3^4*trace[mb] + (320*g3^6*trace[mb])/9 - 
  (8*g1^4*MatMul[Adj[Yb], Yb]*trace[mb])/25 - 
  (16*g1^4*MatMul[Adj[Yt], Yt]*trace[mb])/25 - (398*g1^6*trace[me])/375 - 
  (6*g1^4*g2^2*trace[me])/25 - (32*g1^4*g3^2*trace[me])/75 - 
  (24*g1^4*MatMul[Adj[Yb], Yb]*trace[me])/25 - 
  (48*g1^4*MatMul[Adj[Yt], Yt]*trace[me])/25 - (199*g1^6*trace[ml])/375 - 
  (3*g1^4*g2^2*trace[ml])/25 - (g1^2*g2^4*trace[ml])/5 - 3*g2^6*trace[ml] - 
  (16*g1^4*g3^2*trace[ml])/75 - 16*g2^4*g3^2*trace[ml] - 
  (12*g1^4*MatMul[Adj[Yb], Yb]*trace[ml])/25 - 
  (24*g1^4*MatMul[Adj[Yt], Yt]*trace[ml])/25 - (199*g1^6*trace[mq])/1125 - 
  (g1^4*g2^2*trace[mq])/25 - (3*g1^2*g2^4*trace[mq])/5 - 9*g2^6*trace[mq] - 
  (16*g1^4*g3^2*trace[mq])/225 - 48*g2^4*g3^2*trace[mq] - 
  (32*g1^2*g3^4*trace[mq])/45 - 32*g2^2*g3^4*trace[mq] + 
  (640*g3^6*trace[mq])/9 - (4*g1^4*MatMul[Adj[Yb], Yb]*trace[mq])/25 - 
  (8*g1^4*MatMul[Adj[Yt], Yt]*trace[mq])/25 - (1592*g1^6*trace[mt])/1125 - 
  (8*g1^4*g2^2*trace[mt])/25 - (128*g1^4*g3^2*trace[mt])/225 - 
  (16*g1^2*g3^4*trace[mt])/45 - 16*g2^2*g3^4*trace[mt] + 
  (320*g3^6*trace[mt])/9 - (32*g1^4*MatMul[Adj[Yb], Yb]*trace[mt])/25 - 
  (64*g1^4*MatMul[Adj[Yt], Yt]*trace[mt])/25 - (14*g1^4*trace[hb, Adj[hb]])/
   25 - 54*g2^4*trace[hb, Adj[hb]] - 64*g3^4*trace[hb, Adj[hb]] + 
  (32*g1^2*MatMul[Adj[Yb], Yb]*trace[hb, Adj[hb]])/5 + 
  36*g2^2*MatMul[Adj[Yb], Yb]*trace[hb, Adj[hb]] - 
  16*g3^2*MatMul[Adj[Yb], Yb]*trace[hb, Adj[hb]] - 
  (48*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[hb]])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[hb]] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[hb, Adj[hb]] + 
  (32*g1^2*MatMul[Adj[hb], Yb]*trace[hb, Adj[Yb]])/5 + 
  36*g2^2*MatMul[Adj[hb], Yb]*trace[hb, Adj[Yb]] - 
  16*g3^2*MatMul[Adj[hb], Yb]*trace[hb, Adj[Yb]] - 
  (48*g1^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]]*trace[hb, Adj[Yb]])/5 + 
  96*g3^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]]*trace[hb, Adj[Yb]] + 
  12*MatMul[Adj[hb], Yb, Adj[Yb], Yb]*trace[hb, Adj[Yb]] + 
  12*MatMul[Adj[Yb], Yb, Adj[hb], Yb]*trace[hb, Adj[Yb]] - 
  (18*g1^4*trace[he, Adj[he]])/25 - 18*g2^4*trace[he, Adj[he]] - 
  (16*g1^2*MatMul[Adj[Yb], Yb]*trace[he, Adj[he]])/5 + 
  12*g2^2*MatMul[Adj[Yb], Yb]*trace[he, Adj[he]] + 
  16*g3^2*MatMul[Adj[Yb], Yb]*trace[he, Adj[he]] + 
  (24*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[he, Adj[he]])/5 + 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[he, Adj[he]] - 
  (16*g1^2*MatMul[Adj[hb], Yb]*trace[he, Adj[Ye]])/5 + 
  12*g2^2*MatMul[Adj[hb], Yb]*trace[he, Adj[Ye]] + 
  16*g3^2*MatMul[Adj[hb], Yb]*trace[he, Adj[Ye]] + 
  (24*g1^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]]*trace[he, Adj[Ye]])/5 + 
  4*MatMul[Adj[hb], Yb, Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
  4*MatMul[Adj[Yb], Yb, Adj[hb], Yb]*trace[he, Adj[Ye]] - 
  (26*g1^4*trace[ht, Adj[ht]])/25 - 54*g2^4*trace[ht, Adj[ht]] - 
  64*g3^4*trace[ht, Adj[ht]] + 4*g1^2*MatMul[Adj[Yt], Yt]*
   trace[ht, Adj[ht]] + 36*g2^2*MatMul[Adj[Yt], Yt]*trace[ht, Adj[ht]] - 
  16*g3^2*MatMul[Adj[Yt], Yt]*trace[ht, Adj[ht]] - 
  (48*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[ht]])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[ht]] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*trace[ht, Adj[ht]] + 
  4*g1^2*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]] + 
  36*g2^2*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]] - 
  16*g3^2*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]] - 
  (48*g1^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]]*trace[ht, Adj[Yt]])/5 + 
  96*g3^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]]*trace[ht, Adj[Yt]] + 
  12*MatMul[Adj[ht], Yt, Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
  12*MatMul[Adj[Yt], Yt, Adj[ht], Yt]*trace[ht, Adj[Yt]] + 
  (32*g1^2*MatMul[Adj[Yb], hb]*trace[Adj[hb], Yb])/5 + 
  36*g2^2*MatMul[Adj[Yb], hb]*trace[Adj[hb], Yb] - 
  16*g3^2*MatMul[Adj[Yb], hb]*trace[Adj[hb], Yb] - 
  (48*g1^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]]*trace[Adj[hb], Yb])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]]*trace[Adj[hb], Yb] + 
  12*MatMul[Adj[Yb], hb, Adj[Yb], Yb]*trace[Adj[hb], Yb] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], hb]*trace[Adj[hb], Yb] - 
  36*MatMul[Adj[Yb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  12*MatMul[Adj[Yb], Yb]*trace[he, Adj[Ye]]*trace[Adj[hb], Yb] - 
  (16*g1^2*MatMul[Adj[Yb], hb]*trace[Adj[he], Ye])/5 + 
  12*g2^2*MatMul[Adj[Yb], hb]*trace[Adj[he], Ye] + 
  16*g3^2*MatMul[Adj[Yb], hb]*trace[Adj[he], Ye] + 
  (24*g1^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]]*trace[Adj[he], Ye])/5 + 
  4*MatMul[Adj[Yb], hb, Adj[Yb], Yb]*trace[Adj[he], Ye] + 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], hb]*trace[Adj[he], Ye] - 
  12*MatMul[Adj[Yb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[he], Ye] - 
  4*MatMul[Adj[Yb], Yb]*trace[he, Adj[Ye]]*trace[Adj[he], Ye] + 
  M2*(-36*g2^2*MatMul[Adj[Yb], Yb]*trace[hb, Adj[Yb]] - 
    12*g2^2*MatMul[Adj[Yb], Yb]*trace[he, Adj[Ye]] - 
    36*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[hb], Yb] - 
    12*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[he], Ye]) + 
  M3*(16*g3^2*MatMul[Adj[Yb], Yb]*trace[hb, Adj[Yb]] - 
    96*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]] - 
    16*g3^2*MatMul[Adj[Yb], Yb]*trace[he, Adj[Ye]] + 
    16*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[hb], Yb] - 
    96*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[hb], Yb] - 
    16*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[he], Ye]) + 
  M1*((-32*g1^2*MatMul[Adj[Yb], Yb]*trace[hb, Adj[Yb]])/5 + 
    (48*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[hb, Adj[Yb]])/5 + 
    (16*g1^2*MatMul[Adj[Yb], Yb]*trace[he, Adj[Ye]])/5 - 
    (24*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[he, Adj[Ye]])/5 - 
    (32*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[hb], Yb])/5 + 
    (48*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[hb], Yb])/5 + 
    (16*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[he], Ye])/5 - 
    (24*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[he], Ye])/5) + 
  4*g1^2*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt] + 
  36*g2^2*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt] - 
  16*g3^2*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt] - 
  (48*g1^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]]*trace[Adj[ht], Yt])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]]*trace[Adj[ht], Yt] + 
  12*MatMul[Adj[Yt], ht, Adj[Yt], Yt]*trace[Adj[ht], Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], ht]*trace[Adj[ht], Yt] - 
  36*MatMul[Adj[Yt], Yt]*trace[ht, Adj[Yt]]*trace[Adj[ht], Yt] + 
  M1*((14*g1^4*trace[hb, Adj[Yb]])/15 + (6*g1^4*trace[he, Adj[Ye]])/5 + 
    (26*g1^4*trace[ht, Adj[Yt]])/15 + (14*g1^4*trace[Adj[hb], Yb])/15 + 
    (6*g1^4*trace[Adj[he], Ye])/5 + (26*g1^4*trace[Adj[ht], Yt])/15) + 
  M2*(90*g2^4*trace[hb, Adj[Yb]] + 30*g2^4*trace[he, Adj[Ye]] + 
    90*g2^4*trace[ht, Adj[Yt]] + 90*g2^4*trace[Adj[hb], Yb] + 
    30*g2^4*trace[Adj[he], Ye] + 90*g2^4*trace[Adj[ht], Yt]) + 
  M3*((320*g3^4*trace[hb, Adj[Yb]])/3 + (320*g3^4*trace[ht, Adj[Yt]])/3 + 
    (320*g3^4*trace[Adj[hb], Yb])/3 + (320*g3^4*trace[Adj[ht], Yt])/3) + 
  M2*(-36*g2^2*MatMul[Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
    36*g2^2*MatMul[Adj[Yt], Yt]*trace[Adj[ht], Yt]) + 
  M1*(-4*g1^2*MatMul[Adj[Yt], Yt]*trace[ht, Adj[Yt]] + 
    (48*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]])/5 - 
    4*g1^2*MatMul[Adj[Yt], Yt]*trace[Adj[ht], Yt] + 
    (48*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[ht], Yt])/5) + 
  M3*(16*g3^2*MatMul[Adj[Yt], Yt]*trace[ht, Adj[Yt]] - 
    96*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[ht, Adj[Yt]] + 
    16*g3^2*MatMul[Adj[Yt], Yt]*trace[Adj[ht], Yt] - 
    96*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[ht], Yt]) - 
  (14*g1^4*mh1*trace[Adj[Yb], Yb])/25 - 54*g2^4*mh1*trace[Adj[Yb], Yb] - 
  64*g3^4*mh1*trace[Adj[Yb], Yb] + 
  (32*g1^2*MatMul[Adj[hb], hb]*trace[Adj[Yb], Yb])/5 + 
  36*g2^2*MatMul[Adj[hb], hb]*trace[Adj[Yb], Yb] - 
  16*g3^2*MatMul[Adj[hb], hb]*trace[Adj[Yb], Yb] + 
  (64*g1^2*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 + 
  72*g2^2*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  32*g3^2*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  (48*g1^2*MatMul[MatMul[Adj[hb], hb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  96*g3^2*MatMul[MatMul[Adj[hb], hb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (96*g1^2*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  192*g3^2*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (24*g1^2*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  48*g3^2*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (48*g1^2*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  (24*g1^2*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
  48*g3^2*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (16*g1^2*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 + 
  18*g2^2*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb] - 
  8*g3^2*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  (32*g1^2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb])/5 + 
  36*g2^2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb] - 
  16*g3^2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb] + 
  (16*g1^2*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb])/5 + 
  18*g2^2*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb] - 
  8*g3^2*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[hb], hb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[hb], Yb, Adj[Yb], hb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Yb], hb, Adj[hb], Yb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Yb], Yb, Adj[hb], hb]*trace[Adj[Yb], Yb] + 
  36*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  6*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb] + 
  6*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb] - 
  36*MatMul[Adj[Yb], Yb]*trace[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  36*MatMul[Adj[hb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Adj[Yb], Yb]*trace[he, Adj[he]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Adj[hb], Yb]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  36*MatMul[Adj[Yb], hb]*trace[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 
  12*MatMul[Adj[Yb], hb]*trace[Adj[he], Ye]*trace[Adj[Yb], Yb] - 
  18*MatMul[Adj[hb], hb]*trace[Adj[Yb], Yb]^2 - 
  54*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]^2 - 
  9*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb]^2 - 
  9*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb]^2 - 
  (18*g1^4*mh1*trace[Adj[Ye], Ye])/25 - 18*g2^4*mh1*trace[Adj[Ye], Ye] - 
  (16*g1^2*MatMul[Adj[hb], hb]*trace[Adj[Ye], Ye])/5 + 
  12*g2^2*MatMul[Adj[hb], hb]*trace[Adj[Ye], Ye] + 
  16*g3^2*MatMul[Adj[hb], hb]*trace[Adj[Ye], Ye] - 
  (32*g1^2*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
  24*g2^2*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  32*g3^2*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  (24*g1^2*MatMul[MatMul[Adj[hb], hb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
  (48*g1^2*mh1*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
  (12*g1^2*MatMul[MatMul[mq, Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
  (24*g1^2*MatMul[MatMul[Adj[Yb], mb, Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5 + 
  (12*g1^2*MatMul[MatMul[Adj[Yb], Yb, mq], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  (8*g1^2*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
  6*g2^2*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  8*g3^2*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  (16*g1^2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Ye], Ye])/5 + 
  12*g2^2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Ye], Ye] + 
  16*g3^2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Ye], Ye] - 
  (8*g1^2*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye])/5 + 
  6*g2^2*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye] + 
  8*g3^2*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[hb], hb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[hb], Yb, Adj[Yb], hb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Yb], hb, Adj[hb], Yb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Yb], Yb, Adj[hb], hb]*trace[Adj[Ye], Ye] + 
  12*mh1*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  2*MatMul[mq, Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Yb], mb, Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Yb], Yb, mq, Adj[Yb], Yb]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], mb, Yb]*trace[Adj[Ye], Ye] + 
  2*MatMul[Adj[Yb], Yb, Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[Yb], Yb]*trace[hb, Adj[hb]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[hb], Yb]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Adj[Yb], Yb]*trace[he, Adj[he]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Adj[hb], Yb]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[Yb], hb]*trace[Adj[hb], Yb]*trace[Adj[Ye], Ye] - 
  4*MatMul[Adj[Yb], hb]*trace[Adj[he], Ye]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[hb], hb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  36*mh1*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  6*MatMul[mq, Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  6*MatMul[Adj[Yb], Yb, mq]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  2*MatMul[Adj[hb], hb]*trace[Adj[Ye], Ye]^2 - 6*mh1*MatMul[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye]^2 - MatMul[mq, Adj[Yb], Yb]*trace[Adj[Ye], Ye]^2 - 
  2*MatMul[Adj[Yb], mb, Yb]*trace[Adj[Ye], Ye]^2 - 
  MatMul[Adj[Yb], Yb, mq]*trace[Adj[Ye], Ye]^2 + 
  M2*(-36*g2^2*MatMul[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 
    12*g2^2*MatMul[Adj[hb], Yb]*trace[Adj[Ye], Ye]) + 
  M3*(16*g3^2*MatMul[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 
    96*g3^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] - 
    16*g3^2*MatMul[Adj[hb], Yb]*trace[Adj[Ye], Ye]) + 
  M2*(-36*g2^2*MatMul[Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
    12*g2^2*MatMul[Adj[Yb], hb]*trace[Adj[Ye], Ye]) + 
  M3*(16*g3^2*MatMul[Adj[Yb], hb]*trace[Adj[Yb], Yb] - 
    96*g3^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb] - 
    16*g3^2*MatMul[Adj[Yb], hb]*trace[Adj[Ye], Ye]) + 
  M2^2*(72*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
    24*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]) + 
  M3^2*(-32*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb] + 
    192*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb] + 
    32*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]) + 
  M1*((-32*g1^2*MatMul[Adj[hb], Yb]*trace[Adj[Yb], Yb])/5 + 
    (48*g1^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
    (16*g1^2*MatMul[Adj[hb], Yb]*trace[Adj[Ye], Ye])/5 - 
    (24*g1^2*MatMul[MatMul[Adj[hb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M1*((-32*g1^2*MatMul[Adj[Yb], hb]*trace[Adj[Yb], Yb])/5 + 
    (48*g1^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
    (16*g1^2*MatMul[Adj[Yb], hb]*trace[Adj[Ye], Ye])/5 - 
    (24*g1^2*MatMul[MatMul[Adj[Yb], hb], Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M1^2*((64*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb])/5 - 
    (96*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
    (32*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye])/5 + 
    (48*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye])/5) - 
  (26*g1^4*mh2*trace[Adj[Yt], Yt])/25 - 54*g2^4*mh2*trace[Adj[Yt], Yt] - 
  64*g3^4*mh2*trace[Adj[Yt], Yt] + 4*g1^2*MatMul[Adj[ht], ht]*
   trace[Adj[Yt], Yt] + 36*g2^2*MatMul[Adj[ht], ht]*trace[Adj[Yt], Yt] - 
  16*g3^2*MatMul[Adj[ht], ht]*trace[Adj[Yt], Yt] - 
  36*g2^2*M2*MatMul[Adj[ht], Yt]*trace[Adj[Yt], Yt] - 
  36*g2^2*M2*MatMul[Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  72*g2^2*M2^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  8*g1^2*mh2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  72*g2^2*mh2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  32*g3^2*mh2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  (48*g1^2*MatMul[MatMul[Adj[ht], ht], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  96*g3^2*MatMul[MatMul[Adj[ht], ht], Zeta[3]]*trace[Adj[Yt], Yt] - 
  (96*g1^2*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  192*g3^2*mh2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
  (24*g1^2*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  48*g3^2*MatMul[MatMul[mq, Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
  (48*g1^2*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yt], mt, Yt], Zeta[3]]*trace[Adj[Yt], Yt] - 
  (24*g1^2*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]]*trace[Adj[Yt], Yt])/5 + 
  48*g3^2*MatMul[MatMul[Adj[Yt], Yt, mq], Zeta[3]]*trace[Adj[Yt], Yt] + 
  2*g1^2*MatMul[mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
  8*g3^2*MatMul[mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  4*g1^2*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] + 
  36*g2^2*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] - 
  16*g3^2*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] + 
  2*g1^2*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] + 
  18*g2^2*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] - 
  8*g3^2*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[ht], ht, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[ht], Yt, Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], ht, Adj[ht], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[ht], ht]*trace[Adj[Yt], Yt] + 
  36*mh2*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  6*MatMul[mq, Adj[Yt], Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], mt, Yt, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Adj[Yt], Yt, Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt] + 
  6*MatMul[Adj[Yt], Yt, Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt] - 
  36*MatMul[Adj[Yt], Yt]*trace[ht, Adj[ht]]*trace[Adj[Yt], Yt] - 
  36*MatMul[Adj[ht], Yt]*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  36*MatMul[Adj[Yt], ht]*trace[Adj[ht], Yt]*trace[Adj[Yt], Yt] - 
  18*MatMul[Adj[ht], ht]*trace[Adj[Yt], Yt]^2 - 
  54*mh2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt]^2 - 
  9*MatMul[mq, Adj[Yt], Yt]*trace[Adj[Yt], Yt]^2 - 
  18*MatMul[Adj[Yt], mt, Yt]*trace[Adj[Yt], Yt]^2 - 
  9*MatMul[Adj[Yt], Yt, mq]*trace[Adj[Yt], Yt]^2 + 
  M1^2*((-14*g1^4*trace[Adj[Yb], Yb])/5 - (18*g1^4*trace[Adj[Ye], Ye])/5 - 
    (26*g1^4*trace[Adj[Yt], Yt])/5) + 
  M2^2*(-270*g2^4*trace[Adj[Yb], Yb] - 90*g2^4*trace[Adj[Ye], Ye] - 
    270*g2^4*trace[Adj[Yt], Yt]) + M3^2*(-320*g3^4*trace[Adj[Yb], Yb] - 
    320*g3^4*trace[Adj[Yt], Yt]) + 
  M1*(-4*g1^2*MatMul[Adj[ht], Yt]*trace[Adj[Yt], Yt] + 
    (48*g1^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M3*(16*g3^2*MatMul[Adj[ht], Yt]*trace[Adj[Yt], Yt] - 
    96*g3^2*MatMul[MatMul[Adj[ht], Yt], Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M1*(-4*g1^2*MatMul[Adj[Yt], ht]*trace[Adj[Yt], Yt] + 
    (48*g1^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M3*(16*g3^2*MatMul[Adj[Yt], ht]*trace[Adj[Yt], Yt] - 
    96*g3^2*MatMul[MatMul[Adj[Yt], ht], Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M1^2*(8*g1^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] - 
    (96*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M3^2*(-32*g3^2*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt] + 
    192*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*trace[Adj[Yt], Yt]) - 
  (14*g1^4*trace[Yb, Adj[Yb], mb])/25 - 54*g2^4*trace[Yb, Adj[Yb], mb] - 
  64*g3^4*trace[Yb, Adj[Yb], mb] + 
  (32*g1^2*MatMul[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb])/5 + 
  36*g2^2*MatMul[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  16*g3^2*MatMul[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  (48*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Yb, Adj[Yb], mb])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Yb, Adj[Yb], mb] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  36*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  12*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
  (18*g1^4*trace[Ye, Adj[Ye], me])/25 - 18*g2^4*trace[Ye, Adj[Ye], me] - 
  (16*g1^2*MatMul[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me])/5 + 
  12*g2^2*MatMul[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] + 
  16*g3^2*MatMul[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] + 
  (24*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Ye, Adj[Ye], me])/5 + 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
  12*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
  4*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]*trace[Ye, Adj[Ye], me] - 
  (26*g1^4*trace[Yt, Adj[Yt], mt])/25 - 54*g2^4*trace[Yt, Adj[Yt], mt] - 
  64*g3^4*trace[Yt, Adj[Yt], mt] + 4*g1^2*MatMul[Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], mt] + 36*g2^2*MatMul[Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], mt] - 16*g3^2*MatMul[Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], mt] - (48*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*
    trace[Yt, Adj[Yt], mt])/5 + 96*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*
   trace[Yt, Adj[Yt], mt] + 12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], mt] - 36*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], mt] - (14*g1^4*trace[Adj[Yb], Yb, mq])/25 - 
  54*g2^4*trace[Adj[Yb], Yb, mq] - 64*g3^4*trace[Adj[Yb], Yb, mq] + 
  (32*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq])/5 + 
  36*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
  16*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
  (48*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb, mq])/5 + 
  96*g3^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Yb], Yb, mq] + 
  12*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
  36*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
  12*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, mq] - 
  (18*g1^4*trace[Adj[Ye], Ye, ml])/25 - 18*g2^4*trace[Adj[Ye], Ye, ml] - 
  (16*g1^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml])/5 + 
  12*g2^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] + 
  16*g3^2*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] + 
  (24*g1^2*MatMul[MatMul[Adj[Yb], Yb], Zeta[3]]*trace[Adj[Ye], Ye, ml])/5 + 
  4*MatMul[Adj[Yb], Yb, Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
  12*MatMul[Adj[Yb], Yb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
  4*MatMul[Adj[Yb], Yb]*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] - 
  (26*g1^4*trace[Adj[Yt], Yt, mq])/25 - 54*g2^4*trace[Adj[Yt], Yt, mq] - 
  64*g3^4*trace[Adj[Yt], Yt, mq] + 4*g1^2*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, mq] + 36*g2^2*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, mq] - 16*g3^2*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, mq] - (48*g1^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*
    trace[Adj[Yt], Yt, mq])/5 + 96*g3^2*MatMul[MatMul[Adj[Yt], Yt], Zeta[3]]*
   trace[Adj[Yt], Yt, mq] + 12*MatMul[Adj[Yt], Yt, Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, mq] - 36*MatMul[Adj[Yt], Yt]*trace[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, mq] + 12*MatMul[Adj[Yb], Yb]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 12*MatMul[Adj[Yt], Yt]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 12*MatMul[Adj[hb], Yb]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 12*MatMul[Adj[ht], Yt]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 72*MatMul[Adj[Yb], Yb]*
   trace[Adj[hb], hb, Adj[Yb], Yb] + 12*MatMul[Adj[Yb], Yb]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 12*MatMul[Adj[Yt], Yt]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 24*MatMul[Adj[Yb], Yb]*
   trace[Adj[he], he, Adj[Ye], Ye] + 12*MatMul[Adj[Yb], Yb]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 12*MatMul[Adj[Yt], Yt]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 72*MatMul[Adj[Yt], Yt]*
   trace[Adj[ht], ht, Adj[Yt], Yt] + 12*MatMul[Adj[Yb], hb]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 12*MatMul[Adj[Yt], ht]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 72*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yb], hb, Adj[hb], Yb] + 72*MatMul[Adj[hb], Yb]*
   trace[Adj[Yb], hb, Adj[Yb], Yb] + 72*MatMul[Adj[Yb], hb]*
   trace[Adj[Yb], Yb, Adj[hb], Yb] + 36*MatMul[Adj[hb], hb]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 108*mh1*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 18*MatMul[mq, Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 36*MatMul[Adj[Yb], mb, Yb]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 18*MatMul[Adj[Yb], Yb, mq]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 24*MatMul[Adj[Yb], Yb]*
   trace[Adj[Ye], he, Adj[he], Ye] + 24*MatMul[Adj[hb], Yb]*
   trace[Adj[Ye], he, Adj[Ye], Ye] + 24*MatMul[Adj[Yb], hb]*
   trace[Adj[Ye], Ye, Adj[he], Ye] + 12*MatMul[Adj[hb], hb]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 36*mh1*MatMul[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*MatMul[mq, Adj[Yb], Yb]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 12*MatMul[Adj[Yb], mb, Yb]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*MatMul[Adj[Yb], Yb, mq]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 12*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 12*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 72*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], ht, Adj[ht], Yt] + 12*MatMul[Adj[hb], Yb]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 12*MatMul[Adj[ht], Yt]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 72*MatMul[Adj[ht], Yt]*
   trace[Adj[Yt], ht, Adj[Yt], Yt] + 12*MatMul[Adj[Yb], hb]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] + 12*MatMul[Adj[Yt], ht]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] + 72*MatMul[Adj[Yt], ht]*
   trace[Adj[Yt], Yt, Adj[ht], Yt] + 12*MatMul[Adj[hb], hb]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*MatMul[Adj[ht], ht]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 24*mh1*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*mh2*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*mh1*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 24*mh2*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*MatMul[mq, Adj[Yb], Yb]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*MatMul[mq, Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*MatMul[Adj[Yb], mb, Yb]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*MatMul[Adj[Yb], Yb, mq]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*MatMul[Adj[Yt], mt, Yt]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*MatMul[Adj[Yt], Yt, mq]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 36*MatMul[Adj[ht], ht]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 108*mh2*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 18*MatMul[mq, Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 36*MatMul[Adj[Yt], mt, Yt]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 18*MatMul[Adj[Yt], Yt, mq]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt] + 72*MatMul[Adj[Yb], Yb]*
   trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 12*MatMul[Adj[Yb], Yb]*
   trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 12*MatMul[Adj[Yt], Yt]*
   trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 24*MatMul[Adj[Yb], Yb]*
   trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 12*MatMul[Adj[Yb], Yb]*
   trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 12*MatMul[Adj[Yt], Yt]*
   trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 72*MatMul[Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 72*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 12*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 12*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 24*MatMul[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 12*MatMul[Adj[Yb], Yb]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 12*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 72*MatMul[Adj[Yt], Yt]*
   trace[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  M1^2*((57511*g1^6)/1125 - (2388*g1^6*Zeta[3])/125) + 
  M1^2*((33*g1^4*g2^2)/25 - (162*g1^4*g2^2*Zeta[3])/25) + 
  M1*M2*((22*g1^4*g2^2)/25 - (108*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((4*g1^4*g2^2)/5 - (54*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((379*g1^2*g2^4)/5 - (486*g1^2*g2^4*Zeta[3])/5) + 
  M1*M2*(50*g1^2*g2^4 - (324*g1^2*g2^4*Zeta[3])/5) + 
  M1^2*((152*g1^2*g2^4)/5 - (162*g1^2*g2^4*Zeta[3])/5) + 
  M2^2*(2157*g2^6 + 3780*g2^6*Zeta[3]) + 
  M1^2*((776*g1^4*g3^2)/75 - (528*g1^4*g3^2*Zeta[3])/25) + 
  M1*M3*((1552*g1^4*g3^2)/225 - (352*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((208*g1^4*g3^2)/45 - (176*g1^4*g3^2*Zeta[3])/25) + 
  M2^2*(664*g2^4*g3^2 - 1296*g2^4*g3^2*Zeta[3]) + 
  M2*M3*(400*g2^4*g3^2 - 864*g2^4*g3^2*Zeta[3]) + 
  M3^2*(272*g2^4*g3^2 - 432*g2^4*g3^2*Zeta[3]) + 
  M3^2*((2464*g1^2*g3^4)/15 - (1056*g1^2*g3^4*Zeta[3])/5) + 
  M1*M3*((4864*g1^2*g3^4)/45 - (704*g1^2*g3^4*Zeta[3])/5) + 
  M1^2*((592*g1^2*g3^4)/9 - (352*g1^2*g3^4*Zeta[3])/5) + 
  M3^2*(192*g2^2*g3^4 - 864*g2^2*g3^4*Zeta[3]) + 
  M2*M3*(64*g2^2*g3^4 - 576*g2^2*g3^4*Zeta[3]) + 
  M2^2*(80*g2^2*g3^4 - 288*g2^2*g3^4*Zeta[3]) + 
  M3^2*((20512*g3^6)/9 + 7680*g3^6*Zeta[3])}
