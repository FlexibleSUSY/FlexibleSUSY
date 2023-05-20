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

{(-8*g1^2*M1^2)/15 - (32*g3^2*M3^2)/3 + 4*MatMul[hb, Adj[hb]] + 
  4*mh1*MatMul[Yb, Adj[Yb]] + 2*MatMul[mb, Yb, Adj[Yb]] + 
  4*MatMul[Yb, mq, Adj[Yb]] + 2*MatMul[Yb, Adj[Yb], mb], 
 (808*g1^4*M1^2)/75 + (128*g1^2*g3^2*M1^2)/45 + (128*g1^2*g3^2*M1*M3)/45 + 
  (128*g1^2*g3^2*M3^2)/45 - (128*g3^4*M3^2)/3 + (4*g1^4*mh1)/25 + 
  (4*g1^4*mh2)/25 + (4*g1^2*MatMul[hb, Adj[hb]])/5 + 
  12*g2^2*MatMul[hb, Adj[hb]] - (4*g1^2*M1*MatMul[hb, Adj[Yb]])/5 - 
  12*g2^2*M2*MatMul[hb, Adj[Yb]] - (4*g1^2*M1*MatMul[Yb, Adj[hb]])/5 - 
  12*g2^2*M2*MatMul[Yb, Adj[hb]] + (8*g1^2*M1^2*MatMul[Yb, Adj[Yb]])/5 + 
  24*g2^2*M2^2*MatMul[Yb, Adj[Yb]] + (4*g1^2*mh1*MatMul[Yb, Adj[Yb]])/5 + 
  12*g2^2*mh1*MatMul[Yb, Adj[Yb]] - 
  4*MatMul[MatMul[hb, Adj[Yb], Yb], Adj[hb]] - 
  4*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]] - 
  4*MatMul[MatMul[Yb, Adj[Yb], hb], Adj[hb]] - 
  4*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]] + 
  (2*g1^2*MatMul[mb, Yb, Adj[Yb]])/5 + 6*g2^2*MatMul[mb, Yb, Adj[Yb]] + 
  (4*g1^2*MatMul[Yb, mq, Adj[Yb]])/5 + 12*g2^2*MatMul[Yb, mq, Adj[Yb]] + 
  (2*g1^2*MatMul[Yb, Adj[Yb], mb])/5 + 6*g2^2*MatMul[Yb, Adj[Yb], mb] - 
  4*MatMul[hb, Adj[hb], Yb, Adj[Yb]] - 4*MatMul[hb, Adj[ht], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[hb], hb, Adj[Yb]] - 4*MatMul[Yb, Adj[ht], ht, Adj[Yb]] - 
  8*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]] - 
  4*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] - 
  2*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]] - 
  2*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]] - 
  4*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb] - 
  4*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb] + (8*g1^4*trace[mb])/75 + 
  (16*g3^4*trace[mb])/3 + (8*g1^4*trace[me])/25 + (4*g1^4*trace[ml])/25 + 
  (4*g1^4*trace[mq])/75 + (32*g3^4*trace[mq])/3 + (32*g1^4*trace[mt])/75 + 
  (16*g3^4*trace[mt])/3 - 12*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[hb]] - 
  12*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yb]] - 4*MatMul[Yb, Adj[Yb]]*
   trace[he, Adj[he]] - 4*MatMul[Yb, Adj[hb]]*trace[he, Adj[Ye]] - 
  12*MatMul[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 4*MatMul[hb, Adj[Yb]]*
   trace[Adj[he], Ye] - 12*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  24*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  6*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  6*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb] - 
  4*MatMul[hb, Adj[hb]]*trace[Adj[Ye], Ye] - 8*mh1*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Ye], Ye] - 2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye] - 
  12*MatMul[Yb, Adj[Yb]]*trace[Yb, Adj[Yb], mb] - 
  4*MatMul[Yb, Adj[Yb]]*trace[Ye, Adj[Ye], me] - 
  12*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb, mq] - 
  4*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, ml], 
 (-808*g1^6*mh1)/375 - (64*g1^4*g3^2*mh1)/75 - (808*g1^6*mh2)/375 - 
  (64*g1^4*g3^2*mh2)/75 - (337*g1^4*MatMul[hb, Adj[hb]])/15 - 
  (86*g1^2*g2^2*MatMul[hb, Adj[hb]])/5 - 87*g2^4*MatMul[hb, Adj[hb]] - 
  (48*g1^2*g3^2*MatMul[hb, Adj[hb]])/5 - 176*g2^2*g3^2*MatMul[hb, Adj[hb]] + 
  (32*g3^4*MatMul[hb, Adj[hb]])/3 - 
  (14*g1^4*MatMul[hb, MatMul[Adj[hb], Zeta[3]]])/75 + 
  (84*g1^2*g2^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]])/5 - 
  18*g2^4*MatMul[hb, MatMul[Adj[hb], Zeta[3]]] + 
  (64*g1^2*g3^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]])/3 + 
  192*g2^2*g3^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]] - 
  (1088*g3^4*MatMul[hb, MatMul[Adj[hb], Zeta[3]]])/3 - 
  (1721*g1^4*mh1*MatMul[Yb, Adj[Yb]])/75 - 
  (86*g1^2*g2^2*mh1*MatMul[Yb, Adj[Yb]])/5 - 
  99*g2^4*mh1*MatMul[Yb, Adj[Yb]] - (48*g1^2*g3^2*mh1*MatMul[Yb, Adj[Yb]])/
   5 - 176*g2^2*g3^2*mh1*MatMul[Yb, Adj[Yb]] + 
  (32*g3^4*mh1*MatMul[Yb, Adj[Yb]])/3 - (12*g1^4*mh2*MatMul[Yb, Adj[Yb]])/
   25 - 12*g2^4*mh2*MatMul[Yb, Adj[Yb]] + 
  M1*((674*g1^4*MatMul[Yb, Adj[hb]])/15 + 
    (28*g1^4*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]])/75) + 
  M1*((86*g1^2*g2^2*MatMul[Yb, Adj[hb]])/5 - 
    (84*g1^2*g2^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]])/5) + 
  M2*((86*g1^2*g2^2*MatMul[Yb, Adj[hb]])/5 - 
    (84*g1^2*g2^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]])/5) + 
  M2*(174*g2^4*MatMul[Yb, Adj[hb]] + 
    36*g2^4*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]) + 
  M1*((48*g1^2*g3^2*MatMul[Yb, Adj[hb]])/5 - 
    (64*g1^2*g3^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]])/3) + 
  M3*((48*g1^2*g3^2*MatMul[Yb, Adj[hb]])/5 - 
    (64*g1^2*g3^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]])/3) + 
  M2*(176*g2^2*g3^2*MatMul[Yb, Adj[hb]] - 192*g2^2*g3^2*
     MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]) + 
  M3*(176*g2^2*g3^2*MatMul[Yb, Adj[hb]] - 192*g2^2*g3^2*
     MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]) + 
  M3*((-64*g3^4*MatMul[Yb, Adj[hb]])/3 + 
    (2176*g3^4*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]])/3) + 
  M1*((674*g1^4*MatMul[hb, Adj[Yb]])/15 + 
    (28*g1^4*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]])/75) + 
  M1*((86*g1^2*g2^2*MatMul[hb, Adj[Yb]])/5 - 
    (84*g1^2*g2^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]])/5) + 
  M2*((86*g1^2*g2^2*MatMul[hb, Adj[Yb]])/5 - 
    (84*g1^2*g2^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]])/5) + 
  M2*(174*g2^4*MatMul[hb, Adj[Yb]] + 36*g2^4*MatMul[MatMul[hb, Adj[Yb]], 
      Zeta[3]]) + M1*((48*g1^2*g3^2*MatMul[hb, Adj[Yb]])/5 - 
    (64*g1^2*g3^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]])/3) + 
  M3*((48*g1^2*g3^2*MatMul[hb, Adj[Yb]])/5 - 
    (64*g1^2*g3^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]])/3) + 
  M2*(176*g2^2*g3^2*MatMul[hb, Adj[Yb]] - 192*g2^2*g3^2*
     MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]) + 
  M3*(176*g2^2*g3^2*MatMul[hb, Adj[Yb]] - 192*g2^2*g3^2*
     MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]) + 
  M3*((-64*g3^4*MatMul[hb, Adj[Yb]])/3 + 
    (2176*g3^4*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]])/3) - 
  (14*g1^4*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/75 + 
  (84*g1^2*g2^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/5 - 
  18*g2^4*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]] + 
  (64*g1^2*g3^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3 + 
  192*g2^2*g3^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]] - 
  (1088*g3^4*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3 + 
  M1^2*((-674*g1^4*MatMul[Yb, Adj[Yb]])/5 - 
    (28*g1^4*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/25) + 
  M1^2*((-172*g1^2*g2^2*MatMul[Yb, Adj[Yb]])/5 + 
    (168*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/5) + 
  M1*M2*((-172*g1^2*g2^2*MatMul[Yb, Adj[Yb]])/5 + 
    (168*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/5) + 
  M2^2*((-172*g1^2*g2^2*MatMul[Yb, Adj[Yb]])/5 + 
    (168*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/5) + 
  M2^2*(-474*g2^4*MatMul[Yb, Adj[Yb]] - 108*g2^4*MatMul[MatMul[Yb, Adj[Yb]], 
      Zeta[3]]) + M1^2*((-96*g1^2*g3^2*MatMul[Yb, Adj[Yb]])/5 + 
    (128*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3) + 
  M1*M3*((-96*g1^2*g3^2*MatMul[Yb, Adj[Yb]])/5 + 
    (128*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3) + 
  M3^2*((-96*g1^2*g3^2*MatMul[Yb, Adj[Yb]])/5 + 
    (128*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]])/3) + 
  M2^2*(-352*g2^2*g3^2*MatMul[Yb, Adj[Yb]] + 
    384*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]) + 
  M2*M3*(-352*g2^2*g3^2*MatMul[Yb, Adj[Yb]] + 
    384*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]) + 
  M3^2*(-352*g2^2*g3^2*MatMul[Yb, Adj[Yb]] + 
    384*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]) + 
  M3^2*(64*g3^4*MatMul[Yb, Adj[Yb]] - 2176*g3^4*MatMul[MatMul[Yb, Adj[Yb]], 
      Zeta[3]]) - (2*g1^2*MatMul[MatMul[hb, Adj[Yb], Yb], Adj[hb]])/3 + 
  18*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb], Adj[hb]] + 
  (128*g3^2*MatMul[MatMul[hb, Adj[Yb], Yb], Adj[hb]])/3 + 
  (12*g1^2*MatMul[MatMul[hb, Adj[Yb], Yb], MatMul[Adj[hb], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb], MatMul[Adj[hb], Zeta[3]]] - 
  (58*g1^2*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]])/15 + 
  18*g2^2*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]] + 
  (128*g3^2*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]])/3 + 
  (36*g1^2*MatMul[MatMul[hb, Adj[Yt], Yt], MatMul[Adj[hb], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[hb, Adj[Yt], Yt], MatMul[Adj[hb], Zeta[3]]] - 
  (7*g1^4*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]])/75 + 
  (42*g1^2*g2^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]] + 
  (32*g1^2*g3^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]])/3 + 
  96*g2^2*g3^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]] - 
  (544*g3^4*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]])/3 - 
  (14*g1^4*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]])/75 + 
  (84*g1^2*g2^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]])/5 - 
  18*g2^4*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]] + 
  (64*g1^2*g3^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]])/3 + 
  192*g2^2*g3^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]] - 
  (1088*g3^4*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]])/3 - 
  (2*g1^2*MatMul[MatMul[Yb, Adj[Yb], hb], Adj[hb]])/3 + 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], hb], Adj[hb]] + 
  (128*g3^2*MatMul[MatMul[Yb, Adj[Yb], hb], Adj[hb]])/3 + 
  (12*g1^2*MatMul[MatMul[Yb, Adj[Yb], hb], MatMul[Adj[hb], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yb], hb], MatMul[Adj[hb], Zeta[3]]] - 
  (7*g1^4*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]])/75 + 
  (42*g1^2*g2^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]] + 
  (32*g1^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]])/3 + 
  96*g2^2*g3^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]] - 
  (544*g3^4*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]])/3 - 
  (128*g3^2*M3*MatMul[MatMul[Yb, Adj[Yb], Yb], Adj[hb]])/3 + 
  M1*((2*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Adj[hb]])/3 - 
    (12*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb], MatMul[Adj[hb], Zeta[3]]])/5) + 
  M2*(-18*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], Adj[hb]] + 
    36*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb], MatMul[Adj[hb], Zeta[3]]]) - 
  (58*g1^2*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]])/15 + 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]] + 
  (128*g3^2*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]])/3 + 
  (36*g1^2*MatMul[MatMul[Yb, Adj[Yt], ht], MatMul[Adj[hb], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yt], ht], MatMul[Adj[hb], Zeta[3]]] - 
  (128*g3^2*M3*MatMul[MatMul[Yb, Adj[Yt], Yt], Adj[hb]])/3 + 
  M1*((58*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Adj[hb]])/15 - 
    (36*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt], MatMul[Adj[hb], Zeta[3]]])/5) + 
  M2*(-18*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt], Adj[hb]] + 
    36*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt], MatMul[Adj[hb], Zeta[3]]]) + 
  (12*g1^2*MatMul[MatMul[hb, Adj[hb], Yb, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[hb, Adj[hb], Yb, Adj[Yb]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[hb, Adj[ht], Yt, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[hb, Adj[ht], Yt, Adj[Yb]], Zeta[3]] + 
  (12*g1^2*MatMul[MatMul[Yb, Adj[hb], hb, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[hb], hb, Adj[Yb]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yb, Adj[ht], ht, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[ht], ht, Adj[Yb]], Zeta[3]] + 
  (24*g1^2*mh1*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]])/5 - 
  72*g2^2*mh1*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  (36*g1^2*mh1*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*mh1*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] + 
  (36*g1^2*mh2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*mh2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] + 
  12*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb], Adj[hb]] + 
  24*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb], Yb], MatMul[Adj[hb], Zeta[3]]] - 
  4*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yt], Yt], Adj[hb]] - 
  4*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yb], Yb], Adj[hb]] + 
  12*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yt], Yt], Adj[hb]] + 
  (6*g1^2*MatMul[MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  (18*g1^2*MatMul[MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] + 
  (12*g1^2*MatMul[MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb], Adj[hb]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[Yb], Yb], MatMul[Adj[hb], Zeta[3]]] - 
  4*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[Yt], Yt], Adj[hb]] + 
  (12*g1^2*MatMul[MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]], Zeta[3]] + 
  (12*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb], Adj[hb]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], hb], MatMul[Adj[hb], Zeta[3]]] + 
  (6*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb], Zeta[3]] - 
  4*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yt], ht], Adj[hb]] - 
  4*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yb], Yb], Adj[hb]] + 
  12*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yt], Yt], Adj[hb]] + 
  (36*g1^2*MatMul[MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]], Zeta[3]] - 
  4*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], hb], Adj[hb]] + 
  (18*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yt], ht], Adj[hb]] + 
  24*MatMul[MatMul[hb, Adj[hb], Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[hb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[hb], hb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[hb], Yb, Adj[Yb], hb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[hb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[hb], hb, Adj[Yb]], Zeta[3]] + 
  72*mh1*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  12*MatMul[MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb, Yb, Adj[Yb]], Zeta[3]] + 
  24*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, mq, Adj[Yb]], Zeta[3]] + 
  12*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], mb], Zeta[3]] - 
  (337*g1^4*MatMul[mb, Yb, Adj[Yb]])/30 - 
  (43*g1^2*g2^2*MatMul[mb, Yb, Adj[Yb]])/5 - 
  (87*g2^4*MatMul[mb, Yb, Adj[Yb]])/2 - 
  (24*g1^2*g3^2*MatMul[mb, Yb, Adj[Yb]])/5 - 
  88*g2^2*g3^2*MatMul[mb, Yb, Adj[Yb]] + (16*g3^4*MatMul[mb, Yb, Adj[Yb]])/
   3 - (337*g1^4*MatMul[Yb, mq, Adj[Yb]])/15 - 
  (86*g1^2*g2^2*MatMul[Yb, mq, Adj[Yb]])/5 - 
  87*g2^4*MatMul[Yb, mq, Adj[Yb]] - (48*g1^2*g3^2*MatMul[Yb, mq, Adj[Yb]])/
   5 - 176*g2^2*g3^2*MatMul[Yb, mq, Adj[Yb]] + 
  (32*g3^4*MatMul[Yb, mq, Adj[Yb]])/3 - (337*g1^4*MatMul[Yb, Adj[Yb], mb])/
   30 - (43*g1^2*g2^2*MatMul[Yb, Adj[Yb], mb])/5 - 
  (87*g2^4*MatMul[Yb, Adj[Yb], mb])/2 - 
  (24*g1^2*g3^2*MatMul[Yb, Adj[Yb], mb])/5 - 
  88*g2^2*g3^2*MatMul[Yb, Adj[Yb], mb] + (16*g3^4*MatMul[Yb, Adj[Yb], mb])/
   3 - (2*g1^2*MatMul[hb, Adj[hb], Yb, Adj[Yb]])/3 + 
  18*g2^2*MatMul[hb, Adj[hb], Yb, Adj[Yb]] + 
  (128*g3^2*MatMul[hb, Adj[hb], Yb, Adj[Yb]])/3 - 
  (58*g1^2*MatMul[hb, Adj[ht], Yt, Adj[Yb]])/15 + 
  18*g2^2*MatMul[hb, Adj[ht], Yt, Adj[Yb]] + 
  (128*g3^2*MatMul[hb, Adj[ht], Yt, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[hb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  M1*((-12*g1^2*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]])/5 + 
    (2*g1^2*MatMul[hb, Adj[Yb], Yb, Adj[Yb]])/3) + 
  M2*(36*g2^2*MatMul[MatMul[hb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] - 
    18*g2^2*MatMul[hb, Adj[Yb], Yb, Adj[Yb]]) - 
  (128*g3^2*M3*MatMul[hb, Adj[Yt], Yt, Adj[Yb]])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 + 
    (58*g1^2*MatMul[hb, Adj[Yt], Yt, Adj[Yb]])/15) + 
  M2*(36*g2^2*MatMul[MatMul[hb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] - 
    18*g2^2*MatMul[hb, Adj[Yt], Yt, Adj[Yb]]) - 
  (2*g1^2*MatMul[Yb, Adj[hb], hb, Adj[Yb]])/3 + 
  18*g2^2*MatMul[Yb, Adj[hb], hb, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, Adj[hb], hb, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[Yb, Adj[hb], Yb, Adj[Yb]])/3 + 
  M1*((-12*g1^2*MatMul[MatMul[Yb, Adj[hb], Yb, Adj[Yb]], Zeta[3]])/5 + 
    (2*g1^2*MatMul[Yb, Adj[hb], Yb, Adj[Yb]])/3) + 
  M2*(36*g2^2*MatMul[MatMul[Yb, Adj[hb], Yb, Adj[Yb]], Zeta[3]] - 
    18*g2^2*MatMul[Yb, Adj[hb], Yb, Adj[Yb]]) - 
  (58*g1^2*MatMul[Yb, Adj[ht], ht, Adj[Yb]])/15 + 
  18*g2^2*MatMul[Yb, Adj[ht], ht, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, Adj[ht], ht, Adj[Yb]])/3 - 
  (128*g3^2*M3*MatMul[Yb, Adj[ht], Yt, Adj[Yb]])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[Yb, Adj[ht], Yt, Adj[Yb]], Zeta[3]])/5 + 
    (58*g1^2*MatMul[Yb, Adj[ht], Yt, Adj[Yb]])/15) + 
  M2*(36*g2^2*MatMul[MatMul[Yb, Adj[ht], Yt, Adj[Yb]], Zeta[3]] - 
    18*g2^2*MatMul[Yb, Adj[ht], Yt, Adj[Yb]]) - 
  (128*g3^2*M3*MatMul[Yb, Adj[Yb], hb, Adj[Yb]])/3 + 
  M1*((-12*g1^2*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[Yb]], Zeta[3]])/5 + 
    (2*g1^2*MatMul[Yb, Adj[Yb], hb, Adj[Yb]])/3) + 
  M2*(36*g2^2*MatMul[MatMul[Yb, Adj[Yb], hb, Adj[Yb]], Zeta[3]] - 
    18*g2^2*MatMul[Yb, Adj[Yb], hb, Adj[Yb]]) + 
  (256*g3^2*M3^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 - 
  (4*g1^2*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  36*g2^2*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]] + 
  (256*g3^2*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  M1^2*((24*g1^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]])/5 - 
    (4*g1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]])/3) + 
  M2^2*(-72*g2^2*MatMul[MatMul[Yb, Adj[Yb], Yb, Adj[Yb]], Zeta[3]] + 
    36*g2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]) - 
  (128*g3^2*M3*MatMul[Yb, Adj[Yt], ht, Adj[Yb]])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yb]], Zeta[3]])/5 + 
    (58*g1^2*MatMul[Yb, Adj[Yt], ht, Adj[Yb]])/15) + 
  M2*(36*g2^2*MatMul[MatMul[Yb, Adj[Yt], ht, Adj[Yb]], Zeta[3]] - 
    18*g2^2*MatMul[Yb, Adj[Yt], ht, Adj[Yb]]) + 
  (256*g3^2*M3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/3 - 
  (58*g1^2*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/15 + 
  18*g2^2*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] + 
  (128*g3^2*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/3 - 
  (58*g1^2*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/15 + 
  18*g2^2*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]] + 
  (128*g3^2*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/3 + 
  M1^2*((72*g1^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]])/5 - 
    (116*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]])/15) + 
  M2^2*(-72*g2^2*MatMul[MatMul[Yb, Adj[Yt], Yt, Adj[Yb]], Zeta[3]] + 
    36*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]) - 
  (g1^2*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]])/3 + 
  9*g2^2*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]] + 
  (64*g3^2*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]])/3 - 
  (29*g1^2*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]])/15 + 
  9*g2^2*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]] + 
  (64*g3^2*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]])/3 - 
  (2*g1^2*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]])/3 + 
  18*g2^2*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]])/3 - 
  (58*g1^2*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]])/15 + 
  18*g2^2*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]])/3 - 
  (2*g1^2*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]])/3 + 
  18*g2^2*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]])/3 - 
  (2*g1^2*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]])/3 + 
  18*g2^2*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]])/3 - 
  (g1^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb])/3 + 
  9*g2^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  (64*g3^2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb])/3 - 
  (58*g1^2*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]])/15 + 
  18*g2^2*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]])/3 - 
  (58*g1^2*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]])/15 + 
  18*g2^2*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]] + 
  (128*g3^2*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]])/3 - 
  (29*g1^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb])/15 + 
  9*g2^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  (64*g3^2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb])/3 + 
  12*MatMul[hb, Adj[hb], Yb, Adj[Yb], Yb, Adj[Yb]] - 
  4*MatMul[hb, Adj[hb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[hb, Adj[ht], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[hb, Adj[ht], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[hb, Adj[Yb], Yb, Adj[hb], Yb, Adj[Yb]] - 
  4*MatMul[hb, Adj[Yb], Yb, Adj[ht], Yt, Adj[Yb]] - 
  4*MatMul[hb, Adj[Yt], Yt, Adj[hb], Yb, Adj[Yb]] + 
  12*MatMul[hb, Adj[Yt], Yt, Adj[ht], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[hb], hb, Adj[Yb], Yb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[hb], hb, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[hb], Yb, Adj[Yb], hb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[hb], Yb, Adj[Yt], ht, Adj[Yb]] - 
  4*MatMul[Yb, Adj[ht], ht, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[ht], ht, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[ht], Yt, Adj[Yb], hb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[ht], Yt, Adj[Yt], ht, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yb], hb, Adj[hb], Yb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], hb, Adj[ht], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[hb], hb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[ht], ht, Adj[Yb]] + 
  36*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]] - 
  8*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*mh2*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], ht, Adj[hb], Yb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], ht, Adj[ht], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[hb], hb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[ht], ht, Adj[Yb]] - 
  8*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] - 
  4*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  24*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  6*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]] - 
  2*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  2*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  6*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb]] - 
  4*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb], Yb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb], Yb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yt], Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb, Yb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, mq, Adj[Yb]] + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yb], mb] - 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], mt, Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, mq, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yb], mb] - 
  4*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb], Yb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yt], Yt, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb, Yb, Adj[Yb]] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, mq, Adj[Yb]] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yb], mb] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], mt, Yt, Adj[Yb]] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, mq, Adj[Yb]] + 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yb], mb] - 
  (1616*g1^6*trace[mb])/1125 - (128*g1^4*g3^2*trace[mb])/225 - 
  (64*g1^2*g3^4*trace[mb])/45 + (320*g3^6*trace[mb])/9 - 
  (8*g1^4*MatMul[Yb, Adj[Yb]]*trace[mb])/25 - (1616*g1^6*trace[me])/375 - 
  (128*g1^4*g3^2*trace[me])/75 - (24*g1^4*MatMul[Yb, Adj[Yb]]*trace[me])/25 - 
  (808*g1^6*trace[ml])/375 - (64*g1^4*g3^2*trace[ml])/75 - 
  (12*g1^4*MatMul[Yb, Adj[Yb]]*trace[ml])/25 - 12*g2^4*MatMul[Yb, Adj[Yb]]*
   trace[ml] - (808*g1^6*trace[mq])/1125 - (64*g1^4*g3^2*trace[mq])/225 - 
  (128*g1^2*g3^4*trace[mq])/45 + (640*g3^6*trace[mq])/9 - 
  (4*g1^4*MatMul[Yb, Adj[Yb]]*trace[mq])/25 - 36*g2^4*MatMul[Yb, Adj[Yb]]*
   trace[mq] - (6464*g1^6*trace[mt])/1125 - (512*g1^4*g3^2*trace[mt])/225 - 
  (64*g1^2*g3^4*trace[mt])/45 + (320*g3^6*trace[mt])/9 - 
  (32*g1^4*MatMul[Yb, Adj[Yb]]*trace[mt])/25 - (56*g1^4*trace[hb, Adj[hb]])/
   25 - 64*g3^4*trace[hb, Adj[hb]] + 14*g1^2*MatMul[Yb, Adj[Yb]]*
   trace[hb, Adj[hb]] + 54*g2^2*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[hb]] - 
  32*g3^2*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[hb]] - 
  12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[hb, Adj[hb]] - 
  108*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[hb, Adj[hb]] + 
  192*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[hb, Adj[hb]] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[hb, Adj[hb]] - 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[hb, Adj[hb]] + 
  14*g1^2*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yb]] + 
  54*g2^2*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yb]] - 
  32*g3^2*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yb]] - 
  12*g1^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[hb, Adj[Yb]] - 
  108*g2^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[hb, Adj[Yb]] + 
  192*g3^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[hb, Adj[Yb]] + 
  12*MatMul[MatMul[Yb, Adj[Yb], Yb], Adj[hb]]*trace[hb, Adj[Yb]] - 
  12*MatMul[MatMul[Yb, Adj[Yt], Yt], Adj[hb]]*trace[hb, Adj[Yb]] + 
  12*MatMul[Yb, Adj[hb], Yb, Adj[Yb]]*trace[hb, Adj[Yb]] - 
  12*MatMul[Yb, Adj[ht], Yt, Adj[Yb]]*trace[hb, Adj[Yb]] - 
  (72*g1^4*trace[he, Adj[he]])/25 - 6*g1^2*MatMul[Yb, Adj[Yb]]*
   trace[he, Adj[he]] + 18*g2^2*MatMul[Yb, Adj[Yb]]*trace[he, Adj[he]] + 
  32*g3^2*MatMul[Yb, Adj[Yb]]*trace[he, Adj[he]] + 
  12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[he, Adj[he]] - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[he, Adj[he]] + 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[he, Adj[he]] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[he, Adj[he]] - 
  6*g1^2*MatMul[Yb, Adj[hb]]*trace[he, Adj[Ye]] + 
  18*g2^2*MatMul[Yb, Adj[hb]]*trace[he, Adj[Ye]] + 
  32*g3^2*MatMul[Yb, Adj[hb]]*trace[he, Adj[Ye]] + 
  12*g1^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[he, Adj[Ye]] - 
  36*g2^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[he, Adj[Ye]] + 
  4*MatMul[MatMul[Yb, Adj[Yb], Yb], Adj[hb]]*trace[he, Adj[Ye]] - 
  4*MatMul[MatMul[Yb, Adj[Yt], Yt], Adj[hb]]*trace[he, Adj[Ye]] + 
  4*MatMul[Yb, Adj[hb], Yb, Adj[Yb]]*trace[he, Adj[Ye]] - 
  4*MatMul[Yb, Adj[ht], Yt, Adj[Yb]]*trace[he, Adj[Ye]] - 
  (104*g1^4*trace[ht, Adj[ht]])/25 - 64*g3^4*trace[ht, Adj[ht]] + 
  24*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[ht, Adj[ht]] + 
  24*MatMul[MatMul[Yb, Adj[Yt], Yt], Adj[hb]]*trace[ht, Adj[Yt]] + 
  24*MatMul[Yb, Adj[ht], Yt, Adj[Yb]]*trace[ht, Adj[Yt]] + 
  14*g1^2*MatMul[hb, Adj[Yb]]*trace[Adj[hb], Yb] + 
  54*g2^2*MatMul[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  32*g3^2*MatMul[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  12*g1^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[hb], Yb] - 
  108*g2^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[hb], Yb] + 
  192*g3^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[hb], Yb] + 
  12*MatMul[hb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  12*MatMul[hb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[hb], Yb] + 
  12*MatMul[Yb, Adj[Yb], hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  12*MatMul[Yb, Adj[Yt], ht, Adj[Yb]]*trace[Adj[hb], Yb] - 
  72*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  24*MatMul[Yb, Adj[Yb]]*trace[he, Adj[Ye]]*trace[Adj[hb], Yb] - 
  6*g1^2*MatMul[hb, Adj[Yb]]*trace[Adj[he], Ye] + 
  18*g2^2*MatMul[hb, Adj[Yb]]*trace[Adj[he], Ye] + 
  32*g3^2*MatMul[hb, Adj[Yb]]*trace[Adj[he], Ye] + 
  12*g1^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[he], Ye] - 
  36*g2^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[he], Ye] + 
  4*MatMul[hb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[he], Ye] - 
  4*MatMul[hb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[he], Ye] + 
  4*MatMul[Yb, Adj[Yb], hb, Adj[Yb]]*trace[Adj[he], Ye] - 
  4*MatMul[Yb, Adj[Yt], ht, Adj[Yb]]*trace[Adj[he], Ye] - 
  24*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[Yb]]*trace[Adj[he], Ye] - 
  8*MatMul[Yb, Adj[Yb]]*trace[he, Adj[Ye]]*trace[Adj[he], Ye] + 
  M3*(32*g3^2*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[Yb]] - 
    192*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[hb, Adj[Yb]] - 
    32*g3^2*MatMul[Yb, Adj[Yb]]*trace[he, Adj[Ye]] + 
    32*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[hb], Yb] - 
    192*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[hb], Yb] - 
    32*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[he], Ye]) + 
  M1*(-14*g1^2*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[Yb]] + 
    12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[hb, Adj[Yb]] + 
    6*g1^2*MatMul[Yb, Adj[Yb]]*trace[he, Adj[Ye]] - 
    12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[he, Adj[Ye]] - 
    14*g1^2*MatMul[Yb, Adj[Yb]]*trace[Adj[hb], Yb] + 
    12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[hb], Yb] + 
    6*g1^2*MatMul[Yb, Adj[Yb]]*trace[Adj[he], Ye] - 
    12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[he], Ye]) + 
  M2*(-54*g2^2*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[Yb]] + 
    108*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[hb, Adj[Yb]] - 
    18*g2^2*MatMul[Yb, Adj[Yb]]*trace[he, Adj[Ye]] + 
    36*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[he, Adj[Ye]] - 
    54*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[hb], Yb] + 
    108*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[hb], Yb] - 
    18*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[he], Ye] + 
    36*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[he], Ye]) + 
  24*MatMul[hb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[ht], Yt] + 
  24*MatMul[Yb, Adj[Yt], ht, Adj[Yb]]*trace[Adj[ht], Yt] + 
  M1*((56*g1^4*trace[hb, Adj[Yb]])/15 + (24*g1^4*trace[he, Adj[Ye]])/5 + 
    (104*g1^4*trace[ht, Adj[Yt]])/15 + (56*g1^4*trace[Adj[hb], Yb])/15 + 
    (24*g1^4*trace[Adj[he], Ye])/5 + (104*g1^4*trace[Adj[ht], Yt])/15) + 
  M3*((320*g3^4*trace[hb, Adj[Yb]])/3 + (320*g3^4*trace[ht, Adj[Yt]])/3 + 
    (320*g3^4*trace[Adj[hb], Yb])/3 + (320*g3^4*trace[Adj[ht], Yt])/3) - 
  (56*g1^4*mh1*trace[Adj[Yb], Yb])/25 - 64*g3^4*mh1*trace[Adj[Yb], Yb] + 
  14*g1^2*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb] + 
  54*g2^2*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  32*g3^2*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  12*g1^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Yb], Yb] - 
  108*g2^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Yb], Yb] + 
  192*g3^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Yb], Yb] + 
  28*g1^2*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  108*g2^2*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  64*g3^2*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  24*g1^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
  216*g2^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  384*g3^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  12*MatMul[MatMul[hb, Adj[Yb], Yb], Adj[hb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]]*trace[Adj[Yb], Yb] - 
  6*g1^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
  54*g2^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
  12*g1^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
  108*g2^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  192*g3^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  12*MatMul[MatMul[Yb, Adj[Yb], hb], Adj[hb]]*trace[Adj[Yb], Yb] - 
  6*g1^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  54*g2^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]]*trace[Adj[Yb], Yb] - 
  12*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]]*trace[Adj[Yb], Yb] + 
  7*g1^2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  27*g2^2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  16*g3^2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  14*g1^2*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  54*g2^2*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  32*g3^2*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  7*g1^2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb] + 
  27*g2^2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb] - 
  16*g3^2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb] + 
  12*MatMul[hb, Adj[hb], Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[hb, Adj[ht], Yt, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[hb], hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yb, Adj[ht], ht, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  36*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  24*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  6*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  6*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb] + 
  6*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  6*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb]*trace[Adj[Yb], Yb] - 
  72*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  72*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  24*MatMul[Yb, Adj[Yb]]*trace[he, Adj[he]]*trace[Adj[Yb], Yb] - 
  24*MatMul[Yb, Adj[hb]]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  72*MatMul[hb, Adj[Yb]]*trace[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 
  24*MatMul[hb, Adj[Yb]]*trace[Adj[he], Ye]*trace[Adj[Yb], Yb] - 
  36*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb]^2 - 
  108*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb]^2 - 
  36*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb]^2 - 
  (72*g1^4*mh1*trace[Adj[Ye], Ye])/25 - 6*g1^2*MatMul[hb, Adj[hb]]*
   trace[Adj[Ye], Ye] + 18*g2^2*MatMul[hb, Adj[hb]]*trace[Adj[Ye], Ye] + 
  32*g3^2*MatMul[hb, Adj[hb]]*trace[Adj[Ye], Ye] + 
  12*g1^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Ye], Ye] - 
  36*g2^2*MatMul[hb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Ye], Ye] - 
  12*g1^2*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  36*g2^2*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  64*g3^2*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  24*g1^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] - 
  72*g2^2*mh1*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  4*MatMul[MatMul[hb, Adj[Yb], Yb], Adj[hb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]]*trace[Adj[Ye], Ye] + 
  6*g1^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] - 
  18*g2^2*MatMul[MatMul[mb, Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  12*g1^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] - 
  36*g2^2*MatMul[MatMul[Yb, mq, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  4*MatMul[MatMul[Yb, Adj[Yb], hb], Adj[hb]]*trace[Adj[Ye], Ye] + 
  6*g1^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]]*trace[Adj[Ye], Ye] - 
  18*g2^2*MatMul[MatMul[Yb, Adj[Yb], mb], Zeta[3]]*trace[Adj[Ye], Ye] - 
  4*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]]*trace[Adj[Ye], Ye] - 
  3*g1^2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  9*g2^2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  16*g3^2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  6*g1^2*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  18*g2^2*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  32*g3^2*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  3*g1^2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye] + 
  9*g2^2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye] + 
  16*g3^2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye] + 
  4*MatMul[hb, Adj[hb], Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[hb, Adj[ht], Yt, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yb, Adj[hb], hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yb, Adj[ht], ht, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  12*mh1*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  8*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  2*MatMul[mb, Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  2*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yb, mq, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yb, Adj[Yb], mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yb, Adj[Yb], Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye] + 
  2*MatMul[Yb, Adj[Yb], Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb]*trace[Adj[Ye], Ye] - 
  24*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[hb]]*trace[Adj[Ye], Ye] - 
  24*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  8*MatMul[Yb, Adj[Yb]]*trace[he, Adj[he]]*trace[Adj[Ye], Ye] - 
  8*MatMul[Yb, Adj[hb]]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  24*MatMul[hb, Adj[Yb]]*trace[Adj[hb], Yb]*trace[Adj[Ye], Ye] - 
  8*MatMul[hb, Adj[Yb]]*trace[Adj[he], Ye]*trace[Adj[Ye], Ye] - 
  24*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  72*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  24*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  4*MatMul[hb, Adj[hb]]*trace[Adj[Ye], Ye]^2 - 12*mh1*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Ye], Ye]^2 - 2*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye]^2 - 
  4*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye]^2 - 
  2*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye]^2 + 
  M3*(32*g3^2*MatMul[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    192*g3^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
    32*g3^2*MatMul[hb, Adj[Yb]]*trace[Adj[Ye], Ye]) + 
  M3*(32*g3^2*MatMul[Yb, Adj[hb]]*trace[Adj[Yb], Yb] - 
    192*g3^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Yb], Yb] - 
    32*g3^2*MatMul[Yb, Adj[hb]]*trace[Adj[Ye], Ye]) + 
  M3^2*(-64*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
    384*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
    64*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye]) + 
  M1*(-14*g1^2*MatMul[Yb, Adj[hb]]*trace[Adj[Yb], Yb] + 
    12*g1^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Yb], Yb] + 
    6*g1^2*MatMul[Yb, Adj[hb]]*trace[Adj[Ye], Ye] - 
    12*g1^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Ye], Ye]) + 
  M2*(-54*g2^2*MatMul[Yb, Adj[hb]]*trace[Adj[Yb], Yb] + 
    108*g2^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Yb], Yb] - 
    18*g2^2*MatMul[Yb, Adj[hb]]*trace[Adj[Ye], Ye] + 
    36*g2^2*MatMul[Yb, MatMul[Adj[hb], Zeta[3]]]*trace[Adj[Ye], Ye]) + 
  M1*(-14*g1^2*MatMul[hb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
    12*g1^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
    6*g1^2*MatMul[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
    12*g1^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye]) + 
  M2*(-54*g2^2*MatMul[hb, Adj[Yb]]*trace[Adj[Yb], Yb] + 
    108*g2^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
    18*g2^2*MatMul[hb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
    36*g2^2*MatMul[MatMul[hb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye]) + 
  M1^2*(28*g1^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    24*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] - 
    12*g1^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] + 
    24*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye]) + 
  M2^2*(108*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
    216*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Yb], Yb] + 
    36*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
    72*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye]) - 
  (104*g1^4*mh2*trace[Adj[Yt], Yt])/25 - 64*g3^4*mh2*trace[Adj[Yt], Yt] + 
  24*MatMul[MatMul[hb, Adj[Yt], Yt], Adj[hb]]*trace[Adj[Yt], Yt] + 
  24*MatMul[MatMul[Yb, Adj[Yt], ht], Adj[hb]]*trace[Adj[Yt], Yt] + 
  24*MatMul[hb, Adj[ht], Yt, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  24*MatMul[Yb, Adj[ht], ht, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  24*mh1*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  48*mh2*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  12*MatMul[mb, Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  24*MatMul[Yb, mq, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  24*MatMul[Yb, Adj[Yt], mt, Yt, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  24*MatMul[Yb, Adj[Yt], Yt, mq, Adj[Yb]]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb], mb]*trace[Adj[Yt], Yt] + 
  M1^2*((-56*g1^4*trace[Adj[Yb], Yb])/5 - (72*g1^4*trace[Adj[Ye], Ye])/5 - 
    (104*g1^4*trace[Adj[Yt], Yt])/5) + 
  M3^2*(-320*g3^4*trace[Adj[Yb], Yb] - 320*g3^4*trace[Adj[Yt], Yt]) - 
  (56*g1^4*trace[Yb, Adj[Yb], mb])/25 - 64*g3^4*trace[Yb, Adj[Yb], mb] + 
  14*g1^2*MatMul[Yb, Adj[Yb]]*trace[Yb, Adj[Yb], mb] + 
  54*g2^2*MatMul[Yb, Adj[Yb]]*trace[Yb, Adj[Yb], mb] - 
  32*g3^2*MatMul[Yb, Adj[Yb]]*trace[Yb, Adj[Yb], mb] - 
  12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Yb, Adj[Yb], mb] - 
  108*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Yb, Adj[Yb], mb] + 
  192*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Yb, Adj[Yb], mb] + 
  12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Yb, Adj[Yb], mb] - 
  12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Yb, Adj[Yb], mb] - 
  72*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
  (72*g1^4*trace[Ye, Adj[Ye], me])/25 - 6*g1^2*MatMul[Yb, Adj[Yb]]*
   trace[Ye, Adj[Ye], me] + 18*g2^2*MatMul[Yb, Adj[Yb]]*
   trace[Ye, Adj[Ye], me] + 32*g3^2*MatMul[Yb, Adj[Yb]]*
   trace[Ye, Adj[Ye], me] + 12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*
   trace[Ye, Adj[Ye], me] - 36*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*
   trace[Ye, Adj[Ye], me] + 4*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*
   trace[Ye, Adj[Ye], me] - 4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   trace[Ye, Adj[Ye], me] - 24*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*
   trace[Ye, Adj[Ye], me] - 8*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye]*
   trace[Ye, Adj[Ye], me] - (104*g1^4*trace[Yt, Adj[Yt], mt])/25 - 
  64*g3^4*trace[Yt, Adj[Yt], mt] + 24*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   trace[Yt, Adj[Yt], mt] - (56*g1^4*trace[Adj[Yb], Yb, mq])/25 - 
  64*g3^4*trace[Adj[Yb], Yb, mq] + 14*g1^2*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Yb], Yb, mq] + 54*g2^2*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Yb], Yb, mq] - 32*g3^2*MatMul[Yb, Adj[Yb]]*
   trace[Adj[Yb], Yb, mq] - 12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*
   trace[Adj[Yb], Yb, mq] - 108*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*
   trace[Adj[Yb], Yb, mq] + 192*g3^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*
   trace[Adj[Yb], Yb, mq] + 12*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*
   trace[Adj[Yb], Yb, mq] - 12*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*
   trace[Adj[Yb], Yb, mq] - 72*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, mq] - 24*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, mq] - (72*g1^4*trace[Adj[Ye], Ye, ml])/25 - 
  6*g1^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, ml] + 
  18*g2^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, ml] + 
  32*g3^2*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, ml] + 
  12*g1^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye, ml] - 
  36*g2^2*MatMul[MatMul[Yb, Adj[Yb]], Zeta[3]]*trace[Adj[Ye], Ye, ml] + 
  4*MatMul[Yb, Adj[Yb], Yb, Adj[Yb]]*trace[Adj[Ye], Ye, ml] - 
  4*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Ye], Ye, ml] - 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
  8*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] - 
  (104*g1^4*trace[Adj[Yt], Yt, mq])/25 - 64*g3^4*trace[Adj[Yt], Yt, mq] + 
  24*MatMul[Yb, Adj[Yt], Yt, Adj[Yb]]*trace[Adj[Yt], Yt, mq] + 
  24*MatMul[Yb, Adj[Yb]]*trace[hb, Adj[ht], Yt, Adj[Yb]] + 
  24*MatMul[Yb, Adj[hb]]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  144*MatMul[Yb, Adj[Yb]]*trace[Adj[hb], hb, Adj[Yb], Yb] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[hb], hb, Adj[Yt], Yt] + 
  48*MatMul[Yb, Adj[Yb]]*trace[Adj[he], he, Adj[Ye], Ye] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[ht], ht, Adj[Yb], Yb] + 
  24*MatMul[hb, Adj[Yb]]*trace[Adj[ht], Yt, Adj[Yb], Yb] + 
  144*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], hb, Adj[hb], Yb] + 
  144*MatMul[Yb, Adj[hb]]*trace[Adj[Yb], hb, Adj[Yb], Yb] + 
  144*MatMul[hb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[hb], Yb] + 
  72*MatMul[hb, Adj[hb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  216*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  72*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  36*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yb], Yb, Adj[Yb], Yb] + 
  48*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], he, Adj[he], Ye] + 
  48*MatMul[Yb, Adj[hb]]*trace[Adj[Ye], he, Adj[Ye], Ye] + 
  48*MatMul[hb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[he], Ye] + 
  24*MatMul[hb, Adj[hb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  72*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  24*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  12*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Ye], Ye, Adj[Ye], Ye] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[Yt], ht, Adj[hb], Yb] + 
  24*MatMul[Yb, Adj[hb]]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  24*MatMul[hb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[hb], Yb] + 
  24*MatMul[hb, Adj[hb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  48*mh1*MatMul[Yb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*mh2*MatMul[Yb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[mb, Yb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*MatMul[Yb, mq, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[Yb, Adj[Yb], mb]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  144*MatMul[Yb, Adj[Yb]]*trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  48*MatMul[Yb, Adj[Yb]]*trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  144*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  48*MatMul[Yb, Adj[Yb]]*trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 
  24*MatMul[Yb, Adj[Yb]]*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  M1^2*((227548*g1^6)/1125 - (9552*g1^6*Zeta[3])/125) + 
  M1^2*((108*g1^4*g2^2)/5 - (648*g1^4*g2^2*Zeta[3])/25) + 
  M1*M2*((72*g1^4*g2^2)/5 - (432*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((216*g1^4*g2^2)/25 - (216*g1^4*g2^2*Zeta[3])/25) + 
  M1^2*((2912*g1^4*g3^2)/75 - (2112*g1^4*g3^2*Zeta[3])/25) + 
  M1*M3*((5824*g1^4*g3^2)/225 - (1408*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((3968*g1^4*g3^2)/225 - (704*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((1936*g1^2*g3^4)/15 - (1056*g1^2*g3^4*Zeta[3])/5) + 
  M1*M3*((3616*g1^2*g3^4)/45 - (704*g1^2*g3^4*Zeta[3])/5) + 
  M1^2*((2336*g1^2*g3^4)/45 - (352*g1^2*g3^4*Zeta[3])/5) + 
  M3^2*(720*g2^2*g3^4 - 864*g2^2*g3^4*Zeta[3]) + 
  M2*M3*(480*g2^2*g3^4 - 576*g2^2*g3^4*Zeta[3]) + 
  M2^2*(288*g2^2*g3^4 - 288*g2^2*g3^4*Zeta[3]) + 
  M3^2*((20512*g3^6)/9 + 7680*g3^6*Zeta[3])}
