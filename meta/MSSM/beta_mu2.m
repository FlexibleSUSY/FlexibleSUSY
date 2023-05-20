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

{(-32*g1^2*M1^2)/15 - (32*g3^2*M3^2)/3 + 4*MatMul[ht, Adj[ht]] + 
  4*mh2*MatMul[Yt, Adj[Yt]] + 2*MatMul[mt, Yt, Adj[Yt]] + 
  4*MatMul[Yt, mq, Adj[Yt]] + 2*MatMul[Yt, Adj[Yt], mt], 
 (3424*g1^4*M1^2)/75 + (512*g1^2*g3^2*M1^2)/45 + (512*g1^2*g3^2*M1*M3)/45 + 
  (512*g1^2*g3^2*M3^2)/45 - (128*g3^4*M3^2)/3 + (16*g1^4*mh1)/25 + 
  (16*g1^4*mh2)/25 - (4*g1^2*MatMul[ht, Adj[ht]])/5 + 
  12*g2^2*MatMul[ht, Adj[ht]] + (4*g1^2*M1*MatMul[ht, Adj[Yt]])/5 - 
  12*g2^2*M2*MatMul[ht, Adj[Yt]] + (4*g1^2*M1*MatMul[Yt, Adj[ht]])/5 - 
  12*g2^2*M2*MatMul[Yt, Adj[ht]] - (8*g1^2*M1^2*MatMul[Yt, Adj[Yt]])/5 + 
  24*g2^2*M2^2*MatMul[Yt, Adj[Yt]] - (4*g1^2*mh2*MatMul[Yt, Adj[Yt]])/5 + 
  12*g2^2*mh2*MatMul[Yt, Adj[Yt]] - 
  4*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]] - 
  4*MatMul[MatMul[ht, Adj[Yt], Yt], Adj[ht]] - 
  4*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]] - 
  4*MatMul[MatMul[Yt, Adj[Yt], ht], Adj[ht]] - 
  (2*g1^2*MatMul[mt, Yt, Adj[Yt]])/5 + 6*g2^2*MatMul[mt, Yt, Adj[Yt]] - 
  (4*g1^2*MatMul[Yt, mq, Adj[Yt]])/5 + 12*g2^2*MatMul[Yt, mq, Adj[Yt]] - 
  (2*g1^2*MatMul[Yt, Adj[Yt], mt])/5 + 6*g2^2*MatMul[Yt, Adj[Yt], mt] - 
  4*MatMul[ht, Adj[hb], Yb, Adj[Yt]] - 4*MatMul[ht, Adj[ht], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[hb], hb, Adj[Yt]] - 4*MatMul[Yt, Adj[ht], ht, Adj[Yt]] - 
  4*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] - 
  4*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] - 
  8*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]] - 
  2*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]] - 
  2*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt] - 
  4*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt] + (32*g1^4*trace[mb])/75 + 
  (16*g3^4*trace[mb])/3 + (32*g1^4*trace[me])/25 + (16*g1^4*trace[ml])/25 + 
  (16*g1^4*trace[mq])/75 + (32*g3^4*trace[mq])/3 + (128*g1^4*trace[mt])/75 + 
  (16*g3^4*trace[mt])/3 - 12*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[ht]] - 
  12*MatMul[Yt, Adj[ht]]*trace[ht, Adj[Yt]] - 12*MatMul[ht, Adj[Yt]]*
   trace[Adj[ht], Yt] - 12*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt] - 
  24*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  6*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  6*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, Adj[Yt]]*trace[Yt, Adj[Yt], mt] - 
  12*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, mq], 
 (-3424*g1^6*mh1)/375 - (256*g1^4*g3^2*mh1)/75 - (3424*g1^6*mh2)/375 - 
  (256*g1^4*g3^2*mh2)/75 - (799*g1^4*MatMul[ht, Adj[ht]])/25 - 
  (134*g1^2*g2^2*MatMul[ht, Adj[ht]])/5 - 87*g2^4*MatMul[ht, Adj[ht]] - 
  (16*g1^2*g3^2*MatMul[ht, Adj[ht]])/15 - 176*g2^2*g3^2*MatMul[ht, Adj[ht]] + 
  (32*g3^4*MatMul[ht, Adj[ht]])/3 - 
  (494*g1^4*MatMul[ht, MatMul[Adj[ht], Zeta[3]]])/75 + 
  (156*g1^2*g2^2*MatMul[ht, MatMul[Adj[ht], Zeta[3]]])/5 - 
  18*g2^4*MatMul[ht, MatMul[Adj[ht], Zeta[3]]] - 
  (448*g1^2*g3^2*MatMul[ht, MatMul[Adj[ht], Zeta[3]]])/15 + 
  192*g2^2*g3^2*MatMul[ht, MatMul[Adj[ht], Zeta[3]]] - 
  (1088*g3^4*MatMul[ht, MatMul[Adj[ht], Zeta[3]]])/3 + 
  (12*g1^4*mh1*MatMul[Yt, Adj[Yt]])/25 - 12*g2^4*mh1*MatMul[Yt, Adj[Yt]] - 
  (787*g1^4*mh2*MatMul[Yt, Adj[Yt]])/25 - 
  (134*g1^2*g2^2*mh2*MatMul[Yt, Adj[Yt]])/5 - 
  99*g2^4*mh2*MatMul[Yt, Adj[Yt]] - (16*g1^2*g3^2*mh2*MatMul[Yt, Adj[Yt]])/
   15 - 176*g2^2*g3^2*mh2*MatMul[Yt, Adj[Yt]] + 
  (32*g3^4*mh2*MatMul[Yt, Adj[Yt]])/3 + 
  M1*((1598*g1^4*MatMul[Yt, Adj[ht]])/25 + 
    (988*g1^4*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]])/75) + 
  M1*((134*g1^2*g2^2*MatMul[Yt, Adj[ht]])/5 - 
    (156*g1^2*g2^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]])/5) + 
  M2*((134*g1^2*g2^2*MatMul[Yt, Adj[ht]])/5 - 
    (156*g1^2*g2^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]])/5) + 
  M2*(174*g2^4*MatMul[Yt, Adj[ht]] + 
    36*g2^4*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]) + 
  M1*((16*g1^2*g3^2*MatMul[Yt, Adj[ht]])/15 + 
    (448*g1^2*g3^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]])/15) + 
  M3*((16*g1^2*g3^2*MatMul[Yt, Adj[ht]])/15 + 
    (448*g1^2*g3^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]])/15) + 
  M2*(176*g2^2*g3^2*MatMul[Yt, Adj[ht]] - 192*g2^2*g3^2*
     MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]) + 
  M3*(176*g2^2*g3^2*MatMul[Yt, Adj[ht]] - 192*g2^2*g3^2*
     MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]) + 
  M3*((-64*g3^4*MatMul[Yt, Adj[ht]])/3 + 
    (2176*g3^4*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]])/3) + 
  M1*((1598*g1^4*MatMul[ht, Adj[Yt]])/25 + 
    (988*g1^4*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]])/75) + 
  M1*((134*g1^2*g2^2*MatMul[ht, Adj[Yt]])/5 - 
    (156*g1^2*g2^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]])/5) + 
  M2*((134*g1^2*g2^2*MatMul[ht, Adj[Yt]])/5 - 
    (156*g1^2*g2^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]])/5) + 
  M2*(174*g2^4*MatMul[ht, Adj[Yt]] + 36*g2^4*MatMul[MatMul[ht, Adj[Yt]], 
      Zeta[3]]) + M1*((16*g1^2*g3^2*MatMul[ht, Adj[Yt]])/15 + 
    (448*g1^2*g3^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]])/15) + 
  M3*((16*g1^2*g3^2*MatMul[ht, Adj[Yt]])/15 + 
    (448*g1^2*g3^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]])/15) + 
  M2*(176*g2^2*g3^2*MatMul[ht, Adj[Yt]] - 192*g2^2*g3^2*
     MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]) + 
  M3*(176*g2^2*g3^2*MatMul[ht, Adj[Yt]] - 192*g2^2*g3^2*
     MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]) + 
  M3*((-64*g3^4*MatMul[ht, Adj[Yt]])/3 + 
    (2176*g3^4*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]])/3) - 
  (494*g1^4*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/75 + 
  (156*g1^2*g2^2*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/5 - 
  18*g2^4*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]] - 
  (448*g1^2*g3^2*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/15 + 
  192*g2^2*g3^2*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]] - 
  (1088*g3^4*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/3 + 
  M1^2*((-4794*g1^4*MatMul[Yt, Adj[Yt]])/25 - 
    (988*g1^4*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/25) + 
  M1^2*((-268*g1^2*g2^2*MatMul[Yt, Adj[Yt]])/5 + 
    (312*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/5) + 
  M1*M2*((-268*g1^2*g2^2*MatMul[Yt, Adj[Yt]])/5 + 
    (312*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/5) + 
  M2^2*((-268*g1^2*g2^2*MatMul[Yt, Adj[Yt]])/5 + 
    (312*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/5) + 
  M2^2*(-474*g2^4*MatMul[Yt, Adj[Yt]] - 108*g2^4*MatMul[MatMul[Yt, Adj[Yt]], 
      Zeta[3]]) + M1^2*((-32*g1^2*g3^2*MatMul[Yt, Adj[Yt]])/15 - 
    (896*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/15) + 
  M1*M3*((-32*g1^2*g3^2*MatMul[Yt, Adj[Yt]])/15 - 
    (896*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/15) + 
  M3^2*((-32*g1^2*g3^2*MatMul[Yt, Adj[Yt]])/15 - 
    (896*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]])/15) + 
  M2^2*(-352*g2^2*g3^2*MatMul[Yt, Adj[Yt]] + 
    384*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]) + 
  M2*M3*(-352*g2^2*g3^2*MatMul[Yt, Adj[Yt]] + 
    384*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]) + 
  M3^2*(-352*g2^2*g3^2*MatMul[Yt, Adj[Yt]] + 
    384*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]) + 
  M3^2*(64*g3^4*MatMul[Yt, Adj[Yt]] - 2176*g3^4*MatMul[MatMul[Yt, Adj[Yt]], 
      Zeta[3]]) + (38*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]])/15 + 
  18*g2^2*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]] + 
  (128*g3^2*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]])/3 + 
  (36*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb], MatMul[Adj[ht], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[ht, Adj[Yb], Yb], MatMul[Adj[ht], Zeta[3]]] - 
  (2*g1^2*MatMul[MatMul[ht, Adj[Yt], Yt], Adj[ht]])/3 + 
  18*g2^2*MatMul[MatMul[ht, Adj[Yt], Yt], Adj[ht]] + 
  (128*g3^2*MatMul[MatMul[ht, Adj[Yt], Yt], Adj[ht]])/3 + 
  12*g1^2*MatMul[MatMul[ht, Adj[Yt], Yt], MatMul[Adj[ht], Zeta[3]]] - 
  36*g2^2*MatMul[MatMul[ht, Adj[Yt], Yt], MatMul[Adj[ht], Zeta[3]]] - 
  (247*g1^4*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]])/75 + 
  (78*g1^2*g2^2*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]] - 
  (224*g1^2*g3^2*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]])/15 + 
  96*g2^2*g3^2*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]] - 
  (544*g3^4*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]])/3 - 
  (494*g1^4*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]])/75 + 
  (156*g1^2*g2^2*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]])/5 - 
  18*g2^4*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]] - 
  (448*g1^2*g3^2*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]])/15 + 
  192*g2^2*g3^2*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]] - 
  (1088*g3^4*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]])/3 + 
  (38*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]])/15 + 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]] + 
  (128*g3^2*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]])/3 + 
  (36*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb], MatMul[Adj[ht], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yb], hb], MatMul[Adj[ht], Zeta[3]]] - 
  (128*g3^2*M3*MatMul[MatMul[Yt, Adj[Yb], Yb], Adj[ht]])/3 + 
  M1*((-38*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Adj[ht]])/15 - 
    (36*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb], MatMul[Adj[ht], Zeta[3]]])/5) + 
  M2*(-18*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb], Adj[ht]] + 
    36*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb], MatMul[Adj[ht], Zeta[3]]]) - 
  (2*g1^2*MatMul[MatMul[Yt, Adj[Yt], ht], Adj[ht]])/3 + 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yt], ht], Adj[ht]] + 
  (128*g3^2*MatMul[MatMul[Yt, Adj[Yt], ht], Adj[ht]])/3 + 
  12*g1^2*MatMul[MatMul[Yt, Adj[Yt], ht], MatMul[Adj[ht], Zeta[3]]] - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yt], ht], MatMul[Adj[ht], Zeta[3]]] - 
  (247*g1^4*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]])/75 + 
  (78*g1^2*g2^2*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]] - 
  (224*g1^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]])/15 + 
  96*g2^2*g3^2*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]] - 
  (544*g3^4*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]])/3 - 
  (128*g3^2*M3*MatMul[MatMul[Yt, Adj[Yt], Yt], Adj[ht]])/3 + 
  M1*((2*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Adj[ht]])/3 - 
    12*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt], MatMul[Adj[ht], Zeta[3]]]) + 
  M2*(-18*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt], Adj[ht]] + 
    36*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt], MatMul[Adj[ht], Zeta[3]]]) + 
  (36*g1^2*MatMul[MatMul[ht, Adj[hb], Yb, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[ht, Adj[hb], Yb, Adj[Yt]], Zeta[3]] + 
  12*g1^2*MatMul[MatMul[ht, Adj[ht], Yt, Adj[Yt]], Zeta[3]] - 
  36*g2^2*MatMul[MatMul[ht, Adj[ht], Yt, Adj[Yt]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yt, Adj[hb], hb, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yt, Adj[hb], hb, Adj[Yt]], Zeta[3]] + 
  12*g1^2*MatMul[MatMul[Yt, Adj[ht], ht, Adj[Yt]], Zeta[3]] - 
  36*g2^2*MatMul[MatMul[Yt, Adj[ht], ht, Adj[Yt]], Zeta[3]] + 
  (36*g1^2*mh1*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*mh1*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] + 
  (36*g1^2*mh2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*mh2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] + 
  24*g1^2*mh2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] - 
  72*g2^2*mh2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  12*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yb], Yb], Adj[ht]] - 
  4*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yt], Yt], Adj[ht]] - 
  4*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yb], Yb], Adj[ht]] + 
  12*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt], Adj[ht]] + 
  24*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt], Yt], MatMul[Adj[ht], Zeta[3]]] + 
  (18*g1^2*MatMul[MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] + 
  6*g1^2*MatMul[MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] - 
  18*g2^2*MatMul[MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] + 
  12*g1^2*MatMul[MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] - 
  36*g2^2*MatMul[MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yb], Yb], Adj[ht]] - 
  4*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yt], Yt], Adj[ht]] + 
  (36*g1^2*MatMul[MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]], Zeta[3]] + 
  (36*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yb], hb], Adj[ht]] - 
  4*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], ht], Adj[ht]] + 
  (18*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt], Zeta[3]] - 
  4*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[Yb], Yb], Adj[ht]] + 
  12*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt], Adj[ht]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[Yt], Yt], MatMul[Adj[ht], Zeta[3]]] + 
  12*g1^2*MatMul[MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]], Zeta[3]] - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]], Zeta[3]] + 
  12*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]], Zeta[3]] - 
  36*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]], Zeta[3]] - 
  4*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yb], hb], Adj[ht]] + 
  12*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht], Adj[ht]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], ht], MatMul[Adj[ht], Zeta[3]]] + 
  6*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt], Zeta[3]] - 
  18*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt], Zeta[3]] + 
  24*MatMul[MatMul[ht, Adj[ht], Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[ht], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[ht], ht, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[ht], Yt, Adj[Yt], ht, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[ht], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[ht], ht, Adj[Yt]], Zeta[3]] + 
  72*mh2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  12*MatMul[MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt, Yt, Adj[Yt]], Zeta[3]] + 
  24*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, mq, Adj[Yt]], Zeta[3]] + 
  12*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], mt], Zeta[3]] - 
  (799*g1^4*MatMul[mt, Yt, Adj[Yt]])/50 - 
  (67*g1^2*g2^2*MatMul[mt, Yt, Adj[Yt]])/5 - 
  (87*g2^4*MatMul[mt, Yt, Adj[Yt]])/2 - (8*g1^2*g3^2*MatMul[mt, Yt, Adj[Yt]])/
   15 - 88*g2^2*g3^2*MatMul[mt, Yt, Adj[Yt]] + 
  (16*g3^4*MatMul[mt, Yt, Adj[Yt]])/3 - (799*g1^4*MatMul[Yt, mq, Adj[Yt]])/
   25 - (134*g1^2*g2^2*MatMul[Yt, mq, Adj[Yt]])/5 - 
  87*g2^4*MatMul[Yt, mq, Adj[Yt]] - (16*g1^2*g3^2*MatMul[Yt, mq, Adj[Yt]])/
   15 - 176*g2^2*g3^2*MatMul[Yt, mq, Adj[Yt]] + 
  (32*g3^4*MatMul[Yt, mq, Adj[Yt]])/3 - (799*g1^4*MatMul[Yt, Adj[Yt], mt])/
   50 - (67*g1^2*g2^2*MatMul[Yt, Adj[Yt], mt])/5 - 
  (87*g2^4*MatMul[Yt, Adj[Yt], mt])/2 - (8*g1^2*g3^2*MatMul[Yt, Adj[Yt], mt])/
   15 - 88*g2^2*g3^2*MatMul[Yt, Adj[Yt], mt] + 
  (16*g3^4*MatMul[Yt, Adj[Yt], mt])/3 + 
  (38*g1^2*MatMul[ht, Adj[hb], Yb, Adj[Yt]])/15 + 
  18*g2^2*MatMul[ht, Adj[hb], Yb, Adj[Yt]] + 
  (128*g3^2*MatMul[ht, Adj[hb], Yb, Adj[Yt]])/3 - 
  (2*g1^2*MatMul[ht, Adj[ht], Yt, Adj[Yt]])/3 + 
  18*g2^2*MatMul[ht, Adj[ht], Yt, Adj[Yt]] + 
  (128*g3^2*MatMul[ht, Adj[ht], Yt, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[ht, Adj[Yb], Yb, Adj[Yt]])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 - 
    (38*g1^2*MatMul[ht, Adj[Yb], Yb, Adj[Yt]])/15) + 
  M2*(36*g2^2*MatMul[MatMul[ht, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] - 
    18*g2^2*MatMul[ht, Adj[Yb], Yb, Adj[Yt]]) - 
  (128*g3^2*M3*MatMul[ht, Adj[Yt], Yt, Adj[Yt]])/3 + 
  M1*(-12*g1^2*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
    (2*g1^2*MatMul[ht, Adj[Yt], Yt, Adj[Yt]])/3) + 
  M2*(36*g2^2*MatMul[MatMul[ht, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] - 
    18*g2^2*MatMul[ht, Adj[Yt], Yt, Adj[Yt]]) + 
  (38*g1^2*MatMul[Yt, Adj[hb], hb, Adj[Yt]])/15 + 
  18*g2^2*MatMul[Yt, Adj[hb], hb, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, Adj[hb], hb, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[Yt, Adj[hb], Yb, Adj[Yt]])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[Yt, Adj[hb], Yb, Adj[Yt]], Zeta[3]])/5 - 
    (38*g1^2*MatMul[Yt, Adj[hb], Yb, Adj[Yt]])/15) + 
  M2*(36*g2^2*MatMul[MatMul[Yt, Adj[hb], Yb, Adj[Yt]], Zeta[3]] - 
    18*g2^2*MatMul[Yt, Adj[hb], Yb, Adj[Yt]]) - 
  (2*g1^2*MatMul[Yt, Adj[ht], ht, Adj[Yt]])/3 + 
  18*g2^2*MatMul[Yt, Adj[ht], ht, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, Adj[ht], ht, Adj[Yt]])/3 - 
  (128*g3^2*M3*MatMul[Yt, Adj[ht], Yt, Adj[Yt]])/3 + 
  M1*(-12*g1^2*MatMul[MatMul[Yt, Adj[ht], Yt, Adj[Yt]], Zeta[3]] + 
    (2*g1^2*MatMul[Yt, Adj[ht], Yt, Adj[Yt]])/3) + 
  M2*(36*g2^2*MatMul[MatMul[Yt, Adj[ht], Yt, Adj[Yt]], Zeta[3]] - 
    18*g2^2*MatMul[Yt, Adj[ht], Yt, Adj[Yt]]) - 
  (128*g3^2*M3*MatMul[Yt, Adj[Yb], hb, Adj[Yt]])/3 + 
  M1*((-36*g1^2*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yt]], Zeta[3]])/5 - 
    (38*g1^2*MatMul[Yt, Adj[Yb], hb, Adj[Yt]])/15) + 
  M2*(36*g2^2*MatMul[MatMul[Yt, Adj[Yb], hb, Adj[Yt]], Zeta[3]] - 
    18*g2^2*MatMul[Yt, Adj[Yb], hb, Adj[Yt]]) + 
  (256*g3^2*M3^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/3 + 
  (38*g1^2*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/15 + 
  18*g2^2*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] + 
  (128*g3^2*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/3 + 
  (38*g1^2*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/15 + 
  18*g2^2*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]] + 
  (128*g3^2*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/3 + 
  M1^2*((72*g1^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]])/5 + 
    (76*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]])/15) + 
  M2^2*(-72*g2^2*MatMul[MatMul[Yt, Adj[Yb], Yb, Adj[Yt]], Zeta[3]] + 
    36*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]) - 
  (128*g3^2*M3*MatMul[Yt, Adj[Yt], ht, Adj[Yt]])/3 + 
  M1*(-12*g1^2*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[Yt]], Zeta[3]] + 
    (2*g1^2*MatMul[Yt, Adj[Yt], ht, Adj[Yt]])/3) + 
  M2*(36*g2^2*MatMul[MatMul[Yt, Adj[Yt], ht, Adj[Yt]], Zeta[3]] - 
    18*g2^2*MatMul[Yt, Adj[Yt], ht, Adj[Yt]]) + 
  (256*g3^2*M3^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 - 
  (4*g1^2*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  36*g2^2*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]] + 
  (256*g3^2*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  M1^2*(24*g1^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] - 
    (4*g1^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]])/3) + 
  M2^2*(-72*g2^2*MatMul[MatMul[Yt, Adj[Yt], Yt, Adj[Yt]], Zeta[3]] + 
    36*g2^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]) + 
  (19*g1^2*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]])/15 + 
  9*g2^2*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]] + 
  (64*g3^2*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]])/3 - 
  (g1^2*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  9*g2^2*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]] + 
  (64*g3^2*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]])/3 + 
  (38*g1^2*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]])/15 + 
  18*g2^2*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]])/3 - 
  (2*g1^2*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]])/3 + 
  18*g2^2*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]])/3 + 
  (38*g1^2*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]])/15 + 
  18*g2^2*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]])/3 + 
  (38*g1^2*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]])/15 + 
  18*g2^2*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]])/3 + 
  (19*g1^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt])/15 + 
  9*g2^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  (64*g3^2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt])/3 - 
  (2*g1^2*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]])/3 + 
  18*g2^2*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]])/3 - 
  (2*g1^2*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]])/3 + 
  18*g2^2*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]] + 
  (128*g3^2*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]])/3 - 
  (g1^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt])/3 + 
  9*g2^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
  (64*g3^2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt])/3 + 
  12*MatMul[ht, Adj[hb], Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[ht, Adj[hb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[ht, Adj[ht], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[ht, Adj[ht], Yt, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[ht, Adj[Yb], Yb, Adj[hb], Yb, Adj[Yt]] - 
  4*MatMul[ht, Adj[Yb], Yb, Adj[ht], Yt, Adj[Yt]] - 
  4*MatMul[ht, Adj[Yt], Yt, Adj[hb], Yb, Adj[Yt]] + 
  12*MatMul[ht, Adj[Yt], Yt, Adj[ht], Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[hb], hb, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[hb], hb, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[hb], Yb, Adj[Yb], hb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[hb], Yb, Adj[Yt], ht, Adj[Yt]] - 
  4*MatMul[Yt, Adj[ht], ht, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[ht], ht, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[ht], Yt, Adj[Yb], hb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[ht], Yt, Adj[Yt], ht, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], hb, Adj[hb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], hb, Adj[ht], Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[hb], hb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[ht], ht, Adj[Yt]] + 
  24*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] + 
  12*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  8*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], ht, Adj[hb], Yb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yt], ht, Adj[ht], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], Yt, Adj[hb], hb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yt], Yt, Adj[ht], ht, Adj[Yt]] - 
  4*mh1*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] - 
  8*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  36*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]] + 
  6*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] - 
  2*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  2*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  6*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yb], Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt], Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], mb, Yb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, mq, Adj[Yt]] + 
  6*MatMul[Yt, Adj[Yb], Yb, Adj[Yb], Yb, Adj[Yt], mt] - 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt, Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, mq, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], Yt, Adj[Yt], mt] - 
  4*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yb], Yb, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt], Yt, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], mb, Yb, Adj[Yt]] - 
  4*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, mq, Adj[Yt]] - 
  2*MatMul[Yt, Adj[Yt], Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  12*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt, Yt, Adj[Yt]] + 
  12*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, mq, Adj[Yt]] + 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], Yt, Adj[Yt], mt] - 
  (6848*g1^6*trace[mb])/1125 - (512*g1^4*g3^2*trace[mb])/225 - 
  (256*g1^2*g3^4*trace[mb])/45 + (320*g3^6*trace[mb])/9 + 
  (8*g1^4*MatMul[Yt, Adj[Yt]]*trace[mb])/25 - (6848*g1^6*trace[me])/375 - 
  (512*g1^4*g3^2*trace[me])/75 + (24*g1^4*MatMul[Yt, Adj[Yt]]*trace[me])/25 - 
  (3424*g1^6*trace[ml])/375 - (256*g1^4*g3^2*trace[ml])/75 + 
  (12*g1^4*MatMul[Yt, Adj[Yt]]*trace[ml])/25 - 12*g2^4*MatMul[Yt, Adj[Yt]]*
   trace[ml] - (3424*g1^6*trace[mq])/1125 - (256*g1^4*g3^2*trace[mq])/225 - 
  (512*g1^2*g3^4*trace[mq])/45 + (640*g3^6*trace[mq])/9 + 
  (4*g1^4*MatMul[Yt, Adj[Yt]]*trace[mq])/25 - 36*g2^4*MatMul[Yt, Adj[Yt]]*
   trace[mq] - (27392*g1^6*trace[mt])/1125 - (2048*g1^4*g3^2*trace[mt])/225 - 
  (256*g1^2*g3^4*trace[mt])/45 + (320*g3^6*trace[mt])/9 + 
  (32*g1^4*MatMul[Yt, Adj[Yt]]*trace[mt])/25 - (224*g1^4*trace[hb, Adj[hb]])/
   25 - 64*g3^4*trace[hb, Adj[hb]] + 24*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*
   trace[hb, Adj[hb]] + 24*MatMul[MatMul[Yt, Adj[Yb], Yb], Adj[ht]]*
   trace[hb, Adj[Yb]] + 24*MatMul[Yt, Adj[hb], Yb, Adj[Yt]]*
   trace[hb, Adj[Yb]] - (288*g1^4*trace[he, Adj[he]])/25 + 
  8*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[he, Adj[he]] + 
  8*MatMul[MatMul[Yt, Adj[Yb], Yb], Adj[ht]]*trace[he, Adj[Ye]] + 
  8*MatMul[Yt, Adj[hb], Yb, Adj[Yt]]*trace[he, Adj[Ye]] - 
  (416*g1^4*trace[ht, Adj[ht]])/25 - 64*g3^4*trace[ht, Adj[ht]] + 
  14*g1^2*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[ht]] + 
  54*g2^2*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[ht]] - 
  32*g3^2*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[ht]] + 
  (84*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[ht, Adj[ht]])/5 - 
  108*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[ht, Adj[ht]] + 
  192*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[ht, Adj[ht]] - 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[ht, Adj[ht]] + 
  12*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*trace[ht, Adj[ht]] + 
  14*g1^2*MatMul[Yt, Adj[ht]]*trace[ht, Adj[Yt]] + 
  54*g2^2*MatMul[Yt, Adj[ht]]*trace[ht, Adj[Yt]] - 
  32*g3^2*MatMul[Yt, Adj[ht]]*trace[ht, Adj[Yt]] + 
  (84*g1^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]*trace[ht, Adj[Yt]])/5 - 
  108*g2^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]*trace[ht, Adj[Yt]] + 
  192*g3^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]*trace[ht, Adj[Yt]] - 
  12*MatMul[MatMul[Yt, Adj[Yb], Yb], Adj[ht]]*trace[ht, Adj[Yt]] + 
  12*MatMul[MatMul[Yt, Adj[Yt], Yt], Adj[ht]]*trace[ht, Adj[Yt]] - 
  12*MatMul[Yt, Adj[hb], Yb, Adj[Yt]]*trace[ht, Adj[Yt]] + 
  12*MatMul[Yt, Adj[ht], Yt, Adj[Yt]]*trace[ht, Adj[Yt]] + 
  24*MatMul[ht, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[hb], Yb] + 
  24*MatMul[Yt, Adj[Yb], hb, Adj[Yt]]*trace[Adj[hb], Yb] + 
  8*MatMul[ht, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[he], Ye] + 
  8*MatMul[Yt, Adj[Yb], hb, Adj[Yt]]*trace[Adj[he], Ye] + 
  14*g1^2*MatMul[ht, Adj[Yt]]*trace[Adj[ht], Yt] + 
  54*g2^2*MatMul[ht, Adj[Yt]]*trace[Adj[ht], Yt] - 
  32*g3^2*MatMul[ht, Adj[Yt]]*trace[Adj[ht], Yt] + 
  (84*g1^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]*trace[Adj[ht], Yt])/5 - 
  108*g2^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]*trace[Adj[ht], Yt] + 
  192*g3^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]*trace[Adj[ht], Yt] - 
  12*MatMul[ht, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[ht], Yt] + 
  12*MatMul[ht, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[ht], Yt] - 
  12*MatMul[Yt, Adj[Yb], hb, Adj[Yt]]*trace[Adj[ht], Yt] + 
  12*MatMul[Yt, Adj[Yt], ht, Adj[Yt]]*trace[Adj[ht], Yt] - 
  72*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[Yt]]*trace[Adj[ht], Yt] + 
  M1*((224*g1^4*trace[hb, Adj[Yb]])/15 + (96*g1^4*trace[he, Adj[Ye]])/5 + 
    (416*g1^4*trace[ht, Adj[Yt]])/15 + (224*g1^4*trace[Adj[hb], Yb])/15 + 
    (96*g1^4*trace[Adj[he], Ye])/5 + (416*g1^4*trace[Adj[ht], Yt])/15) + 
  M3*((320*g3^4*trace[hb, Adj[Yb]])/3 + (320*g3^4*trace[ht, Adj[Yt]])/3 + 
    (320*g3^4*trace[Adj[hb], Yb])/3 + (320*g3^4*trace[Adj[ht], Yt])/3) + 
  M1*(-14*g1^2*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[Yt]] - 
    (84*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[ht, Adj[Yt]])/5 - 
    14*g1^2*MatMul[Yt, Adj[Yt]]*trace[Adj[ht], Yt] - 
    (84*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[ht], Yt])/5) + 
  M2*(-54*g2^2*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[Yt]] + 
    108*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[ht, Adj[Yt]] - 
    54*g2^2*MatMul[Yt, Adj[Yt]]*trace[Adj[ht], Yt] + 
    108*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[ht], Yt]) + 
  M3*(32*g3^2*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[Yt]] - 
    192*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[ht, Adj[Yt]] + 
    32*g3^2*MatMul[Yt, Adj[Yt]]*trace[Adj[ht], Yt] - 
    192*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[ht], Yt]) - 
  (224*g1^4*mh1*trace[Adj[Yb], Yb])/25 - 64*g3^4*mh1*trace[Adj[Yb], Yb] + 
  24*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]]*trace[Adj[Yb], Yb] + 
  24*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]]*trace[Adj[Yb], Yb] + 
  24*MatMul[ht, Adj[hb], Yb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[hb], hb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  48*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  24*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  12*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  24*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt]*trace[Adj[Yb], Yb] - 
  (288*g1^4*mh1*trace[Adj[Ye], Ye])/25 + 
  8*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]]*trace[Adj[Ye], Ye] + 
  8*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]]*trace[Adj[Ye], Ye] + 
  8*MatMul[ht, Adj[hb], Yb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  8*MatMul[Yt, Adj[hb], hb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  16*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  8*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  4*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  8*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  8*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  8*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt]*trace[Adj[Ye], Ye] - 
  (416*g1^4*mh2*trace[Adj[Yt], Yt])/25 - 64*g3^4*mh2*trace[Adj[Yt], Yt] + 
  14*g1^2*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt] + 
  54*g2^2*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt] - 
  32*g3^2*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt] + 
  (84*g1^2*MatMul[ht, MatMul[Adj[ht], Zeta[3]]]*trace[Adj[Yt], Yt])/5 - 
  108*g2^2*MatMul[ht, MatMul[Adj[ht], Zeta[3]]]*trace[Adj[Yt], Yt] + 
  192*g3^2*MatMul[ht, MatMul[Adj[ht], Zeta[3]]]*trace[Adj[Yt], Yt] + 
  28*g1^2*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  108*g2^2*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  64*g3^2*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  (168*g1^2*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  216*g2^2*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] + 
  384*g3^2*mh2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] - 
  12*MatMul[MatMul[ht, Adj[Yb], Yb], Adj[ht]]*trace[Adj[Yt], Yt] + 
  12*MatMul[MatMul[ht, Adj[Yt], Yt], Adj[ht]]*trace[Adj[Yt], Yt] + 
  (42*g1^2*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  54*g2^2*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] + 
  96*g3^2*MatMul[MatMul[mt, Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] + 
  (84*g1^2*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  108*g2^2*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] + 
  192*g3^2*MatMul[MatMul[Yt, mq, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt] - 
  12*MatMul[MatMul[Yt, Adj[Yb], hb], Adj[ht]]*trace[Adj[Yt], Yt] + 
  12*MatMul[MatMul[Yt, Adj[Yt], ht], Adj[ht]]*trace[Adj[Yt], Yt] + 
  (42*g1^2*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]]*trace[Adj[Yt], Yt])/5 - 
  54*g2^2*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]]*trace[Adj[Yt], Yt] + 
  96*g3^2*MatMul[MatMul[Yt, Adj[Yt], mt], Zeta[3]]*trace[Adj[Yt], Yt] + 
  7*g1^2*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  27*g2^2*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  16*g3^2*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  14*g1^2*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  54*g2^2*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  32*g3^2*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  7*g1^2*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt] + 
  27*g2^2*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt] - 
  16*g3^2*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt] - 
  12*MatMul[ht, Adj[hb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[ht, Adj[ht], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, Adj[hb], hb, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yt, Adj[ht], ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  12*mh1*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  24*mh2*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  36*mh2*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  6*MatMul[mt, Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  6*MatMul[mt, Yt, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, mq, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yt, mq, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, Adj[Yb], mb, Yb, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  12*MatMul[Yt, Adj[Yb], Yb, mq, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  6*MatMul[Yt, Adj[Yb], Yb, Adj[Yt], mt]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yt, Adj[Yt], mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  12*MatMul[Yt, Adj[Yt], Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt] + 
  6*MatMul[Yt, Adj[Yt], Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt] - 
  72*MatMul[Yt, Adj[Yt]]*trace[ht, Adj[ht]]*trace[Adj[Yt], Yt] - 
  72*MatMul[Yt, Adj[ht]]*trace[ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
  72*MatMul[ht, Adj[Yt]]*trace[Adj[ht], Yt]*trace[Adj[Yt], Yt] - 
  36*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt]^2 - 
  108*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt]^2 - 
  18*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt]^2 - 
  36*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt]^2 - 
  18*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt]^2 + 
  M1^2*((-224*g1^4*trace[Adj[Yb], Yb])/5 - (288*g1^4*trace[Adj[Ye], Ye])/5 - 
    (416*g1^4*trace[Adj[Yt], Yt])/5) + 
  M3^2*(-320*g3^4*trace[Adj[Yb], Yb] - 320*g3^4*trace[Adj[Yt], Yt]) + 
  M1*(-14*g1^2*MatMul[Yt, Adj[ht]]*trace[Adj[Yt], Yt] - 
    (84*g1^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]*trace[Adj[Yt], Yt])/5) + 
  M2*(-54*g2^2*MatMul[Yt, Adj[ht]]*trace[Adj[Yt], Yt] + 
    108*g2^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]*trace[Adj[Yt], Yt]) + 
  M3*(32*g3^2*MatMul[Yt, Adj[ht]]*trace[Adj[Yt], Yt] - 
    192*g3^2*MatMul[Yt, MatMul[Adj[ht], Zeta[3]]]*trace[Adj[Yt], Yt]) + 
  M1*(-14*g1^2*MatMul[ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
    (84*g1^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M2*(-54*g2^2*MatMul[ht, Adj[Yt]]*trace[Adj[Yt], Yt] + 
    108*g2^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M3*(32*g3^2*MatMul[ht, Adj[Yt]]*trace[Adj[Yt], Yt] - 
    192*g3^2*MatMul[MatMul[ht, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M1^2*(28*g1^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
    (168*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt])/5) + 
  M2^2*(108*g2^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] - 
    216*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt]) + 
  M3^2*(-64*g3^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt] + 
    384*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt]) - 
  (224*g1^4*trace[Yb, Adj[Yb], mb])/25 - 64*g3^4*trace[Yb, Adj[Yb], mb] + 
  24*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Yb, Adj[Yb], mb] - 
  (288*g1^4*trace[Ye, Adj[Ye], me])/25 + 8*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*
   trace[Ye, Adj[Ye], me] - (416*g1^4*trace[Yt, Adj[Yt], mt])/25 - 
  64*g3^4*trace[Yt, Adj[Yt], mt] + 14*g1^2*MatMul[Yt, Adj[Yt]]*
   trace[Yt, Adj[Yt], mt] + 54*g2^2*MatMul[Yt, Adj[Yt]]*
   trace[Yt, Adj[Yt], mt] - 32*g3^2*MatMul[Yt, Adj[Yt]]*
   trace[Yt, Adj[Yt], mt] + (84*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*
    trace[Yt, Adj[Yt], mt])/5 - 108*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*
   trace[Yt, Adj[Yt], mt] + 192*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*
   trace[Yt, Adj[Yt], mt] - 12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*
   trace[Yt, Adj[Yt], mt] + 12*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*
   trace[Yt, Adj[Yt], mt] - 72*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt]*
   trace[Yt, Adj[Yt], mt] - (224*g1^4*trace[Adj[Yb], Yb, mq])/25 - 
  64*g3^4*trace[Adj[Yb], Yb, mq] + 24*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*
   trace[Adj[Yb], Yb, mq] - (288*g1^4*trace[Adj[Ye], Ye, ml])/25 + 
  8*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Ye], Ye, ml] - 
  (416*g1^4*trace[Adj[Yt], Yt, mq])/25 - 64*g3^4*trace[Adj[Yt], Yt, mq] + 
  14*g1^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, mq] + 
  54*g2^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, mq] - 
  32*g3^2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, mq] + 
  (84*g1^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt, mq])/5 - 
  108*g2^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt, mq] + 
  192*g3^2*MatMul[MatMul[Yt, Adj[Yt]], Zeta[3]]*trace[Adj[Yt], Yt, mq] - 
  12*MatMul[Yt, Adj[Yb], Yb, Adj[Yt]]*trace[Adj[Yt], Yt, mq] + 
  12*MatMul[Yt, Adj[Yt], Yt, Adj[Yt]]*trace[Adj[Yt], Yt, mq] - 
  72*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt]*trace[Adj[Yt], Yt, mq] + 
  24*MatMul[Yt, Adj[Yt]]*trace[hb, Adj[ht], Yt, Adj[Yb]] + 
  24*MatMul[Yt, Adj[ht]]*trace[hb, Adj[Yt], Yt, Adj[Yb]] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Adj[hb], hb, Adj[Yt], Yt] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Adj[ht], ht, Adj[Yb], Yb] + 
  144*MatMul[Yt, Adj[Yt]]*trace[Adj[ht], ht, Adj[Yt], Yt] + 
  24*MatMul[ht, Adj[Yt]]*trace[Adj[ht], Yt, Adj[Yb], Yb] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], ht, Adj[hb], Yb] + 
  144*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], ht, Adj[ht], Yt] + 
  24*MatMul[Yt, Adj[ht]]*trace[Adj[Yt], ht, Adj[Yb], Yb] + 
  144*MatMul[Yt, Adj[ht]]*trace[Adj[Yt], ht, Adj[Yt], Yt] + 
  24*MatMul[ht, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[hb], Yb] + 
  144*MatMul[ht, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[ht], Yt] + 
  24*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*mh1*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  48*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  24*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  12*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt, Adj[Yb], Yb] + 
  72*MatMul[ht, Adj[ht]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  216*mh2*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  36*MatMul[mt, Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  72*MatMul[Yt, mq, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  36*MatMul[Yt, Adj[Yt], mt]*trace[Adj[Yt], Yt, Adj[Yt], Yt] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 
  144*MatMul[Yt, Adj[Yt]]*trace[Yt, Adj[Yt], Yt, Adj[Yt], mt] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 
  24*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  144*MatMul[Yt, Adj[Yt]]*trace[Adj[Yt], Yt, Adj[Yt], Yt, mq] + 
  M1^2*((864496*g1^6)/1125 - (38208*g1^6*Zeta[3])/125) + 
  M1^2*((432*g1^4*g2^2)/5 - (2592*g1^4*g2^2*Zeta[3])/25) + 
  M1*M2*((288*g1^4*g2^2)/5 - (1728*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((864*g1^4*g2^2)/25 - (864*g1^4*g2^2*Zeta[3])/25) + 
  M1^2*((8576*g1^4*g3^2)/75 - (8448*g1^4*g3^2*Zeta[3])/25) + 
  M1*M3*((17152*g1^4*g3^2)/225 - (5632*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((512*g1^4*g3^2)/9 - (2816*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((-176*g1^2*g3^4)/15 - (1056*g1^2*g3^4*Zeta[3])/5) + 
  M1*M3*((-1376*g1^2*g3^4)/45 - (704*g1^2*g3^4*Zeta[3])/5) + 
  M1^2*((-32*g1^2*g3^4)/9 - (352*g1^2*g3^4*Zeta[3])/5) + 
  M3^2*(720*g2^2*g3^4 - 864*g2^2*g3^4*Zeta[3]) + 
  M2*M3*(480*g2^2*g3^4 - 576*g2^2*g3^4*Zeta[3]) + 
  M2^2*(288*g2^2*g3^4 - 288*g2^2*g3^4*Zeta[3]) + 
  M3^2*((20512*g3^6)/9 + 7680*g3^6*Zeta[3])}
