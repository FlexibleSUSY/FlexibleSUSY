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

{(-24*g1^2*M1^2)/5 + 4*MatMul[he, Adj[he]] + 4*mh1*MatMul[Ye, Adj[Ye]] + 
  2*MatMul[me, Ye, Adj[Ye]] + 4*MatMul[Ye, ml, Adj[Ye]] + 
  2*MatMul[Ye, Adj[Ye], me], (2808*g1^4*M1^2)/25 + (36*g1^4*mh1)/25 + 
  (36*g1^4*mh2)/25 - (12*g1^2*MatMul[he, Adj[he]])/5 + 
  12*g2^2*MatMul[he, Adj[he]] + (12*g1^2*M1*MatMul[he, Adj[Ye]])/5 - 
  12*g2^2*M2*MatMul[he, Adj[Ye]] + (12*g1^2*M1*MatMul[Ye, Adj[he]])/5 - 
  12*g2^2*M2*MatMul[Ye, Adj[he]] - (24*g1^2*M1^2*MatMul[Ye, Adj[Ye]])/5 + 
  24*g2^2*M2^2*MatMul[Ye, Adj[Ye]] - (12*g1^2*mh1*MatMul[Ye, Adj[Ye]])/5 + 
  12*g2^2*mh1*MatMul[Ye, Adj[Ye]] - 
  4*MatMul[MatMul[he, Adj[Ye], Ye], Adj[he]] - 
  4*MatMul[MatMul[Ye, Adj[Ye], he], Adj[he]] - 
  (6*g1^2*MatMul[me, Ye, Adj[Ye]])/5 + 6*g2^2*MatMul[me, Ye, Adj[Ye]] - 
  (12*g1^2*MatMul[Ye, ml, Adj[Ye]])/5 + 12*g2^2*MatMul[Ye, ml, Adj[Ye]] - 
  (6*g1^2*MatMul[Ye, Adj[Ye], me])/5 + 6*g2^2*MatMul[Ye, Adj[Ye], me] - 
  4*MatMul[he, Adj[he], Ye, Adj[Ye]] - 4*MatMul[Ye, Adj[he], he, Adj[Ye]] - 
  8*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]] - 
  2*MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]] - 
  4*MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]] - 
  4*MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]] - 
  4*MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]] - 
  2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me] + (24*g1^4*trace[mb])/25 + 
  (72*g1^4*trace[me])/25 + (36*g1^4*trace[ml])/25 + (12*g1^4*trace[mq])/25 + 
  (96*g1^4*trace[mt])/25 - 12*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[hb]] - 
  12*MatMul[Ye, Adj[he]]*trace[hb, Adj[Yb]] - 4*MatMul[Ye, Adj[Ye]]*
   trace[he, Adj[he]] - 4*MatMul[Ye, Adj[he]]*trace[he, Adj[Ye]] - 
  12*MatMul[he, Adj[Ye]]*trace[Adj[hb], Yb] - 4*MatMul[he, Adj[Ye]]*
   trace[Adj[he], Ye] - 12*MatMul[he, Adj[he]]*trace[Adj[Yb], Yb] - 
  24*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  6*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  6*MatMul[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb] - 
  4*MatMul[he, Adj[he]]*trace[Adj[Ye], Ye] - 8*mh1*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Ye], Ye] - 2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Ye], Ye] - 
  12*MatMul[Ye, Adj[Ye]]*trace[Yb, Adj[Yb], mb] - 
  4*MatMul[Ye, Adj[Ye]]*trace[Ye, Adj[Ye], me] - 
  12*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb, mq] - 
  4*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye, ml], 
 (-2808*g1^6*mh1)/125 - (2808*g1^6*mh2)/125 - (1503*g1^4*MatMul[he, Adj[he]])/
   25 - 54*g1^2*g2^2*MatMul[he, Adj[he]] - 87*g2^4*MatMul[he, Adj[he]] - 
  (162*g1^4*MatMul[he, MatMul[Adj[he], Zeta[3]]])/5 + 
  (324*g1^2*g2^2*MatMul[he, MatMul[Adj[he], Zeta[3]]])/5 - 
  18*g2^4*MatMul[he, MatMul[Adj[he], Zeta[3]]] - 
  (1467*g1^4*mh1*MatMul[Ye, Adj[Ye]])/25 - 
  54*g1^2*g2^2*mh1*MatMul[Ye, Adj[Ye]] - 99*g2^4*mh1*MatMul[Ye, Adj[Ye]] + 
  (36*g1^4*mh2*MatMul[Ye, Adj[Ye]])/25 - 12*g2^4*mh2*MatMul[Ye, Adj[Ye]] + 
  M1*((3006*g1^4*MatMul[Ye, Adj[he]])/25 + 
    (324*g1^4*MatMul[Ye, MatMul[Adj[he], Zeta[3]]])/5) + 
  M1*(54*g1^2*g2^2*MatMul[Ye, Adj[he]] - 
    (324*g1^2*g2^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]])/5) + 
  M2*(54*g1^2*g2^2*MatMul[Ye, Adj[he]] - 
    (324*g1^2*g2^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]])/5) + 
  M2*(174*g2^4*MatMul[Ye, Adj[he]] + 
    36*g2^4*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]) + 
  M1*((3006*g1^4*MatMul[he, Adj[Ye]])/25 + 
    (324*g1^4*MatMul[MatMul[he, Adj[Ye]], Zeta[3]])/5) + 
  M1*(54*g1^2*g2^2*MatMul[he, Adj[Ye]] - 
    (324*g1^2*g2^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]])/5) + 
  M2*(54*g1^2*g2^2*MatMul[he, Adj[Ye]] - 
    (324*g1^2*g2^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]])/5) + 
  M2*(174*g2^4*MatMul[he, Adj[Ye]] + 36*g2^4*MatMul[MatMul[he, Adj[Ye]], 
      Zeta[3]]) - (162*g1^4*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]])/5 + 
  (324*g1^2*g2^2*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]])/5 - 
  18*g2^4*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]] + 
  M1^2*((-9018*g1^4*MatMul[Ye, Adj[Ye]])/25 - 
    (972*g1^4*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]])/5) + 
  M1^2*(-108*g1^2*g2^2*MatMul[Ye, Adj[Ye]] + 
    (648*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]])/5) + 
  M1*M2*(-108*g1^2*g2^2*MatMul[Ye, Adj[Ye]] + 
    (648*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]])/5) + 
  M2^2*(-108*g1^2*g2^2*MatMul[Ye, Adj[Ye]] + 
    (648*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]])/5) + 
  M2^2*(-474*g2^4*MatMul[Ye, Adj[Ye]] - 108*g2^4*MatMul[MatMul[Ye, Adj[Ye]], 
      Zeta[3]]) + (18*g1^2*MatMul[MatMul[he, Adj[Ye], Ye], Adj[he]])/5 + 
  18*g2^2*MatMul[MatMul[he, Adj[Ye], Ye], Adj[he]] + 
  (108*g1^2*MatMul[MatMul[he, Adj[Ye], Ye], MatMul[Adj[he], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[he, Adj[Ye], Ye], MatMul[Adj[he], Zeta[3]]] - 
  (81*g1^4*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]])/5 + 
  (162*g1^2*g2^2*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]] - 
  (162*g1^4*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]])/5 + 
  (324*g1^2*g2^2*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]])/5 - 
  18*g2^4*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]] + 
  (18*g1^2*MatMul[MatMul[Ye, Adj[Ye], he], Adj[he]])/5 + 
  18*g2^2*MatMul[MatMul[Ye, Adj[Ye], he], Adj[he]] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye], he], MatMul[Adj[he], Zeta[3]]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye], he], MatMul[Adj[he], Zeta[3]]] - 
  (81*g1^4*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]])/5 + 
  (162*g1^2*g2^2*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]])/5 - 
  9*g2^4*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]] + 
  M1*((-18*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Adj[he]])/5 - 
    (108*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye], MatMul[Adj[he], Zeta[3]]])/5) + 
  M2*(-18*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], Adj[he]] + 
    36*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye], MatMul[Adj[he], Zeta[3]]]) + 
  (108*g1^2*MatMul[MatMul[he, Adj[he], Ye, Adj[Ye]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[he, Adj[he], Ye, Adj[Ye]], Zeta[3]] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[he], he, Adj[Ye]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[he], he, Adj[Ye]], Zeta[3]] + 
  (216*g1^2*mh1*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]])/5 - 
  72*g2^2*mh1*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  12*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye], Adj[he]] + 
  24*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye], Ye], MatMul[Adj[he], Zeta[3]]] + 
  (54*g1^2*MatMul[MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  (108*g1^2*MatMul[MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  12*MatMul[MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye], Adj[he]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], he, Adj[Ye], Ye], MatMul[Adj[he], Zeta[3]]] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]], Zeta[3]] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]], Zeta[3]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]], Zeta[3]] + 
  12*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he], Adj[he]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], he], MatMul[Adj[he], Zeta[3]]] + 
  (54*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me], Zeta[3]])/5 - 
  18*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me], Zeta[3]] + 
  24*MatMul[MatMul[he, Adj[he], Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[he, Adj[Ye], Ye, Adj[he], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[he], he, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[he], Ye, Adj[Ye], he, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], he, Adj[he], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[he], he, Adj[Ye]], Zeta[3]] + 
  72*mh1*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  12*MatMul[MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me, Ye, Adj[Ye]], Zeta[3]] + 
  24*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, ml, Adj[Ye]], Zeta[3]] + 
  12*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], me], Zeta[3]] - 
  (1503*g1^4*MatMul[me, Ye, Adj[Ye]])/50 - 
  27*g1^2*g2^2*MatMul[me, Ye, Adj[Ye]] - (87*g2^4*MatMul[me, Ye, Adj[Ye]])/
   2 - (1503*g1^4*MatMul[Ye, ml, Adj[Ye]])/25 - 
  54*g1^2*g2^2*MatMul[Ye, ml, Adj[Ye]] - 87*g2^4*MatMul[Ye, ml, Adj[Ye]] - 
  (1503*g1^4*MatMul[Ye, Adj[Ye], me])/50 - 
  27*g1^2*g2^2*MatMul[Ye, Adj[Ye], me] - (87*g2^4*MatMul[Ye, Adj[Ye], me])/
   2 + (18*g1^2*MatMul[he, Adj[he], Ye, Adj[Ye]])/5 + 
  18*g2^2*MatMul[he, Adj[he], Ye, Adj[Ye]] + 
  M1*((-108*g1^2*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye]], Zeta[3]])/5 - 
    (18*g1^2*MatMul[he, Adj[Ye], Ye, Adj[Ye]])/5) + 
  M2*(36*g2^2*MatMul[MatMul[he, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] - 
    18*g2^2*MatMul[he, Adj[Ye], Ye, Adj[Ye]]) + 
  (18*g1^2*MatMul[Ye, Adj[he], he, Adj[Ye]])/5 + 
  18*g2^2*MatMul[Ye, Adj[he], he, Adj[Ye]] + 
  M1*((-108*g1^2*MatMul[MatMul[Ye, Adj[he], Ye, Adj[Ye]], Zeta[3]])/5 - 
    (18*g1^2*MatMul[Ye, Adj[he], Ye, Adj[Ye]])/5) + 
  M2*(36*g2^2*MatMul[MatMul[Ye, Adj[he], Ye, Adj[Ye]], Zeta[3]] - 
    18*g2^2*MatMul[Ye, Adj[he], Ye, Adj[Ye]]) + 
  M1*((-108*g1^2*MatMul[MatMul[Ye, Adj[Ye], he, Adj[Ye]], Zeta[3]])/5 - 
    (18*g1^2*MatMul[Ye, Adj[Ye], he, Adj[Ye]])/5) + 
  M2*(36*g2^2*MatMul[MatMul[Ye, Adj[Ye], he, Adj[Ye]], Zeta[3]] - 
    18*g2^2*MatMul[Ye, Adj[Ye], he, Adj[Ye]]) + 
  (36*g1^2*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]])/5 + 
  36*g2^2*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]] + 
  M1^2*((216*g1^2*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]])/5 + 
    (36*g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]])/5) + 
  M2^2*(-72*g2^2*MatMul[MatMul[Ye, Adj[Ye], Ye, Adj[Ye]], Zeta[3]] + 
    36*g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]) + 
  (9*g1^2*MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]])/5 + 
  9*g2^2*MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]] + 
  (18*g1^2*MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]])/5 + 
  18*g2^2*MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]] + 
  (18*g1^2*MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]])/5 + 
  18*g2^2*MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]] + 
  (18*g1^2*MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]])/5 + 
  18*g2^2*MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]] + 
  (9*g1^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me])/5 + 
  9*g2^2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me] + 
  12*MatMul[he, Adj[he], Ye, Adj[Ye], Ye, Adj[Ye]] + 
  12*MatMul[he, Adj[Ye], Ye, Adj[he], Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[he], he, Adj[Ye], Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[he], Ye, Adj[Ye], he, Adj[Ye]] + 
  12*MatMul[Ye, Adj[Ye], he, Adj[he], Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[he], he, Adj[Ye]] + 
  36*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]] + 
  6*MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]] + 
  12*MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye], Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye], Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me, Ye, Adj[Ye]] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, ml, Adj[Ye]] + 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], me] - 
  (1872*g1^6*trace[mb])/125 + (24*g1^4*MatMul[Ye, Adj[Ye]]*trace[mb])/25 - 
  (5616*g1^6*trace[me])/125 + (72*g1^4*MatMul[Ye, Adj[Ye]]*trace[me])/25 - 
  (2808*g1^6*trace[ml])/125 + (36*g1^4*MatMul[Ye, Adj[Ye]]*trace[ml])/25 - 
  12*g2^4*MatMul[Ye, Adj[Ye]]*trace[ml] - (936*g1^6*trace[mq])/125 + 
  (12*g1^4*MatMul[Ye, Adj[Ye]]*trace[mq])/25 - 36*g2^4*MatMul[Ye, Adj[Ye]]*
   trace[mq] - (7488*g1^6*trace[mt])/125 + 
  (96*g1^4*MatMul[Ye, Adj[Ye]]*trace[mt])/25 - (504*g1^4*trace[hb, Adj[hb]])/
   25 + (214*g1^2*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[hb]])/5 + 
  54*g2^2*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[hb]] - 
  128*g3^2*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[hb]] + 
  (84*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[hb, Adj[hb]])/5 - 
  108*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[hb, Adj[hb]] + 
  192*g3^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[hb, Adj[hb]] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[hb, Adj[hb]] + 
  (214*g1^2*MatMul[Ye, Adj[he]]*trace[hb, Adj[Yb]])/5 + 
  54*g2^2*MatMul[Ye, Adj[he]]*trace[hb, Adj[Yb]] - 
  128*g3^2*MatMul[Ye, Adj[he]]*trace[hb, Adj[Yb]] + 
  (84*g1^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[hb, Adj[Yb]])/5 - 
  108*g2^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[hb, Adj[Yb]] + 
  192*g3^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[hb, Adj[Yb]] + 
  12*MatMul[MatMul[Ye, Adj[Ye], Ye], Adj[he]]*trace[hb, Adj[Yb]] + 
  12*MatMul[Ye, Adj[he], Ye, Adj[Ye]]*trace[hb, Adj[Yb]] - 
  (648*g1^4*trace[he, Adj[he]])/25 + 
  (18*g1^2*MatMul[Ye, Adj[Ye]]*trace[he, Adj[he]])/5 + 
  18*g2^2*MatMul[Ye, Adj[Ye]]*trace[he, Adj[he]] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[he, Adj[he]])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[he, Adj[he]] + 
  4*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[he, Adj[he]] + 
  (18*g1^2*MatMul[Ye, Adj[he]]*trace[he, Adj[Ye]])/5 + 
  18*g2^2*MatMul[Ye, Adj[he]]*trace[he, Adj[Ye]] + 
  (108*g1^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[he, Adj[Ye]])/5 - 
  36*g2^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[he, Adj[Ye]] + 
  4*MatMul[MatMul[Ye, Adj[Ye], Ye], Adj[he]]*trace[he, Adj[Ye]] + 
  4*MatMul[Ye, Adj[he], Ye, Adj[Ye]]*trace[he, Adj[Ye]] - 
  (936*g1^4*trace[ht, Adj[ht]])/25 + 
  (214*g1^2*MatMul[he, Adj[Ye]]*trace[Adj[hb], Yb])/5 + 
  54*g2^2*MatMul[he, Adj[Ye]]*trace[Adj[hb], Yb] - 
  128*g3^2*MatMul[he, Adj[Ye]]*trace[Adj[hb], Yb] + 
  (84*g1^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[hb], Yb])/5 - 
  108*g2^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[hb], Yb] + 
  192*g3^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[hb], Yb] + 
  12*MatMul[he, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[hb], Yb] + 
  12*MatMul[Ye, Adj[Ye], he, Adj[Ye]]*trace[Adj[hb], Yb] - 
  72*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  24*MatMul[Ye, Adj[Ye]]*trace[he, Adj[Ye]]*trace[Adj[hb], Yb] + 
  M3*(128*g3^2*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[Yb]] - 
    192*g3^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[hb, Adj[Yb]] + 
    128*g3^2*MatMul[Ye, Adj[Ye]]*trace[Adj[hb], Yb] - 
    192*g3^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[hb], Yb]) + 
  (18*g1^2*MatMul[he, Adj[Ye]]*trace[Adj[he], Ye])/5 + 
  18*g2^2*MatMul[he, Adj[Ye]]*trace[Adj[he], Ye] + 
  (108*g1^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[he], Ye])/5 - 
  36*g2^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[he], Ye] + 
  4*MatMul[he, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[he], Ye] + 
  4*MatMul[Ye, Adj[Ye], he, Adj[Ye]]*trace[Adj[he], Ye] - 
  24*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[Yb]]*trace[Adj[he], Ye] - 
  8*MatMul[Ye, Adj[Ye]]*trace[he, Adj[Ye]]*trace[Adj[he], Ye] + 
  M1*((-214*g1^2*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[Yb]])/5 - 
    (84*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[hb, Adj[Yb]])/5 - 
    (18*g1^2*MatMul[Ye, Adj[Ye]]*trace[he, Adj[Ye]])/5 - 
    (108*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[he, Adj[Ye]])/5 - 
    (214*g1^2*MatMul[Ye, Adj[Ye]]*trace[Adj[hb], Yb])/5 - 
    (84*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[hb], Yb])/5 - 
    (18*g1^2*MatMul[Ye, Adj[Ye]]*trace[Adj[he], Ye])/5 - 
    (108*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[he], Ye])/5) + 
  M2*(-54*g2^2*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[Yb]] + 
    108*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[hb, Adj[Yb]] - 
    18*g2^2*MatMul[Ye, Adj[Ye]]*trace[he, Adj[Ye]] + 
    36*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[he, Adj[Ye]] - 
    54*g2^2*MatMul[Ye, Adj[Ye]]*trace[Adj[hb], Yb] + 
    108*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[hb], Yb] - 
    18*g2^2*MatMul[Ye, Adj[Ye]]*trace[Adj[he], Ye] + 
    36*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[he], Ye]) + 
  M1*((168*g1^4*trace[hb, Adj[Yb]])/5 + (216*g1^4*trace[he, Adj[Ye]])/5 + 
    (312*g1^4*trace[ht, Adj[Yt]])/5 + (168*g1^4*trace[Adj[hb], Yb])/5 + 
    (216*g1^4*trace[Adj[he], Ye])/5 + (312*g1^4*trace[Adj[ht], Yt])/5) - 
  (504*g1^4*mh1*trace[Adj[Yb], Yb])/25 + 
  (214*g1^2*MatMul[he, Adj[he]]*trace[Adj[Yb], Yb])/5 + 
  54*g2^2*MatMul[he, Adj[he]]*trace[Adj[Yb], Yb] - 
  128*g3^2*MatMul[he, Adj[he]]*trace[Adj[Yb], Yb] + 
  (84*g1^2*MatMul[he, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Yb], Yb])/5 - 
  108*g2^2*MatMul[he, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Yb], Yb] + 
  192*g3^2*MatMul[he, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Yb], Yb] + 
  (428*g1^2*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb])/5 + 
  108*g2^2*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  256*g3^2*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  (168*g1^2*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  216*g2^2*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  384*g3^2*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  12*MatMul[MatMul[he, Adj[Ye], Ye], Adj[he]]*trace[Adj[Yb], Yb] + 
  (42*g1^2*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  54*g2^2*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (84*g1^2*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  108*g2^2*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  192*g3^2*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
  12*MatMul[MatMul[Ye, Adj[Ye], he], Adj[he]]*trace[Adj[Yb], Yb] + 
  (42*g1^2*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
  54*g2^2*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]]*trace[Adj[Yb], Yb] + 
  (107*g1^2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb])/5 + 
  27*g2^2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  64*g3^2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  (214*g1^2*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb])/5 + 
  54*g2^2*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  128*g3^2*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  (107*g1^2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb])/5 + 
  27*g2^2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb] - 
  64*g3^2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb] + 
  12*MatMul[he, Adj[he], Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Ye, Adj[he], he, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  36*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  6*MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  12*MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb] + 
  6*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me]*trace[Adj[Yb], Yb] - 
  72*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  72*MatMul[Ye, Adj[he]]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  24*MatMul[Ye, Adj[Ye]]*trace[he, Adj[he]]*trace[Adj[Yb], Yb] - 
  24*MatMul[Ye, Adj[he]]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  72*MatMul[he, Adj[Ye]]*trace[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 
  24*MatMul[he, Adj[Ye]]*trace[Adj[he], Ye]*trace[Adj[Yb], Yb] - 
  36*MatMul[he, Adj[he]]*trace[Adj[Yb], Yb]^2 - 
  108*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb]^2 - 
  36*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb]^2 + 
  M3*(128*g3^2*MatMul[Ye, Adj[he]]*trace[Adj[Yb], Yb] - 
    192*g3^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Yb], Yb]) + 
  M3*(128*g3^2*MatMul[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
    192*g3^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M3^2*(-256*g3^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb] + 
    384*g3^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb]) - 
  (648*g1^4*mh1*trace[Adj[Ye], Ye])/25 + 
  (18*g1^2*MatMul[he, Adj[he]]*trace[Adj[Ye], Ye])/5 + 
  18*g2^2*MatMul[he, Adj[he]]*trace[Adj[Ye], Ye] + 
  (108*g1^2*MatMul[he, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Ye], Ye])/5 - 
  36*g2^2*MatMul[he, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Ye], Ye] + 
  (36*g1^2*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye])/5 + 
  36*g2^2*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  (216*g1^2*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  72*g2^2*mh1*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  4*MatMul[MatMul[he, Adj[Ye], Ye], Adj[he]]*trace[Adj[Ye], Ye] + 
  (54*g1^2*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  18*g2^2*MatMul[MatMul[me, Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  (108*g1^2*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  36*g2^2*MatMul[MatMul[Ye, ml, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye] + 
  4*MatMul[MatMul[Ye, Adj[Ye], he], Adj[he]]*trace[Adj[Ye], Ye] + 
  (54*g1^2*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]]*trace[Adj[Ye], Ye])/5 - 
  18*g2^2*MatMul[MatMul[Ye, Adj[Ye], me], Zeta[3]]*trace[Adj[Ye], Ye] + 
  (9*g1^2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Ye], Ye])/5 + 
  9*g2^2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  (18*g1^2*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Ye], Ye])/5 + 
  18*g2^2*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  (9*g1^2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Ye], Ye])/5 + 
  9*g2^2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Ye], Ye] + 
  4*MatMul[he, Adj[he], Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Ye, Adj[he], he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  12*mh1*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  2*MatMul[me, Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Ye, ml, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Ye, Adj[Ye], me, Ye, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  4*MatMul[Ye, Adj[Ye], Ye, ml, Adj[Ye]]*trace[Adj[Ye], Ye] + 
  2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me]*trace[Adj[Ye], Ye] - 
  24*MatMul[Ye, Adj[Ye]]*trace[hb, Adj[hb]]*trace[Adj[Ye], Ye] - 
  24*MatMul[Ye, Adj[he]]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  8*MatMul[Ye, Adj[Ye]]*trace[he, Adj[he]]*trace[Adj[Ye], Ye] - 
  8*MatMul[Ye, Adj[he]]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  24*MatMul[he, Adj[Ye]]*trace[Adj[hb], Yb]*trace[Adj[Ye], Ye] - 
  8*MatMul[he, Adj[Ye]]*trace[Adj[he], Ye]*trace[Adj[Ye], Ye] - 
  24*MatMul[he, Adj[he]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  72*mh1*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  24*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[Ye, Adj[Ye], me]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  4*MatMul[he, Adj[he]]*trace[Adj[Ye], Ye]^2 - 12*mh1*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Ye], Ye]^2 - 2*MatMul[me, Ye, Adj[Ye]]*trace[Adj[Ye], Ye]^2 - 
  4*MatMul[Ye, ml, Adj[Ye]]*trace[Adj[Ye], Ye]^2 - 
  2*MatMul[Ye, Adj[Ye], me]*trace[Adj[Ye], Ye]^2 + 
  M1*((-214*g1^2*MatMul[Ye, Adj[he]]*trace[Adj[Yb], Yb])/5 - 
    (84*g1^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Yb], Yb])/5 - 
    (18*g1^2*MatMul[Ye, Adj[he]]*trace[Adj[Ye], Ye])/5 - 
    (108*g1^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Ye], Ye])/5) + 
  M2*(-54*g2^2*MatMul[Ye, Adj[he]]*trace[Adj[Yb], Yb] + 
    108*g2^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Yb], Yb] - 
    18*g2^2*MatMul[Ye, Adj[he]]*trace[Adj[Ye], Ye] + 
    36*g2^2*MatMul[Ye, MatMul[Adj[he], Zeta[3]]]*trace[Adj[Ye], Ye]) + 
  M1*((-214*g1^2*MatMul[he, Adj[Ye]]*trace[Adj[Yb], Yb])/5 - 
    (84*g1^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb])/5 - 
    (18*g1^2*MatMul[he, Adj[Ye]]*trace[Adj[Ye], Ye])/5 - 
    (108*g1^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M2*(-54*g2^2*MatMul[he, Adj[Ye]]*trace[Adj[Yb], Yb] + 
    108*g2^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] - 
    18*g2^2*MatMul[he, Adj[Ye]]*trace[Adj[Ye], Ye] + 
    36*g2^2*MatMul[MatMul[he, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye]) + 
  M1^2*((428*g1^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb])/5 + 
    (168*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb])/5 + 
    (36*g1^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye])/5 + 
    (216*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye])/5) + 
  M2^2*(108*g2^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb] - 
    216*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb] + 
    36*g2^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye] - 
    72*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye]) - 
  (936*g1^4*mh2*trace[Adj[Yt], Yt])/25 + 
  M1^2*((-504*g1^4*trace[Adj[Yb], Yb])/5 - (648*g1^4*trace[Adj[Ye], Ye])/5 - 
    (936*g1^4*trace[Adj[Yt], Yt])/5) - (504*g1^4*trace[Yb, Adj[Yb], mb])/25 + 
  (214*g1^2*MatMul[Ye, Adj[Ye]]*trace[Yb, Adj[Yb], mb])/5 + 
  54*g2^2*MatMul[Ye, Adj[Ye]]*trace[Yb, Adj[Yb], mb] - 
  128*g3^2*MatMul[Ye, Adj[Ye]]*trace[Yb, Adj[Yb], mb] + 
  (84*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Yb, Adj[Yb], mb])/5 - 
  108*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Yb, Adj[Yb], mb] + 
  192*g3^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Yb, Adj[Yb], mb] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Yb, Adj[Yb], mb] - 
  72*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Yb, Adj[Yb], mb] - 
  24*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
  (648*g1^4*trace[Ye, Adj[Ye], me])/25 + 
  (18*g1^2*MatMul[Ye, Adj[Ye]]*trace[Ye, Adj[Ye], me])/5 + 
  18*g2^2*MatMul[Ye, Adj[Ye]]*trace[Ye, Adj[Ye], me] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Ye, Adj[Ye], me])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Ye, Adj[Ye], me] + 
  4*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Ye, Adj[Ye], me] - 
  24*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Ye, Adj[Ye], me] - 
  8*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye]*trace[Ye, Adj[Ye], me] - 
  (936*g1^4*trace[Yt, Adj[Yt], mt])/25 - (504*g1^4*trace[Adj[Yb], Yb, mq])/
   25 + (214*g1^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb, mq])/5 + 
  54*g2^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb, mq] - 
  128*g3^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb, mq] + 
  (84*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb, mq])/5 - 
  108*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb, mq] + 
  192*g3^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Yb], Yb, mq] + 
  12*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Yb], Yb, mq] - 
  72*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Adj[Yb], Yb, mq] - 
  24*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye]*trace[Adj[Yb], Yb, mq] - 
  (648*g1^4*trace[Adj[Ye], Ye, ml])/25 + 
  (18*g1^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye, ml])/5 + 
  18*g2^2*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye, ml] + 
  (108*g1^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye, ml])/5 - 
  36*g2^2*MatMul[MatMul[Ye, Adj[Ye]], Zeta[3]]*trace[Adj[Ye], Ye, ml] + 
  4*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]]*trace[Adj[Ye], Ye, ml] - 
  24*MatMul[Ye, Adj[Ye]]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye, ml] - 
  8*MatMul[Ye, Adj[Ye]]*trace[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml] - 
  (936*g1^4*trace[Adj[Yt], Yt, mq])/25 + 24*MatMul[Ye, Adj[Ye]]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 24*MatMul[Ye, Adj[he]]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 144*MatMul[Ye, Adj[Ye]]*
   trace[Adj[hb], hb, Adj[Yb], Yb] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 48*MatMul[Ye, Adj[Ye]]*
   trace[Adj[he], he, Adj[Ye], Ye] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 24*MatMul[he, Adj[Ye]]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 144*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yb], hb, Adj[hb], Yb] + 144*MatMul[Ye, Adj[he]]*
   trace[Adj[Yb], hb, Adj[Yb], Yb] + 144*MatMul[he, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[hb], Yb] + 72*MatMul[he, Adj[he]]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 216*mh1*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 36*MatMul[me, Ye, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 72*MatMul[Ye, ml, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 36*MatMul[Ye, Adj[Ye], me]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 48*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Ye], he, Adj[he], Ye] + 48*MatMul[Ye, Adj[he]]*
   trace[Adj[Ye], he, Adj[Ye], Ye] + 48*MatMul[he, Adj[Ye]]*
   trace[Adj[Ye], Ye, Adj[he], Ye] + 24*MatMul[he, Adj[he]]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 72*mh1*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 12*MatMul[me, Ye, Adj[Ye]]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 24*MatMul[Ye, ml, Adj[Ye]]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 12*MatMul[Ye, Adj[Ye], me]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 24*MatMul[Ye, Adj[he]]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 24*MatMul[he, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] + 24*MatMul[he, Adj[he]]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 48*mh1*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 24*mh2*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*MatMul[me, Ye, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 24*MatMul[Ye, ml, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*MatMul[Ye, Adj[Ye], me]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 144*MatMul[Ye, Adj[Ye]]*
   trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 48*MatMul[Ye, Adj[Ye]]*
   trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 144*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 48*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 24*MatMul[Ye, Adj[Ye]]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
  M1^2*((191964*g1^6)/125 - (85968*g1^6*Zeta[3])/125) + 
  M1^2*((972*g1^4*g2^2)/5 - (5832*g1^4*g2^2*Zeta[3])/25) + 
  M1*M2*((648*g1^4*g2^2)/5 - (3888*g1^4*g2^2*Zeta[3])/25) + 
  M2^2*((1944*g1^4*g2^2)/25 - (1944*g1^4*g2^2*Zeta[3])/25) + 
  M1^2*((3168*g1^4*g3^2)/5 - (19008*g1^4*g3^2*Zeta[3])/25) + 
  M1*M3*((2112*g1^4*g3^2)/5 - (12672*g1^4*g3^2*Zeta[3])/25) + 
  M3^2*((6336*g1^4*g3^2)/25 - (6336*g1^4*g3^2*Zeta[3])/25)}
