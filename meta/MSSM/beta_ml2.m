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

{(-6*g1^2*M1^2)/5 - 6*g2^2*M2^2 + 2*MatMul[Adj[he], he] + 
  2*mh1*MatMul[Adj[Ye], Ye] + MatMul[ml, Adj[Ye], Ye] + 
  2*MatMul[Adj[Ye], me, Ye] + MatMul[Adj[Ye], Ye, ml], 
 (621*g1^4*M1^2)/25 + (18*g1^2*g2^2*M1^2)/5 + (18*g1^2*g2^2*M1*M2)/5 + 
  (18*g1^2*g2^2*M2^2)/5 + 33*g2^4*M2^2 + (9*g1^4*mh1)/25 + 3*g2^4*mh1 + 
  (9*g1^4*mh2)/25 + 3*g2^4*mh2 + (12*g1^2*MatMul[Adj[he], he])/5 - 
  (12*g1^2*M1*MatMul[Adj[he], Ye])/5 - (12*g1^2*M1*MatMul[Adj[Ye], he])/5 + 
  (24*g1^2*M1^2*MatMul[Adj[Ye], Ye])/5 + (12*g1^2*mh1*MatMul[Adj[Ye], Ye])/
   5 + (6*g1^2*MatMul[ml, Adj[Ye], Ye])/5 + (12*g1^2*MatMul[Adj[Ye], me, Ye])/
   5 + (6*g1^2*MatMul[Adj[Ye], Ye, ml])/5 - 
  4*MatMul[Adj[he], he, Adj[Ye], Ye] - 4*MatMul[Adj[he], Ye, Adj[Ye], he] - 
  4*MatMul[Adj[Ye], he, Adj[he], Ye] - 4*MatMul[Adj[Ye], Ye, Adj[he], he] - 
  8*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye] - 
  2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye] - 
  2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml] + (6*g1^4*trace[mb])/25 + 
  (18*g1^4*trace[me])/25 + (9*g1^4*trace[ml])/25 + 3*g2^4*trace[ml] + 
  (3*g1^4*trace[mq])/25 + 9*g2^4*trace[mq] + (24*g1^4*trace[mt])/25 - 
  6*MatMul[Adj[Ye], Ye]*trace[hb, Adj[hb]] - 6*MatMul[Adj[he], Ye]*
   trace[hb, Adj[Yb]] - 2*MatMul[Adj[Ye], Ye]*trace[he, Adj[he]] - 
  2*MatMul[Adj[he], Ye]*trace[he, Adj[Ye]] - 6*MatMul[Adj[Ye], he]*
   trace[Adj[hb], Yb] - 2*MatMul[Adj[Ye], he]*trace[Adj[he], Ye] - 
  6*MatMul[Adj[he], he]*trace[Adj[Yb], Yb] - 12*mh1*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb] - 3*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  6*MatMul[Adj[Ye], me, Ye]*trace[Adj[Yb], Yb] - 
  3*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb] - 
  2*MatMul[Adj[he], he]*trace[Adj[Ye], Ye] - 4*mh1*MatMul[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye] - MatMul[ml, Adj[Ye], Ye]*trace[Adj[Ye], Ye] - 
  2*MatMul[Adj[Ye], me, Ye]*trace[Adj[Ye], Ye] - 
  MatMul[Adj[Ye], Ye, ml]*trace[Adj[Ye], Ye] - 
  6*MatMul[Adj[Ye], Ye]*trace[Yb, Adj[Yb], mb] - 
  2*MatMul[Adj[Ye], Ye]*trace[Ye, Adj[Ye], me] - 
  6*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb, mq] - 
  2*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye, ml], 
 (-621*g1^6*mh1)/125 - (27*g1^4*g2^2*mh1)/25 - (9*g1^2*g2^4*mh1)/5 - 
  3*g2^6*mh1 - (621*g1^6*mh2)/125 - (27*g1^4*g2^2*mh2)/25 - 
  (9*g1^2*g2^4*mh2)/5 - 3*g2^6*mh2 - (549*g1^4*MatMul[Adj[he], he])/10 - 
  (81*g1^2*g2^2*MatMul[Adj[he], he])/5 - (45*g2^4*MatMul[Adj[he], he])/2 - 
  (2817*g1^4*mh1*MatMul[Adj[Ye], Ye])/50 - 
  (81*g1^2*g2^2*mh1*MatMul[Adj[Ye], Ye])/5 - 
  (45*g2^4*mh1*MatMul[Adj[Ye], Ye])/2 - (36*g1^4*mh2*MatMul[Adj[Ye], Ye])/
   25 + (81*g1^4*MatMul[MatMul[Adj[he], he], Zeta[3]])/25 + 
  (162*g1^2*g2^2*MatMul[MatMul[Adj[he], he], Zeta[3]])/5 - 
  63*g2^4*MatMul[MatMul[Adj[he], he], Zeta[3]] + 
  M1*((549*g1^4*MatMul[Adj[he], Ye])/5 - 
    (162*g1^4*MatMul[MatMul[Adj[he], Ye], Zeta[3]])/25) + 
  M1*((81*g1^2*g2^2*MatMul[Adj[he], Ye])/5 - 
    (162*g1^2*g2^2*MatMul[MatMul[Adj[he], Ye], Zeta[3]])/5) + 
  M2*((81*g1^2*g2^2*MatMul[Adj[he], Ye])/5 - 
    (162*g1^2*g2^2*MatMul[MatMul[Adj[he], Ye], Zeta[3]])/5) + 
  M2*(45*g2^4*MatMul[Adj[he], Ye] + 126*g2^4*MatMul[MatMul[Adj[he], Ye], 
      Zeta[3]]) + M1*((549*g1^4*MatMul[Adj[Ye], he])/5 - 
    (162*g1^4*MatMul[MatMul[Adj[Ye], he], Zeta[3]])/25) + 
  M1*((81*g1^2*g2^2*MatMul[Adj[Ye], he])/5 - 
    (162*g1^2*g2^2*MatMul[MatMul[Adj[Ye], he], Zeta[3]])/5) + 
  M2*((81*g1^2*g2^2*MatMul[Adj[Ye], he])/5 - 
    (162*g1^2*g2^2*MatMul[MatMul[Adj[Ye], he], Zeta[3]])/5) + 
  M2*(45*g2^4*MatMul[Adj[Ye], he] + 126*g2^4*MatMul[MatMul[Adj[Ye], he], 
      Zeta[3]]) + (81*g1^4*mh1*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]])/25 + 
  (162*g1^2*g2^2*mh1*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]])/5 - 
  63*g2^4*mh1*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]] + 
  M1^2*((-1647*g1^4*MatMul[Adj[Ye], Ye])/5 + 
    (486*g1^4*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]])/25) + 
  M1^2*((-162*g1^2*g2^2*MatMul[Adj[Ye], Ye])/5 + 
    (324*g1^2*g2^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]])/5) + 
  M1*M2*((-162*g1^2*g2^2*MatMul[Adj[Ye], Ye])/5 + 
    (324*g1^2*g2^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]])/5) + 
  M2^2*((-162*g1^2*g2^2*MatMul[Adj[Ye], Ye])/5 + 
    (324*g1^2*g2^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]])/5) + 
  M2^2*(-135*g2^4*MatMul[Adj[Ye], Ye] - 378*g2^4*MatMul[MatMul[Adj[Ye], Ye], 
      Zeta[3]]) + (81*g1^4*MatMul[MatMul[ml, Adj[Ye], Ye], Zeta[3]])/50 + 
  (81*g1^2*g2^2*MatMul[MatMul[ml, Adj[Ye], Ye], Zeta[3]])/5 - 
  (63*g2^4*MatMul[MatMul[ml, Adj[Ye], Ye], Zeta[3]])/2 + 
  (81*g1^4*MatMul[MatMul[Adj[Ye], me, Ye], Zeta[3]])/25 + 
  (162*g1^2*g2^2*MatMul[MatMul[Adj[Ye], me, Ye], Zeta[3]])/5 - 
  63*g2^4*MatMul[MatMul[Adj[Ye], me, Ye], Zeta[3]] + 
  (81*g1^4*MatMul[MatMul[Adj[Ye], Ye, ml], Zeta[3]])/50 + 
  (81*g1^2*g2^2*MatMul[MatMul[Adj[Ye], Ye, ml], Zeta[3]])/5 - 
  (63*g2^4*MatMul[MatMul[Adj[Ye], Ye, ml], Zeta[3]])/2 - 
  (108*g1^2*MatMul[MatMul[Adj[he], he, Adj[Ye], Ye], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[he], he, Adj[Ye], Ye], Zeta[3]] - 
  (108*g1^2*MatMul[MatMul[Adj[he], Ye, Adj[Ye], he], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[he], Ye, Adj[Ye], he], Zeta[3]] - 
  (108*g1^2*MatMul[MatMul[Adj[Ye], he, Adj[he], Ye], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Ye], he, Adj[he], Ye], Zeta[3]] - 
  (108*g1^2*MatMul[MatMul[Adj[Ye], Ye, Adj[he], he], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Ye], Ye, Adj[he], he], Zeta[3]] - 
  (216*g1^2*mh1*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]])/5 + 
  72*g2^2*mh1*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] - 
  (54*g1^2*MatMul[MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] - 
  (108*g1^2*MatMul[MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye], Zeta[3]] - 
  (108*g1^2*MatMul[MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye], Zeta[3]] - 
  (108*g1^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye], Zeta[3]])/5 + 
  36*g2^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye], Zeta[3]] - 
  (54*g1^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml], Zeta[3]])/5 + 
  18*g2^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml], Zeta[3]] + 
  12*MatMul[MatMul[Adj[he], he, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[he], Ye, Adj[Ye], he, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[he], Ye, Adj[Ye], Ye, Adj[Ye], he], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], he, Adj[he], Ye, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], he, Adj[Ye], Ye, Adj[he], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[he], he, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[he], Ye, Adj[Ye], he], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], he, Adj[he], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[he], he], Zeta[3]] + 
  36*mh1*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  6*MatMul[MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml, Adj[Ye], Ye], Zeta[3]] + 
  12*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], me, Ye], Zeta[3]] + 
  6*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye, Adj[Ye], Ye, ml], Zeta[3]] - 
  (549*g1^4*MatMul[ml, Adj[Ye], Ye])/20 - 
  (81*g1^2*g2^2*MatMul[ml, Adj[Ye], Ye])/10 - 
  (45*g2^4*MatMul[ml, Adj[Ye], Ye])/4 - (549*g1^4*MatMul[Adj[Ye], me, Ye])/
   10 - (81*g1^2*g2^2*MatMul[Adj[Ye], me, Ye])/5 - 
  (45*g2^4*MatMul[Adj[Ye], me, Ye])/2 - (549*g1^4*MatMul[Adj[Ye], Ye, ml])/
   20 - (81*g1^2*g2^2*MatMul[Adj[Ye], Ye, ml])/10 - 
  (45*g2^4*MatMul[Adj[Ye], Ye, ml])/4 + 
  18*g1^2*MatMul[Adj[he], he, Adj[Ye], Ye] - 
  6*g2^2*MatMul[Adj[he], he, Adj[Ye], Ye] + 
  18*g1^2*MatMul[Adj[he], Ye, Adj[Ye], he] - 
  6*g2^2*MatMul[Adj[he], Ye, Adj[Ye], he] + 
  M1*((108*g1^2*MatMul[MatMul[Adj[he], Ye, Adj[Ye], Ye], Zeta[3]])/5 - 
    18*g1^2*MatMul[Adj[he], Ye, Adj[Ye], Ye]) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[he], Ye, Adj[Ye], Ye], Zeta[3]] + 
    6*g2^2*MatMul[Adj[he], Ye, Adj[Ye], Ye]) + 
  18*g1^2*MatMul[Adj[Ye], he, Adj[he], Ye] - 
  6*g2^2*MatMul[Adj[Ye], he, Adj[he], Ye] + 
  M1*((108*g1^2*MatMul[MatMul[Adj[Ye], he, Adj[Ye], Ye], Zeta[3]])/5 - 
    18*g1^2*MatMul[Adj[Ye], he, Adj[Ye], Ye]) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Ye], he, Adj[Ye], Ye], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Ye], he, Adj[Ye], Ye]) + 
  18*g1^2*MatMul[Adj[Ye], Ye, Adj[he], he] - 
  6*g2^2*MatMul[Adj[Ye], Ye, Adj[he], he] + 
  M1*((108*g1^2*MatMul[MatMul[Adj[Ye], Ye, Adj[he], Ye], Zeta[3]])/5 - 
    18*g1^2*MatMul[Adj[Ye], Ye, Adj[he], Ye]) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Ye], Ye, Adj[he], Ye], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Ye], Ye, Adj[he], Ye]) + 
  M1*((108*g1^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], he], Zeta[3]])/5 - 
    18*g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], he]) + 
  M2*(-36*g2^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], he], Zeta[3]] + 
    6*g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], he]) + 
  36*g1^2*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye] - 
  12*g2^2*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye] + 
  M1^2*((-216*g1^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]])/5 + 
    36*g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  M2^2*(72*g2^2*MatMul[MatMul[Adj[Ye], Ye, Adj[Ye], Ye], Zeta[3]] - 
    12*g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]) + 
  9*g1^2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye] - 
  3*g2^2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye] + 
  18*g1^2*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye] - 
  6*g2^2*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye] + 
  18*g1^2*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye] - 
  6*g2^2*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye] + 
  18*g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye] - 
  6*g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye] + 
  9*g1^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml] - 
  3*g2^2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml] - (414*g1^6*trace[mb])/125 - 
  (18*g1^4*g2^2*trace[mb])/25 - (24*g1^4*MatMul[Adj[Ye], Ye]*trace[mb])/25 - 
  (1242*g1^6*trace[me])/125 - (54*g1^4*g2^2*trace[me])/25 - 
  (72*g1^4*MatMul[Adj[Ye], Ye]*trace[me])/25 - (621*g1^6*trace[ml])/125 - 
  (27*g1^4*g2^2*trace[ml])/25 - (9*g1^2*g2^4*trace[ml])/5 - 
  3*g2^6*trace[ml] - (36*g1^4*MatMul[Adj[Ye], Ye]*trace[ml])/25 - 
  (207*g1^6*trace[mq])/125 - (9*g1^4*g2^2*trace[mq])/25 - 
  (27*g1^2*g2^4*trace[mq])/5 - 9*g2^6*trace[mq] - 
  (12*g1^4*MatMul[Adj[Ye], Ye]*trace[mq])/25 - (1656*g1^6*trace[mt])/125 - 
  (72*g1^4*g2^2*trace[mt])/25 - (96*g1^4*MatMul[Adj[Ye], Ye]*trace[mt])/25 - 
  (126*g1^4*trace[hb, Adj[hb]])/25 - 54*g2^4*trace[hb, Adj[hb]] + 
  16*g1^2*MatMul[Adj[Ye], Ye]*trace[hb, Adj[hb]] + 
  36*g2^2*MatMul[Adj[Ye], Ye]*trace[hb, Adj[hb]] - 
  64*g3^2*MatMul[Adj[Ye], Ye]*trace[hb, Adj[hb]] - 
  24*g1^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[hb]] + 
  96*g3^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[hb]] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*trace[hb, Adj[hb]] + 
  16*g1^2*MatMul[Adj[he], Ye]*trace[hb, Adj[Yb]] + 
  36*g2^2*MatMul[Adj[he], Ye]*trace[hb, Adj[Yb]] - 
  64*g3^2*MatMul[Adj[he], Ye]*trace[hb, Adj[Yb]] - 
  24*g1^2*MatMul[MatMul[Adj[he], Ye], Zeta[3]]*trace[hb, Adj[Yb]] + 
  96*g3^2*MatMul[MatMul[Adj[he], Ye], Zeta[3]]*trace[hb, Adj[Yb]] + 
  12*MatMul[Adj[he], Ye, Adj[Ye], Ye]*trace[hb, Adj[Yb]] + 
  12*MatMul[Adj[Ye], Ye, Adj[he], Ye]*trace[hb, Adj[Yb]] - 
  (162*g1^4*trace[he, Adj[he]])/25 - 18*g2^4*trace[he, Adj[he]] + 
  12*g2^2*MatMul[Adj[Ye], Ye]*trace[he, Adj[he]] + 
  4*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*trace[he, Adj[he]] + 
  12*g2^2*MatMul[Adj[he], Ye]*trace[he, Adj[Ye]] + 
  4*MatMul[Adj[he], Ye, Adj[Ye], Ye]*trace[he, Adj[Ye]] + 
  4*MatMul[Adj[Ye], Ye, Adj[he], Ye]*trace[he, Adj[Ye]] - 
  (234*g1^4*trace[ht, Adj[ht]])/25 - 54*g2^4*trace[ht, Adj[ht]] + 
  16*g1^2*MatMul[Adj[Ye], he]*trace[Adj[hb], Yb] + 
  36*g2^2*MatMul[Adj[Ye], he]*trace[Adj[hb], Yb] - 
  64*g3^2*MatMul[Adj[Ye], he]*trace[Adj[hb], Yb] - 
  24*g1^2*MatMul[MatMul[Adj[Ye], he], Zeta[3]]*trace[Adj[hb], Yb] + 
  96*g3^2*MatMul[MatMul[Adj[Ye], he], Zeta[3]]*trace[Adj[hb], Yb] + 
  12*MatMul[Adj[Ye], he, Adj[Ye], Ye]*trace[Adj[hb], Yb] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], he]*trace[Adj[hb], Yb] - 
  36*MatMul[Adj[Ye], Ye]*trace[hb, Adj[Yb]]*trace[Adj[hb], Yb] - 
  12*MatMul[Adj[Ye], Ye]*trace[he, Adj[Ye]]*trace[Adj[hb], Yb] + 
  M1*(-16*g1^2*MatMul[Adj[Ye], Ye]*trace[hb, Adj[Yb]] + 
    24*g1^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[Yb]] - 
    16*g1^2*MatMul[Adj[Ye], Ye]*trace[Adj[hb], Yb] + 
    24*g1^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[Adj[hb], Yb]) + 
  M3*(64*g3^2*MatMul[Adj[Ye], Ye]*trace[hb, Adj[Yb]] - 
    96*g3^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[hb, Adj[Yb]] + 
    64*g3^2*MatMul[Adj[Ye], Ye]*trace[Adj[hb], Yb] - 
    96*g3^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[Adj[hb], Yb]) + 
  12*g2^2*MatMul[Adj[Ye], he]*trace[Adj[he], Ye] + 
  4*MatMul[Adj[Ye], he, Adj[Ye], Ye]*trace[Adj[he], Ye] + 
  4*MatMul[Adj[Ye], Ye, Adj[Ye], he]*trace[Adj[he], Ye] - 
  12*MatMul[Adj[Ye], Ye]*trace[hb, Adj[Yb]]*trace[Adj[he], Ye] - 
  4*MatMul[Adj[Ye], Ye]*trace[he, Adj[Ye]]*trace[Adj[he], Ye] + 
  M2*(-36*g2^2*MatMul[Adj[Ye], Ye]*trace[hb, Adj[Yb]] - 
    12*g2^2*MatMul[Adj[Ye], Ye]*trace[he, Adj[Ye]] - 
    36*g2^2*MatMul[Adj[Ye], Ye]*trace[Adj[hb], Yb] - 
    12*g2^2*MatMul[Adj[Ye], Ye]*trace[Adj[he], Ye]) + 
  M1*((42*g1^4*trace[hb, Adj[Yb]])/5 + (54*g1^4*trace[he, Adj[Ye]])/5 + 
    (78*g1^4*trace[ht, Adj[Yt]])/5 + (42*g1^4*trace[Adj[hb], Yb])/5 + 
    (54*g1^4*trace[Adj[he], Ye])/5 + (78*g1^4*trace[Adj[ht], Yt])/5) + 
  M2*(90*g2^4*trace[hb, Adj[Yb]] + 30*g2^4*trace[he, Adj[Ye]] + 
    90*g2^4*trace[ht, Adj[Yt]] + 90*g2^4*trace[Adj[hb], Yb] + 
    30*g2^4*trace[Adj[he], Ye] + 90*g2^4*trace[Adj[ht], Yt]) - 
  (126*g1^4*mh1*trace[Adj[Yb], Yb])/25 - 54*g2^4*mh1*trace[Adj[Yb], Yb] + 
  16*g1^2*MatMul[Adj[he], he]*trace[Adj[Yb], Yb] + 
  36*g2^2*MatMul[Adj[he], he]*trace[Adj[Yb], Yb] - 
  64*g3^2*MatMul[Adj[he], he]*trace[Adj[Yb], Yb] + 
  32*g1^2*mh1*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  72*g2^2*mh1*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  128*g3^2*mh1*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  24*g1^2*MatMul[MatMul[Adj[he], he], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[Adj[he], he], Zeta[3]]*trace[Adj[Yb], Yb] - 
  48*g1^2*mh1*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
  192*g3^2*mh1*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] - 
  12*g1^2*MatMul[MatMul[ml, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
  48*g3^2*MatMul[MatMul[ml, Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb] - 
  24*g1^2*MatMul[MatMul[Adj[Ye], me, Ye], Zeta[3]]*trace[Adj[Yb], Yb] + 
  96*g3^2*MatMul[MatMul[Adj[Ye], me, Ye], Zeta[3]]*trace[Adj[Yb], Yb] - 
  12*g1^2*MatMul[MatMul[Adj[Ye], Ye, ml], Zeta[3]]*trace[Adj[Yb], Yb] + 
  48*g3^2*MatMul[MatMul[Adj[Ye], Ye, ml], Zeta[3]]*trace[Adj[Yb], Yb] + 
  8*g1^2*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  18*g2^2*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
  32*g3^2*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  16*g1^2*MatMul[Adj[Ye], me, Ye]*trace[Adj[Yb], Yb] + 
  36*g2^2*MatMul[Adj[Ye], me, Ye]*trace[Adj[Yb], Yb] - 
  64*g3^2*MatMul[Adj[Ye], me, Ye]*trace[Adj[Yb], Yb] + 
  8*g1^2*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb] + 
  18*g2^2*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb] - 
  32*g3^2*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[he], he, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[he], Ye, Adj[Ye], he]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Ye], he, Adj[he], Ye]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Ye], Ye, Adj[he], he]*trace[Adj[Yb], Yb] + 
  36*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  6*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
  12*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye]*trace[Adj[Yb], Yb] + 
  6*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb] - 
  36*MatMul[Adj[Ye], Ye]*trace[hb, Adj[hb]]*trace[Adj[Yb], Yb] - 
  36*MatMul[Adj[he], Ye]*trace[hb, Adj[Yb]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Adj[Ye], Ye]*trace[he, Adj[he]]*trace[Adj[Yb], Yb] - 
  12*MatMul[Adj[he], Ye]*trace[he, Adj[Ye]]*trace[Adj[Yb], Yb] - 
  36*MatMul[Adj[Ye], he]*trace[Adj[hb], Yb]*trace[Adj[Yb], Yb] - 
  12*MatMul[Adj[Ye], he]*trace[Adj[he], Ye]*trace[Adj[Yb], Yb] - 
  18*MatMul[Adj[he], he]*trace[Adj[Yb], Yb]^2 - 
  54*mh1*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb]^2 - 
  9*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb]^2 - 
  18*MatMul[Adj[Ye], me, Ye]*trace[Adj[Yb], Yb]^2 - 
  9*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb]^2 + 
  M1*(-16*g1^2*MatMul[Adj[he], Ye]*trace[Adj[Yb], Yb] + 
    24*g1^2*MatMul[MatMul[Adj[he], Ye], Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M3*(64*g3^2*MatMul[Adj[he], Ye]*trace[Adj[Yb], Yb] - 
    96*g3^2*MatMul[MatMul[Adj[he], Ye], Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M1*(-16*g1^2*MatMul[Adj[Ye], he]*trace[Adj[Yb], Yb] + 
    24*g1^2*MatMul[MatMul[Adj[Ye], he], Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M3*(64*g3^2*MatMul[Adj[Ye], he]*trace[Adj[Yb], Yb] - 
    96*g3^2*MatMul[MatMul[Adj[Ye], he], Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M1^2*(32*g1^2*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb] - 
    48*g1^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb]) + 
  M3^2*(-128*g3^2*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
    192*g3^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*trace[Adj[Yb], Yb]) - 
  (162*g1^4*mh1*trace[Adj[Ye], Ye])/25 - 18*g2^4*mh1*trace[Adj[Ye], Ye] + 
  12*g2^2*MatMul[Adj[he], he]*trace[Adj[Ye], Ye] + 
  24*g2^2*mh1*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  6*g2^2*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  12*g2^2*MatMul[Adj[Ye], me, Ye]*trace[Adj[Ye], Ye] + 
  6*g2^2*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[he], he, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[he], Ye, Adj[Ye], he]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Ye], he, Adj[he], Ye]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Ye], Ye, Adj[he], he]*trace[Adj[Ye], Ye] + 
  12*mh1*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  2*MatMul[ml, Adj[Ye], Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Ye], me, Ye, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Ye], Ye, ml, Adj[Ye], Ye]*trace[Adj[Ye], Ye] + 
  4*MatMul[Adj[Ye], Ye, Adj[Ye], me, Ye]*trace[Adj[Ye], Ye] + 
  2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[Ye], Ye]*trace[hb, Adj[hb]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[he], Ye]*trace[hb, Adj[Yb]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], Ye]*trace[he, Adj[he]]*trace[Adj[Ye], Ye] - 
  4*MatMul[Adj[he], Ye]*trace[he, Adj[Ye]]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[Ye], he]*trace[Adj[hb], Yb]*trace[Adj[Ye], Ye] - 
  4*MatMul[Adj[Ye], he]*trace[Adj[he], Ye]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[he], he]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  36*mh1*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  6*MatMul[ml, Adj[Ye], Ye]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  12*MatMul[Adj[Ye], me, Ye]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  6*MatMul[Adj[Ye], Ye, ml]*trace[Adj[Yb], Yb]*trace[Adj[Ye], Ye] - 
  2*MatMul[Adj[he], he]*trace[Adj[Ye], Ye]^2 - 6*mh1*MatMul[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye]^2 - MatMul[ml, Adj[Ye], Ye]*trace[Adj[Ye], Ye]^2 - 
  2*MatMul[Adj[Ye], me, Ye]*trace[Adj[Ye], Ye]^2 - 
  MatMul[Adj[Ye], Ye, ml]*trace[Adj[Ye], Ye]^2 + 
  M2*(-36*g2^2*MatMul[Adj[he], Ye]*trace[Adj[Yb], Yb] - 
    12*g2^2*MatMul[Adj[he], Ye]*trace[Adj[Ye], Ye]) + 
  M2*(-36*g2^2*MatMul[Adj[Ye], he]*trace[Adj[Yb], Yb] - 
    12*g2^2*MatMul[Adj[Ye], he]*trace[Adj[Ye], Ye]) + 
  M2^2*(72*g2^2*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb] + 
    24*g2^2*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye]) - 
  (234*g1^4*mh2*trace[Adj[Yt], Yt])/25 - 54*g2^4*mh2*trace[Adj[Yt], Yt] + 
  M1^2*((-126*g1^4*trace[Adj[Yb], Yb])/5 - (162*g1^4*trace[Adj[Ye], Ye])/5 - 
    (234*g1^4*trace[Adj[Yt], Yt])/5) + 
  M2^2*(-270*g2^4*trace[Adj[Yb], Yb] - 90*g2^4*trace[Adj[Ye], Ye] - 
    270*g2^4*trace[Adj[Yt], Yt]) - (126*g1^4*trace[Yb, Adj[Yb], mb])/25 - 
  54*g2^4*trace[Yb, Adj[Yb], mb] + 16*g1^2*MatMul[Adj[Ye], Ye]*
   trace[Yb, Adj[Yb], mb] + 36*g2^2*MatMul[Adj[Ye], Ye]*
   trace[Yb, Adj[Yb], mb] - 64*g3^2*MatMul[Adj[Ye], Ye]*
   trace[Yb, Adj[Yb], mb] - 24*g1^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*
   trace[Yb, Adj[Yb], mb] + 96*g3^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*
   trace[Yb, Adj[Yb], mb] + 12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*
   trace[Yb, Adj[Yb], mb] - 36*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb]*
   trace[Yb, Adj[Yb], mb] - 12*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye]*
   trace[Yb, Adj[Yb], mb] - (162*g1^4*trace[Ye, Adj[Ye], me])/25 - 
  18*g2^4*trace[Ye, Adj[Ye], me] + 12*g2^2*MatMul[Adj[Ye], Ye]*
   trace[Ye, Adj[Ye], me] + 4*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*
   trace[Ye, Adj[Ye], me] - 12*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb]*
   trace[Ye, Adj[Ye], me] - 4*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye]*
   trace[Ye, Adj[Ye], me] - (234*g1^4*trace[Yt, Adj[Yt], mt])/25 - 
  54*g2^4*trace[Yt, Adj[Yt], mt] - (126*g1^4*trace[Adj[Yb], Yb, mq])/25 - 
  54*g2^4*trace[Adj[Yb], Yb, mq] + 16*g1^2*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, mq] + 36*g2^2*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, mq] - 64*g3^2*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, mq] - 24*g1^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*
   trace[Adj[Yb], Yb, mq] + 96*g3^2*MatMul[MatMul[Adj[Ye], Ye], Zeta[3]]*
   trace[Adj[Yb], Yb, mq] + 12*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, mq] - 36*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb]*
   trace[Adj[Yb], Yb, mq] - 12*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, mq] - (162*g1^4*trace[Adj[Ye], Ye, ml])/25 - 
  18*g2^4*trace[Adj[Ye], Ye, ml] + 12*g2^2*MatMul[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, ml] + 4*MatMul[Adj[Ye], Ye, Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, ml] - 12*MatMul[Adj[Ye], Ye]*trace[Adj[Yb], Yb]*
   trace[Adj[Ye], Ye, ml] - 4*MatMul[Adj[Ye], Ye]*trace[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, ml] - (234*g1^4*trace[Adj[Yt], Yt, mq])/25 - 
  54*g2^4*trace[Adj[Yt], Yt, mq] + 12*MatMul[Adj[Ye], Ye]*
   trace[hb, Adj[ht], Yt, Adj[Yb]] + 12*MatMul[Adj[he], Ye]*
   trace[hb, Adj[Yt], Yt, Adj[Yb]] + 72*MatMul[Adj[Ye], Ye]*
   trace[Adj[hb], hb, Adj[Yb], Yb] + 12*MatMul[Adj[Ye], Ye]*
   trace[Adj[hb], hb, Adj[Yt], Yt] + 24*MatMul[Adj[Ye], Ye]*
   trace[Adj[he], he, Adj[Ye], Ye] + 12*MatMul[Adj[Ye], Ye]*
   trace[Adj[ht], ht, Adj[Yb], Yb] + 12*MatMul[Adj[Ye], he]*
   trace[Adj[ht], Yt, Adj[Yb], Yb] + 72*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], hb, Adj[hb], Yb] + 72*MatMul[Adj[he], Ye]*
   trace[Adj[Yb], hb, Adj[Yb], Yb] + 72*MatMul[Adj[Ye], he]*
   trace[Adj[Yb], Yb, Adj[hb], Yb] + 36*MatMul[Adj[he], he]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 108*mh1*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 18*MatMul[ml, Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 36*MatMul[Adj[Ye], me, Ye]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 18*MatMul[Adj[Ye], Ye, ml]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb] + 24*MatMul[Adj[Ye], Ye]*
   trace[Adj[Ye], he, Adj[he], Ye] + 24*MatMul[Adj[he], Ye]*
   trace[Adj[Ye], he, Adj[Ye], Ye] + 24*MatMul[Adj[Ye], he]*
   trace[Adj[Ye], Ye, Adj[he], Ye] + 12*MatMul[Adj[he], he]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 36*mh1*MatMul[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*MatMul[ml, Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 12*MatMul[Adj[Ye], me, Ye]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 6*MatMul[Adj[Ye], Ye, ml]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye] + 12*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yt], ht, Adj[hb], Yb] + 12*MatMul[Adj[he], Ye]*
   trace[Adj[Yt], ht, Adj[Yb], Yb] + 12*MatMul[Adj[Ye], he]*
   trace[Adj[Yt], Yt, Adj[hb], Yb] + 12*MatMul[Adj[he], he]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 24*mh1*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*mh2*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*MatMul[ml, Adj[Ye], Ye]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 12*MatMul[Adj[Ye], me, Ye]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 6*MatMul[Adj[Ye], Ye, ml]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb] + 72*MatMul[Adj[Ye], Ye]*
   trace[Yb, Adj[Yb], Yb, Adj[Yb], mb] + 12*MatMul[Adj[Ye], Ye]*
   trace[Yb, Adj[Yt], Yt, Adj[Yb], mb] + 24*MatMul[Adj[Ye], Ye]*
   trace[Ye, Adj[Ye], Ye, Adj[Ye], me] + 12*MatMul[Adj[Ye], Ye]*
   trace[Yt, Adj[Yb], Yb, Adj[Yt], mt] + 72*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, Adj[Yb], Yb, mq] + 12*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yb], Yb, Adj[Yt], Yt, mq] + 24*MatMul[Adj[Ye], Ye]*
   trace[Adj[Ye], Ye, Adj[Ye], Ye, ml] + 12*MatMul[Adj[Ye], Ye]*
   trace[Adj[Yt], Yt, Adj[Yb], Yb, mq] + 
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
  M3^2*(432*g2^4*g3^2 - 432*g2^4*g3^2*Zeta[3])}
