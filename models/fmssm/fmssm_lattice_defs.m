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

<< meta/writeRGE.m

norm[z_] := z Conjugate[z];
H[m_] := ConjugateTranspose[m];
TR[m_] := Re[Tr[m]];

mapVars[{
real[t],
real[g1, g2, g3],
cmat[3,3, Yu, Yd, Ye],
comp[mu],
comp[b],
comp[M1, M2, M3],
real[m2Hu, m2Hd],
hmat[3, m2Q, m2U, m2D, m2L, m2E],
cmat[3,3, TAu, TAd, TAe]
}];
Yn := 0 U;
m2N := 0 U;
TAn := 0 U;

U := IdentityMatrix[3];

lsf := 1/(16 Pi^2);

b1 := 33/5;
b2 := 1;
b3 := -3;

Nq :=  Yu.H[Yu] + Yd.H[Yd] - U(g3^2 8/3 + g2^2 3/2 + g1^2 1/30);
Nu := 2 H[Yu].Yu           - U(g3^2 8/3            + g1^2 8/15);
Nd := 2 H[Yd].Yd           - U(g3^2 8/3            + g1^2 2/15);
Nl :=  Ye.H[Ye] + Yn.H[Yn] - U(           g2^2 3/2 + g1^2 3/10);
Ne := 2 H[Ye].Ye           - U(                      g1^2 6/5 );
Nn := 2 H[Yn].Yn;
NHu:= 3 TR[H[Yu].Yu] + TR[H[Yn].Yn] - 	 (g2^2 3/2 + g1^2 3/10);
NHd:= 3 TR[H[Yd].Yd] + TR[H[Ye].Ye] - 	 (g2^2 3/2 + g1^2 3/10);

Pq := U(g3^2 M3 8/3+g2^2 M2 3/2+g1^2 M1  /30) + TAu.H[Yu] + TAd.H[Yd];
Pu := U(g3^2 M3 8/3            +g1^2 M1 8/15) + 2 H[Yu].TAu;
Pd := U(g3^2 M3 8/3	       +g1^2 M1 2/15) + 2 H[Yd].TAd;
Pl := U(            g2^2 M2 3/2+g1^2 M1 3/10) + TAe.H[Ye] + TAn.H[Yn];
Pe := U(                        g1^2 M1 6/5 ) + 2 H[Ye].TAe;
Pn :=                                           2 H[Yn].TAn;
PHu:= g2^2 M2 3/2 + g1^2 M1 3/10 + 3 Tr[H[Yu].TAu] + Tr[H[Yn].TAn];
PHd:= g2^2 M2 3/2 + g1^2 M1 3/10 + 3 Tr[H[Yd].TAd] + Tr[H[Ye].TAe];

S := m2Hu - m2Hd + TR[m2Q] - TR[m2L] - 2 TR[m2U] + TR[m2D] + TR[m2E];

gMs3 := g3^2 norm[M3];
gMs2 := g2^2 norm[M2];
gMs1 := g1^2 norm[M1];
gsS  := g1^2 S;
