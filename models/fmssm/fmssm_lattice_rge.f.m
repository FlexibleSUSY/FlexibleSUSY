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

<< models/fmssm/fmssm_lattice_defs.m

writeRGEs[filename, "Fmssm",
    "a,\n",
    "      double precision a\n", {
{t, 1},
{g1, lsf g1^3 b1},
{g2, lsf g2^3 b2},
{g3, lsf g3^3 b3},
{Yu, a lsf (Nq.Yu + Yu.Nu + NHu Yu)},
{Yd, a lsf (Nq.Yd + Yd.Nd + NHd Yd)},
(* {Yn, a lsf (Nl.Yn + Yn.Nn + NHu Yn)}, *)
{Ye, a lsf (Nl.Ye + Ye.Ne + NHd Ye)},
{mu, a lsf (NHu + NHd) mu},
{b, a lsf ((NHu + NHd) b + 2 (PHu + PHd) mu)},
{M1, a 2 lsf g1^2 b1 M1},
{M2, a 2 lsf g2^2 b2 M2},
{M3, a 2 lsf g3^2 b3 M3},
{m2Hu, a 2 lsf (
       -2 (gMs2 3/2 + gMs1 3/10 - gsS 3/10)
       +3 (TR[Yu.m2Q.H[Yu]]+TR[Yu.m2U.H[Yu]]+m2Hu TR[Yu.H[Yu]]+TR[TAu.H[TAu]])
       +   TR[Yn.m2L.H[Yn]]+TR[Yn.m2N.H[Yn]]+m2Hu TR[Yn.H[Yn]]+TR[TAn.H[TAn]]
       )},
{m2Hd, a 2 lsf (
       -2 (gMs2 3/2 + gMs1 3/10 + gsS 3/10)
       +3 (TR[Yd.m2Q.H[Yd]]+TR[Yd.m2D.H[Yd]]+m2Hd TR[Yd.H[Yd]]+TR[TAd.H[TAd]])
       +   TR[Ye.m2L.H[Ye]]+TR[Ye.m2E.H[Ye]]+m2Hd TR[Ye.H[Ye]]+TR[TAe.H[TAe]]
       )},
{m2Q, a 2 lsf (
      -2 U(gMs3 8/3 + gMs2 3/2 + gMs1/30 - gsS/10)
      +Yu.H[Yu].m2Q/2+m2Q.Yu.H[Yu]/2+Yu.m2U.H[Yu]+m2Hu Yu.H[Yu]+TAu.H[TAu]
      +Yd.H[Yd].m2Q/2+m2Q.Yd.H[Yd]/2+Yd.m2D.H[Yd]+m2Hd Yd.H[Yd]+TAd.H[TAd]
      )},
{m2U, a 2 lsf (
      -2 U(gMs3 8/3 + gMs1 8/15 + gsS 2/5)
      +2 (H[Yu].Yu.m2U/2+m2U.H[Yu].Yu/2+H[Yu].m2Q.Yu+m2Hu H[Yu].Yu+H[TAu].TAu)
      )},
{m2D, a 2 lsf (
      -2 U(gMs3 8/3 + gMs1 2/15 - gsS/5)
      +2 (H[Yd].Yd.m2D/2+m2D.H[Yd].Yd/2+H[Yd].m2Q.Yd+m2Hd H[Yd].Yd+H[TAd].TAd)
      )},
{m2L, a 2 lsf (
      -2 U(gMs2 3/2 + gMs1 3/10 + gsS 3/10)
      +Ye.H[Ye].m2L/2+m2L.Ye.H[Ye]/2+Ye.m2E.H[Ye]+m2Hd Ye.H[Ye]+TAe.H[TAe]
      +Yn.H[Yn].m2L/2+m2L.Yn.H[Yn]/2+Yn.m2N.H[Yn]+m2Hu Yn.H[Yn]+TAn.H[TAn]
      )},
{m2E, a 2 lsf (
      -2 U(gMs1 6/5 - gsS 3/5)
      +2 (H[Ye].Ye.m2E/2+m2E.H[Ye].Ye/2+H[Ye].m2L.Ye+m2Hd H[Ye].Ye+H[TAe].TAe)
      )},
(* {m2N, a 2 lsf (
      2 (H[Yn].Yn.m2N/2+m2N.H[Yn].Yn/2+H[Yn].m2L.Yn+m2Hu H[Yn].Yn+H[TAn].TAn)
      )}, *)
{TAu, a lsf (Nq.TAu + TAu.Nu + NHu TAu + 2 Pq.Yu + 2 Yu.Pu + 2 PHu Yu)},
{TAd, a lsf (Nq.TAd + TAd.Nd + NHd TAd + 2 Pq.Yd + 2 Yd.Pd + 2 PHd Yd)},
(* {TAn, a lsf (Nl.TAn + TAn.Nn + NHu TAn + 2 Pl.Yn + 2 Yn.Pn + 2 PHu Yn)}, *)
{TAe, a lsf (Nl.TAe + TAe.Ne + NHd TAe + 2 Pl.Ye + 2 Ye.Pe + 2 PHd Ye)}
}];

Close[filename];
