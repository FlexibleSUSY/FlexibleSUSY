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

<< models/fmssmn/fmssmn_lattice_defs.m
<< meta/writeNRGE.m

gY := Sqrt[3/5] g1;

m2stop := {
   {m2Q[[3,3]] + vu^2 norm[Yu[[3,3]]] + 1/12 (3 g2^2 - gY^2) (vd^2-vu^2),
    vu Conjugate[TAu[[3,3]]] - vd Conjugate[Yu[[3,3]]] mu},
   {vu TAu[[3,3]] - vd Yu[[3,3]] Conjugate[mu],
    m2U[[3,3]] + vu^2 norm[Yu[[3,3]]] + 1/3 gY^2 (vd^2-vu^2)}};

hTable[mat_, nr_] := Table[Which[
    i===j, Re[mat[i,j]],
    i < j, mat[i,j],
    i > j, Conjugate[mat[j,i]]], {i,nr}, {j,nr}];

hTableT[mat_, nr_] := Transpose[hTable[mat, nr]];

writeNBCs[filename,
"g1i,g2i,g3i,\n" <>
"     &Yui,Ydi,Yni,Yei,\n" <>
"     &m2Hui,m2Hdi,m2Qi,m2Ui,m2Di,m2Li,m2Ni,m2Ei,\n" <>
"     &Aui,Adi,Ani,Aei,\n" <>
"     &M1i,M2i,M3i,\n" <>
"     &vu,vd,\n",
"      double precision g1i,g2i,g3i\n" <>
"      double complex   Yui(3,3),Ydi(3,3),Yni(3,3),Yei(3,3)\n" <>
"      double precision m2Hui,m2Hdi\n" <>
"      double complex   m2Qi(3,3),m2Ui(3,3),m2Di(3,3),\n" <>
"     &                 m2Li(3,3),m2Ni(3,3),m2Ei(3,3)\n" <>
"      double complex   Aui(3,3),Adi(3,3),Ani(3,3),Aei(3,3)\n" <>
"      double complex   M1i,M2i,M3i\n" <>
"      double precision vu,vd\n",
{
    Yui, Ydi, Yni, Yei,
    m2Qi, m2Ui, m2Di, m2Li, m2Ni, m2Ei,
    Aui, Adi, Ani, Aei
},{
    M1i, M2i, M3i,
    _m2Qi, _m2Ui, _m2Di, _m2Li, _m2Ei, _m2Ni,
    _Aui, _Adi, _Ani, _Aei,
    _Yui, _Ydi, _Yni, _Yei
},{
{"Fmssmn_mx", g1 - g2},
{"Fmssmn_higgs_masses", {m2Hu - m2Hui, m2Hd - m2Hdi}},
{"Fmssmn_gaugino_masses", {M1 - M1i, M2 - M2i, M3 - M3i}},
{"Fmssmn_sfermion_masses", {
    m2Q - hTable[m2Qi, 3],
    m2U - hTable[m2Ui, 3],
    m2D - hTable[m2Di, 3],
    m2L - hTable[m2Li, 3],
    m2E - hTable[m2Ei, 3],
    m2N - hTable[m2Ni, 3]}},
{"Fmssmn_trilinears", {
    TAu - Yu Table[Aui[i,j], {i,3}, {j,3}],
    TAd - Yd Table[Adi[i,j], {i,3}, {j,3}],
    TAn - Yn Table[Ani[i,j], {i,3}, {j,3}],
    TAe - Ye Table[Aei[i,j], {i,3}, {j,3}]}},
{"Fmssmn_ms", Det[m2stop] - (scale0 Exp[t])^4},
{"Fmssmn_gauge_couplings", {g1 - g1i, g2 - g2i, g3 - g3i}},
{"Fmssmn_yude", {
    Yu - Table[Yui[i,j], {i,3}, {j,3}],
    Yd - Table[Ydi[i,j], {i,3}, {j,3}],
    Ye - Table[Yei[i,j], {i,3}, {j,3}]}},
{"Fmssmn_yn", Yn - Table[Yni[i,j], {i,3}, {j,3}]},
{"Fmssmn_ewsb", {
    b vd/vu + (g2^2+gY^2) (vd^4-vu^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hu),
    b vu/vd + (g2^2+gY^2) (vu^4-vd^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hd),
    Im[mu], Im[b]}}
},{}];

Close[filename];
