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

headername =
    StringReplace[filename, "_dependence.cpp" ~~ EndOfString -> ".hpp"];

WriteString[filename,
  "#include <vector>\n",
  "#include ",InputForm[headername],"\n",
  "\n",
  "\nnamespace flexiblesusy {",
  "\n"
];

writeNBCDeps[filename, {
    M1i, M2i, M3i, A0i,
    _m2Qi, _m2Ui, _m2Di, _m2Li, _m2Ei, _m2Ni,
    _Aui, _Adi, _Ani, _Aei,
    _Yui, _Ydi, _Yni, _Yei
},{
{"Fmssm_constraint_on_ms_n::dependence", Det[m2stop] - (scale0 Exp[t])^4},
{"Fmssm_constraint_on_gauge_couplings_n_::dependence", {g1 - g1i, g2 - g2i, g3 - g3i}},
{"Fmssm_constraint_on_yukawas_n_::dependence", {
    Yu - Table[Yui[i,j], {i,3}, {j,3}],
    Yd - Table[Ydi[i,j], {i,3}, {j,3}],
    Ye - Table[Yei[i,j], {i,3}, {j,3}]}},
{"Fmssm_constraint_on_ewsb_n_::dependence", {
    b vd/vu + (g2^2+gY^2) (vd^4-vu^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hu),
    b vu/vd + (g2^2+gY^2) (vu^4-vd^4)/(4(vd^2+vu^2)) - (norm[mu]+m2Hd),
    Im[mu], Im[b]}}
},{TAu, TAd, TAe}];

WriteString[filename,
  "\n",
  "\n}",
  "\n"
];

Close[filename];
