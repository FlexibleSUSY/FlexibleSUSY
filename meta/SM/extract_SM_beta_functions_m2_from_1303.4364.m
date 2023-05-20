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

(*
   This script extracts the SM beta-functions from 1303.4364v2.

   Author: Alexander Voigt

   Run it like this:

   math -run "<< extract_SM_beta_functions_m2_from_1303.4364.m; Quit[]"
*)

Get["beta_m2_1303.4364.m"];

rules = {
    al1 -> g1^2,
    al2 -> g2^2,
    al3 -> g3^2,
    at -> gt^2,
    ab -> gb^2,
    atau -> g\[Tau]^2,
    lam -> \[Lambda]/4,
    NR -> 3,
    cR -> 4/3,
    NG -> 3,
    z3 -> N[Zeta[3]],
    TF -> 1/2,
    cA -> 3
};

beta = {
    Simplify[2 m2 Coefficient[bms, h, 1] /. rules],
    Simplify[2 m2 Coefficient[bms, h, 2] /. rules],
    Simplify[2 m2 Coefficient[bms, h, 3] /. rules]
};

outputDir = ".";
filename = FileNameJoin[{outputDir, "beta_m2.m"}];

Put[beta, filename];
