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

(* loop functions from Appendix A of arxiv:1407.4081 *)

F1[0] := 0;
F1[1] := 1;
F1[x_] := x Log[x^2] / (x^2 - 1);

F2[0] := 0;
F2[1] := 1;
F2[x_] := 6 x^2 (2 - 2 x^2 + (1 + x^2) Log[x^2]) / (x^2 - 1)^3;

F3[0] := 0;
F3[1] := 1;
F3[x_] := 2 x (5 (1 - x^2) + (1 + 4 x^2) Log[x^2]) / (3 (x^2 - 1)^2);

F4[0] := 0;
F4[1] := 1;
F4[x_] := 2 x (x^2 - 1 - Log[x^2]) / (x^2 - 1)^2;

F5[0] := 0;
F5[1] := 1;
F5[x_] := 3 x (1 - x^4 + 2 x^2 Log[x^2]) / (1 - x^2)^3;

F6[0] := -3/4;
F6[1] := 0;
F6[x_] := (x^2 - 3)/(4 (1 - x^2)) + x^2 (x^2 - 2) Log[x^2] / (2 (1 - x^2)^2);

F7[0] := -3/2;
F7[1] := 1;
F7[x_] := -3 (x^4 - 6 x^2 + 1)/(2 (x^2 - 1)^2) + 3 x^4 (x^2 - 3) Log[x^2] / (x^2 - 1)^3;

F8[0, 0] := -2;
F8[0, 1] := 0;
F8[1, 0] := 0;
F8[1, 1] := 1;
F8[0, y_] := -2 + (2*y^2*Log[y^2])/(-1 + y^2);
F8[x_, 0] := F8[0, x];
F8[1, y_] := -2 + (2*(1 - y^2 + y^4*Log[y^2]))/(-1 + y^2)^2;
F8[x_, 1] := F8[1, x];
F8[x_, x_] := (2*(-1 + x^2 - 2*x^2*Log[x^2] + x^4*Log[x^2]))/(-1 + x^2)^2;
F8[x1_, x2_] := -2 + 2/(x1^2 - x2^2) (x1^4 Log[x1^2]/(x1^2 - 1) - x2^4 Log[x2^2]/(x2^2 - 1));

F9[0, 1] := 2;
F9[1, 0] := 2;
F9[1, 1] := 1;
F9[0, y_] := (2*Log[y^2])/(-1 + y^2);
F9[x_, 0] := F9[0, x];
F9[1, y_] := (2*(1 - y^2 + y^2*Log[y^2]))/(-1 + y^2)^2;
F9[x_, 1] := F9[1, x];
F9[x_, x_] := (2*(-1 + x^2 - Log[x^2]))/(-1 + x^2)^2;
F9[x1_, x2_] := 2/(x1^2 - x2^2) (x1^2 Log[x1^2]/(x1^2 - 1) - x2^2 Log[x2^2]/(x2^2 - 1));

f[r_] := F5[r];

g[r_] := F7[r];

f1[0] := 0;
f1[1] := 1;
f1[r_] := 6 (r^2 + 3) r^2 / (7 (r^2 - 1)^2) + 6 (r^2 - 5) r^4 Log[r^2] / (7 (r^2 - 1)^3);

f2[0] := 0;
f2[1] := 1;
f2[r_] := 2 (r^2 + 11) r^2 / (9 (r^2 - 1)^2) + 2 (5 r^2 - 17) r^4 Log[r^2] / (9 (r^2 - 1)^3);

f3[0] := 4/3;
f3[1] := 1;
f3[r_] := 2 (r^4 + 9 r^2 + 2) / (3 (r^2 - 1)^2) + 2 (r^4 - 7 r^2 - 6) r^2 Log[r^2] / (3 (r^2 - 1)^3);

f4[0] := 12/7;
f4[1] := 1;
f4[r_] := 2 (5 r^4 + 25 r^2 + 6) / (7 (r^2 - 1)^2) + 2 (r^4 - 19 r^2 - 18) r^2 Log[r^2] / (7 (r^2 - 1)^3);

f5[0, 0] := 3/4;
f5[0, 1] := 3/4;
f5[1, 0] := 3/4;
f5[1, 1] := 1;
f5[0, y_] := (3*(1 + y^2)*(1 - y^2 + y^2*Log[y^2]))/(4*(-1 + y^2)^2);
f5[x_, 0] := f5[0, x];
f5[1, y_] := (3*(-1 + y + 2*y^2 - y^4 - y^5 + (y^3 + y^5)*Log[y^2]))/(4*(-1 + y)^3*(1 + y)^2);
f5[x_, 1] := f5[1, x];
f5[x_, x_] := (3*(-1 - 5*x^2 + 5*x^4 + x^6 + x^2*(-3 - 6*x^2 + x^4)*Log[x^2]))/(4*(-1 + x^2)^3);
f5[r1_, r2_] := (3/4) (
    (1 + (r1 + r2)^2 - r1^2 r2^2)/((r1^2 - 1) (r2^2 - 1))
    + (r1^3 (r1^2 + 1) Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - (r2^3 (r2^2 + 1) Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

f6[0, 0] := 0;
f6[0, 1] := 3/7;
f6[1, 0] := 3/7;
f6[1, 1] := 1;
f6[0, y_] := (6*(y^2 - y^4 + y^4*Log[y^2]))/(7*(-1 + y^2)^2);
f6[x_, 0] := f6[0, x];
f6[1, y_] := (3*(-1 + 2*y^2 + 2*y^3 - y^4 - 2*y^5 + 2*y^5*Log[y^2]))/(7*(-1 + y)^3*(1 + y)^2);
f6[x_, 1] := f6[1, x];
f6[x_, x_] := (6*x^2*(-3 + 2*x^2 + x^4 + x^2*(-5 + x^2)*Log[x^2]))/(7*(-1 + x^2)^3);
f6[r1_, r2_] := (6/7) (
    (r1^2 + r2^2 + r1 r2 - r1^2 r2^2)/((r1^2 - 1) (r2^2 - 1))
    + ((r1)^5 Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - ((r2)^5 Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

f7[0, 0] := 6;
f7[0, 1] := 3;
f7[1, 0] := 3;
f7[1, 1] := 1;
f7[0, y_] := (6*(1 - y^2 + y^2*Log[y^2]))/(-1 + y^2)^2;
f7[x_, 0] := f7[0, x];
f7[1, y_] := (-3*(1 - 2*y - 2*y^2 + 2*y^3 + y^4 - 2*y^3*Log[y^2]))/((-1 + y)^3*(1 + y)^2);
f7[x_, 1] := f7[1, x];
f7[x_, x_] := (-6*(1 + 2*x^2 - 3*x^4 + x^2*(3 + x^2)*Log[x^2]))/(-1 + x^2)^3;
f7[r1_, r2_] := 6 (
    (1 + r1 r2)/((r1^2 - 1) (r2^2 - 1))
    + ((r1)^3 Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - ((r2)^3 Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

f8[0, 0] := 0;
f8[0, 1] := 3/4;
f8[1, 0] := 3/4;
f8[1, 1] := 1;
f8[0, y_] := (3*(y - y^3 + y^3*Log[y^2]))/(2*(-1 + y^2)^2);
f8[x_, 0] := f8[0, x];
f8[1, y_] := (3*(-1 + 4*y^2 - 3*y^4 + 2*y^4*Log[y^2]))/(4*(-1 + y)^3*(1 + y)^2);
f8[x_, 1] := f8[1, x];
f8[x_, x_] := (3*x*(-1 + x^4 - 2*x^2*Log[x^2]))/(-1 + x^2)^3;
f8[r1_, r2_] := (3/2) (
    (r1 + r2)/((r1^2 - 1) (r2^2 - 1))
    + ((r1)^4 Log[r1^2])/((r1^2 - 1)^2 (r1 - r2))
    - ((r2)^4 Log[r2^2])/((r1 - r2) (r2^2 - 1)^2)
);

(* Delta and Phi function for 2-loop threshold corrections *)

Cl2[x_] := Im[Li2[Exp[I x]]];

Li2[x_] := PolyLog[2, x];

Delta[x_, y_, z_] := x^2 + y^2 + z^2 - 2 (x y + x z + y z);

Lambda2[u_, v_] := (1 - u - v)^2 - 4 u v;

PhiPos[u_, v_] :=
    Module[{lam = Sqrt[Lambda2[u, v]]},
           1/lam (
               - Log[u] Log[v]
               + 2 Log[(1 - lam + u - v)/2] Log[(1 - lam - u + v)/2]
               - 2 Li2[(1 - lam + u - v)/2]
               - 2 Li2[(1 - lam - u + v)/2]
               + Pi^2/3
           )
          ];

PhiNeg[u_, v_] :=
    Module[{lam = Sqrt[-Lambda2[u, v]]},
           2/lam (
               + Cl2[2 ArcCos[(1 + u - v)/(2 Sqrt[u])]]
               + Cl2[2 ArcCos[(1 - u + v)/(2 Sqrt[v])]]
               + Cl2[2 ArcCos[(-1 + u + v)/(2 Sqrt[u v])]]
           )
          ];

PhiUV[u_, v_] :=
    Module[{lam = Lambda2[u, v]},
           If[lam > 0,
              Which[
                  u <= 1 && v <= 1, PhiPos[u, v],
                  u >= 1 && v/u <= 1, PhiPos[1/u, v/u]/u,
                  True, PhiPos[1/v, u/v]/v
                   ],
              PhiNeg[u, v]
             ]
          ];

(* Phi[x,y,z] function from
 * Davydychev and Tausk, Nucl. Phys. B397 (1993) 23.
 * Arguments are interpreted as squared mass parameters.
 *)
Phi[x_, y_, z_] := PhiUV[x/z, y/z];
