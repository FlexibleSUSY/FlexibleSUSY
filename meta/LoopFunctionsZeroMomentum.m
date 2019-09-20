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

BeginPackage["LoopFunctionsZeroMomentum`"]

loopFunctionsZeroMomentum::usage =
"Passarino-Veltman 1-loop functions with squared mass arguments and
vanishing external momenta:

  A0[x,q]
  B0[x,y,q]
  D1B0[x,y,q]
  C0[x,y,z,q]
  D0[x,y,z,u,q],
  E0[x,y,z,u,v,q]
  F0[x,y,z,u,v,w,q]

where

  x, y, z, u, v, w: squared masses
  q: squared renormalization scale
"

sortLoopFunctionsZeroMomentum::usage =
"List of replacement rules to bring arguments of loop functions into a
canonical order."

eps::usage =
"\[Epsilon] = 4 - d, with d being the space-time dimension in
dimensional regularization (DREG) or dimensional reduction (DRED),
respectively."

epsIR::usage = "infrared regulator"

L::usage = "L[a,b] = Log[a/b]"

L[m_, m_] := 0

Derivative[1, 0][L][x_, y_] := 1/x

Derivative[0, 1][L][x_, y_] := -1/y

loopFunctionsZeroMomentum = {
    (* A0 *)
    A0[0, q_] :>
      0,
    A0[x_, q_] :>
      x (2/eps + 1 - L[x, q]),
    (* B0 *)
    B0[0, 0, q_] :>
      2/eps - 2/epsIR,
    B0[x_, x_, q_] :>
      A0[x, q]/x - 1,
    B0[x_, y_, q_] :>
      (A0[x, q] - A0[y, q])/(x - y),
    (* derivative of B0 w.r.t. p^2 *)
    D1B0[x_, y_, q_] :>
      -(C0[x, y, y, q] + x D0[x, x, y, y, q])/2,
    (* C0 *)
    C0[0, 0, 0, q_] ->
      0,
    C0[x_, x_, x_, q_] :>
      -1/(2 x),
    C0[y_, y_, z_, q_] :>
      C0[y, z, y, q],
    C0[x_, y_, z_, q_] :>
      (B0[x, z, q] - B0[y, z, q])/(x - y),
    (* D0 *)
    D0[0, 0, 0, 0, q_] ->
      0,
    D0[x_, x_, x_, x_, q_] :>
      1/(6 x^2),
    D0[y_, y_, z_, u_, q_] :>
      D0[y, z, u, y, q],
    D0[x_, y_, z_, u_, q_] :>
      (C0[x, z, u, q] - C0[y, z, u, q])/(x - y),
    (* E0 *)
    E0[0, 0, 0, 0, 0, q_] :>
      0,
    E0[x_, x_, x_, x_, x_, q_] :>
      -1/(12 x^3),
    E0[y_, y_, z_, u_, v_, q_] :>
      E0[y, z, u, v, y, q],
    E0[x_, y_, z_, u_, v_, q_] :>
      (D0[x, z, u, v, q] - D0[y, z, u, v, q])/(x - y),
    (* F0 *)
    F0[0, 0, 0, 0, 0, 0, q_] :>
      0,
    F0[x_, x_, x_, x_, x_, x_, q_] :>
      1/(20 x^4),
    F0[y_, y_, z_, u_, v_, w_, q_] :>
      F0[y, z, u, v, w, y, q],
    F0[x_, y_, z_, u_, v_, w_, q_] :>
      (E0[x, z, u, v, w, q] - E0[y, z, u, v, w, q])/(x - y)
}

sortLoopFunctionsZeroMomentum = {
    B0[x_, y_, q_] :> B0[Sequence @@ Sort[{x, y}], q],
    D1B0[x_, y_, q_] :> D1B0[Sequence @@ Sort[{x, y}], q],
    C0[x_, y_, z_, q_] :> C0[Sequence @@ Sort[{x, y, z}], q],
    D0[x_, y_, z_, u_, q_] :> D0[Sequence @@ Sort[{x, y, z, u}], q],
    E0[x_, y_, z_, u_, v_, q_] :> E0[Sequence @@ Sort[{x, y, z, u, v}], q],
    F0[x_, y_, z_, u_, v_, w_, q_] :> F0[Sequence @@ Sort[{x, y, z, u, v, w}], q]
}

EndPackage[]
