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

Needs["Utils`"];
Needs["TreeMasses`"];

BeginPackage["Wrappers`"];
Begin["`Private`"];

warnVVVV[fields:{_, _, _, _}, vertex:{__}] :=
   "Input fields "<> ToString@fields<> " lead to the following combination "<>
   "of Lorentz indices in the vertex "<> ToString@First@vertex<> ". "<>
   "They should be equal to {lt1, lt2, lt3, lt4} up to any permutation! "<>
   "Most likely the vertex is wrong. Please, cross-check the input!";

Module[{status = True, res},
   SARAH`Vertex[fields:{Repeated[_?TreeMasses`IsVector, {4}]},
                opts___]/; status :=
   Module[{naked},
      status = False;
      res = SARAH`Vertex[fields, opts];
      naked = Sort@Cases[First@res, {___, lt_}:> lt, Infinity];
      If[Not@MatchQ[naked, {SARAH`lt1, SARAH`lt2, SARAH`lt3, SARAH`lt4}],
         Utils`FSFancyWarning@warnVVVV[fields, res];
      ];
      status = True;
      res
   ];
];

End[];
Block[{$ContextPath}, EndPackage[]];
