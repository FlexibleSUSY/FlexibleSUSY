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

Module[{status = True, res},
   SARAH`Vertex[fields:{Repeated[_?TreeMasses`IsVector, {4}]},
                opts___]/; status :=
   Module[{naked},
      status = False;
      res = SARAH`Vertex[fields, opts];
      naked = Sort@Cases[First@res, {___, lt_}:> lt, Infinity];
      If[Not@MatchQ[naked, {SARAH`lt1, SARAH`lt2, SARAH`lt3, SARAH`lt4}],
         Utils`FSFancyWarning["Input fields are "<>ToString@fields<>
            ". It leads to the combination of Lorentz indices in the vertex "<>
            ToString@First@res<>", which is not equal to {lt1, lt2, lt3, lt4}"<>
            " up to a permutation. Most likely, it will give an incorrect form"<>
            " of the vertex. Please, cross-check the input!"
         ];
      ];
      status = True;
      res
   ];
];

End[];
Block[{$ContextPath}, EndPackage[]];
