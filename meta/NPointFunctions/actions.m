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

BeginPackage["NPointFunctions`"];
Begin["NPointFunctions`internal`"];

define[getActions,
   {keepProcesses:`type`keepProcesses, settings:{}} :> {},

   {keepProcesses:`type`keepProcesses, settings:{Rule[_,{{___},{___}}]..}} :>
   Module[{
         positiveRules = settings /. Rule[s:_, {p:_, _}] :> Rule[s, p],
         negativeRules = settings /. Rule[s:_, {_, n:_}] :> Rule[s, n],
         discardProcesses = Complement[settings[[All, 1]], keepProcesses],
         actions
      },
      actions = DeleteDuplicates@Join[
         DeleteDuplicates@Flatten[keepProcesses /. positiveRules, 1],
         DeleteDuplicates@Flatten[discardProcesses /. negativeRules, 1]
      ]
   ]
];

getTruePositions::usage = "
@brief Converts a list with a boolean variables to the list of positions for
       all true entries.
@param list A list of booleans.
@returns A list of integers (of an empty one)."
define[getTruePositions, {list:{(True|False)...}} :> Flatten@Position[list, True]];

Module[{template, func, delete, append, restrict},

define[applyAction,
   {
      {diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet},
      {text:_String, topologyQ:_, {function:UnsameQ, name:_, value:_}}
   } :>
   Module[{
         daPairs = getAmplitudeNumbers[diagrams, topologyQ],
         amplitudeNumbers, saveClassRules, viPairs, insertions, res
      },
      amplitudeNumbers = Cases[daPairs, Rule[True, {e:__}] :> e];
      saveClassRules = Table[
         viPairs = getClassRules@amplitudes[[i]];
         insertions = Cases[viPairs, (name -> {e:__}) :> e, Infinity ];
         i -> getTruePositions[function[#, value] &/@ insertions] /. {} -> All,
         {i, amplitudeNumbers}
      ];
      res = {
         deleteClasses[diagrams, daPairs, saveClassRules],
         deleteClasses[amplitudes, daPairs, saveClassRules]
      };
      Print@text;
      printDiagramsInfo@res[[1]];
      res
   ],

   {
      diagrams:`type`diagramSet,
      {topologyQ:_, function:_, crit:_, text:_String}
   } :>
   Module[{
         d = diagrams
      },
      d = If[topologyQ@#[[1]], Part[#,1] -> function[Part[#,2], crit], #] &/@ d;
      d = removeTopologiesWithoutInsertions@d;
      Print@text;
      printDiagramsInfo@d;
      d
   ],

   {
      d:`type`diagramSet,
      {s:_String, t:_, {n:_Integer, f:_}}
   } :> Sequence@@foreach[template[s, t, restrict[f, n]], d],

   {
      d:`type`diagramSet,
      {s:_String, t:_, {Hold, e:_Integer}}
   } :> Sequence@@foreach[template[s, t, delete@e], d],

   {
      d:`type`diagramSet,
      {s:_String, t:_, {Append, (f:`type`field)[n:_Integer] :> e:_Integer}}
   } :> Sequence@@foreach[template[s, t, append[f, n, e]], d]
];

delete[e_] := With[{
      pos = {{2*e}, {2*e-1}}
   },
   Delete[#, pos]&
];

append[field_, number:_, e_] := With[{
      rhs = (First /@ getZeroMassRules[])[[2*e-1]]
   },
   Append[#, genericMass[field, number] :> rhs]&
];

restrict[field_, number_] := {number -> Or[field, -field]};

template[text:_, topologyQ:_, realization:_] := If[topologyQ@getTopology@#2,
   Print@text;
   getTopology@#2 -> foreach[realization&, getInsertions@#2],

   (##&)[]
]&;


];

End[];
EndPackage[];
