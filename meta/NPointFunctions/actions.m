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

getActions[keepProcesses:`type`keepProcesses, settings:{}] := {};
getActions[keepProcesses:`type`keepProcesses, settings:{Rule[_,{{___},{___}}]..}] :=
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
];
getActions // Utils`MakeUnknownInputDefinition;
getActions ~ SetAttributes ~ {Protected,Locked};

applyAction[
   {diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet},
   {topologyQ:_, {name:_, function:_, value:_}, text:_String}
] :=
Module[
   {
      daPairs = getAmplitudeRules[diagrams, topologyQ],
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
];
applyAction[
   diagrams:`type`diagramSet,
   {topologyQ:_, function:_, crit:_, text:_String}
] :=
Module[{
      d = diagrams
   },
   d = If[topologyQ@#[[1]], Part[#,1] -> function[Part[#,2], crit], #] &/@ d;
   d = removeTopologiesWithoutInsertions@d;
   Print@text;
   printDiagramsInfo@d;
   d
];
applyAction[
   diagrams:`type`diagramSet,
   {text:_String, topologyQ:_, {n:_Integer, f:_}}
] :=
Module[{
   },
   Do[
      If[topologyQ@getTopology@d,
         Print@text;
         Return[
            getTopology@d -> Table[{n -> Or[f, -f]}, {Length@getInsertions@d}]
         ]
      ];,
      {d, List@@diagrams}
   ]
];

applyAction[
   diagrams:`type`diagramSet,
   {
      text:_String,
      topologyQ:_,
      {Append, RuleDelayed[(t:`type`fa`field)[n:_Integer], e:_Integer]}
   }
] :=
Module[{
      ext = First /@ `get`zeroMassRules[]
   },
   Do[
      If[topologyQ@getTopology@d,
         Print@text;
         Return[
            getTopology@d -> Table[Append[#, genericMass[t, n] :> ext[[2*e-1]]]&,
               {Length@getInsertions@d}
            ]
         ]
      ];,
      {d, List@@diagrams}
   ]
];

applyAction[
   diagrams:`type`diagramSet,
   {
      text:_String,
      topologyQ:_,
      {Hold, e:_Integer}
   }
] :=
Module[{
   },
   Do[
      If[topologyQ@getTopology@d,
         Print@text;
         Return[
            getTopology@d -> Table[Delete[#, {{2*e}, {2*e-1}}]&,
               {Length@getInsertions@d}
            ]
         ]
      ];,
      {d, List@@diagrams}
   ]
];
applyAction // Utils`MakeUnknownInputDefinition;
applyAction ~ SetAttributes ~ {Protected, Locked};

End[];
EndPackage[];
