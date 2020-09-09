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

Module[{
      template, func
   },

template[{text:_, topologyQ:_}, realization:_] := If[topologyQ@getTopology@#2,
   Print@text;
   getTopology@#2 -> Table[realization, {Length@getInsertions@#2}],

   Null
]&;

makeApplication::usage = "
@brief This metafunction allows to define applyAction functions without the
       repetition of the code. One need to specify the pattern of setting to
       be used with applyAction and the realization of this setting. Patterns
       for names text and topologyQ are reserved for information text to be
       printed and the topologyQ-function to parse the correct topology.
@param setting A pattern, replesenting some setting.
@param realization An expression, which represent something, which should be
       done for a setting.
@returns Null.";
makeApplication[{setting:_, realization:_}] := (
   applyAction[diagrams:`type`diagramSet, setting] := (
      func = template[{text, topologyQ}, realization];
      Sequence@@(foreach[func, List@@diagrams] /. Null -> Sequence[])
   );
);
makeApplication // Utils`MakeUnknownInputDefinition;
makeApplication ~ SetAttributes ~ {Protected, Locked};

];

makeApplication@{
   {
      text:_String,
      topologyQ:_,
      {n:_Integer, f:_}
   },
   {n -> Or[f, -f]}
};

makeApplication@{
   {
      text:_String,
      topologyQ:_,
      {Hold, e:_Integer}
   },
   Delete[#, {{2*e}, {2*e-1}}]&
};

makeApplication@{
   {
      text:_String,
      topologyQ:_,
      {Append, (t:`type`fa`field)[n:_Integer] :> e:_Integer}
   },
   Append[#, genericMass[t, n] :> (First /@ `get`zeroMassRules[])[[2*e-1]]]&
};

applyAction // Utils`MakeUnknownInputDefinition;
applyAction ~ SetAttributes ~ {Protected, Locked};

End[];
EndPackage[];
