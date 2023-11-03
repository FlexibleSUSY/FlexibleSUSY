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

BeginPackage@"NPointFunctions`";
Begin@"`Private`";

IsVertex[e_] := MatchQ[e, FeynArts`Vertex[_Integer][_Integer]];
IsMass[e_] := MatchQ[e, Repeated[FeynArts`Loop|FeynArts`Internal, {0, 1}]];
IsPropagator[e_] := MatchQ[e, FeynArts`Propagator[ FeynArts`External|FeynArts`Incoming|
   FeynArts`Outgoing|FeynArts`Internal|FeynArts`Loop[_Integer] ][
      _?IsVertex, _?IsVertex, Repeated[FeynArts`Field[_Integer],{0,1}]]
];
IsAmplitude[e_] := MatchQ[e,
   FeynArts`FeynAmp[
      FeynArts`GraphID[FeynArts`Topology==_Integer,Generic==_Integer],
      (* For 1loop. *)Integral[FeynArts`FourMomentum[FeynArts`Internal,_Integer]]|
      (* For 0loop. *)Integral[],
      _,
      {__}->FeynArts`Insertions[FeynArts`Classes][{__}..]
   ]
];

IsTopology[e_] := MatchQ[e, FeynArts`Topology[_Integer][__?IsPropagator]];
IsDiagram[e_] := MatchQ[e, Rule[_?IsTopology, FeynArts`Insertions[Generic][__]]];
IsGeneric[e_] := MatchQ[e, FeynArts`Insertions[Generic][__]];
IsClasses[e_] := MatchQ[e, FeynArts`Insertions[FeynArts`Classes][__]];
IsTreeHead[e_] := MatchQ[e, {_FeynArts`TopologyList, ___}];
IsTree[e_] := MatchQ[e,
   node[_?IsTreeHead, node[_?IsTopology, node[_?IsGeneric, node[_?IsClasses]..]..]..]
];

IsTopologyListHead[e_] := MatchQ[e, FeynArts`TopologyList[_]];
IsDiagramSet[e_] := MatchQ[e, FeynArts`TopologyList[_][__?IsDiagram]];
IsAmplitudeSet[e_] := MatchQ[e, FeynArts`FeynAmpList[__][__?IsAmplitude]];

IsColorIndex[e_] := MatchQ[e, FeynArts`Index[Global`Colour, _Integer]];
IsGluonIndex[e_] := MatchQ[e, FeynArts`Index[Global`Gluon, _Integer]];
IsGenerationIndex[e_] := !IsColorIndex[e] && !IsGluonIndex[e] && MatchQ[e, FeynArts`Index[_Symbol, _Integer]];

IsFeynArtsField[e_] := MatchQ[e, FeynArts`S|FeynArts`F|FeynArts`V|FeynArts`U];
IsGenericField[e_] := MatchQ[e, (_?IsFeynArtsField)[FeynArts`Index[Generic, _Integer]]];

IsFormCalcField[e_] := MatchQ[e, (_?IsFeynArtsField)[_Integer, Repeated[{_Symbol}, {0, 1}]]];
IsFormCalcMass[e_] := MatchQ[e, 0|_Symbol|_Symbol@_Symbol];
IsFormCalcExternal[e_] := MatchQ[e, {_?IsFormCalcField|(-1*_?IsFormCalcField), FormCalc`k@_Integer, _?IsFormCalcMass, {}}];
IsFormCalcProcess[e_] := MatchQ[e, {__?IsFormCalcExternal} -> {__?IsFormCalcExternal}];
IsFormCalcAmplitude[e_] := MatchQ[e, FormCalc`Amp[_?IsFormCalcProcess][_]];
IsFormCalcSet[e_] := MatchQ[e, {__?IsFormCalcAmplitude}];

End[];
EndPackage[];
