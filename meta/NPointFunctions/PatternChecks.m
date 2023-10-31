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

BeginPackage@"type`";

vertex = FeynArts`Vertex[_Integer][_Integer];

propagator = FeynArts`Propagator[ FeynArts`External|FeynArts`Incoming|
   FeynArts`Outgoing|FeynArts`Internal|FeynArts`Loop[_Integer] ][
      vertex, vertex, Repeated[FeynArts`Field[_Integer],{0,1}]];



amplitude = FeynArts`FeynAmp[
   FeynArts`GraphID[FeynArts`Topology==_Integer,Generic==_Integer],
   (* For 1loop. *)Integral[FeynArts`FourMomentum[FeynArts`Internal,_Integer]]|
   (* For 0loop. *)Integral[],
   _,
   {__}->FeynArts`Insertions[FeynArts`Classes][{__}..]];

colorIndex = FeynArts`Index[Global`Colour, _Integer];
gluonIndex = FeynArts`Index[Global`Gluon, _Integer];
generationIndex =
   FeynArts`Index[
      `Private`s_Symbol /; !MatchQ[`Private`s, Global`Gluon|Global`Colour],
      _Integer];
genericIndex = FeynArts`Index[Generic, _Integer];

field = FeynArts`S|FeynArts`F|FeynArts`V|FeynArts`U;
genericField = field@genericIndex;
mass = Repeated[FeynArts`Loop|FeynArts`Internal, {0, 1}];

`fc`particle = field[_Integer, Repeated[{_Symbol}, {0, 1}]];
`fc`mass = 0|_Symbol|_Symbol@_Symbol;
`fc`external = {`fc`particle|-`fc`particle,
   FormCalc`k@_Integer, `fc`mass, {}};
`fc`process = {`fc`external..} -> {`fc`external..};
`fc`amplitude = FormCalc`Amp[`fc`process][_];

head = {_FeynArts`TopologyList, ___};
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "type`"];

Begin@"NPointFunctions`Private`";

IsTopology[e_] := MatchQ[e, FeynArts`Topology[_Integer][type`propagator..]];
IsDiagram[e_] := MatchQ[e, Rule[_?IsTopology, FeynArts`Insertions[Generic][__]]];
IsGeneric[e_] := MatchQ[e, FeynArts`Insertions[Generic][__]];
IsClasses[e_] := MatchQ[e, FeynArts`Insertions[FeynArts`Classes][__]];
IsTree[e_] := MatchQ[e,
   node[type`head, node[_?IsTopology, node[_?IsGeneric, node[_?IsClasses]..]..]..]
];

IsTopologyListHead[e_] := MatchQ[e, FeynArts`TopologyList[_]];
IsDiagramSet[e_] := MatchQ[e, FeynArts`TopologyList[_][__?IsDiagram]];
IsAmplitudeSet[e_] := MatchQ[e, FeynArts`FeynAmpList[__][type`amplitude..]];
IsFormCalcSet[e_] := MatchQ[e, {type`fc`amplitude..}];

End[];

