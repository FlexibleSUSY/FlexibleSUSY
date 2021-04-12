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
topology = FeynArts`Topology[_Integer][propagator..];

diagram = Rule[topology, FeynArts`Insertions[Generic][__]];
diagramSet = FeynArts`TopologyList[_][diagram..];

amplitude = FeynArts`FeynAmp[
   FeynArts`GraphID[FeynArts`Topology==_Integer,Generic==_Integer],
   Integral[FeynArts`FourMomentum[FeynArts`Internal,_Integer]],
   _,
   {__}->FeynArts`Insertions[FeynArts`Classes][{__}..]];
amplitudeSet = FeynArts`FeynAmpList[__][amplitude..];

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
genericMass =
   FeynArts`Mass[field@genericIndex, mass];

`fc`particle = field[_Integer, Repeated[{_Symbol}, {0, 1}]];
`fc`mass = 0|_Symbol|_Symbol@_Symbol;
`fc`external = {`fc`particle|-`fc`particle,
   FormCalc`k@_Integer, `fc`mass, {}};
`fc`process = {`fc`external..} -> {`fc`external..};
`fc`amplitude = FormCalc`Amp[`fc`process][_];
`fc`amplitudeSet = {`fc`amplitude..};

head =
   {_FeynArts`TopologyList, ___};
generic =
   FeynArts`Insertions[Generic]@__;
classes =
   FeynArts`Insertions[FeynArts`Classes]@__;
With[{node = NPointFunctions`Private`node},
   tree =
      node[head,
         node[topology,
            node[generic,
               node[classes]..]..]..];];

EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "type`"];
