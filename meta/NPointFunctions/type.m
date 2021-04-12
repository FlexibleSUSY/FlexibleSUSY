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

`type`vertex = FeynArts`Vertex[_Integer][_Integer];
`type`propagator = FeynArts`Propagator[ FeynArts`External|FeynArts`Incoming|
   FeynArts`Outgoing|FeynArts`Internal|FeynArts`Loop[_Integer] ][
      `type`vertex, `type`vertex, Repeated[FeynArts`Field[_Integer],{0,1}]];
`type`topology = FeynArts`Topology[_Integer][`type`propagator..];

`type`diagram = Rule[`type`topology, FeynArts`Insertions[Generic][__]];
`type`diagramSet = FeynArts`TopologyList[_][`type`diagram..];

`type`amplitude = FeynArts`FeynAmp[
   FeynArts`GraphID[FeynArts`Topology==_Integer,Generic==_Integer],
   Integral[FeynArts`FourMomentum[FeynArts`Internal,_Integer]],
   _,
   {__}->FeynArts`Insertions[FeynArts`Classes][{__}..]];
`type`amplitudeSet = FeynArts`FeynAmpList[__][`type`amplitude..];

`type`colorIndex = FeynArts`Index[Global`Colour, _Integer];
`type`gluonIndex = FeynArts`Index[Global`Gluon, _Integer];
`type`generationIndex =
   FeynArts`Index[
      s_Symbol /; !MatchQ[s, Global`Gluon|Global`Colour],
      _Integer];
`type`genericIndex = FeynArts`Index[Generic, _Integer];

`type`field = FeynArts`S|FeynArts`F|FeynArts`V|FeynArts`U;
`type`genericField = `type`field[`type`genericIndex];
`type`mass = Repeated[FeynArts`Loop|FeynArts`Internal, {0, 1}];
`type`genericMass =
   FeynArts`Mass[`type`field[`type`genericIndex], `type`mass];

`type`fc`particle = `type`field[_Integer, Repeated[{_Symbol}, {0, 1}]];
`type`fc`mass = 0|_Symbol|_Symbol@_Symbol;
`type`fc`external = {`type`fc`particle|-`type`fc`particle,
   FormCalc`k@_Integer, `type`fc`mass, {}};
`type`fc`process = {`type`fc`external..} -> {`type`fc`external..};
`type`fc`amplitude = FormCalc`Amp[`type`fc`process][_];
`type`fc`amplitudeSet = {`type`fc`amplitude..};

`type`head =
   {_FeynArts`TopologyList, ___};
`type`generic =
   FeynArts`Insertions[Generic]@__;
`type`classes =
   FeynArts`Insertions[FeynArts`Classes]@__;
`type`tree =
   node[`type`head,
      node[`type`topology,
         node[`type`generic,
            node[`type`classes]..]..]..];

End[];
EndPackage[];
