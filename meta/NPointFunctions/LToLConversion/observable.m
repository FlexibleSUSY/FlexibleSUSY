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

Begin["Observables`Private`"];

`args`LToLConversion = Sequence[
   (lIn:_?TreeMasses`IsLepton)[gIn:_Integer] ->
   (lOut_?TreeMasses`IsLepton)[gOut:_Integer],
   nucleus:_,
   contribution:Alternatives[All,
      NPointFunctions`noScalars,
      NPointFunctions`Penguins,
      NPointFunctions`FourFermionScalarPenguins,
      NPointFunctions`FourFermionMassiveVectorPenguins,
      NPointFunctions`FourFermionFlavourChangingBoxes]
];

GetObservableName@FlexibleSUSYObservable`LToLConversion@`args`LToLConversion :=
   StringJoin[
      #@lIn, #@gIn, "_to_", #@lOut, #@gOut, "_conversion_in_", #@nucleus,
      "_for_", #@contribution
   ] &@ CConversion`ToValidCSymbolString;

GetObservableDescription@FlexibleSUSYObservable`LToLConversion@
`args`LToLConversion :=
   StringJoin[
      #@lIn, "(", #@gIn, ") to ", #@lOut, "(", #@gOut, ") conversion in ",
      #@nucleus, " for ",#@contribution
   ] &@ CConversion`ToValidCSymbolString;

GetObservableType@FlexibleSUSYObservable`LToLConversion@`args`LToLConversion :=
   CConversion`ArrayType[CConversion`complexScalarCType, 13];

CalculateObservable[
   FlexibleSUSYObservable`LToLConversion@`args`LToLConversion,
   structName:_String
] :=
   StringJoin[
      structName, ".",
      GetObservableName@FlexibleSUSYObservable`LToLConversion[
         lIn@gIn->lOut@gOut, nucleus, contribution],
      " = ", #2, "::calculate_", #1@lIn, "_to_", #1@lOut, "_for_",
      #1@contribution, "(", #1@gIn, ", ", #1@gOut, ", ", #2, "::Nucleus::",
      #1@nucleus, ", MODEL, qedqcd);"
   ] & [CConversion`ToValidCSymbolString,
      FlexibleSUSY`FSModelName<>"_l_to_l_conversion"
   ];

End[];
