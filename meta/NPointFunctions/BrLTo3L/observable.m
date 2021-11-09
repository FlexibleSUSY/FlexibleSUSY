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

Utils`DynamicInclude@"type.m";

Begin@"Observables`Private`";
With[{args = BrLTo3L`arguments[lep, nI -> {nO, nA}, proc, loopN],
      obs = FlexibleSUSYObservable`BrLTo3L,
      cxx = CConversion`ToValidCSymbolString},

   AppendTo[FlexibleSUSYObservable`FSObservables, obs];

   GetObservableName@obs@args := StringJoin[
      cxx@lep, cxx@nI, "_to_", cxx@lep, cxx@nO, cxx@lep, cxx@nA,
      cxx@SARAH`bar@lep, cxx@nA, "_for_", SymbolName@proc,
      "_", ToString@loopN, "loop"];

   GetObservableDescription@obs@args := StringJoin[
      cxx@lep, "(", cxx@nI, ") to ", cxx@lep, "(", cxx@nO, ")", cxx@lep, "(",
      cxx@nA, ")", cxx@SARAH`bar@lep, "(", cxx@nA, ")",
      " for ", SymbolName@proc,
      " ", ToString@loopN, " loop"];

   GetObservableType@obs@args :=
      CConversion`ArrayType[CConversion`complexScalarCType, 13];

   CalculateObservable[obs@args, structName:_String] := StringJoin[
      structName, ".", GetObservableName@obs@##, " = ",
      BrLTo3L`calculate@obs@##, ";"]&
         [lep@nI -> {lep@nO, lep@nA, SARAH`bar@lep@nA}, proc, loopN];
];
End[];
