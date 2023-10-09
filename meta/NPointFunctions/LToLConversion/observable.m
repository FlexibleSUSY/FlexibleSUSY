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

Begin@"Observables`Private`";
Block[{inF, inG, outF, outG, nucl, contr, loopN},
With[{
   obs = FlexibleSUSYObservable`LToLConversion,
   lhs = Sequence[inF_@inG_ -> outF_@outG_, nucl_, contr_, loopN_],
   rhs = Sequence[inF, inG, outF, outG, nucl, contr, loopN],
   cxx = CConversion`ToValidCSymbolString},

AppendTo[FlexibleSUSYObservable`FSObservables, obs];

GetObservableName@obs@lhs := StringTemplate[
   "`1``2`to`3``4`conversion_`5`_`6`_`7`loop",
   InsertionFunction -> cxx]@rhs;

GetObservableDescription@obs@lhs := StringTemplate[
   "`1`(`2`) to `3`(`4`) conversion in `5` for `6` at `7` loop",
   InsertionFunction -> cxx]@rhs;

GetObservableType@obs@lhs := CConversion`ArrayType[CConversion`complexScalarCType, 13];

CalculateObservable[obs@lhs, structName:_String] := StringTemplate[
   "`9`.`1``2`to`3``4`conversion_`5`_`6`_`7`loop = `8`_ltolconversion::" <>
   "calculate_`1``3`_for`6`_`7`loop"<>
   "(`2`, `4`, `8`_ltolconversion::Nucleus::`5`, MODEL, ltolconversion_settings, qedqcd);",
   InsertionFunction -> cxx][rhs, FlexibleSUSY`FSModelName, structName];

]; (* With *)
]; (* Block *)
End[];
