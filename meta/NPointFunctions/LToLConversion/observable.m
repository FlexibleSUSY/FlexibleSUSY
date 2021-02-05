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

Get@FileNameJoin@{DirectoryName@$Input, "type.m"};

Begin["Observables`Private`"];
With[{
   args = LToLConversion`arguments[in@inN, out@outN, nucl, proc, massless],
   obs = FlexibleSUSYObservable`LToLConversion,
   cxx = CConversion`ToValidCSymbolString
},

GetObservableName@obs@args := StringJoin[
   cxx@in, cxx@inN, "_to_", cxx@out, cxx@outN, "_conversion_in_",
   cxx@nucl, "_for_", cxx@proc, cxx@massless];

GetObservableDescription@obs@args := StringJoin[
   cxx@in, "(", cxx@inN, ") to ", cxx@out, "(", cxx@outN, ") conversion in ",
   cxx@nucl, " for ",cxx@proc];

GetObservableType@obs@args :=
   CConversion`ArrayType[CConversion`complexScalarCType, 13];

CalculateObservable[obs@args, structName:_String] := StringJoin[
   structName, ".",
   GetObservableName@obs[in@inN->out@outN, nucl, proc, massless],
   " = ", #, "::calculate_", cxx@in, "_to_", cxx@out, "_for_",
   cxx@proc, cxx@massless, "(", cxx@inN, ", ", cxx@outN, ", ", #, "::Nucleus::",
   cxx@nucl, ", MODEL, qedqcd);"
] & [FlexibleSUSY`FSModelName<>"_l_to_l_conversion"];

];
End[];
