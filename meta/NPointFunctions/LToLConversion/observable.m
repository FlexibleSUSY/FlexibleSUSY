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

Begin@"Observables`Private`";
With[{args = LToLConversion`arguments[in@inN, out@outN, nucl, proc, loopN],
      obs = FlexibleSUSYObservable`LToLConversion,
      cxx = CConversion`ToValidCSymbolString,
      namespace = LToLConversion`namespace@C},

   GetObservableName@obs@args := StringJoin[
      cxx@in, cxx@inN, "to", cxx@out, cxx@outN, "conversion_in",
      cxx@nucl, "_for", SymbolName@proc, "_", cxx[loopN+0], "loop"];

   GetObservableDescription@obs@args := StringJoin[
      cxx@in, "(", cxx@inN, ") to ", cxx@out, "(", cxx@outN, ") conversion in ",
      cxx@nucl, " for ", SymbolName@proc, " ", cxx[loopN+0], " loop"];

   GetObservableType@obs@args :=
      CConversion`ArrayType[CConversion`complexScalarCType, 13];

   CalculateObservable[obs@args, structName:_String] := StringJoin[
      structName, ".",
      GetObservableName@obs[in@inN -> out@outN, nucl, proc, loopN], " = ",
      namespace, "calculate_", cxx@in, cxx@out, "_for",
      SymbolName@proc, "_", cxx[loopN+0], "loop(", cxx@inN, ", ", cxx@outN, ", ",
      namespace, "Nucleus::", cxx@nucl, ", MODEL, ",
      "ltolconversion_settings, qedqcd);"];];

End[];
