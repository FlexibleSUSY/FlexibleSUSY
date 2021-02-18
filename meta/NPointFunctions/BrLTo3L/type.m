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

BeginPackage@"BrLTo3L`";Quiet[

BrLTo3L`namespace::usage = "
@brief Returns a namespace for C++ code of a given observable";

BrLTo3L`calculate::usage = "
@brief Creates a C++ calculate call or prototype for the observable.
@param obs The observable.
@param . Defines, whether prototype or call should be generated.
@returns A C++ code for prototype of call of calculate function.";

BrLTo3L`arguments::usage = "
@brief Generates a named pattern sequence (by inserting given symbols in
       correct places), which is set of arguments for observable.
@param lepton A symbol for lepton field name.
@param nI A symbol for number of decaying lepton.
@param nO A symbol for number of outgoing lepton.
@param nA A symbol for number of outgoing lepton-antilepton pair.
@param contribution A symbol for a desired contribution.
@return A sequence with named patterns, which is set of arguments for
        observable.";

];Begin@"`Private`";

namespace[] := "l_to_3l";
namespace // Utils`MakeUnknownInputDefinition;
namespace ~ SetAttributes ~ {Protected, Locked};

Off@RuleDelayed::rhs;
arguments[lepton_Symbol, nI_Symbol -> {nO_Symbol, nA_Symbol},
   contribution_Symbol] :=
Sequence[
   lepton_Symbol?TreeMasses`IsLepton[nI_Integer] ->
      {  lepton_[nO_Integer],
         lepton_[nA_Integer],
         SARAH`bar[lepton_[nA_Integer]]},
   contribution:_Symbol|{__Symbol}];
On@RuleDelayed::rhs;
arguments // Utils`MakeUnknownInputDefinition;
arguments ~ SetAttributes ~ {Protected, Locked};

`type`observable = FlexibleSUSYObservable`BrLTo3L@
   arguments[lep, nI -> {nO, nA}, proc];
`type`observable ~ SetAttributes ~ {Protected, Locked};

With[{i = TextFormatting`IndentText, m = FlexibleSUSY`FSModelName,
      cxx = CConversion`ToValidCSymbolString},
   calculate[obs:`type`observable] := StringJoin[
      m, "_", namespace[], "::calculate_",
      cxx@lep, "_to_", cxx@lep, cxx@lep, cxx@SARAH`bar@lep,
      "_for_", SymbolName@proc, "(", cxx@nI, ", ", cxx@nO, ", ", cxx@nA,
      ", MODEL, qedqcd)"];
   calculate[obs:`type`observable, Head] := StringJoin["calculate_",
      cxx@lep, "_to_", cxx@lep, cxx@lep, cxx@SARAH`bar@lep,
      "_for_", SymbolName@proc, "(\n", i@"int nI, int nO, int nA,\n",
      i@"const ", m, "_mass_eigenstates& model,\n",
      i@"const softsusy::QedQcd& qedqcd)"];];
prototype // Utils`MakeUnknownInputDefinition;
prototype ~ SetAttributes ~ {Protected, Locked};

End[];EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "BrLTo3L`"];
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "BrLTo3L`"];
Protect@$Packages;
