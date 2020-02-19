(* ::Package:: *)

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

BeginPackage@"LToLConversion`";

create::usage = "";

Begin["`internal`"];

getPrototype[lIn_->lOut_,contribution_] :=
   "Eigen::Array<std::complex<double>,10,1>"<>
   " calculate_"<>#@lIn<>"_to_"<>#@lOut<>"_for_"<>#@contribution<>"(\n"<>
   "   int generationIndex1,\n"<>
   "   int generationIndex2,\n"<>
   "   const " <> FlexibleSUSY`FSModelName <>
      "_l_to_l_conversion::Nucleus nucleus,\n" <>
   "   const " <> FlexibleSUSY`FSModelName <>
      "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)" &@
      CConversion`ToValidCSymbolString;

getPrototype // Utils`MakeUnknownInputDefinition;

create[
   FlexibleSUSYObservable`LToLConversion[lIn_[gIn_]->lOut_[gOut_],_,contribution_,_->True]
] :=
Module[
   {
      inName = CConversion`ToValidCSymbolString@lIn,
      outName = CConversion`ToValidCSymbolString@lOut,
      npfHeader = "",
      npfDefinition = "",
      calculatePrototype = getPrototype[lIn->lOut,contribution],
      calculateDefinition, n=2
   },
   calculateDefinition = calculatePrototype <> " {\n" <>
   "   return Eigen::Array<std::complex<double>,10,1>::Zero();\n"<>
   "}\n";

   {
      {},
      {npfHeader,npfDefinition},
      {calculatePrototype <> ";", calculateDefinition}
   }
];

create[list:{__}] :=
   {
      DeleteDuplicates[ Join@@#[[All,1]] ],
      {
         #[[1,2,1]],
         StringJoin@Riffle[#[[All,2,2]], "\n\n"]
      },
      {
         StringJoin@Riffle[#[[All,3,1]], "\n\n"],
         StringJoin@Riffle[#[[All,3,2]], "\n\n"]
      }
   } & [create /@ list];

create[{}] := "";

create // Utils`MakeUnknownInputDefinition;

End[];
EndPackage[];
