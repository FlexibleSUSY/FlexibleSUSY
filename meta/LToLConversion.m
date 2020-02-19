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

getPrototype[inFermion_ -> outFermion_] :=
Module[
   {
      inName = CXXDiagrams`CXXNameOfField@inFermion,
      outName = CXXDiagrams`CXXNameOfField@outFermion
   },
   "Eigen::Array<std::complex<double>,10,1>"<>
   " calculate_"<>inName<>"_to_"<>outName<>"(\n"<>
      "   int generationIndex1,\n"<>
      "   int generationIndex2,\n"<>
      "   const " <> FlexibleSUSY`FSModelName <> "_l_to_l_conversion_in_nucleus_wilson::Nucleus nucleus,\n" <>
      "   const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)"
];

getPrototype[list:{Rule[_,_]..}] :=
   StringJoin@Riffle[ getPrototype /@ list, "\n\n" ];

getPrototype[{}] := "";

getPrototype // Utils`MakeUnknownInputDefinition;

create[inFermion_ -> outFermion_] :=
Module[
   {
      inName = CXXNameOfField@inFermion,
      outName = CXXNameOfField@outFermion,
      npfHeader = "",
      npfDefinition = "",
      calculatePrototype = getPrototype[inFermion -> outFermion],
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

create[list:{Rule[_,_]..}] :=
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
