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
Utils`DynamicInclude@"main.m";

Begin@"FlexibleSUSY`Private`";

WriteBrLTo3LClass::usage = "
@brief Takes corresponding to the observable .hpp.in and .cpp.in files from
       template directory. Then inserts parts of code, generated with help of
       main.m and puts resulting files inside directory of configured model.
@param blocks A list of extra slha output blocks. It should contain calls for
       the observable in order to calculate it.
@param files A list of input-output file names to work with.
@returns A list with two entries: 1) a list of all external fields of the
         observable, 2) a list of all vertices, required by the observable.";
WriteBrLTo3LClass[blocks:_List, files:{{_?FileExistsQ, _String}..}] :=
Module[{obs, ffvFields = {},
      fermions = {}, ffvV = {}, npfV = {},
      calcProto = "", npfHead = "", calcDef = "", npfDef = ""},
   obs = DeleteDuplicates@Cases[
      Observables`GetRequestedObservables@blocks,
      _FlexibleSUSYObservable`BrLTo3L];
   If[And[obs =!= {},
         FlexibleSUSY`FSFeynArtsAvailable,
         FlexibleSUSY`FSFormCalcAvailable],
      Print@"\nCreating BrLTo3L class ...";
      fermions = DeleteDuplicates@Cases[obs, {_, f_, bf_} :> {bf, f},
         Infinity] /. f_[_Integer]:>f;
      ffvFields = DeleteDuplicates@Cases[obs,
         Rule[in_, {out, __}] :> {in, out}, Infinity] /. f_[_Integer]:>f;

      ffvV = Flatten/@Tuples@{fermions, {TreeMasses`GetPhoton[]}};

      {npfV, npfHead, npfDef, calcProto, calcDef} = BrLTo3L`create@obs;];
   WriteOut`ReplaceInFiles[files,
      {
         "@npf_headers@" -> npfHead,
         "@npf_definitions@" -> npfDef,
         "@calc_prototypes@" -> calcProto,
         "@calc_definitions@" -> calcDef,
         "@get_MSUSY@" -> TextFormatting`IndentText@
            TextFormatting`WrapLines@AMM`AMMGetMSUSY[],
         Sequence@@GeneralReplacementRules[]
      }
   ];
   {ffvFields, Join[ffvV, npfV], {}}
];
WriteBrLTo3LClass // Utils`MakeUnknownInputDefinition;
WriteBrLTo3LClass // Protect;

End[];
