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

Utils`DynamicInclude@"NPointFunctions.m";
Begin@"FlexibleSUSY`Private`";

FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`BrLTo3L, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      exportFields = {},
      npfDefinitions = "", obsPrototypes = "", obsDefinitions = "",
      ffvV = {}, npfV = {}, fermions = {}
   },

   If[observables =!= {} && FlexibleSUSY`FSFeynArtsAvailable && FlexibleSUSY`FSFormCalcAvailable,
      Print["Creating ", SymbolName@obs, " class ..."];

      exportFields = Cases[observables, {f_@_, __} :> {f, f}, Infinity];
      fermions =     Cases[observables, {f_@_, __} :> {SARAH`bar@f, f}, Infinity];
      ffvV = Flatten/@Tuples@{fermions, {TreeMasses`GetPhoton[]}};

      {npfV, npfDefinitions, obsPrototypes, obsDefinitions} = BrLTo3L`create@observables;
   ];

   WriteOut`ReplaceInFiles[files,
      {
         "@npointfunctions_headers@" -> NPointFunctions`CreateCXXHeaders[],
         "@npointfunctions_definitions@" -> npfDefinitions,
         "@calculate_prototypes@" -> obsPrototypes,
         "@calculate_definitions@" -> obsDefinitions,
         "@include_guard@" -> SymbolName@obs,
         "@namespace@" -> Observables`GetObservableNamespace@obs,
         "@filename@" -> Observables`GetObservableFileName@obs,
         "@get_MSUSY@" -> TextFormatting`IndentText@TextFormatting`WrapLines@AMM`AMMGetMSUSY[],
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];
   {
      "FFV fields" -> DeleteDuplicates@exportFields,
      "C++ vertices" -> DeleteDuplicates@Join[ffvV, npfV]
   }
];

End[];
