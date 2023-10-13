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

BeginPackage["BrLTo3L`"];
Begin["`Private`"];
`type`observable = FlexibleSUSYObservable`BrLTo3L[lep_@nI_ -> {lep_@nO_, lep_@nA_, SARAH`bar@lep_@nA_}, proc_, loopN_];
End[];
Block[{$ContextPath}, EndPackage[]];

Utils`DynamicInclude@"main.m";

Begin@"FlexibleSUSY`Private`";

WriteClass[obs:FlexibleSUSYObservable`BrLTo3L, blocks_List, files_] :=
Module[{observables, ffvFields = {},
      fermions = {}, ffvV = {}, npfV = {},
      calcProto = "", npfHead = "", calcDef = "", npfDef = ""},
   observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@blocks, _obs];
   If[And[observables =!= {},
         FlexibleSUSY`FSFeynArtsAvailable,
         FlexibleSUSY`FSFormCalcAvailable],
      Print@"\nCreating BrLTo3L class ...";
      fermions = DeleteDuplicates@Cases[observables, {_, f_, bf_} :> {bf, f},
         Infinity] /. f_[_Integer]:>f;
      ffvFields = DeleteDuplicates@Cases[observables,
         Rule[in_, {out, __}] :> {in, out}, Infinity] /. f_[_Integer]:>f;

      ffvV = Flatten/@Tuples@{fermions, {TreeMasses`GetPhoton[]}};

      {npfV, npfHead, npfDef, calcProto, calcDef} = BrLTo3L`create@observables;];
   WriteOut`ReplaceInFiles[files,
      {
         "@npf_headers@" -> npfHead,
         "@npf_definitions@" -> npfDef,
         "@calculate_prototypes@" -> calcProto,
         "@calculate_definitions@" -> calcDef,
         "@include_guard@" -> SymbolName@obs,
         "@namespace@" -> Observables`GetObservableNamespace@obs,
         "@filename@" -> Observables`GetObservableFileName@obs,
         "@get_MSUSY@" -> TextFormatting`IndentText@
            TextFormatting`WrapLines@AMM`AMMGetMSUSY[],
         Sequence@@GeneralReplacementRules[]
      }
   ];
   {ffvFields, Join[ffvV, npfV], {}}
];

End[];
