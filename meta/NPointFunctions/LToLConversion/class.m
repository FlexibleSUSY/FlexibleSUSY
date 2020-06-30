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

Begin["FlexibleSUSY`Private`"];

WriteLToLConversionClass[
   extraSLHAOutputBlocks:_List,
   files:{{_?FileExistsQ,_String}..}
] :=
Module[
   {
      observables = DeleteDuplicates@Cases[
         Observables`GetRequestedObservables@extraSLHAOutputBlocks,
         FlexibleSUSYObservable`LToLConversion[pIn_[_]->pOut_[_],__]
      ],
      fields = {}, vertices = {}, additionalVertices = {},
      prototypes = "", npfHeaders = "", definitions = "", npfDefinitions = "",
      masslessNeutralVectorBosons
   },

   If[observables =!= {},
      Print["Creating LToLConversion class ..."];

      fields = DeleteDuplicates[Head/@#&/@observables[[All,1]]/.Rule->List];

      (* additional vertices needed for the calculation *)
      masslessNeutralVectorBosons = Select[
         TreeMasses`GetVectorBosons[],
         And[TreeMasses`IsMassless@#,
            !TreeMasses`IsElectricallyCharged@#,
            !TreeMasses`ColorChargedQ@#
         ]&
      ];

      vertices = Flatten /@ Tuples[
         {
            {CXXDiagrams`LorentzConjugate@#, #} &/@ Flatten@Join[
               TreeMasses`GetSMQuarks[], vertices],
            masslessNeutralVectorBosons
         }
      ];

      {additionalVertices,{npfHeaders,npfDefinitions},{prototypes,definitions}} =
         LToLConversion`create@observables;
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@npf_headers@" -> npfHeaders,
         "@npf_definitions@" -> npfDefinitions,
         "@calculate_prototypes@" -> prototypes,
         "@calculate_definitions@" -> definitions,
         Sequence@@GeneralReplacementRules[]
      }
   ];
   {
      fields,
      DeleteDuplicates@Join[vertices,additionalVertices]
   }
];
WriteLToLConversionClass // Utils`MakeUnknownInputDefinition;

End[];
