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

Begin@"FlexibleSUSY`Private`";

With[{main = FileNameJoin@{DirectoryName@$Input, "main.m"}},
WriteLToLConversionClass[
   extraSLHAOutputBlocks:_List,
   files:{{_?FileExistsQ,_String}..}] :=
Module[
   {
      observables = DeleteDuplicates@Cases[
         Observables`GetRequestedObservables@extraSLHAOutputBlocks,
         FlexibleSUSYObservable`LToLConversion[pIn_[_]->pOut_[_],__]
      ],
      fields = {}, vertices = {}, additionalVertices = {},
      prototypes = "", npfHeaders = "", definitions = "", npfDefinitions = "",
      masslessNeutralVectorBosons, newRules
   },

   newRules = {
      "@LToLConversion_fill@" -> "slha_io.fill(ltolconversion_settings);",
      "@LToLConversion_init@" -> "LToLConversion_settings ltolconversion_settings;",
      "@LToLConversion_class_name@" -> "ltolconversion_settings,",
      "@LToLConversion_class_declaration@" -> "class LToLConversion_settings;",
      "@LToLConversion_named_argument@" -> "const LToLConversion_settings& ltolconversion_settings,",
      (* Live inside librarylink. *)
      "@LToLConversion_private@" -> "LToLConversion_settings ltolconversion_settings{}; ///< LToLConversion settings",
      "@LToLConversion_setter@" -> "void set_ltolconversion_settings(const LToLConversion_settings& s) { ltolconversion_settings = s; }",
      "@LToLConversion_set_data@" -> "data.set_ltolconversion_settings(ltolconversion_settings);",
      "@LToLConversion_set_slha@" -> "slha_io.set_LToLConversion_settings(ltolconversion_settings);",
      "@LToLConversion_reset@" -> "ltolconversion_settings.reset();"};

   If[Or[observables === {},
         Not@FlexibleSUSY`FSFeynArtsAvailable,
         Not@FlexibleSUSY`FSFormCalcAvailable],
      newRules = newRules /. Rule[x_, _]:> Rule[x, ""];];
   If[And[observables =!= {},
         FlexibleSUSY`FSFeynArtsAvailable,
         FlexibleSUSY`FSFormCalcAvailable],
      Print["Creating LToLConversion class ..."];
      Get@main;
      fields = DeleteDuplicates[Head/@#&/@observables[[All,1]]/.Rule->List];

      (* additional vertices needed for the 1 loop calculation *)
      masslessNeutralVectorBosons =
         Select[TreeMasses`GetVectorBosons[],
            And[TreeMasses`IsMassless@#,
               !TreeMasses`IsElectricallyCharged@#,
               !TreeMasses`ColorChargedQ@#]&];

      vertices = Flatten/@ Tuples[
         {  {CXXDiagrams`LorentzConjugate@#, #}&/@
               Flatten@Join[TreeMasses`GetSMQuarks[], vertices],
            masslessNeutralVectorBosons}];

      {additionalVertices,
         {npfHeaders, npfDefinitions},
         {prototypes, definitions}} = LToLConversion`create@observables;];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@npf_headers@" -> npfHeaders,
         "@npf_definitions@" -> npfDefinitions,
         "@calculate_prototypes@" -> prototypes,
         "@calculate_definitions@" -> definitions,
         "@get_MSUSY@" -> TextFormatting`IndentText@
            TextFormatting`WrapLines@AMuon`AMuonGetMSUSY[],
         Sequence@@GeneralReplacementRules[]
      }
   ];
   {
      fields,
      DeleteDuplicates@Join[vertices, additionalVertices],
      newRules
   }
];
];
WriteLToLConversionClass // Utils`MakeUnknownInputDefinition;
WriteLToLConversionClass // Protect;

End[];
