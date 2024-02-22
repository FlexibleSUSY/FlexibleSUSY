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

BeginPackage@"NPointFunctions`";

VerticesForNPointFunction;
CreateCXXHeaders;
CreateCXXFunctions;
NPFDefinitions;
GetProcess;
GetSubexpressions;
GetGenericSums;
GetParticleName;
NPointFunction;
IsNPointFunction;
ApplySubexpressions;

SetAttributes[
   {
      LoopLevel, Regularize, UseCache, ZeroExternalMomenta, OnShellFlag,
      OperatorsOnly, ExceptLoops, KeepProcesses, Irreducible, Triangles,
      GenericS, GenericF, GenericV, GenericU, GenericSum, GenericIndex,
      LorentzIndex, Mat, DiracChain, Observable
   },
   {Protected}
];

Begin@"`Private`";

IsWilsonBasis[list_] := MatchQ[list, {Rule[_String,_]..}];

IsParticle[p_] := Or[TreeMasses`IsScalar@p, TreeMasses`IsFermion@p, TreeMasses`IsVector@p, TreeMasses`IsGhost@p];

IsGenericParticle[field_] := MatchQ[field,
   (GenericS | GenericF | GenericV | GenericU)[GenericIndex[_Integer]] |
   SARAH`bar[(GenericF | GenericU)[GenericIndex[_Integer]]] |
   Susyno`LieGroups`conj[(GenericS | GenericV)[GenericIndex[_Integer]]]
];

IsSubexpressions[list_] := MatchQ[list, {Rule[_Symbol,_]...}];

IsSummation[list_] := MatchQ[list, {{_?IsGenericParticle, _}..}];

IsGenericSum[sum_] := MatchQ[sum, GenericSum[{__}, _?IsSummation]];

IsClassFields[list_] := MatchQ[list, {{__}..}];

IsCombinatoricalFactors[list_] := MatchQ[list, {__Integer}];

IsColorProjector[function_] := MatchQ[function, Identity|SARAH`Delta];

IsColorFactors[list_] := MatchQ[list, {__}];

IsNPointFunction[e_] := MatchQ[e, {{{__}, {__}}, {{{__?IsGenericSum}, {__?IsClassFields}, {__?IsCombinatoricalFactors}, {__?IsColorFactors}}, _?IsSubexpressions}}];

secure::usage = "Make sure that imput of important functions is secured.";
secure[sym_Symbol] := Protect@Evaluate@Utils`MakeUnknownInputDefinition@sym;
secure // secure;

GetIndex[obj_?IsGenericParticle] := First@Level[obj, {-1}];
GetIndex // secure;

GetLength[obj_?IsWilsonBasis] := ToString@Length@obj;
GetLength // secure;

GetProcess[obj_?IsNPointFunction] := obj[[1]];
GetProcess // secure;

GetExternalMomenta[obj_?IsNPointFunction] :=
DeleteDuplicates@Cases[
   {GetGenericSums@obj, GetSubexpressions@obj},
   HoldPattern@SARAH`Mom[_Integer, ___],
   Infinity
];
GetExternalMomenta // secure;

GetExternalIndices[obj_?IsNPointFunction] := DeleteDuplicates@Flatten@Level[GetProcess@obj, {4, 5}];
GetExternalIndices // secure;

GetIndexRange[obj_List] :=  {1, Length@obj};
GetIndexRange // secure;

GetClassFields[obj_?IsNPointFunction] := obj[[2, 1, 2]];
GetClassFields // secure;

GetCombinatoricalFactors[obj_?IsNPointFunction] := obj[[2, 1, 3]];
GetCombinatoricalFactors // secure;

GetColorFactors[obj_?IsNPointFunction] := obj[[2, 1, 4]];
GetColorFactors // secure;

GetSubexpressions[obj_?IsNPointFunction] := obj[[2, 2]];
GetSubexpressions // secure;

GetGenericSums[obj_?IsNPointFunction] := obj[[2, 1, 1]];
GetGenericSums // secure;

GetParticleName[particle_?IsParticle] :=
   CXXDiagrams`RemoveLorentzConjugation[particle] /. {_@n_@_List :> n, n_@_List :> n, _@n_ :> n};
GetParticleName[particle_?IsGenericParticle] :=
   CXXDiagrams`RemoveLorentzConjugation[particle] /. {n_@_ :> n};
GetParticleName // secure;

GetGenericFields[obj_?IsSummation] := First/@obj;
GetGenericFields[obj_?IsGenericSum] := First/@Last[obj];
GetGenericFields[list:{__?IsGenericSum}] := GetGenericFields/@list;
GetGenericFields // secure;

GetSumExpression[obj_?IsGenericSum] := First@obj;
GetSumExpression // secure;

GetSumParticles[obj_?IsGenericSum] := Last@obj;
GetSumParticles // secure;

GetClassFieldRules[obj_?IsNPointFunction] :=
MapThread[
   Function[fields,MapThread[Rule,{#1,fields}]]/@#2&,
   {GetGenericFields@GetGenericSums@obj, GetClassFields@obj}
];
GetClassFieldRules // secure;

SetSubexpressions[old_?IsNPointFunction, new_?IsSubexpressions] := ReplacePart[old, {2, 2} -> new];
SetSubexpressions // secure;

ApplySubexpressions[obj_?IsNPointFunction] :=
Module[{result},
   If[{} === GetSubexpressions@obj, Return@obj];
   WriteString[$Output, "Applying subexpressions ... "];
   result = SetSubexpressions[
      ReplacePart[
         obj,
         {2, 1, 1} -> ReplaceRepeated[GetGenericSums@obj, GetSubexpressions@obj]
      ],
      {}
   ];
   WriteString[$Output, "done\n"];
   result
];
ApplySubexpressions // secure;

subWrite[str_String] := WriteString[$Output, str];
subWrite // secure;

OutputPaths[] := Module[{outputPath, eigenstates = ToString@FlexibleSUSY`FSEigenstates},
   outputPath = {SARAH`$sarahCurrentOutputMainDir, eigenstates};
   FileNameJoin@*Flatten /@
   {
      {outputPath, "NPointFunctions"},
      {outputPath, "FormCalc"},
      {outputPath, "FeynArts", SARAH`ModelName<>eigenstates},
      {outputPath, "FeynArts", "ParticleNamesFeynArts.dat"},
      {outputPath, "FeynArts", "ParticleNamespaces.m"}
   }
];

CheckOptionValues[opts___] := (
   Cases[{opts},
      Rule[First@#, v_] :>
         Utils`AssertOrQuit[MatchQ[v, Alternatives@@Last@#], NPointFunction::errOption, v, First@#, Last@#]
   ] &/@
   {
      OnShellFlag         -> {True, False},
      UseCache            -> {True, False},
      ZeroExternalMomenta -> {True, False, OperatorsOnly, ExceptLoops},
      KeepProcesses       -> {{__Symbol}},
      LoopLevel           -> {0, 1},
      Observable          -> {None, _[___]},
      Regularize          -> {FlexibleSUSY`MSbar, FlexibleSUSY`DRbar}
   };
   True
);

NPointFunction::errOption = "Value `1` for `2` option must match any from `3`.";
NPointFunction[
   inFields: {__?(TreeMasses`IsScalar@# || TreeMasses`IsFermion@# || TreeMasses`IsVector@# &)},
   outFields:{__?(TreeMasses`IsScalar@# || TreeMasses`IsFermion@# || TreeMasses`IsVector@# &)},
   opts:___ /; CheckOptionValues@opts
] :=
Module[{
      npfDir, fcDir, faModel, particleNamesFile,
      particleNamespaceFile, kernel, currentDirectory,
      npf, metaInfo = {{inFields, outFields}, {opts}}
   },
   {npfDir, fcDir, faModel, particleNamesFile, particleNamespaceFile} = OutputPaths[];
   If[DirectoryQ@npfDir === False, CreateDirectory@npfDir];
   If[Utils`FSGetOption[{opts}, UseCache],
      npf = GetCache[inFields, outFields, npfDir, metaInfo];
      If[npf =!= Null, Return@npf];
   ];
   If[FileExistsQ[faModel <> ".mod"] === False,
      {kernel} = LaunchKernels@1;
      WriteModelFile@kernel;
      WriteNamespaceFile@particleNamespaceFile;
      CloseKernels@kernel;
   ];

   {kernel} = LaunchKernels@1;
   SetSharedFunction[subWrite, Print];
   With[{path = $Path,
         data = {fcDir, faModel, particleNamesFile, particleNamespaceFile,
            ToFeynArtsFields[inFields, particleNamesFile],
            ToFeynArtsFields[outFields, particleNamesFile]
         },
         options = {Utils`FSGetOption[{opts}, Observable],
            Utils`FSGetOption[{opts}, LoopLevel],
            SymbolName/@If[List=!=Head@#, {#}, #]&@Utils`FSGetOption[{opts}, KeepProcesses],
            Utils`FSGetOption[{opts}, ZeroExternalMomenta],
            Utils`FSGetOption[{opts}, OnShellFlag],
            Utils`FSGetOption[{opts}, Regularize]
         },
         meta = FlexibleSUSY`$flexiblesusyMetaDir
      },
      npf = DeleteZeroGenericSum@ParallelEvaluate[
         $Path = path;
         Get@FileNameJoin@{meta, "NPointFunctions", "NPointFunction.m"};
         NPointFunction[data, options],
         kernel,
         DistributedContexts -> None
      ];
   ];
   CloseKernels@kernel;
   UnsetShared[subWrite, Print];

   WriteCache[npf, npfDir, metaInfo];
   npf
];
NPointFunction // secure;

VerticesForNPointFunction[obj_?IsNPointFunction] :=
Module[{v, getVertex},
   getVertex[vertGen_, rules_] := vertGen/.#&/@rules;
   v = DeleteDuplicates[
      Cases[#, SARAH`Cp[f__] :> {f}, Infinity, Heads->True] &/@ GetGenericSums@obj
   ];
   DeleteDuplicates[
      Vertices`StripFieldIndices/@ #&/@ Flatten[MapThread[getVertex, {v, GetClassFieldRules@obj}], 2]
   ]
];
VerticesForNPointFunction // secure;

With[{dir = DirectoryName@$InputFileName},
CacheNameForMeta[meta:{__}] :=
Module[{inFields, outFields, dirMod = "", obsDir},
   inFields  = meta[[1, 1]];
   outFields = meta[[1, 2]];
   opts = SortBy[Cases[meta[[2]], _Rule], First];
   obsDir = FileNameJoin@{dir, "Observables", SymbolName@Head@Utils`FSGetOption[opts, Observable]};
   If[DirectoryQ@obsDir,
      dirMod = "modtime" <> DateString@FileDate[obsDir, "Modification"];
   ];
   StringReplace[
      Utils`StringJoinWithSeparator[
         Flatten@{inFields, "to", outFields, (#/.x_Symbol :> SymbolName@x)&/@Last/@opts , dirMod},
         "_",
         ToString
      ] <> ".m",
      {" " -> "", ":" -> "-", "Blank[]" -> "_"}
   ]
];
CacheNameForMeta // secure;
];

WriteCache[new_, dir_, meta_] :=
Module[{path, file, npfs, iosition},
   path = FileNameJoin@{dir, CacheNameForMeta@meta};
   If[FileExistsQ@path,
      npfs = Get@path;,
      npfs = {};
   ];
   position = Position[npfs[[All,1]], new[[1]]];
   If[Length@position === 1,
      npfs[[position[[1]]]] = new;,
      AppendTo[npfs, new];
   ];
   file = OpenWrite@path;
   Write[file, npfs];
   Close@file;
];
WriteCache // secure;

GetCache[inFields_, outFields_, cacheDir_, meta_] :=
Module[{nPointFunctionsFile, nPointFunctions, position},
   nPointFunctionsFile = FileNameJoin@{cacheDir, CacheNameForMeta@meta};
   If[FileExistsQ@nPointFunctionsFile === False,
      Return@Null
   ];
   nPointFunctions = Get@nPointFunctionsFile;
   position = Position[
      Vertices`StripFieldIndices@nPointFunctions[[All,1]],
      {inFields, outFields}
   ];
   If[Length@position == 1, nPointFunctions[[position[[1,1]]]], Null]
];
GetCache // secure;

WriteModelFile[kernel_] :=
With[{currentPath = $Path,
      currentDir = Directory[],
      metaDir = FlexibleSUSY`$flexiblesusyMetaDir,
      sarahInputDirs = SARAH`SARAH@SARAH`InputDirectories,
      sarahOutputDir = SARAH`SARAH@SARAH`OutputDirectory,
      SARAHModelName = If[SARAH`submodeldir =!= False,
         SARAH`modelDir <> "-" <> SARAH`submodeldir,
         SARAH`modelDir
      ],
      eigenstates = FlexibleSUSY`FSEigenstates
   },
   Print["Generating FeynArts model file ..."];
   SetSharedFunction[Print];
   ParallelEvaluate[
      $Path = currentPath;
      SetDirectory@currentDir;
      Get@FileNameJoin@{metaDir, "NPointFunctions", "MakeModelFile.m"};
      NPointFunctions`MakeModelFile[sarahInputDirs,sarahOutputDir, SARAHModelName, eigenstates];,
      kernel,
      DistributedContexts -> None
   ];
   UnsetShared[Print];
   Print["Generating FeynArts model file ... done"];
];

WriteNamespaceFile[fileName_String] :=
Module[{fileHandle = OpenWrite@fileName},
   Write[fileHandle, {ToString@#, Context@#}&/@ TreeMasses`GetParticles[]];
   Close@fileHandle;
];

ToFeynArtsFields[fields_, namesFile_String] :=
Module[{unique, faNames},
   unique = DeleteDuplicates[CXXDiagrams`RemoveLorentzConjugation/@fields];
   faNames = Flatten[
      StringCases[Utils`ReadLinesInFile@namesFile,
         SymbolName@# ~~ ": " ~~ x__ ~~ "]" ~~ ___ :> "FeynArts`" <> x <> "]"
      ]&/@ unique
   ];
   fields /. MapThread[Rule, {unique, faNames}] /.
      {
         SARAH`bar@field_String :> "-" <> field,
         Susyno`LieGroups`conj@field_String :> "-" <> field
      }
];

DeleteZeroGenericSum[npf_?IsNPointFunction] := npf;
DeleteZeroGenericSum[{fields_, {{sums_, rules_, comb_, col_}, subs_}}]:=
Module[{pos = Position[sums, GenericSum[{0}, {}]]},
   Print["Removing zero GenericSum"];
   {fields, {Delete[#,pos]&/@{sums,rules,comb,col}, subs}}
];
DeleteZeroGenericSum // secure;

LorentzConjugation[obj_?IsGenericParticle] :=
Switch[Head@obj,
   SARAH`bar | Susyno`LieGroups`conj, CXXDiagrams`RemoveLorentzConjugation@obj,
   GenericS | GenericV,               Susyno`LieGroups`conj@obj,
   GenericF | GenericU,               SARAH`bar@obj
];
LorentzConjugation[obj_?TreeMasses`IsScalar|_?TreeMasses`IsVector] := Susyno`LieGroups`conj@obj;
LorentzConjugation[obj_?TreeMasses`IsFermion] := SARAH`bar@obj;
LorentzConjugation // secure;

CXXIndex[obj_?IsGenericParticle] :=
   "i"<>StringTake[SymbolName[obj[[0]]], -1]<>ToString[obj[[1, 1]]] &@ LorentzConjugation@obj;
CXXIndex[obj:_?IsParticle] :=
Module[{nakedField},
   nakedField = CXXDiagrams`RemoveLorentzConjugation@obj;
   Switch[nakedField,
   _Symbol,              "",
   (_Symbol)[{_Symbol}], StringDrop[ToString[nakedField[[1, 1]]], 2]]
];
CXXIndex // secure;

CreateCXXHeaders[] :=
TextFormatting`ReplaceCXXTokens["
   #include \"loop_libraries/loop_library.hpp\"
   #include \"cxx_qft/@ModelName@_npointfunctions_wilsoncoeffs.hpp\"
   #include \"concatenate.hpp\"
   #include <limits>
   #include <type_traits>
   #include <boost/fusion/include/at_key.hpp>
   #include \"wrappers.hpp\"
   ",
   {"@ModelName@" -> FlexibleSUSY`FSModelName}
];
CreateCXXHeaders // secure;

CreateCXXFunctions[
   npf_?IsNPointFunction,
   functionName_String,
   colorProjector_?IsColorProjector,
   wilsonBasis:_?IsWilsonBasis:{"value"->"dummy string"}] :=
Module[{GetFunctionPrototype, prototype, definition},
   GetFunctionPrototype = StringReplace[
      "std::array<std::complex<double>, L> name(arg)",
      {"L" -> GetLength@wilsonBasis, "name" -> functionName, "arg" -> #}
   ]&;
   SetHelperClassName@npf;
   SetBasis@wilsonBasis;
   SetColorProjector@colorProjector;
   prototype = GetFunctionPrototype@GetCXXArguments[npf, Default]<>";";
   definition = TextFormatting`ReplaceCXXTokens["
      @class_body@

      @npf_function@ {
         @class_name@ helper{model, indices, momenta};
         return helper.calculate();
      }",
      {
         "@class_name@" -> GetHelperClassName[],
         "@class_body@" -> CreateClassBody@ApplySubexpressions@npf,
         "@npf_function@" -> GetFunctionPrototype@GetCXXArguments@npf
      }
   ];
   {prototype, definition}
];

CreateCXXFunctions // secure;

NPFDefinitions[first__, basis:{__String}] :=
   CreateCXXFunctions[first, Rule[#, "dummy"]&/@basis][[2]];
NPFDefinitions // secure;

GetCXXArguments[npf_?IsNPointFunction, control_:Null] :=
TextFormatting`ReplaceCXXTokens["
   const @eigenstates@ &model,
   const std::array<int, @indices_length@> &indices,
   const std::array<Eigen::Vector4d, @momenta_length@> &momenta
   @default@",
   {
      "@eigenstates@" -> FlexibleSUSY`FSModelName<>"_mass_eigenstates",
      "@indices_length@" -> ToString@Length@GetExternalIndices@npf,
      "@momenta_length@" -> ToString@Length@GetExternalMomenta@npf,
      "@default@" -> If[control === Default,
            "= { "<> Utils`StringJoinWithReplacement@Array["Eigen::Vector4d::Zero()", Length@GetExternalMomenta@npf]<>" }", ""
         ]
   }
];
GetCXXArguments // secure;

Module[{helperClassName, basis, projector},
SetHelperClassName[obj_?IsNPointFunction] :=
Module[{fieldNames = Vertices`StripFieldIndices/@Join@@GetProcess[obj]},
   helperClassName = "nPoint" <>
      StringJoin@Map[ToString, fieldNames/.a_[b_] :> Sequence@@{a,b}] <> "_" <>
      ToString@Ceiling[10^6*AbsoluteTime[]];
];
SetHelperClassName // secure;
GetHelperClassName[] := helperClassName;

SetBasis[obj_?IsWilsonBasis] := (basis = obj);
SetBasis // secure;
GetBasis[] := basis;

SetColorProjector[colorProjector_?IsColorProjector] := (projector = colorProjector);
SetColorProjector // secure;
GetColorProjector[] := projector;
];

CreateClassBody[npf_?IsNPointFunction] :=
Module[{genSums, extIndices, numberOfMomenta, genFields, genSumNames},
   genSums = GetGenericSums@npf;
   extIndices = GetExternalIndices@npf;
   numberOfMomenta = Length[GetExternalMomenta@npf];
   genFields = DeleteDuplicates[Flatten@GetClassFieldRules@npf /. Rule[x_,_] :> x];
   genSumNames = Array["genericSum"<>ToString@#&, Length@genSums];
   SetCXXRules[extIndices, genFields];
   TextFormatting`ReplaceCXXTokens["
      class @ClassName@ : public @Context@
      {
         using generic_sum_base = @Context@;

         template<class GenericFieldMap>
         struct subexpression_base :
         generic_sum_base, index_map_interface<GenericFieldMap>
         {
            subexpression_base( const subexpression_base & ) = default;

            subexpression_base( const generic_sum_base &gsb,
               const typename field_index_map<GenericFieldMap>::type &fim ) :
            generic_sum_base( gsb ), index_map_interface<GenericFieldMap>( fim )
            {}
         }; // End of subexpression_base<GenericFieldMap>

         @KeyStructsInitialization@

         @GenericSums@

         public:
         @ClassName@( @Arguments@ ) :
         @Context@ { model, indices, momenta }
         {}

         @CalculateFunction@
      }; // End of @ClassName@"
      ,
      {
         "@ClassName@" -> GetHelperClassName[],
         "@Context@" -> "correlation_function_context<"<>ToString@Length@extIndices<>","<>ToString@numberOfMomenta<>">",
         "@KeyStructsInitialization@" -> InitializeKeyStructs@genFields,
         "@GenericSums@" -> CXXGenericSum[npf, genSumNames],
         "@Arguments@" -> GetCXXArguments@npf,
         "@CalculateFunction@" -> CXXCalculateFunction@genSumNames
      }
   ]
];
CreateClassBody // secure;

InitializeKeyStructs[fields:{__?IsGenericParticle}]:=
Utils`StringJoinWithSeparator[
   "struct "<>#<>" {};" &/@ CXXGenericKey/@fields,
   "\n"
];
InitializeKeyStructs // secure;

Module[{rules = {{}}},
SetCXXRules[extIdx:{___Symbol}, genericFields:{__?IsGenericParticle}] :=
Module[{externalIndexRules, short, index, genericRules, massRules, couplingRules},
   externalIndexRules = MapIndexed[#1 -> "i"<>ToString@@#2&, extIdx];
   genericRules = Flatten[
      {
         LorentzConjugation@# -> CXXFullFieldName[LorentzConjugation@#][CXXIndex@#],
         # -> CXXFullFieldName[#][CXXIndex@#]
      } &/@ genericFields
   ];

   short[fields__] := Utils`StringJoinWithSeparator[CXXShortFieldName/@{fields}, ", "];
   index[fields__] := Utils`StringJoinWithSeparator[CXXFieldIndices/@{fields}, ", "];
   couplingRules =
   {
      SARAH`Cp[fields__][1] :>                                 StringJoin["NPF_S(", short@fields, ") NPF_I(", index@fields, ")"],
      SARAH`Cp[fields__][SARAH`PL] :>                          StringJoin["NPF_L(", short@fields, ") NPF_I(", index@fields, ")"],
      SARAH`Cp[fields__][SARAH`PR] :>                          StringJoin["NPF_R(", short@fields, ") NPF_I(", index@fields, ")"],
      SARAH`Cp[s1_, s2_, v_][SARAH`Mom@s1_-SARAH`Mom@s2_] :>   StringJoin["NPF_SSV(", short[s1, s2, v], ") NPF_D(0, 1) NPF_I(", index[s1, s2, v], ")"],
      SARAH`Cp[s1_, s2_, v_][SARAH`Mom@s2_-SARAH`Mom@s1_] :>   StringJoin["NPF_SSV(", short[s1, s2, v], ") NPF_D(1, 0) NPF_I(", index[s1, s2, v], ")"],
      SARAH`Cp[fields__][SARAH`g[_,_]] :>                      StringJoin["NPF_G(", short@fields, ") NPF_I(", index@fields, ")"],
      SARAH`Cp[fields__][
         (SARAH`Mom[f2_,_]-SARAH`Mom[f1_,_])*SARAH`g[_,_],
         (SARAH`Mom[f1_,_]-SARAH`Mom[f3_,_])*SARAH`g[_,_],
         (SARAH`Mom[f3_,_]-SARAH`Mom[f2_,_])*SARAH`g[_,_]] :>  StringJoin["NPF_T(", short@fields, ") NPF_I(", index@fields, ")"]
   };

   massRules =
   {
      SARAH`Mass[f_] :> StringJoin["context.mass<", CXXShortFieldName@f, ">(", CXXFieldIndices@f, ")"]
   };

   rules = {externalIndexRules, massRules, couplingRules, genericRules};
];
SetCXXRules // secure;

ApplyCXXRules[obj_] :=
StringReplace[
   Parameters`ExpressionToString[Fold[ReplaceAll, obj, rules]],
   "\"" -> ""
   ];
];

Module[{sandwich, n},

sandwich[f_] = Switch[Head@f,
   SARAH`bar, "typename cxx_diagrams::bar<" <> # <> ">::type",
   Susyno`LieGroups`conj, "typename cxx_diagrams::conj<" <> # <> ">::type",
   _, #]&;

CXXFullFieldName[f_?IsParticle] := (
   n = CXXDiagrams`RemoveLorentzConjugation@f;
   sandwich[f]["fields::" <> ToString@Switch[n, _Symbol, n, _, Head@n]]
);

CXXFullFieldName[f_?IsGenericParticle] := (
   n = CXXDiagrams`RemoveLorentzConjugation@f;
   sandwich[f]["g"<>StringTake[ToString@Head@n, -1] <> ToString@GetIndex@n]
);

CXXFullFieldName // secure;
];

CXXFieldIndices[SARAH`bar[field_]] :=
   CXXFieldIndices[SARAH`bar[field]] =
      CXXFieldIndices@field;

CXXFieldIndices[Susyno`LieGroups`conj[field_]] :=
   CXXFieldIndices[Susyno`LieGroups`conj[field]] =
      CXXFieldIndices@field;

CXXFieldIndices[head_[GenericIndex[index_Integer]]] :=
   CXXFieldIndices[head[GenericIndex[index]]] =
      "i"<>StringTake[SymbolName@head, -1]<>ToString@index;

CXXFieldIndices[field_] := CXXFieldIndices[field] =
   If[Length@field === 0, "i0", field[[1, 1]]];

CXXFieldIndices // Utils`MakeUnknownInputDefinition;

CXXGenericKey[fields:{__}] := Utils`StringJoinWithSeparator[fields, ", ", CXXGenericKey];

CXXGenericKey[head_[GenericIndex[index_Integer]]] := SymbolName@head<>ToString@index<>"Key";

CXXGenericKey // secure;

RemoveNumbers[expr_] := expr /. {_Integer-> 1, _Rational-> 1, _Complex-> 1, Pi-> 1};

ExtractColorFactor[colorfactors:{__?IsColorFactors}, projection_?IsColorProjector] :=
Module[{projectedFactors},
   projectedFactors = Switch[projection,
      Identity,
         colorfactors,
      _,
         Utils`AssertOrQuit[
            And@@Flatten[MatchQ[#, _projection]&/@#&/@ RemoveNumbers@colorfactors],
            ExtractColorFactor::errMultipleStructures,
            colorfactors,
            projection
         ];
         colorfactors /. _projection -> 1
   ];
   Utils`AssertOrQuit[
      And@@Flatten[NumericQ/@#&/@ projectedFactors],
      ExtractColorFactor::errNotNumber,
      projectedFactors
   ];
   projectedFactors
];
ExtractColorFactor::errMultipleStructures = "Some color factors in `1` have no head `2`.";
ExtractColorFactor::errNotNumber = "Some color factors in `1` are not numbers.";
ExtractColorFactor // secure;

CXXGenericSum[npf_?IsNPointFunction, genSumNames:{__String}] :=
Utils`StringJoinWithSeparator[
   MapThread[
      CXXGenericSum[##, GetSubexpressions@npf, npf]&,
      {
         GetGenericSums@npf,
         GetClassFields@npf,
         GetCombinatoricalFactors@npf,
         ExtractColorFactor[GetColorFactors@npf, GetColorProjector[]],
         genSumNames
      }
   ],
   "\n\n"
];

CXXGenericSum[
   sum_?IsGenericSum,
   genericInsertions_?IsClassFields,
   combinatorialFactors_?IsCombinatoricalFactors,
   colorFactors:{__?NumericQ},
   genSumName_String,
   subexpressions_?IsSubexpressions,
   npf_?IsNPointFunction
] :=
TextFormatting`ReplaceCXXTokens["
   template<class GenericFieldMap>
   struct @GenericSum_NAME@_impl : generic_sum_base {
      @GenericSum_NAME@_impl( const generic_sum_base &base ) :
      generic_sum_base( base ) {
      } // End of constructor @GenericSum_NAME@_impl

      std::array<std::complex<double>,@WilsonBasisLength@> operator()( void ) {
         using boost::mpl::at;
         using boost::fusion::at_key;
         @GenericFieldShortNames@

         @ExternalIndices@
         typename field_index_map<GenericFieldMap>::type index_map;
         const context_with_vertices &context = *this;
         @InitializeOutputVars@

         @SummationOverGenericFields@

         return @ReturnOutputVars@;
      } // End of operator()( void )
   }; // End of struct @GenericSum_NAME@_impl<GenericFieldMap>

   std::array<std::complex<double>,@WilsonBasisLength@> @GenericSum_NAME@( void ) {
      using GenericKeys = boost::mpl::vector< @GenericKeys@ >;
      using GenericInsertions = boost::mpl::vector<
         @ClassInsertions@
         >;
      using combinatorial_factors = boost::mpl::vector<
         @CombinatoricalFactors@
         >;
      using color_factors = boost::mpl::vector<
         @ColorFactors@
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         color_factors,
         boost::mpl::int_<@WilsonBasisLength@>,
         @GenericSum_NAME@_impl
         >( *this );
   } // End of function @GenericSum_NAME@()",
   {
      "@GenericSum_NAME@" -> genSumName,
      "@GenericFieldShortNames@" -> CXXGenericAbbreviations@GetGenericFields@sum,
      "@ExternalIndices@" -> InitializeExternalIndices@npf,
      "@InitializeOutputVars@" -> Utils`StringJoinWithSeparator["std::complex<double> "<>#<>" = 0.0;"&/@First/@GetBasis[],"\n"],
      "@SummationOverGenericFields@" -> CXXSumOverGeneric[GetSumParticles@sum, GetSumExpression@sum],
      "@ReturnOutputVars@" -> ToString[First/@GetBasis[]],
      "@GenericKeys@" -> CXXGenericKey@GetGenericFields@sum,
      "@ClassInsertions@" -> CXXInsertFields@genericInsertions,
      "@CombinatoricalFactors@" -> CXXInsertFactors@combinatorialFactors,
      "@ColorFactors@" -> CXXInsertColors@colorFactors,
      "@WilsonBasisLength@" -> GetLength@GetBasis[]
   }
];
CXXGenericSum // secure;
CXXGenericSum::errColors = "Color factor is not a number after projection: `1`";

InitializeExternalIndices[npf_?IsNPointFunction] :=
Module[{extIndices = GetExternalIndices@npf},
   indices = Array[
      "std::array<int, 1> i" <> ToString@# <>
         " {this->external_indices("<>ToString[#-1]<>")};\n"&,
      Length@extIndices];
   If[Length@extIndices < Length[Flatten@GetProcess@npf],
      PrependTo[indices, "std::array<int, 0> i0 {};\n"];];
   StringJoin@indices];
InitializeExternalIndices // secure;

CXXSumOverGeneric[summation_?IsSummation, expr:{__}] :=
Module[{modifiedExpr = expr,
      masses,massDefine,codeMass,massRules,
      couplings,couplingDefine,codeCoupling,couplingRules,
      loopRules,loopArrayDefine,loopArraySet,
      cxxExpr,updatingVars, loopIdentifiers
   },
   masses = Cases[modifiedExpr, _SARAH`Mass, Infinity, Heads -> True];
   {massDefine, codeMass, massRules} = NameMasses@masses;
   modifiedExpr = modifiedExpr /. massRules;

   couplings = Cases[modifiedExpr, SARAH`Cp[__][___], Infinity, Heads -> True];
   {couplingDefine, codeCoupling, couplingRules} = NameCouplings@couplings;
   modifiedExpr = modifiedExpr /. couplingRules;

   modifiedExpr = modifiedExpr /. LoopTools`A0[args__] :> LoopTools`A0i[LoopTools`aa0, args];
   {loopRules, loopArrayDefine, loopArraySet, loopIdentifiers} = CreateLoopFunctions@modifiedExpr;
   modifiedExpr = modifiedExpr /. loopRules;

   cxxExpr = ApplyCXXRules/@modifiedExpr;
   updatingVars = MapThread[#1<>" += "<>#2<>";"&, {First/@GetBasis[], cxxExpr}];

   TextFormatting`ReplaceCXXTokens["
      // Shorter aliases for large types
      @fieldAliases@
      // Aliases for loop function identifiers
      @loopIdentifiers@

      @defineMasses@
      @defineCouplings@
      @defineLoopFunctions@
      // Start of summation over generic fields.
      @BeginSum@

         @setMasses@
         @setCouplings@
         @skipZeroAmplitude@
         @setLoopFunctions@

         @ChangeOutputValues@

      @EndSum@",
      {
         "@loopIdentifiers@" -> loopIdentifiers,
         "@fieldAliases@" -> CXXSetFieldAliases@couplings,
         "@defineMasses@" -> massDefine,
         "@defineCouplings@" -> couplingDefine,
         "@defineLoopFunctions@" -> loopArrayDefine,
         "@BeginSum@" -> CXXGenericSumBegin@summation,
         "@setMasses@" -> codeMass,
         "@setCouplings@" -> codeCoupling,
         "@skipZeroAmplitude@" -> CXXSkipZeroAmplitude[modifiedExpr, loopRules, massRules],
         "@setLoopFunctions@" -> loopArraySet,
         "@ChangeOutputValues@" -> Utils`StringJoinWithSeparator[updatingVars,"\n"],
         "@EndSum@" -> CXXGenericSumEnd@GetGenericFields@summation
      }
   ]
];
CXXSumOverGeneric::errUnimplementedLoops = "Loop functions `1` are not supported.";
CXXSumOverGeneric // secure;

CXXSetFieldAliases[couplings:{SARAH`Cp[__][___]..}] :=
Module[{l = Sort@DeleteDuplicates[(Sequence@@Head@#&)/@couplings], func},
   func = If[#1 =!= #2, "using " <> #1 <> " = " <> #2 <> ";\n", ""]&;
   StringJoin@MapThread[func, {CXXShortFieldName/@l, CXXFullFieldName/@l}]
];
CXXSetFieldAliases // secure;

Module[{c, name},
c[f_] := Switch[Head@f, SARAH`bar|Susyno`LieGroups`conj, "_"<>#, _, #]&;
name[f_] := ToString@Switch[f, _Symbol, f, _, Head@f];

CXXShortFieldName[f_?IsParticle] :=
   c[f][name@CXXDiagrams`RemoveLorentzConjugation@f];
CXXShortFieldName[f_?IsGenericParticle] :=
   c[f][CXXFullFieldName@CXXDiagrams`RemoveLorentzConjugation@f];
CXXShortFieldName // secure;
];

CXXLoopArray[num_Integer] := "l"<>ToString@num;
CXXLoopArray // secure;

CXXLoopType[sym_Symbol] :=
sym /. {
   LoopTools`A0i -> "looplibrary::Acoeff_t",
   LoopTools`B0i -> "looplibrary::Bcoeff_t",
   LoopTools`C0i -> "looplibrary::Ccoeff_t",
   LoopTools`D0i -> "looplibrary::Dcoeff_t"
};
CXXLoopType // secure;

CXXLoopLibrary[sym_Symbol] :=
sym /. {
   LoopTools`A0i -> "A",
   LoopTools`B0i -> "B",
   LoopTools`C0i -> "C",
   LoopTools`D0i -> "D"
};
CXXLoopLibrary // secure;

CreateLoopFunctions[modifiedExpr:{__}] :=
Module[{
      all, gathered, using, loopIdentifiers, abbreviate, abbreviations, loopRules,
      define, loopArrayDefine, set, loopArraySet
   },
   all = DeleteDuplicates@Cases[modifiedExpr, h_@__ /; Context@h == "LoopTools`", Infinity];
   gathered = GatherBy[all, # /. _[_, a___] :> {a}&];

   using[i_List] := i /. h_[id_, __] /; h =!= List :> "using looplibrary::" <> SymbolName@id <> ";";
   loopIdentifiers = StringRiffle[DeleteDuplicates@Flatten@Map[using, gathered], "\n"];

   abbreviate[i_List, {num_}] := i /. h_[id_, __] /; h =!= List :> CXXLoopArray@num <> "[" <> SymbolName@id <> "]";
   abbreviations = MapIndexed[abbreviate, gathered];
   loopRules = MapThread[Rule, Flatten/@{gathered, abbreviations}];

   define[{h_@__, ___}, {num_}] := CXXLoopType@h <> " " <> CXXLoopArray@num <> " {};";
   loopArrayDefine = StringRiffle[MapIndexed[define, gathered], "\n"];

   set[{h_[_, args__], ___}, {num_}] := Parameters`ExpressionToString[StringJoin["Loop_library::get().", CXXLoopLibrary@h][CXXLoopArray@num, args, "Sqr(context.scale())"]] <> ";";
   loopArraySet = Utils`StringJoinWithReplacement[MapIndexed[set, gathered], "\n", "\""->""];

   {loopRules, loopArrayDefine, loopArraySet, loopIdentifiers}
];
CreateLoopFunctions // secure;

CXXSkipZeroAmplitude[
   modifiedExpr:{__},
   loopRules:{Rule[_,_]...},
   massRules:{Rule[_,_]..}
] :=
Module[{massesToOne = Rule[#,1] & /@ massRules[[All,2]],
      loopsToOne = {},
      result, func = "z[" <> # <> "]"&
   },
   If[0 < Length@loopRules,
      loopsToOne = Rule[#,1] & /@ loopRules[[All,2]]
   ];
   result = RemoveNumbers[ExpandAll[modifiedExpr]] /.Plus->List/.
      massesToOne/.loopsToOne;
   result = DeleteDuplicates[Flatten@DeleteCases[result,1]];
   If[2 === LeafCount@result,
      result = func@@result,
      result = Plus@@(result/.HoldPattern[Times[x__]]:>Times@@(func@#&/@{x}))
   ];
   result = "if( "<>
      StringReplace[Parameters`ExpressionToString@result,"\""->""]<>
      " == 0 ) continue;";
   StringReplace[result,
      RegularExpression["g(\\d+)"]:> ToString[ToExpression@"$1"-1]
   ]
];
CXXSkipZeroAmplitude // secure;

CXXMassAbbreviation[(m_/; m === SARAH`Mass)[obj_?IsGenericParticle]] :=
Switch[GetParticleName@obj,
   GenericS, "mS",
   GenericF, "mF",
   GenericV, "mV",
   GenericU, "mU"
]<>ToString@GetIndex@obj;

CXXMassAbbreviation[(m_/; m === SARAH`Mass)[obj_?IsParticle]] := "m"<>ToString@GetParticleName@obj<>CXXIndex@obj;

CXXMassAbbreviation // secure;

NameMasses[masses:{__}] :=
Module[{info, i, r, sortedMasses, abbr, full},
   sortedMasses = Sort@DeleteDuplicates@masses;
   abbr = CXXMassAbbreviation /@ sortedMasses;
   full = ApplyCXXRules /@ sortedMasses;
   {
      "double "<>Utils`StringJoinWithSeparator[abbr, ", "]<>";",
      StringJoin@MapThread[#1<>" = "<>#2<>";\n"&, {abbr, full}],
      MapThread[Rule, {sortedMasses, abbr}]
   }
];
NameMasses // secure;

NameCouplings[couplings:{__}] :=
Module[{sortedCpl, abbr, full, initCpl, initZero},
   sortedCpl = Sort@DeleteDuplicates@couplings;
   abbr = Table["g"<>ToString@n, {n, Length@sortedCpl}];
   full = ApplyCXXRules/@sortedCpl;

   initCpl = StringJoin@MapThread[#1<>" = "<>#2<>";\n"&, {abbr, full}];
   initZero = Table["z[" <> ToString[n-1] <> "] = (std::abs(g" <> ToString@n <> ") < std::numeric_limits<double>::epsilon())?0:1;\n", {n, Length@sortedCpl}];


   info = {"g"<>ToString@#, ApplyCXXRules@sort[[#,1]], sort[[#,1]],
      sort[[#,2]]}&;
   {
      "std::complex<double> "<>Utils`StringJoinWithSeparator[abbr, ", "]<>";\nstd::array<int, "<> ToString@Length@sortedCpl <> "> z{};",
      initCpl <> "\n" <> initZero,
      MapThread[#1 -> I*#2&, {sortedCpl, abbr}]
   }
];
NameCouplings // secure;

CXXGenericAbbreviations[genFields:{__?IsGenericParticle}] :=
Utils`StringJoinWithSeparator[
   Apply[
      "using "<>#1<>" = typename at<GenericFieldMap,"<>ToString@#2<>">::type;"&,
      {CXXFullFieldName@#, CXXGenericKey@#} &/@ genFields,
      {1}
   ],
   "\n"
];
CXXGenericAbbreviations // secure;

CXXGenericSumBegin[summation_?IsSummation]:=
Module[{beginsOfFor},
   beginsOfFor = "for( const auto &"<>CXXIndex[#[[1]]]<>" : "<>
      "index_range<"<> CXXFullFieldName[#[[1]]]<>">() ) {\n"<>
      "at_key<"<>CXXGenericKey@#[[1]]<>">( index_map ) = "<>
      CXXIndex[#[[1]]]<>";"<>parseRestrictionRule[#] &/@summation;
   Utils`StringJoinWithSeparator[beginsOfFor, "\n"]
];
CXXGenericSumBegin // secure;

parseRestrictionRule[{genericField_?IsGenericParticle,rule_}] :=
Module[{f1,f2,GetIndexOfExternalField,OrTwoDifferent},
   GetIndexOfExternalField[_[_[{ind_}]]] := ApplyCXXRules@ind;
   GetIndexOfExternalField[_[{ind_}]] := ApplyCXXRules@ind;
   GetIndexOfExternalField[_] := "i0";
   OrTwoDifferent[] :=
   Module[{type1 = CXXFullFieldName@First@rule,
         type2 = CXXFullFieldName@Last@rule,
         ind = GetIndexOfExternalField@First@rule,
         typeGen = CXXFullFieldName@genericField,
         indGen = CXXIndex@genericField},
      "\nif( (std::is_same_v<"<>typeGen<>","<>type1<>
         "> || std::is_same_v<"<>typeGen<>","<>type2<>
         ">) && "<>indGen<>" == "<>ind<>" ) continue;"];
   Switch[rule,
      Or[f1_,f2_],
         OrTwoDifferent[],
      False,
         "",
      _,
         "@todo This rule is not implemented yet!";Quit[1]]];
parseRestrictionRule // secure;

CXXGenericSumEnd[genFields:{__?IsGenericParticle}] :=
StringJoin[
   Array["}"&,Length@genFields],
   " // End of summation over generic fields"
];
CXXGenericSumEnd // secure;

CXXCalculateFunction[genSumNames:{__String}] :=
Module[{
      varName = "genericsum" (* Feel free to change me to another C++ name *),
      varNames,initVars,sumOfSums},
   varNames = Array[varName<>ToString@#&,Length@genSumNames];
   initVars = Utils`StringJoinWithSeparator[
      MapThread["const auto "<>#1<>" = "<>#2<>"();"&, {varNames, genSumNames}],
      "\n"
   ];
   sumOfSums = Utils`StringJoinWithSeparator[#<>"[i]"&/@varNames, "+"];
   TextFormatting`ReplaceCXXTokens["
      std::array<std::complex<double>,@BasisLength@> calculate( void ) {
         std::array<std::complex<double>,@BasisLength@> genericSummation;
         constexpr int coeffsLength = genericSummation.size();
         @InitializeVariablesWhichStoreGenericSumsOutput@

         for ( std::size_t i=0; i<coeffsLength; i++ ) {
            genericSummation[i] += @SumOfVariables@;
         }
         return genericSummation;
      } // End of calculate()",
      {
         "@BasisLength@"->GetLength@GetBasis[],
         "@InitializeVariablesWhichStoreGenericSumsOutput@"->initVars,
         "@SumOfVariables@"->sumOfSums
      }
   ]
];
CXXCalculateFunction // secure;

CXXInsertFields[genInsertions_?IsClassFields] :=
Utils`StringJoinWithSeparator["boost::mpl::vector<"<>
   Utils`StringJoinWithSeparator[CXXFullFieldName@#&/@#, ", "]<>
   ">"&/@genInsertions,
   ",\n"
];
CXXInsertFields // secure;

CXXInsertFactors[combinatorialFactors_?IsCombinatoricalFactors] :=
Utils`StringJoinWithSeparator["boost::mpl::int_<"<>ToString@#<>
   ">"&/@combinatorialFactors,
   ",\n"
];
CXXInsertFactors // secure;

CXXInsertColors[colorFactors:{__?NumberQ}] :=
Module[{
      ReRatioColorFactors = {Numerator@#, Denominator@#} &/@ Re@colorFactors,
      ImRatioColorFactors = {Numerator@#, Denominator@#} &/@ Im@colorFactors
   },
   Utils`StringJoinWithSeparator[
      StringReplace[
         MapThread[
            "detail::complex_helper<"<>
            "detail::ratio_helper<"<>ToString@#1<>">,"<>
            "detail::ratio_helper<"<>ToString@#2 <> ">>"&,
            {ReRatioColorFactors, ImRatioColorFactors}
         ],
         {"{" -> "", "}" -> ""}
      ],
   ",\n"
   ]
];
CXXInsertColors // secure;

End[];
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "NPointFunctions`"];
