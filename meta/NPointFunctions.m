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

(* Functions *)
VerticesForNPointFunction;
CreateCXXHeaders;
CreateCXXFunctions;
GetProcess;
GetSubexpressions;
GetGenericSums;
GetParticleName;
NPointFunction;
IsNPointFunction;

(* Symbols *)
SetAttributes[
   {
      LoopLevel, Regularize, UseCache, ZeroExternalMomenta, OnShellFlag,
      OperatorsOnly, ExceptLoops, KeepProcesses, Irreducible, Triangles,
      GenericS, GenericF, GenericV, GenericU, GenericSum, GenericIndex,
      LorentzIndex, Mat, DiracChain
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

GetCXXIndex[obj_?IsGenericParticle] :=
   "i"<>StringTake[SymbolName[obj[[0]]],-1]<>ToString[obj[[1,1]]] &@ conj@obj;
GetCXXIndex[obj:_?IsParticle] :=
Module[{nakedField},
   nakedField = obj /. {SARAH`bar->Identity, Susyno`LieGroups`conj->Identity};
   Switch[nakedField,
   _Symbol,              "",
   (_Symbol)[{_Symbol}], StringDrop[ToString[nakedField[[1, 1]]],2]]
];
GetCXXIndex // secure;

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

GetParticleName[particle_?IsParticle] := particle /. {_@n_@_List :> n, n_@_List :> n, _@n_ :> n};
GetParticleName[particle_?IsGenericParticle] := particle /. {_@n_@_@_ :> n, n_@_ :> n};
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
      {}];
   WriteString[$Output, "done\n"];
   result
];
ApplySubexpressions // secure;

conj[obj_?IsGenericParticle] :=
Switch[Head@obj,
   SARAH`bar | Susyno`LieGroups`conj, First@obj,
   GenericS | GenericV,               Susyno`LieGroups`conj@obj,
   GenericF | GenericU,               SARAH`bar@obj
];
conj[obj_?TreeMasses`IsScalar|_?TreeMasses`IsVector] := Susyno`LieGroups`conj@obj;
conj[obj_?TreeMasses`IsFermion] := SARAH`bar@obj;
conj // secure;

OutputPaths[] := Module[{outputPath, eigenstates = ToString@FlexibleSUSY`FSEigenstates},
   outputPath = {SARAH`$sarahCurrentOutputMainDir, eigenstates};
   FileNameJoin@*Flatten /@
   {
      {outputPath, "NPointFunctionis"},
      {outputPath, "FormCalc"},
      {outputPath, "FeynArts", SARAH`ModelName<>eigenstates},
      {outputPath, "FeynArts", "ParticleNamesFeynArts.dat"},
      {outputPath, "FeynArts", "ParticleNamespaces.m"}
   }
];

Options[NPointFunction] =
{
   LoopLevel -> 1,
   Regularize -> FlexibleSUSY`FSRenormalizationScheme,
   UseCache -> True,
   ZeroExternalMomenta -> True,
   OnShellFlag -> True,
   KeepProcesses -> {},
   Observable -> None
};

NPointFunction[inFields_List, outFields_List, opts:OptionsPattern[]] :=
Module[{
      npfDir, fcDir, faModel, particleNamesFile,
      particleNamespaceFile, subKernel, currentDirectory,
      nPointFunction
   },
   {npfDir, fcDir, faModel, particleNamesFile, particleNamespaceFile} = OutputPaths[];
   If[!DirectoryQ@npfDir,
      CreateDirectory@npfDir
   ];
   If[OptionValue@UseCache,
      nPointFunction = CachedNPointFunction[
         inFields,outFields,npfDir,
         {Join[inFields, outFields], OptionValue@KeepProcesses}];
      If[nPointFunction =!= Null, Return@nPointFunction];
   ];
   If[!FileExistsQ[faModel <> ".mod"],
      subKernel = LaunchSubkernelFor@"creation of FeynArts model file";
      makeModelFile@subKernel;
      writeNamespaceFile@particleNamespaceFile;
      CloseKernels@subKernel;
   ];

   subKernel = LaunchSubkernelFor@"FormCalc code generation";
   SetSharedFunction[subWrite, Print];
   With[{path = $Path,
         data = {fcDir, faModel, particleNamesFile,
            particleNamespaceFile,
            renameFields[inFields, particleNamesFile, "SARAH"->"FeynArts"],
            renameFields[outFields, particleNamesFile, "SARAH"->"FeynArts"]},
         options = {OptionValue@Observable,
            OptionValue@LoopLevel,
            SymbolName/@If[List=!=Head@#, {#}, #]&@OptionValue@KeepProcesses,
            OptionValue@ZeroExternalMomenta,
            OptionValue@OnShellFlag,
            OptionValue@Regularize},
         meta = FlexibleSUSY`$flexiblesusyMetaDir},
      nPointFunction = RemoveEmptyGenSums@ParallelEvaluate[
         $Path = path;
         Get@FileNameJoin@{meta, "NPointFunctions", "internal.m"};
         NPointFunction[data, options],
         subKernel,
         DistributedContexts -> None];];
   CloseKernels@subKernel;

   UnsetShared[subWrite, Print];

   cache[nPointFunction, npfDir,
      {Join[inFields, outFields], OptionValue@KeepProcesses}];
   nPointFunction] /; inputCheck[inFields,outFields,opts];
NPointFunction::errFields = "
The field '`1`' must belong to (bar or conj can be applied):
   `2`.";
NPointFunction::errOptions = "
Unknown option(s):
   `1`.
Currently supported options are:
   `2`.";
NPointFunction::errLoopLevel = "
Currently loop levels <= 1 are supported.";
NPointFunction::errRegularize = "
Unknown regularization scheme `1`.
Supported schemes:
   FlexibleSUSY`.`DRbar,
   FlexibleSUSY`.`MSbar.";
NPointFunction::errUseCache=
"UseCache must be either True or False.";
NPointFunction::errZeroExternalMomenta=
"ZeroExternalMomenta must be True, False, ExceptLoops, OperatorsOnly.";
NPointFunction::errOnShellFlag=
"OnShellFlag must be either True or False.";
NPointFunction // secure;

subWrite::usage = "
@brief Prints a string.
@note ``SetSharedFunction`` does not cause names leaking.";
subWrite[str_String] := WriteString[$Output, str];
subWrite // secure;

inputCheck[inFields:{__},outFields:{__},opts___] :=
Module[{aoq, ip, allowedParticles, options, unknown},
   aoq = Utils`AssertOrQuit;
   ip = TreeMasses`IsParticle;
   (* TODO add |_?TreeMasses`IsVector *)
   allowedParticles = Cases[
      TreeMasses`GetParticles[],
      _?TreeMasses`IsScalar|_?TreeMasses`IsFermion];
   aoq[ip@#, NPointFunction::errFields, #, allowedParticles] &/@ inFields;
   aoq[ip@#, NPointFunction::errFields, #, allowedParticles] &/@ outFields;
   options = Options[NPointFunction][[All,1]];
   unknown = FilterRules[{opts}, Except@options];
   aoq[unknown === {}, NPointFunction::errOptions, unknown, options];
   (*Now we know that all options are iside allowed list.*)
   Cases[{opts}, Rule[LoopLevel, x_] :>
      aoq[0<=x<=1, NPointFunction::errLoopLevel]];
   Cases[{opts}, Rule[Regularize, x_] :>
      aoq[MemberQ[{FlexibleSUSY`MSbar, FlexibleSUSY`DRbar}, x],
         NPointFunction::errRegularize, x]];
   Cases[{opts}, Rule[UseCache, x_] :>
      aoq[x===True || x===False, NPointFunction::errUseCache]];
   Cases[{opts}, Rule[ZeroExternalMomenta,x_] :>
      aoq[MemberQ[{True, False, OperatorsOnly, ExceptLoops}, x],
         NPointFunction::errZeroExternalMomenta]];
   Cases[{opts}, Rule[OnShellFlag, x_] :>
      aoq[x===True || x===False, NPointFunction::errOnShellFlag]];
   True];
inputCheck // secure;

VerticesForNPointFunction::usage = "
@brief Return a list of all vertices needed to calculate a given
       n-point correlation function.
@param obj The n-point correlation function.
@returns A list of all vertices needed for the calculation.";
VerticesForNPointFunction[obj_?IsNPointFunction] :=
Module[{v, getVertex},
   getVertex[vertGen_, rules_] := vertGen/.#&/@rules;
   v = DeleteDuplicates@
      Cases[#, SARAH`Cp[f__] :> {f}, Infinity, Heads->True] &/@
         GetGenericSums@obj;
   DeleteDuplicates[
      Vertices`StripFieldIndices/@#&/@
         Flatten[MapThread[getVertex, {v, GetClassFieldRules@obj}],2]]];
VerticesForNPointFunction // secure;

GetSARAHModelName::usage = "
@returns The ``SARAH`` model name as to be passed to ``SARAH`Start``.";
GetSARAHModelName[] :=
If[SARAH`submodeldir =!= False,
      SARAH`modelDir <> "-" <> SARAH`submodeldir,
      SARAH`modelDir];
GetSARAHModelName // secure;

LaunchSubkernelFor::usage = "
@brief Tries to launch a subkernel without errors.
@param message String with a reason to launch a subkernel.
@returns A subkernels name.";
LaunchSubkernelFor[message_String] :=
Module[{kernelName},
   kernelName = LaunchKernels@1;
   If[Head@kernelName === List, kernelName[[1]], kernelName]];
LaunchSubkernelFor // secure;

CacheNameForMeta::usage = "
@param meta The given meta information.
@returns The name of the cache file for given meta information.";
CacheNameForMeta[meta:{__}] :=
   Utils`StringJoinWithSeparator[Flatten@meta, "", SymbolName]<>".m";
CacheNameForMeta // secure;

cache::usage = "
@brief Writes a given n-point function to the cache.
@param new The given n-point function.
@param dir The directory to save cache.
@param meta The meta information about the given n-point function.";
cache[new_, dir_, meta:{__}] :=
Module[{path, file, npfs, position},
   path = FileNameJoin@{dir, CacheNameForMeta@meta};
   If[FileExistsQ@path,
      npfs = Get@path;,
      npfs = {};];
   position = Position[npfs[[All,1]], new[[1]]];
   If[Length@position === 1,
      npfs[[position[[1]]]] = new;,
      AppendTo[npfs, new];];
   file = OpenWrite@path;
   Write[file, npfs];
   Close@file;];
cache // secure;

CachedNPointFunction::usage = "
@brief Retrieve an n-point correlation function from the cache.
@param inFields the incoming fields of the n-point correlation function.
@param outFields the outgoing fields of the n-point correlation function.
@param cacheDir the directory to save cache.
@param nPointMeta the meta information of the n-point correlation function.
@returns The n-point function from cache or ``Null`` if the function is missing.";
CachedNPointFunction[inFields_,outFields_,cacheDir_,nPointMeta:{__}] :=
Module[{nPointFunctionsFile, nPointFunctions, position},
   nPointFunctionsFile = FileNameJoin@{cacheDir,CacheNameForMeta@nPointMeta};
   If[!FileExistsQ@nPointFunctionsFile,Return@Null];
   nPointFunctions = Get@nPointFunctionsFile;
   position = Position[Vertices`StripFieldIndices[ nPointFunctions[[All,1]] ],
      {inFields, outFields}];
   If[Length@position == 1,nPointFunctions[[ position[[1,1]] ]],Null]];
CachedNPointFunction // secure;

makeModelFile::usage = "
@brief Generate the ``FeynArts`` model file on a given subkernel.";
makeModelFile[kernel:_Parallel`Kernels`kernel|_KernelObject] :=
With[{currentPath = $Path,
      currentDir = Directory[],
      fsMetaDir = FlexibleSUSY`$flexiblesusyMetaDir,
      sarahInputDirs = SARAH`SARAH@SARAH`InputDirectories,
      sarahOutputDir = SARAH`SARAH@SARAH`OutputDirectory,
      SARAHModelName = GetSARAHModelName[],
      eigenstates = FlexibleSUSY`FSEigenstates},
   Print["Generating FeynArts model file ..."];
   SetSharedFunction[Print];
   ParallelEvaluate[
      $Path = currentPath;
      SetDirectory@currentDir;
      Get@FileNameJoin@{fsMetaDir, "NPointFunctions", "MakeModelFile.m"};
      NPointFunctions`MakeModelFile[sarahInputDirs,sarahOutputDir,
         SARAHModelName, eigenstates];,
      kernel, DistributedContexts -> None];
   UnsetShared[Print];
   Print["Generating FeynArts model file ... done"];];
makeModelFile // secure;

writeNamespaceFile::usage = "
@brief Write a file containing all field names and the contexts in which they
       live in ``Mathematica``.
@note This is necessary because ``SARAH`` puts fields into different
      contexts.";
writeNamespaceFile[fileName_String] :=
Module[{fileHandle = OpenWrite@fileName},
   Write[fileHandle, {ToString@#, Context@#} & /@ TreeMasses`GetParticles[]];
   Close@fileHandle;];
writeNamespaceFile // secure;

renameFields[fields_, particles_String, "SARAH" -> "FeynArts"] :=
   Module[{unique, faNames},
      unique = DeleteDuplicates[
         CXXDiagrams`RemoveLorentzConjugation[#]&/@ fields];
      faNames = Flatten[
         StringCases[Utils`ReadLinesInFile@particles,
            ToString@# ~~ ": " ~~ x__ ~~ "]" ~~ ___ :>
               "FeynArts`" <> x <> "]"]&/@ unique];
      fields /. MapThread[Rule, {unique, faNames}] /.
         {  SARAH`bar@field_String :> "-" <> field,
            Susyno`LieGroups`conj@field_String :> "-" <> field}];
renameFields // secure;

RemoveEmptyGenSums::usage = "
@brief Somehing went wrong in a subkernel. Most likely, it is dead.
@param npfObject n-point function object to clean.
@returns Cleaned from empty GenericSums npfObject.";
RemoveEmptyGenSums[npfObject_?IsNPointFunction]:=npfObject;
RemoveEmptyGenSums[
   {  fields:{{__},{__}},
      {  {  sums:{GenericSum[_,{___}]..},
            rules:{{{__}..}..},
            comb:{{__Integer}..},
            col:{{__}..}},
         subs:{Rule[_,_]...}}}]:=
Module[{poss=Position[sums,GenericSum[{0},{}]]},
   Print["Removing zero GenericSum at positions ",
      Utils`StringJoinWithSeparator[Flatten@poss,", "],"."];
   {fields,{Delete[#,poss]&/@{sums,rules,comb,col},subs}}];
RemoveEmptyGenSums // secure;

CreateCXXHeaders::usage = "
@brief Create the C++ code for the necessary headers.
@returns The C++ code for the necessary headers.";
CreateCXXHeaders[] :=
TextFormatting`ReplaceCXXTokens["
   #include \"loop_libraries/loop_library.hpp\"
   #include \"cxx_qft/@ModelName@_npointfunctions_wilsoncoeffs.hpp\"
   #include \"concatenate.hpp\"
   #include <limits>
   #include <boost/fusion/include/at_key.hpp>
   #include <boost/core/is_same.hpp>",
   {"@ModelName@"->FlexibleSUSY`FSModelName}];
CreateCXXHeaders // secure;

CreateCXXFunctions::usage = "
@brief Creates a C++ code for the n-point correllation functions.
@param npf n-point correlation function object.
@param name function name
@param colourProjector A function that represents color structure of ampitudes.
@param wilsonBasis Basis for matching.
@returns A list with C++ prototypes and definitions.";
CreateCXXFunctions[
   npf_?IsNPointFunction,
   name:_String,
   colourProjector:_?IsColorProjector,
   wilsonBasis:_?IsWilsonBasis:{"value"->"dummy string"}] :=
Module[{mainFunction, prototype, definition},
   mainFunction = "std::array<std::complex<double>,"<>
      GetLength@wilsonBasis<>"> "<>name<>"("<>#<>")"&;
   setHelperClassName@npf;
   setBasis@wilsonBasis;
   prototype = mainFunction@`cxx`arguments[npf,Default]<>";";
   definition =
      `cxx`npfClass[ApplySubexpressions@npf,colourProjector] <> "\n\n" <>
      mainFunction@`cxx`arguments@npf <>
      "{\n   "<>$helperClassName<>
      " helper{ model, indices, momenta };\n   return helper.calculate();\n}";
   {prototype, definition}] /; And@@
      (Utils`AssertOrQuit[Length@wilsonBasis === Length@#,
         CreateCXXFunctions::errNoMatch]&/@GetSumExpression/@GetGenericSums@npf);
CreateCXXFunctions::errNoMatch = "
Length of basis and the given n-point function one does not match."
CreateCXXFunctions // secure;

`cxx`arguments::usage = "
@brief Returns the C++ arguments that the C++ version of the given n-point
       correlation function shall take.
       ``Default`` value of zero for all external momenta is chosen if the
       second parameter is ``Default``.
@param npf The given n-point correlation function
@param control String that sets up the type of argument string
@returns The C++ arguments that the C++ version of the given n-point
         correlation function shall take.";
`cxx`arguments[npf_?IsNPointFunction,control:Null|Default:Null] :=
   "const "<>#1<>" &model,"<>
   " const std::array<int,"<>#2<>"> &indices,"<>
   " const std::array<Eigen::Vector4d,"<>#3<>"> &momenta"<>
   If[control === Default," = { "<>#4<>" }",""]&[
      FlexibleSUSY`FSModelName<>"_mass_eigenstates",
      ToString@Length@GetExternalIndices@npf,
      ToString@#,
      Utils`StringJoinWithReplacement@Array["Eigen::Vector4d::Zero()",#]] &@
         Length@GetExternalMomenta@npf;
`cxx`arguments // secure;

$basis = {"value"->"dummy string"};
$basis ~ SetAttributes ~ {Protected};

setBasis[obj:_?IsWilsonBasis] :=
(  Unprotect@$basis;
   $basis = obj;
   Protect@$basis;);
setBasis // secure;

$helperClassName = "";
$helperClassName ~ SetAttributes ~ {Protected};

setHelperClassName::usage = "
@brief Sets the C++ name for the helper class of a n-point function.
@param obj n-point function object.";
setHelperClassName[obj_?IsNPointFunction] :=
Module[{fieldNames = Vertices`StripFieldIndices/@Join@@GetProcess[obj]},
   Unprotect@$helperClassName;
   $helperClassName = "nPoint" <>
      StringJoin@Map[ToString,fieldNames/.a_[b_]:>Sequence@@{a,b}] <> "_" <>
      ToString@Ceiling[10^6*AbsoluteTime[]];
   Protect@$helperClassName;];
setHelperClassName // secure;

`cxx`npfClass::usage = "
@param npf The given n-point correlation function.
@param projCol The colour factor projection to be applied for the
       given n-point correlation function.
@returns The C++ code for the helper class of n-point function.";
`cxx`npfClass[npf_?IsNPointFunction, projCol:_?IsColorProjector] :=
Module[{genSums, extIndices, numberOfMomenta, genFields, genSumNames},
   genSums = GetGenericSums@npf;
   extIndices = GetExternalIndices@npf;
   numberOfMomenta = Length[GetExternalMomenta@npf];
   genFields = DeleteDuplicates[Flatten@GetClassFieldRules@npf /. Rule[x_,_] :> x];
   code = "
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
   }; // End of @ClassName@";
   genSumNames = Array["genericSum"<>ToString@#&,Length@genSums];
   `cxx`setRules[extIndices,genFields];
   TextFormatting`ReplaceCXXTokens[code,
      {  "@ClassName@"->$helperClassName,
         "@Context@"->"correlation_function_context<"<>
            ToString@Length@extIndices<>","<>ToString@numberOfMomenta<>">",
         "@KeyStructsInitialization@"->`cxx`initializeKeyStructs@genFields,
         "@GenericSums@"->`cxx`genericSum[npf,projCol,genSumNames],
         "@Arguments@"->`cxx`arguments@npf,
         "@CalculateFunction@"->`cxx`calculateFunction@genSumNames}]];
`cxx`npfClass // secure;

`cxx`initializeKeyStructs::usage = "
@brief Generates required C++ code for key structs initialization.
@param fields List of generic fields.
@returns C++ code for subexpression if generic fields present there.";
`cxx`initializeKeyStructs[fields:{__?IsGenericParticle}]:=
   Utils`StringJoinWithSeparator[
      "struct "<>#<>" {};"&/@`cxx`genericFieldKey/@fields,
      "\n"];
`cxx`initializeKeyStructs // secure;

Module[{rules = {{}}},

`cxx`setRules::usage = "
@brief Generate a list of rules for translating Mathematica expressions to
       C++ ones.
@param extIndices The external indices of an n-point correlation function.
@param genericFields The generic fields appearing in an n-point correlation
       function.
@returns A list of rules for translating Mathematica expressions to C++ ones.
@note All couplings have to be multiplied by ``I``.";
`cxx`setRules[extIndices:{___Symbol},genericFields:{__?IsGenericParticle}] :=
Module[{externalIndexRules, wrap, index, genericRules, massRules, couplingRules},
   externalIndexRules = MapThread[Rule,
      {  extIndices,
         Array["i"<>ToString[#]&,Length@extIndices]}];
   genericRules=Flatten[Thread@Rule[
      {conj[#],#},
      {  `cxx`fieldName[conj@#][GetCXXIndex@#],
         `cxx`fieldName[#][GetCXXIndex@#]}] &/@ genericFields];
   wrap[fields__] := Utils`StringJoinWithSeparator[`cxx`fieldAlias/@{fields},", "];
   index[fields__] := Utils`StringJoinWithSeparator[`cxx`fieldIndices/@{fields},", "];

   couplingRules =
   {  SARAH`Cp[fields__][1] :>
         ("NPF_S(" <> wrap@fields <>") NPF_I(" <> index@fields <> ")"),
      SARAH`Cp[fields__][SARAH`PL] :>
         ("NPF_L(" <> wrap@fields <>") NPF_I(" <> index@fields <> ")"),
      SARAH`Cp[fields__][SARAH`PR] :>
         ("NPF_R(" <> wrap@fields <>") NPF_I(" <> index@fields <> ")"),
      (* On C++ level this vertex is such, that _D should have 0 and 1.
       * The first position corresponds to the momenta with + sign.
       * The second position corresponds to the momenta with - sign.
       * Example: Cp[S1, S2, V][Mom@S1 - Mom@S2] leads to
       *          NPF_MD(S1, S2, V) NPF_D(0, 1) <indices>
       * Example: Cp[S1, S2, V][Mom@S2 - Mom@S1] leads to
       *          NPF_MD(S1, S2, V) NPF_D(1, 0) <indices>*)
      SARAH`Cp[fields__][SARAH`Mom[f1_] - SARAH`Mom[f2_]] :>
      StringReplace["NPF_MD(`1`) NPF_D(`2`, `3`) NPF_I(`4`)",
         {  "`1`"->ToString@wrap@fields,
            "`2`"->ToString[First@@Position[{fields},f1,{1}]-1],
            "`3`"->ToString[First@@Position[{fields},f2,{1}]-1],
            "`4`"->ToString@index@fields}],
      SARAH`Cp[fields__][SARAH`g[_, _]] :>
         ("NPF_G(" <> wrap@fields <>") NPF_I(" <> index@fields <> ")"),
      SARAH`Cp[fields__][(SARAH`Mom[f2_, _]-SARAH`Mom[f1_, _])*SARAH`g[_,_],
         (SARAH`Mom[f1_,_]-SARAH`Mom[f3_,_])*SARAH`g[_,_],
         (SARAH`Mom[f3_,_]-SARAH`Mom[f2_,_])*SARAH`g[_,_]] :>
         ("NPF_T(" <> wrap@fields <>") NPF_I(" <> index@fields <> ")")};

   massRules = {SARAH`Mass[f_] :>
      StringJoin["context.mass<",`cxx`fieldAlias@f,">(",
         `cxx`fieldIndices@f,")"]};

   rules = {externalIndexRules, massRules, couplingRules, genericRules};];
`cxx`setRules // secure;

`cxx`applyRules[obj_] :=
StringReplace[
   Parameters`ExpressionToString[Fold[ReplaceAll, obj, rules]],
   "\"" -> ""];
`cxx`applyRules // secure;];

Module[{strip, sandwich, n},

strip[f_] := f /. {SARAH`bar -> Identity, Susyno`LieGroups`conj -> Identity};

sandwich[f_] = Switch[Head@f,
   SARAH`bar, "typename cxx_diagrams::bar<" <> # <> ">::type",
   Susyno`LieGroups`conj, "typename cxx_diagrams::conj<" <> # <> ">::type",
   _, #]&;

`cxx`fieldName::usage = "
@brief Returns a C++ representation for a field expression or a field name.
@param f A generic or external field, or an explicit field name.
@returns A C++ representation for a field.";
`cxx`fieldName[f_?IsParticle] :=
(  n = strip@f;
   sandwich[f]["fields::" <> ToString@Switch[n, _Symbol, n, _, Head@n]]);

`cxx`fieldName[f_?IsGenericParticle] :=
(  n = strip@f;
   sandwich[f]["g"<>StringTake[ToString@Head@n, -1] <> ToString@Part[n, 1, 1]]);

`cxx`fieldName // secure;];

`cxx`fieldIndices::usage = "
@brief Return the C++ expression for the given field.
@param field The given field.
@returns The C++ expression for the given field.
@note Saves its previous calls to improve the speed.";
`cxx`fieldIndices[SARAH`bar[field_]] := `cxx`fieldIndices[SARAH`bar[field]] =
   `cxx`fieldIndices@field;
`cxx`fieldIndices[Susyno`LieGroups`conj[field_]] :=
   `cxx`fieldIndices[Susyno`LieGroups`conj[field]] =
      `cxx`fieldIndices@field;
`cxx`fieldIndices[head_[GenericIndex[index_Integer]]] :=
   `cxx`fieldIndices[head[GenericIndex[index]]] =
      "i"<>StringTake[SymbolName@head,-1]<>ToString@index;
`cxx`fieldIndices[field_] := `cxx`fieldIndices[field] =
   If[Length@field === 0, "i0", field[[1, 1]]];
`cxx`fieldIndices // Utils`MakeUnknownInputDefinition;
`cxx`fieldIndices ~ SetAttributes ~ {Locked};

`cxx`genericFieldKey::usage = "
@brief Determines C++ key type for a generic field.
@param genericField Given generic field.
@returns C++ key type of a generic field(s).";
`cxx`genericFieldKey[fields:{__}] :=
   Utils`StringJoinWithSeparator[fields, ", ", `cxx`genericFieldKey];
`cxx`genericFieldKey[head_[GenericIndex[index:_Integer]]] :=
   SymbolName@head<>ToString@index<>"Key";
`cxx`genericFieldKey // secure;

removeNumbers[expr_] :=
   expr /. {_Integer-> 1, _Rational-> 1, _Complex-> 1, Pi-> 1};
removeNumbers // secure;

colorFactor::usage = "
@brief Extracts the colour factor for a given colour structure.";
colorFactor[colourfactors:{__?IsColorFactors},
            projection:_?IsColorProjector] :=
Module[{projectedFactors},
   projectedFactors =
      Switch[projection,
         Identity,
            colourfactors,
         _,
            Utils`AssertOrQuit[
               And@@Flatten[MatchQ[#, _projection]&/@#&/@ removeNumbers@colourfactors],
               colorFactor::errMultipleColourStructures,
               colourfactors,
               projection];
            colourfactors /. _projection -> 1];
   Utils`AssertOrQuit[
      And@@Flatten[NumericQ[#]&/@#&/@ projectedFactors],
      colorFactor::errNotNumber,
      projectedFactors];
   projectedFactors];
colorFactor // secure;
colorFactor::errMultipleColourStructures = "
Not all colour factors in
   `1`
have head `2`.";
colorFactor::errNotNumber = "
Not all colour factors in
   `1`
are numbers.";

`cxx`genericSum::usage = "
@brief Create the C++ code form of a generic sums.
@param obj n-point function object.
@param colourProjector An expression, representing the projector.
@param genSumNames Set of names for generic sums.";
`cxx`genericSum[obj_?IsNPointFunction, colourProjector:_?IsColorProjector,
   genSumNames:{__String}] :=
Utils`StringJoinWithSeparator[
   MapThread[
      `cxx`genericSum[##,GetSubexpressions@obj,obj]&,
      {  GetGenericSums@obj,
         GetClassFields@obj,
         GetCombinatoricalFactors@obj,
         colorFactor[GetColorFactors@obj,colourProjector],
         genSumNames}],
   "\n\n"];
`cxx`genericSum[
   sum_?IsGenericSum,
   genericInsertions_?IsClassFields,
   combinatorialFactors_?IsCombinatoricalFactors,
   colourFactors:{__?NumericQ},
   genSumName_String,
   subexpressions_?IsSubexpressions,
   npf_?IsNPointFunction] :=
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
      using colour_factors = boost::mpl::vector<
         @ColorFactors@
         >;
      return accumulate_generic<
         GenericKeys,
         GenericInsertions,
         combinatorial_factors,
         colour_factors,
         boost::mpl::int_<@WilsonBasisLength@>,
         @GenericSum_NAME@_impl
         >( *this );
   } // End of function @GenericSum_NAME@()",
   {  "@GenericSum_NAME@"->genSumName,
      "@GenericFieldShortNames@"->`cxx`shortNames@GetGenericFields@sum,
      "@ExternalIndices@"->`cxx`initializeExternalIndices@npf,
      "@InitializeOutputVars@"->Utils`StringJoinWithSeparator[
         "std::complex<double> "<>#<>" = 0.0;"&/@First/@$basis,"\n"],
      "@SummationOverGenericFields@"->`cxx`changeGenericExpressions[
         GetSumParticles@sum,GetSumExpression@sum],
      "@ReturnOutputVars@"->ToString[First/@$basis],
      "@GenericKeys@"->`cxx`genericFieldKey@GetGenericFields@sum,
      "@ClassInsertions@"->`cxx`insertFields@genericInsertions,
      "@CombinatoricalFactors@"->`cxx`insertFactors@combinatorialFactors,
      "@ColorFactors@"->`cxx`insertColours@colourFactors,
      "@WilsonBasisLength@"->GetLength@$basis}];
`cxx`genericSum // secure;
`cxx`genericSum::errColours = "
Colour factor is not a number after projection: `1`";

`cxx`initializeExternalIndices[npf_?IsNPointFunction] :=
Module[{extIndices = GetExternalIndices@npf},
   indices = Array[
      "std::array<int, 1> i" <> ToString@# <>
         " {this->external_indices("<>ToString[#-1]<>")};\n"&,
      Length@extIndices];
   If[Length@extIndices < Length[Flatten@GetProcess@npf],
      PrependTo[indices, "std::array<int, 0> i0 {};\n"];];
   StringJoin@indices];
`cxx`initializeExternalIndices // secure;

`cxx`changeGenericExpressions::usage = "
@brief Generates C++ code for output value updating inside generic sum.
@param summation Summation structure of generic index.
@param expr List of expressions to be converted into C++ code.
@returns C++ code for output value initializations inside generic sum.";
`cxx`changeGenericExpressions::errUnimplementedLoops =
"Unsupported loop functions
   `1`
were detected.";
`cxx`changeGenericExpressions[summation_?IsSummation, expr:{__}] :=
Module[{
      code = "
      // Shorter aliases for large types
      @fieldAliases@
      // These definitions are repeated multiple times.
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
      modifiedExpr = expr,
      masses,massDefine,codeMass,massRules,
      couplings,couplingDefine,codeCoupling,couplingRules,
      loopRules,loopArrayDefine,loopArraySet,
      cxxExpr,updatingVars},
   masses = Tally@Cases[modifiedExpr, _SARAH`Mass, Infinity, Heads->True];
   {massDefine, codeMass, massRules} = `cxx`nameMasses@masses;
   modifiedExpr = modifiedExpr /. massRules;

   couplings = Tally@Cases[modifiedExpr,SARAH`Cp[__][___],Infinity,Heads->True];
   {couplingDefine, codeCoupling, couplingRules} = `cxx`nameCouplings@couplings;
   modifiedExpr = modifiedExpr /. couplingRules;

   {loopRules,loopArrayDefine,loopArraySet,modifiedExpr} =
      createLoopFunctions@modifiedExpr;

   cxxExpr = `cxx`applyRules/@modifiedExpr;
   updatingVars = MapThread[#1<>" += "<>#2<>";"&, {First/@$basis, cxxExpr}];

   TextFormatting`ReplaceCXXTokens[code,{
      "@fieldAliases@"->`cxx`setFieldAliases@couplings,
      "@defineMasses@"->massDefine,
      "@defineCouplings@"->couplingDefine,
      "@defineLoopFunctions@"->loopArrayDefine,
      "@BeginSum@"->`cxx`beginSum@summation,
      "@setMasses@"->codeMass,
      "@setCouplings@"->codeCoupling,
      "@skipZeroAmplitude@"->`cxx`skipZeroAmplitude[
         modifiedExpr,loopRules,massRules],
      "@setLoopFunctions@"->loopArraySet,
      "@ChangeOutputValues@"->Utils`StringJoinWithSeparator[updatingVars,"\n"],
      "@EndSum@"->`cxx`endSum@GetGenericFields@summation
   }]];
`cxx`changeGenericExpressions // secure;

Module[{func},

func = If[#1 =!= #2, "using " <> #1 <> " = " <> #2 <> ";\n", ""]&;

`cxx`setFieldAliases::usage = "
@brief C++ names for fields with bar and conj are cumbersome. In order to
       improve readaability and simplify checks some aliases are generated.
@param couplings A list of tuples with couplings.
@returns A string with C++ code of alias definitions.";
`cxx`setFieldAliases[couplings:{{SARAH`Cp[__][___], _Integer}..}] :=
Module[{l = Sort@DeleteDuplicates[(Sequence@@Head@First@#&)/@couplings]},
   StringJoin@MapThread[func, {`cxx`fieldAlias/@l, `cxx`fieldName/@l}]];
`cxx`setFieldAliases // secure;

];

Module[{c, name, strip},

c[f_] := Switch[Head@f, SARAH`bar|Susyno`LieGroups`conj, "_"<>#, _, #]&;
name[f_] := ToString@Switch[f, _Symbol, f, _, Head@f];
strip[f_] := f /. {SARAH`bar -> Identity, Susyno`LieGroups`conj -> Identity};

`cxx`fieldAlias::usage = "
@brief Generates a short C++ name for a field, whether conjugated or not.
@param f A external or generic field.
@returns A C++ name for a field.";
`cxx`fieldAlias[f:_?IsParticle] := c[f][name@strip@f];
`cxx`fieldAlias[f_?IsGenericParticle] := c[f][`cxx`fieldName@strip@f];
`cxx`fieldAlias // secure;];

createLoopFunctions[modifiedExpr:{__}] :=
Module[{onePoint,
      onePointTemplate =
         {  LoopTools`A0@@#2 -> "a"<>#1<>"[0]",
            LoopTools`A0i[LoopTools`aa0,Sequence@@#2] -> "a"<>#1<>"[0]"}&,
      twoPoint,
      twoPointTemplate =
         {  LoopTools`B0i[LoopTools`bb0,Sequence@@#2] -> "b"<>#1<>"[0]",
            LoopTools`B0i[LoopTools`bb1,Sequence@@#2] -> "b"<>#1<>"[1]",
            LoopTools`B0i[LoopTools`bb00,Sequence@@#2] -> "b"<>#1<>"[2]"}&,
      threePoint,
      threePointTemplate =
         {  LoopTools`C0i[LoopTools`cc0,Sequence@@#2] -> "c"<>#1<>"[0]",
            LoopTools`C0i[LoopTools`cc1,Sequence@@#2] -> "c"<>#1<>"[1]",
            LoopTools`C0i[LoopTools`cc2,Sequence@@#2] -> "c"<>#1<>"[2]",
            LoopTools`C0i[LoopTools`cc00,Sequence@@#2] -> "c"<>#1<>"[3]",
            LoopTools`C0i[LoopTools`cc11,Sequence@@#2] -> "c"<>#1<>"[4]",
            LoopTools`C0i[LoopTools`cc12,Sequence@@#2] -> "c"<>#1<>"[5]",
            LoopTools`C0i[LoopTools`cc22,Sequence@@#2] -> "c"<>#1<>"[6]"}&,
      fourPoint,
      fourPointTemplate =
         {  LoopTools`D0i[LoopTools`dd0,Sequence@@#2] -> "d"<>#1<>"[0]",
            LoopTools`D0i[LoopTools`dd1,Sequence@@#2] -> "d"<>#1<>"[1]",
            LoopTools`D0i[LoopTools`dd2,Sequence@@#2] -> "d"<>#1<>"[2]",
            LoopTools`D0i[LoopTools`dd3,Sequence@@#2] -> "d"<>#1<>"[3]",
            LoopTools`D0i[LoopTools`dd00,Sequence@@#2] -> "d"<>#1<>"[4]",
            LoopTools`D0i[LoopTools`dd11,Sequence@@#2] -> "d"<>#1<>"[5]",
            LoopTools`D0i[LoopTools`dd12,Sequence@@#2] -> "d"<>#1<>"[6]",
            LoopTools`D0i[LoopTools`dd13,Sequence@@#2] -> "d"<>#1<>"[7]",
            LoopTools`D0i[LoopTools`dd22,Sequence@@#2] -> "d"<>#1<>"[8]",
            LoopTools`D0i[LoopTools`dd23,Sequence@@#2] -> "d"<>#1<>"[9]",
            LoopTools`D0i[LoopTools`dd33,Sequence@@#2] -> "d"<>#1<>"[10]"}&,
      append,
      loopRules = {},
      loopArrayDefine = {},
      loopArraySet = {}},
   append[loopFunctions_List,function_,functionName_String,arrayName_String] :=
   If[loopFunctions=!={},
      AppendTo[loopArrayDefine,Array[
         "looplibrary::"<>functionName<>"coeff_t "<>arrayName<>ToString@#<>"{}"&,
         Length@loopFunctions]];
      AppendTo[loopArraySet,Array[
         Parameters`ExpressionToString[
            StringJoin["Loop_library::get().",functionName][
               arrayName<>ToString@#,
               Sequence@@loopFunctions[[#,1]],
               "Sqr(context.scale())"]]<>"; // "<>
            ToString@loopFunctions[[#,2]]<>" copies."&,
         Length@loopFunctions]];
      AppendTo[loopRules,
         Join@@function@@@Array[{ToString@#,loopFunctions[[#,1]]}&,
         Length@loopFunctions]];];
   onePoint = Tally@Join[
      Cases[modifiedExpr,LoopTools`A0i[_,args:__]:>{args},Infinity],
      Cases[modifiedExpr,LoopTools`A0[args:__]:>{args},Infinity]];
   twoPoint = Tally@Join[
      Cases[modifiedExpr,LoopTools`B0i[_,args:__]:>{args},Infinity],
      Cases[modifiedExpr,(LoopTools`B0|LoopTools`B0)[args:__]:>{args},
         Infinity]];
   threePoint = Tally@Join[
      Cases[modifiedExpr,LoopTools`C0i[_,args:__]:>{args},Infinity],
      Cases[modifiedExpr,(LoopTools`C0)[args:__]:>{args},Infinity]];
   fourPoint = Tally@Join[
      Cases[modifiedExpr,LoopTools`D0i[_,args:__]:>{args},Infinity],
      Cases[modifiedExpr,(LoopTools`D0)[args:__]:>{args},Infinity]];
   append[onePoint,onePointTemplate,"A","a"];
   append[twoPoint,twoPointTemplate,"B","b"];
   append[threePoint,threePointTemplate,"C","c"];
   append[fourPoint,fourPointTemplate,"D","d"];
   {  Flatten@loopRules,
      Utils`StringJoinWithSeparator[Join@@loopArrayDefine,";\n"]<>";\n",
      Utils`StringJoinWithReplacement[Join@@loopArraySet,"\n","\""->""],
      modifiedExpr/.Flatten@loopRules}];
createLoopFunctions // secure;

`cxx`skipZeroAmplitude::usage = "
@brief If some combination of couplings in the amplitude is zero, then the full
       amplitude is zero as well, so that we can skip it.
@note How this is implemented? Any amplitude has the following form: several
      couplings times some masses times loop integral. We can get all different
      combinations of this coupling coefficients and if all of them are zero,
      then amplitude is zero as well.";
`cxx`skipZeroAmplitude[
   modifiedExpr:{__},
   loopRules:{Rule[_,_]...},
   massRules:{Rule[_,_]..}] :=
Module[{massesToOne = Rule[#,1] & /@ massRules[[All,2]],
      loopsToOne = {},
      result, func = "z[" <> # <> "]"&},
   If[0 < Length@loopRules,
      loopsToOne = Rule[#,1] & /@ loopRules[[All,2]]
   ];
   result = removeNumbers[ExpandAll[modifiedExpr]] /.Plus->List/.
      massesToOne/.loopsToOne;
   result = DeleteDuplicates[Flatten@DeleteCases[result,1]];
   If[ 2 === LeafCount@result,
      result = func@@result,
      result = Plus@@(result/.HoldPattern[Times[x__]]:>Times@@(func@#&/@{x}))];
   result = "if( "<>
      StringReplace[Parameters`ExpressionToString@result,"\""->""]<>
      " == 0 ) continue;";
   StringReplace[result,
      RegularExpression["g(\\d+)"]:> ToString[ToExpression@"$1"-1]]];
`cxx`skipZeroAmplitude // secure;

`cxx`getVariableName[SARAH`Mass[obj_?IsGenericParticle]] :=
Switch[GetParticleName@obj,
   GenericS, "mS",
   GenericF, "mF",
   GenericV, "mV",
   GenericU, "mU"]<>ToString@GetIndex@obj;
`cxx`getVariableName[SARAH`Mass[obj:_?IsParticle]] :=
   "m"<>ToString@GetParticleName@obj<>GetCXXIndex@obj;
`cxx`getVariableName // secure;

Module[{
      info = {`cxx`getVariableName@#1, `cxx`applyRules@#1, #1, #2}&@@#&,
      d = "double " <> Utils`StringJoinWithSeparator[#, ", "] <> ";"&,
      i = #1 <> " = " <> #2 <> "; // " <> ToString@#4 <>" copies.\n"&@@#&,
      r = Rule@@@#[[All, {3, 1}]]&},

`cxx`nameMasses::usage = "
@brief Generates names for masses to be used inside generic sums and then
       creates:

       1. a C++ code for definition and initialisation for masses,
       2. a ``Mathematica`` to C++ rules for generated names of masses.
@param masses A list of tallies with Mathematica expressions for mass and the
       number of repetition of it.
@returns A list of:

         1. a C++ string with definitions,
         2. a C++ string with initialisations,
         3. a Mathematica list of rules for mass convertion
            to C++ code.";
`cxx`nameMasses[masses:{{_, _Integer}..}] :=
   {d[First/@#], StringJoin[i/@#], r@#} &@ Table[info@m, {m, Sort@masses}];
`cxx`nameMasses // secure;

`cxx`nameCouplings::usage = "
@brief Generates names for couplings to be used inside generic sums and then
       creates:

       1. C++ code for definition and initialisation for couplings,
       2. ``Mathematica`` to C++ rules for generated names of couplings.
@param masses A list of tallies with Mathematica expressions for coupling
       and the number of repetition of it.
@returns A list of:

         1. a C++ string with definitions,
         2. a C++ string with initialisations,
         3. a ``Mathematica`` list of rules for coupling convertion
            to C++ code.
@note All couplings have to be multiplied by ``I``, as it is done in these
      rule replacements.";
`cxx`nameCouplings[couplings:{{_,_Integer}..}] :=
Module[{
      sort = Sort@couplings, info, isZero,
      d = "std::complex<double> " <> Utils`StringJoinWithSeparator[#, ", "] <>
         ";"&,
      timesI = Rule[#1, I*#2]&@@#&,
      setZero},
   isZero = "std::array<int, "<> ToString@Length@sort <> "> z{};";
   setZero = "z[" <> ToString[#-1] <> "] = (std::abs(g" <> ToString@# <>
      ") < std::numeric_limits<double>::epsilon()) ? 0: 1;\n"&;
   info = {"g"<>ToString@#, `cxx`applyRules@sort[[#,1]], sort[[#,1]],
      sort[[#,2]]}&;
   {  d[First/@#]<>"\n"<>isZero,
      StringJoin[i/@#]<>"\n"<>StringJoin[setZero/@Range@Length@sort],
      timesI/@r@#} &@ Array[info, Length@sort]];
`cxx`nameCouplings // secure;];

`cxx`shortNames::usage = "
@brief Generates C++ code for type abbreviations stored in
       ``GenericFieldMap`` (Associative Sequence) at ``Key`` positions.
@param genFields List of generic fields.
@returns String C++ code for type abbreviations stored in GenericFieldMap
         (Associative Sequence) at Key positions.";
`cxx`shortNames[genFields:{__?IsGenericParticle}] :=
   Utils`StringJoinWithSeparator[Apply[
      "using "<>#1<>" = typename at<GenericFieldMap,"<>ToString@#2<>">::type;"&,
      {`cxx`fieldName@#,`cxx`genericFieldKey@#}&/@genFields,
      {1}],"\n"];
`cxx`shortNames // secure;

`cxx`beginSum::usage = "
@brief Generates C++ code for sum beginning used inside GenericSum.
@param summation List of generic index restriction rules pares, which,
       if are true should lead to a skip of summation.
@returns String C++ code for sum beginning used inside generic sums.";
`cxx`beginSum[summation_?IsSummation]:=
Module[{beginsOfFor},
   beginsOfFor = "for( const auto &"<>GetCXXIndex[#[[1]]]<>" : "<>
      "index_range<"<> `cxx`fieldName[#[[1]]]<>">() ) {\n"<>
      "at_key<"<>`cxx`genericFieldKey@#[[1]]<>">( index_map ) = "<>
      GetCXXIndex[#[[1]]]<>";"<>parseRestrictionRule[#] &/@summation;
   Utils`StringJoinWithSeparator[beginsOfFor,"\n"]];
`cxx`beginSum // secure;

parseRestrictionRule[{genericField_?IsGenericParticle,rule_}] :=
Module[{f1,f2,GetIndexOfExternalField,OrTwoDifferent},
   GetIndexOfExternalField[_[_[{ind_}]]] := `cxx`applyRules@ind;
   GetIndexOfExternalField[_[{ind_}]] := `cxx`applyRules@ind;
   GetIndexOfExternalField[_] := "i0";
   OrTwoDifferent[] :=
   Module[{type1 = `cxx`fieldName@First@rule,
         type2 = `cxx`fieldName@Last@rule,
         ind = GetIndexOfExternalField@First@rule,
         typeGen = `cxx`fieldName@genericField,
         indGen = GetCXXIndex@genericField},
      "\nif( (boost::core::is_same<"<>typeGen<>","<>type1<>
         ">::value || boost::core::is_same<"<>typeGen<>","<>type2<>
         ">::value) && "<>indGen<>" == "<>ind<>" ) continue;"];
   Switch[rule,
      Or[f1_,f2_],
         OrTwoDifferent[],
      False,
         "",
      _,
         "@todo This rule is not implemented yet!";Quit[1]]];
parseRestrictionRule // secure;

`cxx`endSum::usage = "
@brief Generates C++ code for end of sum over generic fields inside
       GenericSum.
@param genFields List of generic fields.
@returns String C++ code for end of sum over generic fields inside
         GenericSum.";
`cxx`endSum[genFields:{__?IsGenericParticle}] :=
   StringJoin[
      Array["}"&,Length@genFields],
      " // End of summation over generic fields"];
`cxx`endSum // secure;

`cxx`calculateFunction::usage = "
@brief Generates C++ code for functions which return result of generic sum
       calculation.
@param genSumNames list of strings with names of generic sums.
@returns Generates C++ code for functions which return result of
         generic sum calculation.";
`cxx`calculateFunction[genSumNames:{__String}] :=
Module[{
      varName = "genericsum" (* Feel free to change me to another C++ name *),
      varNames,initVars,sumOfSums},
   varNames = Array[varName<>ToString@#&,Length@genSumNames];
   initVars = Utils`StringJoinWithSeparator[
      MapThread["const auto "<>#1<>" = "<>#2<>"();"&,{varNames,genSumNames}],
      "\n"];
   sumOfSums = Utils`StringJoinWithSeparator[#<>"[i]"&/@varNames,"+"];
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
      {  "@BasisLength@"->GetLength@$basis,
         "@InitializeVariablesWhichStoreGenericSumsOutput@"->initVars,
         "@SumOfVariables@"->sumOfSums}]];
`cxx`calculateFunction // secure;

`cxx`insertFields::usage = "
@brief Generates C++ code for class insertions inside GenericSum.
@param genInsertions list of list with SARAH particle names.
@returns String C++ code for class insertions inside GenericSum.";
`cxx`insertFields[genInsertions_?IsClassFields] :=
   Utils`StringJoinWithSeparator["boost::mpl::vector<"<>
      Utils`StringJoinWithSeparator[`cxx`fieldName@#&/@#,", "]<>
      ">"&/@genInsertions,",\n"];
`cxx`insertFields // secure;

`cxx`insertFactors::usage = "
@brief Generates C++ code for combinatorical factor insertions inside GenericSum.
@param combinatorialFactors List of integers.
@returns String C++ code for combinatorical factor insertions inside GenericSum.";
`cxx`insertFactors[combinatorialFactors_?IsCombinatoricalFactors] :=
   Utils`StringJoinWithSeparator["boost::mpl::int_<"<>ToString@#<>
      ">"&/@combinatorialFactors,",\n"];
`cxx`insertFactors // secure;

`cxx`insertColours::usage = "
@brief Generates C++ code for colour factor insertions inside GenericSum.
@param colourFactors list of numbers.
@returns String C++ code for colour factor insertions inside GenericSum.";
`cxx`insertColours[colourFactors:{__?NumberQ}] :=
Module[{
      ReRatioColourFactors = {Numerator@#,Denominator@#} &/@ Re@colourFactors,
      ImRatioColourFactors = {Numerator@#,Denominator@#} &/@ Im@colourFactors},
   Utils`StringJoinWithSeparator[
      StringReplace[
         MapThread[
            "detail::complex_helper<"<>
            "detail::ratio_helper<"<>ToString@#1<>">,"<>
            "detail::ratio_helper<"<>ToString@#2 <> ">>"&,
            {ReRatioColourFactors, ImRatioColourFactors}],
         {"{" -> "", "}" -> ""}],
   ",\n"]];
`cxx`insertColours // secure;

End[];
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "NPointFunctions`"];
