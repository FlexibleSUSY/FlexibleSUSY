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

Needs["Utils`", FileNameJoin@{ParentDirectory@DirectoryName@$InputFileName, "Utils.m"}];

Block[{Format},
   Needs@"FeynArts`";
   Needs@"FormCalc`";
   Print[];
];

BeginPackage@"NPointFunctions`";

(* Start reserving names *)

Off[General::shdw];
   DiracChain;
   Mat;
On[General::shdw];

NPointFunction;
SetAttributes[
   {
      LorentzIndex, GenericSum, GenericIndex, OperatorsOnly, ExceptLoops,
      GenericS, GenericF, GenericV, GenericU
   },
   Protected
];

(* End reserving names *)

Begin@"`Private`";

secure[sym_Symbol] :=
   Protect@Evaluate@Utils`MakeUnknownInputDefinition@sym;
secure // secure;

subWrite // Protect;

Utils`DynamicInclude/@{
   "PatternChecks.m",
   "rules.m",
   "settings.m",
   "chains.m",
   "topologies.m",
   "tree.m",
   "mass.m"
};

NPointFunction[
   {
      FCOutputDir:_String,
      FAModelName:_String,
      particleNamesFile_?FileExistsQ,
      particleNamespaceFile_?FileExistsQ,
      FAIncomingFields:{__String},
      FAOutgoingFields:{__String}
   },
   {
      observable:None | _?(Context@Evaluate@Head[#] == "FlexibleSUSYObservable`"&),
      loops_,
      processes:{___String},
      momenta:_Symbol,
      onShell:_Symbol,
      mainRegularization:_Symbol
   }
] :=
Module[{tree},
   BeginPackage["NPointFunction`"];
   Begin["`Private`"];
      $particleNamesFile = particleNamesFile;
      $particleNamespaceFile = particleNamespaceFile;

      $observableName = SymbolName@Head@observable;
      $externalFieldNumbers = {Length@FAIncomingFields, Length@FAOutgoingFields};
      $loopNumber = loops;
      $expressionsToDerive = processes;
      $zeroExternalMomenta = momenta;
      $onShell = onShell;
      $regularizationScheme =  Switch[mainRegularization, FlexibleSUSY`DRbar, 4, FlexibleSUSY`MSbar, D];
   End[];
   EndPackage[];

   FeynArts`$FAVerbose = 0;
   FeynArts`InitializeModel@FAModelName;
   SetOptions[FeynArts`InsertFields,
      FeynArts`Model -> FAModelName,
      FeynArts`InsertionLevel -> FeynArts`Classes
   ];
   FormCalc`$FCVerbose = 0;
   If[!DirectoryQ@FCOutputDir, CreateDirectory@FCOutputDir];
   SetDirectory@FCOutputDir;

   DefineAllowedTopologies[];
   LoadAllSettings[];
   tree = settings[TreeFromDiagrams[FAIncomingFields, FAOutgoingFields], diagrams];
   tree = settings[plant@tree, amplitudes];
   picture@tree;
   {`rules`fields@fields@tree, calculateAmplitudes@tree}
];
NPointFunction // secure;

genericIndex[index:_Integer] := FeynArts`Index[Generic, index];
genericIndex // secure;

process[set:_?IsDiagramSet|_?IsAmplitudeSet] :=
   FirstCase[Head@set, (FeynArts`Process -> e_) :> e];

process[set_?IsFormCalcSet] :=
   Part[Head@Part[set, 1], 1];

process // secure;

getField[set:_?IsDiagramSet, i:_Integer] :=
   Flatten[List@@process@set, 1][[i]] /; 0<i<=Plus@@(Length/@process@set);
getField // secure;

fieldInsertions::usage = "
@brief Finds insertions, related to fields.
@param tree A ``tree`` object.
@param diag A single diagram.
@param graph A ``FeynmanGraph[__][__]`` object.
@param insert A ``Insertions[Classes][__]`` object.
@param keepNumQ Responsible for the type of output field names.
       ``FeynmanGraph`` on a generic level contains
       ``Field[num] -> <generic particle>``.
       ``FeynmanGraph`` on a classes level contains
       ``Field[num] -> <classes particle>``.

       * ``True``: then ``Field[_] -> <classes particle>`` is created
         for a diagram.
       * ``False``: ``<generic particle> -> <classes particle>`` is created
         for a diagram.
@returns * For a single diagram returns 1) ``List`` for topology level of
           2) ``List`` for generic level of 3) ``List`` for classes level of
           field insertion rules::

              1) 2) 3)
              {  {  {Rule[<expr>, <class field>]..}..}..}

         * For a set of diagrams only <class field> is taken instead of the
           whole ``Rule``.
@note All indices in rhs. of rules are removed.";
fieldInsertions[tree:_?IsTree] :=
   Map[Last, #, {3}] &@ Flatten[ fieldInsertions /@ List@@diagrams@tree, 1];
fieldInsertions[diag_?IsDiagram, keepNumQ:True|False:False] :=
   fieldInsertions[#, keepNumQ] &/@ Apply[List, Last@diag, {0, 1}];
fieldInsertions[{graph_, insert_}, keepNumQ_] :=
Module[{toGenericIndexConventionRules, fieldsGen, genericInsertions},
   toGenericIndexConventionRules = Cases[graph,
      Rule[FeynArts`Field[index_Integer],type_Symbol] :>
      Rule[FeynArts`Field@index, type[FeynArts`Index[Generic,index]]]];
   fieldsGen = toGenericIndexConventionRules[[All,1]];
   genericInsertions = Cases[#,
      Rule[genericField_,classesField_] /; MemberQ[fieldsGen, genericField] :>
      Rule[genericField, removeParticleIndices@classesField]] &/@ insert;
   SortBy[#,First]&/@ If[keepNumQ,
      List @@ genericInsertions,
      List @@ genericInsertions /. toGenericIndexConventionRules]];
fieldInsertions // secure;

removeParticleIndices[Times[-1, field_]] := -removeParticleIndices@field;
removeParticleIndices[name_[class_, ___]] := name@class;
removeParticleIndices // secure;

calculateAmplitudes::usage = "
@brief Applies ``FormCalc`` routines to amplitude set, simplifies the result.
@param tree A set of data in the form of a ``tree`` object.
@returns The main part of n-point function object, containing:

         * generic amplitudes,
         * class specific insertions,
         * subexpressions.";
calculateAmplitudes[tree:_?IsTree] :=
Module[{proc, ampsGen, feynAmps, generic, chains, subs, zeroedRules},
   proc = process@amplitudes@tree;
   ampsGen = FeynArts`PickLevel[Generic][amplitudes@tree];
   If[$zeroExternalMomenta,
      ampsGen = FormCalc`OffShell[ampsGen,
         Sequence@@Array[#->0&, Plus@@Length/@proc]
      ]
   ];
   feynAmps = mapThread[
      FormCalc`CalcFeynAmp[Head[ampsGen][#1],
         FormCalc`Dimension -> #2,
         FormCalc`OnShell -> $onShell,
         FormCalc`FermionChains -> FormCalc`Chiral,
         FormCalc`FermionOrder -> settings@order,
         FormCalc`Invariants -> False,
         FormCalc`MomElim -> #3]&,
      {ampsGen, settings[tree, regularization], settings[tree, momenta]},
      "Amplitude calculation"
   ] //. FormCalc`GenericList[];

   generic = MapThread[getGenericSum, {feynAmps, settings[tree, sum]}];

   {generic, chains, subs} = proceedChains[tree, generic];

   mass`rules[tree, feynAmps];
   {generic, chains, subs} = mass`modify[{generic, chains, subs},
      tree,
      $zeroExternalMomenta
   ];

   convertToFS[
      {  generic,
         fieldInsertions@tree,
         combinatoricalFactors@tree,
         colorFactors@tree},
      chains,
      subs] /. `rules`externalMomenta[tree, $zeroExternalMomenta]];
calculatedAmplitudes // secure;

mapThread::usage = "
@brief Behaves like ``MapThread``, but also prints a progress bar.
@param func A function to apply to set of data.
@param exprs A ``List`` of listable sets with data.
@param text A string to be printed.
@todo Add check for equality of length for exprs.";
mapThread[func_, exprs:{__}, text_String] :=
   Module[{printed = 0, delta, out, tot, print, def = 70},
      tot = Length@First@exprs;
      print[i_] :=
      (  delta = Floor[(def-StringLength[text]-4)*i/tot] - printed;
         subWrite[StringJoin@@Array["."&, delta]];
         printed += delta;);
      subWrite[text<>": ["];
      out = Table[print@i; func@@exprs[[All, i]], {i, tot}];
      subWrite@"]\n";
      out];
mapThread // secure;

getGenericFields::usage = "
@brief Generates a list of unique sorted generic fields in expression.
@param expr An expression, where to search.
@returns A list of unique sorted generic fields.";
getGenericFields[expr:_] :=
   Sort@DeleteDuplicates[Cases[expr, type`genericField, Infinity]];
getGenericFields // secure;

getGenericSum::usage= "
@brief Converts ``FormCalc`Amp`` into ``NPointFunctions`GenericSum`` object
       using restriction rules for generic fields.
@param amplitude ``FormCalc`Amp`` expression.
@param sumRules A set of rules, restricting the summation.
@returns A ``NPointFunctions`GenericSum`` object.";
getGenericSum[amplitude:type`fc`amplitude, sumRules:{Rule[_Integer, _]...}] :=
Module[{sort, rules},
   sort = getGenericFields@amplitude;
   rules = Append[sumRules, _Integer -> False];
   GenericSum[
      List@@amplitude,
      sort /. f_[_[_,i_]] :> {f@GenericIndex@i, i /. rules}]];
getGenericSum // secure;

convertToFS::usage = "
@brief Translate a list of ``FormCalc`` amplitudes, abbreviations and
       subexpressions into ``FlexibleSUSY`` language.
@param amplitudes A ``List`` of amplitudes.
@param abbreviations A ``List`` of abbreviations.
@param subexpressions A ``List`` of subexpressions.
@returns A list of amplitudes and joined abbreviations and subexpressions.";
convertToFS[amplitudes_, abbreviations_, subexpressions_] :=
   {  `rules`amplitude@amplitudes,
      `rules`subexpressions/@Join[abbreviations, subexpressions]};
convertToFS // secure;

End[];
Block[{$ContextPath}, EndPackage[]];
