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
      GenericS, GenericF, GenericV, GenericU, Present, Absent
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
   "ModifyMasses.m"
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
   tree = GenerateDiagrams[FAIncomingFields, FAOutgoingFields];
   tree = ApplyObservableSetting[tree, diagrams];
   tree = GenerateColorlessAmplitudes[tree];
   tree = ApplyObservableSetting[tree, amplitudes];
   ExportFeynArtsPaint@tree;
   {FieldRules@GetFields@tree, CalculateAmplitudes@tree}
];
NPointFunction // secure;

CalculateAmplitudes[tree_?IsTree] :=
Module[{ampsGen, feynAmps, generic, chains, subs, zeroedRules},
   ampsGen = FeynArts`PickLevel[Generic][GetAmplitudes@tree];
   If[$zeroExternalMomenta,
      ampsGen = FormCalc`OffShell[ampsGen, Sequence@@Array[#->0&, Tr@$externalFieldNumbers]]
   ];
   feynAmps = MapThreadWithBar[
      FormCalc`CalcFeynAmp[
         Head[ampsGen][#1],
         FormCalc`Dimension -> #2,
         FormCalc`OnShell -> $onShell,
         FormCalc`FermionChains -> FormCalc`Chiral,
         FormCalc`FermionOrder -> GetObservableSetting@order,
         FormCalc`Invariants -> False,
         FormCalc`MomElim -> #3
      ]&,
      {
         ampsGen,
         ApplyObservableSetting[tree, regularization],
         ApplyObservableSetting[tree, momenta]
      },
      "Amplitude calculation"
   ] //. FormCalc`GenericList[];
   generic = MapThread[getGenericSum, {feynAmps, ApplyObservableSetting[tree, sum]}];
   {generic, chains, subs} = ProceedChains[tree, generic];

   MassRules[tree, feynAmps];
   {generic, chains, subs} = ModifyMasses[{generic, chains, subs}, tree, $zeroExternalMomenta];

   convertToFS[
      {
         generic,
         GetFieldInsertions@tree,
         combinatoricalFactors@tree,
         colorFactors@tree
      },
      chains,
      subs
   ] /. `rules`externalMomenta[tree, $zeroExternalMomenta]
];
CalculateAmplitudes // secure;

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

GetFieldInsertions[tree_?IsTree] := GetFieldInsertions[diagrams@tree, False];
GetFieldInsertions[diagrams_, withRules_:False] :=
Module[{genericFields = {}, classesFields = {}, genericRules, finalRule, removeIndices, res},
   removeIndices = f_[n_, {_}] :> f@n;
   Cases[
      diagrams,
      Rule[_[_, Generic == _]@g__, _@c__] :> (
         AppendTo[genericFields, {g}];
         AppendTo[classesFields, {c} /. _[_, FeynArts`Classes == _] -> List]
      ),
      Infinity
   ];
   genericRules = Switch[withRules,
      True,  Rule[f_[n_], t_Symbol] :> (n -> f[n]),
      False, Rule[f_[n_], t_Symbol] :> (t[FeynArts`Index[Generic, n]] -> f[n])
   ];
   finalRule = Switch[withRules,
      True,  Rule[n_, rhs_] :> Rule[FeynArts`Field@n, rhs],
      False, Rule[_, rhs_] :> rhs
   ];
   Table[
      res = Cases[genericFields[[i]], genericRules] /. classesFields[[i]];
      res = SortBy[#, First] &/@ res;
      res = res /. removeIndices /. finalRule,
      {i, Length@genericFields}
   ]
];
GetFieldInsertions // secure;

removeParticleIndices[Times[-1, field_]] := -removeParticleIndices@field;
removeParticleIndices[name_[class_, ___]] := name@class;
removeParticleIndices // secure;

MapThreadWithBar[func_, exprs:{__}, text_String] :=
Module[{printed = 0, delta, out, tot, print, def = 70},
   tot = Length@First@exprs;
   print[i_] :=
   (
      delta = Floor[(def-StringLength[text]-4)*i/tot] - printed;
      subWrite[StringJoin@@Array["."&, delta]];
      printed += delta;
   );
   subWrite[text<>": ["];
   out = Table[print@i; func@@exprs[[All, i]], {i, tot}];
   subWrite@"]\n";
   out
];
MapThreadWithBar // secure;

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
