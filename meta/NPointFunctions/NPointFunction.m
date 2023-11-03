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

Begin@"`Private`";

secure[sym_Symbol] :=
   Protect@Evaluate@Utils`MakeUnknownInputDefinition@sym;
secure // secure;

subWrite // Protect;

Utils`DynamicInclude/@{
   "PatternChecks.m",
   "Rules.m",
   "Settings.m",
   "chains.m",
   "Topologies.m",
   "TreeWrapper.m",
   "ModifyMasses.m"
};

NPointFunction[
   {
      FCOutputDir_String,
      FAModelName_String,
      particleNamesFile_?FileExistsQ,
      particleNamespaceFile_?FileExistsQ,
      FAIncomingFields:{__String},
      FAOutgoingFields:{__String}
   },
   {
      observable:None | _?(Context@Evaluate@Head@# == "FlexibleSUSYObservable`"&),
      loops_,
      processes:{___String},
      momenta_Symbol,
      onShell_Symbol,
      mainRegularization_Symbol
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
   ampsGen = FeynArts`PickLevel[Generic][ExtractAmplitudes@tree];
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
   generic = MapThread[GetGenericSum, {feynAmps, ApplyObservableSetting[tree, sum]}];
   {generic, chains, subs} = ProceedChains[tree, generic];

   MassRules[tree, feynAmps];
   {generic, chains, subs} = ModifyMasses[{generic, chains, subs}, tree, $zeroExternalMomenta];

   ConvertToFS[
      {
         generic,
         GetFieldInsertions@tree,
         CombinatoricalFactors@tree,
         ColorFactors@tree
      },
      chains,
      subs
   ] /. ExternalMomentaRules[tree, $zeroExternalMomenta]
];
CalculateAmplitudes // secure;

GetProcess[set:_?IsDiagramSet|_?IsAmplitudeSet] :=
   FirstCase[Head@set, (FeynArts`Process -> e_) :> e];
GetProcess[set_?IsFormCalcSet] :=
   Part[Head@Part[set, 1], 1];
GetProcess // secure;

GetFieldInsertions[tree_?IsTree] := GetFieldInsertions[ExtractDiagrams@tree, False];
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

GetGenericFields[expr_] := Sort@DeleteDuplicates[Cases[expr, _?IsGenericField, Infinity]];
GetGenericFields // secure;

GetGenericSum[amplitude_?IsFormCalcAmplitude, sumRules:{Rule[_Integer, _]...}] :=
Module[{sort, rules},
   sort = GetGenericFields@amplitude;
   rules = Append[sumRules, _Integer -> False];
   GenericSum[
      List@@amplitude,
      sort /. f_[_[_,i_]] :> {f@GenericIndex@i, i /. rules}]];
GetGenericSum // secure;

ConvertToFS[amplitudes_, abbreviations_, subexpressions_] :=
{
   AmplitudeRules@amplitudes,
   SubexpressionRules /@ Join[abbreviations, subexpressions]
};
ConvertToFS // secure;

End[];
Block[{$ContextPath}, EndPackage[]];
