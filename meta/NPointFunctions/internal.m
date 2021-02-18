(* ::Package:: *)

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

(* There is a problem with Global`args which comes from mathematica paclets.*)
Quiet[Needs@"FeynArts`", {FeynArts`args::shdw}];
FeynArts`$FAVerbose = 0;
FeynArts`InitializeModel@NPointFunctions`$FAMod;
SetOptions[FeynArts`InsertFields,
   FeynArts`Model -> NPointFunctions`$FAMod,
   FeynArts`InsertionLevel -> FeynArts`Classes];

Needs@"FormCalc`";Print[];
FormCalc`$FCVerbose = 0;
(If[!DirectoryQ@#, CreateDirectory@#]; SetDirectory@#;)&@NPointFunctions`$FCDir;

(* Next Format makes some pattern generate mistakes. *)
Format[FormCalc`DiracChain[FormCalc`Private`s1_FormCalc`Spinor,FormCalc`Private`om_,FormCalc`Private`g___,FormCalc`Private`s2_FormCalc`Spinor]] =.;
Needs@"Utils`";

BeginPackage@"NPointFunctions`";

{  SetInitialValues, NPointFunctionFAFC};
{  LorentzIndex, GenericSum, GenericIndex,
   GenericS, GenericF, GenericV, GenericU,
   DimensionalReduction, DimensionalRegularization,
   OperatorsOnly, ExceptLoops} ~ SetAttributes ~ {Protected};

Begin@"`internal`";

secure[sym:_Symbol] :=
Protect@Evaluate@Utils`MakeUnknownInputDefinition@sym;
secure // secure;

Module[{impl},

impl[s:_Symbol, RuleDelayed[{p:___}, d:_]] :=
   SetDelayed[s[p], d];
impl ~ SetAttributes ~ {HoldAllComplete};

define::usage = "
@brief Defines a set of function with a given name in a safe way.
@param s A symbol, which represent a function name.
@param e A sequence of delayed rules. On lhs there is a list with pattern for a
       new function, on rhs there is function body.
@param a A set of attributes to be applied to a new function. Default is a
       protected-locked combination.";
define[s:_Symbol, e:RuleDelayed[{___},_].., a:{__Symbol}:{Protected, Locked}
] := (
   impl[s, ##] &@@@ Hold /@ {e};
   s // Utils`MakeUnknownInputDefinition;
   s ~ SetAttributes ~ a;
);
define // Utils`MakeUnknownInputDefinition;
define ~ SetAttributes ~ {HoldAllComplete, Protected, Locked};

];

createGetSetOnce::usage = "
@brief Defines safe getter and setter for a given symbol by prepending \"set\"
       and \"get\" in front of it. This variable can be set only once.
@param sym A symbol, which serves as a root name for a new functions, i.e. for
       MySymbol symbol setMySymbol and getMySymbol will be created.
@param pattern An expression, which will be used for a pattern for a set
       function, i.e. it will have the form setMySymbol[new:pattern].";
createGetSetOnce[{sym:_Symbol, pattern:_}] :=
Module[{set, get, once, value},
   set = Symbol["set"<>SymbolName@sym];
   get = Symbol["get"<>SymbolName@sym];
   set::errOnce = "The value can be set only once.";
   set[new:pattern] := (
      Utils`AssertOrQuit[!TrueQ@once, set::errOnce];
      value = new;
      once = True;);
   secure@Evaluate@set;

   get::errNotSet = "The value should be set first.";
   get[] := (
      Utils`AssertOrQuit[TrueQ@once, get::errNotSet];
      value);
   secure@Evaluate@get];
createGetSetOnce // Utils`MakeUnknownInputDefinition;
createGetSetOnce ~ SetAttributes ~ {Protected, Locked};

createGetSetOnce /@ {
   {InternalDirectory, _String},
   {SubexpressionRules, {__}},
   {AmplitudeRules, {__}},
   {SettingsFile, {Rule[{__Symbol}, _String]..}}};

setInternalDirectory@DirectoryName@$Input;

`type`vertex = FeynArts`Vertex[_Integer][_Integer];
`type`vertex ~ SetAttributes ~ {Protected, Locked};

`type`propagator = FeynArts`Propagator[ FeynArts`External|FeynArts`Incoming|FeynArts`Outgoing|FeynArts`Internal|FeynArts`Loop[_Integer] ][`type`vertex,`type`vertex,Repeated[FeynArts`Field[_Integer],{0,1}]];
`type`propagator ~ SetAttributes ~ {Protected, Locked};

`type`topology = FeynArts`Topology[_Integer][`type`propagator..];
`type`topology // Protect;

`type`diagram = Rule[`type`topology, FeynArts`Insertions[Generic][__]];

`type`diagramSet = FeynArts`TopologyList[_][`type`diagram..];

`type`amplitude = FeynArts`FeynAmp[
   FeynArts`GraphID[FeynArts`Topology==_Integer,Generic==_Integer],
   Integral[FeynArts`FourMomentum[FeynArts`Internal,_Integer]],
   _,
   {__}->FeynArts`Insertions[FeynArts`Classes][{__}..]
];
`type`amplitudeSet = FeynArts`FeynAmpList[__][`type`amplitude..];

`type`indexCol = FeynArts`Index[Global`Colour,_Integer];
`type`indexGlu = FeynArts`Index[Global`Gluon,_Integer];
`type`indexGeneric = FeynArts`Index[Generic, _Integer];
Module[{filter, rules, indices},
   filter = (FeynArts`Indices -> e_) :> e;
   rules = {  FeynArts`Index -> Identity,
              Global`Colour :> {},
              Global`Gluon :> {}};
   indices = Cases[FeynArts`M$ClassesDescription, filter, Infinity];
   indices = DeleteDuplicates@Flatten[indices //. rules];
   `type`indexGeneration = FeynArts`Index[Alternatives@@indices, _Integer];];

`type`field = FeynArts`S|FeynArts`F|FeynArts`V|FeynArts`U;
`type`genericField = `type`field[`type`indexGeneric];

`type`fc`particle = `type`field[_Integer, Repeated[{_Symbol}, {0, 1}]];
`type`fc`mass = 0|_Symbol|_Symbol@_Symbol;
`type`fc`external = {`type`fc`particle|-`type`fc`particle,
   FormCalc`k@_Integer, `type`fc`mass, {}};
`type`fc`process = {`type`fc`external..} -> {`type`fc`external..};
`type`fc`amplitude = FormCalc`Amp[`type`fc`process][_];
`type`fc`amplitudeSet = {`type`fc`amplitude..};

`type`pickTopoAmp = {Rule[True | False,{__Integer}]..};
`type`saveAmpClass = {Rule[_Integer,{__Integer} | All]..};

`type`massType = Repeated[Alternatives[FeynArts`Loop, FeynArts`Internal], {0, 1}];

`type`genericMass = FeynArts`Mass[
   `type`field[`type`indexGeneric], `type`massType];
`type`specificMass = FeynArts`Mass[
   `type`field[_Integer,
      {Alternatives[`type`indexCol, `type`indexGlu, `type`indexGeneration]..}]];

Get@FileNameJoin@{getInternalDirectory[], "actions.m"};
Get@FileNameJoin@{getInternalDirectory[], "chains.m"};
Get@FileNameJoin@{getInternalDirectory[], "time.m"};
Get@FileNameJoin@{getInternalDirectory[], "topologies.m"};

define[getTopology, {d:`type`diagram} :> First@d];

define[getInsertions, {d:`type`diagram} :> Last@d];

define[indexGeneric, {index:_Integer} :> FeynArts`Index[Generic, index]];

define[genericMass,
   {field:`type`field, index:_Integer} :>
   FeynArts`Mass[field@indexGeneric@index, `type`massType],

   {field:`type`field} :>
   FeynArts`Mass[field@`type`indexGeneric, `type`massType]
];

define[getClassVariables, {amp:`type`amplitude} :> amp[[4, 1]]];

define[getClassInsertions, {amp:`type`amplitude} :> amp[[4, 2]]];

define[getClassRules, {amp:`type`amplitude} :>
   Thread[getClassVariables@amp -> Transpose[List@@getClassInsertions@amp]]
];

define[getProcess,
   {set:`type`diagramSet|`type`amplitudeSet} :>
   Cases[Head@set, (FeynArts`Process -> e:_) :> e][[1]],

   {set:`type`fc`amplitudeSet} :>
   Part[Head@Part[set, 1], 1]
];

getExternalMasses[set:`type`fc`amplitudeSet] :=
   Flatten[List@@getProcess@set, 1][[All, 3]];
getExternalMasses[set:`type`amplitudeSet] :=
   FeynArts`Mass[# /. -1 -> 1] &/@ getField[set, All];
getExternalMasses // secure;

getField[set:`type`diagramSet, i:_Integer] :=
   Flatten[List@@getProcess@set, 1][[i]] /; 0<i<=Plus@@(Length/@getProcess@set);
getField[set:`type`amplitudeSet, In] :=
   First /@ getProcess[set][[1]];
getField[set:`type`amplitudeSet, Out] :=
   First /@ getProcess[set][[2]];
getField[set:`type`amplitudeSet, All] :=
   First /@ Flatten[List @@ getProcess[set], 1];
getField // secure;

fieldPattern[d:Head@`type`diagramSet, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]] /.
      `type`indexGeneration :> Blank[];
fieldPattern[d:Head@`type`diagramSet, a:HoldPattern@Alternatives@__] :=
   fieldPattern[d, #] &/@ a;
fieldPattern // secure;

SetInitialValues[] :=
Module[{fieldNames},
   fieldNames = getFieldNames[];
   setFieldRules@fieldNames;
   setSubexpressionRules@Join[
      getMassRules@fieldNames,
      getFieldRules[],
      getCouplingRules[],
      getGeneralRules[]];
   setAmplitudeRules@Join[
      getSubexpressionRules[],
      getSumRules[],
      {FeynArts`IndexSum -> Sum}];];
SetInitialValues // secure;

getFieldNames::usage = "
@brief Loads and connects particle definitions for SARAH` and FeynArts`.
@returns List of Lists, each contains four String entries:
         1) SARAH` context of particle;
         2) SARAH` name of particle;
         3) FeynArts` type of particle;
         4) FeynArts` integer number of particle.";
define[getFieldNames, {} :>
   Module[{
         regex = "(\\w+): ([SFVU])\\[(\\d+)\\]",
         lines = Utils`ReadLinesInFile@$ParticleFile,
         namespaceRules = Rule[First@#, Sequence[Last@#, First@#]] &/@ Get@$ContextFile,
         names
      },
      names = StringCases[lines, RegularExpression@regex :> {"$1","$2","$3"}] ~ Flatten ~ 1;
      names /. namespaceRules
   ]
];

getMassRules::usage = "
@brief Generates mass replacement rules for a given set of particles.
@param fieldNames Set of field names.
@returns A List of replacements rules for masses.";
define[getMassRules, {fieldNames:{{_, _, _, _}..}} :>
   Module[{
         faMasses = Symbol["Mass" <> #[[2]]] &/@ fieldNames,
         sarahNames = Symbol[#[[1]] <> #[[2]]] &/@ fieldNames,
         massRules
      },
      massRules = MapThread[
         {
            #1[index_] :> SARAH`Mass@#2@{Symbol["SARAH`gt" <> StringTake[SymbolName@index, -1]]},
            #1[indices__] :> SARAH`Mass@#2@indices,
            #1 :> SARAH`Mass@#2
         } &,
         {faMasses, sarahNames}
      ];
      Append[Flatten@massRules,
         FeynArts`Mass[field_, _ : Null] :> SARAH`Mass[field]
      ]
   ]
];

Module[{once, rules},

setFieldRules::errOnce = "The value can be set only once.";
setFieldRules::usage = "
@brief Generates field replacement rules for a given set of particles.
@param fieldNames Set of field names.
@returns A List of replacements rules for fields.";
define[setFieldRules, {fieldNames:{{_, _, _, _}..}} :>
   Module[{
         fullNames = Map[ToExpression, {#[[1]] <> #[[2]], #[[3]], #[[4]]} &/@ fieldNames, 2],
         bose = FeynArts`S|FeynArts`V, fermi = FeynArts`U|FeynArts`F
      },
      Utils`AssertOrQuit[Head@once === Symbol, setFieldRules::errOnce];
      rules = Join[
         # /. {n_, t_, i_} :> Rule[t@i, n],
         # /. {n_, t_, i_} :> RuleDelayed[t[i, {ind__}], n@{ind}],
         # /.
         {
            {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n], Susyno`LieGroups`conj@n],
            {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n], SARAH`bar@n]
         },
         # /.
         {
            {n_, t:bose,  _} :> RuleDelayed[Times[-1,f:n@{ind__}], Susyno`LieGroups`conj@n@{ind}],
            {n_, t:fermi, _} :> RuleDelayed[Times[-1,f:n@{ind__}], SARAH`bar@n@{ind}]
         },
         {
            ind:`type`indexGeneration :> Symbol["SARAH`gt" <> ToString@Last@ind],
            ind:`type`indexCol :> Symbol["SARAH`ct" <> ToString@Last@ind],
            ind:`type`indexGlu :> (Print["Warning: check indexRules of internal.m"];Symbol["SARAH`ct" <> ToString@Last@ind])
         },
         {
            FeynArts`S -> GenericS,
            FeynArts`F -> GenericF,
            FeynArts`V -> GenericV,
            FeynArts`U -> GenericU
         },
         {
            Times[-1,field:_GenericS|_GenericV] :> Susyno`LieGroups`conj@field,
            Times[-1,field:_GenericF|_GenericU] :> SARAH`bar@field
         }
      ] &@ fullNames;
      once = {};
   ]
];

getFieldRules::errNotSet = "The value should be set first.";
define[getFieldRules, {} :> (
   Utils`AssertOrQuit[Head@once =!= Symbol, getFieldRules::errNotSet];
   rules)
];

];

define[getExternalMomentumRules,
   {option:True|False|OperatorsOnly|ExceptLoops, amplitudes:`type`amplitudeSet} :>
   Module[{
         fsFields
      },
      Switch[option,
         True,
            {SARAH`Mom[_Integer,_] :> 0},
         False|OperatorsOnly|ExceptLoops,
            (
               fsFields = getField[amplitudes, All] //. getFieldRules[];
               {SARAH`Mom[i_Integer, lorIndex_] :> SARAH`Mom[fsFields[[i]], lorIndex]}
            )
      ]
   ]
];

define[getSumRules, {} :>
   {
      FeynArts`SumOver[_,_,FeynArts`External] :> Sequence[],
      Times[e:_, FeynArts`SumOver[i:_Symbol, max:_Integer]] :>
         SARAH`sum[i, 1, max, e],
      Times[expr:_, FeynArts`SumOver[index:_Symbol, {min:_Integer, max:_Integer}]] :>
         SARAH`sum[index, min, max, expr],
      SARAH`sum[i:_Symbol, _Integer, max:_Integer, FeynArts`SumOver[_Symbol, max2:_Integer]] :>
         SARAH`sum[i, 1, max, max2],
      SARAH`sum[i:_Symbol, _Integer, max:_Integer, FeynArts`SumOver[_, {min2:_Integer, max2:_Integer}]] :>
         SARAH`sum[i, 1, max, max2-min2]
   }
];

getGeneralRules::usage = "
@brief General translation rules from FeynArts/FormCalc to FlexibleSUSY
       language.
@note See sec. 4.4. of FormCalc manual for details.";
define[getGeneralRules, {} :>
   ({
      FormCalc`Finite -> 1,
      FormCalc`Den[a:_,b:_] :> 1/(a-b),
      FormCalc`Pair[a:_,b:_] :> SARAH`sum[#, 1, 4, SARAH`g[#, #]*Append[a, #]*Append[b, #]],
      f:`type`genericField :> Head[f][GenericIndex@Last@Last@f],
      FormCalc`k[i:_Integer, pairIndex:___] :> SARAH`Mom[i, pairIndex]
   } &@ Unique@"SARAH`lt")
];

Module[{PL, PR, MT, FV, g, md},

   {PL, PR} = Alternatives[
      FeynArts`NonCommutative@Global`ChiralityProjector@#,

      FeynArts`NonCommutative[
         Global`DiracMatrix@FeynArts`KI1@3,
         Global`ChiralityProjector@#
      ]
   ] &/@ {-1, 1};

   Quiet[
      MT[i1:_Symbol, i2:_Symbol] :=
         Global`MetricTensor[FeynArts`KI1[i1:_Integer], FeynArts`KI1[i2:_Integer]];
      If[FormCalc`$FormCalc < 9.7,

         FV[i1:_Symbol, i2:_Symbol, Repeated[_, {0, 1}]] :=
            FeynArts`Mom[i1:_Integer] - FeynArts`Mom[i2:_Integer];,

         FV[i1:_Symbol, i2:_Symbol] := Global`FourVector[
            FeynArts`Mom[i1:_Integer] - FeynArts`Mom[i2:_Integer],
            FeynArts`KI1[3]
         ];
         FV[i1:_Symbol, i2:_Symbol, i3:_Symbol] := Global`FourVector[
            FeynArts`Mom[i1:_Integer] - FeynArts`Mom[i2:_Integer],
            FeynArts`KI1[i3:_Integer]
         ];
      ],
      RuleDelayed::rhs
   ];

   g[f:{__}, i1:_Integer, i2:_Integer] :=
      SARAH`g[LorentzIndex@Part[f, i1], LorentzIndex@Part[f, i2]];

   md[fields:{__}, i1:_Integer, i2:_Integer] :=
      SARAH`Mom@Part[fields, i1] - SARAH`Mom@Part[fields, i2];

   md[fields:{__}, i1:_Integer, i2:_Integer, i3:_Integer] :=
      SARAH`Mom[Part[fields, i1], LorentzIndex@Part[fields, i3]] -
      SARAH`Mom[Part[fields, i2], LorentzIndex@Part[fields, i3]];

   With[{
         f = FeynArts`G[_][0][fields__], s = SARAH`Cp[fields]
      },
      define[getCouplingRules, {} :>
         {
            f@1 :> s@1,
            f@PL :> s@SARAH`PL,
            f@PR :> s@SARAH`PR,
            f@MT[i1, i2] :> s@g[{fields}, i1, i2],
            f@FV[i1, i2] :> s@md[{fields}, i1, i2],

            f[
               FV[i2, i1, i3] * MT[i1, i2] +
               FV[i1, i3, i2] * MT[i1, i3] +
               FV[i3, i2, i1] * MT[i2, i3]
            ] :>
            s[
               md[{fields}, i2, i1, i3] * g[{fields}, i1, i2],
               md[{fields}, i1, i3, i2] * g[{fields}, i1, i3],
               md[{fields}, i3, i2, i1] * g[{fields}, i2, i3]
            ]
         }
      ];
   ];
];

getSettings::usage = "
@brief Loads the file with process-specific settings. If there is no process
       file to load, defines default settings.
@todo Check which settings exist and apply defauld definitions for non-existing
      ones.
@todo Define default first, then redefine during load (pure functions?)";
getSettings[] := Module[{file},
   file = FileNameJoin@{getInternalDirectory[], SymbolName@Head@$Observable,
      "settings.m"};
   BeginPackage@"NPointFunctions`";
   Begin@"`internal`";
   `settings`topology = Default;
   `settings`diagrams = Default;
   `settings`amplitudes = Default;
   `settings`sum = Default;
   `settings`massless = Default;
   `settings`momenta = {};
   `settings`regularization = {};
   `settings`order = Default;
   `settings`chains = Default;
   If[FileExistsQ@file, Get@file;];
   Protect@Evaluate[Context[]<>"settings`*"];
   End[];
   EndPackage[];];
getSettings // secure;

emptyQ::usage = "
@brief Checks, whether the lenght of expression is zero or not. In the latter
       case returns in, otherwise writes an error message and stops
       the evaluation.
@param input Any expression.
@returns An input in case of non-zero length.";
emptyQ[input:_] := With[{sym = Head@Unevaluated@input},
   If[Length@input =!= 0,
      input,
      sym::errEmpty = "The result is empty.";
      Utils`AssertOrQuit[_, sym::errEmpty];]];
emptyQ // Utils`MakeUnknownInputDefinition;
emptyQ ~ SetAttributes ~ {Protected, HoldFirst};

NPointFunctionFAFC::usage = "
@brief Applies FeynArts` routines for a given process, preparing it for
       FormCalc`.
@returns A structure, representing NPF object.
@todo If topologies are not generated, then check and return.";
NPointFunctionFAFC[inFields_, outFields_] :=
Module[{topologies, diagrams, amplitudes},
   getSettings[];
   topologies = emptyQ@FeynArts`CreateTopologies[
      $LoopLevel,
      Length@inFields -> Length@outFields,
      FeynArts`ExcludeTopologies -> getExcludeTopologies[]];
   diagrams = emptyQ@FeynArts`InsertFields[topologies, inFields -> outFields];
   diagrams = modify@diagrams;

   amplitudes = FeynArts`CreateFeynAmp@diagrams;
   {diagrams, amplitudes} = modify[diagrams, amplitudes];

   debugMakePictures[
      diagrams,
      StringJoin[ToString /@ (Join[getField[amplitudes, In],getField[amplitudes, Out]] //. getFieldRules[] /. e_[{_}]:>e)]
   ];

   {
      {getField[amplitudes, In], getField[amplitudes, Out]} //. getFieldRules[],
      calculateAmplitudes[diagrams, amplitudes]
   }
];

getSumSettings::usage = "
@brief Some topologies can lead to physically incorrect summation on C++ level.
       One can use `settings`sum to specify, which fields to skip.
@param ds A set of diagrams.
@returns A set of restrictions for generic sums.";
getSumSettings[ds:`type`diagramSet] :=
Module[{repl = {}, res},
   sum = If[#===Default, {}, #]&@`settings`sum;
   Do[repl = Join[repl, applyAction[ds, ac]];, {ac, sum}];
   res = List@@(ds /. repl);
   res = If[MatchQ[#, `type`diagram],
            Table[{}, {Length@getInsertions@#}],
            First@#] &/@ res;
   Flatten[res, 1]];
getSumSettings // secure;

collectSame::usage = "
@brief Finds the same keys in the list of rules and for them collects RHSs
       into one list, i.e.:
          collectSame@{a->{1, 1}, b->{2, 2, 2}, a->{3, 3}} leads to
          {a->{{1, 3}, {1, 3}}, b->{{2, 2, 2}}}.
@param list A list of rules.
@return A list of rules.";
define[collectSame, {list:{Rule[_, {__}]...}} :>
   Module[{
         tally = Tally[First /@ list]
      },
      If[#2 == 1,
         {#1/.list},
         #1 -> Transpose@Cases[list, Rule[#1, el_] :> el]
      ]&@@#&/@tally
   ]
];

getMasslessSettings::usage = "
@brief In some topologies field insertions can lead to physically incorrect
       simplifications if they are done naively. This function provides
       rules in order to prevent this.
@param diagrams A set of diagrams.
@returns A set of rules for amplitudes.";
getMasslessSettings[diagrams:`type`diagramSet] :=
Module[{parse, rules, set},
   set = If[Default===#, {}, #]&@`settings`massless;
   parse = If[MatchQ[#1, `type`diagram],
      List@@MapIndexed[{}&, getInsertions@#1],
      First@#1]&;
   rules = MapIndexed[applyAction[diagrams, #1]&, set];
   Flatten[List@@MapIndexed[parse, diagrams /. collectSame@rules], 1]];
getMasslessSettings // secure;

getRegularizationSettings::usage = "
@brief Some amplitudes are calculated incorrectly in some schemes (like box
       diagrams in CDR). For handling this, one can overwrite used scheme for
       some topologies.
@param diagrams A set of diagrams.
@returns A set of settings for FormCalc`Dimension.
@todo Unify with getMomSettings.";
getRegularizationSettings::errOverlap = "
Topology rules in `.`settings`.`regularization overlap.";
getRegularizationSettings[diagrams:`type`diagramSet] :=
Module[{scheme, replacements, f},
   scheme = Switch[$Scheme, DimensionalReduction, 4,
                            DimensionalRegularization, D];
   f = (getAmplitudeNumbers[diagrams, First@#] /. x:_Integer :> Last@#) &;
   If[`settings`regularization === {},
      Array[scheme&, getClassAmount@diagrams],
      replacements = Transpose[f /@ `settings`regularization];
      Flatten[Switch[ Count[First/@#,True],
         0, Array[scheme&, Length[False /. #]],
         1, True /. #,
         _, Utils`AssertOrQuit[False, getRegularizationSettings::errOverlap]
         ] &/@ replacements]]];
getRegularizationSettings // secure;

getMomSettings::usage = "
@brief Uses settings to eliminale specific momenta in specific topologies.
@param diagrams A set of topologies with class insertions.
@returns A List of option values for FormCalc`MomElim for every generic
         amplitude.";
getMomSettings::errOverlap = "
Topology rules in `.`settings`.`momenta overlap.";
define[getMomSettings, {diagrams:`type`diagramSet} :>
   Module[{
         replacements,
         f = (getAmplitudeNumbers[diagrams, First@#] /. x:_Integer :> Last@#) &
      },
      If[`settings`momenta === {},
         Array[Automatic&, getClassAmount@diagrams],
         replacements = Transpose[f /@ `settings`momenta];
         Flatten[Switch[ Count[First/@#,True],
            0, Array[Automatic&, Length[False /. #]],
            1, True /. #,
            _, Utils`AssertOrQuit[False, getMomSettings::errOverlap]
            ] &/@ replacements
         ]
      ]
   ]
];

modify::usage = "
@brief Changes amplitudes and diagrams according to `settings`diagrams and
       `settings`amplitudes.
@param diagrams A set of diagrams.
@param amplitudes A set of amplitudes.
@returns A modified set of diagrams and amplitudes.";
modify[diagrams:`type`diagramSet] :=
Module[{d = diagrams},
   Do[d = applyAction[d, ac];, {ac, getActions@`settings`diagrams}];
   d];
modify[diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet] :=
Module[{d = diagrams, a = removeColours@amplitudes},
   Do[{d, a} = applyAction[{d, a}, ac];, {ac, getActions@`settings`amplitudes}];
   {d, a}];
modify // secure;

removeColours[i:`type`amplitudeSet] :=
Delete[i, Position[i, FeynArts`Index[Global`Colour, _Integer]]];
removeColours // secure;

getAmplitudeNumbers::usage = "
@brief Gives numbers of amplitudes, accepted by a criterion on topology.
@param diagrams A set of diagrams.
@param <one argument function> critFunction Function for topology selection. If
       critFunction[<topology>] gives True, then topology is accepted.
@returns {<Rule>} List of rules of the form <boolean>->{<integer>..}. LHS
         stands for the topology, RHS gives numbers of classes (and the numbers
         of amplitudes the same time).";
define[getAmplitudeNumbers, {diagrams:`type`diagramSet, critFunction_} :>
   Module[
      {
         topologies = List@@First/@diagrams,
         genNums = Length/@(List@@Last/@diagrams),
         numRegions,takeOrNot
      },
      numRegions = Array[Range[Plus@@genNums[[1;;#-1]]+1,Plus@@genNums[[1;;#]]]&,Length@genNums];
      takeOrNot = Array[TrueQ@critFunction@Part[topologies,#]&,Length@topologies];
      MapThread[#1->#2&,{takeOrNot,numRegions}]
   ]
];

define[printDiagramsInfo, {diagrams:`type`diagramSet, where_String:" "} :>
   Module[{
         nGeneric = Length@Cases[diagrams,Generic==_Integer:>1,Infinity,Heads -> True],
         nClasses = getClassAmount@diagrams
      },
      Print[where,"in total: ",nGeneric," Generic, ",nClasses," Classes insertions"];
   ]
];

define[getClassAmount, {set:`type`diagramSet} :>
   Length@Cases[set, FeynArts`Classes==_Integer:>1, Infinity, Heads -> True]
];

define[debugMakePictures, {diagrams:`type`diagramSet, name_String:"classes"} :>
   Module[{
         out = {}, directory = FileNameJoin[Most[FileNameSplit@@FeynArts`$Model]]
      },
      DeleteFile@FileNames@"*class_*.jpg";
      Export["class.jpg",FeynArts`Paint[diagrams,
         FeynArts`PaintLevel -> {FeynArts`Classes},
         FeynArts`ColumnsXRows -> 1,
         FeynArts`SheetHeader -> None,
         FeynArts`Numbering -> FeynArts`Simple,
         DisplayFunction :> (AppendTo[out, #] &/@ Render[##, "JPG"] &)
      ]];

      Put[out, FileNameJoin@{directory, name<>".m"}];
   ]
];

getFieldInsertions::usage = "
@brief Applies FindGenericInsertions[] to a set of diagrams or one.
@param set A set of diagrams.
@param diag A single diagram.
@param numQ Responsible for the type of output field names.
@returns For a single diagram returns List (for a given topology) of List
         (for all generic fields) of List (for all class fields) of rules
         {{{x->y,..},..},..}. For a set of diagrams, this construct is further
         transformed.";
define[getFieldInsertions,
   {set:`type`diagramSet} :>
   (Map[Last, #, {3}] &@ Flatten[ getFieldInsertions /@ (List @@ set), 1]),

   {diag:`type`diagram, numQ:True|False:False} :>
   (FindGenericInsertions[#, numQ] &/@ Apply[List, getInsertions@diag, {0, 1}])
];

FindGenericInsertions::usage = "
@brief generic FeynmanGraph has rules Field[num]->particleType,
        class FeynmanGraph has rules Field[num]->particleClass.
        This function gives pairs particleType[gen,num]->particleClass, avoiding
        Field[_] mediator (if keepFieldNum==True then Field[_]->particleClass
        is given).
@param 1st argument is of the form
       {FeynmanGraph[__][__],Insertions[Classes][__]}.
@param 2nd argument changes the type of output field names
       True gives Field[_] names, False gives particleClass names.
@returns list (sorted; for all generic fields) of list (for all class fields)
         of rules {{x->y,..},..}.
@note this function is called by GenericInsertionsForDiagram[].
@note this function doesn't look at external particles.
@note all indices in rhs. of rules are removed.";
define[FindGenericInsertions, {{graphGen_,insertCl_}, keepFieldNum_} :>
   Module[{
         toGenericIndexConventionRules = Cases[graphGen,
            Rule[FeynArts`Field[index_Integer],type_Symbol] :>
            Rule[FeynArts`Field@index, type[FeynArts`Index[Generic,index]]]
         ],
         fieldsGen, genericInsertions
      },
      fieldsGen = toGenericIndexConventionRules[[All,1]];
      genericInsertions = Cases[#,
         Rule[genericField_,classesField_] /; MemberQ[fieldsGen, genericField] :>
         Rule[genericField, StripParticleIndices@classesField]] &/@ insertCl;
      SortBy[#,First]&/@ If[keepFieldNum,
         List @@ genericInsertions,
         List @@ genericInsertions /. toGenericIndexConventionRules
      ]
   ]
];

StripParticleIndices::usage = "
@brief Removes particle indices from a given (possibley generic) field.
@param field the given field.
@returns The given field with all indices removed.";
define[StripParticleIndices,
   {Times[-1,field_]} :>
   Times[-1, StripParticleIndices[field]],

   {genericType_[classIndex_, ___]} :>
   genericType[classIndex]
];

getColourFactors::usage = "
@brief Creates colour factors for a given diagram.
@param ds A diagram set.
@param diagram A diagram to work with.
@returns List (for a given topology) of lists (for all generic fields) of
         colour factors.
@note During generation of genericDiagram at 1-loop level the ii-type loop
      propagators have the largest number because of FeynArts.
@note In seqProp numbers of the first vertices inside propagators are sorted
      by FeynArts.
@note External fields always come at first places in adjacency matrix.
@note This function doesn't know anything about CXXDiagrams`.` context.";
define[getColourFactors,
   {ds:`type`diagramSet} :>
   (Flatten[getColourFactors /@ (List @@ ds), 1] //. getFieldRules[]),

   {diagram:(_[_][seqProp__]->_[_][_[__][rulesFields__]->_,___])} :>
   Module[{
         propPatt,adjacencyMatrix,externalRules,genericDiagram,genericInsertions
      },
      propPatt[i_, j_, f_] := _[_][_[_][i], _[_][j], f];

      adjacencyMatrix = Module[
         {adjs = Tally[{seqProp}/.propPatt[i_,j_,_]:>{{i,j},{j,i}}] },
         Normal@SparseArray@Flatten[{#[[1,1]]->#[[2]],#[[1,2]]->#[[2]]} &/@ adjs]];

      externalRules = Cases[{rulesFields}, HoldPattern[_[_]->_Symbol[__]]];

      genericDiagram = Module[
         {fld = Flatten[{seqProp}/.propPatt[i_,j_,f_]:>{{j,i,-f},{i,j,f}}, 1] },
         GatherBy[SortBy[fld,First],First] /. {_Integer, _Integer, f_} :> f
         ] /. Join[ {#} -> # &/@ externalRules[[All, 1]]];

      genericInsertions = getFieldInsertions[diagram, True];

      Map[CXXDiagrams`ColourFactorForIndexedDiagramFromGraph[
         CXXDiagrams`IndexDiagramFromGraph[
            genericDiagram /. externalRules /. #, adjacencyMatrix],
         adjacencyMatrix] &,
         genericInsertions,
         {2}
      ]
   ]
];

getFermionOrder::usage = "
@brief Returns the order of fermions for FormCalc`FermionOrder option. Default
       is a reversed one. Can be overwritten by `settings`order.
@param expression A set of diagrams.
@returns A list of integers, representing an order of fermions.";
getFermionOrder[expression:`type`diagramSet] :=
Switch[`settings`order,
   Default, Reverse@Range[Plus@@Length/@getProcess@expression],
   _, `settings`order];
getFermionOrder // secure;

calculateAmplitudes::usage = "
@brief Applies FormCalc` routines to amplitude set, simplifies the result.
@param diagrams A set of diagrams.
@param amplitudes A of amplitudes (without colours).
@returns The main part of NPF object, containing: generic amplitudes,
         class specific insertions, subexpressions.";
calculateAmplitudes[diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet] :=
Module[{
      proc = getProcess@amplitudes,
      genericInsertions = getFieldInsertions@diagrams,

      combinatorialFactors = CombinatorialFactorsForClasses /@ List@@amplitudes,
      ampsGen = FeynArts`PickLevel[Generic][amplitudes],
      feynAmps, generic, chains, subs,
      zeroedRules
   },
   If[$ZeroMomenta,
      ampsGen = FormCalc`OffShell[ampsGen,
         Sequence@@Array[#->0&, Plus@@Length/@proc]]];

   feynAmps = mapThread[
      FormCalc`CalcFeynAmp[Head[ampsGen][#1],
         FormCalc`Dimension -> #2,
         FormCalc`OnShell -> $OnShell,
         FormCalc`FermionChains -> FormCalc`Chiral,
         FormCalc`FermionOrder -> getFermionOrder@diagrams,
         FormCalc`Invariants -> False,
         FormCalc`MomElim -> #3]&,
      {  ampsGen,
         getRegularizationSettings@#,
         getMomSettings@#} &@ diagrams,
      "Amplitude calculation"] //. FormCalc`GenericList[];
   generic = MapThread[getGenericSum, {feynAmps, getSumSettings@diagrams}];

   {generic, chains, subs} = proceedChains[diagrams, amplitudes, generic];

   setZeroMassRules@{amplitudes, feynAmps};
   {generic, chains, subs} = makeMassesZero[
      {generic, chains, subs},
      diagrams,
      $ZeroMomenta];

   FCAmplitudesToFSConvention[
      {  generic,
         genericInsertions,
         combinatorialFactors,
         getColourFactors@diagrams},
      chains,
      subs] /. getExternalMomentumRules[$ZeroMomenta, amplitudes]];
calculatedAmplitudes // secure;

setZeroMassRules::usage = "
@brief For a given sets of FeynArts` amd FormCalc` amplitudes creates rules to
       nullify masses of external particles.
@param fa A set of FeynArts` amplitudes.
@param fc A set of FormCalc` amplitudes.
@returns Null.
@note Both of sets are required, because, unfortunately, FormCalc` introduces
      new abbreviations, which mix with FeynArts` ones.
@note Amplitudes are taken, because they do not have colour structures already.
@note Explicit names for masses are expected only for external particles.";
getZeroMassRules::errNotSet = "
Call setZeroMassRules to set up rules first.";
getZeroMassRules::usage = "
@brief Returns a set of rules to nullify masses of external particles.
@return A list of rules to nullify masses of external particles.
@note Rules of external particle #i are under numbers (2*#i) and (2*#i-1).";
Module[{rules},
   setZeroMassRules[{fa:`type`amplitudeSet, fc:`type`fc`amplitudeSet}] :=
      rules = RuleDelayed[#, 0] &/@
         Riffle[getExternalMasses@fa, getExternalMasses@fc];
   setZeroMassRules // secure;
   getZeroMassRules[] := (
      Utils`AssertOrQuit[Head@rules =!= Symbol, getZeroMassRules::errNotSet];
      rules);
   getZeroMassRules // secure;];

makeMassesZero::usage = "
@brief Sets the masses of external particles to zero everywhere, except loop
       integrals, applies subexpressions.
@param generic An expression to modify.
@param chains A lis with fermionic chains.
@param subs A list of subexpressions.
@param diagrams A set of diagrams.
@returns A list with modified expression, chains and empty non-applied subexpressions.";
makeMassesZero[
   {generic_, chains_, subs_}, diagrams:`type`diagramSet, ExceptLoops] :=
Module[{funcs, names, pattern, uniqueIntegrals, hideInt, showInt, rules, new},
   funcs = getMasslessSettings@diagrams;
   names = ToExpression/@Names@RegularExpression@"LoopTools`[ABCD]\\d+i*";
   new = generic //. subs;
   pattern = Alternatives @@ ( #[__] &/@ names );
   uniqueIntegrals = DeleteDuplicates@Cases[new, pattern, Infinity];
   hideInt = Rule[#, Unique@"loopIntegral"] &/@ uniqueIntegrals;
   showInt = hideInt /. Rule[x_, y_] -> Rule[y, x];
   rules = List@@MapIndexed[Composition[Sequence@@funcs[[#2[[1]]]]]@
      getZeroMassRules[]&, new];
   {  List@@MapIndexed[#1 //. rules[[#2[[1]]]] /. showInt&, new /. hideInt /.
         FormCalc`Pair[_,_] -> 0],
      setZeroExternalMomentaInChains@chains /. getZeroMassRules[],
      {}}];
makeMassesZero[{generic_, chains_, subs_}, diagrams:`type`diagramSet, True] :=
Module[{
      zeroedRules = Cases[subs, Rule[_, pair:FormCalc`Pair[_, _]] :> (pair->0)], new
   },
   {new, zeroedRules} = ZeroRules[subs, zeroedRules];
   {  generic /. zeroedRules,
      setZeroExternalMomentaInChains@chains,
      new
   }];
makeMassesZero[{expr_, chains_, subs_}, diagrams:`type`diagramSet, _] :=
   {expr, chains, subs};
makeMassesZero // secure;

mapThread::usage = "
@brief Maps a function onto multiple sets of equal length, accompanying it by
       printing a progress bar.
@param func A function to apply to set of data.
@param exprs A list of listable sets with data.
@param text A string to be printed.
@todo Add check for equality of length for exprs.";
define[mapThread, {func_, exprs:{__}, text:_String:""} :>
   Module[{
         sr, print, percent, init, end, bar, def = 70, dots,
         tot = Length@First@exprs, result
      },
      `time`set[];
      subWrite["\n"<>text<>" ...\n"];

      sr[str:_String, num:_Integer] := StringJoin@@Array[str&, num];

      print = (
         percent = #/tot;
         init = "["<>ToString@now<>"/"<>ToString@tot<>"] [";
         end = "] "<>ToString@Floor[100*percent]<>"%";
         bar = def - StringLength[init<>end];
         dots = Floor[bar*percent];
         subWrite@StringJoin[init, sr[".", dots], sr[" ", bar-dots],end,"\r"];
      )&;

      result = Table[(print@now; func@@(#[[now]]&/@ exprs)), {now, tot}];

      subWrite@"\033[K\033[A";
      subWrite[text<>" ... done in "<>`time`get[]<>" seconds.\n"];
      result
   ]
];

CombinatorialFactorsForClasses::usage="
@brief Takes generic amplitude and finds numerical combinatirical factors
       which arise at class level.
@returns list of combinatorical factors for a given generic amplitude
@param FeynArts`.`FeynAmp[__]";
define[CombinatorialFactorsForClasses, {FeynArts`FeynAmp[_,_,_,rules_->_[_][classReplacements__]]} :>
   ({classReplacements}[[ All,#[[1,1]] ]] /.
      {
         FeynArts`IndexDelta[___] -> 1,
         FeynArts`SumOver[__] -> 1
      } &@ Position[rules, FeynArts`RelativeCF])
];

getGenericFields::usage = "
@brief Generates a list of unique sorted generic fields in expression.
@param expr An expression, where to search.
@returns A list of unique sorted generic fields.";
define[getGenericFields, {expr:_} :>
   Sort@DeleteDuplicates[Cases[expr, `type`genericField, Infinity]]
];

getGenericSum::usage= "
@brief Converts FormCalc`Amp into NPointFunctions`GenericSum object using
       restriction rules for generic fields.
@param amplitude FormCalc`Amp expression.
@param sumRules A set of rules, restricting the summation.
@returns A NPointFunctions`GenericSum object.";
getGenericSum[amplitude:`type`fc`amplitude, sumRules:{Rule[_Integer, _]...}] :=
Module[{sort, rules},
   sort = getGenericFields@amplitude;
   rules = Append[sumRules, _Integer -> False];
   GenericSum[
      List@@amplitude,
      sort /. f_[_[_,i_]] :> {f@GenericIndex@i, i /. rules}]];
getGenericSum // Utils`MakeUnknownInputDefinition;
getGenericSum // Protect;

ZeroRules::usage = "
@brief Given a set of rules that map to zero and a set that does
       not map to zero, apply the zero rules to the non-zero ones
       recursively until the non-zero rules do not change anymore.
@param nonzeroRules The list of nonzero rules.
@param zeroRules The list of zero rules.
@returns a list of rules that map the same expressions as the initial rules.";
define[ZeroRules, {nonzeroRules:{Rule[_,_]...}, zeroRules:{Rule[_,0]...}} :>
   Module[{newNonzero, newZeroRules},
      newNonzero = Thread[
         Rule[nonzeroRules[[All,1]],nonzeroRules[[All,2]] /. zeroRules]];

      If[newNonzero === nonzeroRules, Return[{nonzeroRules, zeroRules}]];

      newZeroRules = Cases[newNonzero,HoldPattern[_->0]];
      newNonzero = Complement[newNonzero, newZeroRules];

      ZeroRules[newNonzero, Join[zeroRules,newZeroRules]]
   ]
];

FCAmplitudesToFSConvention::usage=
"@brief Translate a list of FormCalc amplitudes and their abbreviations and
        subexpressions into FlexibleSUSY language.
@param amplitudes The given list of amplitudes
@param abbreviations A list of abbreviations
@param subexpressions A list of subexpressions
@returns A list of amplitudes and joined abbreviations and subexpressions.";
define[FCAmplitudesToFSConvention,
   {amplitudes_, abbreviations_, subexpressions_} :>
   Module[{
         fsAmplitudes = amplitudes //. getAmplitudeRules[],
         fsAbbreviations = abbreviations //. getSubexpressionRules[] //. {
            FormCalc`DiracChain -> NPointFunctions`internal`dc,
            FormCalc`Spinor -> SARAH`DiracSpinor,
            FormCalc`Lor -> SARAH`Lorentz
         },
         fsSubexpressions = subexpressions //. getSubexpressionRules[]
      },
      {fsAmplitudes, Join[fsAbbreviations,fsSubexpressions]}
   ]
];

End[];
EndPackage[];
