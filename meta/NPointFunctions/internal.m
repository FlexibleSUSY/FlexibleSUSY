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

Block[{Format},
   Needs@"FeynArts`";
   FeynArts`$FAVerbose = 0;
   FeynArts`InitializeModel@NPointFunctions`$FAMod;
   SetOptions[FeynArts`InsertFields,
      FeynArts`Model -> NPointFunctions`$FAMod,
      FeynArts`InsertionLevel -> FeynArts`Classes];
   Needs@"FormCalc`";
   Print[];
   FormCalc`$FCVerbose = 0;
   (  If[!DirectoryQ@#, CreateDirectory@#];
      SetDirectory@#;)&@NPointFunctions`$FCDir;];

Needs@"Utils`";

BeginPackage@"NPointFunctions`";

Off[General::shdw];
{  DiracChain, Mat};
On[General::shdw];

{  NPointFunctionFAFC};
{  LorentzIndex, GenericSum, GenericIndex,
   GenericS, GenericF, GenericV, GenericU,
   OperatorsOnly, ExceptLoops} ~ SetAttributes ~ {Protected};

Begin@"`Private`";

secure[sym:_Symbol] :=
Protect@Evaluate@Utils`MakeUnknownInputDefinition@sym;
secure // secure;

$InternalDirectory = DirectoryName@$Input;
Get@FileNameJoin@{$InternalDirectory, #<>".m"}&/@
   {"type", "rules", "actions", "chains", "topologies"};

getTopology[d:`type`diagram] := First@d;
getTopology // secure;

getInsertions[d:`type`diagram] := Last@d;
getInsertions // secure;

indexGeneric[index:_Integer] := FeynArts`Index[Generic, index];
indexGeneric // secure;

genericMass[field:`type`field, index:_Integer] :=
   FeynArts`Mass[field@indexGeneric@index, `type`massType];
genericMass[field:`type`field] :=
   FeynArts`Mass[field@`type`indexGeneric, `type`massType];
genericMass // secure;

getProcess[set:`type`diagramSet|`type`amplitudeSet] :=
   Cases[Head@set, (FeynArts`Process -> e:_) :> e][[1]];
getProcess[set:`type`fc`amplitudeSet] :=
   Part[Head@Part[set, 1], 1];
getProcess // secure;

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

getExternalMomentumRules[option:True|False|OperatorsOnly|ExceptLoops,
   amplitudes:`type`amplitudeSet] :=
Module[{fsFields},
   Switch[option,
      True,
         {SARAH`Mom[_Integer,_] :> 0},
      False|OperatorsOnly|ExceptLoops,
         fsFields = getField[amplitudes, All] //. $FieldRules;
         {SARAH`Mom[i_Integer, lorIndex_] :>
            SARAH`Mom[fsFields[[i]], lorIndex]}]];
getExternalMomentumRules // secure;

getSettings::usage = "
@brief Loads the file with process-specific settings. If there is no process
       file to load, defines default settings.";
getSettings[] :=
(  BeginPackage@"NPointFunctions`";
   Begin@"`Private`";
   `settings`topology = Default;
   `settings`diagrams = Default;
   `settings`amplitudes = Default;
   `settings`sum = Default;
   `settings`massless = Default;
   `settings`momenta = Default;
   `settings`regularization = Default;
   `settings`order = Default;
   `settings`chains = Default;
   If[FileExistsQ@#, Get@#;]&@FileNameJoin@
      {$InternalDirectory, SymbolName@Head@$Observable, "settings.m"};
   Protect@Evaluate[Context[]<>"settings`*"];
   End[];
   EndPackage[];);
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
@brief Applies ``FeynArts`` routines for a given process, preparing it for
       ``FormCalc``.
@returns A structure, representing n-point function object.
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

   debugMakePictures[diagrams, amplitudes];

   {  {  getField[amplitudes, In], getField[amplitudes, Out]} //. $FieldRules,
         calculateAmplitudes[diagrams, amplitudes]}];
NPointFunctionFAFC // secure;

getSumSettings::usage = "
@brief Some topologies can lead to physically incorrect summation on
       C++ level.
       One can use ``settings`sum`` to specify, which fields to skip.
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
       into one list, i.e.::

          collectSame@{a->{1, 1}, b->{2, 2, 2}, a->{3, 3}}

       leads to::

          {a->{{1, 3}, {1, 3}}, b->{{2, 2, 2}}}.
@param list A list of rules.
@returns A list of rules.";
collectSame[list:{Rule[_, {__}]...}] :=
Module[{tally = Tally[First /@ list]},
   If[#2 == 1,
      {#1/.list},
      #1 -> Transpose@Cases[list, Rule[#1, el_] :> el]]&@@#&/@tally];
collectSame // secure;

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
@returns A set of settings for ``FormCalc`Dimension``.
@todo Unify with ``getMomSettings``.";
getRegularizationSettings::errOverlap = "
Topology rules in `.`settings`.`regularization overlap.";
getRegularizationSettings[diagrams:`type`diagramSet] :=
Module[{scheme, replacements, f},
   scheme = Switch[$Scheme, FlexibleSUSY`DRbar, 4,
                            FlexibleSUSY`MSbar, D];
   f = (getAmplitudeNumbers[diagrams, First@#] /. x:_Integer :> Last@#) &;
   If[`settings`regularization === Default,
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
@returns A List of option values for ``FormCalc`MomElim`` for every generic
         amplitude.";
getMomSettings::errOverlap = "
Topology rules in `.`settings`.`momenta overlap.";
getMomSettings[diagrams:`type`diagramSet] :=
Module[{replacements, f},
   f = (getAmplitudeNumbers[diagrams, First@#] /. x:_Integer :> Last@#) &;
   If[`settings`momenta === Default,
      Array[Automatic&, getClassAmount@diagrams],
      replacements = Transpose[f /@ `settings`momenta];
      Flatten[Switch[ Count[First/@#,True],
         0, Array[Automatic&, Length[False /. #]],
         1, True /. #,
         _, Utils`AssertOrQuit[False, getMomSettings::errOverlap]
         ] &/@ replacements]]];
getMomSettings // secure;

modify::usage = "
@brief Changes amplitudes and diagrams according to ```settings`diagrams``
       and ```settings`amplitudes``.
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
@param critFunction A function of one argumentfor topology selection.
       If critFunction gives True, then topology is accepted.
@returns List of rules of the form::

            <boolean> -> {<integer>..}

         LHS stands for the topology, RHS gives numbers of classes
         (and the numbers of amplitudes the same time).";
getAmplitudeNumbers[diagrams:`type`diagramSet, critFunction_] :=
Module[{topologies, genNums, numRegions, takeOrNot},
   topologies = List@@First/@diagrams;
   genNums = Length/@(List@@Last/@diagrams);
   numRegions = Array[Range[Plus@@genNums[[1;;#-1]]+1,Plus@@genNums[[1;;#]]]&,Length@genNums];
   takeOrNot = Array[TrueQ@critFunction@Part[topologies,#]&,Length@topologies];
   MapThread[#1->#2&,{takeOrNot,numRegions}]];
getAmplitudeNumbers // secure;

printDiagramsInfo[diagrams:`type`diagramSet, where_String:" "] :=
Module[{nGeneric, nClasses},
   nGeneric = Length@Cases[diagrams,Generic==_Integer:>1,Infinity,Heads -> True];
   nClasses = getClassAmount@diagrams;
   Print[where,"in total: ",nGeneric," Generic, ",nClasses," Classes insertions"];];
printDiagramsInfo // secure;

getClassAmount[set:`type`diagramSet] :=
   Length@Cases[set, FeynArts`Classes==_Integer:>1, Infinity, Heads -> True];
getClassAmount // secure;

debugMakePictures[diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet] :=
Module[{out = {}, directory, name},
   name = StringJoin[ToString /@ (Join[getField[amplitudes, All],
      $Processes] //. $FieldRules /. e_[{_}]:>e)];
   directory = DirectoryName[FeynArts`$Model<>".mod"];
   FeynArts`Paint[diagrams,
      FeynArts`PaintLevel -> {Generic},
      FeynArts`ColumnsXRows -> 1,
      FeynArts`FieldNumbers -> True,
      FeynArts`SheetHeader -> None,
      FeynArts`Numbering -> FeynArts`Simple,
      DisplayFunction :> (AppendTo[out, #] &/@ Render[##, "JPG"] &)];
   Put[out, FileNameJoin@{directory, name<>".m"}]];
debugMakePictures // secure;

getFieldInsertions::usage = "
@brief Applies ``FindGenericInsertions`` to a set of diagrams or one.
@param set A set of diagrams.
@param diag A single diagram.
@param numQ Responsible for the type of output field names.
@returns For a single diagram returns ``List`` (for a given topology) of
         ``List`` (for all generic fields) of ``List``
         (for all class fields) of rules
         ``{{{x->y,..},..},..}``. For a set of diagrams, this construct
         is further transformed.";
getFieldInsertions[set:`type`diagramSet] :=
   Map[Last, #, {3}] &@ Flatten[ getFieldInsertions /@ (List @@ set), 1];
getFieldInsertions[diag:`type`diagram, numQ:True|False:False] :=
   FindGenericInsertions[#, numQ] &/@ Apply[List, getInsertions@diag, {0, 1}];
getFieldInsertions // secure;

FindGenericInsertions::usage = "
@brief generic ``FeynmanGraph`` has rules ``Field[num]->particleType``,
       class ``FeynmanGraph`` has rules ``Field[num]->particleClass``.
       This function gives pairs ``particleType[gen,num]->particleClass``,
       avoiding ``Field[_]`` mediator (if ``keepFieldNum==True``, then
       ``Field[_]->particleClass`` is given).
@param graphGen ``FeynmanGraph[__][__]``.
@param insertCl ``Insertions[Classes][__]``.
@param keepFieldNum Changes the type of output field names.
       ``True`` gives ``Field[_]`` names, ``False`` gives ``particleClass``
       names.
@returns A sorted ``List`` for all generic fields of ``List`` for all
         class fields of rules ``{{x->y,..},..}``.
@note This function doesn't look at external particles.
@note All indices in rhs. of rules are removed.";
FindGenericInsertions[{graphGen_,insertCl_}, keepFieldNum_] :=
Module[{toGenericIndexConventionRules, fieldsGen, genericInsertions},
   toGenericIndexConventionRules = Cases[graphGen,
      Rule[FeynArts`Field[index_Integer],type_Symbol] :>
      Rule[FeynArts`Field@index, type[FeynArts`Index[Generic,index]]]];
   fieldsGen = toGenericIndexConventionRules[[All,1]];
   genericInsertions = Cases[#,
      Rule[genericField_,classesField_] /; MemberQ[fieldsGen, genericField] :>
      Rule[genericField, StripParticleIndices@classesField]] &/@ insertCl;
   SortBy[#,First]&/@ If[keepFieldNum,
      List @@ genericInsertions,
      List @@ genericInsertions /. toGenericIndexConventionRules]];
FindGenericInsertions // secure;

StripParticleIndices::usage = "
@brief Removes particle indices from a given (possibley generic) field.
@param field the given field.
@returns The given field with all indices removed.";
StripParticleIndices[Times[-1,field_]] :=
   Times[-1, StripParticleIndices[field]];
StripParticleIndices[genericType_[classIndex_, ___]] :=
   genericType@classIndex;
StripParticleIndices // secure;

getColourFactors::usage = "
@brief Creates colour factors for a given diagram.
@param ds A diagram set.
@param diagram A diagram to work with.
@returns List (for a given topology) of lists (for all generic fields) of
         colour factors.
@note During generation of genericDiagram at 1-loop level the ii-type loop
      propagators have the largest number because of ``FeynArts``.
@note In seqProp numbers of the first vertices inside propagators are sorted
      by ``FeynArts``.
@note External fields always come at first places in adjacency matrix.
@note This function doesn't know anything about ``CXXDiagrams`` context.";
getColourFactors[ds:`type`diagramSet] :=
   Flatten[getColourFactors /@ (List @@ ds), 1] //. $FieldRules;
getColourFactors[diagram:(_[_][seqProp__]->_[_][_[__][rulesFields__]->_,___])] :=
Module[{propPatt, adjacencyMatrix, externalRules, genericDiagram,
      genericInsertions},
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
      {2}]];
getColourFactors // secure;

getFermionOrder::usage = "
@brief Returns the order of fermions for ``FermionOrder`` option.
       Default is a reversed one. Can be overwritten by ``settings`order``.
@param expression A set of diagrams.
@returns A list of integers, representing an order of fermions.";
getFermionOrder[expression:`type`diagramSet] :=
Switch[`settings`order,
   Default, Reverse@Range[Plus@@Length/@getProcess@expression],
   _, `settings`order];
getFermionOrder // secure;

calculateAmplitudes::usage = "
@brief Applies ``FormCalc`` routines to amplitude set, simplifies the result.
@param diagrams A set of diagrams.
@param amplitudes A of amplitudes (without colours).
@returns The main part of n-point function object, containing:

         * generic amplitudes,
         * class specific insertions,
         * subexpressions.";
calculateAmplitudes[diagrams:`type`diagramSet, amplitudes:`type`amplitudeSet] :=
Module[{
      proc = getProcess@amplitudes,
      genericInsertions = getFieldInsertions@diagrams,
      combinatorialFactors = CombinatorialFactorsForClasses /@ List@@amplitudes,
      ampsGen = FeynArts`PickLevel[Generic][amplitudes],
      feynAmps, generic, chains, subs, zeroedRules},
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
      {ampsGen, getRegularizationSettings@#, getMomSettings@#} &@ diagrams,
      "Amplitude calculation"] //. FormCalc`GenericList[];
   generic = MapThread[getGenericSum, {feynAmps, getSumSettings@diagrams}];
   {generic, chains, subs} = proceedChains[diagrams, amplitudes, generic];
   setZeroMassRules@{amplitudes, feynAmps};
   {generic, chains, subs} = makeMassesZero[
      {generic, chains, subs}, diagrams, $ZeroMomenta];
   FCAmplitudesToFSConvention[
      {  generic,
         genericInsertions,
         combinatorialFactors,
         getColourFactors@diagrams},
      chains,
      subs] /. getExternalMomentumRules[$ZeroMomenta, amplitudes]];
calculatedAmplitudes // secure;

setZeroMassRules::usage = "
@brief For a given sets of ``FeynArts`` amd ``FormCalc`` amplitudes creates
       rules to nullify masses of external particles.
@param fa A set of ``FeynArts`` amplitudes.
@param fc A set of ``FormCalc`` amplitudes.
@returns Null.
@note Both of sets are required, because, unfortunately, ``FormCalc``
      introduces new abbreviations, which mix with ``FeynArts`` ones.
@note Amplitudes are taken, because they do not have colour structures already.
@note Explicit names for masses are expected only for external particles.
@returns A list of rules to nullify masses of external particles.
@note Rules of external particle ``#i`` are under numbers
      ``2*#i`` and ``2*#i-1``.";
Module[{rules},
   setZeroMassRules[{fa:`type`amplitudeSet, fc:`type`fc`amplitudeSet}] :=
      rules = RuleDelayed[#, 0] &/@
         Riffle[getExternalMasses@fa, getExternalMasses@fc];
   setZeroMassRules // secure;
   getZeroMassRules[] := (
      Utils`AssertOrQuit[Head@rules =!= Symbol, getZeroMassRules::errNotSet];
      rules);];
getZeroMassRules::errNotSet = "
Call setZeroMassRules to set up rules first.";
getZeroMassRules // secure;

makeMassesZero::usage = "
@brief Sets the masses of external particles to zero everywhere, except loop
       integrals, applies subexpressions.
@param generic An expression to modify.
@param chains A lis with fermionic chains.
@param subs A list of subexpressions.
@param diagrams A set of diagrams.
@returns A list with modified expression, chains and empty non-applied
         subexpressions.";
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
Module[{zeroedRules, new},
   zeroedRules = Cases[subs, Rule[_, pair:FormCalc`Pair[_, _]] :> (pair->0)];
   {new, zeroedRules} = ZeroRules[subs, zeroedRules];
   {  generic /. zeroedRules,
      setZeroExternalMomentaInChains@chains,
      new}];
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
mapThread[func_, exprs:{__}, text:_String:""] :=
Module[{sr, print, percent, init, end, bar, def = 70, dots,
      tot = Length@First@exprs, result},
   subWrite["\n"<>text<>" ...\n"];
   sr[str:_String, num:_Integer] := StringJoin@@Array[str&, num];
   print = (
      percent = #/tot;
      init = "["<>ToString@now<>"/"<>ToString@tot<>"] [";
      end = "] "<>ToString@Floor[100*percent]<>"%";
      bar = def - StringLength[init<>end];
      dots = Floor[bar*percent];
      subWrite@StringJoin[init, sr[".", dots], sr[" ", bar-dots],end,"\r"];)&;
   result = Table[(print@now; func@@(#[[now]]&/@ exprs)), {now, tot}];
   subWrite@"\033[K\033[A";
   subWrite[text<>" ... done\n"];
   result];
mapThread // secure;

CombinatorialFactorsForClasses::usage="
@brief Takes generic amplitude and finds numerical combinatirical factors
       which arise at class level.
@returns List of combinatorical factors for a given generic amplitude.
@param amp ``FeynArts`FeynAmp``.";
CombinatorialFactorsForClasses[amp:FeynArts`FeynAmp[_,_,_,rules_->_[_][classReplacements__]]] :=
{classReplacements}[[ All,#[[1,1]] ]] /.
   {  FeynArts`IndexDelta[___] -> 1,
      FeynArts`SumOver[__] -> 1} &@ Position[rules, FeynArts`RelativeCF];
CombinatorialFactorsForClasses // secure;

getGenericFields::usage = "
@brief Generates a list of unique sorted generic fields in expression.
@param expr An expression, where to search.
@returns A list of unique sorted generic fields.";
getGenericFields[expr:_] :=
   Sort@DeleteDuplicates[Cases[expr, `type`genericField, Infinity]];
getGenericFields // secure;

getGenericSum::usage= "
@brief Converts ``FormCalc`Amp`` into ``NPointFunctions`GenericSum`` object
       using restriction rules for generic fields.
@param amplitude ``FormCalc`Amp`` expression.
@param sumRules A set of rules, restricting the summation.
@returns A ``NPointFunctions`GenericSum`` object.";
getGenericSum[amplitude:`type`fc`amplitude, sumRules:{Rule[_Integer, _]...}] :=
Module[{sort, rules},
   sort = getGenericFields@amplitude;
   rules = Append[sumRules, _Integer -> False];
   GenericSum[
      List@@amplitude,
      sort /. f_[_[_,i_]] :> {f@GenericIndex@i, i /. rules}]];
getGenericSum // secure;

ZeroRules::usage = "
@brief Given a set of rules that map to zero and a set that does
       not map to zero, apply the zero rules to the non-zero ones
       recursively until the non-zero rules do not change anymore.
@param nonzeroRules The list of nonzero rules.
@param zeroRules The list of zero rules.
@returns a list of rules that map the same expressions as the initial rules.";
ZeroRules[nonzeroRules:{Rule[_,_]...}, zeroRules:{Rule[_,0]...}] :=
Module[{newNonzero, newZeroRules},
   newNonzero = Thread[
      Rule[nonzeroRules[[All,1]],nonzeroRules[[All,2]] /. zeroRules]];
   If[newNonzero === nonzeroRules, Return[{nonzeroRules, zeroRules}]];
   newZeroRules = Cases[newNonzero,HoldPattern[_->0]];
   newNonzero = Complement[newNonzero, newZeroRules];
   ZeroRules[newNonzero, Join[zeroRules,newZeroRules]]];
ZeroRules // secure;

FCAmplitudesToFSConvention::usage = "
@brief Translate a list of ``FormCalc`` amplitudes and their abbreviations and
       subexpressions into ``FlexibleSUSY`` language.
@param amplitudes The given list of amplitudes.
@param abbreviations A list of abbreviations.
@param subexpressions A list of subexpressions.
@returns A list of amplitudes and joined abbreviations and subexpressions.";
FCAmplitudesToFSConvention[amplitudes_, abbreviations_, subexpressions_] :=
Module[{fsAmplitudes, fsAbbreviations, fsSubexpressions},
   fsSubexpressions = subexpressions //. $SubexpressionRules;
   fsAmplitudes = amplitudes //. $AmplitudeRules;
   fsAbbreviations = abbreviations //.  $SubexpressionRules //.
      {  FormCalc`Spinor -> SARAH`DiracSpinor,
         FormCalc`Lor -> SARAH`Lorentz};
   {fsAmplitudes, Join[fsAbbreviations,fsSubexpressions]}];
FCAmplitudesToFSConvention // secure;

End[];
EndPackage[];
