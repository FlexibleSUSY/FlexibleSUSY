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

Needs["Utils`",
   FileNameJoin@{ParentDirectory@DirectoryName@$InputFileName, "Utils.m"}];

Block[{Format},
   Needs@"FeynArts`";
   Needs@"FormCalc`";
   Print[];];

BeginPackage@"NPointFunctions`";

(* Reserving names *)
Off[General::shdw];
{  DiracChain, Mat};
On[General::shdw];

{  NPointFunction};
{  LorentzIndex, GenericSum, GenericIndex,
   GenericS, GenericF, GenericV, GenericU,
   OperatorsOnly, ExceptLoops} ~ SetAttributes ~ {Protected};

Begin@"`Private`";

secure[sym:_Symbol] :=
   Protect@Evaluate@Utils`MakeUnknownInputDefinition@sym;
secure // secure;

With[{dir = DirectoryName@$InputFileName},
   Get@FileNameJoin@{dir, #<>".m"}&/@
      {"type", "rules", "settings", "chains", "topologies", "tree"};];

NPointFunction::usage = "
@brief The entry point of the calculation.
       Sets files and options up, calculates diagrams and amplitudes using
       ``FeynArts`` and ``FormCalc``.
@param formcalc A ``String`` directory of ``FormCalc`` output.
@param model A ``String`` name of ``FeynArts`` model file.
@param particles A ``String`` name of a file, which is created by ``SARAH``
       and contains ``FeynArts`` particle names.
@param contexts A ``String`` name of a file, which contains information
       about particle contexts in ``SARAH`` conventions.
@param in An expression, containing ``FeynArts`` fields (each as ``String``),
       required to be *incoming* fields for the amplitude.
@param out An expression, containing ``FeynArts`` fields (each as ``String``),
       required to be *outgoing* fields for the amplitude.
@param observable An observable name in the form
       ``FlexibleSUSYObservable`<Name>[<Arguments>]``.
@param loops The loop level of calculation (as ``Integer``).
@param processes A ``List`` of one or many ``Symbol`` of processes
       to calculate.
@param momenta A ``Symbol``, which defines how to treat momenta of external
       particles.
@param onShell A flag, corresponding to the \"on-shellness\" of external
       particles.
@param scheme The main regularization scheme to be used.
       Can be overriden by ```settings`regularization``.
@returns An object of the n-point function in ``FlexibleSUSY`` conventions."
NPointFunction[
   {formcalc_, model_, particles_, contexts_, in_, out_},
   {observable_, loops_, processes_, momenta_, onShell_, scheme_}] :=
   Module[{tree},
      BeginPackage["NPointFunction`"];
      Begin["`Private`"];
      `file`particles[] := particles;
      `file`contexts[] := contexts;
      `options`observable[] := SymbolName@Head@observable;
      `options`loops[] := loops;
      `options`processes[] := processes;
      `options`momenta[] := momenta;
      `options`onShell[] := onShell;
      `options`scheme[] :=
         Switch[scheme, FlexibleSUSY`DRbar, 4, FlexibleSUSY`MSbar, D];
      End[];
      EndPackage[];

      FeynArts`$FAVerbose = 0;
      FeynArts`InitializeModel@model;
      SetOptions[FeynArts`InsertFields,
         FeynArts`Model -> model,
         FeynArts`InsertionLevel -> FeynArts`Classes];
      FormCalc`$FCVerbose = 0;
      If[!DirectoryQ@formcalc, CreateDirectory@formcalc];
      SetDirectory@formcalc;

      settings[];
      tree = cut[plant[in, out], `settings`diagrams];
      tree = cut[plant@tree, `settings`amplitudes];
      picture@tree;
      {`rules`fields@fields@tree, calculateAmplitudes@tree}];
NPointFunction // secure;

genericIndex[index:_Integer] := FeynArts`Index[Generic, index];
genericIndex // secure;

genericMass[field:`type`field, index:_Integer] :=
   FeynArts`Mass[field@genericIndex@index, `type`mass];
genericMass[field:`type`field] :=
   FeynArts`Mass[field@`type`genericIndex, `type`mass];
genericMass // secure;

process[set:`type`diagramSet|`type`amplitudeSet] :=
   Cases[Head@set, (FeynArts`Process -> e:_) :> e][[1]];
process[set:`type`fc`amplitudeSet] :=
   Part[Head@Part[set, 1], 1];
process // secure;

externalMasses[set:`type`fc`amplitudeSet] :=
   Flatten[List@@process@set, 1][[All, 3]];
externalMasses[tree:`type`tree] :=
   FeynArts`Mass[# /. -1 -> 1] &/@ fields[tree, Flatten];
externalMasses // secure;

getField[set:`type`diagramSet, i:_Integer] :=
   Flatten[List@@process@set, 1][[i]] /; 0<i<=Plus@@(Length/@process@set);
getField // secure;

Field[d:Head@`type`diagramSet, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]];
Field // secure;

FieldPattern[d:Head@`type`diagramSet, i_Integer] :=
   Flatten[List@@(FeynArts`Process /. List@@d), 1][[i]] /.
      `type`generationIndex :> Blank[];
FieldPattern[d:Head@`type`diagramSet, a:HoldPattern@Alternatives@__] :=
   FieldPattern[d, #] &/@ a;
FieldPattern // secure;

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
fieldInsertions[tree:`type`tree] :=
   Map[Last, #, {3}] &@ Flatten[ fieldInsertions /@ List@@diagrams@tree, 1];
fieldInsertions[diag:`type`diagram, keepNumQ:True|False:False] :=
   fieldInsertions[#, keepNumQ] &/@ Apply[List, Last@diag, {0, 1}];
fieldInsertions[{graph_, insert_}, keepNumQ_] :=
Module[{toGenericIndexConventionRules, fieldsGen, genericInsertions},
   toGenericIndexConventionRules = Cases[graph,
      Rule[FeynArts`Field[index_Integer],type_Symbol] :>
      Rule[FeynArts`Field@index, type[FeynArts`Index[Generic,index]]]];
   fieldsGen = toGenericIndexConventionRules[[All,1]];
   genericInsertions = Cases[#,
      Rule[genericField_,classesField_] /; MemberQ[fieldsGen, genericField] :>
      Rule[genericField, stripIndices@classesField]] &/@ insert;
   SortBy[#,First]&/@ If[keepNumQ,
      List @@ genericInsertions,
      List @@ genericInsertions /. toGenericIndexConventionRules]];
fieldInsertions // secure;

stripIndices::usage = "
@brief Removes particle indices from a given field.
@param field the given field.
@returns The given field with all indices removed.";
stripIndices[Times[-1, field_]] :=
   -stripIndices@field;
stripIndices[name_[class_, ___]] :=
   name@class;
stripIndices // secure;

calculateAmplitudes::usage = "
@brief Applies ``FormCalc`` routines to amplitude set, simplifies the result.
@param tree A set of data in the form of a ``tree`` object.
@returns The main part of n-point function object, containing:

         * generic amplitudes,
         * class specific insertions,
         * subexpressions.";
calculateAmplitudes[tree:`type`tree] :=
Module[{
      proc = process@amplitudes@tree,
      ampsGen = FeynArts`PickLevel[Generic][amplitudes@tree],
      feynAmps, generic, chains, subs, zeroedRules},
   If[`options`momenta[],
      ampsGen = FormCalc`OffShell[ampsGen,
         Sequence@@Array[#->0&, Plus@@Length/@proc]]];
   feynAmps = mapThread[
      FormCalc`CalcFeynAmp[Head[ampsGen][#1],
         FormCalc`Dimension -> #2,
         FormCalc`OnShell -> `options`onShell[],
         FormCalc`FermionChains -> FormCalc`Chiral,
         FormCalc`FermionOrder -> fermionOrder@tree,
         FormCalc`Invariants -> False,
         FormCalc`MomElim -> #3]&,
      {  ampsGen,
         settings[tree, `settings`regularization, `options`scheme[]],
         settings[tree, `settings`momenta, Automatic]},
      "Amplitude calculation"] //. FormCalc`GenericList[];
   generic = MapThread[getGenericSum,
      {feynAmps, settings[tree, `settings`sum, {}]}];
   {generic, chains, subs} = proceedChains[tree, generic];
   setZeroMassRules@{tree, feynAmps};
   {generic, chains, subs} = makeMassesZero[
      {generic, chains, subs}, tree, `options`momenta[]];
   convertToFS[
      {  generic,
         fieldInsertions@tree,
         combinatoricalFactors@tree,
         colorFactors@tree},
      chains,
      subs] /. `rules`externalMomenta[tree, `options`momenta[]]];
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
   setZeroMassRules[{tree:`type`tree, fc:`type`fc`amplitudeSet}] :=
      rules = RuleDelayed[#, 0] &/@
         Riffle[externalMasses@tree, externalMasses@fc];
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
   {generic_, chains_, subs_}, tree:`type`tree, ExceptLoops] :=
Module[{funcs, names, pattern, uniqueIntegrals, hideInt, showInt, rules, new},
   funcs = settings[tree, `settings`massless, {}, Identity];
   names = ToExpression/@Names@RegularExpression@"LoopTools`[ABCD]\\d+i*";
   subWrite@"Applying subexpressions ... ";
   new = generic //. subs;
   subWrite@"done\n";
   pattern = Alternatives @@ ( #[__] &/@ names );
   uniqueIntegrals = DeleteDuplicates@Cases[new, pattern, Infinity];
   hideInt = Rule[#, Unique@"loopIntegral"] &/@ uniqueIntegrals;
   showInt = hideInt /. Rule[x_, y_] -> Rule[y, x];
   rules = List@@MapIndexed[Composition[Sequence@@funcs[[#2[[1]]]]]@
      getZeroMassRules[]&, new];
   {  List@@MapIndexed[#1 //. rules[[#2[[1]]]] /. showInt&, new /. hideInt /.
         FormCalc`Pair[_,_] -> 0],
      zeroMomenta@chains /. getZeroMassRules[],
      {}}];
makeMassesZero[{generic_, chains_, subs_}, _, True] :=
Module[{zeroedRules, new},
   zeroedRules = Cases[subs, Rule[_, pair:FormCalc`Pair[_, _]] :> (pair->0)];
   {new, zeroedRules} = ZeroRules[subs, zeroedRules];
   {  generic /. zeroedRules,
      zeroMomenta@chains,
      new}];
makeMassesZero[e:{_, _, _}, __] := e;
makeMassesZero // secure;

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
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "NPointFunctions`"];
