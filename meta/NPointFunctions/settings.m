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
Begin@"`Private`";

settings::usage = "
@brief Defines default settings. Loads all settings from file, corresponding
       to an observable. Can be used to parse the following settings:

       ```settings`regularization``
          Some amplitudes are calculated incorrectly in CDR.
          For handling this, one can override used scheme for some topologies.
          Uses rules one-by-one for replacing if allowed by topology.
       ```settings`momenta``
          Eliminate specific momenta in some topologies.
          Uses rules one-by-one for replacing if allowed by topology.
       ```settings`sum``
          Skip summation over some indices of particles.
       ```settings`massless``
          Do not put some masses to zero.
@param tree An object to work with.
@param settings A list of data, which specifies replacements for each topology.
@param default A default value to be used after the ones, given by
       ``settings``.
@param head Several settings can act on the same topology. All results are
       wrapped by this function. Default is ``First``.
@returns A set of settings for ``FormCalc`Dimension``.";
With[{dir = DirectoryName@$InputFileName},
   settings[] :=
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
            {dir, `options`observable[], "settings.m"};
         Protect@Evaluate[Context[]<>"settings`*"];
         End[];
         EndPackage[];);];
settings[tree:type`tree,
   settings:{__}|Default, default_, head:_:First] :=
   Module[{res = {tree}},
      If[settings =!= Default,
         AppendTo[res, applySetting[tree, #]]&/@ settings];
      res = res /.
         node[type`generic, __] -> default /.
         node[type`topology, rest__] :> rest /.
         node[type`head, rest__] :> {rest};
      DeleteDuplicates/@Transpose@res /.
         {default, rest__} :> head@{rest} /. {default} -> default];
settings // secure;

parseSettings::usage = "
@brief Converts ```settings`diagrams`` or ```settings`amplitudes`` to a set
       of operations, which define what to remove from the ``tree`` object.
       The structure of these settings is the following::

          {  anySymbol -> {  {  inProcesses...},
                             {  notInProcesses...}}..}

       anySymbol
          Any ``Symbol``, which defines the action to do with the ``tree``.
       inProcesses
          Is applied, if ``anySymbol`` is inside ```options`processes[]``.
          Has the syntax, which obeys the pattern::

             {  str, tQ, fun}
       notInProcesses
          Is applied, if ``anySymbol`` is **not** inside
          ```options`processes[]``. See the syntax above.
       str
          Any ``String``, which is printed, if the action is applied.
       tQ
          A criterion, which decides, whether to apply action or not.
          It is applied if ``tQ[topology]`` returns ``True``.
       fun
          A function, which is used to remove generic or classes levels.
          It is applied to 3 arguments: 1)current expression (any node at
          classes or generic level), 2) topology, 3) and the head of diagram
          set. If it returns ``True``, then the node is kept.
@param settings An entry to define actions.
@returns A set of settings.";
parseSettings[settings:Default] := {};
parseSettings[settings:{Rule[_, {{___}, {___}}]..}] :=
Module[{positiveRules, negativeRules, discardProcesses, clean, parsed},
   parsed = Rule[SymbolName@First@#, Last@#]&/@settings;
   positiveRules = parsed /. Rule[s:_, {p:_, _}] :> Rule[s, p];
   negativeRules = parsed /. Rule[s:_, {_, n:_}] :> Rule[s, n];
   discardProcesses = Complement[parsed[[All, 1]], `options`processes[]];
   clean = RuleDelayed[#, Sequence[]] &/@ `options`processes[];
   DeleteDuplicates@Join[
      DeleteDuplicates@Flatten[`options`processes[] /. positiveRules, 1],
      DeleteDuplicates@Flatten[discardProcesses /. negativeRules, 1]] /.
         clean];
parseSettings // secure;

applySetting[tree:type`tree, {str_String, tQ_, fun_}] :=
   info[cut[tree, tQ, fun], str];

makeApply[pattern_, function:_Symbol] :=
   (  Off@RuleDelayed::rhs;
      applySetting[tree:type`tree, pattern] :=
         Module[{once},
            once[arg_] := once@arg =
               If[# =!= {}, Print@@#]&@
                  Cases[pattern, _String, Infinity, Heads -> True];
            tree /. node[t:type`topology /; tQ@t, rest__] :>
               (once@_; node[t, rest] /. node[g:type`generic, __] :>
                  function[fun, g, t, head@tree])];
      On@RuleDelayed::rhs;);

makeApply[tQ_ -> fun_, value];
makeApply[{str_String, tQ_, fun:{_Integer, _}}, restrict];
makeApply[{str_String, tQ_, fun:{Append, _}}, append];
makeApply[{str_String, tQ_, fun:{Hold, _}}, hold];
applySetting // secure;

value[val_, ___] := val;
restrict[{int_, fun_}, __, head_] :=
   {int -> Or[fun[_, _, head], -fun[_, _, head]]};
append[{Append, (f:type`field)[n_Integer] :> e_Integer}, ___] :=
   With[{rhs = (First /@ getZeroMassRules[])[[2*e-1]]},
      Append[#, genericMass[f, n] :> rhs]&];
hold[{Hold, e_}, ___] :=
   With[{pos = {{2*e}, {2*e-1}}}, Delete[#, pos]&];

End[];
EndPackage[];
