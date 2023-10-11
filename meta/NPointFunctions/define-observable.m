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

Begin@"Observables`Private`";

Options@Observables`DefineObservable = {
   InsertionFunction -> CConversion`ToValidCSymbolString,
   GetObservableName -> Unset,
   GetObservableDescription -> Unset,
   GetObservableType -> Unset,
   CalculateObservable -> Unset,
   Context -> Unset
};

Observables`DefineObservable[obs_@pattern___, OptionsPattern[]] :=
Module[{stringPattern, patternNames, uniqueNames, lhsRepl, rhsRepl, warn},
   warn := Utils`FSFancyWarning[#," for ", ToString@obs, " might not be specified."]&;

   stringPattern = ToString@FullForm@{pattern};
   patternNames = StringCases[stringPattern, "Pattern[" ~~ Shortest@x__ ~~ "," :> x];
   uniqueNames = Unique[#<>"$"]&/@patternNames;

   lhsRepl = MapThread[Rule, {patternNames, ToString/@uniqueNames}];
   rhsRepl = MapThread[
      RuleDelayed[#1, OptionValue[InsertionFunction]@#2]&,
      {patternNames, uniqueNames}
   ];

   If[OptionValue@Context === Unset,
      warn@"Context";,
      AppendTo[rhsRepl, "context" :> OptionValue[InsertionFunction][
         FlexibleSUSY`FSModelName <> "_" <> OptionValue@Context]];
   ];

   With[{args = Sequence@@ToExpression@StringReplace[stringPattern, lhsRepl],
         repl = rhsRepl,
         name = OptionValue@GetObservableName,
         description = OptionValue@GetObservableDescription,
         type = OptionValue@GetObservableType,
         calculate = OptionValue@CalculateObservable
      },
      AppendTo[FlexibleSUSYObservable`FSObservables, obs];
      If[name === Unset,
         warn@"GetObservableName";,
         Observables`GetObservableName@obs@args := StringReplace[name, repl];
      ];
      If[description === Unset,
         warn@"GetObservableDescription";,
         Observables`GetObservableDescription@obs@args := StringReplace[description, repl];
      ];
      If[type === Unset,
         warn@"GetObservableType";,
         Observables`GetObservableType@obs@args := type;
      ];
      If[calculate === Unset,
         warn@"CalculateObservable";,
         Observables`CalculateObservable[obs@args, structName:_String] :=
            structName <> "." <> StringReplace[name, repl] <>
            StringReplace[" = context::"<>calculate<>";", repl];
      ];
   ];
];
Utils`MakeUnknownInputDefinition@Observables`DefineObservable;

End[];
