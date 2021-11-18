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

BeginPackage["Utils`"];

AppendOrReplaceInList::usage="Replaces existing element in list,
or appends it if not already present.";

ApplyAndConcatenate::usage = "Applies a function to a list and
concatenates the resulting list.";

InputFormOfNonStrings::usage = "apply InputForm to non-strings";

MaxRelDiff::usage="Returns maximum relative difference between numbers
 in a given list.  The second argument defines the numeric underflow.

In[]:= MaxRelDiff[{1, 1.1, 0.9}]
Out[]= 0.181818

In[]:= MaxRelDiff[{0, 10^(-100)}]
Out[]= 0

In[]:= MaxRelDiff[{0, 10^(-100)}, 10^(-101)]
Out[]= 1

In[]:= MaxRelDiff[{1, -1}]
Out[]= 2
";

StringJoinWithSeparator::usage = "Joins a list of strings with a given separator string";

StringJoinWithReplacement::usage =
"Joins a list of strings with a given separator string, making string replacement afterwards";

Zip::usage = "Combines two lists to a list of touples.
Example:

In[]:= Zip[{a,b,c},{d,e,f}]
Out[]= {{a, d}, {b, e}, {c, f}}
";

StringZip::usage = "Combines lists of strings to a list of concated
strings.  Example:

In[]:= StringZip[{\"a\",\"b\"},{\"d\",\"e\"}]
Out[]= {\"ad\", \"be\"}}
";

StringZipWithSeparator::usage = "Combines lists of strings to a list
of concated strings using a separator.  Example:

In[]:= StringZipWithSeparator[{\"a\",\"b\"},{\"d\",\"e\"}, \"_\"]
Out[]= {\"a_d\", \"b_e\"}}
";

SplitList::usage = "split list into list of sub-lists of given maximum
 size";

ForceJoin::usage = "Joins the given arguments if they are lists.";

FSGetOption::usage = "Returns the value of an option from a list of
options.

Example:

In[]:= opts = {o1 -> x, o2 -> y}

In[]:= FSGetOption[opts, o1]
Out[]= x

In[]:= FSGetOption[opts, o2]
Out[]= y

In[]:= FSGetOption[opts, o3]
Error: option o3 not found
";

FSSetOption::usage = "Returns given list of options, where the
occurrence of the given rule is replaced (if it exists) or added (if
it does not exist).

Example:

In[]:= opts = {o1 -> x, o2 -> y}

In[]:= FSSetOption[opts, o1 -> z]
Out[]= {o2 -> y, o1 -> z}

In[]:= FSSetOption[opts, o3 -> z]
Out[]= {o1 -> x, o2 -> y, o3 -> z}
";

FSImportString::usage = "Returns the content of a file in form of a string.  If the file does not exist, \"unknown\" is returned.";

FSStringPadLeft::usage = "StringPadLeft[] for Mathematica 9 and below.";

StartProgressBar::usage = "Starts progress indicator.

Example:

   StartProgressBar[Dynamic[k], 100];
   For[k = 1, k <= 100, k++,
       UpdateProgressBar[k, 100];
       DoSomething[];
      ];
   StopProgressBar[100];
";

UpdateProgressBar::usage = "updates progress indicator.";

StopProgressBar::usage = "stops progress bar.";

FSColor::usage = "Default FlexibleSUSY color";

FSFancyPrint::usage = "Print text in fancy headline style";

FSFancyLine::usage = "Print separator line in command line mode";

FSFancyWarning::usage = "Print a warning with a style."

PrintHeadline::usage = "Print fancy head line";

PrintAndReturn::usage = "Print result and return it";

AssertWithMessage::usage = "AssertWithMessage[assertion_, message_String]:
If assertion does not evaluate to True, print message and Quit[1].";

AssertOrQuit::usage = "
@brief If assertion === True, True is returned.
       If assertion =!= True, an error message is printed and kernel is killed.
@param assertion Any expression.
@param sym::tag A message string.
@param insertions Expressions inserted into the message via StringTemplate.";

MakeUnknownInputDefinition::usage = "
@brief After the definition of any function use it like:

          Utils`MakeUnknownInputDefinition[function];

       It creates a new definition function[args___], which prints an error and
       kills the kernel if used in the code.
@param sym A name of a function.
@note UpValues for symbol[args___] are not cleared.";

DynamicInclude::usage = "
@brief Once gets a file(s), specified by the input string.
@param file A name of a file to load.";

ReadLinesInFile::usage = "ReadLinesInFile[fileName_String]:
Read the entire contents of the file given by fileName and return it
as a list of Strings representing the lines in the file.
Warning: This function may ignore empty lines.";

MathIndexToCPP::usage = "Converts integer-literal index from mathematica to c/c++ convention";
FSPermutationSign::usage = "Returns the sign of a permutation given in a Cycles form";

DumpStart::usage ="
@brief Used to start the model from existing \"kernel snapshot\", i.e. from
       the existing copy of kernel definitions after the first initialization
       of the model via \"normal\" SARAH`Start way.
@param model The name of the model to work with.
@returns Null.";

Begin["`Private`"];

AppendOrReplaceInList[values_List, elem_, test_:SameQ] :=
    Module[{matches, result},
           matches = test[elem, #]& /@ values;
           matches = (If[# =!= True && # =!= False, False, #])& /@ matches;
           If[!Or @@ matches,
              result = Append[values, elem];,
              result = ReplacePart[values, Position[matches, True] -> elem];
             ];
           result
          ];

ApplyAndConcatenate[Func_, l_List] :=
    Module[{result = ""},
           (result = result <> Evaluate[Func[#]])& /@ l;
           result
          ];

ApplyAndConcatenate[Func_, l_] := Evaluate[Func[l]];

SetAttributes[ApplyAndConcatenate, HoldFirst];

StringJoinWithSeparator[list_List, separator_String, transformer_:ToString] :=
    StringJoin[Riffle[transformer /@ list, separator]];

Zip[list1_List, list2_List] :=
    MapThread[List, {list1, list2}];

StringZip[lists___List] :=
    MapThread[StringJoin, {lists}];

StringZipWithSeparator[lists___List, separator_String] :=
    MapThread[StringJoinWithSeparator[{##},separator]&, {lists}];

SplitList[lst_List, 0] := {lst};

SplitList[lst_List, size_Integer] :=
    Module[{result = {}, list = lst, drops},
           While[list =!= {},
                 drops = Min[size,Length[list]];
                 AppendTo[result, Take[list, drops]];
                 list = Drop[list, drops];
                ];
           result
          ];

FSGetOption[opts_List, opt_] :=
    Module[{values},
           values = Cases[opts, (Rule[opt, value_] | RuleDelayed[opt, value_]) :> value];
           Switch[Length[values],
                  0, Print["Error: option ", opt, " not found"];
                     Null,
                  1, values[[1]],
                  _, FSFancyWarning["Option ", opt, " is not unique"];
                     values[[1]]
                 ]
          ];

FSSetOption[opts_List, rule:Rule[opt_, value_]] :=
    Join[
        Cases[opts, Except[Rule[opt,_] | RuleDelayed[opt,_]]],
        {rule}
        ];

FSSetOption[opts_List, rule:RuleDelayed[opt_, value_]] :=
    Join[
        Cases[opts, Except[Rule[opt,_] | RuleDelayed[opt,_]]],
        {rule}
        ];

FSImportString[fileName_String] :=
    Module[{str = Import[fileName, "String"]},
           If[str =!= $Failed,
              str,
              "unknown"
             ]
          ];

FSStringPadLeft[str_String, width_, pad_String] :=
    StringJoin[PadLeft[Characters[str], width, pad]];

ForceJoin[elem___] :=
    Join[Sequence @@ Select[{elem}, (Head[#] === List)&]];

InputFormOfNonStrings[a_String] := a;
InputFormOfNonStrings[a_] := InputForm[a];

MaxRelDiff[{}, _] := 0;

MaxRelDiff[{a_, b_}, underflow_:10^(-16)] :=
    If[Max[Abs[{a,b}]] < underflow,
       0,
       Abs[(a-b)/Max[Abs[{a,b}]]]
      ];

MaxRelDiff[numbers_List, underflow_:10^(-16)] :=
    Max[MaxRelDiff[#,underflow]& /@ Tuples[numbers, 2]];

StartProgressBar[dyn:Dynamic[x_], total_, len_:50] :=
    If[$Notebooks,
       PrintTemporary @ Row[{ProgressIndicator[dyn, {0, total}], " Total: ", total}]
      ];

UpdateProgressBar[x_, total_, len_:50] :=
    Module[{i},
           If[!$Notebooks,
              WriteString["stdout", "[" <> StringJoin[
                  Join[Table[".",{i,Round[len*x/total]}],
                       Table[" ",{i,Round[len*(1-x/total)]}]]
                 ] <> "] " <> ToString[x] <> "/" <> ToString[total] <> "\r"];
             ];
          ];

StopProgressBar[total_, len_:50] :=
    If[!$Notebooks,
       WriteString["stdout", "[" <> StringJoin[
           Table[".",{i,Round[len]}]
       ] <> "] " <> ToString[total] <> "/" <> ToString[total] <> "\n"];
      ];

FSColor = Blue;

FSFancyPrint[text_, level_:1] :=
    Print[Style[text, "Section", FontSize->14 - 2 level, FSColor]]

FSFancyLine[type_:"-", style__:Bold] :=
    If[!$Notebooks, Print[Style[StringJoin[Array[type&, 70]], style]]];

PrintHeadline[text__] :=
    Block[{},
          Print[""];
          FSFancyLine[];
          FSFancyPrint[text];
          FSFancyLine[];
         ];

FSFancyWarning[string_String, len_Integer:70] :=
Module[{warning, chopped},
   warning = If[!$Notebooks,
      If[FlexibleSUSY`FSEnableColors,
         "\033[1;36mWarning\033[1;0m: ",
         "Warning: "
      ],
      Style["Warning: ", Cyan]
   ];
   chopped = InsertLinebreaks[StringReplace[string, "\n"-> " "], len-9];
   chopped = StringReplace[chopped, "\n"-> "\n         "];
   WriteString[$Output, warning <> chopped <> "\n"]
   If[!$Notebooks,
      WriteString[$Output, warning <> chopped <> "\n"];,
      Print[warning, chopped];
   ];
];

PrintAndReturn[e___] := (Print[e]; e)

AssertWithMessage[assertion_, message_String] :=
	If[assertion =!= True, Print[message]; Quit[1]];

AssertOrQuit[assertion_,
             HoldPattern@MessageName[sym_, tag_],
             insertions___] :=
Module[{RedString, WriteOut},
   If[TrueQ@assertion, Return@True];
   RedString[str_] := If[$Notebooks,
      Style[str, Red],
      If[TrueQ@FlexibleSUSY`FSEnableColors,
         "\033[1;31m"<>str<>"\033[1;0m",
         str
      ]
   ];
   WriteOut[str__] := If[$Notebooks,
      Print[str],
      WriteString[$Output, StringJoin[str]<>"\n"]
   ];
   WriteOut[Context@sym, SymbolName@sym, ": ", RedString@tag, ":"];
   WriteOut@StringTemplate[MessageName[sym, tag]][insertions];
   WriteOut["Wolfram Language kernel session ", RedString@"terminated", "."];
   Quit[1];
];

AssertOrQuit[x___] := AssertOrQuit[False, AssertOrQuit::errInput, {x}];
AssertOrQuit ~ SetAttributes ~ {HoldAll, Protected};
AssertOrQuit::errInput = "Input '`1`' is not supported.";

Options[FSFancyWarning] = {
   PageWidth-> 70
};

FSFancyWarning[input__, OptionsPattern[]] :=
Module[{warning, chopped, string},
   string = StringReplace[StringJoin[ToString/@ {input}], "\n"-> " "];
   warning = If[!$Notebooks,
      If[TrueQ@FlexibleSUSY`FSEnableColors,
         "\033[1;36mWarning:\033[1;0m ",
         "Warning: "
      ],
      Style["Warning: ", Cyan]
   ];
   chopped = InsertLinebreaks[string, OptionValue[PageWidth]-9];
   chopped = StringReplace[chopped, "\n"-> "\n         "];
   If[!$Notebooks,
      WriteString[$Output, warning <> chopped <> "\n"];,
      Print[warning, chopped];
   ];
];

MakeUnknownInputDefinition[sym_Symbol] :=
Module[{up, down, all, usageString, infoString, simplify},
   Off[Unset::norep];
      sym[args___] =.;
   On[Unset::norep];

   If[MatchQ[sym::usage, _String],
      usageString = StringJoin["Usage:\n", sym::usage, "\n\n"];,
      usageString = "";
   ];

   simplify[expr_] := StringReplace[ToString@expr,
      StartOfString ~~ "HoldPattern[" ~~ x___ ~~ "]" ~~ EndOfString :> x];

   up = simplify/@ First/@ UpValues@ sym;
   down = simplify/@ First/@ DownValues@ sym;
   all = Join[up, down];

   If[all === {},
      infoString = "",
      infoString = Array[(ToString[#]<>") "<>all[[#]])&, Length@all];
      infoString = StringRiffle[infoString, "\n"];
      infoString = "The behavior for case(s):\n"<>infoString<>
         "\nis defined only.\n\n";
   ];

   With[{name = ToString@sym},
      sym::errUnknownInput = "`1``2`Call\n`3`[`4`]\nis not supported.";
      sym[args___] := AssertOrQuit[False,
         sym::errUnknownInput,
         usageString,
         infoString,
         name,
         StringJoinWithSeparator[{args},", "]
      ];
   ];
];
MakeUnknownInputDefinition@MakeUnknownInputDefinition;
MakeUnknownInputDefinition // Protect;

DynamicInclude[path_String] :=
Module[{separatedQ, wildQ, files},
   separatedQ = StringContainsQ[path, $PathnameSeparator];
   wildQ = StringContainsQ[path, "*"];

   Switch[{wildQ, separatedQ},
      {False, False},
         Return@singleInclude@FileNameJoin@{DirectoryName@$Input, path};,
      {False, True},
         Return@singleInclude@path;,
      {True, False},
         files = FileNames@FileNameJoin@{DirectoryName@$Input, path};
         Return[singleInclude/@files];,
      {True, True},
         Return[singleInclude/@FileNames@path];
   ];
];
DynamicInclude // MakeUnknownInputDefinition;
DynamicInclude // Protect;

singleInclude[path_String] :=
With[{inserted = path},
   AssertOrQuit[And[FileExistsQ@path, Not@DirectoryQ@path],
      DynamicInclude::errNoFile,
      path
   ];
   Once@Get@inserted;
   Return@inserted;
];
singleInclude // MakeUnknownInputDefinition;
singleInclude // Protect;

DynamicInclude::errNoFile = "\nFile\n`1`\ndoes not exist and can't be loaded!";

abbreviateLongString[expr_] :=
Module[{str, b},
   b[s_] := If[!$Notebooks,"\033[1;36m"<>s<>"\033[0m", s];
   str = ToString@expr;
   If[StringLength@str > 1000, b["<"]<>StringTake[str, 1000]<>b[">"], str]];
abbreviateLongString // MakeUnknownInputDefinition;
abbreviateLongString // Protect;

StringJoinWithReplacement[
   list_List,
   separator:_String:", ",
   replacement:Rule[_String,_String]:Rule["`","`.`"],
   transformer_:ToString
] :=
StringReplace[StringJoinWithSeparator[list,separator,transformer],replacement];
StringJoinWithReplacement // MakeUnknownInputDefinition;
StringJoinWithReplacement ~ SetAttributes ~ {Locked,Protected};

ReadLinesInFile[fileName_String] :=
	Module[{fileHandle, lines = {}, line},
		fileHandle = OpenRead[fileName, BinaryFormat -> True];

		While[(line = Read[fileHandle, String]) =!= EndOfFile,
			AssertWithMessage[line =!= $Failed,
				"Utils`ReadLinesInFile[]: Unable to read line from file '" <>
				fileName <> "'"];
			AppendTo[lines, line];
			];

    Close[fileHandle];
    lines
	]

(* MathIndexToCPP *)

MathIndexToCPP[i_Integer /; i>0] := i-1;

MathIndexToCPP::wrongInt =
"Cannot convert index of value \"`1`\". Index value cannot be smaller than \"1\".";
MathIndexToCPP[i_Integer] := AssertOrQuit[False, MathIndexToCPP::wrongInt, StringJoin@@Riffle[ToString/@{i},", "]];

MathIndexToCPP::nonIntInput =
"Cannot convert a non integer index \"`1`\".";
MathIndexToCPP[i___] := AssertOrQuit[False, MathIndexToCPP::nonIntInput, StringJoin@@Riffle[ToString/@{i},", "]];

(* FSPermutationSign *)

(* from https://reference.wolfram.com/language/tutorial/Permutations.html *)
FSPermutationSign[perm_?PermutationCyclesQ] :=
    Apply[Times, (-1)^(Length /@ First[perm] - 1)];
FSPermutationSign[perm___] :=
    (Print[perm, " is not a permutation in disjoint cyclic form."];Quit[1]);

End[];

EndPackage[];
