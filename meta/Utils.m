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

PrintHeadline::usage = "Print fancy head line";

PrintAndReturn::usage = "Print result and return it";

AssertWithMessage::usage = "AssertWithMessage[assertion_, message_String]:
If assertion does not evaluate to True, print message and Quit[1].";

AssertOrQuit::usage =
"@brief AssertOrQuit[assertion, sym::tag, insertions...]:
If assertion evaluate to True, returns True.
If assertion evaluate to False, print message with sequence of insertions
and Quit[1].
@param assertion Some expression which one want to check being True or something
else.
@param sym::tag Controlling MessageName String. It is assumed that symbol \"\"
appears only as controllling one, i.e. in the form \"1\", \"25\" etc.
@param insertions... None or more expressions which are inserted inside sym::tag
controlling String. Its length should be equal or more then the maximal
controlling number in sym::tag.
@note First check controlling sym::tag.
@note __~~\"Private\"~~str__:>str is done because Mathematica prints the names
of Private` variables in a quite weird way.
@note \"\\n->\"dummy_new_line\" because Mathematica's StringForm incorrectly
parses \\n symbol.
@note hard-coded width of output colourless text is 70.";

EvaluateOrQuit::usage =
"@brief EvaluateOrQuit[expression, sym::tag, insertions...]:
Evaluates expression and returns the result if Wolfram messages aren't
generated.
If any Message is generated it
   stops evaluation,
   prints message with sequence of insertions,
   prints content of Message and
   Quit[1].
@param expression Some expression which one want to evaluate without Wolfram
messages.
@param sym::tag Controlling MessageName String. It is assumed that symbol \"`\"
appears only as controllling one, i.e. in the form \"`1`\", \"`25`\" etc.
@param insertions... None or more expressions which are inserted inside sym::tag
controlling String. Its length should be equal or more then the maximal
controlling number in sym::tag.
@note First check controlling sym::tag.
@note hard-coded width of output colourless text is 70.
@note It seems that function is not working for messages generated in subkernels
if they are running inside expression.
@note __~~\"Private`\"~~str__:>str is done because Mathematica prints the names
of Private` variables in a quite weird way.
@note \"\\n->\"dummy_new_line\" because Mathematica's StringForm incorrectly
parses \\n symbol.
@note Internal`HandlerBlock is not documented, this is why System`Dump` prefix
is used for function arguments to avoid any unexpected behavior.";

MakeUnknownInputDefinition::usage =
"@brief Creates definition for a given symbol for a case when input is not defined
explicitly, i.e. creates definition for the pattern symbol[args___].
@example Step 1) Make all desired definitions for a function (here: foo), like
foo[a_Integer] := Module[{<...>},<...>];
foo[c:{__Integer}] := Module[{<...>},<...>];
Step 2) Secure the usage of foo for specified above cases simply by writing
foo // Utils`MakeUnknownInputDefinition;
or
Utils`MakeUnknownInputDefinition[foo];
or
Utils`MakeUnknownInputDefinition@foo;
somewhere in the scope of the package, where foo is defined.
@param <Symbol> sym Symbol to make definition for.
@returns None.
@note UpValues for symbol[args___] are not cleared.";

ReadLinesInFile::usage = "ReadLinesInFile[fileName_String]:
Read the entire contents of the file given by fileName and return it
as a list of Strings representing the lines in the file.
Warning: This function may ignore empty lines.";

FSReIm::usage = "FS replacement for the mathematica's function ReIm";
FSBooleanQ::usage = "FS replacement for the mathematica's function BooleanQ";
MathIndexToCPP::usage = "Converts integer-literal index from mathematica to c/c++ convention";

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
                  _, Print["Warning: option ", opt, " not unique"];
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

PrintAndReturn[e___] := (Print[e]; e)

AssertWithMessage[assertion_, message_String] :=
	If[assertion =!= True, Print[message]; Quit[1]];

AssertOrQuit::errNotDefined =
"Error message \"`1`\" is not defined in the code.";
AssertOrQuit::errStrokes =
"Even number of `.` symbols in sym::tag should be given in:
\"`1`\"";
AssertOrQuit::errControl =
"Only control symbols \"`.`int Number`.`\" and \"`.`.`.`\" are allowed in:
\"`1`\"";
AssertOrQuit::errInsertions =
"The length of insertions
`1`
should large or equal to the max control number `2` in:
\"`3`\"";
AssertOrQuit::errInput =
"Input should be of the following form:
AssertOrQuit[assertion, sym::tag, insertions...] and not
AssertOrQuit@@`1`.

Read AssertOrQuit::usage for more information.";
AssertOrQuit[assertion_,HoldPattern@MessageName[sym_, tag_],insertions___] :=
   internalAssertOrQuit[assertion,MessageName[sym,tag],insertions] /;
   internalOrQuitInputCheck[AssertOrQuit,MessageName[sym,tag],insertions];
AssertOrQuit[x___] :=
   AssertOrQuit[False,AssertOrQuit::errInput,{x}];
If[!$Notebooks,
   internalAssertOrQuit[assertion_,HoldPattern@MessageName[sym_, tag_],insertions___] :=
   Module[{RedString,WriteOut,MultilineToDummy,replacedMessage},
      If[assertion === True,Return@True];

      RedString[str_] := "\033[1;31m"<>str<>"\033[1;0m";
      WriteOut[str__] := WriteString["stdout"~OutputStream~1,StringJoin@str];
      MultilineToDummy[args___] := Sequence@@(StringReplace[ToString@#,"\n"->"dummy_n"]&/@{args});
      replacedMessage = StringReplace[sym~MessageName~tag,"\n"->"dummy_n"];

      Utils`FSFancyLine[];
      WriteOut[Context@sym,StringReplace[ToString@sym,__~~"`"~~str__:>str],": ",RedString@tag,":\n"];
      WriteOut@StringReplace[ToString@StringForm[replacedMessage,MultilineToDummy@insertions],"dummy_n"->"\n"];
      WriteOut["\nWolfram Language kernel session ",RedString@"terminated",".\n"];
      Utils`FSFancyLine[];
      Quit[1];
   ];,
   (* Else *)
   internalAssertOrQuit[assertion_,HoldPattern@MessageName[sym_, tag_],insertions___] :=
   Module[{WriteColourless,MultilineToDummy,replacedMessage},
      If[assertion === True,Return@True];

      MultilineToDummy[args___] := Sequence@@(StringReplace[ToString@#,"\n"->"dummy_n"]&/@{args});
      replacedMessage = StringReplace[sym~MessageName~tag,"\n"->"dummy_n"];

      Print[Context@sym,StringReplace[ToString@sym,__~~"`"~~str__:>str],": ",Style[tag,Red],":\n",
         StringReplace[ToString@StringForm[replacedMessage,MultilineToDummy@insertions],"dummy_n"->"\n"],
         "\nWolfram Language kernel session ","terminated"~Style~Red,"."];
      Quit[1];
   ];
];
SetAttributes[{AssertOrQuit,internalAssertOrQuit},{HoldAll,Locked,Protected}];

EvaluateOrQuit::errNotDefined = AssertOrQuit::errNotDefined;
EvaluateOrQuit::errStrokes = AssertOrQuit::errStrokes;
EvaluateOrQuit::errControl = AssertOrQuit::errControl;
EvaluateOrQuit::errInsertions = AssertOrQuit::errInsertions;
EvaluateOrQuit::errInput =
"Input should be of the following form:
EvaluateOrQuit[expression, sym::tag, insertions...] and not
EvaluateOrQuit@@`1`.
Read EvaluateOrQuit::usage for more information.";
EvaluateOrQuit[expression_,HoldPattern@MessageName[sym_, tag_],insertions___] :=
   internalEvaluateOrQuit[expression,MessageName[sym,tag],insertions] /;
   internalOrQuitInputCheck[EvaluateOrQuit,MessageName[sym,tag],insertions];
EvaluateOrQuit[x___] :=
   AssertOrQuit[False,EvaluateOrQuit::errInput,{x}];
internalEvaluateOrQuit[
   expression_,
   HoldPattern@MessageName[sym_, tag_],
   insertions___
] :=
Module[
   {
      ctrlRed=If[!$Notebooks,"\033[1;31m",""],
      ctrlBack=If[!$Notebooks,"\033[1;0m",""],
      CutString=If[(!$Notebooks)&&MemberQ[$Packages,"TextFormatting`"],
         TextFormatting`WrapLines[#,70,""]&,#&],
      WriteOut,WriteColourless,Filter
   },
   WriteOut[string__] := WriteString[OutputStream["stdout",1],StringJoin@string];
   WriteColourless[string__] := WriteOut@CutString@StringJoin@string;
   Filter[
      System`Dump`str_,
      Hold[MessageName[System`Dump`s_, System`Dump`t_]],
      Hold[Message[_, System`Dump`args___]]
   ] :=
   (
      Utils`FSFancyLine[];
      WriteOut[Context@sym,StringReplace[ToString@sym,__~~"`"~~str__:>str],
         ": ",ctrlRed,tag,ctrlBack,":\n"];
      WriteColourless[#,"\n"]&/@StringSplit[ToString@StringForm[
         StringReplace[MessageName[sym, tag],"\n"->"dummy_n"],insertions],
         "dummy_n"];
      WriteColourless[ToString@System`Dump`s,"::",System`Dump`t," ",
         ToString@StringForm[System`Dump`str,System`Dump`args]];
      WriteOut["\nWolfram Language kernel session ",ctrlRed,"terminated",ctrlBack,".\n"];
      Utils`FSFancyLine[];
      Quit[1]
   );
   Internal`HandlerBlock[{"MessageTextFilter", Filter}, expression]
];
SetAttributes[{EvaluateOrQuit,internalEvaluateOrQuit},{HoldAll,Locked,Protected}];

internalOrQuitInputCheck[func_,message_,insertions___] :=
Module[{nStrokes,controlSubstrings},
   internalAssertOrQuit[StringQ@message,
      func::errNotDefined,message];
   nStrokes = StringCount[message,"`"];
   internalAssertOrQuit[EvenQ@nStrokes,
      func::errStrokes,message];

   If[nStrokes===0,Return@True];

   controlSubstrings=DeleteDuplicates@StringCases[message,{
      "`.`":>0,(* Ok *)
      "`"~~num:DigitCharacter..~~"`":>FromDigits@num,(* Ok *)
      "`"~~___~~"`":>-1(* Something bad *)
      }];
   internalAssertOrQuit[FreeQ[controlSubstrings,-1],
      func::errControl,message];
   internalAssertOrQuit[TrueQ[Max@controlSubstrings<=Length@{insertions}],
      func::errInsertions,{insertions},Max@checkedControl,message]
];
SetAttributes[internalOrQuitInputCheck,{HoldFirst,Locked,Protected}];

MakeUnknownInputDefinition[sym_Symbol] :=
Module[{usageString,info,parsedInfo,infoString,symbolAsString},
   (* Clean existing definitions if they exist for required pattern.. *)
   Off[Unset::norep];
   sym[args___] =.;
   On[Unset::norep];
   (* Maybe some useful definitions already exist*)
   If[MatchQ[sym::usage,_String],usageString="Usage:\n"<>sym::usage<>"\n\n",usageString=""];
   info = MakeBoxes@Definition@sym;
   If[MatchQ[info,InterpretationBox["Null",__]],(* True - No, there is no definitions. *)
      infoString="",
      parsedInfo = Flatten@# &/@ (Cases[info[[1,1]],GridBox[{x:{_}..},__]:>Cases[{x},{_RowBox},1],2]~Flatten~2 //. {RowBox[x_]:>x,StyleBox[x_,_]:>x});
      parsedInfo = MapThread[parsedInfo[[##]]&,{Range@Length@#,First/@#}] &@ (Range@(First@#-1) &@ Position[#,"="|":="|"^="|"^:="] &/@ parsedInfo);
      parsedInfo = DeleteCases[DeleteDuplicates@parsedInfo,{"Options",__}|{"Attributes",__}];
      parsedInfo = Array[Join[{ToString@#,") "},parsedInfo[[#]]]&,Length@parsedInfo];
      infoString = StringJoin@Riffle[StringJoin @@ # & /@ parsedInfo, "\n"];
      infoString = "The behavior for case"<>If[Length@parsedInfo===1,"\n","s\n"]<>infoString<>"\nis defined only.\n\n";
   ];
   symbolAsString=StringReplace[ToString@sym,"`"->"`.`"];
   sym::errUnknownInput = "`1``2`Call\n"<>symbolAsString<>"[`3`]\nis not supported.";
   (* Define a new pattern. *)
   sym[args___] := AssertOrQuit[False,sym::errUnknownInput,usageString,infoString,StringJoinWithSeparator[{args},", "]];
];
MakeUnknownInputDefinition@MakeUnknownInputDefinition;
SetAttributes[MakeUnknownInputDefinition,{Locked,Protected}];

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

FSReIm[z_] := If[$VersionNumber >= 10.1,
   ReIm[z],
   {Re[z], Im[z]}
];

FSBooleanQ[b_] :=
   If[$VersionNumber >= 10.0,
      BooleanQ[b],
      If[b === True || b === False, True, False]
   ];

(* MathIndexToCPP *)

MathIndexToCPP[i_Integer /; i>0] := i-1;

MathIndexToCPP::wrongInt =
"Cannot convert index of value \"`1`\". Index value cannot be smaller than \"1\".";
MathIndexToCPP[i_Integer] := AssertOrQuit[False, MathIndexToCPP::wrongInt, StringJoin@@Riffle[ToString/@{i},", "]];

MathIndexToCPP::nonIntInput =
"Cannot convert a non integer index \"`1`\".";
MathIndexToCPP[i___] := AssertOrQuit[False, MathIndexToCPP::nonIntInput, StringJoin@@Riffle[ToString/@{i},", "]];

End[];

EndPackage[];
