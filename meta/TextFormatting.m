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

$flexiblesusyCSrcChunkSize = 2^21;

BeginPackage["TextFormatting`"];

WrapLines::usage="Breaks long text lines.  This function is a wrapper around WrapText[].";
WrapText::usage="WrapText[text, maxWidth, indentation] breaks text lines.  It tries to wrap a line at a blank or a special character that maximizes the line length within maxWidth characters.";
IndentText::usage="indents text by a given number of spaces";

Begin["`Private`"];

WrapLines[text_String, maxWidth_:79, offset_:"   "] :=
    WrapText[text, maxWidth, StringLength[offset]];

IndentText[text_String, spaces_Integer:3] :=
    Module[{i, whiteSpace = ""},
           If[text == "", Return[text];];
           For[i = 0, i < spaces, i++, whiteSpace = whiteSpace <> " "];
           Return[StringReplace[whiteSpace <> text,
                                { x:("\n" ~~ ("\n"..)) :> x <> whiteSpace,
                                  "\n" ~~ EndOfString -> "\n",
                                  "\n" -> "\n" <> whiteSpace } ]];
          ];

global$IterationLimit = $IterationLimit;

WrapText[text_String, maxWidth_Integer:79, indentation_Integer:2] := Block[{
    maxLength = maxWidth,
    relSkip = indentation,
    $IterationLimit = Max[Global`$flexiblesusyCSrcChunkSize / 32,
			  global$IterationLimit]
  },
  StringJoin @ Riffle[WrapLine /@ StringSplit[text, "\n", All], "\n"]
];

RemoveTrailingWhitespace[str_] :=
    StringTrim[str, Whitespace.. ~~ EndOfString];

RemoveEmptyLines[l_List] :=
    Select[l, !StringMatchQ[#, StartOfString ~~ Whitespace... ~~ EndOfString]&];

WrapLine[blank_String] := {} /;
  StringMatchQ[blank, RegularExpression["^[[:space:]]*$"]];

WrapLine[line_String] := Block[{
    nLeadingSpaces = NLeadingSpaces[line],
    absSkip,
    lst,
    firstBlankStr,
    otherBlankStr
  },
  absSkip = nLeadingSpaces + relSkip;
  lst = RemoveEmptyLines[RemoveTrailingWhitespace /@ SplitLine[nLeadingSpaces, ProtectCTokens@StringTrim[line]]];
  firstBlankStr = StringJoin@Table[" ", {nLeadingSpaces}];
  otherBlankStr = StringJoin@Table[" ", {absSkip}];
  If[lst === {}, {},
     {firstBlankStr, First[lst], {"\n", otherBlankStr, #}& /@ Rest[lst]}]
];

NLeadingSpaces[line_String] := Module[{
    sp = StringPosition[line, RegularExpression["[^ ]"], 1]
  },
  If[sp === {}, 0, First@First[sp] - 1]
];

ProtectCTokens[line_String] :=
  DeleteCases[StringSplit[
      line,
      s:RegularExpression["\".*?(?<!\\\\)\"|##|<:|:>|<%|%>|(%:){1,2}|[0-9]+\\.[0-9]*|[[:alnum:]_]+|\\.\\.\\.|::|\\.\\*|[-+*/%&|^]=|(<<|>>)=?|[=!<>]=|&&|\\|\\||\\+\\+|--|->\\*?"]
      :> Hold[s]], ""];

SplitLine[fstLen_Integer, strs_List] :=
  StringJoin /@ Last@Reap[Fold[SplitString, {1, fstLen}, strs]];

SplitString[{lineN_, curLen_}, Hold[str_String]] := Module[{
    strLen = StringLength[str],
    nextLen, nextLineN = lineN
  },
  nextLen = curLen + strLen;
  If[nextLen > maxLength, nextLineN++; nextLen = absSkip + strLen];
  Sow[str, nextLineN];
  If[nextLen < maxLength, {nextLineN, nextLen}, {nextLineN + 1, absSkip}]
];

SplitString[{lineN_, curLen_}, str_String] := Block[{
    split,
    strLen, nextLen, nextLineN = lineN,
    remainder
  },
  split = Flatten@Last@Reap@SowStrings[curLen, str];
  If[split === {}, {lineN, curLen},
    nextLen = curLen + StringLength@First[split];
    If[nextLen > maxLength, nextLineN++];
    Sow[#, nextLineN++]& /@ split;
    If[Length[split] > 1, nextLen = absSkip + StringLength@Last[split]];
    If[nextLen < maxLength, {nextLineN - 1, nextLen}, {nextLineN, absSkip}]]
];

splitRegex = "[(){}\\[\\]<>,;:?~!=%^&|*+/-]";

SowStrings[_Integer, ""] := True;

SowStrings[curLen_Integer, str_String] :=
  Which[
    curLen + StringLength[str] <= maxLength, Sow[str],
    Head[StringReplace[str,
		       RegularExpression[
			   "^(.{0," <> ToString[maxLength - curLen - 1] <>
			   "}(?:" <> splitRegex <>
			   "|[^[:space:]]?(?=[[:space:]])))[[:space:]]*(.*)"] :>
		       (Sow["$1"]; remainder = "$2"; True)]] =!= String,
      SowStrings[absSkip, remainder],
    curLen > absSkip, SowStrings[absSkip, str],
    Head[StringReplace[str,
		       RegularExpression[
			   "^(.*?(?:" <> splitRegex <>
			   "|[^[:space:]]?(?=[[:space:]])))[[:space:]]*(.*)"] :>
		       (Sow["$1"]; remainder = "$2"; True)]] =!= String,
      SowStrings[absSkip, remainder],
    True, Sow[str]];

End[];

EndPackage[];
