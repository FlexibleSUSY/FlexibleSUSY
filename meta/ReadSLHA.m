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

BeginPackage["ReadSLHA`"];

ReadSLHAFile::usage="Reads SLHA file and returns list of input and
 output parameters.

 Usage:  ReadSLHAFile[fileName, parameters]

 - fileName - name of SLHA file to be read
 - parameters - list of 3-component lists specifying how
   to interpret the SLHA file content:
   #1 is the parameter name (a symbol)
   #2 is the parameter dimension (a list of positive integers)
   #3 is the {block, entry, ...} in the SLHA input

 Example:

 parameters = {
    {Qin  , {0}   , {EXTPAR, 0}},
    {m0   , {0}   , {MINPAR, 1}},
    {CpHPP, {0}   , {EFFHIGGSCOUPLINGS, 25, 22, 22}},
    {k    , {3}   , KappaIn},
    {Yu   , {3, 3}, YuIN}
 };

 ReadSLHAFile[\"file.slha\", parameters]
";

ReadSLHAStream::usage="Reads stream with SLHA input and returns list
 of input and output parameters.

 Usage:  ReadSLHAStream[streamName, parameters]

 - streamName - name of the stream
 - parameters - see documentation of ReadSLHAFile[].
";

ReadSLHAString::usage="Reads string with SLHA input and returns list
 of input and output parameters.

 Usage:  ReadSLHAString[str, parameters]

 - str - string in SLHA format (set of blocks)
 - parameters - see documentation of ReadSLHAFile[].
";

Begin["`Private`"];

IsDataLine[str_String] :=
    !StringMatchQ[str, StartOfString ~~ "#" ~~ ___];

BlockStarts[str_String, blockName_String] :=
    StringMatchQ[str, StartOfString ~~ "BLOCK" ~~ Whitespace ~~ blockName ~~ WordBoundary ~~ ___,
                 IgnoreCase -> True];

BlockStarts[str_String] :=
    StringMatchQ[str, StartOfString ~~ "BLOCK" ~~ WordBoundary ~~ ___,
                 IgnoreCase -> True];

floatRegex = "[+\-]?(?:[0-9]*)(?:\\.[0-9]*)?(?:[eE][+\-]?[0-9]+)?";

ReadBlock[stream_, blockName_String] :=
    Module[{line = "", block = "", inBlock = False},
           SetStreamPosition[stream, 0];
           While[(line = Read[stream, String]) =!= EndOfFile,
                 If[BlockStarts[line, blockName],
                    inBlock = True,
                    If[BlockStarts[line],
                       inBlock = False
                      ];
                   ];
                 If[!inBlock || !IsDataLine[line], Continue[]];
                 block = block <> line <> "\n";
                ];
           block
          ];

ComposeFixedRegEx[idx__] :=
    StringExpression[StartOfString, Whitespace,
                     Sequence @@ (StringExpression[ToString[#], Whitespace]& /@ {idx})];

ToWolframExpr[str_String] :=
    ToExpression[StringReplace[str, {"E"|"e" -> "*10^", "=" -> ""}]];

ReadIndexFromBlock[block_String, idx___] :=
    Module[{value = 0, values, line = "", val,
            stream = StringToStream[block], patt = ComposeFixedRegEx[idx]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 If[BlockStarts[line], Continue[]];
                 values = StringCases[line, patt ~~ v:RegularExpression[floatRegex] :> v];
                 If[values =!= {},
                    val = ToWolframExpr[First[values]];
                    If[NumberQ[val],
                       value = val;
                      ];
                   ];
                ];
           value
          ];

IsIndex[i_?IntegerQ] := True;
IsIndex[_] := False;
IsIndex[indices_List] := And @@ (IsIndex /@ indices);
IsIndex[indices__] := IsIndex[{indices}];

WithinRange[i_Integer, dim_Integer] := i >= 1 && i <= dim;
WithinRange[idx_List, dim_List] := Length[idx] == Length[dim] && (And @@ MapThread[WithinRange, {idx, dim}]);

ReadMatrixFromBlock[block_String, dims__] :=
    Module[{matrix = Array[0&, {dims}],
            ndims = Length[{dims}],
            line = "", values, indices, value,
            stream = StringToStream[block]},
           While[(line = Read[stream, String]) =!= EndOfFile,
                 If[BlockStarts[line], Continue[]];
                 values = StringSplit[line];
                 If[Length[values] > ndims,
                    indices = ToWolframExpr /@ Take[values, ndims];
                    value = ToWolframExpr[values[[ndims + 1]]];
                    If[IsIndex[indices] && NumberQ[value] && WithinRange[indices, {dims}],
                       matrix[[Sequence @@ indices]] = value;
                      ];
                   ];
                ];
           matrix
          ];

ReadParameter[stream_, par_, {} | {0|1} | 0 | 1, {block_, idx___}] :=
    ReadIndexFromBlock[ReadBlock[stream, ToString[block]], idx];

ReadParameter[stream_, par_, {dims__}, block_] :=
    ReadMatrixFromBlock[ReadBlock[stream, ToString[block]], dims];

ReadParameter[stream_, {par_, type_, block_}] :=
    par -> ReadParameter[stream, par, type, block];

ReadSLHAStream[stream_, parameters_List] :=
    ReadParameter[stream, #]& /@ parameters;

ReadSLHAFile[fileName_String, parameters_List] :=
    Module[{file, values},
           file = OpenRead[fileName];
           values = ReadSLHAStream[file, parameters];
           Close[file];
           values
          ];

ReadSLHAString[slha_String, parameters_List] :=
    ReadSLHAStream[StringToStream[slha], parameters];

End[];

EndPackage[];
