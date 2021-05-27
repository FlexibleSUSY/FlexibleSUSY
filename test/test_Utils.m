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

Global`time = AbsoluteTime[];

(* THERE IS A NEED TO CHANGE BEHAVIOR OF LOCKED FILE **************************)
fileName = FileNameJoin@{Directory[],"test","test_Utils.temp"};
Off[Syntax::sntufn];
subKernel = LaunchKernels[1];
DistributeDefinitions[fileName];
ParallelEvaluate[
   Needs["Utils`",FileNameJoin@{Directory[],"meta","Utils.m"}];
   Put[Definition@Utils`Private`internalAssertOrQuit,fileName];
   code = StringJoin@Riffle[Drop[Utils`ReadLinesInFile@fileName,2],"\n"];
   code=StringReplace[code,
      {
         "\\[Raw"~~Shortest@__~~"m"->"",
         "Quit[1];"->"","FSFancyLine[];"->"",
         "WriteString"~~__~~x:"StringJoin[Utils`Private`str]"~~"]":>
            "(`test`out = `test`out<>"<>x<>")",
         "\\nWolfram Language kernel session "->""
      }
   ];
   fileHandle = OpenWrite@fileName;
   fileHandle ~ WriteString ~ code;
   Close@fileHandle;
];
CloseKernels@subKernel;
On[Syntax::sntufn];
(* LOAD TEST DEFINITIONS ******************************************************)
`test`out="";
Get@fileName;
Attributes[Utils`Private`internalAssertOrQuit] = {HoldAll,Protected,Locked};
DeleteFile@fileName;

Off[SetDelayed::write,Attributes::locked];
Needs["TestSuite`", "TestSuite.m"];
Needs["Utils`", "Utils.m"];
On[SetDelayed::write,Attributes::locked];

(* TEST MULTIDIMENSIONAL CONTEXT***********************************************)
a;
a//Utils`MakeUnknownInputDefinition;
a[1];
TestEquality[`test`out,
"Global`a: errUnknownInput:
Call
a[1]
is not supported.terminated.\n"
];
`test`out="";

`subcontext`b;
`subcontext`b//Utils`MakeUnknownInputDefinition;
`subcontext`b[1];
TestEquality[`test`out,
"Global`subcontext`b: errUnknownInput:
Call
Global`subcontext`b[1]
is not supported.terminated.\n"
];
`test`out="";

`subcontext`subcontext`c;
`subcontext`subcontext`c//Utils`MakeUnknownInputDefinition;
`subcontext`subcontext`c[1];
TestEquality[`test`out,
"Global`subcontext`subcontext`c: errUnknownInput:
Call
Global`subcontext`subcontext`c[1]
is not supported.terminated.\n"
];
`test`out="";

d; d::usage="some text";
d//Utils`MakeUnknownInputDefinition;
d[1];
TestEquality[`test`out,
"Global`d: errUnknownInput:
Usage:
some text

Call
d[1]
is not supported.terminated.\n"
];
`test`out="";

(* TEST INPUT TYPES ***********************************************************)

e[x_,y_] = 3;
e//Utils`MakeUnknownInputDefinition;
e[1];
TestEquality[`test`out,
"Global`e: errUnknownInput:
The behavior for case
1) e[x_,y_]
is defined only.

Call
e[1]
is not supported.terminated.\n"
];
`test`out="";

f[x_Integer,y:{_String}:"`"] := Sin[345*x];
f//Utils`MakeUnknownInputDefinition;
f[];
TestEquality[`test`out,
"Global`f: errUnknownInput:
The behavior for case
1) f[x_Integer,y:{_String}:\"`\"]
is defined only.

Call
f[]
is not supported.terminated.\n"
];
`test`out="";

g[a_,3] := g[a,3] = 17;
g//Utils`MakeUnknownInputDefinition;
g["monkey"];
TestEquality[`test`out,
"Global`g: errUnknownInput:
The behavior for case
1) g[a_,3]
is defined only.

Call
g[monkey]
is not supported.terminated.\n"
];
`test`out="";

h /: Dot[h[2,x_,t:String:"34"], somethingelse] = "works";
h//Utils`MakeUnknownInputDefinition;
h[1,4];
TestEquality[`test`out,
"Global`h: errUnknownInput:
The behavior for case
1) h/:h[2,x_,t:String:\"34\"].somethingelse
is defined only.

Call
h[1, 4]
is not supported.terminated.\n"
];
`test`out="";

i /: Print[i[n_], i[m_]] := "";
i//Utils`MakeUnknownInputDefinition;
Print[i[n_], i[m_,1]];
TestEquality[`test`out,
"Global`i: errUnknownInput:
The behavior for case
1) Print[i[n_],i[m_]]
is defined only.

Call
i[n_]
is not supported.terminated.
Global`i: errUnknownInput:
The behavior for case
1) Print[i[n_],i[m_]]
is defined only.

Call
i[m_, 1]
is not supported.terminated.\n"
];
`test`out="";

TestEquality[FSPermutationSign[Cycles[{{1,2}}]], -1];
TestEquality[FSPermutationSign[Cycles[{{1,3,2}}]], 1];
TestEquality[FSPermutationSign[Cycles[{{1,3,2},{5,6}}]], -1];

Print[StringJoin[">>test>> done in ",ToString@N[AbsoluteTime[]-Global`time,{Infinity,3}]," seconds.\n"]];

PrintTestSummary[];
