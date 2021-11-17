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

Needs["TestSuite`", "TestSuite.m"];
Needs["Utils`", "Utils.m"];

evaluate[call_] :=
Block[{Quit, $Output = {OpenWrite[]}, out},
   call;
   out = ReadString[$Output[[1, 1]]];
   Close[$Output[[1]]];
   out
];
evaluate ~ SetAttributes ~ {HoldAll};

a;
a//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@a[1],
"Global`a: errUnknownInput:
Call
a[1]
is not supported.
Wolfram Language kernel session terminated.\n"
];

`subcontext`b;
`subcontext`b//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@`subcontext`b[1],
"Global`subcontext`b: errUnknownInput:
Call
Global`subcontext`b[1]
is not supported.
Wolfram Language kernel session terminated.\n"
];

`subcontext`subcontext`c;
`subcontext`subcontext`c//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@`subcontext`subcontext`c[1],
"Global`subcontext`subcontext`c: errUnknownInput:
Call
Global`subcontext`subcontext`c[1]
is not supported.
Wolfram Language kernel session terminated.\n"
];

d; d::usage="some text";
d//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@d[1],
"Global`d: errUnknownInput:
Usage:
some text

Call
d[1]
is not supported.
Wolfram Language kernel session terminated.\n"
];

e[x_,y_] = 3;
e//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@e[1],
"Global`e: errUnknownInput:
The behavior for case
1) e[x_,y_]
is defined only.

Call
e[1]
is not supported.
Wolfram Language kernel session terminated.\n"
];

f[x_Integer,y:{_String}:"`"] := Sin[345*x];
f//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@f[],
"Global`f: errUnknownInput:
The behavior for case
1) f[x_Integer,y:{_String}:\"`\"]
is defined only.

Call
f[]
is not supported.
Wolfram Language kernel session terminated.\n"
];

g[a_,3] := g[a,3] = 17;
g//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@g["monkey"],
"Global`g: errUnknownInput:
The behavior for case
1) g[a_,3]
is defined only.

Call
g[monkey]
is not supported.
Wolfram Language kernel session terminated.\n"
];

h /: Dot[h[2,x_,t:String:"34"], somethingelse] = "works";
h//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@h[1, 4],
"Global`h: errUnknownInput:
The behavior for case
1) h/:h[2,x_,t:String:\"34\"].somethingelse
is defined only.

Call
h[1, 4]
is not supported.
Wolfram Language kernel session terminated.\n"
];

i /: FUN[i[n_], i[m_]] := "";
i//Utils`MakeUnknownInputDefinition;
TestEquality[evaluate@FUN[i[n_], i[m_,1]],
"Global`i: errUnknownInput:
The behavior for case
1) FUN[i[n_],i[m_]]
is defined only.

Call
i[n_]
is not supported.
Wolfram Language kernel session terminated.
Global`i: errUnknownInput:
The behavior for case
1) FUN[i[n_],i[m_]]
is defined only.

Call
i[m_, 1]
is not supported.
Wolfram Language kernel session terminated.\n"
];

TestEquality[FSPermutationSign[Cycles[{{1,2}}]], -1];
TestEquality[FSPermutationSign[Cycles[{{1,3,2}}]], 1];
TestEquality[FSPermutationSign[Cycles[{{1,3,2},{5,6}}]], -1];

PrintTestSummary[];
