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



time = AbsoluteTime[];

Module[{
      file = OpenRead@FileNameJoin@{Directory[], "Makefile"},
      cond = False, str,
      regex1 = RegularExpression["MODELS\\s+:=\\s+(models/\\w+\\s*)+"],
      regex2 = RegularExpression@"models/(\\w+)"
   },

   While[!cond,
      str = ReadLine[file];
      cond = StringMatchQ[str, regex1];
      If[cond, models = StringCases[str, regex2 :> "$1"]];
   ];

   Close@file;

];

models = Module[{
      file = OpenRead@FileNameJoin@{Directory[], "models", #, "start.m"},
      cond = False, str, out,
      regex1 = RegularExpression["Start..(\\w+)..;"]
   },

   While[!cond,
      str = ReadLine[file];
      cond = StringMatchQ[str, regex1];
      If[cond, out = StringCases[str, regex1 :> "$1"]];
   ];

   Close@file;
   out[[1]]
]&/@models;

Print@models;

AppendTo[$ContextPath, "SARAH`"];
AppendTo[$ContextPath, "Susyno`LieGroups`"];

Needs["TestSuite`", "TestSuite.m"];

CheckRude[str_String, type:Identity|Clear|Replace|Inverse] :=
Module[{
      model = str, fsmodel = fsstr, dir = Directory[], path = $Path, pos3, pos4, prt,
      time = AbsoluteTime[], jobs, kers, res, ok = True,
      func = Switch[type,
         Identity, Identity,
         Clear, #/.fse_[{__}]:>fse&,
         Replace, #/.{SARAH`gt1:>RandomInteger@{1,3},
                      SARAH`gt2:>RandomInteger@{1,3},
                      SARAH`gt3:>RandomInteger@{1,3}}&,
         Inverse, Prepend[Rest@#, SARAH`AntiField@First@#]&@*(#/.fse_[{__}]:>fse&)
      ]
   },
   Print["test for ", str];

   prt = GetPartitions@str;

   While[ok,

      Check[
         kers = LaunchKernels@4;
         ok = False;,

         CloseKernels@kers;
         ok = True;
         Print@"Relaunching kernels",
         LinkObject::linkd
         ];

   ];

   jobs = Table[
      pos3 = prt[[1, i]];
      pos4 = prt[[2, i]];
      ParallelSubmit[{model, dir, path, func, pos3, pos4},
         SetDirectory[dir];
         $Path = path;
         Block[{Print},
            Needs["TestSuite`", "TestSuite.m"];
            Needs["Cache`", "Cache.m"];
            SARAH`Start@model;
         ];
         SARAH`RXi[_] = 1;
         SA`CurrentStates = EWSB;
         SARAH`$sarahCurrentOutputMainDir =
            FileNameJoin@{Directory[], "Output", model};
         all3 = First/@Get@FileNameJoin@{SARAH`$sarahCurrentOutputMainDir,
            "EWSB", "Vertices", "VertexList3.m"};
         all4 = First/@Get@FileNameJoin@{SARAH`$sarahCurrentOutputMainDir,
            "EWSB", "Vertices", "VertexList4.m"};
         vertices = func/@Join[all3[[pos3]], all4[[pos4]]];
         fsS1 = SA`subUnitaryCondition;
         fsS2 = SA`subUnitaryCondition /. Susyno`LieGroups`b :> Susyno`LieGroups`a;
         Get[Cache`GetVertex];
         Do[
            perm = Permutations@vertices[[vert]];
            Do[
               cache = TrigExpand@Expand[Simplify[Expand@Cache`GetVertex@perm[[subv]]] /. fsS1 /. fsS2];
               sarah = TrigExpand@Expand[Simplify[SARAH`Vertex@perm[[subv]]] /. fsS1 /. fsS2];
               If[MatchQ[sarah,{{_, _,_, _}, {0, 1|SARAH`g[_, _]}}],
                  bad = " 4part case: ";
                  sarah = sarah /. SARAH`g[_, _] :> 1;
                  cache = cache /. SARAH`g[_, _] :> 1;,
                  bad = ": ";
               ];
               Print[vert, "/", Length@vertices, " ", subv, "/", Length@perm, bad, perm[[subv]]];
               If[MatchQ[sarah,{{_, _, _}, {_, HoldPattern[-2*SARAH`Mom[_,_]+2*SARAH`Mom[_,_]]}}],
                  Print["Wrong sarah expression for SSV vertex!"];
                  ,
                  TestEquality[cache, sarah];
               ];,

               {subv, Length@perm}
            ];,

            {vert, Length@vertices}
         ];
         If[GetNumberOfFailedTests[] === 0, Put[Cache`GetVertex]];
         PrintTestSummary[];
         GetNumberOfFailedTests[]
      ],
      {i, 4}
   ];
   res = WaitAll@jobs;
   CloseKernels@kers;

   Print["Time needed: ", ToString[N[(AbsoluteTime[]-time)/60.,{Infinity,3}]], "m.\n"];

   TestEquality[res, {0, 0, 0, 0}];
];

GetPartitions[model_String] :=
With[{dir = FileNameJoin@{Directory[], "Output", model, "EWSB", "Vertices"}},
Module[{
      ker = LaunchKernels@1, par3, par4
   },
   {{par3, par4}} = ParallelEvaluate[
      f[{min_, max_, rest___}] :=
         {Range @@ min, Sequence @@ f[{{Last@min + 1, Last@min + First@max}, rest}]};
      f[{{min_, max_}}] :=
         {Range[min, max]};

      FairSplit[l_, n_] :=
      Module[{subs},
         subs = Table[{Floor[l/n]}, n];
         Do[++subs[[i]], {i, Mod[l, n]}];
         f@subs
      ];
      {
         FairSplit[Length@Get@FileNameJoin@{dir, "VertexList3.m"}, 4],
         FairSplit[Length@Get@FileNameJoin@{dir, "VertexList4.m"}, 4]
      }, ker, DistributedContexts -> None];

   CloseKernels@ker;
   {par3, par4}
]
];


(
   CheckRude[#, Inverse];
   CheckRude[#, Clear];
   CheckRude[#, Replace];
   CheckRude[#, Identity];
)& /@ models ;

Print[StringJoin[">>test>> done in ",ToString[N[(AbsoluteTime[]-Global`time)/60.,{Infinity,3}]]," m.\n"]];

PrintTestSummary[];
