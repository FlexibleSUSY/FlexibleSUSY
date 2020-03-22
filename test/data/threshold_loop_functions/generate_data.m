(* calculate func[args] with precision prec *)
GeneratePoint[func_, prec_, args__] :=
    N[#, prec]& @ { args, func[args] }

(* calculate func[arg] with precision prec on a grid *)
GenerateGridData[func_, narg_, prec_, {min_, max_, step_}] :=
    Module[{vals = Table[x, {x, min, max, step}]},
           GeneratePoint[func, prec, Sequence @@ #]& /@ Tuples[vals, narg]
    ]

(* calculate func with precision prec on a grid close to start *)
GenerateLimitData[func_, narg_, prec_, n_, frac_, dir_, start_] :=
    Module[{vals = Table[start + dir frac^k, {k, 1, n, 1}]},
           GeneratePoint[func, prec, Sequence @@ #]& /@ Tuples[vals, narg]
    ]

ExportData[func_, narg_, prec_, zero_] :=
    Module[{data, filename = ToString[func] <> ".txt", step = 1/20},
           Print["Generating data for " <> ToString[func]];
           data = Join[
               GenerateGridData[func, narg, prec, {zero, 5, step}],
               GenerateLimitData[func, narg, prec, prec, 1/10, -1, 1],
               GenerateLimitData[func, narg, prec, prec, 1/10, +1, 1],
               If[PossibleZeroQ[zero],
                  GenerateLimitData[func, narg, prec, prec, 1/10, +1, 0],
                  {}
               ]
           ];
           Print["Writing data to ", filename];
           Export[filename, data, "Table"];
    ]

$MaxExtraPrecision = 10000;
precision = 15;

ExportData[F1, 1, precision, 0];
ExportData[F2, 1, precision, 0];
ExportData[F3, 1, precision, 0];
ExportData[F4, 1, precision, 0];
ExportData[F5, 1, precision, 0];
ExportData[F6, 1, precision, 0];
ExportData[F7, 1, precision, 0];
ExportData[F8, 2, precision, 0];
ExportData[F9, 2, precision, 1/20];

ExportData[f1, 1, precision, 0];
ExportData[f2, 1, precision, 0];
ExportData[f3, 1, precision, 0];
ExportData[f4, 1, precision, 0];
ExportData[f5, 2, precision, 0];
ExportData[f6, 2, precision, 0];
ExportData[f7, 2, precision, 0];
ExportData[f8, 2, precision, 0];
