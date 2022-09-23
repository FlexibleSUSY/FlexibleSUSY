files =
   DeleteCases[
      If[!FileExistsQ[#], Print["Argument ", #, " is not a file name"], #]& /@ DeleteDuplicates@Drop[$ScriptCommandLine, 1],
      Null
   ];

If[files === {},
   Print["Error: No files specified"];
   Print["Usage: wolframscript -f static_analysis.wl FILE(S)"];
   Quit[1];
];

Needs["CodeInspector`"];

(WriteString[$Output, #]; Print[];) &/@
   Map[
      CodeInspectSummarize[File[#], SeverityExclusions -> {}]&,
      files
   ];

Quit[0];
