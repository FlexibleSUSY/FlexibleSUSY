(* Invalid paths to the Lilith and HiggsTools databases must give an error
   instead of aborting the kernel. *)

Needs["TestSuite`", "TestSuite.m"];

Get["models/SM/SM_librarylink.m"];

handle = FSSMOpenHandle[
   fsSettings -> { calculateStandardModelMasses -> 1 },
   fsModelParameters -> { LambdaIN -> 0.2, Qin -> 173.3, QEWSB -> 173.3 }
];

FSSMCalculateSpectrum[handle];

tooLong = "/tmp/" <> StringRepeat["a", 300];
missing = "/tmp/flexiblesusy_no_such_path";
notADataset = CreateFile[];

If[Head[FSSMCallLilith] === LibraryFunction,
   Print["testing Lilith"];
   Quiet[
      TestEquality[FSSMCallLilith[handle, tooLong], $Failed];
      TestEquality[FSSMCallLilith[handle, missing], $Failed];
   ];
];

If[Head[FSSMCallHiggsTools] === LibraryFunction,
   Print["testing HiggsTools"];
   Quiet[
      TestEquality[FSSMCallHiggsTools[handle, tooLong, tooLong], $Failed];
      TestEquality[FSSMCallHiggsTools[handle, missing, missing], $Failed];
      TestEquality[FSSMCallHiggsTools[handle, notADataset, tooLong], $Failed];
      TestEquality[FSSMCallHiggsTools[handle, tooLong, notADataset], $Failed];
   ];
];

DeleteFile[notADataset];

FSSMCloseHandle[handle];

PrintTestSummary[];
