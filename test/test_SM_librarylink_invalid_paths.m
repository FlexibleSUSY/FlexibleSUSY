(* Invalid paths to the Lilith and HiggsTools databases must give an error
   instead of aborting the kernel.  Every error must give a message with a
   defined text and $Failed. *)

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
   TestEquality[Head[FSSMCallLilith::error], String];
];

If[Head[FSSMCallHiggsTools] === LibraryFunction,
   Print["testing HiggsTools"];
   Quiet[
      TestEquality[FSSMCallHiggsTools[handle, tooLong, tooLong], $Failed];
      TestEquality[FSSMCallHiggsTools[handle, missing, missing], $Failed];
      TestEquality[FSSMCallHiggsTools[handle, notADataset, tooLong], $Failed];
      TestEquality[FSSMCallHiggsTools[handle, tooLong, notADataset], $Failed];
      (* empty HiggsBounds path: flexiblesusy::Error from call_higgstools,
         which used to give LibraryFunctionError without a message *)
      TestEquality[FSSMCallHiggsTools[handle, notADataset, ""], $Failed];
   ];
   TestEquality[Head[FSSMCallHiggsTools::error], String];
   TestEquality[Head[FSSMCallHiggsTools::warning], String];
   TestEquality[
      Quiet[Check[FSSMCallHiggsTools[handle, notADataset, ""], "error",
                  {FSSMCallHiggsTools::error}], {FSSMCallHiggsTools::warning}],
      "error"];
   (* only one of the two paths exists: warning about the other one *)
   TestEquality[
      Quiet[Check[FSSMCallHiggsTools[handle, notADataset, missing], "warning",
                  {FSSMCallHiggsTools::warning}], {FSSMCallHiggsTools::error}],
      "warning"];
];

DeleteFile[notADataset];

FSSMCloseHandle[handle];

PrintTestSummary[];
