(* Invalid paths to the Lilith and HiggsTools databases must give an error
   instead of aborting the kernel.  Every error must give a message with a
   defined text and $Failed.  If Lilith or HiggsTools isn't enabled, the
   interface must load without messages, and a call must give an error. *)

Needs["TestSuite`", "TestSuite.m"];

TestEquality[Check[Get["models/SM/SM_librarylink.m"]; "loaded", "messages"], "loaded"];

configH = Import["config/config.h", "Text"];
enabledQ[tool_String] := StringContainsQ[configH, "#define ENABLE_" <> tool <> " 1"];

TestEquality[Head[FSSMCallLilith], LibraryFunction];
TestEquality[Head[FSSMCallHiggsTools], LibraryFunction];
TestEquality[Head[FSSMCallLilith::error], String];
TestEquality[Head[FSSMCallHiggsTools::error], String];

handle = FSSMOpenHandle[
   fsSettings -> { calculateStandardModelMasses -> 1 },
   fsModelParameters -> { LambdaIN -> 0.2, Qin -> 173.3, QEWSB -> 173.3 }
];

FSSMCalculateSpectrum[handle];

tooLong = "/tmp/" <> StringRepeat["a", 300];
missing = "/tmp/flexiblesusy_no_such_path";
notADataset = CreateFile[];

If[enabledQ["LILITH"],
   Print["testing Lilith"];
   Quiet[
      TestEquality[FSSMCallLilith[handle, tooLong], $Failed];
      TestEquality[FSSMCallLilith[handle, missing], $Failed];
   ];
   TestEquality[Head[FSSMCallLilith::error], String];
   ,
   Print["testing Lilith (not enabled)"];
   TestEquality[
      Quiet[Check[FSSMCallLilith[handle, missing], "error", {FSSMCallLilith::error}],
            {FSSMCallLilith::error}],
      "error"];
   TestEquality[Quiet[FSSMCallLilith[handle, missing]], $Failed];
];

If[enabledQ["HIGGSTOOLS"],
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
   ,
   Print["testing HiggsTools (not enabled)"];
   TestEquality[
      Quiet[Check[FSSMCallHiggsTools[handle, missing, missing], "error",
                  {FSSMCallHiggsTools::error}], {FSSMCallHiggsTools::error}],
      "error"];
   TestEquality[Quiet[FSSMCallHiggsTools[handle, missing, missing]], $Failed];
];

DeleteFile[notADataset];

FSSMCloseHandle[handle];

PrintTestSummary[];
