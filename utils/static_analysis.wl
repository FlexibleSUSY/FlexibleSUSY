Needs["CodeInspector`"];

(WriteString[$Output, #];Print[];) &/@
   Map[
      CodeInspectSummarize[File[#], SeverityExclusions -> {}]&,
      FileNames["meta/*.m"]
   ];

Quit[0];
