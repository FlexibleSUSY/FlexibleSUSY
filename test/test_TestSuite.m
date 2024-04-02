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

failedTests = 0;

AbortAssert::trace = "Assertion `` failed.";
AbortAssert /: On[AbortAssert] := On[AbortAssert::trace];
AbortAssert /: Off[AbortAssert] := Off[AbortAssert::trace];
SetAttributes[AbortAssert, {HoldFirst}];
AbortAssert[test_] := Check[TrueQ[test] || Message[AbortAssert::trace, HoldForm[test]], failedTests++];

(* test TestEquality *)
AbortAssert[TestEquality[0, 0]];
AbortAssert[!TestEquality[0, 1]];
AbortAssert[TestEquality["ab", "ab"]];
AbortAssert[!TestEquality["ab", "ba"]];

(* test TestCloseRel *)
AbortAssert[TestCloseRel[0, 0, 1]];
AbortAssert[TestCloseRel[1, 1, 1]];
AbortAssert[TestCloseRel[1, 1 + 1*^-10, 1*^-10]];
AbortAssert[!TestCloseRel[1, 1 + 1*^-10, 1*^-11]];
AbortAssert[TestCloseRel[1*^-100, 1*^-100 + 1*^-200, 1*^-100]];

(* test TestLowerThan *)
AbortAssert[!TestLowerThan[0, 0]];
AbortAssert[!TestLowerThan[1, 0]];
AbortAssert[TestLowerThan[0, 1]];

(* Summary *)
If[failedTests == 0,
   Print["All tests passed."],
   Print["Number of failed tests: ", failedTests]
];

Quit[failedTests];
