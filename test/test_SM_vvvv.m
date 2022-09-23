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
Block[{Print}, (* We do not test SARAH, so let us reduce output. *)
   Needs["SARAH`"];
   Start@"SM";
];

FlexibleSUSY`FSEnableColors = False; (* Desired by Utils`FSFancyWarning *)
FlexibleSUSY`FSEigenstates = EWSB; (* Required by TreeMasses`Private`IsOfType *)
Needs["Wrappers`", "Wrappers.m"];

check[fields:{_, _, _, _}] := Module[{lhs, rhs, vertex, stream, out},
   Print["\ncheck"];
   out = $Output;
   stream = OpenWrite[];
   $Output = {stream};
   vertex = Vertex@fields;
   lhs = ReadString[stream[[1]]];
   $Output = out;
   Close[stream];

   stream = OpenWrite[];
   AppendTo[$Output, stream];
   Utils`FSFancyWarning@Wrappers`Private`warnVVVV[fields, vertex];
   rhs = ReadString[stream[[1]]];
   $Output = Most[$Output];
   Close[stream];
   {lhs, rhs}
];

Utils`FSFancyLine[];
Print["Wrong input"];
TestEquality@@ check@{conj[VWp[{lt1}]], VZ[{lt2}], VWp[{lt3}], VZ[{1}]};
TestEquality@@ check@{conj[VWp[{lt1}]], VZ[{lt2}], VWp[{lt3}], VZ[{"oops"}]};
TestEquality@@ check@{conj[VWp[{lt1}]], VZ[{lt1}], VWp[{lt3}], VZ};
TestEquality@@ check@{conj[VWp], VZ[{lt1}], VWp[{lt3}], VZ};

Utils`FSFancyLine[];
Print["Correct but weird input"];
TestEquality[First@check@{conj[VWp[{lt1}]], VZ[{lt2}], VWp[{lt3}], VZ},
             EndOfFile];
TestEquality[First@check@{conj[VWp[{lt2}]], VZ[{lt1}], VWp[{lt3}], VZ},
             EndOfFile];

PrintTestSummary[];
