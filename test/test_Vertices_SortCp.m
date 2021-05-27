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
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["MSSM"];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;

Print["testing SortCp[] ..."];

TestEquality[SortCp @ Cp[conj[USd[{gO2}]], VG, Sd[{gI2}]],
		      Cp[Sd[{gI2}], conj[USd[{gO2}]], VG]];

TestEquality[SortCp @ Cp[conj[VWm], Hpm[{gI1}], hh[{gI2}]],
		      Cp[hh[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[conj[VWm], Hpm[{gI1}], Ah[{gI2}]],
		      Cp[Ah[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[conj[VWm], hh[{gI2}], Hpm[{gI1}]],
		      Cp[hh[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[Ah[{gI2}], conj[VWm], Hpm[{gI1}]],
		      Cp[Ah[{gI2}], Hpm[{gI1}], conj[VWm]]];

TestEquality[SortCp @ Cp[VWm, conj[VWm], conj[VWm], VWm][1],
		      Cp[conj[VWm], conj[VWm], VWm, VWm][2]];

TestEquality[SortCp @ Cp[VWm, conj[VWm], conj[VWm], VWm][2],
		      Cp[conj[VWm], conj[VWm], VWm, VWm][3]];

TestEquality[SortCp @ Cp[VWm, conj[VWm], conj[VWm], VWm][3],
		      Cp[conj[VWm], conj[VWm], VWm, VWm][1]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], VWm][PL],
		     -Cp[bar[Cha[{gI1}]], UChi[{gO1}], VWm][PR]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], VWm][PR],
		     -Cp[bar[Cha[{gI1}]], UChi[{gO1}], VWm][PL]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], Hpm][PL],
		      Cp[bar[Cha[{gI1}]], UChi[{gO1}], Hpm][PL]];

TestEquality[SortCp @ Cp[UChi[{gO1}], bar[Cha[{gI1}]], Hpm][PR],
		      Cp[bar[Cha[{gI1}]], UChi[{gO1}], Hpm][PR]];

TestEquality[SortCp @ Cp[conj[USd[{gO2}]], Glu[{1}], Fd[{gI2}]][PL],
		      Cp[Glu, Fd[{gI2}], conj[USd[{gO2}]]][PL]];

TestEquality[SortCp @ Cp[VZ, hh[{gI1}], Ah[{2}]],
		      Cp[Ah[{2}], hh[{gI1}], VZ]];

(* For vertex for fields get coefficient in front of
   SARAH`g[idx[[1]], idx[[2]]], SARAH`g[idx[[3]], idx[[4]]] *)
Get4VectorVertexPart[fields_List, idx_List] :=
   Module[{result},
      result =
         Select[
            Delete[SARAH`Vertex[fields], 1],
            MatchQ[
               Last@#,
               HoldPattern[Times][SARAH`g[idx[[1]], idx[[2]]], SARAH`g[idx[[3]], idx[[4]]]]
            ]&
         ];
      If[Length@result =!= 1, Print["Error! Couldn't identify requested part of 4-vector vertex"], First@First@result]
   ];

(* Check if SortCp works correctly for fields, where fields are
   4-vector bosons *)
CheckSortCpForFieldList[fieldsIn_] /; Length[fieldsIn]===4 :=
   Module[{fieldPermuations, indexPermutations, testInputs},
      fieldPermuations = Permutations[fieldsIn];
      (* all combinations of {lt1, lt2, lt3, lt4} without
         elements which can be obtained by interchanging 1 <-> 2
         and/or 3 <-> 4 *)
      indexPermutations =
         DeleteDuplicatesBy[
            Permutations[{lt1, lt2, lt3, lt4}],
            Flatten[Sort /@ TakeList[#, {2, 2}]]&
         ];
      testInputs = Tuples[{fieldPermuations, indexPermutations}];
      Module[{fields = #1, sortedFields, indices = #2, sortedIndices},
         sortedCoupling =
           SortCp[(Cp @@ fields)[SARAH`g[indices[[1]], indices[[2]]] SARAH`g[indices[[3]], indices[[4]]]]];
         {sortedFields, sortedIndices} =
            sortedCoupling /. SARAH`Cp[y___][x___] :> {{y}, x};
         (* g[] g[] -> {} *)
         sortedIndices =
            sortedIndices /. Times[x_, y_] :> Flatten[List @@@ {x, y}, 1];
         TestEquality[
            Get4VectorVertexPart[fields,       indices],
            Get4VectorVertexPart[sortedFields, sortedIndices]
         ];
      ] & @@@ testInputs
   ];

(* this tests takes ~1h will all combinations, so we take every 10th *)
CheckSortCpForFieldList /@
   Select[Tuples[{VP, VZ, VWm, conj[VWm]}, 4], IsNonZeroVertex][[;; ;; 10]];

Print["testing Vertices`Private`RestoreBarOnMajorana[] ..."];

RBOM := Vertices`Private`RestoreBarOnMajorana;

Block[{SARAH`bar},

TestEquality[RBOM[{UChi[{gO1}], Cha[{gI2}], conj[Hpm[{gI1}]]}, PR],
	     {bar @ UChi[{gO1}], Cha[{gI2}], conj[Hpm[{gI1}]]}];

TestEquality[RBOM[{bar[Cha[{gI1}]], UChi[{gO1}], Hpm}, PR],
	     {bar[Cha[{gI1}]], bar @ UChi[{gO1}], Hpm}];

TestEquality[RBOM[{UChi[{gO2}], Cha[{gI2}], conj[VWm]},
		  LorentzProduct[gamma[lt3], PR]],
	     {bar @ UChi[{gO2}], Cha[{gI2}], conj[VWm]}];

TestEquality[RBOM[{bar[Cha[{gI1}]], UChi[{gO1}], VWm},
		  LorentzProduct[gamma[lt3], PR]],
	     {bar[Cha[{gI1}]], (* w/o bar *) UChi[{gO1}], VWm}];

];

DeleteDirectory[SARAH`SARAH[OutputDirectory], DeleteContents -> True];

PrintTestSummary[];
