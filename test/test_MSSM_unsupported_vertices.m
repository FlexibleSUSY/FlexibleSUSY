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

(* Tests the treatment of vertices that CXXDiagrams can't evaluate,
   such as the four-gluon vertex, in the MSSM. *)

Get["utils/load-FlexibleSUSY.m"];

Needs["TestSuite`", "TestSuite.m"];
Needs["CXXDiagrams`", "CXXDiagrams.m"];
Needs["Decays`", "Decays.m"];

Start["MSSM"];

{susyBetaFunctions, susyBreakingBetaFunctions} =
  FlexibleSUSY`ReadSARAHBetaFunctions[];
allParameters =
  FlexibleSUSY`SetupModelParameters[susyBetaFunctions,
   susyBreakingBetaFunctions];

{massMatrices, Lat$massMatrices} =
  FlexibleSUSY`SetupMassMatrices[allParameters];
FlexibleSUSY`SetupOutputParameters[massMatrices];

Print["testing the Lorentz structure of the four-gluon vertex ..."];

(* The four-gluon vertex has a different colour structure for each pair
   of metric tensors, so it is not of the form colour * Lorentz. *)
TestEquality[
   Check[CXXDiagrams`VertexTypeForFields[{VG, VG, VG, VG}], $Failed],
   CXXDiagrams`Private`CompositeQuadrupleVectorVertex
];
TestEquality[
   MemberQ[CXXDiagrams`VertexTypes[], CXXDiagrams`VertexTypeForFields[{VG, VG, VG, VG}]],
   False
];

(* vertices with a single colour structure are still supported *)
TestEquality[CXXDiagrams`VertexTypeForFields[{VG, VG, VG}], CXXDiagrams`TripleVectorVertex];
TestEquality[
   CXXDiagrams`VertexTypeForFields[{conj[VWm], VP, VP, VWm}],
   CXXDiagrams`QuadrupleVectorVertex
];

Print["testing the warning about discarded diagrams ..."];

(* Su -> Su VG is skipped by GetFinalStatesForParticle (it vanishes),
   but its 1-loop diagrams include tadpoles of a gluon attached to the
   external gluon through the four-gluon vertex. *)
oneLoop = Join @@ Cases[
   Decays`Private`GetAllContributingGraphsForDecay[Su, {Su, VG}, 1],
   {1, toposAndDiags_} :> toposAndDiags];

discarded = Flatten[
   Function[{topo, diags},
      {topo, #}& /@ Select[diags, !Decays`Private`IsSupportedDiagram[#]&]] @@@ oneLoop, 1];

TestEquality[discarded =!= {}, True];
TestEquality[And @@ (Decays`Private`IsMasslessTadpoleDiagram @@@ discarded), True];

(* massless tadpoles vanish: no warning *)
TestEquality[
   Reap[Block[{Utils`FSFancyWarning = Sow[{##}]&},
      Decays`Private`WarnAboutDiscardedDiagrams[Su, {Su, VG}, oneLoop]]][[2]],
   {}
];

(* a discarded diagram that may be non-zero is reported *)
TestEquality[
   Length @ First @ Reap[Block[{Utils`FSFancyWarning = Sow[{##}]&,
                                Decays`Private`IsMasslessTadpoleDiagram = False&},
      Decays`Private`WarnAboutDiscardedDiagrams[Su, {Su, VG}, oneLoop]]][[2]],
   1
];

PrintTestSummary[];
