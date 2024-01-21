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

Get["models/CMSSM/CMSSM_librarylink.m"];

(* Create a handle to a model given the input parameters.
   See Options[FSCMSSMOpenHandle] for all default options. *)
handle = FSCMSSMOpenHandle[
  fsSettings -> { precisionGoal -> 1.*^-4 },
  fsSMParameters -> { Mt -> 173.3 },
  fsModelParameters -> {
      m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
];

FSCMSSMCalculateSpectrum[handle];

unitarity = FSCMSSMCalculateUnitarity[handle];

Print[Abs[(FlexibleSUSYUnitarity`RenormalizationScale /. unitarity[[1,2]]) - 866.8013526495104]/866.8013526495104];
TestCloseRel[FlexibleSUSYUnitarity`RenormalizationScale /. unitarity[[1,2]], 866.8013526495104, 10^-16];
Print[Abs[(FlexibleSUSYUnitarity`MaxAbsReEigen /. unitarity[[1,2]]) - 0.02973890765470515]/0.02973890765470515];
TestCloseRel[FlexibleSUSYUnitarity`MaxAbsReEigen /. unitarity[[1,2]], 0.02973890765470515, 7.*10^-16];
