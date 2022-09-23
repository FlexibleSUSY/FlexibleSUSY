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

Off[FSCMSSM::info];
Off[FSCMSSMCalculateSpectrum::error];
Off[FSCMSSMCalculateSpectrum::warning];
Off[FSCMSSMCalculateObservables::error];
Off[FSCMSSMCalculateObservables::warning];

CalcMh[TB_] :=
    Module[{spec, handle},
           handle = FSCMSSMOpenHandle[
               fsSettings -> {calculateStandardModelMasses -> 1,
                              betaFunctionLoopOrder -> 2},
               fsModelParameters -> {m0 -> 125, m12 -> 500, TanBeta -> TB,
                                     SignMu -> 1, Azero -> 0}];
           spec = FSCMSSMCalculateSpectrum[handle];
           FSCMSSMCloseHandle[handle];
           If[spec === $Failed, 0,
              (Pole[M[hh]] /. spec)[[1]]]
          ];

kernels = LaunchKernels[];
Print["Using ", Length[kernels], " kernels."];

DistributeDefinitions[CalcMh];

r = Sequence[1, 50]
range = Range[r, 1];
resultMap = AbsoluteTiming[Map[CalcMh, range]];
resultParallelMap = AbsoluteTiming[ParallelMap[CalcMh, range, Method -> "CoarsestGrained"]];
resultParallelTab = AbsoluteTiming[ParallelTable[CalcMh[tb], {tb,r}]];

Print["time for sequential runs: ", resultMap[[1]]];
Print["time for parallel map   : ", resultParallelMap[[1]]];
Print["time for parallel table : ", resultParallelTab[[1]]];

TestEquality[resultMap[[2]], resultParallelMap[[2]]];
TestEquality[resultMap[[2]], resultParallelTab[[2]]];

PrintTestSummary[];
