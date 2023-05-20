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

Get["models/SM/SM_librarylink.m"];
Get["models/SMEFTHiggs/SMEFTHiggs_librarylink.m"];

hSMEFTHiggs = FSSMEFTHiggsOpenHandle[betaFunctionLoopOrder -> 3, Mt -> 173.34, LambdaIN -> 0.2, Qin -> 2000];

TestEquality[LambdaIN              /. FSSMEFTHiggsGetInputParameters[hSMEFTHiggs]  , 0.2];
TestEquality[Mt                    /. FSSMEFTHiggsGetSMInputParameters[hSMEFTHiggs], 173.34];
TestEquality[betaFunctionLoopOrder /. FSSMEFTHiggsGetSettings[hSMEFTHiggs]         , 3];

TestEquality[FSSMGetInputParameters[hSMEFTHiggs]  , $Failed];
TestEquality[FSSMGetSMInputParameters[hSMEFTHiggs], $Failed];
TestEquality[FSSMGetSettings[hSMEFTHiggs]         , $Failed];
TestEquality[FSSMCalculateSpectrum[hSMEFTHiggs]   , $Failed];
TestEquality[FSSMCalculateObservables[hSMEFTHiggs], $Failed];

hSM = FSSMOpenHandle[betaFunctionLoopOrder -> 2, Mt -> 173, LambdaIN -> 0.1, Qin -> 1000];

TestEquality[LambdaIN              /. FSSMEFTHiggsGetInputParameters[hSMEFTHiggs]  , 0.2];
TestEquality[Mt                    /. FSSMEFTHiggsGetSMInputParameters[hSMEFTHiggs], 173.34];
TestEquality[betaFunctionLoopOrder /. FSSMEFTHiggsGetSettings[hSMEFTHiggs]         , 3];

TestEquality[LambdaIN              /. FSSMGetInputParameters[hSMEFTHiggs]  , 0.1];
TestEquality[Mt                    /. FSSMGetSMInputParameters[hSMEFTHiggs], 173.];
TestEquality[betaFunctionLoopOrder /. FSSMGetSettings[hSMEFTHiggs]         , 2];

FSSMCloseHandle[hSMEFTHiggs];
FSSMCloseHandle[hSM];

PrintTestSummary[];
