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
Needs["SARAH`"];
Needs["FlexibleSUSY`", "FlexibleSUSY.m"];

ResetFlags[] := (
    FlexibleSUSY`UseHiggs2LoopNMSSM = False;
    FlexibleSUSY`UseHiggs3LoopNMSSM = False;
    SARAH`UseHiggs2LoopMSSM = False;
    FlexibleSUSY`UseHiggs3LoopMSSM = False;
    FlexibleSUSY`UseMSSMYukawa2Loop = False;
    FlexibleSUSY`UseMSSMAlphaS2Loop = False;
    FlexibleSUSY`UseMSSM3LoopRGEs = False;
);

Print["testing that UseHiggs3LoopNMSSM implies the 2-loop corrections ..."];

ResetFlags[];
FlexibleSUSY`UseHiggs3LoopNMSSM = True;
FlexibleSUSY`Private`FSCheckFlags[];

TestEquality[FlexibleSUSY`UseHiggs2LoopNMSSM, True];
TestEquality[FlexibleSUSY`UseMSSMYukawa2Loop, True];
TestEquality[FlexibleSUSY`UseMSSMAlphaS2Loop, True];
TestEquality[FlexibleSUSY`UseMSSM3LoopRGEs, True];

Print["testing that UseHiggs3LoopMSSM implies the 2-loop corrections ..."];

ResetFlags[];
FlexibleSUSY`UseHiggs3LoopMSSM = True;
FlexibleSUSY`Private`FSCheckFlags[];

TestEquality[SARAH`UseHiggs2LoopMSSM, True];
TestEquality[FlexibleSUSY`UseMSSMYukawa2Loop, True];
TestEquality[FlexibleSUSY`UseMSSMAlphaS2Loop, True];
TestEquality[FlexibleSUSY`UseMSSM3LoopRGEs, True];

PrintTestSummary[];
