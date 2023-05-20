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

BeginPackage["NPointFunctions`", {"SARAH`"}];

CreateFAModelFile::usage=
"@brief Create the FeynArts model file using SARAH`.`MakeFeynArts[].
@param sarahInputDirectories the SARAH input directories
@param sarahOutputDirectory the SARAH output directory
@param sarahModelName the model name used by SARAH
@param eigenstates the option to forward to SARAH`.`MakeFeynArts[]
@note Effectively Private function hence there are no checks for input";

Begin["`Private`"];

CreateFAModelFile[sarahInputDirs_, sarahOutputDir_, sarahModelName_, 
   eigenstates_] :=
(
   SARAH`SARAH@SARAH`InputDirectories = sarahInputDirs;
   SARAH`SARAH@SARAH`OutputDirectory = sarahOutputDir;
   SARAH`Start@sarahModelName;
   SA`CurrentStates = eigenstates; 
   SARAH`InitVertexCalculation[eigenstates, False];
   SARAH`partDefinition = SARAH`ParticleDefinitions@eigenstates;
   SARAH`Particles@SARAH`Current = SARAH`Particles@eigenstates;
   SARAH`ReadVertexList[eigenstates, False, False, True];
   SARAH`MakeCouplingLists;
   SARAH`MakeFeynArts[SARAH`Eigenstates -> eigenstates];
);
SetAttributes[CreateFAModelFile,{Protected, Locked}];

End[];
EndPackage[];
