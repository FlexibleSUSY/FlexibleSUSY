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

BeginPackage@"BrLTo3L`";Quiet[

BrLTo3L`write::usage = "
@brief Returns an expression to be printed in the C++ output for observable.
@param block A string with a name of Les Houshes block.
@param obs An observable to parse.
@param idx An index, corresponding to the block.
@param comment A string comment.
@returns A C++ code for observable output.";

];Begin@"`Private`";

flha[{nO_, o1_, nI_}, {nA_, o2_, nA_}, a_, as_, c_, name_] :=
If[nO < nA,
   {nO<>nI<>nA<>nA, o1<>o2, a, as, c, name},
   {nA<>nA<>nO<>nI, o2<>o1, a, as, c, name}];
flha // Utils`MakeUnknownInputDefinition;
flha ~ SetAttributes ~ {Protected, Locked};

getFLHA::usage = "
@brief Returns information of Wilson coefficients in format by 1008.0762.
@param An observable.
@returns List of List of Strings in FLHA format."
getFLHA@`type`observable :=
Module[{fields, args},
   fields = {0->"11", 1->"13", 2->"15"};
   args = Sequence[CConversion`ToValidCSymbolString@proc, "0", "0", "2"];
   Quiet[{{nO<>nI, "4322", ##2, "D_L"},
          {nO<>nI, "4422", ##2, "D_R"},
          flha[{nO, "31", nI}, {nA, "31", nA}, ##2, "S_LL "<>#],
          flha[{nO, "31", nI}, {nA, "32", nA}, ##2, "S_LR "<>#],
          flha[{nO, "32", nI}, {nA, "31", nA}, ##2, "S_RL "<>#],
          flha[{nO, "32", nI}, {nA, "32", nA}, ##2, "S_RR "<>#],
          flha[{nO, "41", nI}, {nA, "41", nA}, ##2, "V_LL "<>#],
          flha[{nO, "41", nI}, {nA, "42", nA}, ##2, "V_LR "<>#],
          flha[{nO, "42", nI}, {nA, "41", nA}, ##2, "V_RL "<>#],
          flha[{nO, "42", nI}, {nA, "42", nA}, ##2, "V_RR "<>#],
          flha[{nO, "43", nI}, {nA, "43", nA}, ##2, "T_LL "<>#],
          flha[{nO, "44", nI}, {nA, "44", nA}, ##2, "T_RR "<>#]}&@args,
         StringJoin::string] /. fields];
getFLHA // Utils`MakeUnknownInputDefinition;
getFLHA ~ SetAttributes ~ {Protected, Locked};

FlexibleSUSYObservable`BrLTo3L::errBlock = "
Observable can be called from: FlexibleSUSYLowEnergy, FWCOEF, IMFWCOEF only.";
write[block:_String, obs:`type`observable, idx:_Integer, comment:_String] :=
Switch[block,
   "FlexibleSUSYLowEnergy",
      WriteOut`Private`WriteSLHABlockEntry[block,
         {"Re(OBSERVABLES."<>Observables`GetObservableName@obs<>"(0))",
            idx},
         Observables`GetObservableDescription@obs],
   "FWCOEF"|"IMFWCOEF",
      WriteOut`Private`WriteSLHABlockEntry[block,
         {"OBSERVABLES."<>Observables`GetObservableName@obs,
            Observables`GetObservableType@obs},
         getFLHA@obs],
   _,
      Utils`AssertOrQuit[_, FlexibleSUSYObservable`BrLTo3L::errBlock]];
write // Utils`MakeUnknownInputDefinition;
write ~ SetAttributes ~ {Protected, Locked};

End[];EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "BrLTo3L`"];
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "BrLTo3L`"];
Protect@$Packages;
