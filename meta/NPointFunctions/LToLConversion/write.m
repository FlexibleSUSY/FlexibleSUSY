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

BeginPackage@"LToLConversion`";Quiet[

LToLConversion`write::usage = "
@brief Returns an expression to be printed in the C++ output for observable.
@param block A string with a name of Les Houshes block.
@param obs An observable to parse.
@param idx An index, corresponding to the block.
@param comment A string comment.
@returns A C++ code for observable output.";

];Begin@"`Private`";

getFLHA::usage = "
@brief Returns information of Wilson coefficients, calculated by this observable
       in the format specified by [arxiv:1008.0762].
@param An observable.
@returns List of List of Strings which are: fermions in basis element, operators
         in basis element, orders of perturbative expansion, type of contribution,
         description.
@note We assume that there is a normal ordering of leptons."
getFLHA@`type`observable :=
Module[{rules = {0->"11", 1->"13", 2->"15"}, leptons,
      quarksU = "0202", quarksD = "0101"},
   leptons = StringJoin[{iOut, iIn} /. rules];
   {
      {leptons, "3122", #2, #3, #4, "D_L"},
      {leptons, "3122", #2, #3, #4, "D_R"},
      {#1, "3131", #2, #3, #4, "S_LL "<>#5},
      {#1, "3132", #2, #3, #4, "S_LR "<>#5},
      {#1, "3231", #2, #3, #4, "S_RL "<>#5},
      {#1, "3232", #2, #3, #4, "S_RR "<>#5},
      {#1, "4141", #2, #3, #4, "V_LL "<>#5},
      {#1, "4142", #2, #3, #4, "V_LR "<>#5},
      {#1, "4241", #2, #3, #4, "V_RL "<>#5},
      {#1, "4242", #2, #3, #4, "V_RR "<>#5},
      {#1, "4343", #2, #3, #4, "T_LL "<>#5},
      {#1, "4444", #2, #3, #4, "T_RR "<>#5}
   } & [leptons<>quarksU, "0", "0", "2", CConversion`ToValidCSymbolString@con]];
getFLHA // Utils`MakeUnknownInputDefinition;
getFLHA ~ SetAttributes ~ {Protected, Locked};

FlexibleSUSYObservable`LToLConversion::errBlock = "
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
         LToLConversion`getFLHA@obs],
   _,
      Utils`AssertOrQuit[_, FlexibleSUSYObservable`LToLConversion::errBlock]];
write // Utils`MakeUnknownInputDefinition;
write ~ SetAttributes ~ {Protected, Locked};

End[];EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "LToLConversion`"];
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "LToLConversion`"];
Protect@$Packages;
