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

Off[LToLConversion`write::shdw];
BeginPackage["LToLConversion`"];

LToLConversion`write::usage = "
@brief Write observable in the C++ Les Houshes block.
       Please, don't change write interface.
@param block A Les Houshes block name.
@param obs An observable to parse.
@param idx An index, corresponding to the block.
@param comment A comment.";

Begin["`Private`"];

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
      {leptons, "4322", #2, #3, #4, "D_L"},
      {leptons, "4422", #2, #3, #4, "D_R"},
      {#1, "3131", #2, #3, #4, "S_LL_"<>#5<>"loop "<>#6},
      {#1, "3132", #2, #3, #4, "S_LR_"<>#5<>"loop "<>#6},
      {#1, "3231", #2, #3, #4, "S_RL_"<>#5<>"loop "<>#6},
      {#1, "3232", #2, #3, #4, "S_RR_"<>#5<>"loop "<>#6},
      {#1, "4141", #2, #3, #4, "V_LL_"<>#5<>"loop "<>#6},
      {#1, "4142", #2, #3, #4, "V_LR_"<>#5<>"loop "<>#6},
      {#1, "4241", #2, #3, #4, "V_RL_"<>#5<>"loop "<>#6},
      {#1, "4242", #2, #3, #4, "V_RR_"<>#5<>"loop "<>#6},
      {#1, "4343", #2, #3, #4, "T_LL_"<>#5<>"loop "<>#6},
      {#1, "4444", #2, #3, #4, "T_RR_"<>#5<>"loop "<>#6}
   } & [leptons<>quarksU, "0", "0", "2", ToString@loopN, SymbolName@con]];
getFLHA // Utils`MakeUnknownInputDefinition;
getFLHA ~ SetAttributes ~ {Protected, Locked};

write[block:_String, obs:`type`observable, idx:_Integer, comment:_String] :=
Module[{configured},
   configured = And[FlexibleSUSY`FSFeynArtsAvailable,
                    FlexibleSUSY`FSFormCalcAvailable];
   Utils`AssertOrQuit[configured,
      FlexibleSUSYObservable`LToLConversion::errConfigured,
      obs,
      FlexibleSUSY`FSOutputDir
   ];
   Switch[block,
      "FlexibleSUSYLowEnergy",
         WriteOut`Private`WriteSLHABlockEntry[block,
            {"Re(OBSERVABLES."<>Observables`GetObservableName@obs<>"(0))", idx},
            Observables`GetObservableDescription@obs
         ],
      "FWCOEF"|"IMFWCOEF",
         WriteOut`Private`WriteSLHABlockEntry[block,
            {"OBSERVABLES."<>Observables`GetObservableName@obs,
               Observables`GetObservableType@obs},
            getFLHA@obs
         ],
      _,
         Utils`AssertOrQuit[_,
            FlexibleSUSYObservable`LToLConversion::errBlock,
            obs,
            block,
            FlexibleSUSY`FSOutputDir,
            {"FlexibleSUSYLowEnergy", "FWCOEF", "IMFWCOEF"}
         ];
   ]
];

FlexibleSUSYObservable`LToLConversion::errConfigured = "
FeynArts and/or FormCalc are not installed.
Solutions:
1) Remove an observable:\n      `1`
   from the file:\n      `2`/FlexibleSUSY.m
2) Install FeynArts and FormCalc.";

FlexibleSUSYObservable`LToLConversion::errBlock = "
An observable:\n   `1`
from the file:\n   `3`/FlexibleSUSY.m
is mentioned in Les Houches block:\n   `2`,
but should only be used inside block(s):\n   `4`";

write // Utils`MakeUnknownInputDefinition;
write // Protected;

End[];
Block[{$ContextPath}, EndPackage[]];
