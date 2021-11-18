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

Off[BrLTo3L`write::shdw];
BeginPackage["BrLTo3L`"];

BrLTo3L`write::usage = "
@brief Write observable in the C++ Les Houshes block.
       Please, don't change write interface.
@param block A Les Houshes block name.
@param obs An observable to parse.
@param idx An index, corresponding to the block.
@param comment A comment.";

Begin["`Private`"];

flha[{nO_, o1_, nI_}, {nA_, o2_, nA_}, a_, as_, c_, name_] :=
If[nO < nA,
   {nO<>nI<>nA<>nA, o1<>o2, a, as, c, name},
   {nA<>nA<>nO<>nI, o2<>o1, a, as, c, name}];
flha // Utils`MakeUnknownInputDefinition;
flha ~ SetAttributes ~ {Protected, Locked};

getFLHA::usage = "
@brief Returns information of Wilson coefficients in format by [1008.0762].
@param An observable.
@returns List of List of Strings in FLHA format."
getFLHA@`type`observable :=
Module[{fields, args},
   fields = {0->"11", 1->"13", 2->"15"};
   args = Sequence[ToString@loopN, SymbolName@proc, "0", "0", "2"];
   Quiet[{{nO<>nI, "4322", ##3, "D_L"},
          {nO<>nI, "4422", ##3, "D_R"},
          flha[{nO, "31", nI}, {nA, "31", nA}, ##3, "S_LL_"<>#1<>"loop "<>#2],
          flha[{nO, "31", nI}, {nA, "32", nA}, ##3, "S_LR_"<>#1<>"loop "<>#2],
          flha[{nO, "32", nI}, {nA, "31", nA}, ##3, "S_RL_"<>#1<>"loop "<>#2],
          flha[{nO, "32", nI}, {nA, "32", nA}, ##3, "S_RR_"<>#1<>"loop "<>#2],
          flha[{nO, "41", nI}, {nA, "41", nA}, ##3, "V_LL_"<>#1<>"loop "<>#2],
          flha[{nO, "41", nI}, {nA, "42", nA}, ##3, "V_LR_"<>#1<>"loop "<>#2],
          flha[{nO, "42", nI}, {nA, "41", nA}, ##3, "V_RL_"<>#1<>"loop "<>#2],
          flha[{nO, "42", nI}, {nA, "42", nA}, ##3, "V_RR_"<>#1<>"loop "<>#2],
          flha[{nO, "43", nI}, {nA, "43", nA}, ##3, "T_LL_"<>#1<>"loop "<>#2],
          flha[{nO, "44", nI}, {nA, "44", nA}, ##3, "T_RR_"<>#1<>"loop "<>#2]}&@args,
         StringJoin::string] /. fields];
getFLHA // Utils`MakeUnknownInputDefinition;
getFLHA ~ SetAttributes ~ {Protected, Locked};

write[block:_String, obs:`type`observable, idx:_Integer, comment:_String] :=
Module[{configured},
   configured = And[FlexibleSUSY`FSFeynArtsAvailable,
                    FlexibleSUSY`FSFormCalcAvailable];
   Utils`AssertOrQuit[configured,
      FlexibleSUSYObservable`BrLTo3L::errConfigured,
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
            FlexibleSUSYObservable`BrLTo3L::errBlock,
            obs,
            block,
            FlexibleSUSY`FSOutputDir,
            {"FlexibleSUSYLowEnergy", "FWCOEF", "IMFWCOEF"}
         ];
   ]
];

FlexibleSUSYObservable`BrLTo3L::errConfigured = "
FeynArts and/or FormCalc are not installed.
Solutions:
1) Remove an observable:\n      `1`
   from the file:\n      `2`/FlexibleSUSY.m
2) Install FeynArts and FormCalc.";

FlexibleSUSYObservable`BrLTo3L::errBlock = "
An observable:\n   `1`
from the file:\n   `3`/FlexibleSUSY.m
is mentioned in Les Houches block:\n   `2`,
but should only be used inside block(s):\n   `4`";

write // Utils`MakeUnknownInputDefinition;
write // Protected;

End[];
Block[{$ContextPath}, EndPackage[]];
