WriteOut`WriteObservable[
   blockName_,
   obs:FlexibleSUSYObservable`BrLTo3L[_@iI_ -> {_@iO_, _@iA_, _}, contr_, loopN_]
] :=
Switch[blockName,
   "FlexibleSUSYLowEnergy",
      "Re(OBSERVABLES." <> Observables`GetObservableName@obs <> "(0))",
   "FWCOEF"|"IMFWCOEF",
      Module[{c1, c2, type},
         c1 = StringJoin[{iO, iI} /. {0 -> "11", 1 -> "13", 2 -> "15"}];
         c2 = StringJoin[{iA, iA} /. {0 -> "11", 1 -> "13", 2 -> "15"}];
         type = blockName /. {"FWCOEF"->"Re", "IMFWCOEF"->"Im"};
         StringReplace[
            {
               "c1, 4322, 0, 0, 2, obs, \"D_L 1 loop\"",
               "c1, 4422, 0, 0, 2, obs, \"D_R 1 loop\"",
               "c1c2, 3131, 0, 0, 2, obs, \"S_LL contr loopN loop\"",
               "c1c2, 3132, 0, 0, 2, obs, \"S_LR contr loopN loop\"",
               "c1c2, 3231, 0, 0, 2, obs, \"S_RL contr loopN loop\"",
               "c1c2, 3232, 0, 0, 2, obs, \"S_RR contr loopN loop\"",
               "c1c2, 4141, 0, 0, 2, obs, \"V_LL contr loopN loop\"",
               "c1c2, 4142, 0, 0, 2, obs, \"V_LR contr loopN loop\"",
               "c1c2, 4241, 0, 0, 2, obs, \"V_RL contr loopN loop\"",
               "c1c2, 4242, 0, 0, 2, obs, \"V_RR contr loopN loop\"",
               "c1c2, 4343, 0, 0, 2, obs, \"T_LL contr loopN loop\"",
               "c1c2, 4444, 0, 0, 2, obs, \"T_RR contr loopN loop\""
            },
            {
               "c1" -> c1, "c2" -> c2,
               "contr" -> SymbolName@contr, "loopN" -> ToString@loopN,
               "obs" -> type<>"(OBSERVABLES."<>Observables`GetObservableName@obs<>")"
            }
         ]
      ]
];
