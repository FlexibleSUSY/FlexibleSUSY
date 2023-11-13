WriteOut`WriteObservable[
   blockName_,
   obs:FlexibleSUSYObservable`LToLConversion[_, iI_ -> iO_, _, contr_, loopN_]
] :=
Switch[blockName,
   "FlexibleSUSYLowEnergy",
      "Re(OBSERVABLES." <> Observables`GetObservableName@obs <> "(0))",
   "FWCOEF"|"IMFWCOEF",
      Module[{lep, type},
         lep = StringJoin[{iO, iI} /. {0 -> "11", 1 -> "13", 2 -> "15"}];
         type = blockName /. {"FWCOEF"->"Re", "IMFWCOEF"->"Im"};
         StringReplace[
            {
               "lep, 4322, 0, 0, 2, obs, \"D_L 1 loop\"",
               "lep, 4422, 0, 0, 2, obs, \"D_R 1 loop\"",
               "lepqu, 3131, 0, 0, 2, obs, \"S_LL contr loopN loop\"",
               "lepqu, 3132, 0, 0, 2, obs, \"S_LR contr loopN loop\"",
               "lepqu, 3231, 0, 0, 2, obs, \"S_RL contr loopN loop\"",
               "lepqu, 3232, 0, 0, 2, obs, \"S_RR contr loopN loop\"",
               "lepqu, 4141, 0, 0, 2, obs, \"V_LL contr loopN loop\"",
               "lepqu, 4142, 0, 0, 2, obs, \"V_LR contr loopN loop\"",
               "lepqu, 4241, 0, 0, 2, obs, \"V_RL contr loopN loop\"",
               "lepqu, 4242, 0, 0, 2, obs, \"V_RR contr loopN loop\"",
               "lepqu, 4343, 0, 0, 2, obs, \"T_LL contr loopN loop\"",
               "lepqu, 4444, 0, 0, 2, obs, \"T_RR contr loopN loop\""
            },
            {
               "lep" -> lep, "qu" -> "0202", "qd" -> "0101",
               "contr" -> SymbolName@contr, "loopN" -> ToString@loopN,
               "obs" -> type<>"(OBSERVABLES."<>Observables`GetObservableName@obs<>")"
            }
         ]
      ]
];
