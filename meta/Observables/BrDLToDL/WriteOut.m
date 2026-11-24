WriteOut`WriteObservable[
   blockName_,
   obs:FlexibleSUSYObservable`BrDLToDL[{_[qI_], _[lI_]} -> {_[qO_], _}, contr_, loopN_]
] :=
Switch[blockName,
   "FWCOEF"|"IMFWCOEF",
      Module[{qd, lep, type, contrName},
         qd = StringJoin[{qO, qI} /. {0 -> "01", 1 -> "03", 2 -> "05"}];
         lep = StringJoin[{lI, lI} /. {0 -> "11", 1 -> "13", 2 -> "15"}];
         type = blockName /. {"FWCOEF" -> "Re", "IMFWCOEF" -> "Im"};
         contrName = StringRiffle[SymbolName /@ Flatten@{contr}, "+"];
         StringReplace[
            {
               "qd, 4422, 0, 0, 1, obs, \"C7 1 loop\"",
               "qd, 4322, 0, 0, 1, obs, \"C7' 1 loop\"",
               "qd, 6421, 0, 0, 1, obs, \"C8 1 loop\"",
               "qd, 6321, 0, 0, 1, obs, \"C8' 1 loop\"",
               "qdlep, 4133, 0, 0, 1, obs, \"C9 contr loopN loop\"",
               "qdlep, 4233, 0, 0, 1, obs, \"C9' contr loopN loop\"",
               "qdlep, 4137, 0, 0, 1, obs, \"C10 contr loopN loop\"",
               "qdlep, 4237, 0, 0, 1, obs, \"C10' contr loopN loop\"",
               "qdlep, 3230, 0, 0, 1, obs, \"CS contr loopN loop\"",
               "qdlep, 3130, 0, 0, 1, obs, \"CS' contr loopN loop\"",
               "qdlep, 3234, 0, 0, 1, obs, \"CP contr loopN loop\"",
               "qdlep, 3134, 0, 0, 1, obs, \"CP' contr loopN loop\""
            },
            {
               "qd" -> qd, "lep" -> lep,
               "contr" -> contrName, "loopN" -> ToString@loopN,
               "obs" -> type <> "(observables." <> Observables`GetObservableName@obs <> ")"
            }
         ]
      ]
];
