(* d -> d' nu nu: C_L, C_R (last two entries of the observable) in the FLHA
   normalization of the operators (qbar gamma^mu P_{L,R} b)(nubar gamma_mu P_L nu)
   (codes 4141, 4241), i.e. 2 C_L and 2 C_R of flavio *)
WriteOut`WriteObservable[
   blockName_,
   obs:FlexibleSUSYObservable`BrDLToDL[{_[qI_], lep_[lI_]} -> {_[qO_], _}, contr_, loopN_]
] /; !TreeMasses`IsElectricallyCharged[lep] :=
Switch[blockName,
   "FWCOEF"|"IMFWCOEF",
      Module[{qd, nu, type, contrName},
         qd = StringJoin[{qO, qI} /. {0 -> "01", 1 -> "03", 2 -> "05"}];
         nu = StringJoin[{lI, lI} /. {0 -> "12", 1 -> "14", 2 -> "16"}];
         type = blockName /. {"FWCOEF" -> "Re", "IMFWCOEF" -> "Im"};
         contrName = StringRiffle[SymbolName /@ Flatten@{contr}, "+"];
         StringReplace[
            {
               "qdnu, 4141, 0, 0, 1, obs, \"2 CL contr loopN loop\"",
               "qdnu, 4241, 0, 0, 1, obs, \"2 CR contr loopN loop\""
            },
            {
               "qd" -> qd, "nu" -> nu,
               "contr" -> contrName, "loopN" -> ToString@loopN,
               "obs" -> "2.*" <> type <> "(observables." <> Observables`GetObservableName@obs <> ")"
            }
         ]
      ]
];

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
