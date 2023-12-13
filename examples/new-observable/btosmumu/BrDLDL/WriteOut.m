WriteOut`WriteObservable[
   "FWCOEF",
   obs:FlexibleSUSYObservable`BrDLDL[{qd_[qI_], lep_[lI_]} -> {_[qO_], _}, __]
] :=
StringReplace[
   {
      "qd,    4322, exp, num_value, \"D_L\"",
      "qd,    4422, exp, num_value, \"D_R\"",
      "qdlep, 3131, exp, num_value, \"S_LL\"",
      "qdlep, 3132, exp, num_value, \"S_LR\"",
      "qdlep, 3231, exp, num_value, \"S_RL\"",
      "qdlep, 3232, exp, num_value, \"S_RR\"",
      "qdlep, 4141, exp, num_value, \"V_LL\"",
      "qdlep, 4142, exp, num_value, \"V_LR\"",
      "qdlep, 4241, exp, num_value, \"V_RL\"",
      "qdlep, 4242, exp, num_value, \"V_RR\"",
      "qdlep, 4343, exp, num_value, \"T_LL\"",
      "qdlep, 4444, exp, num_value, \"T_RR\""
   },
   {
      "num_value" -> "Re(observables." <> Observables`GetObservableName@obs <> ")",
      "exp"       -> "9, 9, 2",
      "qd"        -> StringJoin[{qO, qI} /. {0 -> "01", 1 -> "03", 2 -> "05"}],
      "lep"       -> StringJoin[{lI, lI} /. {0 -> "11", 1 -> "13", 2 -> "15"}]
   }
];
