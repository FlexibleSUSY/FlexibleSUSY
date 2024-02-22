Observables`DefineObservable[
FlexibleSUSYObservable`BrDLToDL[{qd_[qI_], lep_[lI_]} -> {qd_[qO_], lep_[lI_]}, contr_, loopN_],
   GetObservableType        -> {13},
   GetObservableName        -> "qdqIleplI_to_qdqOleplI_contr_loopN",
   GetObservableDescription -> "qd$(qI+1) lep$(lI+1) -> qd$(qO+1) lep$(lI+1) contr at loopN loops",
   CalculateObservable      -> "calculate_qdlep_qdlep_contr_loopN(qI, qO, lI, model, qedqcd)",
   GetObservablePrototype   -> "calculate_qdlep_qdlep_contr_loopN(int qi, int qo, int li, auto model, auto qedqcd)"
];
