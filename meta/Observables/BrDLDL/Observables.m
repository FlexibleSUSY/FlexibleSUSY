Observables`DefineObservable[
FlexibleSUSYObservable`BrDLDL[{qd_[qI_], lep_[lI_]} -> {qd_[qO_], lep_[lI_]}, contr_, loopN_],
   GetObservableName        -> "qdqIleplI_to_qdqOleplI_contr_loopN",
   GetObservableDescription -> "qd(qI) lep(lI) -> qd(qO) lep(lI) contr at loopN loops",
   GetObservableType        -> CConversion`ArrayType[CConversion`complexScalarCType, 13],
   CalculateObservable      -> "calculate_qdlep_qdlep_contr_loopN(qI, qO, lI, MODEL, qedqcd)",
   GetObservablePrototype   -> "calculate_qdlep_qdlep_contr_loopN(int qi, int qo, int li, const eigenstates& model, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace   -> "dltodl"
];
