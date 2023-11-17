DefineObservable[
   FlexibleSUSYObservable`BrLTo3L[lep_, iI_ -> {iJ_, iK_, iK_}, contr_, loopN_],
   GetObservableName ->
      "lepiItoiJiKbariK_contr_loopNloop",
   GetObservableDescription ->
      "lep(iI->iJiKbariK) for contr at loopN loop",
   GetObservableType ->
      CConversion`ArrayType[CConversion`complexScalarCType, 13],
   CalculateObservable ->
      "calculate_lep_to_leplepbarlep_for_contr_loopNloop(iI, iJ, iK, MODEL, qedqcd)",
   GetObservablePrototype ->
      "calculate_lep_to_leplepbarlep_for_contr_loopNloop(int nI, int nO, int nA, const eigenstates& model, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace -> "lto3l"
];
