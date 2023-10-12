DefineObservable[
   FlexibleSUSYObservable`BrLTo3L[pI_@iI_ -> {_@iO_, _@iA_, _}, contr_, loopN_],
   GetObservableName ->
      "pIiItoiOiAbariA_contr_loopNloop",
   GetObservableDescription ->
      "pI(iI->iOiAbariA) for contr at loopN loop",
   GetObservableType ->
      CConversion`ArrayType[CConversion`complexScalarCType, 13],
   CalculateObservable ->
      "calculate_pI_to_pIpIbarpI_for_contr_loopNloop(iI, iO, iA, MODEL, qedqcd)",
   GetObservablePrototype ->
      "calculate_pI_to_pIpIbarpI_for_contr_loopNloop(int nI, int nO, int nA, const eigenstates& model, const softsusy::QedQcd& qedqcd)",
   Context -> "lto3l"
];
