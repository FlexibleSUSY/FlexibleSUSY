DefineObservable[
   FlexibleSUSYObservable`BrLToLGamma[pIn_[idxIn_] -> {pOut_[idxOut_], V_}],
   GetObservableName ->
      "pInidxIn_to_pOutidxOut_V",
   GetObservableDescription ->
      "BR(pIn$(idxIn+1) -> pOut$(idxOut+1) V)",
   GetObservableType ->
      CConversion`ScalarType[CConversion`realScalarCType],
   CalculateObservable ->
      "calculate_pIn_to_pOut_V(idxIn, idxOut, model, qedqcd, physical_input)",
   GetObservablePrototype -> "dummy()"
];

DefineObservable[
   FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, V_}],
   GetObservableName ->
      "pIn_to_pOut_V",
   GetObservableDescription ->
      "BR(pIn -> pOut V)",
   GetObservableType ->
      CConversion`ScalarType[CConversion`realScalarCType],
   CalculateObservable ->
      "calculate_pIn_to_pOut_V(model, qedqcd, physical_input)",
   GetObservablePrototype -> "dummy()"
];
