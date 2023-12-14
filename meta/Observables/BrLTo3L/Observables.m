DefineObservable[
   FlexibleSUSYObservable`BrLTo3L[lep_, iI_ -> {iJ_, iK_, iK_}, contr_, loopN_],
   GetObservableType        -> {13},
   GetObservableName        -> "lepiItoiJiKbariK_contr_loopNloop",
   GetObservableDescription -> "lep(iI)->lep(iJ)lep(iK)barlep(iK) for contr at loopN loop",
   CalculateObservable ->
      "calculate_lep_to_leplepbarlep_for_contr_loopNloop($(iI-1), $(iJ-1), $(iK-1), model, qedqcd)",
   GetObservablePrototype ->
      "calculate_lep_to_leplepbarlep_for_contr_loopNloop(int nI, int nO, int nA, auto model, auto qedqcd)"
];
