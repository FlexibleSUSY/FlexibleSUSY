Observables`DefineObservable[
   FlexibleSUSYObservable`HiggsTo2Gluons[higgs_ -> {gl_, gl_}, contr_],
   GetObservableType        -> CConversion`ArrayType[CConversion`complexScalarCType, 1],
   GetObservableName        -> "higgs_glgl_contr",
   GetObservableDescription -> "higgs -> glgl for contr",
   CalculateObservable      -> "calculate_higgstoglgl_contr(MODEL, qedqcd)",
   GetObservablePrototype   -> "calculate_higgstoglgl_contr(const eigenstates& model, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace   -> "h_to_gg",
   GetObservableFileName    -> "h_to_gg"
];
