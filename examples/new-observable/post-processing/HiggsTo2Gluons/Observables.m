Observables`DefineObservable[
   FlexibleSUSYObservable`HiggsTo2Gluons[higgs_ -> {gl_, gl_}, contr_],
   GetObservableType        -> {1},
   GetObservableName        -> "higgs_glgl_contr",
   GetObservableDescription -> "higgs -> gl gl for contr",
   CalculateObservable      -> "calculate_higgstoglgl_contr(model, qedqcd)",
   GetObservablePrototype   -> "calculate_higgstoglgl_contr(auto model, auto qedqcd)",
   GetObservableNamespace   -> "h_to_gg",
   GetObservableFileName    -> "h_to_gg"
];
