Observables`DefineObservable[
   FlexibleSUSYObservable`ExampleLeptonSE[field_@gen_, contr_],
   GetObservableType        -> CConversion`ArrayType[CConversion`complexScalarCType, 3],
   GetObservableName        -> "example_leptongen_se_contr",
   GetObservableDescription -> "self energy of lepton(gen) for contr",
   CalculateObservable      -> "ex_lepton_se_contr(gen, MODEL)",
   GetObservablePrototype   -> "ex_lepton_se_contr(int idx, const eigenstates& model)",
   GetObservableNamespace   -> "example_lepton_se"
];
