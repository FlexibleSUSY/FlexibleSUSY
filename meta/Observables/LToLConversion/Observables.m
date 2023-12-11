DefineObservable[
   FlexibleSUSYObservable`LToLConversion[lep_, iI_ -> iO_, nucl_, contr_, loopN_],
   GetObservableName ->
      "lepiIlepiOinnucl_contr_loopNloop",
   GetObservableDescription ->
      "lep(iI)->lep(iO) in nucl for contr at loopN loop",
   GetObservableType ->
      CConversion`ArrayType[CConversion`complexScalarCType, 13],
   CalculateObservable ->
      "calculate_leplep_forcontr_loopNloop(iI-1, iO-1, context::Nucleus::nucl, MODEL, ltolconversion_settings, qedqcd)",
   GetObservablePrototype ->
      "calculate_leplep_forcontr_loopNloop(int in, int out, const context::Nucleus n, const eigenstates& model, const LToLConversion_settings& parameters, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace -> "ltolconversion",
   GetObservableFileName -> "ltolconversion"
];
