DefineObservable[
   FlexibleSUSYObservable`LToLConversion[pI_@iI_ -> _@iO_, nucl_, contr_, loopN_],
   GetObservableName ->
      "pIiIpIiOinnucl_contr_loopNloop",
   GetObservableDescription ->
      "pI(iI)->pI(iO) in nucl for contr at loopN loop",
   GetObservableType ->
      CConversion`ArrayType[CConversion`complexScalarCType, 13],
   CalculateObservable ->
      "calculate_pIpI_forcontr_loopNloop(iI, iO, context::Nucleus::nucl, MODEL, ltolconversion_settings, qedqcd)",
   GetObservablePrototype ->
      "calculate_pIpI_forcontr_loopNloop(int in, int out, const context::Nucleus n, const eigenstates& model, const LToLConversion_settings& parameters, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace -> "ltolconversion"
];
