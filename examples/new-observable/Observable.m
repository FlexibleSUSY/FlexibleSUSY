DefineObservable[
   FlexibleSUSYObservable`MyNewObservable[num_],
   GetObservableName ->
      "mynewobservable_num",
   GetObservableDescription ->
      "mynewobservable description",
   GetObservableType ->
      CConversion`ArrayType[CConversion`complexScalarCType, 1],
   CalculateObservable ->
      "calculate_mynewobservable(MODEL, qedqcd)",
   GetObservablePrototype ->
      "calculate_mynewobservable(const eigenstates& model, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace -> "mynewobservable"
];
