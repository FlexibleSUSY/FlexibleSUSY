DefineObservable[
   FlexibleSUSYObservable`ExampleConstantObservable[num_],
   GetObservableType        -> CConversion`ArrayType[CConversion`complexScalarCType, 1],
   GetObservableName        -> "example_constant_observable_num",
   GetObservableDescription -> "example observable returns num",
   CalculateObservable      -> "calculate_example_constant_observable(num)",
   GetObservablePrototype   -> "calculate_example_constant_observable(double con)",
   GetObservableNamespace   -> "example_constant_observable"
];
