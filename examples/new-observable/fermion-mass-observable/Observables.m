DefineObservable[
   FlexibleSUSYObservable`ExampleFermionMass[fermion_@gen_],
   GetObservableType        -> CConversion`ArrayType[CConversion`complexScalarCType, 2],
   GetObservableName        -> "example_fermion_gen_mass",
   GetObservableDescription -> "example mass of fermion_gen",
   CalculateObservable      -> "ex_fermion_mass(gen, MODEL, qedqcd)",
   GetObservablePrototype   -> "ex_fermion_mass(int idx, const eigenstates& model, const softsusy::QedQcd& qedqcd)",
   GetObservableNamespace   -> "example_fermion_mass"
];
