WriteOut`WriteObservable[
   "FlexibleSUSYLowEnergy",
   obs:FlexibleSUSYObservable`ExampleFermionMass[_]
] := "Re(observables." <> Observables`GetObservableName@obs <> "(0))";

WriteOut`WriteObservable[
   "ExampleLeptonMass",
   obs:FlexibleSUSYObservable`ExampleFermionMass[_]
] := "Re(observables." <> Observables`GetObservableName@obs <> "(1))";
