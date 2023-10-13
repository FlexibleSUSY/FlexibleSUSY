WriteOut`WriteObservable[
   blockName_,
   obs:FlexibleSUSYObservable`MyNewObservable[_]
] :=
Switch[blockName,
   "FlexibleSUSYLowEnergy",
      "Re(OBSERVABLES." <> Observables`GetObservableName@obs <> "(0))"
];
