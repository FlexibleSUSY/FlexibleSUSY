WriteOut`WriteObservable[
   "FlexibleSUSYLowEnergy",
   obs:FlexibleSUSYObservable`BrLToLGamma[pIn_[idxIn_] -> {pOut_[idxOut_], V_}]
] := "observables." <> Observables`GetObservableName@obs;

WriteOut`WriteObservable[
   "FlexibleSUSYLowEnergy",
   obs:FlexibleSUSYObservable`BrLToLGamma[pIn_ -> {pOut_, V_}]
] := "observables." <> Observables`GetObservableName@obs;


