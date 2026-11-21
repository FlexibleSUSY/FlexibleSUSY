WriteOut`WriteObservable[
   "FWCOEF",
   obs:FlexibleSUSYObservable`ExampleLeptonSE[_@gen_, _]
] :=
StringReplace[
   {
      "leptons, 31, 0, 0, 2, num_value, \"left\"", (* P_L *)
      "leptons, 32, 0, 0, 2, num_value, \"right\"" (* P_R *)
   },
   {
      "num_value" -> "Re(observables." <> Observables`GetObservableName@obs <> ")",
      "leptons"   -> Switch[gen, 0, "1111", 1, "1313", 2, "1515"]
   }
];
