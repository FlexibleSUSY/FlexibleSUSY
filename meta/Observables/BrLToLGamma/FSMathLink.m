Begin@"FSMathLink`Private`";

PutObservable[FlexibleSUSYObservable`BrLToLGamma[p1_[idx1_Integer]->{p2_[idx2_Integer], V_}], type_, link_String, heads_:{}] /; V === TreeMasses`GetPhoton[] := "
MLPutFunction(link, \"Rule\", 2);
MLPutFunction(link, \"FlexibleSUSYObservable`BrLToLGamma\", 1);
MLPutFunction(link, \"Rule\", 2);
MLPutFunction(link, " <> ObsToStr[p1] <> ", 1);
MLPutInteger(link, " <> ObsToStr[idx1] <> ");
MLPutFunction(link, \"List\", 2);
MLPutFunction(link, " <> ObsToStr[p2] <> ", 1);
MLPutInteger(link, " <> ObsToStr[idx2] <> ");
MLPutSymbol(link, " <> ObsToStr[V] <> ");
MLPutReal(link, OBSERVABLE(" <> ToString[p1] <> ToString[idx1] <> "_to_" <> ToString[p2] <> ToString[idx2] <> "_" <> ToString[V] <> "));"

End[];
