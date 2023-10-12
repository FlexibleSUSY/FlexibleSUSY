Begin@"FSMathLink`Private`";
Module[{args = Sequence[lep_@nI_ -> {lep_@nO_, lep_@nA_, SARAH`bar@lep_@nA_}, proc_, loopN_],
      obs = FlexibleSUSYObservable`BrLTo3L,
      cxx = CConversion`ToValidCSymbolString},
      PutObservable[obs@args, type_, link_String, heads_:{}] := StringJoin[
      "MLPutFunction(link, \"Rule\", 2);",
         "MLPutFunction(link, ", ObsToStr@obs, ", 3);",
            "MLPutFunction(link, \"Rule\", 2);",
               "MLPutFunction(link, \"", cxx@lep, "\", 1);",
                  "MLPutInteger(link, ", ObsToStr@nI, ");",
               "MLPutFunction(link, \"List\", 3);",
                  "MLPutFunction(link, \"", cxx@lep, "\", 1);",
                     "MLPutInteger(link, ", ObsToStr@nO, ");",
                  "MLPutFunction(link, \"", cxx@lep, "\", 1);",
                     "MLPutInteger(link, ", ObsToStr@nA, ");",
                  "MLPutFunction(link, \"SARAH`bar\", 1);",
                     "MLPutFunction(link, \"", cxx@lep, "\", 1);",
                        "MLPutInteger(link, ", ObsToStr@nA, ");",
            "MLPutSymbol(link, \"", SymbolName@proc, "\");",
            "MLPutInteger(link, ", ObsToStr@loopN, ");",
         "MLPutReal(link, Re(OBSERVABLE(", GetObservableName@obs[
            lep@nI -> {lep@nO, lep@nA, SARAH`bar@lep@nA},
            proc,
            loopN], ")(0)));"];];
End[];
