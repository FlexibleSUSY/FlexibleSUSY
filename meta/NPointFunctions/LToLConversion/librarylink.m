Begin@"FSMathLink`Private`";
Module[{args = Sequence[in_@inN_ -> _@outN_, nucl_, proc_, loopN_],
      obs = FlexibleSUSYObservable`LToLConversion,
      cxx = CConversion`ToValidCSymbolString},
      PutObservable[obs@args, type_, link_String, heads_:{}] := StringJoin[
      "MLPutFunction(link, \"Rule\", 2);",
         "MLPutFunction(link, ", ObsToStr@obs, ", 4);",
            "MLPutFunction(link, \"Rule\", 2);",
               "MLPutFunction(link, \"", cxx@in, "\", 1);",
                  "MLPutInteger(link, ", ObsToStr@inN, ");",
               "MLPutFunction(link, \"", cxx@in, "\", 1);",
                  "MLPutInteger(link, ", ObsToStr@outN, ");",
            "MLPutSymbol(link, \"", SymbolName@nucl, "\");",
            "MLPutSymbol(link, \"", SymbolName@proc, "\");",
            "MLPutInteger(link, ", ObsToStr@loopN, ");",
         "MLPutReal(link, Re(OBSERVABLE(", GetObservableName@obs[
            in@inN -> in@outN, nucl, proc, loopN], ")(0)));"];];
End[];
