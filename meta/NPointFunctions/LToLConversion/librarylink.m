Begin@"FSMathLink`Private`";
Module[{args = LToLConversion`arguments[in@inN, out@outN, nucl, proc],
      obs = FlexibleSUSYObservable`LToLConversion,
      cxx = CConversion`ToValidCSymbolString},
      PutObservable[obs@args, type_, link_String, heads_:{}] := StringJoin[
      "MLPutFunction(link, \"Rule\", 2);",
         "MLPutFunction(link, \"", ToString@obs, "\", 3);",
            "MLPutFunction(link, \"Rule\", 2);",
               "MLPutFunction(link, \"", cxx@in, "\", 1);",
                  "MLPutInteger(link, ", cxx@inN, ");",
               "MLPutFunction(link, \"", cxx@out, "\", 1);",
                  "MLPutInteger(link, ", cxx@outN, ");",
            "MLPutSymbol(link, \"", SymbolName@nucl, "\");",
            "MLPutSymbol(link, \"", SymbolName@proc, "\");",
         "MLPutReal(link, Re(OBSERVABLE(", GetObservableName@obs@##&[
            in@inN -> out@outN, nucl, proc], ")(0)));"];];
End[];
