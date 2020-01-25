BeginPackage["OneLoopDecaysUtils`"];

ExtractFormFactors::usage="given a FormCalc Amp object, returns the expression for the
amplitude as a list of the form
    {{matrix element, coefficient}, {matrix element, coefficient}, ...}
where the first element of each pair is a non-trivial Lorentz structure
and the second is its corresponding coefficient in the amplitude";
ToFSConventions::usage="replaces FeynArts/FormCalc symbols in an expression by their
equivalents used by SARAH and FlexibleSUSY";

GetGraphCombinatorialFactor::usage="returns the combinatorial factor associated with a
FeynArts FeynmanGraph object";
GetGraphNumber::usage="returns the graph number of a FeynArts FeymanGraph object";
GetGraphInsertions::usage="returns the insertions associated with a FeynmanGraph object as a list";
GetAdjacencyMatrix::usage="determines the adjacency matrix from a FeynArts Topology object";

CollectLorentzStructures::missterms="Missing terms from form factor decomposition: ``";
ExtractFormFactors::missterms="Could not associate the following terms with a matrix element: ``";
CheckFormCalcVersion::usage="";

Begin["`Private`"];

CheckFormCalcVersion[fcver_String] :=
    With[{version = Read[StringToStream[StringSplit[fcver][[2]]], Number]},
      If[!NumericQ[version],
        Print["Error: Could not identify the version of FormCalc"];
        Quit[1]
      ];
      If[!MemberQ[{9.5, 9.6}, version],
        Print["Error: This script only works with FormCalc version 9.5 or 9.6"];
        Quit[1]
      ];
    ];
CheckFormCalcVersion[___] :=
    (Print["String expected"]; Quit[1]);

GetGenericFieldSymbol[{field_[indices__], properties__}] := field;
GetGenericFieldSymbol[{field_, properties__}] := field;

GetGenericProcess[Rule[inFields_List, outFields_List]] :=
    Rule[GetGenericFieldSymbol /@ inFields, GetGenericFieldSymbol /@ outFields];

GetGraphCombinatorialFactor[FeynArts`FeynmanGraph[s_, level_][insertions__]] := s;

GetGraphNumber[FeynArts`FeynmanGraph[s_, level_ == n_][insertions__]] := n;

GetGraphInsertions[FeynArts`FeynmanGraph[s_, level_][insertions__]] := List[insertions];

IsLoopFunction[LoopTools`A0[__]] := True;
IsLoopFunction[LoopTools`A00[__]] := True;
IsLoopFunction[LoopTools`A0i[__]] := True;
IsLoopFunction[LoopTools`B0[__]] := True;
IsLoopFunction[LoopTools`DB0[__]] := True;
IsLoopFunction[LoopTools`B00[__]] := True;
IsLoopFunction[LoopTools`DB00[__]] := True;
IsLoopFunction[LoopTools`B001[__]] := True;
IsLoopFunction[LoopTools`B0i[__]] := True;
IsLoopFunction[LoopTools`B1[__]] := True;
IsLoopFunction[LoopTools`DB1[__]] := True;
IsLoopFunction[LoopTools`B11[__]] := True;
IsLoopFunction[LoopTools`DB11[__]] := True;
IsLoopFunction[LoopTools`B111[__]] := True;
IsLoopFunction[LoopTools`C0[__]] := True;
IsLoopFunction[LoopTools`C0i[__]] := True;
IsLoopFunction[LoopTools`D0[__]] := True;
IsLoopFunction[LoopTools`D0i[__]] := True;
IsLoopFunction[LoopTools`E0[__]] := True;
IsLoopFunction[LoopTools`E0i[__]] := True;
IsLoopFunction[LoopTools`F0[__]] := True;
IsLoopFunction[LoopTools`F0i[__]] := True;
IsLoopFunction[f_] := False;

IsNonTrivialLorentzStructure[FormCalc`Pair[FormCalc`k[i_Integer], FormCalc`k[j_Integer]]] := False;
IsNonTrivialLorentzStructure[FormCalc`Pair[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`Eps[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`DiracChain[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`WeylChain[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`Mat[__]] := True;
IsNonTrivialLorentzStructure[_] := False;

JoinStructureLists[args__StructureList] := Join[args];

CollectCoefficients[StructureList[pairs__]] :=
    Module[{likeTerms, sumTerms},
           likeTerms = Gather[List[pairs], (First[#1] === First[#2])&];
           sumTerms[terms_List] :=
               Module[{lorentzStructure},
                      lorentzStructure = First[First[terms]];
                      {lorentzStructure, Simplify[Plus @@ (Last /@ terms)]}
                     ];
           StructureList @@ (sumTerms /@ likeTerms)
          ];

MultiplyStructureLists[factors__StructureList] := Distribute[Times[factors], StructureList];

FactorOutLorentzStructure[expr_Times] :=
    Module[{factors, lorentz, scalar},
           factors = List @@ expr;
           factored = FactorOutLorentzStructure /@ factors;
           MultiplyStructureLists @@ factored
          ]

IsPositiveInteger[n_] := IntegerQ[n] && n > 0;

FactorOutLorentzStructure[Power[expr_, exponent_?IsPositiveInteger]] :=
    Module[{i},
           MultiplyStructureLists @@ Table[FactorOutLorentzStructure[expr], {i, 1, exponent}]
          ];

FactorOutLorentzStructure[expr_Plus] :=
   CollectCoefficients[JoinStructureLists @@ (FactorOutLorentzStructure /@ (List @@ expr))];

FactorOutLorentzStructure[expr_?IsNonTrivialLorentzStructure] := StructureList[{expr, 1}];

FactorOutLorentzStructure[expr_] := StructureList[{1, expr}];

ToTermList[expr_Plus] := List @@ expr;
ToTermList[expr_Times] := List[expr];

CollectLorentzStructures[FormCalc`Amp[process_][expr_]] :=
    Module[{couplingRules, loopFuncRules, compactExpr, terms, result, difference},
           couplingRules = Rule[#, Unique["coup"]]& /@ DeleteDuplicates[Cases[expr, G[_][_][__][__]]];
           loopFuncRules = Rule[#, Unique["lpFn"]]& /@ DeleteDuplicates[Cases[expr, f_[args__] /; IsLoopFunction[f[args]]]];
           compactExpr = expr /. couplingRules /. loopFuncRules;
           terms = ToTermList[Expand[compactExpr]];
           result = CollectCoefficients[JoinStructureLists @@ (FactorOutLorentzStructure /@ terms)];
           result = List @@ (result /. (Reverse /@ Join[couplingRules, loopFuncRules]));
           (* Check result is consistent with original expression *)
           difference = Simplify[expr - Plus @@ ((Times @@ #)& /@ result)];
           If[difference =!= 0,
              Message[CollectLorentzStructures::missterms, difference];
             ];
           result
          ];

ExtractFormFactors[FormCalc`Amp[process_][expr_]] :=
    CollectLorentzStructures[FormCalc`Amp[process][expr]];

SARAHLorentzIndex[i_Integer] := Symbol["lt" <> ToString[i]];

CouplingToSARAHCpRules[] :=
    {
     RuleDelayed[FeynArts`G[_][0][fields__][1], SARAH`Cp[fields][1]],
     RuleDelayed[FeynArts`G[_][0][fields__][FeynArts`NonCommutative[Global`ChiralityProjector[1]]], SARAH`Cp[fields][SARAH`PR]],
     RuleDelayed[FeynArts`G[_][0][fields__][FeynArts`NonCommutative[Global`ChiralityProjector[-1]]], SARAH`Cp[fields][SARAH`PL]],
     RuleDelayed[FeynArts`G[_][0][fields__][FormCalc`Private`ga[6]], SARAH`Cp[fields][SARAH`PR]],
     RuleDelayed[FeynArts`G[_][0][fields__][FormCalc`Private`ga[7]], SARAH`Cp[fields][SARAH`PL]],
     RuleDelayed[FeynArts`G[_][0][fields__][Global`MetricTensor[FeynArts`KI1[i1_], FeynArts`KI1[i2_]]], SARAH`Cp[fields][SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i2]]]],
       (*
     RuleDelayed[FeynArts`G[_][0][fields__][
       Global`MetricTensor[FeynArts`KI1[i1_], FeynArts`KI1[i2_]] (-FeynArts`Mom[i1_] + FeynArts`Mom[i2_]) +
           Global`MetricTensor[FeynArts`KI1[i1_], FeynArts`KI1[i3_]] (FeynArts`Mom[i1_] - FeynArts`Mom[i3_]) +
           Global`MetricTensor[FeynArts`KI1[i2_], FeynArts`KI1[i3_]] (-FeynArts`Mom[i2_] + FeynArts`Mom[i3_])
     ],
       SARAH`Cp[fields][
         SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i2]] (-SARAH`Mom[{fields}[[i1]]] + SARAH`Mom[{fields}[[i2]]]) +
         SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i3]] (SARAH`Mom[{fields}[[i1]]] - SARAH`Mom[{fields}[[i3]]]) +
         SARAH`g[SARAHLorentzIndex[i2], SARAHLorentzIndex[i3]] (-SARAH`Mom[{fields}[[i2]]] + SARAH`Mom[{fields}[[i3]]])
      ]
      ],
      *)
     RuleDelayed[
       FeynArts`G[_][0][fields__][
         Global`MetricTensor[FeynArts`KI1[i1_], FeynArts`KI1[i2_]]
         Global`MetricTensor[FeynArts`KI1[i3_], FeynArts`KI1[i4_]]
       ],
       SARAH`Cp[fields][
         SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i2]]
         SARAH`g[SARAHLorentzIndex[i3], SARAHLorentzIndex[i4]]
       ]
     ],
     RuleDelayed[
        FeynArts`G[_][0][fields__][FeynArts`NonCommutative[Global`DiracMatrix[FeynArts`KI1[i1_]], Global`ChiralityProjector[1]]],
        SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[SARAHLorentzIndex[i1]], SARAH`PR]]
     ],
     RuleDelayed[
        FeynArts`G[_][0][fields__][FeynArts`NonCommutative[Global`DiracMatrix[FeynArts`KI1[i1_]], Global`ChiralityProjector[-1]]],
        SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[SARAHLorentzIndex[i1]], SARAH`PL]]
     ],
     RuleDelayed[
        FeynArts`G[_][0][fields__][FeynArts`Mom[i1_] - FeynArts`Mom[i2_]],
        SARAH`Cp[fields][SARAH`Mom[{fields}[[i1]]] - SARAH`Mom[{fields}[[i2]]]]
     ],
     RuleDelayed[
        FeynArts`G[_][0][fields__][FeynArts`Mom[i1_]],
        SARAH`Cp[fields][SARAH`Mom[{fields}[[i1]]]]
     ],
       (*
     RuleDelayed[FeynArts`G[_][0][fields__]["d_"[FeynArts`KI1[i1_], FeynArts`KI1[i2_]]], SARAH`Cp[fields][SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i2]]]],
      RuleDelayed[FeynArts`G[_][0][fields__][NonCommutativeMultiply[FeynArts`KI1[i1_], FormCalc`Private`ga[6]]], SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[SARAHLorentzIndex[i1]], SARAH`PR]]],
     RuleDelayed[FeynArts`G[_][0][fields__][NonCommutativeMultiply[FeynArts`KI1[i1_], FormCalc`Private`ga[7]]], SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[SARAHLorentzIndex[i1]], SARAH`PL]]],
     RuleDelayed[FeynArts`G[_][0][fields__]["d_"[FeynArts`KI1[i1_], FeynArts`KI1[i2_]] (FeynArts`Mom[i2_] - FeynArts`Mom[i1_])
                                            + "d_"[FeynArts`KI1[i1_], FeynArts`KI1[i3_]] (FeynArts`Mom[i1_] - FeynArts`Mom[i3_])
                                            + "d_"[FeynArts`KI1[i2_], FeynArts`KI1[i3_]] (FeynArts`Mom[i3_] - FeynArts`Mom[i2_])],
                 Print[fields];SARAH`Cp[fields][SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i2]] (SARAH`Mom[{fields}[[i2]]] - SARAH`Mom[{fields}[[i1]]])
                                  + SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i3]] (SARAH`Mom[{fields}[[i1]]] - SARAH`Mom[{fields}[[i3]]])
                                  + SARAH`g[SARAHLorentzIndex[i2], SARAHLorentzIndex[i3]] (SARAH`Mom[{fields}[[i3]]] - SARAH`Mom[{fields}[[i2]]])]],
     RuleDelayed[FeynArts`G[_][0][fields__]["d_"[FeynArts`KI1[i1_], FeynArts`KI1[i2_]] "d_"[FeynArts`KI1[i3_], FeynArts`KI1[i4_]]],
                 SARAH`Cp[fields][SARAH`g[SARAHLorentzIndex[i1], SARAHLorentzIndex[i2]] SARAH`g[SARAHLorentzIndex[i3], SARAHLorentzIndex[i4]]]],
*)
     RuleDelayed[FeynArts`G[_][0][fields__][lor__], (Print["Unknown coupling ", fields, FullForm@lor];Quit[1])]
    };

IsFeynArtsField[FeynArts`S] := True;
IsFeynArtsField[FeynArts`F] := True;
IsFeynArtsField[FeynArts`V] := True;
IsFeynArtsField[FeynArts`U] := True;
IsFeynArtsField[sym_] := False;

RemoveFieldLineSpecifiers[] :=
    {
     RuleDelayed[field_[indices__, FeynArts`Loop] /; IsFeynArtsField[field], field[indices]],
     RuleDelayed[field_[indices__, FeynArts`External] /; IsFeynArtsField[field], field[indices]],
     RuleDelayed[field_[indices__, FeynArts`Internal] /; IsFeynArtsField[field], field[indices]]
    };

ToSARAHCouplings[expr_] := expr /. CouplingToSARAHCpRules[] /. RemoveFieldLineSpecifiers[];

LoopFunctionToFSNotationRules[] :=
    {
     RuleDelayed[LoopTools`B0i[LoopTools`bb0, args__], SARAH`B0[args]],
     RuleDelayed[LoopTools`B0i[LoopTools`bb00, args__], SARAH`B00[args]],
     RuleDelayed[LoopTools`B0i[LoopTools`bb1, args__], SARAH`B1[args]],
     RuleDelayed[LoopTools`B0i[LoopTools`bb11, args__], SARAH`B11[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc0, args__], SARAH`C0[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc1, args__], SARAH`C1[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc2, args__], SARAH`C2[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc00, args__], SARAH`C00[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc11, args__], SARAH`C11[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc12, args__], SARAH`C12[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc22, args__], SARAH`C22[args]]
    };

ToFSLoopFunctions[expr_] := expr /. LoopFunctionToFSNotationRules[];

ToFSConventions[expr_] := ToFSLoopFunctions[ToSARAHCouplings[expr]];

(* Converts a list of integer pairs {{from, to}, {from, to}, ...}
   representing a graph to the corresponding adjacency matrix *)
GetAdjacencyMatrix[edgeList_List, undirected_:True] :=
    Module[{edgeTest, edgeCount},
           edgeTest = If[undirected,
                         Function[{e1, e2}, Or[e1[[1]] == e2[[1]] && e1[[2]] == e2[[2]],
                                               e1[[2]] == e2[[1]] && e1[[1]] == e2[[2]]]],
                         Function[{e1, e2}, e1[[1]] == e2[[1]] && e1[[2]] == e2[[2]]]
                        ];
           (* Use convention in which self-edges contribute 2 to the diagonal entries
              for undirected graphs *)
           edgeCount = Tally[edgeList, edgeTest] /. If[undirected, {{{i_, i_}, c_} :> {{i, i}, 2 c}}, {}];
           If[undirected,
              edgeCount = DeleteDuplicates[Join[edgeCount, {Reverse[#[[1]]], #[[2]]}& /@ edgeCount]];
             ];
           SparseArray[Rule[#[[1]], #[[2]]]& /@ edgeCount]
          ];

GetAdjacencyMatrix[FeynArts`Topology[s_][propagators__], undirected_:True] :=
    Module[{propList, edgeList},
           propList = List[propagators];
           edgeList = {#[[1]], #[[2]]}& /@ propList;
           Normal[GetAdjacencyMatrix[edgeList, undirected]]
          ];

End[];

EndPackage[];
