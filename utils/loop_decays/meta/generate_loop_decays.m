status = 0;

baseDir = DirectoryName[$InputFileName];
fsDir = Nest[ParentDirectory, baseDir, 3];

AppendTo[$Path, FileNameJoin[{fsDir, "meta"}]];
AppendTo[$Path, baseDir];

loadFSStatus = Needs["WriteOut`"];
If[loadFSStatus === $Failed,
   Quit[1];
  ];

Clear[SARAH`AntiField];
Clear[Vertex];

loadUtilsStatus = Needs["OneLoopDecaysUtils`"];
If[loadUtilsStatus === $Failed,
   Quit[1];
  ];

ToFieldString[-field_[indices__Integer]] := "c" <> ToFieldString[field[indices]];
ToFieldString[-field_] := "c" <> ToFieldString[field];
ToFieldString[field_[indices__Integer]] := ToString[field] <> StringJoin[ToString /@ List[indices]];
ToFieldString[field_] := ToString[field];

ToGenericFieldString[field_[indices__Integer]] := ToFieldString[field];
ToGenericFieldString[field_] := ToFieldString[field];

IsNonZeroDiagram[diagram_] :=
    Module[{formFactors},
           formFactors = Last[diagram];
           Select[formFactors, (Last[#] =!= 0)&] =!= {}
          ];

CreateProcessName[Rule[{initial_}, finalState_List]] :=
    ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]];

CreateFormFactorsOutputFileName[process_] :=
    "form_factors_" <> CreateProcessName[process] <> ".m";

ReadFormFactorsExprs[process_, dataDir_] :=
    Module[{formFactorsFile},
           formFactorsFile = FileNameJoin[{dataDir, CreateFormFactorsOutputFileName[process]}];
           If[FileExistsQ[formFactorsFile],
              result = Get[formFactorsFile],
              Print["Error: input file ", formFactorsFile, " does not exist."];
              Quit[2];
             ]
          ];

GetInsertionsAsList[FeynmanGraph[props__][insertions__]] := List[insertions];
GetInsertionsAsList[insertions_List] := insertions;

GetAmplitudeCType[Rule[{S}, {S, S}]] := "Decay_amplitude_SSS";
GetAmplitudeCType[Rule[{S}, {F, F}]] := "Decay_amplitude_SFF";
GetAmplitudeCType[Rule[{S}, {S, V}]] := "Decay_amplitude_SSV";
GetAmplitudeCType[Rule[{S}, {V, S}]] := "Decay_amplitude_SSV";
GetAmplitudeCType[Rule[{S}, {V, V}]] := "Decay_amplitude_SVV";
GetAmplitudeCType[Rule[{F}, {F, S}]] := "Decay_amplitude_FFS";
GetAmplitudeCType[Rule[{F}, {S, F}]] := "Decay_amplitude_FFS";
GetAmplitudeCType[Rule[{F}, {F, V}]] := "Decay_amplitude_FFV";
GetAmplitudeCType[Rule[{F}, {V, F}]] := "Decay_amplitude_FFV";

GetFormFactorName[Rule[{S}, {S, S}], 1] := "form_factor";

GetFormFactorName[Rule[{S}, {S, V}], Pair[ec[3], k[1]]] := "form_factor";

GetFormFactorName[Rule[{S}, {V, V}], Pair[ec[2], ec[3]]] := "form_factor_g";
GetFormFactorName[Rule[{S}, {V, V}], Pair[ec[2], k[2]] Pair[ec[3], k[2]]] := "form_factor_11";
GetFormFactorName[Rule[{S}, {V, V}], Pair[ec[2], k[2]] Pair[ec[3], k[3]]] := "form_factor_12";
GetFormFactorName[Rule[{S}, {V, V}], Pair[ec[2], k[3]] Pair[ec[3], k[2]]] := "form_factor_21";
GetFormFactorName[Rule[{S}, {V, V}], Pair[ec[2], k[3]] Pair[ec[3], k[3]]] := "form_factor_22";
GetFormFactorName[Rule[{S}, {V, V}], Eps[ec[2], ec[3], k[2], k[3]]] := "form_factor_eps";

GetFormFactorName[Rule[{S}, {F, F}], Mat[DiracChain[Spinor[k[2], Mass[F[Index[Generic, 2]]], 1], 6,
                                                    Spinor[k[3], Mass[F[Index[Generic, 3]]], -1]]]] := "form_factor_right";
GetFormFactorName[Rule[{S}, {F, F}], Mat[DiracChain[Spinor[k[2], Mass[F[Index[Generic, 2]]], 1], 7,
                                                    Spinor[k[3], Mass[F[Index[Generic, 3]]], -1]]]] := "form_factor_left";

ExternalMomentumSymbol[idx_] := k[idx];
CreateExternalMomentumCString[ExternalMomentumSymbol[idx_]] := "mext" <> ToString[idx];
CreateExternalMomentumCString[Pair[k[i_], k[i_]]] := "k" <> ToString[i] <> "sq";

IsExternalMomentumSquared[Pair[ExternalMomentumSymbol[i_], ExternalMomentumSymbol[i_]]] := True;
IsExternalMomentumSquared[_] := False;

GetExternalMomentumCType[psq_] := CConversion`CreateCType[CConversion`ScalarType[realScalarCType]];

IsLoopMass[Mass[fields_[indices__], Loop]] := True;
IsLoopMass[_] := False;

GetLoopMassCType[msq_] := CConversion`CreateCType[CConversion`ScalarType[realScalarCType]];

CreateLoopMassCString[Mass[field_[Index[Generic, idx_], indices___], Loop]^2] :=
    "m" <> ToGenericFieldString[field] <> ToString[idx] <> "sq";

CreateLoopMassCString[Mass[field_[Index[Generic, idx_], indices___], Loop]] :=
    "m" <> ToGenericFieldString[field] <> ToString[idx];

GetCouplingCType[cp_] := "const " <> CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <> "&";

GetCouplingCType[couplings_List] := GetCouplingCType /@ couplings;

CreateCouplingCString[SARAH`Cp[fields__]] :=
    Module[{edgeLabels},
           edgeLabels = ToFieldString /@ (List[fields] /. (field_[___, Index[Generic, i_], ___] :> field[i]));
           "Cp" <> StringJoin[edgeLabels]
          ];

CreateCouplingCString[SARAH`Cp[fields__][SARAH`PL]] :=
    CreateCouplingCString[SARAH`Cp[fields]] <> "PL";

CreateCouplingCString[SARAH`Cp[fields__][SARAH`PR]] :=
    CreateCouplingCString[SARAH`Cp[fields]] <> "PR";

CreateCouplingCString[SARAH`Cp[fields__][SARAH`LorentzProduct[SARAH`gamma[_], SARAH`PL]]] :=
    CreateCouplingCString[SARAH`Cp[fields]] <> "PL";

CreateCouplingCString[SARAH`Cp[fields__][SARAH`LorentzProduct[SARAH`gamma[_], SARAH`PR]]] :=
    CreateCouplingCString[SARAH`Cp[fields]] <> "PR";

CreateCouplingCString[SARAH`Cp[fields__][lor_]] :=
    CreateCouplingCString[SARAH`Cp[fields]];

CreateCouplingCString[couplings_List] :=
    Module[{i, j, couplingStrs, repeatedCounts, count},
           couplingStrs = CreateCouplingCString /@ couplings;
           repeatedCounts = Select[Tally[couplingStrs], (#[[2]] > 1)&];
           For[i = 1, i <= Length[repeatedCounts], i++,
               count = 1;
               For[j = 1, j <= Length[couplingStrs], j++,
                   If[couplingStrs[[j]] == repeatedCounts[[i,1]],
                      couplingStrs = ReplacePart[couplingStrs, j -> couplingStrs[[j]] <> ToString[count]];
                      count++;
                     ];
                  ];
              ];
           couplingStrs
          ];

GetExternalMomentumIndex[ExternalMomentumSymbol[i_]] := i;
GetExternalMomentumIndex[Pair[ExternalMomentumSymbol[i_], ExternalMomentumSymbol[i_]]] := i;

GetExternalMomenta[process_, diagram_] :=
    Module[{nExt},
           nExt = Length[Join[process[[1]], process[[2]]]];
           Table[ExternalMomentumSymbol[i], {i, 1, nExt}]
          ];

SortMassesByGenericIndex[masses_List] :=
    Module[{orderingFn},
           orderingFn[Mass[field1_[___, Index[Generic, i1_], ___], Loop],
                      Mass[field2_[___, Index[Generic, i2_], ___], Loop]] := i1 <= i2;
           Sort[masses, orderingFn]
          ];

GetLoopMassField[Mass[field_[___, Index[Generic, i_], ___], Loop]] := field[i];

GetLoopMasses[process_, diagram_] :=
    Module[{formFactors, masses},
           formFactors = Last[diagram];
           masses = Cases[formFactors, Mass[field_, Loop], {0, Infinity}];
           masses = DeleteDuplicates[masses];
           SortMassesByGenericIndex[masses]
          ];

GroupByFieldContent[couplings_List] :=
    Module[{fieldsGather, result},
           fieldsGather[SARAH`Cp[fields__]] := List[fields];
           fieldsGather[SARAH`Cp[fields__][lor_]] := List[fields];
           GatherBy[couplings, fieldsGather]
          ];

SortChiralCouplings[couplings_List] :=
    Module[{chiralCouplings, chiralSort},
           chiralCouplings = Cases[couplings, SARAH`Cp[fields__][SARAH`PL|SARAH`PR], {0, Infinity}];
           If[chiralCouplings === {},
              Return[couplings];
             ];
           If[chiralCouplings =!= couplings,
              Print["Error: cannot sort mixture of chiral and non-chiral couplings"];
              Quit[1];
             ];
           chiralSort[SARAH`Cp[fields__][SARAH`PL], SARAH`Cp[fields__][SARAH`PR]] := True;
           chiralSort[SARAH`Cp[fields__][SARAH`PR], SARAH`Cp[fields__][SARAH`PL]] := False;
           Sort[couplings, chiralSort]
          ];

GetCouplings[process_, diagram_] :=
    Module[{formFactors},
           formFactors = Last[diagram];
           couplings = Join[Cases[formFactors, SARAH`Cp[fields__], {0, Infinity}],
                            Cases[formFactors, SARAH`Cp[fields__][lor__], {0, Infinity}]];
           couplings = DeleteDuplicates[couplings];
           SortChiralCouplings /@ GroupByFieldContent[couplings]
          ];

CreateProcessString[Rule[{initial_}, finalState_List]] :=
    Module[{initialField, finalFields},
           initialField = ToFieldString[initial];
           finalFields = ToFieldString /@ finalState;
           initialField <> " -> " <> StringJoin[finalFields]
          ];

CreateGraphIDString[GraphID[Topology == t_, Generic == g_, Number == n_]] :=
    "t" <> ToString[t] <> "g" <> ToString[g] <> "n" <> ToString[n];

(* @todo fix an appropriate naming convention *)
CreateOneLoopDiagramName[process_, graphID_, insertions_] :=
    Module[{processName, idString, nExt, internalFields, internalFieldsOrder},
           processName = CreateProcessName[process];
           idString = CreateGraphIDString[graphID];
           nExt = Length[Join[process[[1]], process[[2]]]];
           internalFields = DeleteCases[GetInsertionsAsList[insertions], (Field[i_] -> field_) /; i <= nExt];
           internalFieldsOrder = Ordering[internalFields /. (Field[i_] -> f_) :> i];
           internalFieldsLabel = StringJoin[ToGenericFieldString /@ ((#[[2]]& /@ internalFields)[[internalFieldsOrder]])];
           "diagram_" <> processName <> "_" <> idString <> "_" <> internalFieldsLabel
          ];

CreateDiagramEvaluatorName[process_, diagram_] :=
    Module[{diagramName},
           diagramName = CreateOneLoopDiagramName[process, diagram[[1]], diagram[[3]]];
           "calculate_" <> diagramName
          ];

CreateOneLoopDiagramDeclaration[process_, diagram_] :=
    Module[{returnType, externalMomenta, loopMasses, couplings,
            formFactors, args = ""},
           returnType = GetAmplitudeCType[process];
           externalMomenta = GetExternalMomenta[process, diagram];
           loopMasses = GetLoopMasses[process, diagram];
           couplings = Flatten[GetCouplings[process, diagram]];
           formFactors = Last[diagram];
           args = StringJoin[Riffle[Join[GetExternalMomentumCType /@ externalMomenta,
                                         GetLoopMassCType /@ loopMasses,
                                         GetCouplingCType /@ couplings], ", "]] <>
                  ", double" <> If[!FreeQ[formFactors, Finite], ", double", ""];
           returnType <> " " <> CreateDiagramEvaluatorName[process, diagram] <> "(" <> args <> ");"
          ];

CreateOneLoopDiagramDeclarations[process_, diagrams_] :=
    StringJoin[Riffle[CreateOneLoopDiagramDeclaration[process, #]& /@ diagrams, "\n\n"]];

GetExternalMomentaVars[process_, diagram_] :=
    Module[{i, externalMomenta, vars, types, fields},
           externalMomenta = GetExternalMomenta[process, diagram];
           vars = CreateExternalMomentumCString /@ externalMomenta;
           types = GetExternalMomentumCType /@ externalMomenta;
           fields = Table[Field[i][Index[Generic, i]], {i, 1, Length[externalMomenta]}] /. (List @@ diagram[[3]]);
           MapThread[{#1, #2, #3, #4}&, {externalMomenta, vars, types, fields}]
          ];

GetLoopMassesVars[process_, diagram_] :=
    Module[{loopMasses, vars, types},
           loopMasses = GetLoopMasses[process, diagram];
           vars = CreateLoopMassCString /@ loopMasses;
           types = GetLoopMassCType /@ loopMasses;
           MapThread[{#1, #2, #3}&, {loopMasses, vars, types}]
          ];

GetCouplingsVars[process_, diagram_] :=
    Module[{couplings, vars, types},
           couplings = GetCouplings[process, diagram];
           vars = CreateCouplingCString /@ couplings;
           types = GetCouplingCType /@ couplings;
           MapThread[{#1, #2, #3}&, Flatten /@ {couplings, vars, types}]
          ];

CreateCouplingDescription[SARAH`Cp[fields__]] :=
    "coupling " <> ToString[SARAH`Cp[fields] /. f_[___, Index[Generic, i_], ___] :> f[i]]

CreateCouplingDescription[SARAH`Cp[fields__][1]] :=
    CreateCouplingDescription[SARAH`Cp[fields]];

CreateCouplingDescription[SARAH`Cp[fields__][lor__]] :=
    "coupling " <> ToString[SARAH`Cp[fields][lor] /. f_[___, Index[Generic, i_], ___] :> f[i]]

FillAmplitudeMasses[Rule[{S}, {S, S}], diagram_, struct_:"result"] :=
    Module[{topology, incomingIndex, outgoingIndices, incomingMass, outgoingOne, outgoingTwo},
           topology = diagram[[2]];
           incomingIndex = Cases[topology, Propagator[Incoming][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[incomingIndex] != 1,
              Print["Error: number of incoming fields is not one."];
              Quit[1];
             ];
           outgoingIndices = Cases[topology, Propagator[Outgoing][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[outgoingIndices] != 2,
              Print["Error: number of outgoing fields is not two."];
              Quit[1];
             ];
           incomingMass = CreateExternalMomentumCString[ExternalMomentumSymbol[First[incomingIndex]]];
           outgoingOne = CreateExternalMomentumCString[ExternalMomentumSymbol[First[outgoingIndices]]];
           outgoingTwo = CreateExternalMomentumCString[ExternalMomentumSymbol[Last[outgoingIndices]]];
           struct <> ".m_decay = " <> incomingMass <> ";\n" <>
           struct <> ".m_scalar_1 = " <> outgoingOne  <> ";\n"<>
           struct <> ".m_scalar_2 = " <> outgoingTwo  <> ";\n"
          ];

FillAmplitudeMasses[Rule[{S}, {V, V}], diagram_, struct_:"result"] :=
    Module[{topology, incomingIndex, outgoingIndices, incomingMass, outgoingOne, outgoingTwo},
           topology = diagram[[2]];
           incomingIndex = Cases[topology, Propagator[Incoming][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[incomingIndex] != 1,
              Print["Error: number of incoming fields is not one."];
              Quit[1];
             ];
           outgoingIndices = Cases[topology, Propagator[Outgoing][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[outgoingIndices] != 2,
              Print["Error: number of outgoing fields is not two."];
              Quit[1];
             ];
           incomingMass = CreateExternalMomentumCString[ExternalMomentumSymbol[First[incomingIndex]]];
           outgoingOne = CreateExternalMomentumCString[ExternalMomentumSymbol[First[outgoingIndices]]];
           outgoingTwo = CreateExternalMomentumCString[ExternalMomentumSymbol[Last[outgoingIndices]]];
           struct <> ".m_decay = " <> incomingMass <> ";\n" <>
           struct <> ".m_vector_1 = " <> outgoingOne  <> ";\n"<>
           struct <> ".m_vector_2 = " <> outgoingTwo  <> ";\n"
          ];

FillAmplitudeMasses[Rule[{S}, {V, S}], diagram_, struct_:"result"] :=
    FillAmplitudeMasses[Rule[{S}, {S, V}] , diagram, struct];

FillAmplitudeMasses[Rule[{S}, {S, V}], diagram_, struct_:"result"] :=
    Module[{topology, insertions, incomingIndex, outgoingIndices,
            scalarIndex, vectorIndex, incomingMass, outgoingScalar, outgoingVector},
           topology = diagram[[2]];
           insertions = List @@ diagram[[3]] /. (Field[i_] -> f_) :> {i, f};
           incomingIndex = Cases[topology, Propagator[Incoming][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[incomingIndex] != 1,
              Print["Error: number of incoming fields is not one."];
              Quit[1];
             ];
           outgoingIndices = Cases[topology, Propagator[Outgoing][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[outgoingIndices] != 2,
              Print["Error: number of outgoing fields is not two."];
              Quit[1];
             ];
           scalarIndex = Cases[insertions, ({i_, S} /; MemberQ[outgoingIndices, i]) :> i];
           vectorIndex = Complement[outgoingIndices, scalarIndex];
           incomingMass = CreateExternalMomentumCString[ExternalMomentumSymbol[First[incomingIndex]]];
           outgoingScalar = CreateExternalMomentumCString[ExternalMomentumSymbol[First[scalarIndex]]];
           outgoingVector = CreateExternalMomentumCString[ExternalMomentumSymbol[First[vectorIndex]]];
           struct <> ".m_decay = " <> incomingMass <> ";\n" <>
           struct <> ".m_scalar = " <> outgoingScalar  <> ";\n"<>
           struct <> ".m_vector = " <> outgoingVector  <> ";\n"
          ];

FillAmplitudeMasses[Rule[{S}, {F, F}], diagram_, struct_:"result"] :=
    Module[{topology, incomingIndex, outgoingIndices, incomingMass, outgoingOne, outgoingTwo},
           topology = diagram[[2]];
           incomingIndex = Cases[topology, Propagator[Incoming][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[incomingIndex] != 1,
              Print["Error: number of incoming fields is not one."];
              Quit[1];
             ];
           outgoingIndices = Cases[topology, Propagator[Outgoing][vertices___, Field[i_]] :> i, {0, Infinity}];
           If[Length[outgoingIndices] != 2,
              Print["Error: number of outgoing fields is not two."];
              Quit[1];
             ];
           incomingMass = CreateExternalMomentumCString[ExternalMomentumSymbol[First[incomingIndex]]];
           outgoingOne = CreateExternalMomentumCString[ExternalMomentumSymbol[First[outgoingIndices]]];
           outgoingTwo = CreateExternalMomentumCString[ExternalMomentumSymbol[Last[outgoingIndices]]];
           struct <> ".m_decay = " <> incomingMass <> ";\n" <>
           struct <> ".m_fermion_1 = " <> outgoingOne  <> ";\n"<>
           struct <> ".m_fermion_2 = " <> outgoingTwo  <> ";\n"
          ];

CreateOneLoopDiagramDocString[process_, diagram_] :=
    Module[{idString, brief, externalMomenta, momentaInfo,
            loopMasses, massesInfo, couplings, couplingsInfo, docString = ""},
           idString = CreateGraphIDString[diagram[[1]]];
           brief = " * @brief Evaluates " <> ToUpperCase[idString] <>
                   " diagram for process " <> CreateProcessString[process] <> "\n";
           docString = docString <> brief <> " *\n";
           externalMomenta = GetExternalMomentaVars[process, diagram];
           momentaInfo = (" * @param[in] " <> #[[2]] <> " mass of external field " <>
                          ToString[#[[4]] /. field_[___, Index[Generic, i_], ___] :> field[i]])& /@ externalMomenta;
           loopMasses = GetLoopMassesVars[process, diagram];
           massesInfo = (" * @param[in] " <> #[[2]] <> " mass of internal field " <>
                          ToString[GetLoopMassField[#[[1]]]])& /@ loopMasses;
           couplings = GetCouplingsVars[process, diagram];
           couplingsInfo = (" * @param[in] " <> #[[2]] <> " " <>
                          CreateCouplingDescription[#[[1]]])& /@ couplings;
           docString = docString <> StringJoin[Riffle[momentaInfo, "\n"]] <> "\n" <>
                       StringJoin[Riffle[massesInfo, "\n"]] <> "\n" <>
                       StringJoin[Riffle[couplingsInfo, "\n"]] <> "\n";
           docString = docString <> " *\n * @return value of the one-loop diagram\n";
           "/**\n" <> docString <> " */\n"
          ];

IsSARAHLoopFunction[SARAH`A0] := True;
IsSARAHLoopFunction[SARAH`B0] := True;
IsSARAHLoopFunction[SARAH`B1] := True;
IsSARAHLoopFunction[SARAH`C0] := True;
IsSARAHLoopFunction[SARAH`C1] := True;
IsSARAHLoopFunction[SARAH`C2] := True;
IsSARAHLoopFunction[SARAH`C00] := True;
IsSARAHLoopFunction[SARAH`C11] := True;
IsSARAHLoopFunction[SARAH`C12] := True;
IsSARAHLoopFunction[SARAH`C22] := True;
IsSARAHLoopFunction[f_] := False;

CreateSavedLoopFunctionName[SARAH`A0[args__]] := "a0tmp";
CreateSavedLoopFunctionName[SARAH`B0[args__]] := "b0tmp";
CreateSavedLoopFunctionName[SARAH`B1[args__]] := "b1tmp";
CreateSavedLoopFunctionName[SARAH`C0[args__]] := "c0tmp";
CreateSavedLoopFunctionName[SARAH`C1[args__]] := "c1tmp";
CreateSavedLoopFunctionName[SARAH`C2[args__]] := "c2tmp";
CreateSavedLoopFunctionName[SARAH`C00[args__]] := "c00tmp";
CreateSavedLoopFunctionName[SARAH`C11[args__]] := "c11tmp";
CreateSavedLoopFunctionName[SARAH`C12[args__]] := "c12tmp";
CreateSavedLoopFunctionName[SARAH`C22[args__]] := "c22tmp";

CreateLoopFunctionArgumentName[Pair[ExternalMomentumSymbol[i_], ExternalMomentumSymbol[j_]]] :=
    CreateExternalMomentumCString[ExternalMomentumSymbol[i]] <> "*" <>
    CreateExternalMomentumCString[ExternalMomentumSymbol[j]];

CreateLoopFunctionArgumentName[arg_^2 /; IsLoopMass[arg]] :=
    CreateLoopMassCString[arg] <> "*" <> CreateLoopMassCString[arg];

CreateLoopFunctionArgumentName[arg_ /; IsLoopMass[arg]] :=
    CreateLoopMassCString[arg];

CallLoopFunction[fn_[args__], renScale_String, namespace_:"passarino_veltman"] :=
    Module[{argsStr},
           argsStr = StringJoin[Riffle[CreateLoopFunctionArgumentName /@ List[args], ", "]];
           namespace <> "::" <> ToString[fn] <> "(" <> argsStr <> ", " <> renScale <> "*" <> renScale <> ")"
          ];

SaveLoopIntegrals[diagram_, renScale_String] :=
    Module[{formFactors, loopFunctions, tmpVars, savedValues, subs},
           formFactors = Last[diagram];
           loopFunctions = Sort[DeleteDuplicates[Cases[formFactors, fn_[args__] /; IsSARAHLoopFunction[fn], {0, Infinity}]]];
           tmpVars = MapIndexed[(CreateSavedLoopFunctionName[#1] <> ToString[First[#2]])&, loopFunctions];
           savedValues = MapThread[("const auto " <> #1 <> " = " <> CallLoopFunction[#2, renScale] <> ";\n")&, {tmpVars, loopFunctions}];
           savedValues = StringJoin[savedValues];
           subs = MapThread[Rule[#1, Symbol[#2]]&, {loopFunctions, tmpVars}];
           {savedValues, subs}
          ];

CreateOneLoopDiagramDefinition[process_, diagram_] :=
    Module[{returnType, externalMomenta, externalMomentaArgs, externalMomentaSubs,
            externalMassSubs, loopMasses, loopMassesArgs, loopMassesSubs,
            couplings, couplingsArgs, couplingsSubs, argSubs,
            saveLoopIntegrals, loopIntegralSubs,
            formFactorExprs, fillExternalMasses, calculateFormFactors, renScale = "scale",
            args = "", body = "", docString = ""},
           returnType = GetAmplitudeCType[process];

           externalMomenta = GetExternalMomentaVars[process, diagram];
           externalMomentaArgs = (#[[3]] <> " " <> #[[2]])& /@ externalMomenta;
           externalMomentaSubs = Flatten[{Rule[Pair[#[[1]], #[[1]]], Symbol[#[[2]]]^2], Rule[#[[1]], Symbol[#[[2]]]]}& /@ externalMomenta];
           externalMassSubs = Table[Rule[Mass[Field[i][Index[Generic, i]]],
                                         Symbol[CreateExternalMomentumCString[ExternalMomentumSymbol[i]]]], {i, 1, 3}] /. (List @@ diagram[[3]]);
           loopMasses = GetLoopMassesVars[process, diagram];
           loopMassesArgs = (#[[3]] <> " " <> #[[2]])& /@ loopMasses;
           loopMassesSubs = Rule[#[[1]], Symbol[#[[2]]]]& /@ loopMasses;
           couplings = GetCouplingsVars[process, diagram];
           couplingsArgs = (#[[3]] <> " " <> #[[2]])& /@ couplings;
           couplingsSubs = Rule[#[[1]], Symbol[#[[2]]]]& /@ couplings;
           argSubs = Join[couplingsSubs, loopMassesSubs, externalMomentaSubs];

           {saveLoopIntegrals, loopIntegralSubs} = SaveLoopIntegrals[diagram, renScale];

           body = body <> saveLoopIntegrals <> "\n";

           formFactorExprs = Last[diagram];
           formFactorValues = {GetFormFactorName[process, #[[1]]], "oneOver16PiSqr*(" <>
                               CConversion`RValueToCFormString[Simplify[(16 Pi^2 #[[2]]) /. loopIntegralSubs /. argSubs /. externalMassSubs]] <> ")"}& /@ formFactorExprs;
           fillExternalMasses = FillAmplitudeMasses[process, diagram, "result"];
           calculateFormFactors = returnType <> " result;\n\n" <> fillExternalMasses <> "\n" <>
                                  StringJoin[("result." <> #[[1]] <> " = " <> #[[2]] <> ";\n")& /@ formFactorValues];
           body = body <> StringReplace[calculateFormFactors, "Finite" -> "finite"] <> "\nreturn result;\n";

           args = StringJoin[Riffle[externalMomentaArgs, ", "]] <> ",\n" <>
                  StringJoin[Riffle[loopMassesArgs, ", "]] <> ",\n" <>
                  StringJoin[Riffle[couplingsArgs, ", "]] <>
                  ",\ndouble " <> renScale <>
                  If[!FreeQ[formFactorExprs, Finite], ", double finite", ""];

           docString = CreateOneLoopDiagramDocString[process, diagram];

           docString <>
           returnType <> " " <> CreateDiagramEvaluatorName[process, diagram] <> "(\n" <>
           TextFormatting`IndentText[args] <> ")\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateOneLoopDiagramDefinitions[process_, diagrams_] :=
    StringJoin[Riffle[CreateOneLoopDiagramDefinition[process, #]& /@ diagrams, "\n\n"]];

CreateDiagramEvaluators[process_, diagrams_] :=
    Module[{decls, defs},
           decls = CreateOneLoopDiagramDeclarations[process, diagrams];
           defs = CreateOneLoopDiagramDefinitions[process, diagrams];
           {decls, defs}
          ];

GetNumberOfVertices[topology_] :=
    Length[DeleteDuplicates[Cases[List @@ topology, Vertex[_][_], {0, Infinity}]]];

GetDirectedAdjacencyMatrix[topology_] :=
    Module[{propagators, edges, nVertices},
           propagators = List @@ topology;
           edges = propagators /. Propagator[type_][Vertex[i_][from_], Vertex[j_][to_], Field[k_]] :> {from, to};
           nVertices = GetNumberOfVertices[topology];
           Normal[SparseArray[Rule @@@ Tally[edges], {nVertices, nVertices}]]
          ];

GetUndirectedAdjacencyMatrix[topology_] :=
    Module[{directedAdjacencyMatrix},
           directedAdjacencyMatrix = GetDirectedAdjacencyMatrix[topology];
           directedAdjacencyMatrix + Transpose[directedAdjacencyMatrix]
          ];

GetEdgeLabels[topology_] :=
    Module[{propagators},
          propagators = List @@ topology;
          propagators /. Propagator[type_][Vertex[i_][from_], Vertex[j_][to_], Field[k_]] :> Rule[{from, to}, Field[k]]
         ];

GetCouplings[expr_] :=
    Module[{allCouplings},
           allCouplings = Join[Cases[expr, SARAH`Cp[fields__], {0, Infinity}],
                               Cases[expr, SARAH`Cp[fields__][lor__], {0, Infinity}]];
           DeleteDuplicates[allCouplings]
          ];

GetEdgeLists[adjacencyMatrix_List, vertexLabels_List] :=
    Module[{i, k},
           edgeTriples = Flatten[(Flatten[Position[adjacencyMatrix[[#]], Except[0], Heads -> False]]
                                  /. {i_Integer :> Table[{#, i, k}, {k, 1, adjacencyMatrix[[#,i]]}]}), 1]&
                         /@ vertexLabels
          ];

GetPropagatorType[Propagator[Incoming][info__]] := Incoming;
GetPropagatorType[Propagator[Outgoing][info__]] := Outgoing;
GetPropagatorType[Propagator[Loop[i_]][info__]] := Loop;

GetInternalEdgeLists[topology_] :=
    Module[{propagators, internalPropagators, internalEdges},
           propagators = List @@ topology;
           internalPropagators = Select[propagators, (GetPropagatorType[#] === Loop)&];
           internalEdges = Cases[internalPropagators, Propagator[type_][Vertex[i_][from_], Vertex[j_][to_], Field[k_]] :> Rule[{from, to}, Field[k]], {0, Infinity}];
           internalEdges = internalEdges /. (Rule[{i_, j_}, f_] /; i > j) :> Rule[{j, i}, SARAH`AntiField[f]];
           Flatten @ (MapIndexed[Rule[{#[[1,1]], #[[1,2]], First[#2]}, #[[2]]]&, #]& /@ Gather[internalEdges, (First[#1] === First[#2])&])
          ];

GetCXXDiagramsDiagram[topology_] :=
    Module[{propagators, externalFields, externalFieldRules, vertices, externalVertices, internalVertices,
            adjacencyMatrix, internalFieldCouplings, internalEdges, internalFieldRules, resolvedCouplings},
           propagators = List @@ topology;
           externalFields = Flatten @ Join[Cases[propagators, Propagator[Incoming][Vertex[i_][from_], Vertex[j_][to_], Field[k_]] :> Rule[from, Field[k]], {0, Infinity}],
                                           Cases[propagators, Propagator[Outgoing][Vertex[i_][from_], Vertex[j_][to_], Field[k_]] :> Rule[from, Field[k]], {0, Infinity}]];
           externalVertices = externalFields[[All, 1]];
           externalFieldRules = Flatten @ ({{_,#,_} :> SARAH`AntiField[# /. externalFields],
                                            {#,_,_} :> SARAH`AntiField[# /. externalFields]} & /@ externalVertices);
           vertices = Table[i, {i, 1, GetNumberOfVertices[topology]}];
           internalVertices = Complement[vertices, externalVertices];
           adjacencyMatrix = GetUndirectedAdjacencyMatrix[topology];
           internalFieldCouplings = GetEdgeLists[adjacencyMatrix, internalVertices] /. externalFieldRules;
           internalEdges = GetInternalEdgeLists[topology];
           internalFieldRules = Join[internalEdges, internalEdges /. (Rule[{i_, j_, k_}, f_]) :> Rule[{j, i, k}, SARAH`AntiField[f]]];
           resolvedCouplings = internalFieldCouplings /. internalFieldRules;
           vertices /. externalFields /. Thread[Rule[internalVertices, resolvedCouplings]]
          ];

GetGenericDiagramClass[process_, {graphID_, topology_, insertions_, formFactors_}] :=
    Module[{graphName, adjacencyMatrix, edgeLabels, couplings,
            hasLocalTerms = FreeQ[formFactors, Finite]
            },
           graphName = CreateOneLoopDiagramName[process, graphID, insertions];
           adjacencyMatrix = GetUndirectedAdjacencyMatrix[topology];
           edgeLabels = GetEdgeLabels[topology];
           couplings = GetCouplings[formFactors];
           cxxDiagramsFormat = GetCXXDiagramsDiagram[topology];
           {graphName, adjacencyMatrix, edgeLabels, List @@ insertions,
            couplings, cxxDiagramsFormat, hasLocalTerms}
          ];

GetGenericGraphClassesFileName[] :=
    "generic_loop_decay_diagram_classes.m"

WriteGraphClassesFile[fileName_String, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

genericProcesses = {
    {S} -> {S, S},
    {S} -> {V, V},
    {S} -> {S, V},
    {S} -> {F, F}
};

genericOneLoopDiagramEvaluatorDecls = "";
genericOneLoopDiagramEvaluatorDefs = "";
genericOneLoopDiagramClasses = {};

For[i = 1, i <= Length[genericProcesses], i++,
    process = genericProcesses[[i]];
    Print["Creating C++ code for process: ", process, " ..."];
    Print["... reading input files"];
    diagramExprs = ReadFormFactorsExprs[process, resultsDir];
    diagramExprs = Select[diagramExprs, IsNonZeroDiagram];
    Print["... generating evaluation functions"];
    {decls, defs} = CreateDiagramEvaluators[process, diagramExprs];
    genericOneLoopDiagramEvaluatorDecls = genericOneLoopDiagramEvaluatorDecls <>
                                          If[genericOneLoopDiagramEvaluatorDecls != "", "\n\n", ""] <> decls;
    genericOneLoopDiagramEvaluatorDefs = genericOneLoopDiagramEvaluatorDefs <>
                                         If[genericOneLoopDiagramEvaluatorDefs != "", "\n\n", ""] <> defs;
    genericOneLoopDiagramClasses = Join[genericOneLoopDiagramClasses, GetGenericDiagramClass[process, #]& /@ diagramExprs];
   ];

Print["Writing output files ..."];
decayAmplitudesFiles = {{FileNameJoin[{templatesDir, "one_loop_decay_diagrams.hpp.in"}],
                         FileNameJoin[{resultsDir, "one_loop_decay_diagrams.hpp"}]},
                        {FileNameJoin[{templatesDir, "one_loop_decay_diagrams.cpp.in"}],
                         FileNameJoin[{resultsDir, "one_loop_decay_diagrams.cpp"}]}};
WriteOut`ReplaceInFiles[decayAmplitudesFiles,
                        { "@genericOneLoopDiagramEvaluatorDecls@" -> TextFormatting`WrapLines[genericOneLoopDiagramEvaluatorDecls],
                          "@genericOneLoopDiagramEvaluatorDefs@"  -> TextFormatting`WrapLines[genericOneLoopDiagramEvaluatorDefs]
                        }];

Print["Generating C++ code finished"];

Print["Writing generic graph classes to file ..."];
propertiesFileName = FileNameJoin[{resultsDir, GetGenericGraphClassesFileName[]}];
genericOneLoopDiagramClasses = DeleteDuplicates[genericOneLoopDiagramClasses];
classesOutputStatus = WriteGraphClassesFile[propertiesFileName, genericOneLoopDiagramClasses];
If[classesOutputStatus === $Failed,
   Quit[2];
  ];
