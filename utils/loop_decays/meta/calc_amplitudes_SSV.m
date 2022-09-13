<< "meta/OneLoopDecaysUtils.m";

status = 0;

baseDir = DirectoryName[$InputFileName];
AppendTo[$Path, baseDir];

CheckDirectoryContains[dir_String, files_List] :=
    Module[{workDir = Directory[], result = False},
           If[DirectoryQ[dir],
              result = FileNames[files, {dir}] =!= {};
             ];
           result
          ];

ToFieldString[field_[indices__Integer]] := ToString[field] <> StringJoin[ToString /@ List[indices]];
ToFieldString[field_] := ToString[field];

ToGenericFieldString[field_[indices__Integer]] := ToFieldString[field];
ToGenericFieldString[field_] := ToFieldString[field];

CreateDiagramsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "diagrams_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesExprsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_exprs_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateFormFactorsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "form_factors_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

loadFeynArts = Needs["FeynArts`"];
loadFormCalc = Needs["FormCalc`"];

If[loadFeynArts === $Failed || loadFormCalc === $Failed,
   Quit[1];
  ];
OneLoopDecaysUtils`CheckFormCalcVersion[$FormCalcVersion];

WriteFeynArtsOutputFile[fileName_, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated with " <> $FeynArtsVersion <>
                        " at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

WriteFormCalcOutputFile[fileName_, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated with " <> $FormCalcVersion <>
                        " at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

WriteFormFactorsOutputFile[fileName_, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

CalculateAmplitudes[amplitudeHead_, amplitudeExpr_] :=
    (
     FormCalc`ClearProcess[];
     FormCalc`CalcFeynAmp[amplitudeHead[amplitudeExpr], OnShell -> False] //. Subexpr[] //. Abbr[] //. GenericList[]
    );

(* Calculate all contributing graphs using FeynArts/FormCalc *)
topologies = FeynArts`CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Tadpoles];

process = {S} -> {S, V};

diagramsOutputFile = CreateDiagramsOutputFileName[process];
amplitudesOutputFile = CreateAmplitudesOutputFileName[process];
amplitudesExprsOutputFile = CreateAmplitudesExprsOutputFileName[process];
formFactorsOutputFile = CreateFormFactorsOutputFileName[process];

diags = FeynArts`InsertFields[topologies, process, InsertionLevel -> Generic, Model -> "SM"];
diagsOutputStatus = WriteFeynArtsOutputFile[FileNameJoin[{resultsDir, diagramsOutputFile}], diags];
If[diagsOutputStatus === $Failed,
   status = 2;
  ];

amplitudes = FeynArts`CreateFeynAmp[diags];
ampsOutputStatus = WriteFeynArtsOutputFile[FileNameJoin[{resultsDir, amplitudesOutputFile}], amplitudes];
If[ampsOutputStatus === $Failed,
   status = 2;
  ];

GetGraphID[FeynAmp[GraphID[details__], rest__]] := GraphID[details];
graphIDs = List @@ (GetGraphID /@ amplitudes);

genericAmplitudes = FeynArts`PickLevel[Generic][amplitudes];

amplitudesExprs = CalculateAmplitudes[Head[genericAmplitudes], #]& /@ genericAmplitudes;
exprsOutputStatus = WriteFormCalcOutputFile[FileNameJoin[{resultsDir, amplitudesExprsOutputFile}], amplitudesExprs];
If[exprsOutputStatus === $Failed,
   status = 2;
  ];

loadUtils = Needs["OneLoopDecaysUtils`"];
If[loadUtils === $Failed,
   Quit[1];
  ];

(* Given a list of matrix element and coefficients, apply simplifications assumed
   to hold in SARAH and FlexibleSUSY
*)
CanonicalizeCouplings[formFactors_List] :=
    Module[{countGhosts, countVectors, couplingsUUV, dupCouplingsUUV,
            couplingSubs, result},
           result = formFactors;
           countGhosts[fields_List] := Count[fields, U[__] | -U[__]];
           countVectors[fields_List] := Count[fields, V[__] | -V[__]];
           couplingsUUV = Cases[formFactors, SARAH`Cp[fields__][lor_] /; (countGhosts[List[fields]] == 2 &&
                                                                          countVectors[List[fields]] == 1), {0, Infinity}];
           dupCouplingsUUV = Select[GatherBy[couplingsUUV, #[[0]]&], (Length[#] > 1) &];
           (* Assume that only one component of the kinematic vector for a given U-U-V
              coupling has non-vanishing coefficient *)
           If[dupCouplingsUUV =!= {},
              couplingSubs = DeleteCases[#, SARAH`Cp[U[a___, Index[Generic, i_], b___],
                                                     U[c___, Index[Generic, j_], d___],
                                                     V[x___, Index[Generic, k_], y___]][SARAH`Mom[U[a___, Index[Generic, i_], b___]]]]& /@ dupCouplingsUUV;
              couplingSubs = Rule[#, 0]& /@ Flatten[couplingSubs];
              result = result /. couplingSubs;
             ];
           result
          ];

(* Converts a FeynmanGraph object to a list containing the graph number,
   combinatorial factors, and insertions *)
ConvertFeynmanGraphToList[graph_] :=
    {OneLoopDecaysUtils`GetGraphNumber[graph],
     OneLoopDecaysUtils`GetGraphCombinatorialFactor[graph],
     OneLoopDecaysUtils`GetGraphInsertions[graph]
    };

(* Combines the lists of graph IDs, diagrams, and amplitudes (expressed as a list of matrix
   elements and corresponding coefficients) into a single list where each element is of the form
   {graph ID, topology, insertions, form factors}
*)
CollectDiagramInfo[ids_, diagrams_, formFactors_] :=
    Module[{i, j, nTopologies = Length[diagrams],
            topology, insertions, nGenericInsertions, genericInsertions,
            genericAmps, genericIDs, count = 1},
           First[Last[
               Reap[
                    For[i = 1, i <= nTopologies, i++,
                        topology = diagrams[[i, 1]];
                        insertions = diagrams[[i, 2]];
                        genericInsertions = List @@ insertions;
                        nGenericInsertions = Length[genericInsertions];
                        genericAmps = List @@ formFactors[[count ;; count + nGenericInsertions - 1]];
                        genericIDs = ids[[count ;; count + nGenericInsertions - 1]];
                        count += nGenericInsertions;
                        MapThread[Sow[List[#1, topology, #2, Simplify[#3]]]&, {genericIDs, genericInsertions, genericAmps}];
                       ];
                   ]]]
          ];

Print["Extracting form factors ..."];

formFactors = OneLoopDecaysUtils`ExtractFormFactors /@ amplitudesExprs;

Print["Converting form factors ..."];

formFactors = OneLoopDecaysUtils`ToFSConventions /@ formFactors;
formFactors = CanonicalizeCouplings /@ formFactors;

Print["Combining graph info ..."];

contributions = CollectDiagramInfo[graphIDs, diags, formFactors];

formFactorsOutputStatus = WriteFormFactorsOutputFile[FileNameJoin[{resultsDir, formFactorsOutputFile}], contributions];
If[formFactorsOutputStatus === $Failed,
   status = 3;
  ];
