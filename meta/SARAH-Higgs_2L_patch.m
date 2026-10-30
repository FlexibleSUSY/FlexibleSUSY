Calc2LoopCorrections[states_]:=Block[{sstates},

Print["Calculating 2 loop corrections"];
sstates=ToString[states];
SA`CurrentStates = states;				     
StringReplaceFinal = {};
SPhenoParameters = {};
temp = Select[VertexListNonCC, (#[[-1]] === SSSS || #[[-1]] === SSVV) &];(*extract all 4-vertices with scalars*)
temp = Select[temp, ((Length[Intersection[RE /@ (#[[1, 1]] /. A_[{b__}] -> A)]] < 3) && (Mod[Count[RE /@ (#[[1, 1]] /. A_[{b__}] -> A), RE[(#[[1, 1, 1]] /. A_[{b__}] -> A)]], 2] =!= 1)) &];(*take only those vertices which have 2 pairs of \
		identical fields*)

CouplingUsedForEffPot = True;
dataUnBrokenGaugeGroups = {#, Gauge[[#, 3]], Gauge[[#, 4]], 
     SGauge[[#]] /. A_[{b__}] -> A, getDimFundamental[Gauge[[#, 2]]], 
     getDimAdjoint[
      Gauge[[#, 2]]]} & /@ (Position[SGauge /. A_[{b__}] -> A, #][[1, 
       1]] & /@ 
     Select[SGauge /. A_[{b__}] -> A, 
      FreeQ[Particles[SA`CurrentStates], #] == False &]);
				     
subfourpoint = {};
For[ii = 1, ii <= Length[dataUnBrokenGaugeGroups], ii++, 
  namestub = StringTake[ToString[dataUnBrokenGaugeGroups[[ii, 2]]], 1];
  AppendTo[subfourpoint, 
   ToExpression[namestub <> "t" <> "3"] -> 
    ToExpression[namestub <> "t" <> "1"]];
  AppendTo[subfourpoint, 
   ToExpression[namestub <> "t" <> "4"] -> 
    ToExpression[namestub <> "t" <> "2"]];];

listBrokenGaugeCouplings = BetaGauge[[All, 1]];
For[i = 1, i <= Length[dataUnBrokenGaugeGroups], i++, 
  listBrokenGaugeCouplings = 
   DeleteCases[listBrokenGaugeCouplings, 
    dataUnBrokenGaugeGroups[[i, 3]]]];

subZeroGaugeLess = 
  Table[listBrokenGaugeCouplings[[i]] -> 0, {i, 1, 
    Length[listBrokenGaugeCouplings]}];

AllRelevant = getAllRelevantCouplings[VertexListNonCC];
AllRelevant = 
  Select[AllRelevant /. subZeroGaugeLess, 
   If[Length[#[[1]]] === 
      3, (#[[1, 2, 1]] =!= 0) || (#[[1, 3, 1]] =!= 0), #[[1, 2, 1]] =!=
       0] &];
temp = SPhenoCouplingList[AllRelevant];
SPhenoCouplings3P = temp[[1]];

temp = Select[VertexListNonCC, (#[[-1]] === SSSS) &];
(*extract all 4-vertices with scalars only,drop the ones with vectors*)

temp = Select[temp /. subZeroGaugeLess, #[[1, 2, 1]] =!= 0 &];

specialPOLEvertices = {};

Block[{tt, mi, mj, mk, tttt2, tttt3, POLEstructures, inner, outer}, 
  specialPOLEvertices = 
   Reap[For[mi = 1, mi <= Length[temp], mi++, 
      Sow[Reap[
          Sow[{temp[[mi, 1, 
             1]]}];(*entries start with the particles in the coupling,
          then a list of group structures*)
          For[mj = 1, mj <= Length[dataUnBrokenGaugeGroups], mj++, 
           tt = ExtractStructure[temp[[mi, 1, 2, 1]], 
             dataUnBrokenGaugeGroups[[mj, 2]]];
           tttt2 = Select[tt, #[[2, 1]] =!= 0 &];
           
           POLEstructures = 
            Table[tttt2[[mk, 1]], {mk, 1, Length[tttt2]}];
           Sow[POLEstructures]
           (*AppendTo[specialPOLEvertices,{{temp[[i,1,1]],
           POLEstructures}}];*)];][[2, 
         1]](*end inner reap*)];] (*end for...*)][[2, 
    1]];(*end outer Reap*)
  specialPOLEverticesorg = 
   Table[C @@ specialPOLEvertices[[mi, 1, 1]] /. {A_[{___}] -> A}, {mi, 1, Length[specialPOLEvertices]}];];

CouplingsFor2LPole = True;

temp2 = SPhenoCouplingList4ptPOLE[temp];

CouplingUsedFor2LPole = False;
SPhenoCouplings4Pole = temp2[[1]];

$sarahCurrentOutputDir = ToFileName[{$sarahOutputDir, $sarahModelNameMain, sstates}];
$sarahCurrentSPhenoDir = ToFileName[{$sarahCurrentOutputDir,"SPheno"}];
$sarahStoreTwoLoopDir = ToFileName[{$sarahCurrentOutputDir,"Two-Loop"}];

If[FileExistsQ[$sarahStoreTwoLoopDir]=!=True,
CreateDirectory[$sarahStoreTwoLoopDir]];
				     
ListAllParticles = Particles[states];
AllScalars = Select[ListAllParticles, #[[4]] === S &];
AllFermions = Select[ListAllParticles, #[[4]] === F &];
AllDiracFermions = getDiracFermionList;
DiracList = {#[[1]], #[[6]]} & /@ AllDiracFermions;
generate2LPoleFunctions
];
