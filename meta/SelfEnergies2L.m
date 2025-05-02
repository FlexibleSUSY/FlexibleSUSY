BeginPackage["SelfEnergies2L`", {"SARAH`"}];

ConvertSarah2LDiagramList::usage = "Converts SARAH's list with 2-loop
 self-energy and tadpole diagrams."

Begin["`Private`"];

AppendIndex[p_, idx_, 1] := p;
AppendIndex[p_, idx_, range_] := AppendIndex[p, idx];
AppendIndex[SARAH`bar[p_], idx_] := SARAH`bar[AppendIndex[p, idx]];
AppendIndex[Susyno`LieGroups`conj[p_], idx_] := Susyno`LieGroups`conj[AppendIndex[p, idx]];
AppendIndex[p_, idx_] := p[{idx}];

DistChir[{}, {}] := 1;
DistChir[{0, rest1___}, {c_, rest2___}] := c     DistributeChiralities[{rest1}, {rest2}];
DistChir[{L, rest1___}, {c_, rest2___}] := c[PL] DistributeChiralities[{rest1}, {rest2}];
DistChir[{R, rest1___}, {c_, rest2___}] := c[PR] DistributeChiralities[{rest1}, {rest2}];

DistributeChiralities[chiralities_List, couplings_List] :=
    If[chiralities === (chiralities /. {L -> R, R -> L}),
       DistChir[chiralities, couplings],
       DistChir[chiralities, couplings] + (DistChir[chiralities /. {L -> R, R -> L}, couplings])
      ];

MultLF[{factor_, func_, chiralities_List}, particles_List, expr_] := 
    factor func[Sequence @@ (
        If[#3 === 1, SARAH`Mass2[#1], 
           SARAH`Mass2[#1, #2]]& @@@ particles
    )] DistributeChiralities[chiralities, expr];

MultiplyLoopFunction[loopfuncs_List, particles_List, expr_] := 
    Total[MultLF[#, particles, expr]& /@ loopfuncs];

(* check for ambiguous contraction of indices *)
IsAmbiguousContraction[indices_List] := 
    Or @@ ((#[[2]] > 2)& /@
           Tally[First /@ (DeleteCases[indices, {_, SARAH`gE1 | SARAH`gE2}] /.
                           SARAH`bar -> Identity /. Susyno`LieGroups`conj -> Identity)]);

(* one or more fields with indices appear more once *)
IsAmbiguousIndex[indices_List] := 
    Or @@ ((#[[2]] > 1)& /@ Tally[First /@ indices]);

DistributeIndices[{}, coupling_] := coupling;
DistributeIndices[indices_List, c : C[particles__]] :=
    If[IsAmbiguousIndex[indices],
       (* replace indices in order from indices list *)
       Assert[First /@ indices === {particles}];
       C[Sequence @@ 
         MapThread[(#1 /. #2)&, {{particles}, (Rule[#1, AppendIndex[#1, #2]]& @@@ indices)}]]
       ,
       (* use one common replacement rule for all fields *)
       c /. (Rule[#1, AppendIndex[#1, #2]]& @@@ indices)
    ];

(* multiplies all given couplings *)
MultiplyCouplings[lst_List] := DistributeIndices @@@ lst;

SumOverIndices[{}, expr_] := expr;
SumOverIndices[{{_, _, 1}, rest___}, expr_] :=
    SumOverIndices[{rest}, expr];
SumOverIndices[{{_, idx_, range_}, rest___}, expr_] :=
    SumOverIndices[{rest}, SARAH`sum[idx, 1, range, expr]];

(* combines list of particles/loop functions/couplings to a sum *)
CreateTadpoleDiag[{particles_List, loopFuncs_List, couplings_List}] :=
    SumOverIndices[particles, 
                   MultiplyLoopFunction[loopFuncs, particles, 
                                        MultiplyCouplings[couplings]]];

(* sums all diagrams that for a given type *)
SumTadpoleType[{type_, diags_List}] := Total[CreateTadpoleDiag /@ diags];

ConvertSarah2LDiagramList[tad_List, head_:Total] :=
    head[SumTadpoleType /@ tad] //. {
        (m : (SARAH`Mass | SARAH`Mass2))[(SARAH`bar | Susyno`LieGroups`conj)[p_], idx___] :> m[p, idx],
        (m : (SARAH`Mass | SARAH`Mass2))[p_, idx__] :> m[p[{idx}]],
        C[p__] :> Cp[p]
        (* unrotate to interaction basis *)
        (* SARAH`HiggsBoson[{SARAH`gE1}]   -> Symbol["U" <> ToString[SARAH`HiggsBoson]][{SARAH`gO1}], *)
        (* SARAH`HiggsBoson[{SARAH`gE2}]   -> Symbol["U" <> ToString[SARAH`HiggsBoson]][{SARAH`gO2}], *)
        (* SARAH`PseudoScalar[{SARAH`gE1}] -> Symbol["U" <> ToString[SARAH`PseudoScalar]][{SARAH`gO1}], *)
        (* SARAH`PseudoScalar[{SARAH`gE2}] -> Symbol["U" <> ToString[SARAH`PseudoScalar]][{SARAH`gO2}] *)
    };

End[];

EndPackage[];
