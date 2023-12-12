FlexibleSUSY`WriteClass[obs:FlexibleSUSYObservable`BrDLDL, slha_, files_] :=
Module[
   {
      observables = DeleteDuplicates@Cases[Observables`GetRequestedObservables@slha, _obs],
      prototypes = "", definitions = "", npfHeaders = "", npfDefinitions = "",
      fieldsFFV = {}, verticesFFV = {}, additionalVertices = {}
   },

   If[observables =!= {} && FlexibleSUSY`FSFeynArtsAvailable && FlexibleSUSY`FSFormCalcAvailable,
      Utils`PrintHeadline["Creating " <> SymbolName@obs <> " class ..."];
      observables = DeleteDuplicates[observables /. f_@_Integer -> f@_];
      fieldsFFV   = Cases[observables, Rule[{q_[_], _}, _] :> q, Infinity];
      verticesFFV = Cases[observables, Rule[{_, l_[_]}, _] :> {SARAH`bar@l, l, TreeMasses`GetPhoton[]}, Infinity];
      {additionalVertices, npfDefinitions, prototypes, definitions} = create@observables;
      npfHeaders  = NPointFunctions`CreateCXXHeaders[];
   ];

   WriteOut`ReplaceInFiles[
      files,
      {
         "@calculate_prototypes@"        -> prototypes,
         "@calculate_definitions@"       -> definitions,
         "@npointfunctions_headers@"     -> npfHeaders,
         "@npointfunctions_definitions@" -> npfDefinitions,
         "@include_guard@"               -> SymbolName@obs,
         "@namespace@"                   -> Observables`GetObservableNamespace@obs,
         "@filename@"                    -> Observables`GetObservableFileName@obs,
         Sequence@@FlexibleSUSY`Private`GeneralReplacementRules[]
      }
   ];

   {
      "FFV fields" -> DeleteDuplicates@fieldsFFV,
      "C++ vertices" -> DeleteDuplicates@Join[verticesFFV, additionalVertices]
   }
];

create[manyObservables:{__FlexibleSUSYObservable`BrDLDL}] := {
   DeleteDuplicates[Join@@#[[All,1]]],
   StringRiffle[#[[All, 2]], "\n\n"],
   StringRiffle[#[[All, 3]], "\n\n"],
   StringRiffle[#[[All, 4]], "\n\n"]
}&[create/@DeleteDuplicates[manyObservables /. Rule[_Integer, _Integer] -> Rule[_, _]]];

create[obs:FlexibleSUSYObservable`BrDLDL[{qd_[_], lep_[_]} -> _, contr_, loopN_]] :=
Module[{npfVertices, npfDefinition, calculateDefinition, prototype,
   keep, npf, fields, sp, dc, basis},
   keep = If[Head@contr === List, contr, {contr}];
   Print["Contributions: ", SymbolName/@keep];
   Print["Loop level: ", loopN];
   npf = NPointFunctions`NPointFunction[
      {qd, lep},
      {qd, lep},
      NPointFunctions`UseCache            -> True,
      NPointFunctions`OnShellFlag         -> True,
      NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
      NPointFunctions`LoopLevel           -> loopN,
      NPointFunctions`Regularize          -> FlexibleSUSY`FSRenormalizationScheme,
      NPointFunctions`Observable          -> obs,
      NPointFunctions`KeepProcesses       -> keep
   ];
   npf = npf /. {
      SARAH`sum[__] -> 0,
      LoopTools`B0i[i_, _, mm__] :> LoopTools`B0i[i, 0, mm],
      LoopTools`C0i[i_, Repeated[_, {3}], mm__] :> LoopTools`C0i[i, Sequence@@Array[0&, 3], mm],
      LoopTools`D0i[i_, Repeated[_, {6}], mm__] :> LoopTools`D0i[i, Sequence@@Array[0&, 6], mm]
   };
   fields = Flatten@NPointFunctions`GetProcess@npf;
   sp[i_] := SARAH`DiracSpinor[fields[[i]], 0, 0];
   dc[a_, b__, c_] := NPointFunctions`DiracChain[sp@a, b, sp@c];
   basis = With[{l = SARAH`Lorentz, R = 6, L = 7},
      {
         "S_LL" -> dc[3,L,1] dc[4,L,2],
         "S_LR" -> dc[3,L,1] dc[4,R,2],
         "S_RL" -> dc[3,R,1] dc[4,L,2],
         "S_RR" -> dc[3,R,1] dc[4,R,2],
         "V_LL" -> dc[3,R,l@1,1] dc[4,R,l@1,2],
         "V_LR" -> dc[3,R,l@1,1] dc[4,L,l@1,2],
         "V_RL" -> dc[3,L,l@1,1] dc[4,R,l@1,2],
         "V_RR" -> dc[3,L,l@1,1] dc[4,L,l@1,2],
         "T_LL" -> dc[3,-L,l@1,l@2,1] dc[4,-L,l@1,l@2,2],
         "T_RR" -> dc[3,-R,l@1,l@2,1] dc[4,-R,l@1,l@2,2]
      }
   ];
   npf = WilsonCoeffs`InterfaceToMatching[npf, basis];
   npfVertices = NPointFunctions`VerticesForNPointFunction@npf;
   npfDefinition = NPointFunctions`CreateCXXFunctions[npf, "@class@", SARAH`Delta, basis][[2]];
   prototype = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " " <> Observables`GetObservablePrototype@obs;

   {calculateDefinition, npfDefinition} = FixedPoint[
      StringReplace[#,
         {
            "@Qd@"    -> CConversion`ToValidCSymbolString@qd,
            "@Lep@"   -> CConversion`ToValidCSymbolString@lep,
            "@Ph@"    -> CConversion`ToValidCSymbolString@SARAH`Photon,
            "@loopN@" -> CConversion`ToValidCSymbolString@loopN,
            "@contr@" -> CConversion`ToValidCSymbolString@contr,
            "@class@" -> "@Qd@@Lep@_to_@Qd@@Lep@_@contr@@loopN@",
            "@photon_penguin@" -> If[loopN === 1, "calculate_@Qd@_@Qd@_@Ph@_form_factors", "zero"]
         }
      ]&,
      {
         prototype <> " {
         return forge<
            fields::@Qd@, fields::@Lep@, fields::@Ph@,
            @photon_penguin@,
            npointfunctions::@class@
         >(qi, qo, li, model, qedqcd);\n}",
         npfDefinition
      }
   ];

   {npfVertices, npfDefinition, prototype <> ";", calculateDefinition}
];
