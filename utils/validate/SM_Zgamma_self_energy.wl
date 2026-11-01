args=ToExpression/@Rest[$ScriptCommandLine];

<<FeynArts`
<<FormCalc`
Install["/home/wojciech/HEP-software/LoopTools-2.16/bin/LoopTools"];

Model = SM;


If[Length[args] != 15,
   Print["\033[31m\nWARNING: Wrong number of Arguments, this script needs 15 Parameters!\nUsing Fallback Parameters\033[0m\n"];

   scaleSM = 91.18760000000000332;

   (* Physical Masses *)
   mHSM=123;
   mZSM=91.43879126614187669;
   mWSM=80.173623954610150122;
   mLSM={0.00051099999999999995252,0.10499999999999998224,1.7769900000000002915};

   (* SM Tree-Level Masses *)
   mHSMTree=123;
   mZSMTree=91.43879126614187669;
   mWSMTree=80.173623954610150122;
   mLSMTree={0.00051099999999999995252,0.10499999999999998224,1.7769900000000002915};
   mQSMTree={0.17394826817189068535,0.17394826817189068535,165};
   mDSMTree={0.17394826817189068535,0.17394826817189068535,2.8999999999999999112};

   g1U1SUSY=0.4614901524829831958;
   g2SU2SUSY=0.65181808093178983388;
   vevSM = 246;,

   scaleSM=args[[1]];

   (* SM Physical Masses *)
   mHSM=args[[2]];
   mZSM=args[[3]];
   mWSM=args[[4]];
   mLSM=args[[5]];

   (* SM Tree-Level Masses *)
   mHSMTree=args[[6]];
   mZSMTree=args[[7]];
   mWSMTree=args[[8]];
   mLSMTree=args[[9]];
   mvSMTree=args[[10]];
   mQSMTree=args[[11]];
   mDSMTree=args[[12]];

   g1U1SUSY=args[[13]];
   g2SU2SUSY=args[[14]];
   vevSM=args[[15]];
]

g1U1SM=Sqrt[3/5]*g1U1SUSY;
g2SU2SM=g2SU2SUSY;
\[Theta]WSM=ArcCos[mWSMTree/mZSMTree];

mFSMTree = {mLSMTree, mQSMTree, mDSMTree};

SetMudim[scaleSM^2];

(* === Basic Definitions === *)

selfEnergyTop=CreateTopologies[1,1->1,ExcludeTopologies->Reducible];

PaVeRe[expr_] := expr/.A0i[f_, m_] :> Re[A0i[f, m]]/.B0i[f_, k_, m1_, m2_] :> Re[B0i[f, k, m1, m2]];

ClearAll[CleanNumber];
CleanNumber[num_, prec_: 16, thresh_: 10^-16] := Module[
  {preciseNum},

  preciseNum = SetPrecision[num, prec];
  Chop[preciseNum, thresh]
];

(* === Higgs Self-Energy === *)

higgsSelfEnergyDiagram=InsertFields[selfEnergyTop,
                                    S[1]->S[1], 
                                    Model->Model,
                                    InsertionLevel->Classes];

higgsSelfEnergyAmp = CreateFeynAmp[higgsSelfEnergyDiagram,
                                   PreFactor->1 / (16 \[Pi]^4),
                                   Truncated -> True];

higgsSelfEnergyFermion = CalcFeynAmp[higgsSelfEnergyAmp, 
                                     OnShell -> False][[1]]//.SubExpr[]//.Abbr[]/.{
                                     Alfa -> g2^2 * SW2 / (4 * \[Pi])};
higgsSelfEnergyBoson = CalcFeynAmp[higgsSelfEnergyAmp,
                                   OnShell -> False][[2]]//.SubExpr[]//.Abbr[]/.{
                                   Finite->1,
                                   Alfa -> g2^2 * SW2 / (4 * \[Pi])};

higgsSelfEnergyFermionFunc[p_,mf_,mVW_,g1_,g2_,T_] = higgsSelfEnergyFermion/.{
                                                         Mf2[a_,b_] :> (mf[[a - 1, b]])^2, 
                                                         MW2 -> mVW^2,
                                                         SW2 -> Sin[T]^2,
                                                         CW2 -> Cos[T]^2,
                                                         Pair[__] -> p^2}//Quiet;
higgsSelfEnergyBosonFunc[p_, mh_, mVZ_, mVW_, g1_, g2_, T_] = higgsSelfEnergyBoson/.{
      MH2 -> mh^2, 
      MZ2 -> mVZ^2, 
      MW2 -> mVW^2, 
      SW2 -> Sin[T]^2, 
      CW2 -> Cos[T]^2, 
      Pair[__] -> p^2};

higgsSelfEnergyNum[p_] = Quiet[-I * (ExpandSums[higgsSelfEnergyFermionFunc[p, mFSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM]] + higgsSelfEnergyBosonFunc[p, mHSMTree, mZSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM])];

higgsSelfEnergyDerivNum[p_] = D[higgsSelfEnergyNum[k]/.{k^2 -> a}, a]/.{a -> p^2};

(* === Lepton Self Energy === *)

leptonSelfEnergyDiagram = InsertFields[selfEnergyTop, 
                                       F[2, {i}] -> F[2, {j}],
                                       Model -> Model,
                                       InsertionLevel -> {Classes}];

leptonSelfEnergyAmp = CreateFeynAmp[leptonSelfEnergyDiagram,
                                    PreFactor -> 1 / (16 * \[Pi]^4),
                                    Truncated -> True]

ClearProcess[];
leptonSelfEnergy = CalcFeynAmp[leptonSelfEnergyAmp,
                               OnShell -> False][[1]]//.SubExpr[]//.Abbr[]/.{
                               Finite -> 1,
                               Alfa -> g2^2 * SW2 / (4 * \[Pi])};

leptonSelfEnergyFunc[p_, mh_, mf_, mVZ_, mVW_, g1_, g2_, T_] = leptonSelfEnergy/.{
      MH2 -> mh^2,
      Mf[a_,b_] :> mf[[a - 1,b]],
      Mf2[a_,b_] :> mf[[a - 1,b]]^2,
      MZ2 -> mVZ^2, 
      MW2 -> mVW^2, 
      SW2 -> Sin[T]^2,
      CW2 -> Cos[T]^2,
      Pair[__] -> p^2}//Quiet;

leptonSelfEnergyNum[p_, i_, j_] = -I * ExpandSums[leptonSelfEnergyFunc[p , mHSMTree, mFSMTree, mZSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM]]//Expand//Quiet;

leptonSelfEnergyPLNum[p_, i_, j_] = Coefficient[leptonSelfEnergyNum[p, i, j], DiracChain[6, k[1]]]//Quiet;
leptonSelfEnergyPRNum[p_, i_, j_] = Coefficient[leptonSelfEnergyNum[p, i, j], DiracChain[7, k[1]]]//Quiet;
leptonSelfEnergymLNum[p_, i_, j_] = Coefficient[leptonSelfEnergyNum[p, i, j], DiracChain[7]]//Quiet;
leptonSelfEnergymRNum[p_, i_, j_] = Coefficient[leptonSelfEnergyNum[p, i, j], DiracChain[6]]//Quiet;

leptonSelfEnergyDerivPLNum[p_, i_, j_] = Coefficient[D[leptonSelfEnergyNum[k, i, j]/.{k^2 -> a}, a]/.{a -> p^2}, DiracChain[6, k[1]]]//Quiet;
leptonSelfEnergyDerivPRNum[p_, i_, j_] = Coefficient[D[leptonSelfEnergyNum[k, i, j]/.{k^2 -> a}, a]/.{a -> p^2}, DiracChain[7, k[1]]]//Quiet;
leptonSelfEnergyDerivmLNum[p_, i_, j_] = Coefficient[D[leptonSelfEnergyNum[k, i, j]/.{k^2 -> a}, a]/.{a -> p^2}, DiracChain[7]]//Quiet;
leptonSelfEnergyDerivmRNum[p_, i_, j_] = Coefficient[D[leptonSelfEnergyNum[k, i, j]/.{k^2 -> a}, a]/.{a -> p^2}, DiracChain[6]]//Quiet;

(* === Z Self Energy === *)

gaugeZSelfEnergyDiagram = InsertFields[selfEnergyTop,
                                       V[2] -> V[2],
                                       Model -> Model,
                                       InsertionLevel -> {Classes}];

gaugeZSelfEnergyAmp = CreateFeynAmp[gaugeZSelfEnergyDiagram,
                                    PreFactor -> 1 / (16 \[Pi]^4),
                                    Truncated -> True];

ClearProcess[];
gaugeZSelfEnergyLepton = CalcFeynAmp[gaugeZSelfEnergyAmp,
                                     OnShell -> False][[1]]//.SubExpr[]//.Abbr[]/.{
                                       Alfa -> g2^2 * SW2 / (4 * \[Pi])};
gaugeZSelfEnergyBoson = CalcFeynAmp[gaugeZSelfEnergyAmp,
                                    OnShell -> False][[2]]//.SubExpr[]//.Abbr[]/.{
                                       Finite -> 1,
                                       Alfa -> g2^2 * SW2 / (4 * \[Pi])};

gaugeZSelfEnergyLeptonFunc[p_, mf_, g1_, g2_, T_] = gaugeZSelfEnergyLepton/.{
                                    Mf2[a_,b_] :> mf[[a - 1,b]]^2,
                                    SW2 -> Sin[T]^2,
                                    CW2 -> Cos[T]^2,
                                    Pair[__] -> p^2}/.{
                                    k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                    MetricTensor[Lor1, Lor2] -> -B + A}//Quiet;
gaugeZSelfEnergyBosonFunc[p_, mh_, mVZ_, mVW_, g1_, g2_, T_] = gaugeZSelfEnergyBoson/.{
                                    MH2 -> mh^2,
                                    MZ2 -> mVZ^2,
                                    MW2 -> mVW^2,
                                    SW2 -> Sin[T]^2,
                                    CW2 -> Cos[T]^2,
                                    Pair[__] -> p^2}/.{
                                    k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                    MetricTensor[Lor1, Lor2] -> -B + A};

gaugeZSelfEnergyNum[p_] = I * ExpandSums[gaugeZSelfEnergyLeptonFunc[p, mFSMTree, g1U1SM, g2SU2SM, \[Theta]WSM] + gaugeZSelfEnergyBosonFunc[p, mHSMTree, mZSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM]]//Expand//Quiet

gaugeZSelfEnergyDerivNum[p_] = D[gaugeZSelfEnergyNum[k]/.{k^2 -> a}, a]/.{a -> p^2}

(* === W Self Energy === *)

gaugeWSelfEnergyDiagram = InsertFields[selfEnergyTop,
                                       V[3] -> V[3],
                                       Model -> Model,
                                       InsertionLevel -> {Classes}];

gaugeWSelfEnergyAmp = CreateFeynAmp[gaugeWSelfEnergyDiagram,
                                    PreFactor -> 1 / (16 \[Pi]^4),
                                    Truncated -> True];

ClearProcess[];
gaugeWSelfEnergyLepton = CalcFeynAmp[gaugeWSelfEnergyAmp,
                                     OnShell -> False][[1]]//.SubExpr[]//.Abbr[]/. {
                                    Alfa -> g2^2 * SW2 / (4 * \[Pi])};
gaugeWSelfEnergyBoson = CalcFeynAmp[gaugeWSelfEnergyAmp,
                                    OnShell -> False][[2]]//.SubExpr[]//.Abbr[]/.{
                                    Finite -> 1,
                                    Alfa -> g2^2 * SW2 / (4 * \[Pi])};


gaugeWSelfEnergyLeptonFunc[p_, mf_, g1_, g2_, T_] = gaugeWSelfEnergyLepton/.{
                                    Mf2[a_,b_] :> mf[[a - 1, b]]^2,
                                    SW2 -> Sin[T]^2,
                                    CW2 -> Cos[T]^2,
                                    Pair[__] -> p^2}/.{
                                    k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                    MetricTensor[Lor1, Lor2] -> -B + A}//Quiet;

gaugeWSelfEnergyBosonFunc[p_, mh_, mVZ_, mVW_, g1_, g2_, T_] = gaugeWSelfEnergyBoson/.{
                                    MH2 -> mh^2,
                                    MZ2 -> mVZ^2,
                                    MW2 -> mVW^2,
                                    SW2 -> Sin[T]^2,
                                    CW2 -> Cos[T]^2,
                                    Pair[__] -> p^2}/.{
                                    k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                    MetricTensor[Lor1, Lor2] -> -B + A};

gaugeWSelfEnergyNum[p_] = I * (ExpandSums[gaugeWSelfEnergyLeptonFunc[p, mFSMTree, g1U1SM, g2SU2SM, \[Theta]WSM]] + gaugeWSelfEnergyBosonFunc[p, mHSMTree, mZSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM])//Expand//Quiet;

gaugeWSelfEnergyDerivNum[p_] = D[gaugeWSelfEnergyNum[k]/.{k^2 -> a}, a]/.{a -> p^2}; 

(* === P Self Energy === *)

gaugePSelfEnergyDiagram = InsertFields[selfEnergyTop,
                                       V[1] -> V[1],
                                       Model -> Model,
                                       InsertionLevel -> {Classes}];
gaugePSelfEnergyAmp = CreateFeynAmp[gaugePSelfEnergyDiagram,
                                    PreFactor -> 1 / (16 \[Pi]^4),
                                    Truncated -> True];

ClearProcess[];
gaugePSelfEnergyLepton = CalcFeynAmp[gaugePSelfEnergyAmp,
                                     OnShell -> False][[1]]//.SubExpr[]//.Abbr[]/. {
                                    Alfa -> g2^2 * SW2 / (4 * \[Pi])};
gaugePSelfEnergyBoson = CalcFeynAmp[gaugePSelfEnergyAmp,
                                    OnShell -> False][[2]]//.SubExpr[]//.Abbr[]/.{Finite -> 1,
                                    Alfa -> g2^2 * SW2 / (4 * \[Pi])};

gaugePSelfEnergyLeptonFunc[p_, mf_, g1_, g2_, T_] = gaugePSelfEnergyLepton/.{
                                    Mf2[a_, b_] :> mf[[a - 1, b]]^2,
                                    SW2 -> Sin[T]^2,
                                    CW2 -> Cos[T]^2,
                                    Pair[__] -> p^2}/.{
                                    k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                    MetricTensor[Lor1, Lor2] -> -B + A}//Quiet;
gaugePSelfEnergyBosonFunc[p_, mh_, mVZ_, mVW_, g1_, g2_, T_] = gaugePSelfEnergyBoson/.{
                                    MH2 -> mh^2,
                                    MZ2 -> mVZ^2,
                                    MW2 -> mVW^2,
                                    SW2 -> Sin[T]^2,
                                    CW2 -> Cos[T]^2,
                                    Pair[__] -> p^2}/.{
                                    k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                    MetricTensor[Lor1, Lor2] -> -B + A};

gaugePSelfEnergyNum[p_] = I * (ExpandSums[gaugePSelfEnergyLeptonFunc[p, mFSMTree, g1U1SM, g2SU2SM, \[Theta]WSM]] + gaugePSelfEnergyBosonFunc[p, mHSMTree, mZSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM])//Expand//Quiet;


gaugePSelfEnergyDerivNum[p_] = D[gaugePSelfEnergyNum[k]/.{k -> Sqrt[a]}, a]/.{a -> p^2};

(* === P-Z Self Energy === *)

gaugePZSelfEnergyDiagram = InsertFields[selfEnergyTop,
                                        V[1] -> V[2],
                                        Model -> Model,
                                        InsertionLevel -> {Classes}];
gaugePZSelfEnergyAmp = CreateFeynAmp[gaugePZSelfEnergyDiagram,
                                     PreFactor -> 1 / (16 \[Pi]^4),
                                     Truncated -> True];

ClearProcess[];
gaugePZSelfEnergyLepton = CalcFeynAmp[gaugePZSelfEnergyAmp,
                                      OnShell -> False,
                                      PaVeReduce -> False,
                                      SortDen -> False][[1]]//.SubExpr[]//.Abbr[]/.{
                                      Alfa -> g2^2 * SW2 / (4 * \[Pi])};
gaugePZSelfEnergyBoson = CalcFeynAmp[gaugePZSelfEnergyAmp,
                                     OnShell -> False,
                                     PaVeReduce -> False,
                                     SortDen -> False][[2]]//.SubExpr[]//.Abbr[]/.{
                                     Finite -> 1,
                                     Alfa -> g2^2 * SW2 / (4 * \[Pi])};


gaugePZSelfEnergyLeptonFunc[p_, mf_, g1_, g2_, T_] = gaugePZSelfEnergyLepton/.{
                                     Mf2[a_, b_] :> mf[[a - 1, b]]^2,
                                     SW -> Sin[T],
                                     CW -> Cos[T],
                                     SW2 -> Sin[T]^2,
                                     CW2 -> Cos[T]^2,
                                     Pair[__] -> p^2}/.{
                                     k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                     MetricTensor[Lor1, Lor2] -> -B + A}//Expand//Quiet;
gaugePZSelfEnergyBosonFunc[p_, mh_, mVZ_, mVW_, g1_, g2_, T_] = gaugePZSelfEnergyBoson/.{
                                     MH2 -> mh^2,
                                     MZ2 -> mVZ^2,
                                     MW2 -> mVW^2,
                                     SW -> Sin[T],
                                     CW -> Cos[T],
                                     SW2 -> Sin[T]^2,
                                     CW2 -> Cos[T]^2,
                                     Pair[__] -> p^2}/.{
                                     k[1][Lor1] k[1][Lor2] -> -B*p^2,
                                     MetricTensor[Lor1, Lor2] -> -B + A}//Expand;

gaugePZSelfEnergyNum[p_] = I * (ExpandSums[gaugePZSelfEnergyLeptonFunc[p, mFSMTree, g1U1SM, g2SU2SM, \[Theta]WSM]] + gaugePZSelfEnergyBosonFunc[p, mHSMTree, mZSMTree, mWSMTree, g1U1SM, g2SU2SM, \[Theta]WSM])//Expand//Quiet;

gaugePZSelfEnergyDerivNum[p_] = D[gaugePZSelfEnergyNum[k]/.{k -> Sqrt[a]}, a]/.{a -> p^2};


(* === Self Energy Printouts ===  *)

Print["Self-Energies"]
Print[CleanNumber[higgsSelfEnergyNum[mHSM],16]]
Print[CleanNumber[leptonSelfEnergyPLNum[mLSM[[2]], 2, 2] / 2, 16]]
Print[CleanNumber[leptonSelfEnergyPRNum[mLSM[[2]], 2, 2] / 2,16]]
Print[CleanNumber[leptonSelfEnergymLNum[mLSM[[2]], 2, 2],16]]
Print[CleanNumber[leptonSelfEnergymRNum[mLSM[[2]], 2, 2],16]]
Print[CleanNumber[gaugeZSelfEnergyNum[mZSM],16]]
Print[CleanNumber[gaugeWSelfEnergyNum[mWSM],16]]
Print[CleanNumber[gaugePSelfEnergyNum[0],16]]
Print[CleanNumber[gaugePZSelfEnergyNum[0],16]]

Print["\nSelf-Energy Derivatives"]
Print[CleanNumber[higgsSelfEnergyDerivNum[mHSM],16]]
Print[CleanNumber[leptonSelfEnergyDerivPLNum[mLSM[[2]], 2, 2] / 2, 16]]
Print[CleanNumber[leptonSelfEnergyDerivPRNum[mLSM[[2]], 2, 2] / 2,16]]
Print[CleanNumber[leptonSelfEnergyDerivmLNum[mLSM[[2]], 2, 2],16]]
Print[CleanNumber[leptonSelfEnergyDerivmRNum[mLSM[[2]], 2, 2],16]]
Print[CleanNumber[gaugeZSelfEnergyDerivNum[mZSM],16]]
Print[CleanNumber[gaugeWSelfEnergyDerivNum[mWSM],16]]
Print[CleanNumber[gaugePSelfEnergyDerivNum[0],16]]
Print[CleanNumber[gaugePZSelfEnergyDerivNum[0],16]]
