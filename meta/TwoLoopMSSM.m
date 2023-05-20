(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["TwoLoopMSSM`"];
EndPackage[];

GetMSSMCPEvenHiggsLoopMassMatrix::usage = "Returns the
 loop corrections to the CP-even Higgs mass matrix in the
 CP-conserving MSSM.  The loop-corrections are taken from
 arxiv:hep-ph/0105096.

Note: The return value contains the contributions from tadpole
 diagrams.

Note: The sign of the mu parameter in arxiv:hep-ph/0105096 is opposite
 to the one in SARAH.  In order to switch to the SARAH convention,
 pass parameters -> {signMu -> -1} to the function.

Usage: GetMSSMCPEvenHiggsLoopMassMatrix[
           loopOrder -> {1,1,1}, corrections -> {1}, parameters -> {}]

Parameters:

- loopOrder (optional): List of factors multiplied by each loop order.
  #1: tree-level
  #2: 1-loop level
  #3: 2-loop level
  (default: {1,1,1})

- corrections (optional): List of factors multiplied by each correction.
  (default: {1, 1})
  #1: alpha_t * alpha_s (arxiv:hep-ph/0105096)
  #2: alpha_t^2         (arxiv:hep-ph/0112177)

  The expressions for the O(alpha_t^2) corrections
  [arxiv:hep-ph/0112177] have been kindly provided by Pietro Slavich.

- parameters (optional): List of internal replacement rules for
  parameters, useful when certain limits are considered (default: {}).
  Example:  parameters -> {At -> 0, sin2Theta -> 0}

- abbreviations (optional): List of replacement rules to resolve
  abbreviated symbols and functions (default:
  {TwoLoopMSSM`Private`alphaT2Abbreviations}).
  Example: abbreviations -> {} (* leave abbreviations *)

To express the result by the stop mixing parameter Xt, call

  GetMSSMCPEvenHiggsLoopMassMatrix[
      parameters -> {At -> Xt - signMu Mu/TanBeta}]
";

GetMSSMCPOddHiggsLoopMass::usage = "Returns the loop corrections to
 the CP-odd Higgs mass in the CP-conserving MSSM.  The
 loop-corrections are taken from arxiv:hep-ph/0105096.

The function takes the same arguments as
GetMSSMCPEvenHiggsLoopMassMatrix[].
"

ReplaceStopMasses::usage = "Returns list of replacetment rules which
 replace mst1, mst2 and sin2Theta by DR-bar Lagrangian parameters.

Usage:

  ReplaceStopMasses[parameters -> {}]

Parameters:

- parameters (optional): List of internal replacement rules for
  parameters, useful when certain limits are considered (default: {}).
  Example: parameters -> {At -> 0, sin2Theta -> 0}
";

GetDeltaMPoleOverMRunningMSSMSQCDDRbar2LUniversalMSUSY::usage =
 "Returns two-loop SUSY-QCD contributions (from squarks and gluino) to
 Delta M_f/m_f in the MSSM in DR-bar scheme for universal SUSY mass
 parameters.  Taken from hep-ph/0210258, Eq. (62).  (mst[i] = i'th
 running stop mass, mg = gluino DR-bar mass, Q = ren. scale).";

GetDeltaMPoleOverMRunningMSSMSQCDDRbar::usage = "Returns the general
 one- and two-loop SUSY-QCD contributions (from squarks and gluino) to
 Delta M_f/m_f in the MSSM in DR-bar scheme.  Taken from the files
 mt.res and mb.res from hep-ph/0210258.  (mst[i] = i'th running stop
 mass, mg = gluino DR-bar mass, Q = ren. scale).

Note: The 2-loop QCD contribution is not included in this expression!

Parameters:

- loopOrder (optional): List of factors multiplied by each loop order.
  #1: 1-loop level
  #2: 2-loop level
  Default: loopOrder -> {1,1}

- parameters (optional): List of internal replacement rules for
  parameters
  Default: parameters -> { sin2Theta -> 2 mq Xt/(mst[1]^2 - mst[2]^2) }
";

(* DR-bar parameters *)
{ ht, Mu, mt, mb, At, mQ33, mU33, TanBeta, g3, Q, M3, signMu, BMu };

(* Stop mass parameters *)
{ sin2Theta, cos2Theta, sin4Theta, cos4Theta, mst1, mst2 };

(* Sbottom mass parameters *)
{ msb1, msb2 };

(* abbreviated parameters *)
{ mA, Xt, Yt, delta, Li2, lam, xp, xm, phi2, phi };

(* options *)
{ loopOrder, corrections, parameters, abbreviations };

(* DR-bar parameters for Delta M_t / m_t *)
{ mg, mst, msb, msq, MSUSY, den, fin };

Begin["TwoLoopMSSM`Private`"];

(* Eqs. (17), (19) of arxiv:hep-ph/0105096 *)
CalculateMStop2[] :=
    Module[{mst2 = {0,0}, mL2, mR2, Xt, s2t},
           mL2 = mQ33^2 + mt^2;
           mR2 = mU33^2 + mt^2;
           Xt = mt Abs[At + signMu Mu / TanBeta]; (* Eq. (14) and below *)
           (* Eq. (17) *)
           mst2[[1]] = 1/2 ((mL2 + mR2) + Sqrt[(mL2 - mR2)^2 + 4 Abs[Xt]^2]);
           mst2[[2]] = 1/2 ((mL2 + mR2) - Sqrt[(mL2 - mR2)^2 + 4 Abs[Xt]^2]);
           (* Eq. (19) *)
           s2t = 2 mt (At + signMu Mu / TanBeta) / (mst2[[1]] - mst2[[2]]);
           {mst2[[1]], mst2[[2]], s2t}
          ];

(* Eqs. (25)-(27) of arxiv:hep-ph/0105096 *)
CreateMassMatrixCPEven[F1_, F2_, F3_, DeltaF2_, DeltaF3_, parameters_List] :=
    Module[{mm = {{0,0},{0,0}}},
           mm[[1,1]] = 1/2 ht^2 Mu^2 sin2Theta^2 F3;
           mm[[1,2]] = (ht^2 signMu Mu mt sin2Theta F2
                        + 1/2 ht^2 At signMu Mu sin2Theta^2 (F3 + DeltaF3));
           mm[[2,2]] = (2 ht^2 mt^2 F1 + 2 ht^2 At mt sin2Theta (F2 + DeltaF2)
                        + 1/2 ht^2 At^2 sin2Theta^2 (F3 + 2 DeltaF3));
           If[PossibleZeroQ[At /. parameters],
              mm[[1,2]] = Expand[mm[[1,2]]];
              mm[[2,2]] = Expand[mm[[2,2]]];
             ];
           If[PossibleZeroQ[sin2Theta /. parameters],
              mm[[1,1]] = 0;
              mm[[1,2]] = Expand[ht^2 signMu Mu mt sin2Theta F2];
              mm[[2,2]] = Expand[2 ht^2 mt^2 F1
                                 + 2 ht^2 At mt sin2Theta (F2 + DeltaF2)];
             ];
           mm[[2,1]] = mm[[1,2]];
           mm
          ];

(* Eq. (31) of arxiv:hep-ph/0105096 *)
GetMSSMCPEvenHiggsLoopMassMatrix1LAlphaTAlphaS[parameters_List] :=
    Module[{F1, F2, F3 = 0, Nc = 3, h = 1/(16 Pi^2)},
           F1 = Nc h Log[mst1^2 mst2^2 / mt^4];
           F2 = Nc h Log[mst1^2 / mst2^2];
           If[!PossibleZeroQ[mst1 - mst2 /. parameters],
              F3 = Nc h (2 - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]);
             ];
           CreateMassMatrixCPEven[F1, F2, F3, 0, 0, {}]
          ];

(* Eq. (A4) of arxiv:hep-ph/0105096 *)
Phi[x_, y_, z_] :=
    Module[{u = x/z, v = y/z, lambda, xp, xm},
           lambda = Sqrt[(1 - u - v)^2 - 4 u v];
           xp = 1/2 (1 + (u - v) - lambda);
           xm = 1/2 (1 - (u - v) - lambda);
           1/lambda (2 Log[xp] Log[xm] - Log[u] Log[v]
                     - 2 (PolyLog[2,xp] + PolyLog[2,xm]) + Pi^2/3)
          ];

(* Eq. (A1) of arxiv:hep-ph/0105096 *)
f1[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{delta},
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           (
               4 (mt^2 + mg^2 - mg mt s2t)/(msqu) (1 - Log[mg^2/Q^2])
               + 4 Log[mt^2/mg^2]
               - 2 Log[msqu/mg^2]
               + 2 / delta (
                   + 4 mg^4 Log[msqu/mg^2]
                   + (mg^4 - msqu^2 +mt^2 (10 mg^2 + 3 mt^2 + 2 (mt^2 mg^2 - mt^4)/msqu)) Log[mt^2/mg^2]
                 )
               + (2 mg s2t)/mt (Log[msqu/Q^2]^2 + 2 Log[mt^2/Q^2] Log[msqu/Q^2])
               + (4 mg s2t)/(mt delta) (
                   + mg^2 (msqu - mt^2 - mg^2) Log[msqu/mg^2]
                   + mt^2 (msqu - 3 mg^2 - 2 mt^2 - (mt^2 mg^2 - mt^4)/msqu) Log[mt^2/mg^2]
                 )
               + (4 mg^2 (mt^2 + mg^2 - msqu - 2 mg mt s2t) / delta
                  - (4 mg s2t)/mt) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eq. (A2) of arxiv:hep-ph/0105096 *)
f2[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{delta, diff},
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           diff = msqu - msqd;
           (
               4 (mt^2 + mg^2)/msqu
               - (4 mg s2t)/(mt diff) (
                   + 3 msqu - (mt^2 msqd)/(msqu)
                 )
               + (2 mg s2t)/(mt diff) (
                   + (4 mt^2 + 5 msqu + msqd) Log[msqu/Q^2]
                   - 2 (mt^2 msqd)/msqu Log[mg^2/Q^2]
                 )
               - 4 (mg^2 + mt^2)/msqu Log[mg^2/Q^2] 
               - 2 Log[msqu/mg^2]
               + 2/delta (
                   2 mg^2 (mg^2 + mt^2 - msqu) Log[msqu/mg^2]
                   + 2 mt^2 (
                       + 3 mg^2 + 2 mt^2 - msqu 
                       + (mg^2 mt^2-mt^4)/(msqu)
                     ) Log[mt^2/mg^2]
                 )
               - (4 mg mt s2t)/(msqu delta) (
                   2 msqu mg^2 Log[msqu/mg^2]
                   - ((mt^2 - msqu)^2 - mg^2 (mt^2 + msqu)) Log[mt^2/mg^2]
                 )
               - (8 mg mt)/(s2t diff) (
                   Log[msqu/Q^2] - Log[mt^2/Q^2] Log[msqu/Q^2]
                 )
               - (mg s2t)/(mt diff) (
                   (msqu + msqd) Log[msqu/Q^2]^2
                   + (10 mt^2 - 2 mg^2 + msqu + msqd) Log[mt^2/Q^2] Log[msqu/Q^2] 
                   + (2 mg^2 - 2 mt^2 + msqu + msqd) Log[msqu/Q^2] Log[mg^2/Q^2]
                 )
               + ((8 mg^2 mt^2)/delta
                  - (8 mg mt)/(s2t diff)
                  + (2 s2t (4 mt^2  mg^2 - delta))/(mg mt diff)
                  + (s2t (msqu - mg^2 - mt^2)^3)/(mg mt delta)
                 ) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eq. (A3) of arxiv:hep-ph/0105096 *)
f3[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{delta, diff},
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           diff = msqu - msqd;
           (
               - 4 (msqd (mg^2 + mt^2))/(msqu diff)
               + (4 mg mt s2t)/(diff^2) (
                   21 msqu - msqd^2/msqu
                 )
               + 4/diff ( 
                   (mg^2msqd)/(msqu)Log[mg^2/Q^2]
                   - 2 (mt^2 + mg^2) Log[msqu/Q^2]
                 )
               - (24 mg mt s2t (3 msqu + msqd))/(diff^2)Log[msqu/Q^2]
               + (4 mt^2)/(msqu delta) (
                   2 mg^2 msqu Log[msqu/Q^2] 
                   -mg^2 (mg^2 - mt^2 +msqu) Log[mg^2/Q^2]  
                   - ((mt^2 - msqu)^2-mg^2 (mt^2+msqu)) Log[mt^2/Q^2] 
                 )
               - (4 mg mt s2t)/(msqu delta)(
                   mt^2 (mg^2 -mt^2+msqu) Log[mt^2/Q^2] 
                   - mg^2 (mg^2 -mt^2-msqu) Log[mg^2/Q^2] 
                   + msqu (mg^2 +mt^2-msqu) Log[msqu/Q^2]
                 )
               + 2 (2 mg^2 + 2 mt^2 - msqu - msqd)/diff Log[mt^2 mg^2/Q^4]Log[msqu/Q^2]
               + (12 mg mt s2t)/(diff^2)(
                   2 (mg^2 - mt^2) Log[mg^2/mt^2] Log[msqu/Q^2]
                   + (msqu + msqd) Log[mt^2 mg^2/Q^4] Log[msqu/Q^2]
                 )
               + (8 mg mt)/(s2t diff^2) (
                   - 8 msqu + 2 (3 msqu + msqd)Log[msqu/Q^2] 
                   - 2 (mg^2 - mt^2) Log[mg^2/mt^2] Log[msqu/Q^2]
                   - (msqu + msqd) Log[mt^2 mg^2/Q^4] Log[msqu/Q^2]
                 )
               - (
                   (8/s2t - 12 s2t) (mt (2 delta + (mg^2 + mt^2 - msqu) diff))/(mg diff^2)
                   + (4 delta + 8 mg^2 mt^2)/(mg^2 diff)
                   + (2 (mg^2 + mt^2 - msqu))/(mg^2)
                   - (4 mt^2 (mg^2 + mt^2 - msqu - 2 mg mt s2t))/delta
                 ) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eqs. (32)-(34) of arxiv:hep-ph/0105096 *)
GetMSSMCPEvenHiggsLoopMassMatrix2LAlphaTAlphaS[parameters_List] :=
    Module[{unit, CF = 4/3, Nc = 3, h = 1/(16 Pi^2),
            mg, cos2Theta2, F1, F2, F3, F31L, DeltaF2, DeltaF3},
           unit = g3^2 CF Nc h^2;
           cos2Theta2 = 1 - sin2Theta^2;
           mg = Abs[M3];
           (* Eq. (32) *)
           F1 = unit (
               - 6 (1 - Log[mt^2 / Q^2])
               + 5 Log[mst1^2 mst2^2 / mt^4]
               + Log[mst1^2 mst2^2 / mt^4]^2
               + 8 Log[mt^2 / Q^2]^2
               - 4 (Log[mst1^2 / Q^2]^2 + Log[mst2^2 / Q^2]^2)
               - cos2Theta2 (2 - Log[mst1^2 mst2^2 / Q^4] - Log[mst1^2 / mst2^2]^2)
               - sin2Theta^2 (
                   + mst1^2 / mst2^2 (1 - Log[mst1^2 / Q^2])
                   + mst2^2 / mst1^2 (1 - Log[mst2^2 / Q^2])
                 )
               + f1[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
               + f1[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
           );
           (* Eq. (33) *)
           If[PossibleZeroQ[mst1 - mst2 /. parameters],
              If[PossibleZeroQ[sin2Theta /. parameters],
                 F2 = 0;
                 ,
                 F2 = unit (
                     + f2[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
                     - f2[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
                 );
                ];
              ,
              F2 = unit (
                  + 5 Log[mst1^2 / mst2^2]
                  - 3 (Log[mst1^2 / Q^2]^2 - Log[mst2^2 / Q^2]^2)
                  + cos2Theta2 (
                      + 5 Log[mst1^2 / mst2^2]
                      - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]^2
                      - 2 / (mst1^2 - mst2^2) (
                          + mst1^2 Log[mst1^2 / Q^2]
                          - mst2^2 Log[mst2^2 / Q^2]
                                              ) Log[mst1^2 / mst2^2]
                               )
                  + sin2Theta^2 (
                      + mst1^2 / mst2^2 (1 - Log[mst1^2 / Q^2])
                      - mst2^2 / mst1^2 (1 - Log[mst2^2 / Q^2])
                                )
                  + f2[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
                  - f2[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
              );
             ];
           (* Eq. (33) *)
           F31L = Nc h (2 - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]);
           F3 = unit (
               + F31L/(h Nc) (3 + 9 cos2Theta2)
               + 4
               - (3 + 13 cos2Theta2) / (mst1^2 - mst2^2) (
                   + mst1^2 Log[mst1^2 / Q^2]
                   - mst2^2 Log[mst2^2 / Q^2]
                 )
               + 3 (mst1^2 + mst2^2)/(mst1^2 - mst2^2) (
                   + Log[mst1^2 / Q^2]^2 - Log[mst2^2 / Q^2]^2
                 )
               - cos2Theta2 (
                   + 4
                   - ((mst1^2 + mst2^2)/(mst1^2 - mst2^2))^2 Log[mst1^2 / mst2^2]^2
                   - 6 (mst1^2 + mst2^2)/(mst1^2 - mst2^2)^2 (
                       + mst1^2 Log[mst1^2 / Q^2]
                       - mst2^2 Log[mst2^2 / Q^2]
                     ) Log[mst1^2 / mst2^2]
                 )
               - sin2Theta^2 (
                   + mst1^2 / mst2^2
                   + mst2^2 / mst1^2
                   + 2 Log[mst1^2 mst2^2 / Q^4]
                   - mst1^2^2 / (mst2^2 (mst1^2 - mst2^2)) Log[mst1^2 / Q^2]
                   + mst2^2^2 / (mst1^2 (mst1^2 - mst2^2)) Log[mst2^2 / Q^2]
                 )
               + f3[mt, mg, mst1^2, mst2^2, sin2Theta, Q]
               + f3[mt, mg, mst2^2, mst1^2, -sin2Theta, Q]
           );
           (* Eq. (35) *)
           DeltaF2 = unit 2 mg / At (
               + Log[mst2^2 / Q^2]^2
               - Log[mst1^2 / Q^2]^2
           );
           (* Eq. (36) *)
           DeltaF3 = unit mg / At (
               + 8
               - 2 (mst1^2 + mst2^2)/(mst1^2 - mst2^2) (
                   + Log[mst2^2 / Q^2]^2
                   - Log[mst1^2 / Q^2]^2
                 )
               + 8 / (mst1^2 - mst2^2) (
                   + mst2^2 Log[mst2^2 / Q^2]
                   - mst1^2 Log[mst1^2 / Q^2]
                 )
           );
           CreateMassMatrixCPEven[F1, F2, F3, DeltaF2, DeltaF3, parameters]
          ];

Get[FileNameJoin[{"meta", "MSSM", "BDSZHiggs.m"}]];

alphaT2Abbreviations = {
    mA -> Sqrt[- signMu BMu (TanBeta + 1/TanBeta)],
    Xt -> At + signMu Mu / TanBeta,
    Yt -> At - signMu Mu TanBeta,
    delta[x_,y_,z_] :> x^2 + y^2 + z^2 - 2 (x y + x z + y z),
    Li2[x_] :> PolyLog[2,x],
    lam[x_,y_,z_] :> Sqrt[(1 - x/z - y/z)^2 - 4 x y/z^2],
    xp[x_,y_,z_] :> 1/2 (1 + x/z -y/z - lam[x,y,z]),
    xm[x_,y_,z_] :> 1/2 (1 - x/z +y/z - lam[x,y,z]),
    phi2[x_,y_,z_] :> (2 Log[xp[x,y,z]] Log[xm[x,y,z]]
                       - Log[x/z]*Log[y/z]
                       - 2 (Li2[xp[x,y,z]] + Li2[xm[x,y,z]])
                       + Pi^2/3) / lam[x,y,z],
    phi[x_,y_,z_] :> Which[x <= z && y <= z,     phi2[x,y,z],
                           z <= x && y <= x, z/x phi2[z,y,x],
                           z <= y && x <= y, z/y phi2[z,x,y]]
};

(* arxiv:hep-ph/0112177 *)
GetMSSMCPEvenHiggsLoopMassMatrix2LAlphaTAlphaT[parameters_List, abbreviations_List] :=
    Module[{Mh},
           Mh = {{M11, M12},
                 {M12, M22}};
           If[PossibleZeroQ[sin2Theta /. parameters],
              Mh = Expand[Mh];
             ];
           Mh //. {
               Nc -> 3,
               t -> mt^2,
               A0 -> mA^2,
               q -> Q^2,
               mu -> signMu Mu,
               mu2 -> Mu^2,
               T1 -> mst1^2,
               T2 -> mst2^2,
               BL -> msb1^2,
               BR -> msb2^2,
               cb -> Cos[ArcTan[TanBeta]],
               sb -> Sin[ArcTan[TanBeta]],
               c2t -> cos2Theta,
               s2t -> sin2Theta
           } //. abbreviations
          ];

(* Eqs. (C1) of arxiv:hep-ph/0105096 *)
CalculateMassCPOdd[FA_, parameters_List] :=
    Module[{result = 0},
           If[!PossibleZeroQ[mst1 - mst2 /. parameters],
              result = Expand[ht^2 signMu Mu At FA / (mst1^2 - mst2^2) (TanBeta + 1/TanBeta)];
             ];
           result /. parameters
          ];

(* Eq. (C3) of arxiv:hep-ph/0105096 *)
GetMSSMCPOddHiggsLoopMass1LAlphaTAlphaS[parameters_List] :=
    Module[{FA, Nc = 3, h = 1/(16 Pi^2)},
           FA = Nc h (mst1^2 (1 - Log[mst1^2/Q^2]) - mst2^2 (1 - Log[mst2^2/Q^2]));
           CalculateMassCPOdd[FA, parameters]
          ];

(* Eq. (C5) of arxiv:hep-ph/0105096 *)
fA[mt_, mg_, msqu_, msqd_, s2t_, Q_] :=
    Module[{diff, delta},
           diff = msqu - msqd;
           delta = mg^4 + mt^4 + msqu^2 - 2 (mg^2 mt^2 + mg^2 msqu + mt^2 msqu);
           (
               + (16 msqu mg mt s2t)/diff
               - (Pi^2 mg (msqu + msqd))/(6 At)
               - 4 (mg^2 + mt^2) Log[msqu/Q^2]
               - (2 mg)/At (
                   + msqu (6 - 5 Log[msqu/Q^2])
                   + msqd (1 - Log[msqd/Q^2])
                 )
               - (4 mg mt s2t (3 msqu + msqd))/diff Log[msqu/Q^2]
               - mg/At (msqu Log[msqu/Q^2]^2 + msqd Log[msqd/Q^2]^2)
               + 2 (mg^2 + mt^2 - msqu) Log[mg^2 mt^2/Q^4] Log[msqu/Q^2]
               + 2 msqu (1 + mg/At) Log[mg^2/Q^2] Log[mt^2/Q^2]
               - 2 mg / At (
                   + (mg^2 - mt^2) Log[mg^2/mt^2]
                   + msqu Log[mg^2 mt^2/Q^4]
                 ) Log[msqu/Q^2]
               + (2 mg mt s2t)/diff (
                   + 2 (mg^2 - mt^2) Log[mg^2/mt^2]
                   + (msqu + msqd) Log[mg^2 mt^2/Q^4]
                 ) Log[msqu/Q^2]
               - (
                   4 mt^2
                   + 2 delta/mg^2 (1 + mg/At)
                   - 2 mt s2t (2 delta + (mg^2 + mt^2 - msqu) diff)/(mg diff)
                 ) Phi[mt^2, msqu, mg^2]
           )
          ];

(* Eq. (C4) of arxiv:hep-ph/0105096 *)
GetMSSMCPOddHiggsLoopMass2LAlphaTAlphaS[parameters_List] :=
    Module[{FA, FA1L, CF = 4/3, Nc = 3, h = 1/(16 Pi^2)},
           unit = g3^2 CF Nc h^2;
           FA1L = (mst1^2 (1 - Log[mst1^2/Q^2]) - mst2^2 (1 - Log[mst2^2/Q^2]));
           FA = unit (
               FA1L (
                   8 - sin2Theta^2 (
                       2
                       - (mst1^2 + mst2^2)/(mst1^2 - mst2^2) Log[mst1^2 / mst2^2]
                   )
                 )
               + 2 (mst1^2 Log[mst1^2/Q^2]^2 - mst2^2 Log[mst2^2/Q^2]^2)
               + 2/(mst1^2 - mst2^2) (mst1^2 Log[mst1^2/Q^2] - mst2^2 Log[mst2^2/Q^2])^2
               + fA[mt, Abs[M3], mst1^2, mst2^2, sin2Theta, Q]
               - fA[mt, Abs[M3], mst2^2, mst1^2, -sin2Theta, Q]
           );
           CalculateMassCPOdd[FA, parameters]
          ];

Options[GetMSSMCPEvenHiggsLoopMassMatrix] = \
Options[GetMSSMCPOddHiggsLoopMass] = {
    loopOrder -> {1,1,1},
    corrections -> {1,1},
    parameters -> {},
    abbreviations -> alphaT2Abbreviations
};

GetMSSMCPEvenHiggsLoopMassMatrix[OptionsPattern[]] :=
    (
        OptionValue[loopOrder][[1]] {{0,0}, {0,0}} +
        OptionValue[loopOrder][[2]] OptionValue[corrections][[1]] GetMSSMCPEvenHiggsLoopMassMatrix1LAlphaTAlphaS[OptionValue[parameters]] +
        OptionValue[loopOrder][[3]] OptionValue[corrections][[1]] GetMSSMCPEvenHiggsLoopMassMatrix2LAlphaTAlphaS[OptionValue[parameters]] +
        OptionValue[loopOrder][[3]] OptionValue[corrections][[2]] GetMSSMCPEvenHiggsLoopMassMatrix2LAlphaTAlphaT[OptionValue[parameters], OptionValue[abbreviations]]
    ) /. OptionValue[parameters];

GetMSSMCPOddHiggsLoopMass[OptionsPattern[]] :=
    (
        OptionValue[loopOrder][[1]] 0 +
        OptionValue[loopOrder][[2]] OptionValue[corrections][[1]] GetMSSMCPOddHiggsLoopMass1LAlphaTAlphaS[OptionValue[parameters]] +
        OptionValue[loopOrder][[3]] OptionValue[corrections][[1]] GetMSSMCPOddHiggsLoopMass2LAlphaTAlphaS[OptionValue[parameters]]
    ) /. OptionValue[parameters];

Options[ReplaceStopMasses] = {
    parameters -> {}
};

ReplaceStopMasses[OptionsPattern[]] :=
    Module[{mstop = CalculateMStop2[] /. OptionValue[parameters]},
           {
               mst1 -> Sqrt[mstop[[1]]],
               mst2 -> Sqrt[mstop[[2]]],
               sin2Theta -> mstop[[3]]
           }
          ];

(* hep-ph/0210258, Eq. (59) *)
GetDeltaMPoleOverMRunningMSSMSQCDDRbar1L[] :=
    With[{CF = 4/3, as = g3^2/(4 Pi), diff1 = mst[1]^2/(mst[1]^2 - mg^2),
          diff2 = mst[2]^2/(mst[2]^2 - mg^2), mq = mt},
         CF as/(8 Pi) (
             -3 + diff1 + diff2 + 2 Log[mg^2/Q^2]
             + diff1 (2 - diff1 - 2 sin2Theta mg/mq) Log[mst[1]^2/mg^2]
             + diff2 (2 - diff2 + 2 sin2Theta mg/mq) Log[mst[2]^2/mg^2]
         )
        ];

(* hep-ph/0210258, Eq. (62) *)
GetDeltaMPoleOverMRunningMSSMSQCDDRbar2LUniversalMSUSY[] :=
    With[{CF = 4/3, CA = 3, as = g3^2/(4 Pi), mq = mt, aq = Xt, M = MSUSY},
         CF (as/(4 Pi))^2 (
             47/3 + 20 Log[M^2/Q^2] + 6 Log[M^2/Q^2] Log[M^2/mq^2]
             + CF (23/24 - 13/6 Log[M^2/Q^2] + 1/2 Log[M^2/Q^2]^2 -
                   3 Log[M^2/Q^2] Log[mq^2/Q^2])
             + CA (175/72 + 41/6 Log[M^2/Q^2] - 1/2 Log[M^2/Q^2]^2 -
                   2 Log[M^2/Q^2] Log[mq^2/Q^2])
             + aq/M (-4 - 8 Log[M^2/Q^2])
             + CF aq/M (7/3 - 11/3 Log[M^2/Q^2] + 3 Log[mq^2/Q^2])
             + CA aq/M (-8/3 + 4 Log[M^2/Q^2])
         )
        ];

resmt = Get[FileNameJoin[{"meta", "MSSM", "tquark_2loop_strong.m"}]];

GetDeltaMPoleOverMRunningMSSMSQCDDRbar2L[] :=
    With[{as = g3^2/(4 Pi)},
         (as/(4 Pi))^2 resmt //. {
             GS -> 1, colorCF -> 4/3, colorCA -> 3, Tf -> 1/2,
             log2 -> Log[2], MGl -> mgl, MT -> mt,
             zt2 -> Zeta[2], zt3 -> Zeta[3],
             mmt -> mt^2, mmb -> mb^2,
             mgl -> mg, mmgl -> mgl^2,
             mmst1 -> mst[1]^2, mmst2 -> mst[2]^2,
             mmsb1 -> msb[1]^2, mmsb2 -> msb[2]^2,
             mmsusy -> msq^2, TwoLoopMSSM`Private`mmu -> Q^2,
             cst -> ct, snt -> st, csb -> cb, snb -> sb,
             cs2t -> c2t, sn2t -> s2t, cs2b -> c2b, sn2b -> s2b,
             cs4t -> c4t, sn4t -> s4t, cs4b -> c4b, sn4b -> s4b,
             c2t -> cos2Theta, s2t -> sin2Theta,
             c4t -> cos4Theta, s4t -> sin4Theta
         }
        ];

Options[GetDeltaMPoleOverMRunningMSSMSQCDDRbar] = {
    loopOrder -> {1, 1},
    parameters -> {
        sin2Theta -> 2 mt Xt/(mst[1]^2 - mst[2]^2),
        sin4Theta -> 2 sin2Theta cos2Theta,
        cos2Theta -> Sqrt[1 - sin2Theta^2],
        cos4Theta -> 1 - 2 sin2Theta^2,
        den[x_, n_] :> x^(-n)
    }
};

GetDeltaMPoleOverMRunningMSSMSQCDDRbar[OptionsPattern[]] := (
    OptionValue[loopOrder][[1]] GetDeltaMPoleOverMRunningMSSMSQCDDRbar1L[] +
    OptionValue[loopOrder][[2]] GetDeltaMPoleOverMRunningMSSMSQCDDRbar2L[]
    ) //. OptionValue[parameters];

End[];
