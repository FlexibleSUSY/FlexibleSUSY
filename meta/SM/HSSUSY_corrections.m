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

(*** High-scale MSSM boundary conditions to the MSSM ***)
(* Taken from arXiv:1407.4081, arxiv:1504.05200, arxiv:1703.08166 *)

(* abbreviations *)
xQU = Sqrt[Abs[msq2[3,3]/msu2[3,3]]];
xQD = Sqrt[Abs[msq2[3,3]/msd2[3,3]]];
xLE = Sqrt[Abs[msl2[3,3]/mse2[3,3]]];

At = AtInput;
xt = At - MuInput/TanBeta;
yt = At + MuInput*TanBeta;
xtt = xt^2/Sqrt[Abs[msq2[3,3] msu2[3,3]]];

Ab = AbInput;
xb = Ab - MuInput*TanBeta;
yb = Ab + MuInput/TanBeta;
xbb = xb^2/Sqrt[Abs[msq2[3,3] msd2[3,3]]];

Atau = AtauInput;
xtau = Atau - MuInput*TanBeta;
xtaut = xtau^2/Sqrt[Abs[msl2[3,3] mse2[3,3]]];
ytau = Atau + MuInput/TanBeta;

(* arXiv:1407.4081, Eq. (5) *)
gYu = Sqrt[3/5] g1 Sin[ArcTan[TanBeta]];
gYd = Sqrt[3/5] g1 Cos[ArcTan[TanBeta]];
g2u = g2 Sin[ArcTan[TanBeta]];
g2d = g2 Cos[ArcTan[TanBeta]];

(* arXiv:1407.4081, Eq. (3) *)
lambdaTree = 1/4 (g2^2 + 3/5 g1^2) Cos[2 ArcTan[TanBeta]]^2;

(* arXiv:1407.4081, Eq. (9) *)
lambda1LReg = 1/(4 Pi)^2 (
    - 9/100 g1^4 - 3/10 g1^2 g2^2
    - (3/4 - Cos[2 ArcTan[TanBeta]]^2/6) * g2^4
    );

(* arXiv:1407.4081, Eq. (10) *)
lambda1LPhi = 1/(4 Pi)^2 (
    3 Yu[3,3]^2 (Yu[3,3]^2 + 1/2 (g2^2-g1^2/5) Cos[2 ArcTan[TanBeta]]) Log[msq2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[3,3]])
    + 3 Yu[3,3]^2 (Yu[3,3]^2 + 2/5 g1^2 Cos[2 ArcTan[TanBeta]]) Log[msu2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[3,3]])
    + Cos[2 ArcTan[TanBeta]]^2/300 (
        3 (g1^4 + 25 g2^4) (
            + Log[msq2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[1,1]])
            + Log[msq2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[2,2]])
            + Log[msq2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msq2[3,3]])
        )
        + 24 g1^4 (
            + Log[msu2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[1,1]])
            + Log[msu2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[2,2]])
            + Log[msu2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msu2[3,3]])
        )
        + 6 g1^4 (
            + Log[msd2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msd2[1,1]])
            + Log[msd2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msd2[2,2]])
            + Log[msd2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msd2[3,3]])
        )
        + (9 g1^4 + 25 g2^4) (
            + Log[msl2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[msl2[1,1]])
            + Log[msl2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[msl2[2,2]])
            + Log[msl2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[msl2[3,3]])
        )
        + 18 g1^4 (
            + Log[mse2[1,1]/SCALE^2] (1 + DeltaEFT v^2/Abs[mse2[1,1]])
            + Log[mse2[2,2]/SCALE^2] (1 + DeltaEFT v^2/Abs[mse2[2,2]])
            + Log[mse2[3,3]/SCALE^2] (1 + DeltaEFT v^2/Abs[mse2[3,3]])
        )
    )
    + (
       + 1/4800 (261 g1^4 + 630 g1^2 g2^2 + 1325 g2^4
                 -4 Cos[4 ArcTan[TanBeta]] (9 g1^4 + 90 g1^2 g2^2 + 175 g2^4)
                 -9 Cos[8 ArcTan[TanBeta]] (3 g1^2 + 5 g2^2)^2) Log[mAInput^2/SCALE^2]
       - 3/16 (3/5 g1^2 + g2^2)^2 Sin[4 ArcTan[TanBeta]]^2
    ) (1 + DeltaEFT v^2/mAInput^2)
    + (
       + 6 Yu[3,3]^4 xtt (TCF[1][xQU] - xtt/12 TCF[2][xQU])
       + 3/4 Yu[3,3]^2 xtt Cos[2 ArcTan[TanBeta]] (3/5 g1^2 TCF[3][xQU] + g2^2 TCF[4][xQU])
       - 1/4 Yu[3,3]^2 xtt Cos[2 ArcTan[TanBeta]]^2 (3/5 g1^2 + g2^2) TCF[5][xQU]
    ) (1 + DeltaEFT v^2/Min[Abs[msq2[3,3]],Abs[msu2[3,3]]])
    );

(* arXiv:1407.4081, Eq. (11) *)
lambda1LChi1 = With[
    { r1 = M1Input / MuInput, r2 = M2Input / MuInput },
    1/(4 Pi)^2 (
        1/2 betaLambda Log[MuInput^2/SCALE^2]                       (1 + DeltaEFT v^2/MuInput^2)
        - 7/12 TCf[1][r1] (gYd^4 + gYu^4)                           (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        - 9/4 TCf[2][r2] (g2d^4 + g2u^4)                            (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
        - 3/2 TCf[3][r1] gYd^2 gYu^2                                (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        - 7/2 TCf[4][r2] g2d^2 g2u^2                                (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
        - 8/3 TCf[5][r1,r2] gYd gYu g2d g2u                         (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        - 7/6 TCf[6][r1,r2] (gYd^2 g2d^2 + gYu^2 g2u^2)             (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        - 1/6 TCf[7][r1,r2] (gYd^2 g2u^2 + gYu^2 g2d^2)             (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        - 4/3 TCf[8][r1,r2] (gYd g2u + gYu g2d) (gYd g2d + gYu g2u) (1 + DeltaEFT v^2/Min[M1Input^2,M2Input^2,MuInput^2])
        + 2/3 TCf0[r1] gYd gYu (lambdaTree - 2 (gYd^2 + gYu^2))     (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        + 2 TCf0[r2] g2d g2u (lambdaTree - 2 (g2d^2 + g2u^2))       (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
        + 1/3 TCg0[r1] lambdaTree (gYd^2 + gYu^2)                   (1 + DeltaEFT v^2/Min[M1Input^2,MuInput^2])
        + TCg0[r2] lambdaTree (g2d^2 + g2u^2)                       (1 + DeltaEFT v^2/Min[M2Input^2,MuInput^2])
    )
];

(* arXiv:1407.4081, Eq. (13) *)
lambda1LChi2 = 1/(4 Pi)^2 (
    -1/6 Cos[2 ArcTan[TanBeta]]^2 (
        + 2 g2^4 Log[M2Input^2/SCALE^2]             (1 + DeltaEFT v^2/M2Input^2)
        + (9/25 g1^4 + g2^4) Log[MuInput^2/SCALE^2] (1 + DeltaEFT v^2/MuInput^2)
    )
);

(* complete lambda 1-loop threshold correction O(gb^4); as discussed
   in arxiv:1703.08166.

   We express the correction in terms of the bottom yukawa of the MSSM
   to resum tanb enhanced corrections
*)
lambda1Lbottom = With[{
    Nc = 3,
    mQ3  = Sqrt[msq2[3,3]],
    mQ32 = msq2[3,3],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    mD32 = msd2[3,3],
    M3 = M3Input,
    Mu = MuInput,
    Xt = xt,
    gb = Yd[3,3], (* SM Yukawa coupling *)
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gtau = Ye[3,3], (* SM Yukawa coupling *)
    Q = SCALE,
    Q2 = SCALE^2,
    MA = mAInput,
    cosb = Cos[ArcTan[TanBeta]],
    cos2beta = Cos[2*ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    Xb = xb,
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    k = 1/(4 Pi)^2,
    CF = 4/3
    },
      ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
        deltagsb   = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
        deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
        deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
        deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
        deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
        (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
       ];
      k*(6 Xtildeb ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4 (TCF[1][xQD] - 1/12 Xtildeb TCF[2][xQD])
      + 3/4 Xtildeb ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 cos2beta (-(3/10) g1^2 TCF[3][xQD] + (-((3 g1^2)/10) - g2^2) TCF[4][xQD])
      - 1/4 ((3 g1^2)/5 + g2^2) Xtildeb ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 cos2beta^2 TCF[5][xQD]
      + 3*ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*(ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 - 1/2*(g2^2+g1^2/5)*cos2beta)*Log[mQ32/Q2]
      + 3 ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 (ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2 - 1/5 g1^2 cos2beta) Log[mD32/Q2]
      ) (1 + DeltaEFT v^2/Min[Abs[msq2[3,3]], Abs[msu2[3,3]], Abs[msd2[3,3]], M3Input^2, MuInput^2])
];

(* lambda 1-loop threshold correction O(alpha_tau),
   taken from SUSYHD 1.0.2, arxiv:1504.05200 *)
lambda1Ltau = With[{
    ML32 = msl2[3,3],
    ME32 = mse2[3,3],
    Q2 = SCALE^2,
    gtau = Ye[3,3], (* SM Yukawa coupling *)
    XtildeTau = xtaut,
    cos2beta = Cos[2*ArcTan[TanBeta]],
    k = 1/(4*Pi)^2
    },
     k*(2*XtildeTau*gtau^4*(TCF[1][xLE] - (XtildeTau*TCF[2][xLE])/12)
     + (XtildeTau*gtau^2*cos2beta*((-9*g1^2*TCF[3][xLE])/10 + ((3*g1^2)/10 - g2^2)*TCF[4][xLE]))/4
     -  (((3*g1^2)/5 + g2^2)*XtildeTau*gtau^2*cos2beta^2*TCF[5][xLE])/12
     + gtau^2*(gtau^2 - (3*g1^2*cos2beta)/5)*Log[ME32/Q2] + gtau^2*(gtau^2 - 1/2*(g2^2 - 3/5*g1^2)*cos2beta)*Log[ML32/Q2]
     ) (1 + DeltaEFT v^2/Min[Abs[msl2[3,3]], Abs[mse2[3,3]], MuInput^2])
    ];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   arXiv:1407.4081, Eq. (36)
   Valid in the limit M3 = MQ[3,3] = MU[3,3] = MSUSY *)
lambda2LPhiHSSAlphaTAlphaSDegenerate = With[{ r = xt / SCALE },
    (g3^2 Yu[3,3]^4)/(96 Pi^4) (
        -12 r - 6 r^2 + 14 r^3 + 1/2 r^4 - r^5
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   arXiv:1407.4081, provided by Pietro Slavich,
   in the limit M3 = MQ[3,3] = MU[3,3] != MSUSY *)
lambda2LPhiHSSAlphaTAlphaSMQMUM3Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x1 = msq2[3,3]/M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (
      -(as*CF*htSM^4*Nc*(Xt*(24 + 12*Xt - 28*Xt^2 - Xt^3 + 2*Xt^4) -
       2*Xt*(24 - 24*Xt - 4*Xt^2 + Xt^3)*Log[g/q2] +
       36*Log[g/q2]^2))/(192*Pi^3)
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   arXiv:1407.4081, provided by Pietro Slavich,
    General expression. *)
lambda2LPhiHSSAlphaTAlphaSFull = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = (M3Input)^2,
    x1 = msq2[3,3]*(1-0.005)/(M3Input*(1+0.001))^2,
    x2 = msu2[3,3]/(M3Input*(1-0.005))^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },    
    htSM^4*as*CF*Nc/(4*Pi)^3 (
       -12*Log[g/q2]^2 + ((-6 + 8*x1 - 4*x1^2)*Log[x1]^2)/(-1 + x1)^2 +
        Log[g/q2]*(8 - 4/x1 - 4/x2 - 12*Log[x1] - 12*Log[x2]) +
        (2*(-7 + 7*x2 - 3*x2^2 + x1*(6 - 5*x2 + 2*x2^2))*Log[x2])/
         ((-1 + x1)*(-1 + x2)^2) + ((-6 + 8*x2 - 4*x2^2)*Log[x2]^2)/(-1 + x2)^2 +
        Log[x1]*((2*(-7 + x1*(7 - 5*x2) + 6*x2 + x1^2*(-3 + 2*x2)))/
           ((-1 + x1)^2*(-1 + x2)) - (2*((-2 + x2)*x2 + x1*(-2 + 8*x2 - 4*x2^2) +
             x1^2*(1 - 4*x2 + 2*x2^2))*Log[x2])/((-1 + x1)^2*(-1 + x2)^2)) +
        Xt^2*(-8/(x1*x2) - (4*(4*x1^3 - x2 - 2*x1^2*(4 + x2) + x1*(3 + 4*x2))*
            Log[x1]^2)/((-1 + x1)^2*(x1 - x2)^2) -
          (4*(7 + 3*x1*(-2 + x2) - 4*x2)*Log[x2])/((-1 + x1)*(x1 - x2)*(-1 + x2)) +
          (4*(x2*(-3 + 8*x2 - 4*x2^2) + x1*(1 - 4*x2 + 2*x2^2))*Log[x2]^2)/
           ((x1 - x2)^2*(-1 + x2)^2) + Log[g/q2]*(8/(x1*x2) -
            (24*Log[x1])/(x1 - x2) + (24*Log[x2])/(x1 - x2)) +
          Log[x1]*((4*(7 - 6*x2 + x1*(-4 + 3*x2)))/((-1 + x1)*(x1 - x2)*(-1 + x2)) +
            (4*(((-2 + x1)*x1*(x1 - x2))/(-1 + x1)^2 + ((-2 + x2)*x2*(-x1 + x2))/
                (-1 + x2)^2 + 2*(x1 + x2))*Log[x2])/(x1 - x2)^2)) +
        Xt^5*((-8*x1*(x1 + x2)*Log[x1]^2)/((-1 + x1)*(x1 - x2)^4) +
          (16*x2*Log[x2])/((-1 + x2)*(-x1 + x2)^3) - (8*x2*(x1 + x2)*Log[x2]^2)/
           ((x1 - x2)^4*(-1 + x2)) + Log[x1]*((16*x1)/((-1 + x1)*(x1 - x2)^3) +
            (8*(x1 + x2)*(-x2 + x1*(-1 + 2*x2))*Log[x2])/((-1 + x1)*(x1 - x2)^4*
              (-1 + x2)))) + Xt*((8*(-2 + x1)*Log[x1]^2)/((-1 + x1)*(x1 - x2)) +
          (16*Log[x2])/(x1 - x2) - (8*(-2 + x2)*Log[x2]^2)/((x1 - x2)*(-1 + x2)) +
          Log[g/q2]*((16*Log[x1])/(x1 - x2) - (16*Log[x2])/(x1 - x2)) +
          Log[x1]*(-16/(x1 - x2) - (8*Log[x2])/((-1 + x1)*(-1 + x2))) +
          (32*(-((-1 + x2)*PolyLog[2, (-1 + x1)/x1]) +
             (-1 + x1)*PolyLog[2, (-1 + x2)/x2]))/((-1 + x1)*(x1 - x2)*(-1 + x2))) -
        (2*(4*x1*(-1 + x2)^2*x2*PolyLog[2, (-1 + x1)/x1] +
           (-1 + x1)*((-1 + x2)*(2*(-1 + x2)*x2 + x1*(-2 + 9*x2 - 6*x2^2) +
               x1^2*(2 - 6*x2 + 3*x2^2)) + 4*(-1 + x1)*x1*x2*
              PolyLog[2, (-1 + x2)/x2])))/((-1 + x1)^2*x1*(-1 + x2)^2*x2) +
        Xt^3*((8*(-2 + x1^2 - 3*x1*(-1 + x2) + x2)*Log[x1]^2)/
           ((-1 + x1)*(x1 - x2)^3) - (16*(x1 + 3*x2)*Log[x2])/(x1 - x2)^3 +
          (8*(-2 + x1 + 3*x2 - 3*x1*x2 + x2^2)*Log[x2]^2)/((-1 + x2)*(-x1 + x2)^3) +
          Log[g/q2]*(32/(x1 - x2)^2 - (16*(x1 + x2)*Log[x1])/(x1 - x2)^3 +
            (16*(x1 + x2)*Log[x2])/(x1 - x2)^3) +
          Log[x1]*((16*(3*x1 + x2))/(x1 - x2)^3 + (16*(x1 + x2 - 2*x1*x2)*Log[x2])/
             ((-1 + x1)*(x1 - x2)^2*(-1 + x2))) -
          (16*(4*(x1 - x2) + (-2 + x1 + x2)*PolyLog[2, (-1 + x1)/x1] -
             (-2 + x1 + x2)*PolyLog[2, (-1 + x2)/x2]))/(x1 - x2)^3) +
        Xt^4*((2*(5*x1^4 + 4*x1^3*(-2 + x2) - 2*x2 - x1^2*x2*(10 + x2) +
             2*x1*(1 + 4*x2 + x2^2))*Log[x1]^2)/((-1 + x1)^2*(x1 - x2)^4) +
          (2*(-4 - 7*x2 + 29*x2^2 - 16*x2^3 + x1^2*(6 - 7*x2 + 3*x2^2) +
             x1*(-3 + 15*x2 - 31*x2^2 + 15*x2^3))*Log[x2])/
           ((-1 + x1)*(x1 - x2)^3*(-1 + x2)^2) -
          (2*(-2*x2 + x1^2*(-2 + x2)*x2 + 8*x2^3 - 5*x2^4 -
             2*x1*(-1 + x2)^2*(-1 + 2*x2))*Log[x2]^2)/((x1 - x2)^4*(-1 + x2)^2) +
          Log[g/q2]*((-4*(x1 + x2 + 6*x1*x2))/(x1*(x1 - x2)^2*x2) +
            (4*(2 + 3*x1 + 3*x2)*Log[x1])/(x1 - x2)^3 -
            (4*(2 + 3*x1 + 3*x2)*Log[x2])/(x1 - x2)^3) +
          Log[x1]*((2*(4 + x1^3*(16 - 15*x2) + 3*x2 - 6*x2^2 +
               x1^2*(-29 + 31*x2 - 3*x2^2) + x1*(7 - 15*x2 + 7*x2^2)))/
             ((-1 + x1)^2*(x1 - x2)^3*(-1 + x2)) +
            (2*(-4*x1^2 - 8*x1*x2 - 4*x2^2 - ((-2 + x1)*x1*(x1 - x2)*(x1 + x2))/
                (-1 + x1)^2 - ((-2 + x2)*x2*(-x1 + x2)*(x1 + x2))/(-1 + x2)^2)*
              Log[x2])/(x1 - x2)^4) -
          (4*(-x1^2 + x1^3 - 6*x1^2*x2 + 5*x1^3*x2 + x2^2 + 6*x1*x2^2 -
             5*x1^3*x2^2 - x2^3 - 5*x1*x2^3 + 5*x1^2*x2^3 +
             (-1 + x1)*x1*(-1 + x2)*x2*(-2 + x1 + x2)*PolyLog[2, (-1 + x1)/x1] -
             (-1 + x1)*x1*(-1 + x2)*x2*(-2 + x1 + x2)*PolyLog[2, (-1 + x2)/x2]))/
           ((-1 + x1)*x1*(x1 - x2)^3*(-1 + x2)*x2))
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   in the limit MQ[3,3]/M3 = MU[3,3]/M3
 *)
lambda2LPhiHSSAlphaTAlphaSX1X2Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (as*CF*htSM^4*Nc*(24*x2^2 - 78*x2^3 + 72*x2^4 - 18*x2^5 -
   48*x2^2*Xt + 96*x2^3*Xt - 48*x2^4*Xt - 24*x2*Xt^2 +
   108*x2^2*Xt^2 - 96*x2^3*Xt^2 + 12*x2^4*Xt^2 + 8*x2*Xt^3 -
   72*x2^2*Xt^3 + 64*x2^3*Xt^3 + 4*Xt^4 - 13*x2*Xt^4 +
   10*x2^2*Xt^4 - x2^3*Xt^4 + 4*x2*Xt^5 - 4*x2^2*Xt^5 -
   24*x2^2*Log[g] + 72*x2^3*Log[g] - 72*x2^4*Log[g] +
   24*x2^5*Log[g] + 48*x2^2*Xt*Log[g] - 96*x2^3*Xt*Log[g] +
   48*x2^4*Xt*Log[g] + 24*x2*Xt^2*Log[g] - 120*x2^2*Xt^2*Log[g] +
   168*x2^3*Xt^2*Log[g] - 72*x2^4*Xt^2*Log[g] - 8*x2*Xt^3*Log[g] +
   16*x2^2*Xt^3*Log[g] - 8*x2^3*Xt^3*Log[g] - 4*Xt^4*Log[g] +
   14*x2*Xt^4*Log[g] - 16*x2^2*Xt^4*Log[g] + 6*x2^3*Xt^4*Log[g] -
   36*x2^3*Log[g]^2 + 72*x2^4*Log[g]^2 - 36*x2^5*Log[g]^2 +
   24*x2^2*Log[q2] - 72*x2^3*Log[q2] + 72*x2^4*Log[q2] -
   24*x2^5*Log[q2] - 48*x2^2*Xt*Log[q2] + 96*x2^3*Xt*Log[q2] -
   48*x2^4*Xt*Log[q2] - 24*x2*Xt^2*Log[q2] +
   120*x2^2*Xt^2*Log[q2] - 168*x2^3*Xt^2*Log[q2] +
   72*x2^4*Xt^2*Log[q2] + 8*x2*Xt^3*Log[q2] -
   16*x2^2*Xt^3*Log[q2] + 8*x2^3*Xt^3*Log[q2] + 4*Xt^4*Log[q2] -
   14*x2*Xt^4*Log[q2] + 16*x2^2*Xt^4*Log[q2] -
   6*x2^3*Xt^4*Log[q2] + 72*x2^3*Log[g]*Log[q2] -
   144*x2^4*Log[g]*Log[q2] + 72*x2^5*Log[g]*Log[q2] -
   36*x2^3*Log[q2]^2 + 72*x2^4*Log[q2]^2 - 36*x2^5*Log[q2]^2 +
   96*x2^2*Xt*Log[1 - (-1 + x2)/x2] +
   32*Xt^3*Log[1 - (-1 + x2)/x2] -
   72*x2*Xt^3*Log[1 - (-1 + x2)/x2] +
   48*x2^2*Xt^3*Log[1 - (-1 + x2)/x2] +
   8*Xt^4*Log[1 - (-1 + x2)/x2] -
   18*x2*Xt^4*Log[1 - (-1 + x2)/x2] +
   12*x2^2*Xt^4*Log[1 - (-1 + x2)/x2] + 84*x2^3*Log[x2] -
   72*x2^4*Log[x2] + 24*x2^5*Log[x2] + 96*x2^2*Xt*Log[x2] -
   144*x2^3*Xt*Log[x2] + 48*x2^4*Xt*Log[x2] -
   48*x2^2*Xt^2*Log[x2] + 168*x2^3*Xt^2*Log[x2] -
   72*x2^4*Xt^2*Log[x2] + 32*Xt^3*Log[x2] - 72*x2*Xt^3*Log[x2] -
   8*x2^3*Xt^3*Log[x2] + 8*Xt^4*Log[x2] - 12*x2*Xt^4*Log[x2] -
   4*x2^2*Xt^4*Log[x2] + 6*x2^3*Xt^4*Log[x2] + 4*x2*Xt^5*Log[x2] -
   72*x2^3*Log[g]*Log[x2] + 144*x2^4*Log[g]*Log[x2] -
   72*x2^5*Log[g]*Log[x2] + 72*x2^3*Log[q2]*Log[x2] -
   144*x2^4*Log[q2]*Log[x2] + 72*x2^5*Log[q2]*Log[x2] -
   36*x2^3*Log[x2]^2 + 72*x2^4*Log[x2]^2 - 36*x2^5*Log[x2]^2 -
   48*x2^3*PolyLog[2, (-1 + x2)/x2] +
   96*x2^3*Xt*PolyLog[2, (-1 + x2)/x2]))/(192*Pi^3*(-1 + x2)^2*x2^3)
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   in the limit MQ[3,3] = M3 (x1 -> 1), MU[3,3] != M3 (x2 != 1)
 *)
lambda2LPhiHSSAlphaTAlphaSMQM3Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x1 = msq2[3,3]/M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (
       -(as*CF*htSM^4*Nc*(4 - 25*x2 + 63*x2^2 - 82*x2^3 + 58*x2^4 -
       21*x2^5 + 3*x2^6 - 32*x2*Xt + 128*x2^2*Xt - 192*x2^3*Xt +
       128*x2^4*Xt - 32*x2^5*Xt - 8*Xt^2 + 32*x2*Xt^2 -
       48*x2^2*Xt^2 + 32*x2^3*Xt^2 - 8*x2^4*Xt^2 - 64*x2*Xt^3 +
       192*x2^2*Xt^3 - 192*x2^3*Xt^3 + 64*x2^4*Xt^3 + 4*Xt^4 +
       14*x2*Xt^4 - 62*x2^2*Xt^4 + 66*x2^3*Xt^4 - 22*x2^4*Xt^4 +
       16*x2*Xt^5 - 32*x2^2*Xt^5 + 16*x2^3*Xt^5 +
       12*(-1 + x2)^5*x2*Log[g/q2]^2 + 11*x2*Log[x2] -
       41*x2^2*Log[x2] + 60*x2^3*Log[x2] - 44*x2^4*Log[x2] +
       17*x2^5*Log[x2] - 3*x2^6*Log[x2] + 24*x2*Xt*Log[x2] -
       96*x2^2*Xt*Log[x2] + 144*x2^3*Xt*Log[x2] -
       96*x2^4*Xt*Log[x2] + 24*x2^5*Xt*Log[x2] -
       22*x2*Xt^2*Log[x2] + 76*x2^2*Xt^2*Log[x2] -
       96*x2^3*Xt^2*Log[x2] + 52*x2^4*Xt^2*Log[x2] -
       10*x2^5*Xt^2*Log[x2] - 32*x2*Xt^3*Log[x2] +
       32*x2^2*Xt^3*Log[x2] + 32*x2^3*Xt^3*Log[x2] -
       32*x2^4*Xt^3*Log[x2] + 19*x2*Xt^4*Log[x2] -
       x2^2*Xt^4*Log[x2] - 47*x2^3*Xt^4*Log[x2] +
       29*x2^4*Xt^4*Log[x2] + 8*x2*Xt^5*Log[x2] +
       16*x2^2*Xt^5*Log[x2] - 24*x2^3*Xt^5*Log[x2] -
       6*x2*Log[x2]^2 + 26*x2^2*Log[x2]^2 - 46*x2^3*Log[x2]^2 +
       42*x2^4*Log[x2]^2 - 20*x2^5*Log[x2]^2 + 4*x2^6*Log[x2]^2 -
       16*x2*Xt*Log[x2]^2 + 56*x2^2*Xt*Log[x2]^2 -
       72*x2^3*Xt*Log[x2]^2 + 40*x2^4*Xt*Log[x2]^2 -
       8*x2^5*Xt*Log[x2]^2 + 4*x2*Xt^2*Log[x2]^2 -
       32*x2^2*Xt^2*Log[x2]^2 + 68*x2^3*Xt^2*Log[x2]^2 -
       56*x2^4*Xt^2*Log[x2]^2 + 16*x2^5*Xt^2*Log[x2]^2 -
       8*x2*Xt^3*Log[x2]^2 + 8*x2^2*Xt^3*Log[x2]^2 +
       8*x2^3*Xt^3*Log[x2]^2 - 8*x2^4*Xt^3*Log[x2]^2 -
       4*x2*Xt^4*Log[x2]^2 + 20*x2^2*Xt^4*Log[x2]^2 -
       2*x2^3*Xt^4*Log[x2]^2 - 10*x2^4*Xt^4*Log[x2]^2 +
       8*x2^2*Xt^5*Log[x2]^2 + 8*x2^3*Xt^5*Log[x2]^2 +
       4*(-1 + x2)^2*Log[g/q2]*(-x2^4 - 2*x2^3*(-2 + Xt^2) -
         (-1 + Xt^2)^2 + x2*(4 - 6*Xt^2 + 8*Xt^3 - 6*Xt^4) +
         x2^2*(-6 + 6*Xt^2 - 8*Xt^3 + 7*Xt^4) +
         x2*(-3 + 3*x2^3 - 4*Xt + 6*Xt^2 + 4*Xt^3 - 5*Xt^4 +
           x2^2*(-9 - 4*Xt + 6*Xt^2) + x2*(9 + 8*Xt - 12*Xt^2 +
             4*Xt^3 - 3*Xt^4))*Log[x2]) + 4*(-1 + x2)^3*x2*
        (2 + 8*Xt + 4*Xt^3 + Xt^4)*PolyLog[2, (-1 + x2)/x2]))/
    (64*Pi^3*(-1 + x2)^5*x2)
    )
];

(* lambda 2-loop threshold correction O(alpha_s alpha_t^2),
   in the limit MU[3,3] = M3 (x2 -> 1), MQ[3,3] != M3 (x1 != 1)
 *)
lambda2LPhiHSSAlphaTAlphaSMUM3Degenerate = With[{
    htSM = Yu[3,3],
    CF = 4/3,
    Nc = 3,
    as = g3^2 / (4 Pi),
    g = M3Input^2,
    x1 = msq2[3,3]/M3Input^2,
    x2 = msu2[3,3]/M3Input^2,
    q2 = SCALE^2,
    Xt = xt / M3Input
    },
    (
       -(as*CF*htSM^4*Nc*(4 - 25*x1 + 63*x1^2 - 82*x1^3 + 58*x1^4 -
       21*x1^5 + 3*x1^6 - 32*x1*Xt + 128*x1^2*Xt - 192*x1^3*Xt +
       128*x1^4*Xt - 32*x1^5*Xt - 8*Xt^2 + 32*x1*Xt^2 -
       48*x1^2*Xt^2 + 32*x1^3*Xt^2 - 8*x1^4*Xt^2 - 64*x1*Xt^3 +
       192*x1^2*Xt^3 - 192*x1^3*Xt^3 + 64*x1^4*Xt^3 + 4*Xt^4 +
       14*x1*Xt^4 - 62*x1^2*Xt^4 + 66*x1^3*Xt^4 - 22*x1^4*Xt^4 +
       16*x1*Xt^5 - 32*x1^2*Xt^5 + 16*x1^3*Xt^5 +
       12*(-1 + x1)^5*x1*Log[g/q2]^2 + 11*x1*Log[x1] -
       41*x1^2*Log[x1] + 60*x1^3*Log[x1] - 44*x1^4*Log[x1] +
       17*x1^5*Log[x1] - 3*x1^6*Log[x1] + 24*x1*Xt*Log[x1] -
       96*x1^2*Xt*Log[x1] + 144*x1^3*Xt*Log[x1] -
       96*x1^4*Xt*Log[x1] + 24*x1^5*Xt*Log[x1] -
       22*x1*Xt^2*Log[x1] + 76*x1^2*Xt^2*Log[x1] -
       96*x1^3*Xt^2*Log[x1] + 52*x1^4*Xt^2*Log[x1] -
       10*x1^5*Xt^2*Log[x1] - 32*x1*Xt^3*Log[x1] +
       32*x1^2*Xt^3*Log[x1] + 32*x1^3*Xt^3*Log[x1] -
       32*x1^4*Xt^3*Log[x1] + 19*x1*Xt^4*Log[x1] -
       x1^2*Xt^4*Log[x1] - 47*x1^3*Xt^4*Log[x1] +
       29*x1^4*Xt^4*Log[x1] + 8*x1*Xt^5*Log[x1] +
       16*x1^2*Xt^5*Log[x1] - 24*x1^3*Xt^5*Log[x1] -
       6*x1*Log[x1]^2 + 26*x1^2*Log[x1]^2 - 46*x1^3*Log[x1]^2 +
       42*x1^4*Log[x1]^2 - 20*x1^5*Log[x1]^2 + 4*x1^6*Log[x1]^2 -
       16*x1*Xt*Log[x1]^2 + 56*x1^2*Xt*Log[x1]^2 -
       72*x1^3*Xt*Log[x1]^2 + 40*x1^4*Xt*Log[x1]^2 -
       8*x1^5*Xt*Log[x1]^2 + 4*x1*Xt^2*Log[x1]^2 -
       32*x1^2*Xt^2*Log[x1]^2 + 68*x1^3*Xt^2*Log[x1]^2 -
       56*x1^4*Xt^2*Log[x1]^2 + 16*x1^5*Xt^2*Log[x1]^2 -
       8*x1*Xt^3*Log[x1]^2 + 8*x1^2*Xt^3*Log[x1]^2 +
       8*x1^3*Xt^3*Log[x1]^2 - 8*x1^4*Xt^3*Log[x1]^2 -
       4*x1*Xt^4*Log[x1]^2 + 20*x1^2*Xt^4*Log[x1]^2 -
       2*x1^3*Xt^4*Log[x1]^2 - 10*x1^4*Xt^4*Log[x1]^2 +
       8*x1^2*Xt^5*Log[x1]^2 + 8*x1^3*Xt^5*Log[x1]^2 +
       4*(-1 + x1)^2*Log[g/q2]*(-x1^4 - 2*x1^3*(-2 + Xt^2) -
         (-1 + Xt^2)^2 + x1*(4 - 6*Xt^2 + 8*Xt^3 - 6*Xt^4) +
         x1^2*(-6 + 6*Xt^2 - 8*Xt^3 + 7*Xt^4) +
         x1*(-3 + 3*x1^3 - 4*Xt + 6*Xt^2 + 4*Xt^3 - 5*Xt^4 +
           x1^2*(-9 - 4*Xt + 6*Xt^2) + x1*(9 + 8*Xt - 12*Xt^2 +
             4*Xt^3 - 3*Xt^4))*Log[x1]) + 4*(-1 + x1)^3*x1*
        (2 + 8*Xt + 4*Xt^3 + Xt^4)*PolyLog[2, (-1 + x1)/x1]))/
    (64*Pi^3*(-1 + x1)^5*x1)
    )
];

(* arXiv:1407.4081, Eq. (12) *)
betaLambda = (
    2 lambdaTree (gYd^2 + gYu^2 + 3 g2d^2 + 3 g2u^2)
    - gYd^4 - gYu^4 - 5 g2d^4 - 5 g2u^4
    - 4 gYd gYu g2d g2u
    - 2 (gYd^2 + g2u^2) (gYu^2 + g2d^2)
);

(* lambda 2-loop threshold correction O(alpha_t^2),
   arxiv:1504.05200, Eq. (21), in the limit of degenerate masses for the stops and the pseudoscalar *)
lambda2LHSSAlphaT2 = With[{
    tan\[Beta] = TanBeta,
    mst = Sqrt[Sqrt[msq2[3,3]*msu2[3,3]]], (* average stop mass *)
    Q = SCALE,
    \[Mu] = MuInput,
    K = -0.1953256,
    Xt = xt,
    gt = Yu[3,3] (* SM Yukawa coupling *)
    },
    (
(3*(1 + tan\[Beta]^2)*(1/2 - 8.34993159891064/(1 + tan\[Beta]^2) +
   (6*\[Mu]^2)/mst^2 +
   (Xt^6*(-1/2 + (1/2 - Log[mst^2/Q^2]/2)/(1 + tan\[Beta]^2) +
      Log[mst^2/Q^2]/2))/mst^6 - 4*Log[mst^2/Q^2] +
   (13*Log[mst^2/Q^2])/(1 + tan\[Beta]^2) - (6*\[Mu]^2*Log[mst^2/Q^2])/
    mst^2 + 3*Log[mst^2/Q^2]^2 - (3*Log[mst^2/Q^2]^2)/(1 + tan\[Beta]^2) +
   (Xt^3*(Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)*
     (0.8747904000000002/(1 + tan\[Beta]^2) - (2*Log[mst^2/Q^2])/
       (1 + tan\[Beta]^2)))/mst^3 +
   (Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)^2*
    (0.1252095999999998/(1 + tan\[Beta]^2) + (3*Log[mst^2/Q^2])/
      (1 + tan\[Beta]^2)) +
   (Xt*(Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)*
     (0.5008383999999992/(1 + tan\[Beta]^2) + (12*Log[mst^2/Q^2])/
       (1 + tan\[Beta]^2)))/mst +
   4*Piecewise[{{-9/4, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
     Re[((-1 + (2*\[Mu]^2)/mst^2 + (2*\[Mu]^4)/mst^4)*
        (-Pi^2/6 - (\[Mu]^2*Log[\[Mu]^2/mst^2])/mst^2 + Log[\[Mu]^2/mst^2]*
          Log[Abs[1 - \[Mu]^2/mst^2]] + PolyLog[2, \[Mu]^2/mst^2]))/
       Abs[(1 - \[Mu]^2/mst^2)^2]]] -
   8*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
     (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))] -
   (2*\[Mu]^2*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
      (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))])/mst^2 +
   (3*\[Mu]^2*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
      (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
       (1 - \[Mu]^2/mst^2)])/mst^2 +
   (Xt^2*(-7 + 19.6878144/(1 + tan\[Beta]^2) - (6*\[Mu]^2)/mst^2 +
      27*Log[mst^2/Q^2] - (24*Log[mst^2/Q^2])/(1 + tan\[Beta]^2) +
      (6*\[Mu]^2*Log[mst^2/Q^2])/mst^2 +
      (Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)^2*
       (-0.021147733333332752/(1 + tan\[Beta]^2) - (3*Log[mst^2/Q^2])/
         (1 + tan\[Beta]^2)) +
      4*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))] -
      (6*\[Mu]^2*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))])/mst^2 -
      4*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
         (1 - \[Mu]^2/mst^2)] -
      (6*\[Mu]^2*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
          (1 - \[Mu]^2/mst^2)])/mst^2))/mst^2 +
   (Xt^4*(11/2 - 25/(4*(1 + tan\[Beta]^2)) + \[Mu]^2/mst^2 -
      (13*Log[mst^2/Q^2])/2 + (6*Log[mst^2/Q^2])/(1 + tan\[Beta]^2) -
      (\[Mu]^2*Log[mst^2/Q^2])/mst^2 +
      (Xt/mst + (2*\[Mu]*Csc[2*ArcTan[tan\[Beta]]])/mst)^2*
       (-0.020728533333333354/(1 + tan\[Beta]^2) + Log[mst^2/Q^2]/
         (2*(1 + tan\[Beta]^2))) -
      Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))]/2 +
      (\[Mu]^2*Piecewise[{{-1, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2))])/mst^2 +
      Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
        (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
         (1 - \[Mu]^2/mst^2)]/2 +
      (\[Mu]^2*Piecewise[{{1/2, Abs[-1 + \[Mu]^2/mst^2] < 1/100000}},
         (1 + (\[Mu]^2*Log[\[Mu]^2/mst^2])/(mst^2*(1 - \[Mu]^2/mst^2)))/
          (1 - \[Mu]^2/mst^2)])/(2*mst^2)))/mst^4)*gt^6)/
    (256*Pi^4*tan\[Beta]^2)
    ) /. Piecewise[{{val_, cond_}}, default_] :> If[cond, val, default]
    ];


(* lambda 2-loop threshold correction O(alpha_t^2), computation for generic masses *)
lambda2LHSSAlphaT2Generic = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.02),
    U  = msu2[3,3]*(1-0.02),
    mu2 = MuInput^2*(1+0.01),
    q2 = SCALE^2, (* renormalization/matching scale *)
    gt = Yu[3,3], (* SM Yukawa coupling *)
    A0 = mAInput^2*(1-0.01),
    Xt = xt,
    Yt = yt
    },
      gt^6*k^2*(Nc*(5 + 3*Log[Q]^2 + Log[Q]*(-6 - 6*Log[q2]) + 5*Log[q2]^2 +
    Log[q2]*(10 - 4*Log[U]) - 4*Log[U] + 2*Log[U]^2) +
  Nc*Xt^4*((2*(2*Q^2 - 27*Q*U + U^2))/(Q*(Q - U)^2*U) +
    ((-7*Q - 9*U)*Log[Q]^2)/(Q - U)^3 + ((-6*Q^2 - 44*Q*U + 2*U^2)*Log[U])/
     (Q*(Q - U)^3) + ((3*Q + 5*U)*Log[U]^2)/(Q - U)^3 +
    Log[Q]*((-4*(Q^2 - 10*Q*U - 3*U^2))/((Q - U)^3*U) +
      (2*(5*Q + 7*U)*Log[q2])/(Q - U)^3 + (4*(Q + U)*Log[U])/(Q - U)^3) +
    Log[q2]*((2*(2*Q^2 - 15*Q*U + U^2))/(Q*(Q - U)^2*U) -
      (2*(5*Q + 7*U)*Log[U])/(Q - U)^3) - (6*PolyLog[2, 1 - Q/U])/
     (Q - U)^2) + Nc*Xt^2*((-2*Q^2 - 3*Q*U + U^2)/(Q*(Q - U)*U) +
    ((7*Q - 9*U)*Log[Q]^2)/(Q - U)^2 + ((16*Q^2 - 21*Q*U + U^2)*Log[U])/
     (Q*(Q - U)^2) - (9*Log[U]^2)/(Q - U) +
    Log[Q]*((2*Q^2 - 17*Q*U + 19*U^2)/((Q - U)^2*U) -
      (2*(8*Q - 9*U)*Log[q2])/(Q - U)^2 + (2*Q*Log[U])/(Q - U)^2) +
    Log[q2]*((-2*Q^2 - Q*U + U^2)/(Q*(Q - U)*U) + (2*(8*Q - 9*U)*Log[U])/
       (Q - U)^2) - (6*PolyLog[2, 1 - Q/U])/(Q - U)) +
  Nc*Xt^6*((-2*Q^2 - 11*Q*U + U^2)/(Q*(Q - U)^3*U) +
    ((-9*Q - 5*U)*Log[Q]^2)/(Q - U)^4 + ((-10*Q^2 - 3*Q*U + U^2)*Log[U])/
     (Q*(Q - U)^4) + ((-3*Q - 5*U)*Log[U]^2)/(Q - U)^4 +
    Log[q2]*((-2*Q^2 - 5*Q*U + U^2)/(Q*(Q - U)^3*U) -
      (6*Q*Log[U])/(Q - U)^4) +
    Log[Q]*((2*Q^2 + 13*Q*U - 3*U^2)/((Q - U)^4*U) +
      (6*Q*Log[q2])/(Q - U)^4 + (2*(6*Q + 5*U)*Log[U])/(Q - U)^4) -
    (2*PolyLog[2, 1 - Q/U])/(Q - U)^3 + (4*PolyLog[2, (Q - U)/Q])/
     (Q - U)^3) +
  (cbe^2*(Nc*Yt^2*(-Q^(-1) - 2/U +
        ((((Q - U)^4 - A0^3*(Q + U) - A0*(Q - U)^2*(3*Q + 5*U) +
             A0^2*(3*Q^2 + 4*Q*U + 5*U^2))*delta[A0, U, Q] +
           delta[A0, Q, U]*(-(Q*(-A0^2 + (Q - U)^2)*U) +
             (-Q^2 + 3*Q*U - U^2 + A0*(Q + U))*delta[A0, U, Q]))*Log[A0])/
         (2*Q*U^2*delta[A0, Q, U]*delta[A0, U, Q]) +
        (((A0^3 - (Q - U)^3 - 3*A0^2*(Q + U) + 3*A0*(Q^2 - U^2))*
            delta[A0, U, Q] + delta[A0, Q, U]*(-(U*(A0^2 - Q^2 - 2*A0*U +
                U^2)) + (-A0 + Q + 2*U)*delta[A0, U, Q]))*Log[Q])/
         (2*U^2*delta[A0, Q, U]*delta[A0, U, Q]) + (-Q^(-1) - 2/U)*Log[q2] +
        (((A0^3 + (Q - U)^3 - A0^2*(Q + 5*U) - A0*(Q^2 + 4*Q*U - 5*U^2))*
            delta[A0, U, Q] - delta[A0, Q, U]*(2*Q*(A0 + Q - U)*U +
             (A0 + Q - 3*U)*delta[A0, U, Q]))*Log[U])/(2*Q*U*delta[A0, Q, U]*
          delta[A0, U, Q]) + ((A0^4 + (Q - U)^4 - 4*A0^3*(Q + U) -
           4*A0*(Q - U)^2*(Q + U) + 2*A0^2*(3*Q^2 + 2*Q*U + U^2) -
           3*(A0^2 + (Q - U)^2 - 2*A0*(Q + U))*delta[A0, Q, U] +
           2*delta[A0, Q, U]^2)*phi[A0, Q, U])/(2*U^3*delta[A0, Q, U]) +
        ((-(A0 + Q - U)^2 + delta[A0, U, Q])*phi[A0, U, Q])/
         (2*Q*delta[A0, U, Q])) + Nc*(3/2 + Pi^2 - A0/Q - (2*A0)/U +
        3*Log[A0]^2 + 3*Log[Q]^2 + Log[Q]*(-9/2 - 3*Log[q2]) + 2*Log[q2]^2 +
        Log[A0]*(7 + A0*(Q^(-1) + 2/U) - 3*Log[Q] - 3*Log[U]) +
        Log[q2]*(-((A0*(2*Q + U))/(Q*U)) - Log[U]) - (5*Log[U])/2 +
        2*Log[U]^2 + ((A0^2 - 7*A0*Q - delta[A0, Q, Q])*phi[A0, Q, Q])/Q^2 +
        ((A0^2 - 6*A0*U - delta[A0, U, U])*phi[A0, U, U])/U^2) +
      Nc*Xt*Yt*((6*Log[Q]^2)/(Q - U) + (12*Log[U])/(Q - U) +
        (12*Log[q2]*Log[U])/(Q - U) - (4*Log[U]^2)/(Q - U) +
        Log[A0]*((2*Log[Q])/(Q - U) - (2*Log[U])/(Q - U)) +
        Log[Q]*(-12/(Q - U) - (12*Log[q2])/(Q - U) - (2*Log[U])/(Q - U)) -
        (2*(-(A0*(A0 - 7*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
         (Q^2*(Q - U)) + (2*(A0 + Q - U)*phi[A0, U, Q])/(Q*(Q - U)) +
        (2*(-(A0*(A0 - 6*U)) + delta[A0, U, U])*phi[A0, U, U])/
         ((Q - U)*U^2)) + Nc*Xt^3*Yt*(-48/(Q - U)^2 +
        (2*(4*A0 - 5*Q - 3*U)*Log[Q]^2)/(Q - U)^3 - (12*(Q + 3*U)*Log[U])/
         (Q - U)^3 + (4*(-A0 + Q + U)*Log[U]^2)/(Q - U)^3 +
        Log[A0]*((2*(-6*A0 + Q - U)*Log[Q])/(Q - U)^3 +
          (2*(6*A0 - Q + U)*Log[U])/(Q - U)^3) +
        Log[q2]*(-24/(Q - U)^2 - (12*(Q + U)*Log[U])/(Q - U)^3) +
        Log[Q]*((12*(3*Q + U))/(Q - U)^3 + (12*(Q + U)*Log[q2])/(Q - U)^3 +
          (2*(-2*A0 + 3*Q + U)*Log[U])/(Q - U)^3) +
        (2*(A0*(A0 - 7*Q)*(Q - U) + (-5*Q + U)*delta[A0, Q, Q])*
          phi[A0, Q, Q])/(Q^2*(Q - U)^3) -
        (2*((Q - U)*(A0 + Q - U) - 2*delta[A0, U, Q])*phi[A0, U, Q])/
         (Q*(Q - U)^3) + (2*(A0*(A0 - 6*U)*(Q - U) - (Q - 3*U)*
            delta[A0, U, U])*phi[A0, U, U])/((Q - U)^3*U^2)) +
      Xt^2*(Nc*Yt^2*((4*Q - 2*U)/(Q^2*U - Q*U^2) + (3*Log[Q]^2)/(Q - U)^2 +
          ((-((Q - U)*((Q - U)^2*(2*Q - U) + A0^2*(2*Q + U) - 4*A0*Q*
                 (Q + 2*U))*delta[A0, U, Q]) + delta[A0, Q, U]*
              (2*Q*(Q - U)*(A0 + Q - U)*U + (2*Q^2 + Q*U - U^2)*delta[A0, U,
                 Q]))*Log[U])/(Q*(Q - U)^2*U*delta[A0, Q, U]*
            delta[A0, U, Q]) + (2*Log[U]^2)/(Q - U)^2 +
          Log[Q]*(((Q - U)*(-A0^3 + Q*(Q - U)^2 + A0^2*(3*Q + 4*U) -
                A0*(3*Q^2 + 2*Q*U + 7*U^2))*delta[A0, U, Q] +
              delta[A0, Q, U]*(-((Q - U)*U*(-A0^2 + Q^2 + 2*A0*U - U^2)) +
                (-Q^2 + A0*(Q - U) - 2*Q*U + U^2)*delta[A0, U, Q]))/
             ((Q - U)^2*U^2*delta[A0, Q, U]*delta[A0, U, Q]) -
            (2*Log[q2])/(Q - U)^2 - (5*Log[U])/(Q - U)^2) +
          Log[A0]*(((A0^3*Q - (Q - U)^4 + A0*Q*(3*Q^2 - 2*Q*U - U^2) +
                A0^2*(-3*Q^2 - 2*Q*U + U^2))*delta[A0, U, Q] +
              delta[A0, Q, U]*(Q*(-A0^2 + (Q - U)^2)*U + (-(A0*Q) + Q^2 -
                  3*Q*U + U^2)*delta[A0, U, Q]))/(Q*(Q - U)*U^2*
              delta[A0, Q, U]*delta[A0, U, Q]) + Log[Q]/(Q - U)^2 -
            Log[U]/(Q - U)^2) + Log[q2]*((4*Q - 2*U)/(Q^2*U - Q*U^2) +
            (2*Log[U])/(Q - U)^2) + ((A0*(A0 - 7*Q) - delta[A0, Q, Q])*
            phi[A0, Q, Q])/(Q^2*(Q - U)^2) +
          ((-((Q - U)*(A0^4 - 4*A0*Q^2*(Q - U) + (Q - U)^4 - 4*A0^3*(Q + U) +
                A0^2*(6*Q^2 + 4*Q*U + 6*U^2))) + (A0^2*(3*Q - 5*U) +
               (3*Q - 5*U)*(Q - U)^2 + A0*(-6*Q^2 + 4*Q*U + 14*U^2))*
              delta[A0, Q, U] - 2*(Q - 2*U)*delta[A0, Q, U]^2)*phi[A0, Q, U])/
           ((Q - U)^2*U^3*delta[A0, Q, U]) +
          (((Q - U)*(A0 + Q - U)^2 + A0*delta[A0, U, Q])*phi[A0, U, Q])/
           (Q*(Q - U)^2*delta[A0, U, Q]) + ((A0*(A0 - 6*U) - delta[A0, U, U])*
            phi[A0, U, U])/((Q - U)^2*U^2)) +
        Nc*((4*A0*Q - 2*A0*U - 4*Q*U)/(Q^2*U - Q*U^2) +
          ((A0 + Q - 3*U)*Log[Q]^2)/(Q - U)^2 + ((2*A0 + Q - 5*U)*Log[U])/
           (Q - U)^2 - (2*Log[U]^2)/(Q - U) + Log[q2]*
           ((4*A0*Q - 2*A0*U - 2*Q*U)/(Q^2*U - Q*U^2) + (2*(A0 - Q)*Log[U])/
             (Q - U)^2) + Log[Q]*((-2*A0 + Q + 3*U)/(Q - U)^2 -
            (2*(A0 - Q)*Log[q2])/(Q - U)^2 + ((-A0 + Q + U)*Log[U])/
             (Q - U)^2) + Log[A0]*((2*A0*(-2*Q + U))/(Q*(Q - U)*U) +
            ((A0 - 5*Q + 5*U)*Log[Q])/(Q - U)^2 - ((A0 - 5*Q + 5*U)*Log[U])/
             (Q - U)^2) + ((A0*(A0 - 7*Q)*(Q - U) + (-2*Q + U)*
              delta[A0, Q, Q])*phi[A0, Q, Q])/(Q^2*(Q - U)^2) +
          (delta[A0, U, Q]*phi[A0, U, Q])/(Q*(Q - U)^2) +
          ((-(A0*(A0 - 6*U)) + delta[A0, U, U])*phi[A0, U, U])/
           ((Q - U)*U^2))) +
      Xt^4*(Nc*((3*Q*(Q - U)*U + A0*(-2*Q^2 - 5*Q*U + U^2))/(Q*(Q - U)^3*U) +
          Log[Q]*((3*(4*A0*Q - Q^2 + U^2))/(2*(Q - U)^4) +
            (3*(2*A0*Q - Q^2 + U^2)*Log[q2])/(Q - U)^4) +
          (3*(-4*A0*Q + Q^2 - U^2)*Log[U])/(2*(Q - U)^4) +
          Log[q2]*((6*Q*(Q - U)*U + A0*(-2*Q^2 - 5*Q*U + U^2))/
             (Q*(Q - U)^3*U) + (3*(-2*A0*Q + Q^2 - U^2)*Log[U])/(Q - U)^4) +
          Log[A0]*((6*Q*U*(-Q + U) + A0*(2*Q^2 + 5*Q*U - U^2))/
             (Q*(Q - U)^3*U) + (3*(-2*A0*Q + Q^2 - U^2)*Log[Q])/(Q - U)^4 +
            (3*(2*A0*Q - Q^2 + U^2)*Log[U])/(Q - U)^4)) +
        Nc*Yt^2*((-2*Q^2 - 11*Q*U + U^2)/(Q*(Q - U)^3*U) +
          ((7*A0 - 11*Q - 3*U)*Log[Q]^2)/(Q - U)^4 -
          (((Q - U)^2*(A0^3 - 3*Q^3 + 11*Q^2*U - 5*Q*U^2 - 3*U^3 - A0^2*
                (5*Q + 3*U) + A0*(7*Q^2 + 4*Q*U + 5*U^2))*delta[A0, U, Q] +
             delta[A0, Q, U]*(2*Q*(Q - U)^2*(A0 + Q - U)*U + (3*Q^3 -
                 A0*(Q - U)^2 + 7*Q^2*U + 13*Q*U^2 + U^3)*delta[A0, U, Q]))*
            Log[U])/(2*Q*(Q - U)^4*U*delta[A0, Q, U]*delta[A0, U, Q]) -
          (2*(-2*A0 + Q + 3*U)*Log[U]^2)/(Q - U)^4 +
          Log[q2]*((-2*Q^2 - 5*Q*U + U^2)/(Q*(Q - U)^3*U) -
            (6*Q*Log[U])/(Q - U)^4) + Log[A0]*
           (((Q - U)*(-A0^3 + Q^3 - 3*Q^2*U - Q*U^2 + 3*U^3 + 3*A0^2*
                 (Q + U) - A0*(3*Q^2 + 5*U^2))*delta[A0, U, Q] +
              delta[A0, Q, U]*(-(Q*(-A0^2 + (Q - U)^2)*U) +
                (-Q^2 + A0*(Q - U) + 3*Q*U + 3*U^2)*delta[A0, U, Q]))/
             (2*Q*(Q - U)^2*U^2*delta[A0, Q, U]*delta[A0, U, Q]) -
            (3*(A0 - Q + U)*Log[Q])/(Q - U)^4 + (3*(A0 - Q + U)*Log[U])/
             (Q - U)^4) + Log[Q]*(((Q - U)^2*(A0^3 - Q^3 + Q^2*U + 9*Q*U^2 -
                9*U^3 - A0^2*(3*Q + 5*U) + A0*(3*Q^2 + 4*Q*U + 9*U^2))*delta[
                A0, U, Q] + delta[A0, Q, U]*((Q - U)^2*U*(-A0^2 + Q^2 +
                  2*A0*U - U^2) + (Q^3 - A0*(Q - U)^2 + 2*Q^2*U + 17*Q*U^2 +
                  4*U^3)*delta[A0, U, Q]))/(2*(Q - U)^4*U^2*delta[A0, Q, U]*
              delta[A0, U, Q]) + (6*Q*Log[q2])/(Q - U)^4 +
            ((-11*A0 + 13*Q + 9*U)*Log[U])/(Q - U)^4) +
          ((A0*(A0 - 7*Q)*(Q - U) + (-8*Q + U)*delta[A0, Q, Q])*
            phi[A0, Q, Q])/(Q^2*(Q - U)^4) +
          (((Q - U)^2*(A0^4 - 4*A0^3*(Q + U) - 4*A0*(Q - U)^2*(Q + U) +
               (Q - U)^2*(Q^2 - 2*Q*U - 3*U^2) + A0^2*(6*Q^2 + 4*Q*U +
                 6*U^2)) - (Q - U)*(3*Q^3 + 3*A0^2*(Q - 3*U) - 15*Q^2*U + 29*
                Q*U^2 - 17*U^3 - 6*A0*(Q^2 - 2*Q*U - 3*U^2))*
              delta[A0, Q, U] + 2*(Q^2 - 5*Q*U + 12*U^2)*delta[A0, Q, U]^2)*
            phi[A0, Q, U])/(2*(Q - U)^4*U^3*delta[A0, Q, U]) -
          (((Q - U)^2*(A0 + Q - U)^2 + (4*A0 + 3*Q - 3*U)*(Q - U)*
              delta[A0, U, Q] - 6*delta[A0, U, Q]^2)*phi[A0, U, Q])/
           (2*Q*(Q - U)^4*delta[A0, U, Q]) +
          ((A0*(A0 - 6*U)*(-Q + U) + (Q - 5*U)*delta[A0, U, U])*
            phi[A0, U, U])/((Q - U)^4*U^2)))) +
    Nc*((4*mu2^3*(2*Q + U) - Q*U*(4*Q^2 + 3*Q*U + 2*U^2) -
        3*mu2^2*(4*Q^2 + 7*Q*U + 2*U^2) + mu2*(4*Q^3 + 17*Q^2*U + 13*Q*U^2 +
          2*U^3))/(2*(mu2 - Q)*Q*(mu2 - U)*U) +
      ((4*mu2^4 - 4*mu2^3*U - 2*Q^2*U^2 + 2*mu2^2*U*(-4*Q + U) +
         4*mu2*Q*U*(Q + U))*Log[mu2]^2)/((mu2 - Q)^2*(mu2 - U)^2) +
      (2*(2*mu2 - Q)*Q*Log[Q]^2)/(mu2 - Q)^2 - 2*Log[q2]^2 +
      ((-(Q*U^2*(Q + 2*U)) + mu2^3*(13*Q + 2*U) +
         mu2*U*(4*Q^2 + 9*Q*U + 2*U^2) - mu2^2*(9*Q^2 + 14*Q*U + 4*U^2))*
        Log[U])/(2*(mu2 - Q)*Q*(mu2 - U)^2) + ((2*mu2 - U)*U*Log[U]^2)/
       (mu2 - U)^2 + Log[q2]*(-(((-2*mu2 + Q + U)*(2*Q + U))/(Q*U)) +
        Log[U]) + Log[mu2]*((mu2*(-2*Q^2*U^2*(Q + U) - 2*mu2^4*(2*Q + U) +
           mu2*Q*U*(3*Q^2 + 8*Q*U - 2*U^2) + mu2^3*(8*Q^2 + Q*U + 4*U^2) -
           2*mu2^2*(2*Q^3 + 4*Q^2*U - Q*U^2 + U^3)))/((mu2 - Q)^2*Q*
          (mu2 - U)^2*U) + ((-5*mu2^4 - 2*mu2^3*(Q - 4*U) + 2*Q^2*U^2 -
           4*mu2*Q*U*(Q + U) + mu2^2*(Q^2 + 8*Q*U - 4*U^2))*Log[Q])/
         ((mu2 - Q)^2*(mu2 - U)^2) + ((-3*mu2^4 + 2*mu2^3*Q -
           mu2^2*Q*(Q - 8*U) + 2*Q^2*U^2 - 4*mu2*Q*U*(Q + U))*Log[U])/
         ((mu2 - Q)^2*(mu2 - U)^2)) +
      Log[Q]*((-(Q^2*U*(4*Q + 5*U)) + mu2^3*(4*Q + 15*U) +
          mu2*Q*(4*Q^2 + 15*Q*U + 6*U^2) - mu2^2*(8*Q^2 + 14*Q*U + 13*U^2))/
         (2*(mu2 - Q)^2*(mu2 - U)*U) + 3*Log[q2] +
        ((2*mu2^4 - 2*mu2^3*U - Q^2*U^2 + mu2^2*U*(-4*Q + U) +
           2*mu2*Q*U*(Q + U))*Log[U])/((mu2 - Q)^2*(mu2 - U)^2)) +
      (4*mu2^2*PolyLog[2, (mu2 - Q)/mu2])/(mu2 - Q)^2 +
      (2*(mu2^2 + 2*mu2*Q - Q^2)*PolyLog[2, 1 - Q/mu2])/(mu2 - Q)^2 +
      (2*(mu2^2 + 2*mu2*U - U^2)*PolyLog[2, 1 - U/mu2])/(mu2 - U)^2) +
    Xt^2*(Nc^2*((4*Q*Log[Q]^2)/(Q - U)^2 + (4*Log[U])/(Q - U) +
        (4*Log[q2]*Log[U])/(Q - U) + (4*U*Log[U]^2)/(Q - U)^2 +
        Log[Q]*(-4/(Q - U) - (4*Log[q2])/(Q - U) - (4*(Q + U)*Log[U])/
           (Q - U)^2)) + Nc*((-8*mu2*Q + 4*Q^2 + 4*mu2*U + 6*Q*U - 2*U^2)/
         (Q*(Q - U)*U) + (2*Q*(3*mu2^2 - 2*mu2*(Q + 2*U) + Q*(Q + 2*U))*
          Log[Q]^2)/((mu2 - Q)^2*(Q - U)^2) +
        ((-4*mu2^3*Q + Q*U*(7*Q^2 + 3*Q*U - 2*U^2) +
           mu2^2*(Q^2 + 17*Q*U - 2*U^2) - mu2*(Q^3 + 16*Q^2*U + 5*Q*U^2 -
             2*U^3))*Log[U])/((mu2 - Q)*Q*(mu2 - U)*(Q - U)^2) +
        (2*U*(2*mu2^2 - 2*mu2*(Q + U) + U*(Q + U))*Log[U]^2)/
         ((mu2 - U)^2*(Q - U)^2) + Log[q2]*
         ((-8*mu2*Q + 4*Q^2 + 4*mu2*U + 4*Q*U - 2*U^2)/(Q*(Q - U)*U) +
          (2*(-2*mu2 + 2*Q + U)*Log[U])/(Q - U)^2) +
        Log[Q]*((4*mu2^3*U + Q*U*(-4*Q^2 - 7*Q*U + 3*U^2) -
            mu2^2*(4*Q^2 + Q*U + 11*U^2) + mu2*(4*Q^3 + 9*Q^2*U + 2*Q*U^2 +
              5*U^3))/((mu2 - Q)*(mu2 - U)*(Q - U)^2*U) -
          (2*(-2*mu2 + 2*Q + U)*Log[q2])/(Q - U)^2 +
          ((-6*Q)/(Q - U)^2 + (4*Q*(-2*mu2 + Q))/((mu2 - Q)^2*(Q - U)) -
            (4*U)/(Q - U)^2 + (2*U*(-2*mu2 + U))/((mu2 - U)^2*(-Q + U)))*
           Log[U]) + Log[mu2]*((4*mu2*(mu2^2*(2*Q - U) + Q^2*U +
             mu2*(-2*Q^2 - Q*U + U^2)))/((mu2 - Q)*Q*(mu2 - U)*(Q - U)*U) +
          (2*(-2*mu2^5 + 2*Q^3*U^2 - 2*mu2*Q^2*U*(2*Q + 3*U) +
             mu2^4*(3*Q + 7*U) - 2*mu2^3*(2*Q^2 + 5*Q*U + 3*U^2) +
             mu2^2*(Q^3 + 13*Q^2*U + 4*Q*U^2 + 2*U^3))*Log[Q])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^2) +
          (2*(2*mu2^5 - 2*Q^3*U^2 + 2*mu2*Q^2*U*(2*Q + 3*U) -
             mu2^4*(3*Q + 7*U) + 2*mu2^3*(2*Q^2 + 5*Q*U + 3*U^2) -
             mu2^2*(Q^3 + 13*Q^2*U + 4*Q*U^2 + 2*U^3))*Log[U])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^2)) +
        (4*(mu2 - Q)*PolyLog[2, (mu2 - Q)/mu2])/(Q - U)^2 -
        (4*(mu2 - Q)*PolyLog[2, 1 - U/mu2])/(Q - U)^2)) +
    Xt^4*(Nc^2*(-8/(Q - U)^2 - (4*Q*(Q + U)*Log[Q]^2)/(Q - U)^4 -
        (4*(Q + 3*U)*Log[U])/(Q - U)^3 - (4*U*(Q + U)*Log[U]^2)/(Q - U)^4 +
        Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/(Q - U)^3) +
        Log[Q]*((4*(3*Q + U))/(Q - U)^3 + (4*(Q + U)*Log[q2])/(Q - U)^3 +
          (4*(Q + U)^2*Log[U])/(Q - U)^4)) +
      Nc*((2*mu2^3*(2*Q^2 + 8*Q*U - U^2) - mu2^2*(6*Q^3 + 38*Q^2*U +
            13*Q*U^2 - 3*U^3) + Q*U*(-2*Q^3 - 24*Q^2*U + 7*Q*U^2 + U^3) +
          mu2*(2*Q^4 + 28*Q^3*U + 31*Q^2*U^2 - 6*Q*U^3 - U^4))/
         ((mu2 - Q)*Q*(mu2 - U)*(Q - U)^3*U) -
        (2*Q*(mu2^2*(4*Q + 3*U) - 2*mu2*(3*Q^2 + 3*Q*U + U^2) +
           Q*(3*Q^2 + 3*Q*U + U^2))*Log[Q]^2)/((mu2 - Q)^2*(Q - U)^4) +
        ((12*mu2^4*Q*(2*Q + U) + Q*U^2*(25*Q^3 + 34*Q^2*U - 21*Q*U^2 -
             2*U^3) - mu2^3*(39*Q^3 + 102*Q^2*U + 5*Q*U^2 - 2*U^3) +
           mu2^2*(19*Q^4 + 130*Q^3*U + 93*Q^2*U^2 - 22*Q*U^3 - 4*U^4) +
           mu2*U*(-48*Q^4 - 113*Q^3*U - 6*Q^2*U^2 + 21*Q*U^3 + 2*U^4))*
          Log[U])/(2*(mu2 - Q)*Q*(mu2 - U)^2*(Q - U)^4) -
        (U*(Q + U)*(4*mu2^2 - 2*mu2*(Q + 3*U) + U*(Q + 3*U))*Log[U]^2)/
         ((mu2 - U)^2*(Q - U)^4) + Log[q2]*
         ((-2*Q^3 - 13*Q^2*U + 2*Q*U^2 + U^3 + 2*mu2*(2*Q^2 + 5*Q*U - U^2))/
           (Q*(Q - U)^3*U) + (3*(4*mu2*Q - 3*Q^2 - 2*Q*U + U^2)*Log[U])/
           (Q - U)^4) + Log[Q]*((-36*mu2^4*Q*U - Q^2*U*(4*Q^3 + 47*Q^2*U +
              4*Q*U^2 - 19*U^3) + mu2^3*(4*Q^3 + 121*Q^2*U + 24*Q*U^2 -
              5*U^3) + mu2*Q*(4*Q^4 + 53*Q^3*U + 130*Q^2*U^2 - 9*Q*U^3 -
              34*U^4) + mu2^2*(-8*Q^4 - 126*Q^3*U - 131*Q^2*U^2 + 42*Q*U^3 +
              7*U^4))/(2*(mu2 - Q)^2*(mu2 - U)*(Q - U)^4*U) -
          (3*(4*mu2*Q - 3*Q^2 - 2*Q*U + U^2)*Log[q2])/(Q - U)^4 +
          ((10*Q)/(Q - U)^3 - 2/(Q - U)^2 + (22*Q*U)/(Q - U)^4 +
            (6*U)/(-Q + U)^3 + (2*(2*mu2 - Q)*Q*(Q + U))/((mu2 - Q)^2*
              (Q - U)^3) + ((2*mu2 - U)*U*(Q + U))/((mu2 - U)^2*(-Q + U)^3))*
           Log[U]) + Log[mu2]*((-2*mu2*(-3*mu2*Q^3*U*(Q + 3*U) +
             mu2^4*(2*Q^2 + 2*Q*U - U^2) + Q^2*U^2*(2*Q^2 + 2*Q*U - U^2) -
             mu2^3*(4*Q^3 + 5*Q^2*U + 5*Q*U^2 - 2*U^3) +
             mu2^2*(2*Q^4 + 8*Q^3*U + 7*Q^2*U^2 + 2*Q*U^3 - U^4)))/
           ((mu2 - Q)^2*Q*(mu2 - U)^2*(Q - U)^3*U) +
          ((12*mu2^5*Q + 2*Q^2*U^2*(-2*Q^2 - 2*Q*U + U^2) -
             mu2^4*(25*Q^2 + 28*Q*U + U^2) + 6*mu2^3*Q*(3*Q^2 + 10*Q*U + 3*
                U^2) - mu2^2*Q*(3*Q^3 + 44*Q^2*U + 41*Q*U^2 - 4*U^3) +
             4*mu2*Q*U*(2*Q^3 + 7*Q^2*U + Q*U^2 - U^3))*Log[Q])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^4) +
          ((-12*mu2^5*Q + 2*Q^2*U^2*(2*Q^2 + 2*Q*U - U^2) +
             mu2^4*(25*Q^2 + 28*Q*U + U^2) - 6*mu2^3*Q*(3*Q^2 + 10*Q*U + 3*
                U^2) + mu2^2*Q*(3*Q^3 + 44*Q^2*U + 41*Q*U^2 - 4*U^3) -
             4*mu2*Q*U*(2*Q^3 + 7*Q^2*U + Q*U^2 - U^3))*Log[U])/
           ((mu2 - Q)^2*(mu2 - U)^2*(Q - U)^4)) +
        (2*(3*mu2^2 - 6*mu2*Q + 2*Q^2 + 2*Q*U - U^2)*PolyLog[2,
           (mu2 - Q)/mu2])/(Q - U)^4 + (2*(-3*mu2^2 + 6*mu2*Q - 2*Q^2 -
           2*Q*U + U^2)*PolyLog[2, 1 - U/mu2])/(Q - U)^4)))/sbe^2) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
];


    
(* Top and bottom Yukawa lambda 2-loop threshold correction, all masses degenerate and equal to the matching scale *)
lambda2LHSSAlphaTAlphaBAllDegenerate = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    CF = 4/3,
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    gtau = Ye[3,3],
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb,
    M  = Sqrt[msq2[3,3]*msu2[3,3]],
    K = -0.1953256,
    sgn = MuInput/Abs[MuInput],
    (* For the bottom Yukawa coupling *)
    Mu = MuInput,
    M3 = M3Input,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    Q2 = SCALE^2,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];
    k^2*( ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^6*(5*Nc - (14*Nc*Xb^2)/M + ((11*Nc)/(2*M^2) - (2*Nc^2)/(3*M^2))*Xb^4 +
    (-Nc/(2*M^3) + Nc^2/(18*M^3))*Xb^6 +
    (4*Nc - (2*Nc*Xb^2)/M + (Nc*Xb^4)/M^2)/cbe^2 +
    (sbe^2*(Nc*(-3/2 + 60*K + Pi^2) + (-12/M - (64*K)/M)*Nc*Xb*Yb +
       (4/M^2 + (16*K)/M^2)*Nc*Xb^3*Yb + (-3/M - (16*K)/M)*Nc*Yb^2 +
       Xb^4*(-Nc/(2*M^2) + (-19/(12*M^3) - (8*K)/M^3)*Nc*Yb^2) +
       Xb^2*((-2/M - (24*K)/M)*Nc + (14/(3*M^2) + (24*K)/M^2)*Nc*Yb^2)))/
     cbe^2) + ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*gt^2*(Nc*(-1 + 48*K + (4*Pi^2)/3) +
    (-10/M - (32*K)/M)*Nc*Xb^2 + (Nc*Xb^4)/(6*M^2) +
    ((4*Nc)/M + (-(Nc/M^2) - (2*Nc^2)/(3*M^2))*Xb^2 +
      (Nc/(6*M^3) + Nc^2/(18*M^3))*Xb^4)*Xt^2 +
    ((-2*Nc*sgn*Xb^3)/(3*M^(3/2)) + (4*Nc*sgn*Xt)/Sqrt[M])/(cbe*sbe) +
    (-8*Nc - (Nc*Xb^2)/M + (-(Nc/M) + (Nc*Xb^2)/(2*M^2))*Xt^2)/cbe^2 +
    (-4/M - (32*K)/M)*Nc*Xb*Yb + (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xb^3*Yb +
    (sbe^2*(24*K*Nc + (16*K*Nc*Xb*Yb)/M + (8*K*Nc*Yb^2)/M +
       Xt*((16*K*Nc*Xb)/M + (16*K*Nc*Yb)/M + (2*Nc*Xb^2*Yb)/(3*M^2)) +
       Xb^2*((-M^(-1) - (8*K)/M)*Nc + (Nc*Yb^2)/(3*M^2)) +
       Xt^2*((-M^(-1) - (8*K)/M)*Nc + (2*Nc*Xb*Yb)/(3*M^2) +
         (Nc*Yb^2)/(3*M^2) + (-11/(18*M^3) - (8*K)/(3*M^3))*Nc*Xb^2*Yb^2)))/
     cbe^2 + ((-4/M - (32*K)/M)*Nc*Xb + (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xb^3)*
     Yt + Xt*((16*K*Nc*Xb)/M + (16*K*Nc*Yb)/M + (2*Nc*Xb^2*Yb)/(3*M^2) +
      ((16*K*Nc)/M + (2*Nc*Xb^2)/(3*M^2) + (4*Nc*Xb*Yb)/(3*M^2) +
        (2/(27*M^3) + (16*K)/(9*M^3))*Nc*Xb^3*Yb)*Yt) +
    (4*Nc + (Nc*Xb^4)/(6*M^2) + cbe^2*(Nc*(12*K + (-3 + 2*Pi^2)/6) -
        (8*K*Nc*Xb^2)/M - (Nc*Xb^4)/(6*M^2) + ((-4/M - (32*K)/M)*Nc*Xb +
          (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xb^3)*Yt +
        ((-M^(-1) - (8*K)/M)*Nc + (4/(3*M^2) + (8*K)/M^2)*Nc*Xb^2 +
          (-35/(108*M^3) - (16*K)/(9*M^3))*Nc*Xb^4)*Yt^2))/sbe^2) +
     + ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*gt^4*(Nc*(-1 + 48*K + (4*Pi^2)/3) + (4*Nc*Xb^2)/M +
    (-Nc/(2*M^2) + (Nc*Xb^2)/(6*M^3))*Xt^4 + (4*Nc + (Nc*Xt^4)/(6*M^2))/
     cbe^2 + ((4*Nc*sgn*Xb)/Sqrt[M] - (4*Nc*sgn*Xb*Xt^2)/M^(3/2) -
      (2*Nc*sgn*Xt^3)/(3*M^(3/2)) + (Nc*sgn*Xb*Xt^4)/(3*M^(5/2)))/(cbe*sbe) +
    (16*K*Nc*Xb*Yb)/M + (sbe^2*(Nc*(12*K + (-3 + 2*Pi^2)/6) +
       (-4/M - (32*K)/M)*Nc*Xt*Yb + (8/(9*M^2) + (16*K)/(3*M^2))*Nc*Xt^3*Yb +
       (-M^(-1) - (8*K)/M)*Nc*Yb^2 + Xt^4*(-Nc/(4*M^2) +
         (-35/(108*M^3) - (16*K)/(9*M^3))*Nc*Yb^2) +
       Xt^2*((M^(-1) - (8*K)/M)*Nc + (4/(3*M^2) + (8*K)/M^2)*Nc*Yb^2)))/
     cbe^2 + (16*K*Nc*Xb*Yt)/M + Xt^2*((-2/M - (32*K)/M)*Nc - (Nc*Xb^2)/M^2 +
      (2*Nc*Xb*Yb)/(3*M^2) + (2*Nc*Xb*Yt)/(3*M^2)) +
    Xt^3*((8/(9*M^2) + (16*K)/(3*M^2))*Nc*Yb +
      ((8/(9*M^2) + (16*K)/(3*M^2))*Nc + (2/(27*M^3) + (16*K)/(9*M^3))*Nc*Xb*
         Yb)*Yt) + Xt*((16*K*Nc*Xb)/M + (-4/M - (32*K)/M)*Nc*Yb +
      ((-4/M - (32*K)/M)*Nc + (4*Nc*Xb*Yb)/(3*M^2))*Yt) +
    (-8*Nc - (Nc*Xb^2)/M + (-(Nc/M) + (Nc*Xb^2)/(2*M^2))*Xt^2 +
      cbe^2*(24*K*Nc + (-M^(-1) - (8*K)/M)*Nc*Xb^2 + (16*K*Nc*Xb*Yt)/M +
        ((8*K*Nc)/M + (Nc*Xb^2)/(3*M^2))*Yt^2 +
        Xt*((16*K*Nc*Xb)/M + ((16*K*Nc)/M + (2*Nc*Xb^2)/(3*M^2))*Yt) +
        Xt^2*((-M^(-1) - (8*K)/M)*Nc + (2*Nc*Xb*Yt)/(3*M^2) +
          (Nc/(3*M^2) + (-11/(18*M^3) - (8*K)/(3*M^3))*Nc*Xb^2)*Yt^2)))/
     sbe^2))
];

(* Top and bottom Yukawa lambda 2-loop threshold correction, computation with general masses *)
lambda2LHSSAlphaTAlphaBGeneric = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.02),
    U  = msu2[3,3]*(1-0.02),
    mD  = msd2[3,3]*(1-0.04),
    mu2 = MuInput^2*(1+0.03),
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    sgn = MuInput/Abs[MuInput],      
    q2 = SCALE^2, (* renormalization/matching scale *)
    Q2 = SCALE^2,
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3,
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    gtau = Ye[3,3],
    A0 = mAInput^2*(1-0.02),
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
      ];
k^2*(ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^6*(Nc*(5 + 2*Log[mD]^2 + 3*Log[Q]^2 + Log[Q]*(-6 - 6*Log[q2]) + 
      Log[mD]*(-4 - 4*Log[q2]) + 10*Log[q2] + 5*Log[q2]^2) + 
    (sbe^2*(Nc*Yb^2*(-2/mD - Q^(-1) + ((-(A0^3*(mD + 2*Q)) + 
            (mD - Q)^2*(mD^2 - 4*mD*Q + 2*Q^2) + A0^2*(5*mD^2 + 6*mD*Q + 
              6*Q^2) + A0*(-5*mD^3 + 8*mD^2*Q + 3*mD*Q^2 - 6*Q^3) + 
            ((A0 - mD)*mD + 2*(A0 + 2*mD)*Q - 2*Q^2)*delta[A0, Q, mD])*
           Log[A0])/(2*mD^2*Q*delta[A0, Q, mD]) + 
         ((A0^3 + A0^2*(-5*mD + Q) + 5*A0*(mD^2 - 2*mD*Q - Q^2) - 
            (mD - Q)*(mD^2 - 6*mD*Q + 3*Q^2) - (A0 - 3*mD + 3*Q)*
             delta[A0, Q, mD])*Log[mD])/(2*mD*Q*delta[A0, Q, mD]) - 
         ((-2*A0^3 + A0^2*(7*mD + 6*Q) + (mD - Q)*(mD^2 + 3*mD*Q - 2*Q^2) - 
            2*A0*(mD^2 + mD*Q + 3*Q^2) + (2*A0 - 3*mD - 2*Q)*
             delta[A0, Q, mD])*Log[Q])/(2*mD^2*delta[A0, Q, mD]) + 
         (-2/mD - Q^(-1))*Log[q2] + ((2*A0^4 - 8*A0^3*(mD + Q) - 
            2*A0*(mD - Q)*(3*mD^2 - 4*Q^2) + (mD - Q)^2*(mD^2 - 4*mD*Q + 
              2*Q^2) + A0^2*(7*mD^2 + 8*mD*Q + 12*Q^2) + delta[A0, Q, mD]*
             (-6*A0^2 - 5*mD^2 + 12*mD*Q - 6*Q^2 + 12*A0*(mD + Q) + 
              4*delta[A0, Q, mD]))*phi[A0, Q, mD])/
          (2*mD^3*delta[A0, Q, mD])) + Nc*(3/2 - (2*A0)/mD + Pi^2 - A0/Q + 
         3*Log[A0]^2 + 2*Log[mD]^2 + 3*Log[Q]^2 + Log[Q]*(-6 - 6*Log[q2]) + 
         Log[A0]*(7 + A0*(2/mD + Q^(-1)) - 6*Log[q2]) + 
         Log[mD]*(-4 - 4*Log[q2]) + (3 - A0*(2/mD + Q^(-1)))*Log[q2] + 
         8*Log[q2]^2 + ((A0^2 - 6*A0*mD - delta[A0, mD, mD])*phi[A0, mD, mD])/
          mD^2 + ((2*A0^2 - 11*A0*Q - 2*delta[A0, Q, Q])*phi[A0, Q, Q])/
          Q^2) + Nc*Xb*Yb*((4*Log[mD]^2)/(mD - Q) + (6*Log[Q]^2)/(-mD + Q) + 
         Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
         Log[Q]*(12/(mD - Q) + (12*Log[q2])/(mD - Q)) + 
         Log[mD]*(12/(-mD + Q) + (2*Log[Q])/(mD - Q) + (12*Log[q2])/
            (-mD + Q)) - (2*(-(A0*(A0 - 6*mD)) + delta[A0, mD, mD])*
           phi[A0, mD, mD])/(mD^2*(mD - Q)) + 
         (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
            delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)) + 
         ((-4*A0^2 + 22*A0*Q + 4*delta[A0, Q, Q])*phi[A0, Q, Q])/
          ((mD - Q)*Q^2)) + Nc*Xb^3*Yb*(-48/(mD - Q)^2 + 
         (4*(A0 - mD - Q)*Log[mD]^2)/(mD - Q)^3 + 
         (2*(-4*A0 + 3*mD + 5*Q)*Log[Q]^2)/(mD - Q)^3 + 
         Log[A0]*((-2*(6*A0 + mD - Q)*Log[mD])/(mD - Q)^3 + 
           (2*(6*A0 + mD - Q)*Log[Q])/(mD - Q)^3) - (24*Log[q2])/(mD - Q)^2 + 
         Log[Q]*((-12*(mD + 3*Q))/(mD - Q)^3 - (12*(mD + Q)*Log[q2])/
            (mD - Q)^3) + Log[mD]*((12*(3*mD + Q))/(mD - Q)^3 - 
           (2*(-2*A0 + mD + 3*Q)*Log[Q])/(mD - Q)^3 + (12*(mD + Q)*Log[q2])/
            (mD - Q)^3) + (2*(A0*(A0 - 6*mD)*(mD - Q) + (-3*mD + Q)*
             delta[A0, mD, mD])*phi[A0, mD, mD])/(mD^2*(mD - Q)^3) + 
         (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
            (-3*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)^3) + 
         ((2*A0*(2*A0 - 11*Q)*(mD - Q) - 4*(mD - 3*Q)*delta[A0, Q, Q])*
           phi[A0, Q, Q])/((mD - Q)^3*Q^2)) + 
       Xb^4*(Nc*((A0*(-mD^2 + 5*mD*Q + 2*Q^2))/(mD*(mD - Q)^3*Q) + 
           Log[A0]*((A0*(mD^2 - 5*mD*Q - 2*Q^2))/(mD*(mD - Q)^3*Q) + 
             (6*A0*Q*Log[mD])/(mD - Q)^4 - (6*A0*Q*Log[Q])/(mD - Q)^4) + 
           (A0*(-mD^2 + 5*mD*Q + 2*Q^2)*Log[q2])/(mD*(mD - Q)^3*Q) + 
           Log[mD]*((-6*A0*Q)/(mD - Q)^4 - (6*A0*Q*Log[q2])/(mD - Q)^4) + 
           Log[Q]*((6*A0*Q)/(mD - Q)^4 + (6*A0*Q*Log[q2])/(mD - Q)^4)) + 
         Nc*Yb^2*((-mD^2 + 11*mD*Q + 2*Q^2)/(mD*(mD - Q)^3*Q) - 
           (2*(-2*A0 + 3*mD + Q)*Log[mD]^2)/(mD - Q)^4 + 
           ((7*A0 - 3*mD - 11*Q)*Log[Q]^2)/(mD - Q)^4 + 
           Log[A0]*((A0^3*(mD - 2*Q) - (mD - Q)^2*(3*mD^2 + 4*mD*Q - 2*Q^2) + 
               A0*(mD - Q)*(5*mD^2 + mD*Q + 6*Q^2) + A0^2*(-3*mD^2 + 2*mD*Q + 
                 6*Q^2) + (-(A0*mD) + 3*mD^2 + 2*A0*Q + 4*mD*Q - 2*Q^2)*
                delta[A0, Q, mD])/(2*mD^2*(mD - Q)^2*Q*delta[A0, Q, mD]) + 
             (3*(A0 + mD - Q)*Log[mD])/(mD - Q)^4 - (3*(A0 + mD - Q)*Log[Q])/
              (mD - Q)^4) + ((-mD^2 + 5*mD*Q + 2*Q^2)*Log[q2])/
            (mD*(mD - Q)^3*Q) + Log[mD]*(((mD - Q)^2*(-A0^3 + 
                 A0^2*(3*mD + 7*Q) + (mD - Q)*(3*mD^2 + 12*mD*Q - 5*Q^2) - 
                 A0*(5*mD^2 + 10*mD*Q + 11*Q^2)) - (mD^3 - A0*(mD - Q)^2 + 
                 15*mD^2*Q + 3*mD*Q^2 + 5*Q^3)*delta[A0, Q, mD])/
              (2*mD*(mD - Q)^4*Q*delta[A0, Q, mD]) + ((-11*A0 + 9*mD + 13*Q)*
               Log[Q])/(mD - Q)^4 - (6*Q*Log[q2])/(mD - Q)^4) + 
           Log[Q]*(((mD - Q)^2*(2*A0^3 - 3*A0^2*(3*mD + 2*Q) - (mD - Q)*
                  (11*mD^2 + mD*Q - 2*Q^2) + 2*A0*(7*mD^2 + 3*mD*Q + 
                   3*Q^2)) + (5*mD^3 - 2*A0*(mD - Q)^2 + 16*mD^2*Q + mD*Q^2 + 
                 2*Q^3)*delta[A0, Q, mD])/(2*mD^2*(mD - Q)^4*delta[A0, Q, 
                mD]) + (6*Q*Log[q2])/(mD - Q)^4) + 
           ((A0*(A0 - 6*mD)*(mD - Q) + (-5*mD + Q)*delta[A0, mD, mD])*
             phi[A0, mD, mD])/(mD^2*(mD - Q)^4) + 
           (((mD - Q)^2*(2*A0^4 - 8*A0^3*(mD + Q) - 2*A0*(mD - Q)*
                 (3*mD^2 - 4*Q^2) - (mD - Q)^2*(3*mD^2 + 4*mD*Q - 2*Q^2) + 
                A0^2*(11*mD^2 + 8*mD*Q + 12*Q^2)) - (mD - Q)*(2*A0^2*
                 (8*mD - 3*Q) - 4*A0*(9*mD^2 + 5*mD*Q - 3*Q^2) + 
                (mD - Q)*(27*mD^2 - 22*mD*Q + 6*Q^2))*delta[A0, Q, mD] + 
              2*(18*mD^2 - 9*mD*Q + 2*Q^2)*delta[A0, Q, mD]^2)*
             phi[A0, Q, mD])/(2*mD^3*(mD - Q)^4*delta[A0, Q, mD]) + 
           ((A0*(2*A0 - 11*Q)*(-mD + Q) + (2*mD - 9*Q)*delta[A0, Q, Q])*
             phi[A0, Q, Q])/((mD - Q)^4*Q^2))) + 
       Xb^2*(Nc*Yb^2*((2*mD - 4*Q)/(mD^2*Q - mD*Q^2) + (2*Log[mD]^2)/
            (mD - Q)^2 + (3*Log[Q]^2)/(mD - Q)^2 + 
           Log[A0]*(-(((A0 + mD - Q)*(2*A0^2*Q + A0*(mD^2 - 6*mD*Q - 4*Q^2) - 
                  (mD - Q)*(mD^2 - 4*mD*Q + 2*Q^2)) + (mD^2 - 4*mD*Q + 
                  2*Q*(-A0 + Q))*delta[A0, Q, mD])/(mD^2*(mD - Q)*Q*
                delta[A0, Q, mD])) - Log[mD]/(mD - Q)^2 + 
             Log[Q]/(mD - Q)^2) + ((2*mD - 4*Q)*Log[q2])/(mD^2*Q - mD*Q^2) + 
           Log[Q]*((-2*(mD - Q)*(-A0^3 + mD^3 - 2*mD*Q^2 + Q^3 + 
                 A0^2*(4*mD + 3*Q) - A0*(6*mD^2 + 2*mD*Q + 3*Q^2)) - 2*
                (-mD^2 + A0*(mD - Q) + mD*Q + Q^2)*delta[A0, Q, mD])/
              (mD^2*(mD - Q)^2*delta[A0, Q, mD]) - (2*Log[q2])/(mD - Q)^2) + 
           Log[mD]*(((mD - Q)*(A0^2*(mD + 4*Q) - 2*A0*Q*(7*mD + 4*Q) - 
                 (mD - Q)*(mD^2 - 7*mD*Q + 4*Q^2)) - (mD^2 + mD*Q - 4*Q^2)*
                delta[A0, Q, mD])/(mD*(mD - Q)^2*Q*delta[A0, Q, mD]) - 
             (5*Log[Q])/(mD - Q)^2 + (2*Log[q2])/(mD - Q)^2) + 
           ((A0*(A0 - 6*mD) - delta[A0, mD, mD])*phi[A0, mD, mD])/
            (mD^2*(mD - Q)^2) + (((mD - Q)*(2*A0^4 - 8*A0^3*(mD + Q) - 
                2*A0*(mD - 2*Q)*(mD - Q)*(mD + 2*Q) + (mD - Q)^2*
                 (mD^2 - 4*mD*Q + 2*Q^2) + A0^2*(11*mD^2 + 8*mD*Q + 
                  12*Q^2)) + delta[A0, Q, mD]*(-3*(3*mD - 2*Q)*(mD - Q)^2 + 
                A0^2*(-9*mD + 6*Q) + A0*(23*mD^2 + 6*mD*Q - 12*Q^2) + 
                (7*mD - 4*Q)*delta[A0, Q, mD]))*phi[A0, Q, mD])/
            (mD^3*(mD - Q)^2*delta[A0, Q, mD]) + 
           ((A0*(2*A0 - 11*Q) - 2*delta[A0, Q, Q])*phi[A0, Q, Q])/
            ((mD - Q)^2*Q^2)) + Nc*((2*A0*mD - 4*A0*Q + 4*mD*Q)/
            (mD^2*Q - mD*Q^2) + (2*Log[mD]^2)/(mD - Q) + 
           ((A0 - 3*mD + Q)*Log[Q]^2)/(mD - Q)^2 + 
           Log[A0]*((-2*A0*(mD - 2*Q))/(mD*(mD - Q)*Q) - 
             ((A0 - mD + Q)*Log[mD])/(mD - Q)^2 + ((A0 - mD + Q)*Log[Q])/
              (mD - Q)^2) + ((2*A0*mD - 4*A0*Q + 2*mD*Q)*Log[q2])/
            (mD^2*Q - mD*Q^2) + Log[Q]*((-2*(A0 - 3*mD + Q))/(mD - Q)^2 - 
             (2*(A0 - 3*mD + 2*Q)*Log[q2])/(mD - Q)^2) + 
           Log[mD]*((2*(A0 - 4*mD + 2*Q))/(mD - Q)^2 + 
             ((-A0 + mD + Q)*Log[Q])/(mD - Q)^2 + (2*(A0 - 3*mD + 2*Q)*Log[
                q2])/(mD - Q)^2) + ((A0*(A0 - 6*mD) - delta[A0, mD, mD])*
             phi[A0, mD, mD])/(mD^2*(mD - Q)) + (delta[A0, Q, mD]*
             phi[A0, Q, mD])/(mD*(mD - Q)^2) + 
           ((A0*(2*A0 - 11*Q)*(-mD + Q) + (2*mD - 3*Q)*delta[A0, Q, Q])*
             phi[A0, Q, Q])/((mD - Q)^2*Q^2)))))/cbe^2 + 
    Xb^4*(Nc^2*((4*mD*Q*Log[mD]^2)/(mD - Q)^4 + (2*(mD + Q)*Log[Q])/
         (mD - Q)^3 + (4*mD*Q*Log[Q]^2)/(mD - Q)^4 + 
        Log[mD]*((-2*(mD + Q))/(mD - Q)^3 - (8*mD*Q*Log[Q])/(mD - Q)^4)) + 
      Nc*((2*(mD^2 - 27*mD*Q + 2*Q^2))/(mD*(mD - Q)^2*Q) - 
        (2*(mD + 3*Q)*Log[mD]^2)/(mD - Q)^3 + (4*(3*mD + Q)*Log[Q]^2)/
         (mD - Q)^3 + (2*(mD^2 - 15*mD*Q + 2*Q^2)*Log[q2])/
         (mD*(mD - Q)^2*Q) + Log[Q]*((4*(-3*mD^2 - 10*mD*Q + Q^2))/
           (mD*(mD - Q)^3) - (2*(7*mD + 5*Q)*Log[q2])/(mD - Q)^3) + 
        Log[mD]*((-2*mD^2 + 44*mD*Q + 6*Q^2)/((mD - Q)^3*Q) + 
          (2*(-5*mD + Q)*Log[Q])/(mD - Q)^3 + (2*(7*mD + 5*Q)*Log[q2])/
           (mD - Q)^3) + (6*PolyLog[2, 1 - mD/Q])/(mD - Q)^2)) + 
    Xb^6*(Nc^2*((-2*(mD + Q))/(mD - Q)^4 - (2*mD*Q*(mD + Q)*Log[mD]^2)/
         (mD - Q)^6 - ((mD^2 + 6*mD*Q + Q^2)*Log[Q])/(mD - Q)^5 - 
        (2*mD*Q*(mD + Q)*Log[Q]^2)/(mD - Q)^6 + 
        Log[mD]*((mD^2 + 6*mD*Q + Q^2)/(mD - Q)^5 + (4*mD*Q*(mD + Q)*Log[Q])/
           (mD - Q)^6)) + Nc*((-mD^2 + 11*mD*Q + 2*Q^2)/(mD*(mD - Q)^3*Q) - 
        (2*(3*mD + Q)*Log[mD]^2)/(mD - Q)^4 - (2*(3*mD + 4*Q)*Log[Q]^2)/
         (mD - Q)^4 + ((-mD^2 + 5*mD*Q + 2*Q^2)*Log[q2])/(mD*(mD - Q)^3*Q) + 
        Log[mD]*(((mD - 5*Q)*(mD + 2*Q))/((mD - Q)^4*Q) + 
          (2*(6*mD + 5*Q)*Log[Q])/(mD - Q)^4 - (6*Q*Log[q2])/(mD - Q)^4) + 
        Log[Q]*((-3*mD^2 + 13*mD*Q + 2*Q^2)/(mD*(mD - Q)^4) + 
          (6*Q*Log[q2])/(mD - Q)^4) + (6*PolyLog[2, 1 - mD/Q])/
         (-mD + Q)^3)) + Xb^2*(Nc^2*((2*mD*Q*Log[mD]^2)/(mD - Q)^3 - 
        (2*mD*Q*Log[Q]^2)/(mD - Q)^3 + (2*(mD + Q)*Log[q2])/(mD - Q)^2 + 
        Log[mD]*(-((mD + Q)/(mD - Q)^2) - (4*mD*Q*Log[q2])/(mD - Q)^3) + 
        Log[Q]*(-((mD + Q)/(mD - Q)^2) + (4*mD*Q*Log[q2])/(mD - Q)^3)) + 
      Nc*(-2/mD + 4/(mD - Q) - Q^(-1) + (6*Log[mD]^2)/(mD - Q) - 
        (2*(6*mD - 5*Q)*Log[Q]^2)/(mD - Q)^2 + (-2/mD + 2/(mD - Q) - Q^(-1))*
         Log[q2] + Log[mD]*((mD^2 - 21*mD*Q + 16*Q^2)/((mD - Q)^2*Q) + 
          ((6*mD - 4*Q)*Log[Q])/(mD - Q)^2 - (2*(9*mD - 8*Q)*Log[q2])/
           (mD - Q)^2) + Log[Q]*((19*mD^2 - 17*mD*Q + 2*Q^2)/
           (mD*(mD - Q)^2) + (2*(9*mD - 8*Q)*Log[q2])/(mD - Q)^2) + 
        (6*PolyLog[2, 1 - mD/Q])/(-mD + Q))) + 
    (Nc*(-5/2 + mD*((mD - mu2)^(-1) - Q^(-1)) + (4*mu2 - 2*Q)/mD + 
        2*mu2*(Q^(-1) + (-mu2 + Q)^(-1)) + ((mD^2 - 2*mD*mu2 - mu2^2)*
          Log[mD]^2)/(mD - mu2)^2 + (1 - (4*mu2^2)/(mu2 - Q)^2)*Log[Q]^2 + 
        (2 + (4*mu2)/mD + (2*mD)/(-mD + mu2) + (4*mu2)/(mu2 - Q) - 
          (mD - 2*mu2)/Q - (2*Q)/mD)*Log[q2] + 4*Log[q2]^2 + 
        Log[mD]*((-mD^2 + 2*mD*mu2 + 2*mu2^2)/(mD - mu2)^2 + mD/Q + 
          (2*mu2^2*Log[mu2])/(mD - mu2)^2 + 2*Log[Q] - 
          (2*(2*mD^2 - 4*mD*mu2 + mu2^2)*Log[q2])/(mD - mu2)^2) + 
        Log[mu2]*(mu2*((-2*mD^2 + 3*mD*mu2 - 4*mu2^2)/(mD*(mD - mu2)^2) - 
            (6*mu2)/(mu2 - Q)^2 - 2/Q) + (4*mu2^2*Log[Q])/(mu2 - Q)^2 + 
          2*mu2^2*(-(mD - mu2)^(-2) - 2/(mu2 - Q)^2)*Log[q2]) + 
        Log[Q]*((2*Q)/mD + (3*mu2^2 + 2*mu2*Q + Q^2)/(mu2 - Q)^2 + 
          (-4 + (4*mu2^2)/(mu2 - Q)^2)*Log[q2]) + 
        (2*(mD^2 - 2*mD*mu2 - mu2^2)*PolyLog[2, 1 - mu2/mD])/(mD - mu2)^2 + 
        (2 - (8*mu2^2)/(mu2 - Q)^2)*PolyLog[2, 1 - mu2/Q]) + 
      Xb^2*(Nc^2*((4*mD*Log[mD]^2)/(mD - Q)^2 + (4*Q*Log[Q]^2)/(mD - Q)^2 + 
          Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q)) + 
          Log[mD]*(4/(-mD + Q) - (4*(mD + Q)*Log[Q])/(mD - Q)^2 + 
            (4*Log[q2])/(-mD + Q))) + 
        Nc*((2*mD^2 - 4*mD*mu2 - 6*mD*Q + 8*mu2*Q - 4*Q^2)/
           (mD^2*Q - mD*Q^2) + (2*(2*mD + mu2 - Q)*Log[mD]^2)/(mD - Q)^2 + 
          4*mu2*(1/((mD - mu2)*(mu2 - Q)) + (Q^(-1) + (-mD + Q)^(-1))/mD)*
           Log[mu2] + ((-2*mu2 + 8*Q)*Log[Q]^2)/(mD - Q)^2 + 
          (2*(mD^2 + 2*(2*mu2 - Q)*Q - 2*mD*(mu2 + Q))*Log[q2])/
           (mD*(mD - Q)*Q) + Log[Q]*((-2*mD*mu2*(mD + 2*mu2) + 
              2*mD*(3*mD + 5*mu2)*Q + 2*(-5*mD + 2*mu2)*Q^2 - 4*Q^3)/
             (mD*(mD - Q)^2*(-mu2 + Q)) + ((4*(mD + mu2) - 10*Q)*Log[q2])/
             (mD - Q)^2) + Log[mD]*((-4*(mD^2 + mD*(mu2 - 3*Q) - mu2*
                (mu2 - 2*Q)))/((mD - mu2)*(mD - Q)^2) - 2/Q - 
            (2*(2*mD + 3*Q)*Log[Q])/(mD - Q)^2 + ((-4*(mD + mu2) + 10*Q)*
              Log[q2])/(mD - Q)^2) + (4*(mu2 - Q)*PolyLog[2, 1 - mu2/mD])/
           (mD - Q)^2 + (4*(-mu2 + Q)*PolyLog[2, 1 - mu2/Q])/(mD - Q)^2)) + 
      Xb^4*(Nc^2*(-8/(mD - Q)^2 - (4*mD*(mD + Q)*Log[mD]^2)/(mD - Q)^4 - 
          (4*Q*(mD + Q)*Log[Q]^2)/(mD - Q)^4 - (8*Log[q2])/(mD - Q)^2 + 
          Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - (4*(mD + Q)*Log[q2])/
             (mD - Q)^3) + Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 + 
            (4*(mD + Q)^2*Log[Q])/(mD - Q)^4 + (4*(mD + Q)*Log[q2])/
             (mD - Q)^3)) + Nc*((-mD^3 + 2*mD^2*(mu2 - 5*Q) + 
            2*Q^2*(-2*mu2 + Q) + mD*Q*(-16*mu2 + 27*Q))/(mD*(mD - Q)^3*Q) + 
          ((-5*mD^2 + 3*mu2^2 - 2*(mD + 3*mu2)*Q + 2*Q^2)*Log[mD]^2)/
           (mD - Q)^4 + ((mD^2 - 3*mu2^2 - 8*mD*Q + 6*mu2*Q - 10*Q^2)*
            Log[Q]^2)/(mD - Q)^4 + Log[mu2]*((-2*mu2*(mD^2 - 2*mD*Q - 2*Q^2))/
             (mD*(mD - Q)^3*Q) + (6*mu2^2*Log[Q])/(mD - Q)^4) + 
          ((-mD^3 + 2*mD^2*(mu2 - 4*Q) + 2*Q^2*(-2*mu2 + Q) + 
             mD*Q*(-10*mu2 + 19*Q))*Log[q2])/(mD*(mD - Q)^3*Q) + 
          Log[Q]*((-7*mD^3 - 6*mD*(mD + 3*mu2)*Q + 29*mD*Q^2 + 2*Q^3)/
             (mD*(mD - Q)^4) - (6*(mD^2 - mD*Q + 2*(mu2 - Q)*Q)*Log[q2])/
             (mD - Q)^4) + Log[mD]*((mD^3 + 14*mD^2*Q + 3*mD*(2*mu2 - 7*Q)*
               Q + 12*(mu2 - Q)*Q^2)/((mD - Q)^4*Q) - (6*mu2^2*Log[mu2])/
             (mD - Q)^4 + (2*(2*mD^2 + 5*mD*Q + 4*Q^2)*Log[Q])/(mD - Q)^4 + 
            (6*(mD^2 - mD*Q + 2*(mu2 - Q)*Q)*Log[q2])/(mD - Q)^4) + 
          ((-2*mD^2 + 6*mu2^2 + 4*mD*Q + 4*Q*(-3*mu2 + Q))*
            PolyLog[2, 1 - mu2/mD])/(mD - Q)^4 + 
          (2*(mD^2 - 3*mu2^2 - 2*mD*Q + 6*mu2*Q - 2*Q^2)*
            PolyLog[2, 1 - mu2/Q])/(mD - Q)^4)))/cbe^2) + 
  ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^2*gt^4*(Nc*Xb*Yb*((2*Log[mD]*Log[Q])/(mD - Q) + (2*Log[Q]^2)/(-mD + Q) + 
      Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) - 
      (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(mD - Q)) + 
      (2*A0*phi[A0, Q, Q])/((mD - Q)*Q)) + 
    Nc*(-1 + (4*Pi^2)/3 + 4*Log[A0]^2 + Log[Q]^2 - 2*Log[q2] - Log[q2]^2 + 
      Log[Q]*(2 + 2*Log[q2]) + Log[A0]*(-2*Log[mD] - 4*Log[Q] - 2*Log[U]) + 
      2*Log[mD]*Log[U] - (2*A0*phi[A0, Q, Q])/Q - 
      (2*(A0 + mD - U)*phi[A0, U, mD])/mD) + 
    (sbe^2*(Nc*Yb^2*(-Q^(-1) + ((A0^2 - (mD - Q)^2 + delta[A0, Q, mD])*
           Log[A0])/(2*Q*delta[A0, Q, mD]) + 
         ((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])*Log[mD])/
          (2*Q*delta[A0, Q, mD]) - ((A0 + mD - Q)*Log[Q])/delta[A0, Q, mD] - 
         Log[q2]/Q + ((-(A0 + mD - Q)^2 + delta[A0, Q, mD])*phi[A0, Q, mD])/
          (2*mD*delta[A0, Q, mD])) + Nc*(1/2 + Pi^2/3 - A0/Q + Log[A0]^2 + 
         Log[A0]*((A0 + Q)/Q - Log[mD] - Log[Q]) + Log[Q]*(1/2 + Log[q2]) + 
         Log[q2]*(-(A0/Q) - Log[U]) - (3*Log[U])/2 + Log[mD]*Log[U] - 
         ((A0 + mD - U)*phi[A0, U, mD])/mD) + 
       Nc*Xt*Yb*(Log[Q]*(4/(-Q + U) + (4*Log[q2])/(-Q + U)) + 
         (4*Log[U])/(Q - U) + (4*Log[q2]*Log[U])/(Q - U) + 
         Log[A0]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
         Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
         (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(-Q + U)) + 
         (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(Q - U))) + 
       Nc*Xt^3*Yb*(-16/(Q - U)^2 + Log[Q]*((4*(3*Q + U))/(Q - U)^3 + 
           (4*(Q + U)*Log[q2])/(Q - U)^3) - (4*(Q + 3*U)*Log[U])/(Q - U)^3 + 
         Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/(Q - U)^3) + 
         Log[A0]*((-2*(2*A0 - 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
           (2*(2*A0 - 2*mD + Q + U)*Log[U])/(Q - U)^3) + 
         Log[mD]*((-2*(-2*A0 + 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
           (2*(-2*A0 + 2*mD + Q + U)*Log[U])/(Q - U)^3) - 
         (2*((A0 + mD - Q)*(Q - U) + 2*delta[A0, Q, mD])*phi[A0, Q, mD])/
          (mD*(Q - U)^3) + ((2*(A0 + mD - U)*(-Q + U) + 4*delta[A0, U, mD])*
           phi[A0, U, mD])/(mD*(Q - U)^3)) + 
       Xt^4*(Nc*((Q*(Q - U) + A0*(5*Q + U))/(Q*(Q - U)^3) + 
           Log[Q]*((-(Q*(4*A0 + Q)) - 8*A0*U + U^2)/(2*(Q - U)^4) + 
             ((-(Q*(2*A0 + Q)) - 4*A0*U + U^2)*Log[q2])/(Q - U)^4) + 
           ((4*A0*Q + Q^2 + 8*A0*U - U^2)*Log[U])/(2*(Q - U)^4) + 
           Log[q2]*((2*Q*(Q - U) + A0*(5*Q + U))/(Q*(Q - U)^3) + 
             ((2*A0*Q + Q^2 + 4*A0*U - U^2)*Log[U])/(Q - U)^4) + 
           Log[A0]*((2*Q*(-Q + U) - A0*(5*Q + U))/(Q*(Q - U)^3) + 
             ((2*A0*Q + Q^2 + 4*A0*U - U^2)*Log[Q])/(Q - U)^4 + 
             ((-(Q*(2*A0 + Q)) - 4*A0*U + U^2)*Log[U])/(Q - U)^4)) + 
         Nc*Yb^2*((11*Q + U)/(Q*(Q - U)^3) + Log[Q]*
            ((-4*(2*Q + U))/(Q - U)^4 - (A0 + mD - Q)/((Q - U)^2*delta[A0, Q, 
                mD]) - (2*(Q + 2*U)*Log[q2])/(Q - U)^4) + 
           (2*(Q + 5*U)*Log[U])/(Q - U)^4 + Log[q2]*
            ((5*Q + U)/(Q*(Q - U)^3) + (2*(Q + 2*U)*Log[U])/(Q - U)^4) + 
           Log[A0]*((A0^2 - (mD - Q)^2 + delta[A0, Q, mD])/(2*Q*(Q - U)^
                2*delta[A0, Q, mD]) + ((3*A0 - 3*mD + Q + 2*U)*Log[Q])/
              (Q - U)^4 - ((3*A0 - 3*mD + Q + 2*U)*Log[U])/(Q - U)^4) + 
           Log[mD]*((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])/
              (2*Q*(Q - U)^2*delta[A0, Q, mD]) + ((-3*A0 + 3*mD + Q + 2*U)*
               Log[Q])/(Q - U)^4 - ((-3*A0 + 3*mD + Q + 2*U)*Log[U])/
              (Q - U)^4) + ((-((A0 + mD - Q)^2*(Q - U)^2) + (4*A0 + 4*mD - 
                3*Q - U)*(Q - U)*delta[A0, Q, mD] + 6*delta[A0, Q, mD]^2)*
             phi[A0, Q, mD])/(2*mD*(Q - U)^4*delta[A0, Q, mD]) + 
           (((A0 + mD - U)*(Q - U) - 3*delta[A0, U, mD])*phi[A0, U, mD])/
            (mD*(Q - U)^4))) + Xt^2*(Nc*Yb^2*(-2/(Q^2 - Q*U) + 
           Log[Q]*((2*((A0 + mD - Q)*(-Q + U) + delta[A0, Q, mD]))/
              ((Q - U)^2*delta[A0, Q, mD]) + (2*Log[q2])/(Q - U)^2) - 
           (2*Log[U])/(Q - U)^2 + Log[q2]*(-2/(Q^2 - Q*U) - 
             (2*Log[U])/(Q - U)^2) + Log[A0]*((A0^2 - (mD - Q)^2 + delta[A0, 
                Q, mD])/(Q*(Q - U)*delta[A0, Q, mD]) - Log[Q]/(Q - U)^2 + 
             Log[U]/(Q - U)^2) + Log[mD]*((-A0^2 + mD^2 + 2*A0*Q - Q^2 + 
               delta[A0, Q, mD])/(Q*(Q - U)*delta[A0, Q, mD]) - 
             Log[Q]/(Q - U)^2 + Log[U]/(Q - U)^2) + 
           ((-((A0 + mD - Q)^2*(Q - U)) + (A0 + mD - U)*delta[A0, Q, mD])*
             phi[A0, Q, mD])/(mD*(Q - U)^2*delta[A0, Q, mD]) - 
           ((A0 + mD - U)*phi[A0, U, mD])/(mD*(Q - U)^2)) + 
         Nc*((-2*A0 + 4*Q)/(Q*(Q - U)) + Log[Q]*
            (-((-2*A0 + 3*Q + U)/(Q - U)^2) + (2*(A0 - U)*Log[q2])/
              (Q - U)^2) + ((-2*A0 + Q + 3*U)*Log[U])/(Q - U)^2 + 
           Log[mD]*(((-A0 + mD + Q)*Log[Q])/(Q - U)^2 - 
             ((-A0 + mD + Q)*Log[U])/(Q - U)^2) + 
           Log[A0]*((2*A0)/(Q^2 - Q*U) - ((A0 + mD + Q - 2*U)*Log[Q])/
              (Q - U)^2 + ((A0 + mD + Q - 2*U)*Log[U])/(Q - U)^2) + 
           Log[q2]*((2*(-A0 + Q))/(Q*(Q - U)) - (2*(A0 - U)*Log[U])/
              (Q - U)^2) + (delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(Q - U)^2) + 
           (((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
            (mD*(Q - U)^2)))))/cbe^2 + 
    Nc*Xb*Yt*(Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) + 
      (2*Log[mD]*Log[U])/(mD - Q) + (2*Log[Q]*Log[U])/(-mD + Q) - 
      (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)) + 
      (2*(A0 + Q - U)*phi[A0, U, Q])/((mD - Q)*Q)) + 
    Xt*(Nc*Yb*(Log[Q]*(4/(-Q + U) + (4*Log[q2])/(-Q + U)) + 
        (4*Log[U])/(Q - U) + (4*Log[q2]*Log[U])/(Q - U) + 
        Log[A0]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
        Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
        (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(-Q + U)) + 
        (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(Q - U))) + 
      Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
        Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + (-Q + U)^(-1))*
           Log[Q] + (2*Log[U])/(Q - U)) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
          (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
        (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
        (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
        (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))) + 
      Yt*(Nc*((2*Log[Q]^2)/(Q - U) + (4*Log[U])/(Q - U) + 
          (4*Log[q2]*Log[U])/(Q - U) + Log[A0]*((2*Log[Q])/(Q - U) + 
            (2*Log[U])/(-Q + U)) + Log[Q]*(4/(-Q + U) + (4*Log[q2])/
             (-Q + U) + (2*Log[U])/(-Q + U)) - (2*A0*phi[A0, Q, Q])/
           (Q^2 - Q*U) + (2*(A0 + Q - U)*phi[A0, U, Q])/(Q*(Q - U))) + 
        Nc*Xb*Yb*((2*Log[Q]^2)/((mD - Q)*(-Q + U)) + (2*Log[Q]*Log[U])/
           ((mD - Q)*(Q - U)) + Log[mD]*((2*Log[Q])/((mD - Q)*(Q - U)) + 
            (2*Log[U])/((mD - Q)*(-Q + U))) - 
          (2*(A0 + mD - Q)*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) + 
          (2*A0*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) + 
          (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)*(Q - U)) + 
          (2*(A0 + Q - U)*phi[A0, U, Q])/(Q*(-mD + Q)*(Q - U))))) + 
    Xt^3*(Nc*Yb*(-16/(Q - U)^2 + Log[Q]*((4*(3*Q + U))/(Q - U)^3 + 
          (4*(Q + U)*Log[q2])/(Q - U)^3) - (4*(Q + 3*U)*Log[U])/(Q - U)^3 + 
        Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/(Q - U)^3) + 
        Log[A0]*((-2*(2*A0 - 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
          (2*(2*A0 - 2*mD + Q + U)*Log[U])/(Q - U)^3) + 
        Log[mD]*((-2*(-2*A0 + 2*mD + Q + U)*Log[Q])/(Q - U)^3 + 
          (2*(-2*A0 + 2*mD + Q + U)*Log[U])/(Q - U)^3) - 
        (2*((A0 + mD - Q)*(Q - U) + 2*delta[A0, Q, mD])*phi[A0, Q, mD])/
         (mD*(Q - U)^3) + ((2*(A0 + mD - U)*(-Q + U) + 4*delta[A0, U, mD])*
          phi[A0, U, mD])/(mD*(Q - U)^3)) + 
      Yt*(Nc*Xb*Yb*((2*(-2*A0 + 3*Q + U)*Log[Q]^2)/((mD - Q)*(Q - U)^3) - 
          (2*(-2*A0 + 3*Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)^3) + 
          Log[A0]*((4*Log[Q])/(Q - U)^3 + (4*Log[U])/(-Q + U)^3) + 
          Log[mD]*((-2*(-2*A0 + 2*mD + Q + U)*Log[Q])/((mD - Q)*(Q - U)^3) + 
            (2*(-2*A0 + 2*mD + Q + U)*Log[U])/((mD - Q)*(Q - U)^3)) - 
          (2*((A0 + mD - Q)*(Q - U) + 2*delta[A0, Q, mD])*phi[A0, Q, mD])/
           (mD*(mD - Q)*(Q - U)^3) + (2*(A0*(Q - U) + 2*delta[A0, Q, Q])*
            phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)^3) + 
          ((2*(A0 + mD - U)*(-Q + U) + 4*delta[A0, U, mD])*phi[A0, U, mD])/
           (mD*(mD - Q)*(Q - U)^3) + (2*((Q - U)*(A0 + Q - U) - 
             2*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)*Q*(Q - U)^3)) + 
        Nc*(-16/(Q - U)^2 - (2*(-2*A0 + 3*Q + U)*Log[Q]^2)/(Q - U)^3 - 
          (4*(Q + 3*U)*Log[U])/(Q - U)^3 + Log[A0]*
           ((2*(-2*A0 + Q - U)*Log[Q])/(Q - U)^3 + (2*(2*A0 - Q + U)*Log[U])/
             (Q - U)^3) + Log[q2]*(-8/(Q - U)^2 - (4*(Q + U)*Log[U])/
             (Q - U)^3) + Log[Q]*((4*(3*Q + U))/(Q - U)^3 + 
            (4*(Q + U)*Log[q2])/(Q - U)^3 + (2*(-2*A0 + 3*Q + U)*Log[U])/
             (Q - U)^3) - (2*(A0*(Q - U) + 2*delta[A0, Q, Q])*phi[A0, Q, Q])/
           (Q*(Q - U)^3) + ((-2*(Q - U)*(A0 + Q - U) + 4*delta[A0, U, Q])*
            phi[A0, U, Q])/(Q*(Q - U)^3)))) + 
    Nc*Xb^2*((mD - 5*Q)/(Q*(-mD + Q)) - (2*mD*Log[Q]^2)/(mD - Q)^2 + 
      ((mD - 3*Q)*Log[q2])/(Q*(-mD + Q)) + 
      Log[mD]*((mD*(mD - 5*Q))/((mD - Q)^2*Q) + (2*mD*Log[Q])/(mD - Q)^2 - 
        (2*mD*Log[q2])/(mD - Q)^2) + Log[Q]*((mD + 3*Q)/(mD - Q)^2 + 
        (2*mD*Log[q2])/(mD - Q)^2) + (6*PolyLog[2, 1 - mD/Q])/(-mD + Q)) + 
    (Nc*Xb*Xt^2*((8*Sqrt[mu2]*Q*sgn*Log[Q]^2)/((mu2 - Q)*(-mD + Q)*(Q - U)) + 
        (8*Sqrt[mu2]*Q*sgn*Log[Q]*Log[U])/((mu2 - Q)*(-mD + Q)*(-Q + U)) + 
        Log[mD]*((8*mD*Sqrt[mu2]*sgn*Log[Q])/((mD - mu2)*(mD - Q)*(-Q + U)) + 
          (8*mD*Sqrt[mu2]*sgn*Log[U])/((mD - mu2)*(mD - Q)*(Q - U))) + 
        Log[mu2]*((8*mu2^(3/2)*sgn*Log[Q])/((mD - mu2)*(mu2 - Q)*(Q - U)) + 
          (8*mu2^(3/2)*sgn*Log[U])/((-mD + mu2)*(mu2 - Q)*(Q - U)))) + 
      Nc*Xb*Xt^4*((-4*Sqrt[mu2]*Q*sgn*(Q + U)*Log[Q]^2)/((mu2 - Q)*(-mD + Q)*
          (Q - U)^3) + Log[mD]*((-8*mD*Sqrt[mu2]*sgn)/((mD - mu2)*(mD - Q)*
            (Q - U)^2) + (4*mD*Sqrt[mu2]*sgn*(Q + U)*Log[Q])/
           ((mD - mu2)*(mD - Q)*(Q - U)^3) - (4*mD*Sqrt[mu2]*sgn*(Q + U)*
            Log[U])/((mD - mu2)*(mD - Q)*(Q - U)^3)) + 
        Log[mu2]*((8*mu2^(3/2)*sgn)/((mD - mu2)*(mu2 - Q)*(Q - U)^2) - 
          (4*mu2^(3/2)*sgn*(Q + U)*Log[Q])/((mD - mu2)*(mu2 - Q)*(Q - U)^3) + 
          (4*mu2^(3/2)*sgn*(Q + U)*Log[U])/((mD - mu2)*(mu2 - Q)*
            (Q - U)^3)) + Log[Q]*((8*Sqrt[mu2]*Q*sgn)/((mu2 - Q)*(-mD + Q)*
            (Q - U)^2) + (4*Sqrt[mu2]*Q*sgn*(Q + U)*Log[U])/
           ((mu2 - Q)*(-mD + Q)*(Q - U)^3))) + 
      Nc*Xb*((4*Sqrt[mu2]*(mD + mu2)*sgn*Log[mD]^2)/((mD - mu2)*(mD - Q)) + 
        (4*mu2^(3/2)*sgn*Log[Q]^2)/((mD - Q)*(mu2 - Q)) + 
        (4*Sqrt[mu2]*Q*sgn*Log[Q]*Log[U])/((mu2 - Q)*(-mD + Q)) + 
        Log[mu2]*((-4*mu2^(3/2)*(mD - 2*mu2 + Q)*sgn*Log[Q])/
           ((mD - mu2)*(mD - Q)*(mu2 - Q)) + (4*mu2^(3/2)*sgn*Log[U])/
           ((mD - mu2)*(mu2 - Q))) + Log[mD]*((8*mu2^(3/2)*sgn*Log[mu2])/
           ((-mD + mu2)*(mD - Q)) + (4*mD*Sqrt[mu2]*sgn*Log[Q])/
           ((mD - mu2)*(-mD + Q)) + (4*mD*Sqrt[mu2]*sgn*Log[U])/
           ((mD - mu2)*(-mD + Q))) + (8*Sqrt[mu2]*(mD + mu2)*sgn*
          PolyLog[2, 1 - mu2/mD])/((mD - mu2)*(mD - Q)) - 
        (8*Sqrt[mu2]*(mu2 + Q)*sgn*PolyLog[2, 1 - mu2/Q])/
         ((mu2 - Q)*(-mD + Q))) + Nc*Xt*((4*Sqrt[mu2]*sgn*Log[Q]^2)/
         (-Q + U) + Log[Q]*((8*Sqrt[mu2]*sgn)/(Q - U) + 
          (8*Sqrt[mu2]*sgn*Log[q2])/(Q - U)) + (8*Sqrt[mu2]*sgn*Log[U])/
         (-Q + U) + (8*Sqrt[mu2]*sgn*Log[q2]*Log[U])/(-Q + U) + 
        (4*Sqrt[mu2]*sgn*Log[U]^2)/(Q - U) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/Q])/(-Q + U) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/U])/(Q - U)) + 
      Nc*Xt^3*((32*Sqrt[mu2]*sgn)/(Q - U)^2 + 
        (4*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*Log[Q]^2)/(Q - U)^3 + 
        Log[Q]*((-8*Sqrt[mu2]*sgn*(3*Q + U))/(Q - U)^3 - 
          (8*Sqrt[mu2]*sgn*(Q + U)*Log[q2])/(Q - U)^3) + 
        (8*Sqrt[mu2]*sgn*(Q + 3*U)*Log[U])/(Q - U)^3 - 
        (4*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*Log[U]^2)/(Q - U)^3 + 
        Log[mu2]*((16*mu2^(3/2)*sgn*Log[Q])/(Q - U)^3 + 
          (16*mu2^(3/2)*sgn*Log[U])/(-Q + U)^3) + 
        Log[q2]*((16*Sqrt[mu2]*sgn)/(Q - U)^2 + (8*Sqrt[mu2]*sgn*(Q + U)*
            Log[U])/(Q - U)^3) + (8*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*
          PolyLog[2, 1 - mu2/Q])/(Q - U)^3 - 
        (8*Sqrt[mu2]*sgn*(-2*mu2 + Q + U)*PolyLog[2, 1 - mu2/U])/(Q - U)^3))/
     (cbe*sbe) + (Nc*(-3/2 + mD/(mD - mu2) - (mD - 2*mu2)/Q + 
        ((mD^2 - 2*mD*mu2 - mu2^2)*Log[mD]^2)/(mD - mu2)^2 + Log[Q]^2 + 
        Log[Q]*(1/2 + mD/(-mD + mu2) - Log[q2]) + (3/2 + mD/(-mD + mu2))*
         Log[U] + Log[q2]*(-((mD - 2*mu2 + Q)/Q) + Log[U]) + 
        Log[mD]*(mD*((mD + 2*mu2)/(mD - mu2)^2 + Q^(-1)) + 
          (2*mu2^2*Log[mu2])/(mD - mu2)^2 - (mD*(mD - 2*mu2)*Log[Q])/
           (mD - mu2)^2 - (mD*(mD - 2*mu2)*Log[U])/(mD - mu2)^2) + 
        Log[mu2]*(-((mu2*(2*mD + mu2))/(mD - mu2)^2) - (2*mu2)/Q - 
          (mu2^2*Log[Q])/(mD - mu2)^2 - (mu2^2*Log[U])/(mD - mu2)^2) + 
        (2*(mD^2 - 2*mD*mu2 - mu2^2)*PolyLog[2, 1 - mu2/mD])/(mD - mu2)^2 + 
        2*PolyLog[2, 1 - mu2/Q]) + Nc*Xt^2*((2*(mD - 2*mu2 + 3*Q))/
         (Q*(-Q + U)) + (2*(mu2 - U)*Log[Q]^2)/(Q - U)^2 + 
        Log[Q]*((2*(mD - 2*mu2 + 3*Q) + ((mD + mu2)*(-Q + U))/(mD - mu2))/
           (Q - U)^2 + (2*(mD - 2*mu2 + Q + U)*Log[q2])/(Q - U)^2) + 
        ((-2*(mD - 2*mu2 + 3*Q) + ((5*mD - 3*mu2)*(Q - U))/(mD - mu2))*
          Log[U])/(Q - U)^2 - (2*(mu2 - U)*Log[U]^2)/(Q - U)^2 + 
        Log[mD]*((2*mD)/(Q^2 - Q*U) + (2*mD*(-1 + ((mD - 2*mu2)*(-Q + U))/
              (mD - mu2)^2)*Log[Q])/(Q - U)^2 + 
          (2*(mD + (mD*(mD - 2*mu2)*(Q - U))/(mD - mu2)^2)*Log[U])/
           (Q - U)^2) + Log[mu2]*((-4*mu2)/(Q^2 - Q*U) + 
          (2*mu2^2*Log[Q])/((mD - mu2)^2*(-Q + U)) + (2*mu2^2*Log[U])/
           ((mD - mu2)^2*(Q - U))) + Log[q2]*((2*(mD - 2*mu2 + 2*Q))/
           (Q*(-Q + U)) - (2*(mD - 2*mu2 + Q + U)*Log[U])/(Q - U)^2) + 
        (4*(mu2 - U)*PolyLog[2, 1 - mu2/Q])/(Q - U)^2 - 
        (4*(mu2 - U)*PolyLog[2, 1 - mu2/U])/(Q - U)^2) + 
      Nc*Xt^4*((mD^2*(5*Q + U) + 2*mu2*(-3*Q*(Q + U) + mu2*(8*Q + U)) + 
          mD*(-3*mu2*(7*Q + U) + 4*Q*(Q + 2*U)))/((mD - mu2)*Q*(Q - U)^3) + 
        ((mu2 - U)*(3*mu2 - 2*Q - U)*Log[Q]^2)/(Q - U)^4 + 
        Log[Q]*((Q*(-4*mD^2 + 3*mD*(8*mu2 - Q) + 5*mu2*(-4*mu2 + Q)) - 
            4*(mD - mu2)*(2*mD - 4*mu2 + 5*Q)*U - (mD + mu2)*U^2)/
           (2*(mD - mu2)*(Q - U)^4) - ((Q*(2*mD - 4*mu2 + Q) + 
             4*(mD - 2*mu2 + Q)*U + U^2)*Log[q2])/(Q - U)^4) + 
        ((4*mD^2*(Q + 2*U) + mu2*(-3*(Q + U)*(Q + 3*U) + 4*mu2*(2*Q + 7*U)) + 
           mD*(-12*mu2*(Q + 3*U) + (Q + U)*(Q + 11*U)))*Log[U])/
         (2*(mD - mu2)*(Q - U)^4) + ((3*mu2 - 2*Q - U)*(-mu2 + U)*Log[U]^2)/
         (Q - U)^4 + Log[mu2]*((2*mu2*(-3 + (((mD - mu2)^2 + mu2*Q)*(Q - U))/
              ((mD - mu2)^2*Q)))/(-Q + U)^3 + 
          (mu2^2*(-6*(mD - mu2)^2 + (Q - U)*(Q + U))*Log[Q])/
           ((mD - mu2)^2*(Q - U)^4) + (mu2^2*(6*(mD - mu2)^2 - Q^2 + U^2)*
            Log[U])/((mD - mu2)^2*(Q - U)^4)) + 
        Log[q2]*((3*Q*(Q + U) + mD*(5*Q + U) - 2*mu2*(5*Q + U))/
           (Q*(Q - U)^3) + ((Q*(2*mD - 4*mu2 + Q) + 4*(mD - 2*mu2 + Q)*U + 
             U^2)*Log[U])/(Q - U)^4) + Log[mD]*
         ((mD*(6 - ((mD^2 - 2*mD*(mu2 + Q) + mu2*(mu2 + 4*Q))*(Q - U))/
              ((mD - mu2)^2*Q)))/(-Q + U)^3 + 
          (mD*(Q*(2*(mD - mu2)^2 + (mD - 2*mu2)*Q) + 4*(mD - mu2)^2*U - 
             (mD - 2*mu2)*U^2)*Log[Q])/((mD - mu2)^2*(Q - U)^4) + 
          (mD*(Q*(-2*(mD - mu2)^2 - (mD - 2*mu2)*Q) - 4*(mD - mu2)^2*U + 
             (mD - 2*mu2)*U^2)*Log[U])/((mD - mu2)^2*(Q - U)^4)) + 
        (2*(mu2 - U)*(3*mu2 - 2*Q - U)*PolyLog[2, 1 - mu2/Q])/(Q - U)^4 - 
        (2*(-mu2 + U)*(-3*mu2 + 2*Q + U)*PolyLog[2, 1 - mu2/U])/(Q - U)^4))/
     cbe^2 + (cbe^2*(Nc*Xb^2*(4/(mD - Q) - ((-A0 + mD + Q)*Log[Q]^2)/
           (mD - Q)^2 + Log[A0]*(((A0 + mD - Q)*Log[mD])/(mD - Q)^2 - 
            ((A0 + mD - Q)*Log[Q])/(mD - Q)^2) + (2*Log[q2])/(mD - Q) + 
          Log[mD]*((-4*mD)/(mD - Q)^2 + ((-A0 + mD + Q)*Log[Q])/(mD - Q)^2 - 
            (2*mD*Log[q2])/(mD - Q)^2) + Log[Q]*((2*(mD + Q))/(mD - Q)^2 + 
            (2*mD*Log[q2])/(mD - Q)^2) + (delta[A0, Q, mD]*phi[A0, Q, mD])/
           (mD*(mD - Q)^2) + ((A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*
              delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^2*Q^2)) + 
        Nc*(Log[Q]^2 + Log[Q]*(-2 - 2*Log[q2]) + 2*Log[q2]^2 - 2*Log[U] + 
          Log[mD]*Log[U] - 2*Log[q2]*Log[U] + Log[A0]*(4 - Log[mD] + 
            Log[U]) + ((A0^2 - 5*A0*Q - delta[A0, Q, Q])*phi[A0, Q, Q])/Q^2 - 
          ((A0 + mD - U)*phi[A0, U, mD])/mD) + 
        Nc*Xb*Yt*(Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) + 
          (2*Log[mD]*Log[U])/(mD - Q) + (2*Log[Q]*Log[U])/(-mD + Q) - 
          (2*(A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)) + 
          (2*(A0 + Q - U)*phi[A0, U, Q])/((mD - Q)*Q)) + 
        Yt^2*(Nc*((((A0 + Q - U)^2 - delta[A0, U, Q])*Log[A0])/
             (Q*delta[A0, U, Q]) - (2*(A0 + Q - U)*Log[Q])/delta[A0, U, Q] + 
            ((-A0^2 + Q^2 + 2*A0*U - U^2 + delta[A0, U, Q])*Log[U])/
             (Q*delta[A0, U, Q]) + ((A0 - U + ((A0 + Q - U)*(-A0^2 + 
                  (Q - U)*U + A0*(Q + 2*U)))/delta[A0, U, Q])*phi[A0, U, Q])/
             Q^2) + Nc*Xb^2*(Log[A0]*((-(A0 + Q - U)^2 + delta[A0, U, Q])/(
                (mD - Q)*Q*delta[A0, U, Q]) - Log[mD]/(mD - Q)^2 + 
              Log[Q]/(mD - Q)^2) + ((-A0^2 + Q^2 + 2*A0*U - U^2 + delta[A0, 
                U, Q])*Log[U])/(Q*(-mD + Q)*delta[A0, U, Q]) + 
            (Log[mD]*Log[U])/(mD - Q)^2 + Log[Q]*((2*(A0 + Q - U))/((mD - Q)*
                delta[A0, U, Q]) - Log[U]/(mD - Q)^2) - 
            ((A0 + mD - U)*phi[A0, U, mD])/(mD*(mD - Q)^2) + 
            (((mD - Q)*(A0 + Q - U)*(A0^2 + U*(-Q + U) - A0*(Q + 2*U)) + 
               (-(A0*mD) + 2*A0*Q + Q^2 + mD*U - 2*Q*U)*delta[A0, U, Q])*
              phi[A0, U, Q])/((mD - Q)^2*Q^2*delta[A0, U, Q]))) + 
        Xt*(Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
            (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
            Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + 
                (-Q + U)^(-1))*Log[Q] + (2*Log[U])/(Q - U)) + 
            Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
              (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
            (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
            (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
            (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
            (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))) + 
          Yt*(Nc*((2*Log[Q]^2)/(Q - U) + (2*Log[Q]*Log[U])/(-Q + U) + 
              Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) - 
              (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/(Q^2*
                (Q - U)) - (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2 - 
                 delta[A0, U, Q])*phi[A0, U, Q])/(Q^2*(Q - U))) + 
            Nc*Xb^2*((2*(A0 - mD - Q)*Log[Q]^2)/((mD - Q)^2*(Q - U)) + 
              Log[A0]*((-2*Log[mD])/(mD - Q)^2 + (2*Log[Q])/(mD - Q)^2) + 
              (2*(-A0 + mD + U)*Log[Q]*Log[U])/((mD - Q)^2*(Q - U)) + 
              Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)^2*(Q - U)) - 
                (2*(-A0 + mD + U)*Log[U])/((mD - Q)^2*(Q - U))) + 
              (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(mD - Q)^2*(Q - U)) + 
              (2*(A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*delta[A0, Q, Q])*
                phi[A0, Q, Q])/((mD - Q)^2*Q^2*(Q - U)) - (2*delta[A0, U, mD]*
                phi[A0, U, mD])/(mD*(mD - Q)^2*(Q - U)) + 
              (2*((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2) - 
                 (mD - 2*Q)*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)^2*Q^2*
                (Q - U))))) + 
        Xt^2*(Nc*(4/(Q - U) + Log[Q]*((-4*Q)/(Q - U)^2 - (2*Q*Log[q2])/
               (Q - U)^2) + (2*(Q + U)*Log[U])/(Q - U)^2 + 
            Log[q2]*(2/(Q - U) + (2*Q*Log[U])/(Q - U)^2) + 
            Log[A0]*(((A0 - mD + Q)*Log[Q])/(Q - U)^2 - ((A0 - mD + Q)*
                Log[U])/(Q - U)^2) + Log[mD]*(((-A0 + mD + Q)*Log[Q])/(Q - U)^
                2 - ((-A0 + mD + Q)*Log[U])/(Q - U)^2) + 
            (delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(Q - U)^2) + 
            (((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
             (mD*(Q - U)^2)) + Nc*Xb*Yt*((2*(A0 - 2*Q)*Log[Q]^2)/
             ((mD - Q)*(Q - U)^2) + (2*(A0 - 2*Q)*Log[Q]*Log[U])/
             ((-mD + Q)*(Q - U)^2) + Log[A0]*((-2*Log[Q])/(Q - U)^2 + 
              (2*Log[U])/(Q - U)^2) + Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/(
                (mD - Q)*(Q - U)^2) + (2*(A0 - mD - Q)*Log[U])/((mD - Q)*
                (Q - U)^2)) + (2*delta[A0, Q, mD]*phi[A0, Q, mD])/
             (mD*(mD - Q)*(Q - U)^2) - (2*delta[A0, Q, Q]*phi[A0, Q, Q])/
             ((mD - Q)*Q*(Q - U)^2) + (2*((A0 + mD - U)*(Q - U) - delta[A0, 
                U, mD])*phi[A0, U, mD])/(mD*(mD - Q)*(Q - U)^2) + 
            ((-2*(Q - U)*(A0 + Q - U) + 2*delta[A0, U, Q])*phi[A0, U, Q])/
             ((mD - Q)*Q*(Q - U)^2)) + 
          Yt^2*(Nc*(Log[Q]^2/(Q - U)^2 + ((A0^2 - Q^2 - 2*A0*U + U^2 - 
                 delta[A0, U, Q])*Log[U])/(Q*(Q - U)*delta[A0, U, Q]) + 
              Log[Q]*((2*(A0 + Q - U))/((Q - U)*delta[A0, U, Q]) - 
                Log[U]/(Q - U)^2) + Log[A0]*((-(A0 + Q - U)^2 + delta[A0, U, 
                   Q])/(Q*(Q - U)*delta[A0, U, Q]) - Log[Q]/(Q - U)^2 + 
                Log[U]/(Q - U)^2) + ((A0*(A0 - 5*Q) - delta[A0, Q, Q])*
                phi[A0, Q, Q])/(Q^2*(Q - U)^2) + (((Q - U)*(A0 + Q - U)*
                  (A0^2 + U*(-Q + U) - A0*(Q + 2*U)) + delta[A0, U, Q]*
                  (-A0^2 - 2*(Q - U)^2 + A0*(2*Q + 3*U) + delta[A0, U, Q]))*
                phi[A0, U, Q])/(Q^2*(Q - U)^2*delta[A0, U, Q])) + 
            Nc*Xb^2*(((-(A0 + Q - U)^2 + delta[A0, U, Q])*Log[A0])/(Q*
                (-mD + Q)*(Q - U)*delta[A0, U, Q]) - ((-A0 + mD + Q)*
                Log[Q]^2)/((mD - Q)^2*(Q - U)^2) + ((A0^2 - Q^2 - 2*A0*U + 
                 U^2 - delta[A0, U, Q])*Log[U])/(Q*(-mD + Q)*(Q - U)*
                delta[A0, U, Q]) + Log[mD]*(((-A0 + mD + Q)*Log[Q])/
                 ((mD - Q)^2*(Q - U)^2) - ((-A0 + mD + Q)*Log[U])/
                 ((mD - Q)^2*(Q - U)^2)) + Log[Q]*((2*(A0 + Q - U))/
                 ((-mD + Q)*(Q - U)*delta[A0, U, Q]) + ((-A0 + mD + Q)*
                  Log[U])/((mD - Q)^2*(Q - U)^2)) + (delta[A0, Q, mD]*
                phi[A0, Q, mD])/(mD*(mD - Q)^2*(Q - U)^2) + 
              ((A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*delta[A0, Q, Q])*
                phi[A0, Q, Q])/((mD - Q)^2*Q^2*(Q - U)^2) + 
              (((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/(mD*
                (mD - Q)^2*(Q - U)^2) + (((mD - Q)*(Q - U)*(A0 + Q - U)*
                  (-A0^2 + (Q - U)*U + A0*(Q + 2*U)) + delta[A0, U, Q]*
                  (A0^2*(mD - Q) + (2*mD - 3*Q)*(Q - U)^2 + A0*(-2*mD*Q + 
                     Q^2 - 3*mD*U + 4*Q*U) - (mD - 2*Q)*delta[A0, U, Q]))*
                phi[A0, U, Q])/((mD - Q)^2*Q^2*(Q - U)^2*delta[A0, U, 
                 Q]))))) + Nc*Xb^2*(4/(-mD + Q) - (2*mD*Log[mD]^2)/
         (mD - Q)^2 + (4*mu2*Log[mu2])/((mD - Q)*(mu2 - Q)) + 
        (2*mD*Log[Q]^2)/(mD - Q)^2 + (2*Log[q2])/(-mD + Q) + 
        (2*Log[U])/(-mD + Q) + Log[Q]*((4*(mD*mu2 - Q^2))/
           ((mD - Q)^2*(-mu2 + Q)) - (2*mD*Log[q2])/(mD - Q)^2 - 
          (2*mD*Log[U])/(mD - Q)^2) + Log[mD]*((4*mD)/(mD - Q)^2 + 
          (2*mD*Log[q2])/(mD - Q)^2 + (2*mD*Log[U])/(mD - Q)^2) - 
        (4*mD*PolyLog[2, 1 - mu2/mD])/(mD - Q)^2 + 
        (4*mD*PolyLog[2, 1 - mu2/Q])/(mD - Q)^2) + 
      Nc*(4*mu2*((-mu2 + Q)^(-1) + (-mu2 + U)^(-1))*Log[mu2] - 2*Log[Q]^2 - 
        2*Log[q2]^2 + (4*mu2*Log[U])/(mu2 - U) + 2*Log[q2]*Log[U] - 
        2*Log[U]^2 + Log[Q]*((4*mu2)/(mu2 - Q) + 2*Log[q2] + 2*Log[U]) - 
        4*PolyLog[2, 1 - mu2/Q] - 4*PolyLog[2, 1 - mu2/U]) + 
      Xt^2*(Nc*Xb^2*(2/((mD - Q)*(Q - U)) - (2*mD*Q*Log[Q]^2)/
           ((mD - Q)^2*(Q - U)^2) + (2*Q*Log[U])/((mD - Q)*(Q - U)^2) + 
          Log[mD]*((2*mD)/((mD - Q)^2*(-Q + U)) + (2*mD*Q*Log[Q])/
             ((mD - Q)^2*(Q - U)^2) - (2*mD*Q*Log[U])/((mD - Q)^2*
              (Q - U)^2)) + Log[Q]*((2*(Q^2 - mD*U))/((mD - Q)^2*(Q - U)^2) + 
            (2*mD*Q*Log[U])/((mD - Q)^2*(Q - U)^2))) + 
        Nc*(4/(-Q + U) + (4*mu2*Log[mu2])/((mu2 - U)*(Q - U)) + 
          ((-4*mu2*Q + 4*U^2)*Log[U])/((mu2 - U)*(Q - U)^2) + 
          (2*Q*Log[U]^2)/(Q - U)^2 + Log[q2]*(2/(-Q + U) - 
            (2*Q*Log[U])/(Q - U)^2) + Log[Q]*((2*(Q + U))/(Q - U)^2 + 
            (2*Q*Log[q2])/(Q - U)^2 - (2*Q*Log[U])/(Q - U)^2) - 
          (4*Q*PolyLog[2, 1 - mu2/Q])/(Q - U)^2 + (4*Q*PolyLog[2, 1 - mu2/U])/
           (Q - U)^2)))/sbe^2 + 
    Xt^4*(Nc*(24/(Q - U)^2 + ((7*Q + U)*Log[Q]^2)/(Q - U)^3 + 
        (6*(Q + 3*U)*Log[U])/(Q - U)^3 + (3*Log[U]^2)/(Q - U)^2 + 
        Log[Q]*((-2*(7*Q + 5*U))/(Q - U)^3 - (8*(Q + U)*Log[q2])/(Q - U)^3 + 
          (2*(-5*Q + U)*Log[U])/(Q - U)^3) + Log[A0]*(-8/(Q - U)^2 + 
          (4*(Q + U)*Log[Q])/(Q - U)^3 - (4*(Q + U)*Log[U])/(Q - U)^3) + 
        Log[q2]*(16/(Q - U)^2 + (8*(Q + U)*Log[U])/(Q - U)^3) + 
        (6*PolyLog[2, 1 - Q/U])/(Q - U)^2) + 
      Nc*Xb^2*((11*Q + U)/(Q*(Q - U)^3) + 
        ((-3*(mD + Q - 2*U) + (2*(Q - U)*(4*Q^2 + mD*(-5*Q + U)))/(mD - Q)^2)*
          Log[Q]^2)/(Q - U)^4 + (2*(U*(-7*Q + U) + mD*(Q + 5*U))*Log[U])/
         ((mD - Q)*(Q - U)^4) + (3*(mD + Q - 2*U)*Log[U]^2)/(Q - U)^4 + 
        Log[mD]*(mD/((mD - Q)*Q*(Q - U)^2) + (2*mD*(mD - U)*(3*mD - 2*Q - U)*
            Log[Q])/((mD - Q)^2*(Q - U)^4) - (2*mD*(mD - U)*(3*mD - 2*Q - U)*
            Log[U])/((mD - Q)^2*(Q - U)^4)) + 
        Log[q2]*((5*Q + U)/(Q*(Q - U)^3) + (2*(Q + 2*U)*Log[U])/(Q - U)^4) + 
        Log[Q]*((5*Q^2 + 10*Q*U - 3*U^2 - 4*mD*(2*Q + U))/
           ((mD - Q)*(Q - U)^4) - (2*(Q + 2*U)*Log[q2])/(Q - U)^4 + 
          (2*(5*mD*Q - 4*Q^2 - mD*U)*Log[U])/((mD - Q)^2*(Q - U)^3)) - 
        (6*(mD - U)^2*PolyLog[2, 1 - mD/Q])/((mD - Q)*(Q - U)^4) + 
        (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/((mD - Q)*(Q - U)^4) - 
        (6*PolyLog[2, 1 - Q/U])/((mD - Q)*(Q - U)^2))) + 
    Xt^2*(Nc*Xb*Yb*((2*(-A0 + Q + U)*Log[Q]^2)/((mD - Q)*(Q - U)^2) - 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)^2) + 
        Log[A0]*((2*Log[Q])/(Q - U)^2 - (2*Log[U])/(Q - U)^2) + 
        Log[mD]*((-2*(-A0 + mD + U)*Log[Q])/((mD - Q)*(Q - U)^2) + 
          (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U)^2)) - 
        (2*((A0 + mD - Q)*(Q - U) + delta[A0, Q, mD])*phi[A0, Q, mD])/
         (mD*(mD - Q)*(Q - U)^2) + (2*(A0*(Q - U) + delta[A0, Q, Q])*
          phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)^2) + 
        (2*delta[A0, U, mD]*phi[A0, U, mD])/(mD*(mD - Q)*(Q - U)^2) - 
        (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U)^2)) + 
      Nc*Xb*Yt*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)^2) + 
        (2*(A0 - 2*Q)*Log[Q]*Log[U])/((-mD + Q)*(Q - U)^2) + 
        Log[A0]*((-2*Log[Q])/(Q - U)^2 + (2*Log[U])/(Q - U)^2) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)^2) + 
          (2*(A0 - mD - Q)*Log[U])/((mD - Q)*(Q - U)^2)) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(mD - Q)*(Q - U)^2) - 
        (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)^2) + 
        (2*((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
         (mD*(mD - Q)*(Q - U)^2) + ((-2*(Q - U)*(A0 + Q - U) + 
           2*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)*Q*(Q - U)^2)) + 
      Nc*(8/(Q - U) + (2*A0*Log[Q]^2)/(Q - U)^2 + ((-6*Q + 14*U)*Log[U])/
         (Q - U)^2 + (6*Log[U]^2)/(Q - U) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/(Q - U)^2 + 
          (2*(A0 - mD - Q)*Log[U])/(Q - U)^2) + 
        Log[A0]*((-2*(mD + 2*Q - 3*U)*Log[Q])/(Q - U)^2 + 
          (2*(mD + 2*Q - 3*U)*Log[U])/(Q - U)^2) + 
        Log[Q]*((2*(Q - 5*U))/(Q - U)^2 + (4*(2*Q - 3*U)*Log[q2])/(Q - U)^2 - 
          (2*(A0 + 3*Q - 3*U)*Log[U])/(Q - U)^2) + 
        Log[q2]*(4/(Q - U) + ((-8*Q + 12*U)*Log[U])/(Q - U)^2) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(Q - U)^2) - 
        (2*(A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/(Q*(Q - U)^2) + 
        (2*((A0 + mD - U)*(Q - U) - delta[A0, U, mD])*phi[A0, U, mD])/
         (mD*(Q - U)^2) + (2*delta[A0, U, Q]*phi[A0, U, Q])/(Q*(Q - U)^2) + 
        (12*PolyLog[2, 1 - Q/U])/(Q - U)) + 
      Nc*Xb^2*(-2/(Q^2 - Q*U) + ((-3*mD^2 - 2*mD*Q + Q^2 + 4*mD*U)*Log[Q]^2)/
         ((mD - Q)^2*(Q - U)^2) - (2*(mD - U)*Log[U])/((mD - Q)*(Q - U)^2) + 
        (3*Log[U]^2)/(Q - U)^2 + Log[q2]*(-2/(Q^2 - Q*U) - 
          (2*Log[U])/(Q - U)^2) + Log[mD]*((2*mD)/((mD - Q)*Q*(Q - U)) + 
          (4*mD*(mD - U)*Log[Q])/((mD - Q)^2*(Q - U)^2) + 
          (4*mD*(-mD + U)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
        Log[Q]*(2/(Q - U)^2 + (2*Log[q2])/(Q - U)^2 + 
          ((8*mD*Q - 4*Q^2 - 4*mD*U)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
        (6*(mD + Q - 2*U)*PolyLog[2, 1 - mD/Q])/((-mD + Q)*(Q - U)^2) + 
        (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/((mD - Q)^2*(Q - U)^2) + 
        (6*(-2*mD + Q + U)*PolyLog[2, 1 - Q/U])/((mD - Q)^2*(Q - U))))) + 
  ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*gt^2*(Nc*Xb*Yb*((2*Log[Q]^2)/(-mD + Q) + 
      Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
      Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q)) + 
      Log[mD]*(4/(-mD + Q) + (2*Log[Q])/(mD - Q) + (4*Log[q2])/(-mD + Q)) + 
      (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - delta[A0, Q, mD])*
        phi[A0, Q, mD])/(mD^2*(mD - Q)) + 
      (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
       ((mD - Q)*Q^2)) + Nc*Xb^3*Yb*(-16/(mD - Q)^2 + 
      (2*(-2*A0 + mD + 3*Q)*Log[Q]^2)/(mD - Q)^3 + 
      Log[A0]*((-2*(2*A0 + mD - Q)*Log[mD])/(mD - Q)^3 + 
        (2*(2*A0 + mD - Q)*Log[Q])/(mD - Q)^3) - (8*Log[q2])/(mD - Q)^2 + 
      Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - (4*(mD + Q)*Log[q2])/(mD - Q)^3) + 
      Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 - (2*(-2*A0 + mD + 3*Q)*Log[Q])/
         (mD - Q)^3 + (4*(mD + Q)*Log[q2])/(mD - Q)^3) + 
      (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
         (-3*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)^3) + 
      (2*(A0*(A0 - 5*Q)*(mD - Q) - (mD - 3*Q)*delta[A0, Q, Q])*phi[A0, Q, Q])/
       ((mD - Q)^3*Q^2)) + Nc*(-1 + (4*Pi^2)/3 + 4*Log[A0]^2 + Log[Q]^2 + 
      Log[Q]*(-2 - 2*Log[q2]) + 6*Log[q2] + 7*Log[q2]^2 + 
      Log[A0]*(2*Log[mD] - 8*Log[q2] - 2*Log[U]) + 
      Log[mD]*(-4 - 4*Log[q2] + 2*Log[U]) - 
      (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/Q^2 + 
      (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - delta[A0, U, mD])*
        phi[A0, U, mD])/mD^2) + 
    Yt*(Nc*Xb*(Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
        Log[mD]*(4/(-mD + Q) + (4*Log[q2])/(-mD + Q) + (2*Log[U])/(mD - Q)) + 
        Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q) + (2*Log[U])/(-mD + Q)) + 
        (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
           delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(mD - Q)) - 
        (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2 - delta[A0, U, Q])*
          phi[A0, U, Q])/((mD - Q)*Q^2)) + 
      Nc*Xb^3*(-16/(mD - Q)^2 + Log[A0]*((-2*(2*A0 + mD + Q - 2*U)*Log[mD])/
           (mD - Q)^3 + (2*(2*A0 + mD + Q - 2*U)*Log[Q])/(mD - Q)^3) - 
        (8*Log[q2])/(mD - Q)^2 + Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 + 
          (4*(mD + Q)*Log[q2])/(mD - Q)^3 - (2*(-2*A0 + mD + Q + 2*U)*Log[U])/
           (mD - Q)^3) + Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - 
          (4*(mD + Q)*Log[q2])/(mD - Q)^3 + (2*(-2*A0 + mD + Q + 2*U)*Log[U])/
           (mD - Q)^3) + (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
             3*mD*U + U^2) + (-3*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
         (mD^2*(mD - Q)^3) + (2*((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
             3*Q*U + U^2) - (mD - 3*Q)*delta[A0, U, Q])*phi[A0, U, Q])/
         ((mD - Q)^3*Q^2))) + 
    (sbe^2*(Nc*Yb^2*((((A0 + mD - Q)^2 - delta[A0, Q, mD])*Log[A0])/
          (mD*delta[A0, Q, mD]) - (2*(A0 + mD - Q)*Log[mD])/
          delta[A0, Q, mD] + ((-A0^2 + mD^2 + 2*A0*Q - Q^2 + 
            delta[A0, Q, mD])*Log[Q])/(mD*delta[A0, Q, mD]) + 
         ((A0 - Q + ((A0 + mD - Q)*(-A0^2 + (mD - Q)*Q + A0*(mD + 2*Q)))/
             delta[A0, Q, mD])*phi[A0, Q, mD])/mD^2) + 
       Nc*Xb*Yb*((2*Log[mD]*Log[Q])/(mD - Q) + (2*Log[Q]^2)/(-mD + Q) + 
         Log[A0]*((2*Log[mD])/(-mD + Q) + (2*Log[Q])/(mD - Q)) - 
         (2*(A0 + mD - Q)*phi[A0, Q, mD])/(mD*(mD - Q)) + 
         (2*A0*phi[A0, Q, Q])/((mD - Q)*Q)) + 
       Nc*(Log[Q]^2 + Log[Q]*(-2 - 2*Log[q2]) + 2*Log[q2]^2 + 
         Log[A0]*(4 + Log[mD] - Log[U]) + Log[mD]*(-2 - 2*Log[q2] + Log[U]) - 
         (A0*phi[A0, Q, Q])/Q + ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + 
            U^2 - delta[A0, U, mD])*phi[A0, U, mD])/mD^2) + 
       Xb^2*(Nc*Yb^2*(((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])*
             Log[Q])/(mD*(mD - Q)*delta[A0, Q, mD]) + Log[Q]^2/(mD - Q)^2 + 
           Log[mD]*((-2*(A0 + mD - Q))/((mD - Q)*delta[A0, Q, mD]) - 
             Log[Q]/(mD - Q)^2) + Log[A0]*(((A0 + mD - Q)^2 - delta[A0, Q, 
                mD])/(mD*(mD - Q)*delta[A0, Q, mD]) + Log[mD]/(mD - Q)^2 - 
             Log[Q]/(mD - Q)^2) + (((mD - Q)*(A0 + mD - Q)*(-A0^2 + 
                (mD - Q)*Q + A0*(mD + 2*Q)) + ((mD - Q)^2 + A0*(2*mD - Q))*
               delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)^2*
             delta[A0, Q, mD]) - (A0*phi[A0, Q, Q])/((mD - Q)^2*Q)) + 
         Nc*(4/(-mD + Q) + Log[A0]*(-(((A0 + Q - U)*Log[mD])/(mD - Q)^2) + 
             ((A0 + Q - U)*Log[Q])/(mD - Q)^2) + (2*Log[q2])/(-mD + Q) + 
           Log[mD]*((2*(mD + Q))/(mD - Q)^2 + (2*Q*Log[q2])/(mD - Q)^2 - 
             ((-A0 + Q + U)*Log[U])/(mD - Q)^2) + Log[Q]*((-4*Q)/(mD - Q)^2 - 
             (2*Q*Log[q2])/(mD - Q)^2 + ((-A0 + Q + U)*Log[U])/(mD - Q)^2) + 
           (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
              (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
            (mD^2*(mD - Q)^2) + (delta[A0, U, Q]*phi[A0, U, Q])/
            ((mD - Q)^2*Q))) + 
       Xt^2*(Nc*Yb^2*(((-A0^2 + mD^2 + 2*A0*Q - Q^2 + delta[A0, Q, mD])*
             Log[Q])/(mD*(Q - U)*delta[A0, Q, mD]) + 
           Log[A0]*(((A0 + mD - Q)^2 - delta[A0, Q, mD])/(mD*(Q - U)*delta[
                A0, Q, mD]) + Log[Q]/(Q - U)^2 - Log[U]/(Q - U)^2) + 
           Log[mD]*((-2*(A0 + mD - Q))/((Q - U)*delta[A0, Q, mD]) - 
             Log[Q]/(Q - U)^2 + Log[U]/(Q - U)^2) + 
           ((-((A0 + mD - Q)*(A0^2 + Q*(-mD + Q) - A0*(mD + 2*Q))*(Q - U)) + 
              delta[A0, Q, mD]*(-A0^2 - 2*mD^2 + 3*mD*Q + 3*A0*(mD + Q) - 
                A0*U + Q*(-2*Q + U) + delta[A0, Q, mD]))*phi[A0, Q, mD])/
            (mD^2*(Q - U)^2*delta[A0, Q, mD]) + 
           ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
              delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(Q - U)^2)) + 
         Nc*(4/(-Q + U) - ((-A0 + Q + U)*Log[Q]^2)/(Q - U)^2 - 
           (4*U*Log[U])/(Q - U)^2 + Log[q2]*(2/(-Q + U) - (2*U*Log[U])/
              (Q - U)^2) + Log[A0]*(-(((A0 - Q + U)*Log[Q])/(Q - U)^2) + 
             ((A0 - Q + U)*Log[U])/(Q - U)^2) + 
           Log[Q]*((2*(Q + U))/(Q - U)^2 + (2*U*Log[q2])/(Q - U)^2 + 
             ((-A0 + Q + U)*Log[U])/(Q - U)^2) - 
           ((A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/(Q*(Q - U)^2) + 
           (delta[A0, U, Q]*phi[A0, U, Q])/(Q*(Q - U)^2)) + 
         Nc*Xb^2*Yb^2*((((A0 + mD - Q)^2 - delta[A0, Q, mD])*Log[A0])/
            (mD*(mD - Q)*(Q - U)*delta[A0, Q, mD]) - ((-A0 + Q + U)*Log[Q]^2)/
            ((mD - Q)^2*(Q - U)^2) + Log[mD]*((-2*(A0 + mD - Q))/
              ((mD - Q)*(Q - U)*delta[A0, Q, mD]) + ((-A0 + Q + U)*Log[Q])/
              ((mD - Q)^2*(Q - U)^2) - ((-A0 + Q + U)*Log[U])/
              ((mD - Q)^2*(Q - U)^2)) + Log[Q]*((-A0^2 + mD^2 + 2*A0*Q - Q^
                2 + delta[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)*delta[A0, Q, 
                mD]) + ((-A0 + Q + U)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
           (((mD - Q)*(A0 + mD - Q)*(-A0^2 + (mD - Q)*Q + A0*(mD + 2*Q))*(Q - 
                U) + delta[A0, Q, mD]*(A0^2*(-mD + Q) - (mD - Q)^2*
                 (2*mD - 2*Q + U) + A0*(3*mD^2 + mD*(Q - 2*U) + 
                  Q*(-3*Q + U)) + (2*mD - Q)*delta[A0, Q, mD]))*
             phi[A0, Q, mD])/(mD^2*(mD - Q)^2*(Q - U)^2*delta[A0, Q, mD]) - 
           ((A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^2*Q*
             (Q - U)^2) + (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
                3*mD*U + U^2) + (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
            (mD^2*(mD - Q)^2*(Q - U)^2) + (delta[A0, U, Q]*phi[A0, U, Q])/
            ((mD - Q)^2*Q*(Q - U)^2)) + Nc*Xb*Yb*((2*(-A0 + Q + U)*Log[Q]^2)/
            ((mD - Q)*(Q - U)^2) - (2*(-A0 + Q + U)*Log[Q]*Log[U])/
            ((mD - Q)*(Q - U)^2) + Log[A0]*((2*Log[Q])/(Q - U)^2 - 
             (2*Log[U])/(Q - U)^2) + Log[mD]*((-2*(-A0 + mD + U)*Log[Q])/
              ((mD - Q)*(Q - U)^2) + (2*(-A0 + mD + U)*Log[U])/
              ((mD - Q)*(Q - U)^2)) - (2*((A0 + mD - Q)*(Q - U) + 
              delta[A0, Q, mD])*phi[A0, Q, mD])/(mD*(mD - Q)*(Q - U)^2) + 
           (2*(A0*(Q - U) + delta[A0, Q, Q])*phi[A0, Q, Q])/
            ((mD - Q)*Q*(Q - U)^2) + (2*delta[A0, U, mD]*phi[A0, U, mD])/
            (mD*(mD - Q)*(Q - U)^2) - (2*delta[A0, U, Q]*phi[A0, U, Q])/
            ((mD - Q)*Q*(Q - U)^2))) + 
       Xt*(Nc*Yb*(Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) + 
           Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
           (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
              delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(Q - U)) - 
           (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
              delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(Q - U))) + 
         Nc*Xb^2*Yb*(((-2*A0 + 4*Q)*Log[Q]^2)/((mD - Q)^2*(Q - U)) + 
           Log[A0]*((2*Log[mD])/(mD - Q)^2 - (2*Log[Q])/(mD - Q)^2) - 
           (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)^2*(Q - U)) + 
           Log[mD]*((2*(A0 - 2*Q)*Log[Q])/((mD - Q)^2*(Q - U)) + 
             (2*(-A0 + Q + U)*Log[U])/((mD - Q)^2*(Q - U))) + 
           (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
              (-2*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/
            (mD^2*(mD - Q)^2*(Q - U)) + (2*delta[A0, Q, Q]*phi[A0, Q, Q])/
            ((mD - Q)^2*Q*(Q - U)) - (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 
                2*A0*U - 3*mD*U + U^2) + (-2*mD + Q)*delta[A0, U, mD])*
             phi[A0, U, mD])/(mD^2*(mD - Q)^2*(Q - U)) - 
           (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)^2*Q*(Q - U))) + 
         Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
           (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
           Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + (-Q + U)^(-1))*
              Log[Q] + (2*Log[U])/(Q - U)) + Log[mD]*
            ((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
             (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
           (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
           (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
           (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
           (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))))))/cbe^2 + 
    Xt*(Nc*Yb*(Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) + 
        Log[mD]*((2*Log[Q])/(Q - U) + (2*Log[U])/(-Q + U)) + 
        (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
           delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(Q - U)) - 
        (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
           delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(Q - U))) + 
      Nc*Xb^2*Yb*(((-2*A0 + 4*Q)*Log[Q]^2)/((mD - Q)^2*(Q - U)) + 
        Log[A0]*((2*Log[mD])/(mD - Q)^2 - (2*Log[Q])/(mD - Q)^2) - 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)^2*(Q - U)) + 
        Log[mD]*((2*(A0 - 2*Q)*Log[Q])/((mD - Q)^2*(Q - U)) + 
          (2*(-A0 + Q + U)*Log[U])/((mD - Q)^2*(Q - U))) + 
        (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
           (-2*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/
         (mD^2*(mD - Q)^2*(Q - U)) + (2*delta[A0, Q, Q]*phi[A0, Q, Q])/
         ((mD - Q)^2*Q*(Q - U)) - 
        (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
           (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
         (mD^2*(mD - Q)^2*(Q - U)) - (2*delta[A0, U, Q]*phi[A0, U, Q])/
         ((mD - Q)^2*Q*(Q - U))) + 
      Nc*Xb*((2*(A0 - 2*Q)*Log[Q]^2)/((mD - Q)*(Q - U)) + 
        (2*(-A0 + Q + U)*Log[Q]*Log[U])/((mD - Q)*(Q - U)) + 
        Log[A0]*((2*Log[mD])/(-mD + Q) + 2*((mD - Q)^(-1) + (-Q + U)^(-1))*
           Log[Q] + (2*Log[U])/(Q - U)) + 
        Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/((mD - Q)*(Q - U)) - 
          (2*(-A0 + mD + U)*Log[U])/((mD - Q)*(Q - U))) + 
        (2*delta[A0, Q, mD]*phi[A0, Q, mD])/((mD^2 - mD*Q)*(Q - U)) - 
        (2*delta[A0, Q, Q]*phi[A0, Q, Q])/((mD - Q)*Q*(Q - U)) - 
        (2*delta[A0, U, mD]*phi[A0, U, mD])/((mD^2 - mD*Q)*(Q - U)) + 
        (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)*Q*(Q - U))) + 
      Yt*(Nc*((2*Log[Q]^2)/(Q - U) + (2*Log[Q]*Log[U])/(-Q + U) + 
          Log[A0]*((2*Log[Q])/(-Q + U) + (2*Log[U])/(Q - U)) - 
          (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
           (Q^2*(Q - U)) - (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2 - 
             delta[A0, U, Q])*phi[A0, U, Q])/(Q^2*(Q - U))) + 
        Nc*Xb*Yb*((2*Log[Q]^2)/((mD - Q)*(-Q + U)) + (2*Log[Q]*Log[U])/
           ((mD - Q)*(Q - U)) + Log[mD]*((2*Log[Q])/((mD - Q)*(Q - U)) + 
            (2*Log[U])/((mD - Q)*(-Q + U))) + 
          (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2 - 
             delta[A0, Q, mD])*phi[A0, Q, mD])/(mD^2*(mD - Q)*(Q - U)) + 
          (2*(-(A0*(A0 - 5*Q)) + delta[A0, Q, Q])*phi[A0, Q, Q])/
           ((mD - Q)*Q^2*(Q - U)) - (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
             3*mD*U + U^2 - delta[A0, U, mD])*phi[A0, U, mD])/
           (mD^2*(mD - Q)*(Q - U)) + (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
             3*Q*U + U^2 - delta[A0, U, Q])*phi[A0, U, Q])/
           ((mD - Q)*Q^2*(Q - U))) + Nc*Xb^3*Yb*
         ((2*(-2*A0 + mD + 3*Q)*Log[Q]^2)/((mD - Q)^3*(Q - U)) + 
          Log[A0]*((4*Log[mD])/(mD - Q)^3 + (4*Log[Q])/(-mD + Q)^3) - 
          (2*(-2*A0 + mD + Q + 2*U)*Log[Q]*Log[U])/((mD - Q)^3*(Q - U)) + 
          Log[mD]*((2*(-2*A0 + mD + 3*Q)*Log[Q])/((mD - Q)^3*(-Q + U)) + 
            (2*(-2*A0 + mD + Q + 2*U)*Log[U])/((mD - Q)^3*(Q - U))) + 
          (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*Q - 3*mD*Q + Q^2) + 
             (-3*mD + Q)*delta[A0, Q, mD])*phi[A0, Q, mD])/
           (mD^2*(mD - Q)^3*(Q - U)) + (2*(A0*(A0 - 5*Q)*(mD - Q) - 
             (mD - 3*Q)*delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^3*Q^2*
            (Q - U)) - (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*
                U + U^2) + (-3*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
           (mD^2*(mD - Q)^3*(Q - U)) + 
          (2*(-((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2)) + 
             (mD - 3*Q)*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)^3*Q^2*
            (Q - U))) + Nc*Xb^2*((2*(A0 - mD - Q)*Log[Q]^2)/
           ((mD - Q)^2*(Q - U)) + Log[A0]*((-2*Log[mD])/(mD - Q)^2 + 
            (2*Log[Q])/(mD - Q)^2) + (2*(-A0 + mD + U)*Log[Q]*Log[U])/
           ((mD - Q)^2*(Q - U)) + Log[mD]*((2*(-A0 + mD + Q)*Log[Q])/
             ((mD - Q)^2*(Q - U)) - (2*(-A0 + mD + U)*Log[U])/
             ((mD - Q)^2*(Q - U))) + (2*delta[A0, Q, mD]*phi[A0, Q, mD])/
           (mD*(mD - Q)^2*(Q - U)) + (2*(A0*(A0 - 5*Q)*(-mD + Q) + 
             (mD - 2*Q)*delta[A0, Q, Q])*phi[A0, Q, Q])/((mD - Q)^2*Q^2*
            (Q - U)) - (2*delta[A0, U, mD]*phi[A0, U, mD])/
           (mD*(mD - Q)^2*(Q - U)) + 
          (2*((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2) - 
             (mD - 2*Q)*delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)^2*Q^2*
            (Q - U))))) + Nc*Xb^4*(16/(mD - Q)^2 - (4*(mD + Q)*Log[Q]^2)/
       (mD - Q)^3 + (8*Log[q2])/(mD - Q)^2 + 
      Log[mD]*((-2*(7*mD + Q))/(mD - Q)^3 + (4*(mD + Q)*Log[Q])/(mD - Q)^3 - 
        (4*(mD + Q)*Log[q2])/(mD - Q)^3) + 
      Log[Q]*((2*(3*mD + 5*Q))/(mD - Q)^3 + (4*(mD + Q)*Log[q2])/
         (mD - Q)^3) - (6*PolyLog[2, 1 - mD/Q])/(mD - Q)^2) + 
    Nc*Xb^2*(8/(-mD + Q) + (2*(A0 + 3*mD - 3*Q)*Log[Q]^2)/(mD - Q)^2 + 
      Log[A0]*((2*(mD - 2*Q + U)*Log[mD])/(mD - Q)^2 - 
        (2*(mD - 2*Q + U)*Log[Q])/(mD - Q)^2) + (4*Log[q2])/(-mD + Q) + 
      Log[mD]*((2*(3*mD + Q))/(mD - Q)^2 - (2*(A0 + 3*mD - 3*Q)*Log[Q])/
         (mD - Q)^2 + (4*mD*Log[q2])/(mD - Q)^2 - (2*(-A0 + Q + U)*Log[U])/
         (mD - Q)^2) + Log[Q]*((-2*(mD + 3*Q))/(mD - Q)^2 - 
        (4*mD*Log[q2])/(mD - Q)^2 + (2*(-A0 + Q + U)*Log[U])/(mD - Q)^2) + 
      (2*delta[A0, Q, mD]*phi[A0, Q, mD])/(mD*(mD - Q)^2) + 
      (2*(A0*(A0 - 5*Q)*(-mD + Q) + (mD - 2*Q)*delta[A0, Q, Q])*
        phi[A0, Q, Q])/((mD - Q)^2*Q^2) + 
      (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
         (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/(mD^2*(mD - Q)^2) + 
      (2*delta[A0, U, Q]*phi[A0, U, Q])/((mD - Q)^2*Q) + 
      (12*PolyLog[2, 1 - mD/Q])/(mD - Q)) + 
    (cbe^2*(Nc*Xb^4*((A0*(mD + 5*Q))/(Q*(-mD + Q)^3) + 
          Log[A0]*((A0*(mD + 5*Q))/((mD - Q)^3*Q) - (2*A0*(2*mD + Q)*Log[mD])/
             (mD - Q)^4 + (2*A0*(2*mD + Q)*Log[Q])/(mD - Q)^4) + 
          (A0*(mD + 5*Q)*Log[q2])/(Q*(-mD + Q)^3) + 
          Log[Q]*((-2*A0*(2*mD + Q))/(mD - Q)^4 - (2*A0*(2*mD + Q)*Log[q2])/
             (mD - Q)^4) + Log[mD]*((2*A0*(2*mD + Q))/(mD - Q)^4 + 
            (2*A0*(2*mD + Q)*Log[q2])/(mD - Q)^4)) + 
        Nc*(1/2 + Pi^2/3 - A0/Q + Log[A0]^2 + (1 - A0/Q)*Log[q2] + 
          2*Log[q2]^2 + Log[A0]*((A0 + Q)/Q + Log[mD] - 2*Log[q2] - Log[U]) + 
          Log[mD]*(-2 - 2*Log[q2] + Log[U]) + 
          ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2 - 
             delta[A0, U, mD])*phi[A0, U, mD])/mD^2) + 
        Nc*Xb^2*((2*(A0 - 2*Q))/((mD - Q)*Q) + Log[A0]*
           ((2*A0)/(Q*(-mD + Q)) + ((A0 - Q + U)*Log[mD])/(mD - Q)^2 - 
            ((A0 - Q + U)*Log[Q])/(mD - Q)^2) + (2*(A0 - Q)*Log[q2])/
           ((mD - Q)*Q) + Log[mD]*((2*(-A0 + mD + Q))/(mD - Q)^2 - 
            (2*(A0 - Q)*Log[q2])/(mD - Q)^2 - ((-A0 + Q + U)*Log[U])/
             (mD - Q)^2) + Log[Q]*((2*(A0 - 2*Q))/(mD - Q)^2 + 
            (2*(A0 - Q)*Log[q2])/(mD - Q)^2 + ((-A0 + Q + U)*Log[U])/
             (mD - Q)^2) + (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*
                mD*U + U^2) + (-2*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
           (mD^2*(mD - Q)^2) + (delta[A0, U, Q]*phi[A0, U, Q])/
           ((mD - Q)^2*Q)) + 
        Yt*(Nc*Xb*(Log[A0]*((2*Log[mD])/(mD - Q) + (2*Log[Q])/(-mD + Q)) + 
            Log[mD]*(4/(-mD + Q) + (4*Log[q2])/(-mD + Q) + (2*Log[U])/(mD - 
                Q)) + Log[Q]*(4/(mD - Q) + (4*Log[q2])/(mD - Q) + 
              (2*Log[U])/(-mD + Q)) + (2*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 
               3*mD*U + U^2 - delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)) - (2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^
                2 - delta[A0, U, Q])*phi[A0, U, Q])/((mD - Q)*Q^2)) + 
          Nc*Xb^3*(-16/(mD - Q)^2 + Log[A0]*((-2*(2*A0 + mD + Q - 2*U)*
                Log[mD])/(mD - Q)^3 + (2*(2*A0 + mD + Q - 2*U)*Log[Q])/
               (mD - Q)^3) - (8*Log[q2])/(mD - Q)^2 + 
            Log[mD]*((4*(3*mD + Q))/(mD - Q)^3 + (4*(mD + Q)*Log[q2])/
               (mD - Q)^3 - (2*(-2*A0 + mD + Q + 2*U)*Log[U])/(mD - Q)^3) + 
            Log[Q]*((-4*(mD + 3*Q))/(mD - Q)^3 - (4*(mD + Q)*Log[q2])/
               (mD - Q)^3 + (2*(-2*A0 + mD + Q + 2*U)*Log[U])/(mD - Q)^3) + 
            (2*((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
               (-3*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)^3) + (2*(-((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 
                  2*A0*U - 3*Q*U + U^2)) + (mD - 3*Q)*delta[A0, U, Q])*
              phi[A0, U, Q])/(Q^2*(-mD + Q)^3))) + 
        Yt^2*(Nc*Xb^4*((mD + 11*Q)/(Q*(-mD + Q)^3) + 
            Log[A0]*((-((A0 + Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                   3*Q*U + U^2)) + (A0 + 2*Q - U)*delta[A0, U, Q])/(2*
                (mD - Q)^2*Q^2*delta[A0, U, Q]) - ((3*A0 + 2*mD + Q - 3*U)*
                Log[mD])/(mD - Q)^4 + ((3*A0 + 2*mD + Q - 3*U)*Log[Q])/
               (mD - Q)^4) + ((mD + 5*Q)*Log[q2])/(Q*(-mD + Q)^3) + 
            (((A0 - Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2) + 
               (-A0 + 2*Q + U)*delta[A0, U, Q])*Log[U])/(2*(mD - Q)^2*Q^2*
              delta[A0, U, Q]) + Log[mD]*((2*(5*mD + Q))/(mD - Q)^4 + 
              (2*(2*mD + Q)*Log[q2])/(mD - Q)^4 - ((-3*A0 + 2*mD + Q + 3*U)*
                Log[U])/(mD - Q)^4) + Log[Q]*(-((mD^2 + 2*mD*Q + 9*Q^2 - 
                 ((mD - Q)^2*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2))/
                  delta[A0, U, Q])/((mD - Q)^4*Q)) - (2*(2*mD + Q)*Log[q2])/
               (mD - Q)^4 + ((-3*A0 + 2*mD + Q + 3*U)*Log[U])/(mD - Q)^4) + 
            (((mD - Q)*(A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*U + U^2) + 
               (-4*mD + Q)*delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)^4) + ((-((mD - Q)*(A0^2*(3*mD - 7*Q) + 
                  2*(mD - 5*Q)*Q^2 - 6*(mD - 3*Q)*Q*U + (3*mD - 7*Q)*U^2 + 
                  2*A0*(-3*mD*(Q + U) + Q*(9*Q + 7*U)))) + ((mD - Q)^2*
                 (A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + U^2)*(A0^2 + 
                  U*(-Q + U) - A0*(Q + 2*U)))/delta[A0, U, Q] + 2*
                (mD^2 - 4*mD*Q + 6*Q^2)*delta[A0, U, Q])*phi[A0, U, Q])/
             (2*(mD - Q)^4*Q^3)) + Nc*(-Q^(-1) + 
            ((-((A0 + Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + 
                  U^2)) + (A0 + 2*Q - U)*delta[A0, U, Q])*Log[A0])/
             (2*Q^2*delta[A0, U, Q]) + ((A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*
                U + U^2 - delta[A0, U, Q])*Log[Q])/(Q*delta[A0, U, Q]) - 
            Log[q2]/Q + (((A0 - Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                 3*Q*U + U^2) + (-A0 + 2*Q + U)*delta[A0, U, Q])*Log[U])/
             (2*Q^2*delta[A0, U, Q]) + (((A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                 3*Q*U + U^2)*(A0^2 + U*(-Q + U) - A0*(Q + 2*U)) + 
               delta[A0, U, Q]*(-3*A0^2 - 2*Q^2 + 6*Q*U - 3*U^2 + 
                 6*A0*(Q + U) + 2*delta[A0, U, Q]))*phi[A0, U, Q])/
             (2*Q^3*delta[A0, U, Q])) + Nc*Xb^2*(2/((mD - Q)*Q) + 
            Log[A0]*(((A0 + Q - U)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + 
                  U^2) + (-A0 - 2*Q + U)*delta[A0, U, Q])/((mD - Q)*Q^2*
                delta[A0, U, Q]) + Log[mD]/(mD - Q)^2 - Log[Q]/(mD - Q)^2) + 
            (2*Log[q2])/(mD*Q - Q^2) + ((-((A0 - Q - U)*(A0^2 - 3*A0*Q + 
                  2*Q^2 - 2*A0*U - 3*Q*U + U^2)) + (A0 - 2*Q - U)*
                delta[A0, U, Q])*Log[U])/((mD - Q)*Q^2*delta[A0, U, Q]) + 
            Log[Q]*((2*(-((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 3*Q*U + 
                    U^2)) + mD*delta[A0, U, Q]))/((mD - Q)^2*Q*delta[A0, U, 
                 Q]) + (2*Log[q2])/(mD - Q)^2 - Log[U]/(mD - Q)^2) + 
            Log[mD]*(-2/(mD - Q)^2 - (2*Log[q2])/(mD - Q)^2 + 
              Log[U]/(mD - Q)^2) + ((A0^2 - 3*A0*mD + 2*mD^2 - 2*A0*U - 3*mD*
                U + U^2 - delta[A0, U, mD])*phi[A0, U, mD])/
             (mD^2*(mD - Q)^2) + (((mD - Q)*(A0^2 - 3*A0*Q + 2*Q^2 - 2*A0*U - 
                 3*Q*U + U^2)*(-A0^2 + (Q - U)*U + A0*(Q + 2*U)) + 
               delta[A0, U, Q]*(A0^2*(3*mD - 4*Q) + 2*(mD - 2*Q)*Q^2 + 
                 3*Q*(-2*mD + 3*Q)*U + (3*mD - 4*Q)*U^2 + A0*(-6*mD*(Q + U) + 
                   Q*(9*Q + 8*U)) + (-2*mD + 3*Q)*delta[A0, U, Q]))*
              phi[A0, U, Q])/((mD - Q)^2*Q^3*delta[A0, U, Q])))) + 
      Nc*Xb^2*((2*(-2*mu2 + 3*Q + U))/((mD - Q)*Q) + (2*(mD - mu2)*Log[mD]^2)/
         (mD - Q)^2 + (4*mu2*Log[mu2])/(mD*Q - Q^2) - (2*(mD - mu2)*Log[Q]^2)/
         (mD - Q)^2 + (2*(-2*mu2 + 2*Q + U)*Log[q2])/((mD - Q)*Q) + 
        (2*U*Log[U])/(Q*(-mD + Q)) + Log[Q]*((2*(mD - 2*mu2 + 2*Q + U))/
           (mD - Q)^2 + (2*(2*mD - 2*mu2 + U)*Log[q2])/(mD - Q)^2 - 
          (2*U*Log[U])/(mD - Q)^2) + Log[mD]*((-2*(3*mD - 2*mu2 + U))/
           (mD - Q)^2 - (2*(2*mD - 2*mu2 + U)*Log[q2])/(mD - Q)^2 + 
          (2*U*Log[U])/(mD - Q)^2) + (4*(mD - mu2)*PolyLog[2, 1 - mu2/mD])/
         (mD - Q)^2 - (4*(mD - mu2)*PolyLog[2, 1 - mu2/Q])/(mD - Q)^2) + 
      Nc*Xb^4*((mD*(-2*mu2 + 9*Q + U) + Q*(-16*mu2 + 3*Q + 5*U))/
         (Q*(-mD + Q)^3) + ((-mD + mu2)*(mD - 3*mu2 + 2*Q)*Log[mD]^2)/
         (mD - Q)^4 + ((mD - mu2)*(mD - 3*mu2 + 2*Q)*Log[Q]^2)/(mD - Q)^4 + 
        Log[mu2]*((2*mu2*(mD + 2*Q))/(Q*(-mD + Q)^3) - (6*mu2^2*Log[Q])/
           (mD - Q)^4) + ((mD*(-2*mu2 + 5*Q + U) + Q*(-10*mu2 + Q + 5*U))*
          Log[q2])/(Q*(-mD + Q)^3) + ((mD + 5*Q)*U*Log[U])/((mD - Q)^3*Q) + 
        Log[mD]*((2*(3*mD^2 + Q*(-2*mu2 + U) + mD*(-7*mu2 + 3*Q + 2*U)))/
           (mD - Q)^4 + (6*mu2^2*Log[mu2])/(mD - Q)^4 + 
          (2*(mD^2 + Q*(-2*mu2 + U) + 2*mD*(-2*mu2 + Q + U))*Log[q2])/
           (mD - Q)^4 - (2*(2*mD + Q)*U*Log[U])/(mD - Q)^4) + 
        Log[Q]*(-((mD^2 + Q*(-10*mu2 + Q + 2*U) + mD*(-8*mu2 + 10*Q + 4*U))/
            (mD - Q)^4) - (2*(mD^2 + Q*(-2*mu2 + U) + 2*mD*(-2*mu2 + Q + U))*
            Log[q2])/(mD - Q)^4 + (2*(2*mD + Q)*U*Log[U])/(mD - Q)^4) - 
        (2*(mD - mu2)*(mD - 3*mu2 + 2*Q)*PolyLog[2, 1 - mu2/mD])/(mD - Q)^4 + 
        (2*(mD - mu2)*(mD - 3*mu2 + 2*Q)*PolyLog[2, 1 - mu2/Q])/(mD - Q)^4) + 
      Nc*(-1/2 - U/Q + mu2*(2/Q + (-mu2 + U)^(-1)) + Log[Q]^2 + 
        Log[Q]*(-1 - 2*Log[q2]) + 2*Log[q2]^2 + 
        (U*(mu2^2 + 2*mu2*(Q - U) + U*(Q + U))*Log[U])/(Q*(mu2 - U)^2) + 
        (1 - (2*mu2^2)/(mu2 - U)^2)*Log[U]^2 + 
        Log[q2]*((2*mu2*(mu2 + Q) - 3*mu2*U + U^2)/(Q*(mu2 - U)) + 
          (-2 + (2*mu2^2)/(mu2 - U)^2)*Log[U]) + 
        Log[mu2]*(-((mu2*(2*mu2^2 + mu2*(Q - 4*U) + 2*U*(Q + U)))/
            (Q*(mu2 - U)^2)) - (2*mu2^2*Log[q2])/(mu2 - U)^2 + 
          (2*mu2^2*Log[U])/(mu2 - U)^2) + 2*PolyLog[2, 1 - mu2/Q] + 
        (2 - (4*mu2^2)/(mu2 - U)^2)*PolyLog[2, 1 - mu2/U]))/sbe^2 + 
    (Nc*Xb*((4*Sqrt[mu2]*sgn*Log[mD]^2)/(-mD + Q) + 
        (4*Sqrt[mu2]*sgn*Log[Q]^2)/(mD - Q) + 
        Log[mD]*((8*Sqrt[mu2]*sgn)/(mD - Q) + (8*Sqrt[mu2]*sgn*Log[q2])/
           (mD - Q)) + Log[Q]*((8*Sqrt[mu2]*sgn)/(-mD + Q) + 
          (8*Sqrt[mu2]*sgn*Log[q2])/(-mD + Q)) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/mD])/(-mD + Q) + 
        (8*Sqrt[mu2]*sgn*PolyLog[2, 1 - mu2/Q])/(mD - Q)) + 
      Nc*Xb^3*((32*Sqrt[mu2]*sgn)/(mD - Q)^2 + 
        (4*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*Log[mD]^2)/(mD - Q)^3 + 
        (16*mu2^(3/2)*sgn*Log[mu2]*Log[Q])/(-mD + Q)^3 - 
        (4*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*Log[Q]^2)/(mD - Q)^3 + 
        (16*Sqrt[mu2]*sgn*Log[q2])/(mD - Q)^2 + 
        Log[mD]*((8*Sqrt[mu2]*(3*mD + Q)*sgn)/(-mD + Q)^3 + 
          (16*mu2^(3/2)*sgn*Log[mu2])/(mD - Q)^3 - 
          (8*Sqrt[mu2]*(mD + Q)*sgn*Log[q2])/(mD - Q)^3) + 
        Log[Q]*((8*Sqrt[mu2]*(mD + 3*Q)*sgn)/(mD - Q)^3 + 
          (8*Sqrt[mu2]*(mD + Q)*sgn*Log[q2])/(mD - Q)^3) + 
        (8*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*PolyLog[2, 1 - mu2/mD])/
         (mD - Q)^3 - (8*Sqrt[mu2]*(mD - 2*mu2 + Q)*sgn*
          PolyLog[2, 1 - mu2/Q])/(mD - Q)^3) + 
      Nc*Xt*((4*Sqrt[mu2]*(mu2 + Q)*sgn*Log[Q]^2)/((mu2 - Q)*(-Q + U)) + 
        (8*Sqrt[mu2]*Q*sgn*Log[Q]*Log[q2])/((mu2 - Q)*(Q - U)) + 
        (8*Sqrt[mu2]*sgn*U*Log[q2]*Log[U])/((mu2 - U)*(-Q + U)) - 
        (4*Sqrt[mu2]*sgn*(mu2 + U)*Log[U]^2)/((mu2 - U)*(-Q + U)) + 
        Log[mu2]*((8*mu2^(3/2)*sgn*Log[Q])/((mu2 - Q)*(Q - U)) + 
          (8*mu2^(3/2)*sgn*Log[q2])/((mu2 - Q)*(-mu2 + U)) + 
          (8*mu2^(3/2)*sgn*Log[U])/((mu2 - U)*(-Q + U))) + 
        (8*Sqrt[mu2]*(mu2 + Q)*sgn*PolyLog[2, 1 - mu2/Q])/
         ((mu2 - Q)*(-Q + U)) - (8*Sqrt[mu2]*sgn*(mu2 + U)*
          PolyLog[2, 1 - mu2/U])/((mu2 - U)*(-Q + U))))/(cbe*sbe) + 
    (Nc*(-2*Log[mD]^2 + 4*mu2*((mD - mu2)^(-1) + (-mu2 + Q)^(-1))*Log[mu2] - 
        2*Log[Q]^2 - 2*Log[q2]^2 + Log[Q]*((4*mu2)/(mu2 - Q) + 2*Log[q2]) + 
        Log[mD]*((4*mu2)/(-mD + mu2) + 2*Log[Q] + 2*Log[q2]) - 
        4*PolyLog[2, 1 - mu2/mD] - 4*PolyLog[2, 1 - mu2/Q]) + 
      Nc*Xb^2*(4/(mD - Q) + (2*Q*Log[mD]^2)/(mD - Q)^2 + 
        (4*mu2*Log[mu2])/((mD - mu2)*(mD - Q)) + (2*Log[q2])/(mD - Q) + 
        Log[mD]*((-4*(mD^2 - mu2*Q))/((mD - mu2)*(mD - Q)^2) - 
          (2*Q*Log[Q])/(mD - Q)^2 - (2*Q*Log[q2])/(mD - Q)^2) + 
        Log[Q]*((2*(mD + Q))/(mD - Q)^2 + (2*Q*Log[q2])/(mD - Q)^2) + 
        (4*Q*PolyLog[2, 1 - mu2/mD])/(mD - Q)^2 - (4*Q*PolyLog[2, 1 - mu2/Q])/
         (mD - Q)^2) + Xt^2*(Nc*Xb^2*(2/((mD - Q)*(Q - U)) - 
          (2*Q*U*Log[Q]^2)/((mD - Q)^2*(Q - U)^2) + (2*U*Log[U])/
           ((mD - Q)*(Q - U)^2) + Log[mD]*((2*Q)/((mD - Q)^2*(-Q + U)) + 
            (2*Q*U*Log[Q])/((mD - Q)^2*(Q - U)^2) - (2*Q*U*Log[U])/
             ((mD - Q)^2*(Q - U)^2)) + Log[Q]*((2*(Q^2 - mD*U))/
             ((mD - Q)^2*(Q - U)^2) + (2*Q*U*Log[U])/((mD - Q)^2*
              (Q - U)^2))) + Nc*(4/(Q - U) + (4*mu2*Log[mu2])/
           ((mu2 - Q)*(-Q + U)) + (2*U*Log[Q]^2)/(Q - U)^2 + 
          Log[Q]*((4*(Q^2 - mu2*U))/((mu2 - Q)*(Q - U)^2) - 
            (2*U*Log[q2])/(Q - U)^2) + (4*U*Log[U])/(Q - U)^2 - 
          (2*U*Log[U]^2)/(Q - U)^2 + Log[q2]*(2/(Q - U) + 
            (2*U*Log[U])/(Q - U)^2) + Log[mD]*(2/(Q - U) - 
            (2*U*Log[Q])/(Q - U)^2 + (2*U*Log[U])/(Q - U)^2) + 
          (4*U*PolyLog[2, 1 - mu2/Q])/(Q - U)^2 - (4*U*PolyLog[2, 1 - mu2/U])/
           (Q - U)^2)))/cbe^2 + 
    Xt^2*(Nc^2*((2*Q*U*Log[Q]^2)/(Q - U)^3 + Log[mD]*(-((Q + U)/(Q - U)^2) + 
          (2*Q*U*Log[Q])/(Q - U)^3 - (2*Q*U*Log[U])/(Q - U)^3) + 
        Log[Q]*(-((Q + U)/(Q - U)^2) - (4*Q*U*Log[q2])/(Q - U)^3 - 
          (2*Q*U*Log[U])/(Q - U)^3) + Log[q2]*((2*(Q + U))/(Q - U)^2 + 
          (4*Q*U*Log[U])/(Q - U)^3)) + Nc*((-5*Q + U)/(Q^2 - Q*U) + 
        ((-3*Q + U)*Log[Q]^2)/(Q - U)^2 + (U*(-5*Q + U)*Log[U])/
         (Q*(Q - U)^2) + (3*Log[U]^2)/(-Q + U) + 
        Log[Q]*((3*Q + U)/(Q - U)^2 + (2*U*Log[q2])/(Q - U)^2 + 
          ((6*Q - 4*U)*Log[U])/(Q - U)^2) + Log[q2]*((-3*Q + U)/(Q^2 - Q*U) - 
          (2*U*Log[U])/(Q - U)^2) + (6*PolyLog[2, 1 - Q/U])/(-Q + U)) + 
      Xb^4*(Nc^2*((-2*(Q + U))/((mD - Q)^2*(Q - U)^2) + 
          (2*Q*(mD + Q)*U*Log[Q]^2)/((mD - Q)^3*(Q - U)^3) - 
          (4*Q*U*Log[U])/((mD - Q)^2*(Q - U)^3) + 
          Log[mD]*(((mD + Q)*(Q + U))/((mD - Q)^3*(Q - U)^2) + 
            (2*Q*(mD + Q)*U*Log[Q])/((-mD + Q)^3*(Q - U)^3) + 
            (2*Q*(mD + Q)*U*Log[U])/((mD - Q)^3*(Q - U)^3)) + 
          Log[Q]*((-(Q^2*(mD + Q)) + 4*(mD - Q)*Q*U + (mD + Q)*U^2)/
             ((mD - Q)^3*(Q - U)^3) + (2*Q*(mD + Q)*U*Log[U])/
             ((-mD + Q)^3*(Q - U)^3))) + Nc*((mD + 11*Q)/(Q*(-mD + Q)^3) + 
          ((-2*Q*(3*mD + Q) + 4*(mD + Q)*U)*Log[Q]^2)/((mD - Q)^3*
            (Q - U)^2) + ((mD + 5*Q)*Log[q2])/(Q*(-mD + Q)^3) + 
          (U*Log[U])/((mD - Q)^2*Q*(-Q + U)) + 
          Log[Q]*((-4*(mD + 2*Q))/(mD - Q)^4 + 3/((mD - Q)^2*(Q - U)) - 
            (2*(2*mD + Q)*Log[q2])/(mD - Q)^4 - 
            (2*(mD - U)*(-3*mD*Q + 2*mD*U + Q*U)*Log[U])/((mD - Q)^4*
              (Q - U)^2)) + Log[mD]*((-2*(mD^2 - 7*mD*Q + 5*mD*U + Q*U))/
             ((mD - Q)^4*(Q - U)) + ((2*Q*(3*mD + Q) - 4*(mD + Q)*U)*Log[Q])/
             ((mD - Q)^3*(Q - U)^2) + (2*(2*mD + Q)*Log[q2])/(mD - Q)^4 + 
            (2*(mD - U)*(-3*mD*Q + 2*mD*U + Q*U)*Log[U])/((mD - Q)^4*
              (Q - U)^2)) + (6*PolyLog[2, 1 - mD/Q])/((mD - Q)^2*(-Q + U)) + 
          (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/((mD - Q)^4*(Q - U)) - 
          (6*(mD - U)^2*PolyLog[2, 1 - Q/U])/((mD - Q)^4*(Q - U)))) + 
      Xb^2*(Nc^2*((4*Q*U*Log[Q]^2)/((-mD + Q)*(Q - U)^3) + 
          Log[Q]*((2*(Q + U))/((mD - Q)*(Q - U)^2) + (4*Q*U*Log[U])/
             ((mD - Q)*(Q - U)^3)) + Log[mD]*
           ((2*(Q + U))/((-mD + Q)*(Q - U)^2) + (4*Q*U*Log[Q])/
             ((mD - Q)*(Q - U)^3) + (4*Q*U*Log[U])/((-mD + Q)*(Q - U)^3))) + 
        Nc*(2/(mD*Q - Q^2) + ((3*mD^2 + 6*mD*Q - 5*Q^2 - 8*mD*U + 4*Q*U)*
            Log[Q]^2)/((mD - Q)^2*(Q - U)^2) + (2*Log[q2])/(mD*Q - Q^2) + 
          (2*U*Log[U])/((mD - Q)*Q*(Q - U)) - (3*Log[U]^2)/(Q - U)^2 + 
          Log[mD]*((2*(mD - U))/((mD - Q)^2*(-Q + U)) + 
            ((-6*mD^2 + 2*Q*(Q - 2*U) + 8*mD*U)*Log[Q])/((mD - Q)^2*
              (Q - U)^2) - (2*Log[q2])/(mD - Q)^2 + 
            (2*(3*mD^2 - 4*mD*U + U^2)*Log[U])/((mD - Q)^2*(Q - U)^2)) + 
          Log[Q]*(2/(mD - Q)^2 + (2*Log[q2])/(mD - Q)^2 - 
            (2*(6*mD*Q - 3*Q^2 - 4*mD*U + U^2)*Log[U])/((mD - Q)^2*
              (Q - U)^2)) + (6*(mD + Q - 2*U)*PolyLog[2, 1 - mD/Q])/
           ((mD - Q)*(Q - U)^2) - (6*(mD - U)^2*PolyLog[2, 1 - mD/U])/
           ((mD - Q)^2*(Q - U)^2) + ((12*mD - 6*(Q + U))*PolyLog[2, 1 - Q/U])/
           ((mD - Q)^2*(Q - U))))))) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
    ];


(* Tau Yukawa lambda 2-loop threshold correction, computation with general masses *)
lambda2LHSSTau = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q3 = msq2[3,3]*(1+0.03),
    D3 = msd2[3,3]*(1-0.03),
    L3 = msl2[3,3]*(1-0.01),
    R3 = mse2[3,3]*(1+0.01),
    mu2 = MuInput^2*(1-0.05),
    q2 = SCALE^2, (* renormalization/matching scale *)
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    gtau = Ye[3,3],
    A0 = mAInput^2*(1-0.015),
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb,
    Xtau = xtau,
    Ytau = ytau,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    Q2 = SCALE^2,
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];
    k^2*(gtau^6*(5 + 3*Log[L3]^2 + Log[L3]*(-6 - 6*Log[q2]) + 5*Log[q2]^2 + 
    Log[q2]*(10 - 4*Log[R3]) - 4*Log[R3] + 2*Log[R3]^2 + 
    (sbe^2*(3/2 - A0/L3 + Pi^2 - (2*A0)/R3 + 3*Log[A0]^2 + 3*Log[L3]^2 + 
       Log[L3]*(-9/2 - 3*Log[q2]) + 2*Log[q2]^2 + 
       Log[A0]*(7 + A0*(L3^(-1) + 2/R3) - 3*Log[L3] - 3*Log[R3]) + 
       Log[q2]*(-((A0*(2*L3 + R3))/(L3*R3)) - Log[R3]) - (5*Log[R3])/2 + 
       2*Log[R3]^2 + ((A0^2 - 7*A0*L3 - delta[A0, L3, L3])*phi[A0, L3, L3])/
        L3^2 + Ytau^2*(-L3^(-1) - 2/R3 + 
         ((((L3 - R3)^4 - A0^3*(L3 + R3) - A0*(L3 - R3)^2*(3*L3 + 5*R3) + 
              A0^2*(3*L3^2 + 4*L3*R3 + 5*R3^2))*delta[A0, R3, L3] + 
            delta[A0, L3, R3]*(-(L3*(-A0^2 + (L3 - R3)^2)*R3) + 
              (-L3^2 + 3*L3*R3 - R3^2 + A0*(L3 + R3))*delta[A0, R3, L3]))*
           Log[A0])/(2*L3*R3^2*delta[A0, L3, R3]*delta[A0, R3, L3]) + 
         (((A0^3 - (L3 - R3)^3 - 3*A0^2*(L3 + R3) + 3*A0*(L3^2 - R3^2))*
             delta[A0, R3, L3] + delta[A0, L3, R3]*
             (-(R3*(A0^2 - L3^2 - 2*A0*R3 + R3^2)) + (-A0 + L3 + 2*R3)*delta[
                A0, R3, L3]))*Log[L3])/(2*R3^2*delta[A0, L3, R3]*
           delta[A0, R3, L3]) + (-L3^(-1) - 2/R3)*Log[q2] + 
         (((A0^3 + (L3 - R3)^3 - A0^2*(L3 + 5*R3) - A0*(L3^2 + 4*L3*R3 - 
                5*R3^2))*delta[A0, R3, L3] - delta[A0, L3, R3]*
             (2*L3*(A0 + L3 - R3)*R3 + (A0 + L3 - 3*R3)*delta[A0, R3, L3]))*
           Log[R3])/(2*L3*R3*delta[A0, L3, R3]*delta[A0, R3, L3]) + 
         ((A0^4 + (L3 - R3)^4 - 4*A0^3*(L3 + R3) - 4*A0*(L3 - R3)^2*
             (L3 + R3) + 2*A0^2*(3*L3^2 + 2*L3*R3 + R3^2) - 
            3*(A0^2 + (L3 - R3)^2 - 2*A0*(L3 + R3))*delta[A0, L3, R3] + 
            2*delta[A0, L3, R3]^2)*phi[A0, L3, R3])/
          (2*R3^3*delta[A0, L3, R3]) + 
         ((-(A0 + L3 - R3)^2 + delta[A0, R3, L3])*phi[A0, R3, L3])/
          (2*L3*delta[A0, R3, L3])) + ((A0^2 - 6*A0*R3 - delta[A0, R3, R3])*
         phi[A0, R3, R3])/R3^2 + Xtau*Ytau*((6*Log[L3]^2)/(L3 - R3) + 
         (12*Log[R3])/(L3 - R3) + (12*Log[q2]*Log[R3])/(L3 - R3) - 
         (4*Log[R3]^2)/(L3 - R3) + Log[A0]*((2*Log[L3])/(L3 - R3) - 
           (2*Log[R3])/(L3 - R3)) + Log[L3]*(-12/(L3 - R3) - 
           (12*Log[q2])/(L3 - R3) - (2*Log[R3])/(L3 - R3)) - 
         (2*(-(A0*(A0 - 7*L3)) + delta[A0, L3, L3])*phi[A0, L3, L3])/
          (L3^2*(L3 - R3)) + (2*(A0 + L3 - R3)*phi[A0, R3, L3])/
          (L3*(L3 - R3)) + (2*(-(A0*(A0 - 6*R3)) + delta[A0, R3, R3])*
           phi[A0, R3, R3])/((L3 - R3)*R3^2)) + Xtau^3*Ytau*
        (-48/(L3 - R3)^2 + (2*(4*A0 - 5*L3 - 3*R3)*Log[L3]^2)/(L3 - R3)^3 - 
         (12*(L3 + 3*R3)*Log[R3])/(L3 - R3)^3 + (4*(-A0 + L3 + R3)*Log[R3]^2)/
          (L3 - R3)^3 + Log[A0]*((2*(-6*A0 + L3 - R3)*Log[L3])/(L3 - R3)^3 + 
           (2*(6*A0 - L3 + R3)*Log[R3])/(L3 - R3)^3) + 
         Log[q2]*(-24/(L3 - R3)^2 - (12*(L3 + R3)*Log[R3])/(L3 - R3)^3) + 
         Log[L3]*((12*(3*L3 + R3))/(L3 - R3)^3 + (12*(L3 + R3)*Log[q2])/
            (L3 - R3)^3 + (2*(-2*A0 + 3*L3 + R3)*Log[R3])/(L3 - R3)^3) + 
         (2*(A0*(A0 - 7*L3)*(L3 - R3) + (-5*L3 + R3)*delta[A0, L3, L3])*
           phi[A0, L3, L3])/(L3^2*(L3 - R3)^3) - 
         (2*((L3 - R3)*(A0 + L3 - R3) - 2*delta[A0, R3, L3])*phi[A0, R3, L3])/
          (L3*(L3 - R3)^3) + (2*(A0*(A0 - 6*R3)*(L3 - R3) - 
            (L3 - 3*R3)*delta[A0, R3, R3])*phi[A0, R3, R3])/
          ((L3 - R3)^3*R3^2)) + Xtau^2*((4*A0*L3 - 2*A0*R3 - 4*L3*R3)/
          (L3^2*R3 - L3*R3^2) + ((A0 + L3 - 3*R3)*Log[L3]^2)/(L3 - R3)^2 + 
         ((2*A0 + L3 - 5*R3)*Log[R3])/(L3 - R3)^2 - (2*Log[R3]^2)/(L3 - R3) + 
         Log[q2]*((4*A0*L3 - 2*A0*R3 - 2*L3*R3)/(L3^2*R3 - L3*R3^2) + 
           (2*(A0 - L3)*Log[R3])/(L3 - R3)^2) + 
         Log[L3]*((-2*A0 + L3 + 3*R3)/(L3 - R3)^2 - (2*(A0 - L3)*Log[q2])/
            (L3 - R3)^2 + ((-A0 + L3 + R3)*Log[R3])/(L3 - R3)^2) + 
         Log[A0]*((2*A0*(-2*L3 + R3))/(L3*(L3 - R3)*R3) + 
           ((A0 - 5*L3 + 5*R3)*Log[L3])/(L3 - R3)^2 - 
           ((A0 - 5*L3 + 5*R3)*Log[R3])/(L3 - R3)^2) + 
         ((A0*(A0 - 7*L3)*(L3 - R3) + (-2*L3 + R3)*delta[A0, L3, L3])*
           phi[A0, L3, L3])/(L3^2*(L3 - R3)^2) + 
         (delta[A0, R3, L3]*phi[A0, R3, L3])/(L3*(L3 - R3)^2) + 
         ((-(A0*(A0 - 6*R3)) + delta[A0, R3, R3])*phi[A0, R3, R3])/
          ((L3 - R3)*R3^2) + Ytau^2*((4*L3 - 2*R3)/(L3^2*R3 - L3*R3^2) + 
           (3*Log[L3]^2)/(L3 - R3)^2 + 
           ((-((L3 - R3)*((L3 - R3)^2*(2*L3 - R3) + A0^2*(2*L3 + R3) - 
                 4*A0*L3*(L3 + 2*R3))*delta[A0, R3, L3]) + delta[A0, L3, R3]*(
                2*L3*(L3 - R3)*(A0 + L3 - R3)*R3 + (2*L3^2 + L3*R3 - R3^2)*
                 delta[A0, R3, L3]))*Log[R3])/(L3*(L3 - R3)^2*R3*
             delta[A0, L3, R3]*delta[A0, R3, L3]) + (2*Log[R3]^2)/
            (L3 - R3)^2 + Log[L3]*(((L3 - R3)*(-A0^3 + L3*(L3 - R3)^2 + 
                 A0^2*(3*L3 + 4*R3) - A0*(3*L3^2 + 2*L3*R3 + 7*R3^2))*
                delta[A0, R3, L3] + delta[A0, L3, R3]*(-((L3 - R3)*R3*
                   (-A0^2 + L3^2 + 2*A0*R3 - R3^2)) + (-L3^2 + A0*(L3 - R3) - 
                   2*L3*R3 + R3^2)*delta[A0, R3, L3]))/((L3 - R3)^2*R3^2*
               delta[A0, L3, R3]*delta[A0, R3, L3]) - (2*Log[q2])/
              (L3 - R3)^2 - (5*Log[R3])/(L3 - R3)^2) + 
           Log[A0]*(((A0^3*L3 - (L3 - R3)^4 + A0*L3*(3*L3^2 - 2*L3*R3 - 
                   R3^2) + A0^2*(-3*L3^2 - 2*L3*R3 + R3^2))*delta[A0, R3, 
                 L3] + delta[A0, L3, R3]*(L3*(-A0^2 + (L3 - R3)^2)*R3 + 
                 (-(A0*L3) + L3^2 - 3*L3*R3 + R3^2)*delta[A0, R3, L3]))/
              (L3*(L3 - R3)*R3^2*delta[A0, L3, R3]*delta[A0, R3, L3]) + 
             Log[L3]/(L3 - R3)^2 - Log[R3]/(L3 - R3)^2) + 
           Log[q2]*((4*L3 - 2*R3)/(L3^2*R3 - L3*R3^2) + (2*Log[R3])/
              (L3 - R3)^2) + ((A0*(A0 - 7*L3) - delta[A0, L3, L3])*
             phi[A0, L3, L3])/(L3^2*(L3 - R3)^2) + 
           ((-((L3 - R3)*(A0^4 - 4*A0*L3^2*(L3 - R3) + (L3 - R3)^4 - 
                 4*A0^3*(L3 + R3) + A0^2*(6*L3^2 + 4*L3*R3 + 6*R3^2))) + 
              (A0^2*(3*L3 - 5*R3) + (3*L3 - 5*R3)*(L3 - R3)^2 + 
                A0*(-6*L3^2 + 4*L3*R3 + 14*R3^2))*delta[A0, L3, R3] - 
              2*(L3 - 2*R3)*delta[A0, L3, R3]^2)*phi[A0, L3, R3])/
            ((L3 - R3)^2*R3^3*delta[A0, L3, R3]) + 
           (((L3 - R3)*(A0 + L3 - R3)^2 + A0*delta[A0, R3, L3])*
             phi[A0, R3, L3])/(L3*(L3 - R3)^2*delta[A0, R3, L3]) + 
           ((A0*(A0 - 6*R3) - delta[A0, R3, R3])*phi[A0, R3, R3])/
            ((L3 - R3)^2*R3^2))) + 
       Xtau^4*((3*L3*(L3 - R3)*R3 + A0*(-2*L3^2 - 5*L3*R3 + R3^2))/
          (L3*(L3 - R3)^3*R3) + Log[L3]*((3*(4*A0*L3 - L3^2 + R3^2))/
            (2*(L3 - R3)^4) + (3*(2*A0*L3 - L3^2 + R3^2)*Log[q2])/
            (L3 - R3)^4) + (3*(-4*A0*L3 + L3^2 - R3^2)*Log[R3])/
          (2*(L3 - R3)^4) + Log[q2]*((6*L3*(L3 - R3)*R3 + 
             A0*(-2*L3^2 - 5*L3*R3 + R3^2))/(L3*(L3 - R3)^3*R3) + 
           (3*(-2*A0*L3 + L3^2 - R3^2)*Log[R3])/(L3 - R3)^4) + 
         Log[A0]*((6*L3*R3*(-L3 + R3) + A0*(2*L3^2 + 5*L3*R3 - R3^2))/
            (L3*(L3 - R3)^3*R3) + (3*(-2*A0*L3 + L3^2 - R3^2)*Log[L3])/
            (L3 - R3)^4 + (3*(2*A0*L3 - L3^2 + R3^2)*Log[R3])/(L3 - R3)^4) + 
         Ytau^2*((-2*L3^2 - 11*L3*R3 + R3^2)/(L3*(L3 - R3)^3*R3) + 
           ((7*A0 - 11*L3 - 3*R3)*Log[L3]^2)/(L3 - R3)^4 - 
           (((L3 - R3)^2*(A0^3 - 3*L3^3 + 11*L3^2*R3 - 5*L3*R3^2 - 3*R3^3 - 
                A0^2*(5*L3 + 3*R3) + A0*(7*L3^2 + 4*L3*R3 + 5*R3^2))*delta[
                A0, R3, L3] + delta[A0, L3, R3]*(2*L3*(L3 - R3)^2*
                 (A0 + L3 - R3)*R3 + (3*L3^3 - A0*(L3 - R3)^2 + 7*L3^2*R3 + 
                  13*L3*R3^2 + R3^3)*delta[A0, R3, L3]))*Log[R3])/
            (2*L3*(L3 - R3)^4*R3*delta[A0, L3, R3]*delta[A0, R3, L3]) - 
           (2*(-2*A0 + L3 + 3*R3)*Log[R3]^2)/(L3 - R3)^4 + 
           Log[q2]*((-2*L3^2 - 5*L3*R3 + R3^2)/(L3*(L3 - R3)^3*R3) - 
             (6*L3*Log[R3])/(L3 - R3)^4) + Log[A0]*
            (((L3 - R3)*(-A0^3 + L3^3 - 3*L3^2*R3 - L3*R3^2 + 3*R3^3 + 
                 3*A0^2*(L3 + R3) - A0*(3*L3^2 + 5*R3^2))*delta[A0, R3, L3] + 
               delta[A0, L3, R3]*(-(L3*(-A0^2 + (L3 - R3)^2)*R3) + 
                 (-L3^2 + A0*(L3 - R3) + 3*L3*R3 + 3*R3^2)*delta[A0, R3, 
                   L3]))/(2*L3*(L3 - R3)^2*R3^2*delta[A0, L3, R3]*delta[A0, 
                R3, L3]) - (3*(A0 - L3 + R3)*Log[L3])/(L3 - R3)^4 + 
             (3*(A0 - L3 + R3)*Log[R3])/(L3 - R3)^4) + 
           Log[L3]*(((L3 - R3)^2*(A0^3 - L3^3 + L3^2*R3 + 9*L3*R3^2 - 
                 9*R3^3 - A0^2*(3*L3 + 5*R3) + A0*(3*L3^2 + 4*L3*R3 + 
                   9*R3^2))*delta[A0, R3, L3] + delta[A0, L3, R3]*
                ((L3 - R3)^2*R3*(-A0^2 + L3^2 + 2*A0*R3 - R3^2) + 
                 (L3^3 - A0*(L3 - R3)^2 + 2*L3^2*R3 + 17*L3*R3^2 + 4*R3^3)*
                  delta[A0, R3, L3]))/(2*(L3 - R3)^4*R3^2*delta[A0, L3, R3]*
               delta[A0, R3, L3]) + (6*L3*Log[q2])/(L3 - R3)^4 + 
             ((-11*A0 + 13*L3 + 9*R3)*Log[R3])/(L3 - R3)^4) + 
           ((A0*(A0 - 7*L3)*(L3 - R3) + (-8*L3 + R3)*delta[A0, L3, L3])*
             phi[A0, L3, L3])/(L3^2*(L3 - R3)^4) + 
           (((L3 - R3)^2*(A0^4 - 4*A0^3*(L3 + R3) - 4*A0*(L3 - R3)^2*
                 (L3 + R3) + (L3 - R3)^2*(L3^2 - 2*L3*R3 - 3*R3^2) + 
                A0^2*(6*L3^2 + 4*L3*R3 + 6*R3^2)) - (L3 - R3)*(3*L3^3 + 
                3*A0^2*(L3 - 3*R3) - 15*L3^2*R3 + 29*L3*R3^2 - 17*R3^3 - 
                6*A0*(L3^2 - 2*L3*R3 - 3*R3^2))*delta[A0, L3, R3] + 
              2*(L3^2 - 5*L3*R3 + 12*R3^2)*delta[A0, L3, R3]^2)*
             phi[A0, L3, R3])/(2*(L3 - R3)^4*R3^3*delta[A0, L3, R3]) - 
           (((L3 - R3)^2*(A0 + L3 - R3)^2 + (4*A0 + 3*L3 - 3*R3)*(L3 - R3)*
               delta[A0, R3, L3] - 6*delta[A0, R3, L3]^2)*phi[A0, R3, L3])/
            (2*L3*(L3 - R3)^4*delta[A0, R3, L3]) + 
           ((A0*(A0 - 6*R3)*(-L3 + R3) + (L3 - 5*R3)*delta[A0, R3, R3])*
             phi[A0, R3, R3])/((L3 - R3)^4*R3^2)))))/cbe^2 + 
    Xtau^4*((2*(2*L3^2 - 27*L3*R3 + R3^2))/(L3*(L3 - R3)^2*R3) + 
      ((-7*L3 - 9*R3)*Log[L3]^2)/(L3 - R3)^3 + 
      ((-6*L3^2 - 44*L3*R3 + 2*R3^2)*Log[R3])/(L3*(L3 - R3)^3) + 
      ((3*L3 + 5*R3)*Log[R3]^2)/(L3 - R3)^3 + 
      Log[L3]*((-4*(L3^2 - 10*L3*R3 - 3*R3^2))/((L3 - R3)^3*R3) + 
        (2*(5*L3 + 7*R3)*Log[q2])/(L3 - R3)^3 + (4*(L3 + R3)*Log[R3])/
         (L3 - R3)^3) + Log[q2]*((2*(2*L3^2 - 15*L3*R3 + R3^2))/
         (L3*(L3 - R3)^2*R3) - (2*(5*L3 + 7*R3)*Log[R3])/(L3 - R3)^3) - 
      (6*PolyLog[2, 1 - L3/R3])/(L3 - R3)^2) + 
    Xtau^2*((-2*L3^2 - 3*L3*R3 + R3^2)/(L3*(L3 - R3)*R3) + 
      ((7*L3 - 9*R3)*Log[L3]^2)/(L3 - R3)^2 + 
      ((16*L3^2 - 21*L3*R3 + R3^2)*Log[R3])/(L3*(L3 - R3)^2) - 
      (9*Log[R3]^2)/(L3 - R3) + Log[L3]*((2*L3^2 - 17*L3*R3 + 19*R3^2)/
         ((L3 - R3)^2*R3) - (2*(8*L3 - 9*R3)*Log[q2])/(L3 - R3)^2 + 
        (2*L3*Log[R3])/(L3 - R3)^2) + 
      Log[q2]*((-2*L3^2 - L3*R3 + R3^2)/(L3*(L3 - R3)*R3) + 
        (2*(8*L3 - 9*R3)*Log[R3])/(L3 - R3)^2) - (6*PolyLog[2, 1 - L3/R3])/
       (L3 - R3)) + Xtau^6*((-2*L3^2 - 11*L3*R3 + R3^2)/(L3*(L3 - R3)^3*R3) + 
      ((-9*L3 - 5*R3)*Log[L3]^2)/(L3 - R3)^4 + 
      ((-10*L3^2 - 3*L3*R3 + R3^2)*Log[R3])/(L3*(L3 - R3)^4) + 
      ((-3*L3 - 5*R3)*Log[R3]^2)/(L3 - R3)^4 + 
      Log[q2]*((-2*L3^2 - 5*L3*R3 + R3^2)/(L3*(L3 - R3)^3*R3) - 
        (6*L3*Log[R3])/(L3 - R3)^4) + 
      Log[L3]*((2*L3^2 + 13*L3*R3 - 3*R3^2)/((L3 - R3)^4*R3) + 
        (6*L3*Log[q2])/(L3 - R3)^4 + (2*(6*L3 + 5*R3)*Log[R3])/(L3 - R3)^4) - 
      (2*PolyLog[2, 1 - L3/R3])/(L3 - R3)^3 + (4*PolyLog[2, (L3 - R3)/L3])/
       (L3 - R3)^3) + 
    ((4*mu2^3*(2*L3 + R3) - L3*R3*(4*L3^2 + 3*L3*R3 + 2*R3^2) - 
        3*mu2^2*(4*L3^2 + 7*L3*R3 + 2*R3^2) + mu2*(4*L3^3 + 17*L3^2*R3 + 
          13*L3*R3^2 + 2*R3^3))/(2*L3*(-L3 + mu2)*(mu2 - R3)*R3) - 
      (2*L3*(L3 - 2*mu2)*Log[L3]^2)/(L3 - mu2)^2 + 
      ((4*mu2^4 - 4*mu2^3*R3 - 2*L3^2*R3^2 + 2*mu2^2*R3*(-4*L3 + R3) + 
         4*L3*mu2*R3*(L3 + R3))*Log[mu2]^2)/((L3 - mu2)^2*(mu2 - R3)^2) - 
      2*Log[q2]^2 + ((-2*mu2*(mu2 - R3)^2*R3 + L3^2*(9*mu2^2 - 4*mu2*R3 + 
           R3^2) + L3*(-13*mu2^3 + 14*mu2^2*R3 - 9*mu2*R3^2 + 2*R3^3))*
        Log[R3])/(2*L3*(L3 - mu2)*(mu2 - R3)^2) + ((2*mu2 - R3)*R3*Log[R3]^2)/
       (mu2 - R3)^2 + Log[q2]*(-(((2*L3 + R3)*(L3 - 2*mu2 + R3))/(L3*R3)) + 
        Log[R3]) + Log[mu2]*((mu2*(-2*L3^2*R3^2*(L3 + R3) - 
           2*mu2^4*(2*L3 + R3) + L3*mu2*R3*(3*L3^2 + 8*L3*R3 - 2*R3^2) + 
           mu2^3*(8*L3^2 + L3*R3 + 4*R3^2) - 2*mu2^2*(2*L3^3 + 4*L3^2*R3 - 
             L3*R3^2 + R3^3)))/(L3*(L3 - mu2)^2*(mu2 - R3)^2*R3) + 
        ((2*L3*mu2^3 - 3*mu2^4 - L3*mu2^2*(L3 - 8*R3) + 2*L3^2*R3^2 - 
           4*L3*mu2*R3*(L3 + R3))*Log[R3])/((L3 - mu2)^2*(mu2 - R3)^2)) + 
      Log[L3]*((4*L3^3*(mu2 - R3) + mu2^2*(15*mu2 - 13*R3)*R3 + 
          L3^2*(-8*mu2^2 + 15*mu2*R3 - 5*R3^2) + 2*L3*mu2*
           (2*mu2^2 - 7*mu2*R3 + 3*R3^2))/(2*(L3 - mu2)^2*(mu2 - R3)*R3) + 
        ((-5*mu2^4 - 2*mu2^3*(L3 - 4*R3) + 2*L3^2*R3^2 - 
           4*L3*mu2*R3*(L3 + R3) + mu2^2*(L3^2 + 8*L3*R3 - 4*R3^2))*Log[mu2])/
         ((L3 - mu2)^2*(mu2 - R3)^2) + 3*Log[q2] + 
        ((2*mu2^4 - 2*mu2^3*R3 - L3^2*R3^2 + mu2^2*R3*(-4*L3 + R3) + 
           2*L3*mu2*R3*(L3 + R3))*Log[R3])/((L3 - mu2)^2*(mu2 - R3)^2)) + 
      (2*(-L3^2 + 2*L3*mu2 + mu2^2)*PolyLog[2, 1 - L3/mu2])/(L3 - mu2)^2 + 
      (4*mu2^2*PolyLog[2, (-L3 + mu2)/mu2])/(L3 - mu2)^2 + 
      (2*(mu2^2 + 2*mu2*R3 - R3^2)*PolyLog[2, 1 - R3/mu2])/(mu2 - R3)^2 + 
      Xtau^2*((4*L3^2 - 8*L3*mu2 + 6*L3*R3 + 4*mu2*R3 - 2*R3^2)/
         (L3^2*R3 - L3*R3^2) + (2*L3*(3*L3^2 - 6*L3*mu2 + 5*mu2^2 + 2*L3*R3 - 
           4*mu2*R3)*Log[L3]^2)/((L3 - mu2)^2*(L3 - R3)^2) + 
        ((L3^3*(5*mu2 - 11*R3) + 2*mu2*(mu2 - R3)*R3^2 + 
           L3^2*(-5*mu2^2 + 16*mu2*R3 + R3^2) + L3*(4*mu2^3 - 13*mu2^2*R3 + 
             mu2*R3^2 + 2*R3^3))*Log[R3])/(L3*(L3 - mu2)*(L3 - R3)^2*
          (mu2 - R3)) + (2*R3*(4*mu2^2 - 6*mu2*R3 + 3*R3^2 + 
           L3*(-2*mu2 + R3))*Log[R3]^2)/((L3 - R3)^2*(mu2 - R3)^2) + 
        Log[L3]*((-4*L3^3*(mu2 - R3) + L3*R3*(5*mu2^2 - 2*mu2*R3 - 7*R3^2) - 
            mu2*R3*(4*mu2^2 - 7*mu2*R3 + R3^2) + L3^2*(4*mu2^2 - 13*mu2*R3 + 
              11*R3^2))/((L3 - mu2)*(L3 - R3)^2*(mu2 - R3)*R3) + 
          (2*(-2*mu2^5 + 2*L3^3*R3^2 - 2*L3^2*mu2*R3*(2*L3 + 3*R3) + 
             mu2^4*(3*L3 + 7*R3) - 2*mu2^3*(2*L3^2 + 5*L3*R3 + 3*R3^2) + 
             mu2^2*(L3^3 + 13*L3^2*R3 + 4*L3*R3^2 + 2*R3^3))*Log[mu2])/
           ((L3 - mu2)^2*(L3 - R3)^2*(mu2 - R3)^2) + 
          ((-8*L3 + 4*mu2 + 2*R3)*Log[q2])/(L3 - R3)^2 + 
          ((-6*L3)/(L3 - R3)^2 + (4*L3*(L3 - 2*mu2))/((L3 - mu2)^2*
              (L3 - R3)) - (4*R3)/(L3 - R3)^2 - (4*(L3 + R3))/(L3 - R3)^2 + 
            (2*R3*(-2*mu2 + R3))/((mu2 - R3)^2*(-L3 + R3)))*Log[R3]) + 
        Log[q2]*((4*L3^2 - 8*L3*mu2 + 4*L3*R3 + 4*mu2*R3 - 2*R3^2)/
           (L3^2*R3 - L3*R3^2) + ((8*L3 - 2*(2*mu2 + R3))*Log[R3])/
           (L3 - R3)^2) + Log[mu2]*((4*mu2*(L3^2*(2*mu2 - R3) + 
             mu2*(mu2 - R3)*R3 + L3*mu2*(-2*mu2 + R3)))/(L3*(L3 - mu2)*
            (L3 - R3)*(mu2 - R3)*R3) + (2*(2*mu2^5 - 2*L3^3*R3^2 + 
             2*L3^2*mu2*R3*(2*L3 + 3*R3) - mu2^4*(3*L3 + 7*R3) + 
             2*mu2^3*(2*L3^2 + 5*L3*R3 + 3*R3^2) - mu2^2*(L3^3 + 13*L3^2*R3 + 
               4*L3*R3^2 + 2*R3^3))*Log[R3])/((L3 - mu2)^2*(L3 - R3)^2*
            (mu2 - R3)^2)) + (4*(-L3 + mu2)*PolyLog[2, (-L3 + mu2)/mu2])/
         (L3 - R3)^2 + (4*(L3 - mu2)*PolyLog[2, 1 - R3/mu2])/(L3 - R3)^2) + 
      Xtau^4*((-2*L3^4*(mu2 - R3) + mu2*R3^2*(2*mu2^2 - 3*mu2*R3 + R3^2) + 
          L3^3*(6*mu2^2 - 36*mu2*R3 + 32*R3^2) - 
          L3*R3*(16*mu2^3 - 5*mu2^2*R3 - 14*mu2*R3^2 + R3^3) - 
          L3^2*(4*mu2^3 - 46*mu2^2*R3 + 31*mu2*R3^2 + 15*R3^3))/
         (L3*(L3 - mu2)*(L3 - R3)^3*(mu2 - R3)*R3) - 
        (2*L3*(5*L3^3 + mu2*(5*mu2 - 2*R3)*R3 + 5*L3^2*(-2*mu2 + R3) + 
           L3*(6*mu2^2 - 10*mu2*R3 + R3^2))*Log[L3]^2)/
         ((L3 - mu2)^2*(L3 - R3)^4) + ((-2*mu2*(mu2 - R3)^2*R3^3 + 
           L3^4*(-27*mu2^2 + 64*mu2*R3 - 33*R3^2) + 
           L3^3*(47*mu2^3 - 162*mu2^2*R3 + 153*mu2*R3^2 - 50*R3^3) + 
           L3*R3*(-12*mu2^4 - 19*mu2^3*R3 + 70*mu2^2*R3^2 - 45*mu2*R3^3 + 
             2*R3^4) + L3^2*(-24*mu2^4 + 118*mu2^3*R3 - 101*mu2^2*R3^2 - 
             26*mu2*R3^3 + 45*R3^4))*Log[R3])/(2*L3*(L3 - mu2)*(L3 - R3)^4*
          (mu2 - R3)^2) - (R3*(L3 + R3)*(8*mu2^2 - 14*mu2*R3 + 7*R3^2 + 
           L3*(-2*mu2 + R3))*Log[R3]^2)/((L3 - R3)^4*(mu2 - R3)^2) + 
        Log[q2]*((-2*L3^3 + L3^2*(4*mu2 - 21*R3) + R3^2*(-2*mu2 + R3) + 
            10*L3*R3*(mu2 + R3))/(L3*(L3 - R3)^3*R3) + 
          ((-13*L3^2 + 12*L3*mu2 - 6*L3*R3 + 7*R3^2)*Log[R3])/(L3 - R3)^4) + 
        Log[L3]*((4*L3^5*(mu2 - R3) + mu2^2*R3^3*(-13*mu2 + 15*R3) + 
            L3^4*(-8*mu2^2 + 77*mu2*R3 - 71*R3^2) + 2*L3*mu2*R3*
             (-18*mu2^3 + 4*mu2^2*R3 + 37*mu2*R3^2 - 25*R3^3) + 
            2*L3^3*(2*mu2^3 - 87*mu2^2*R3 + 81*mu2*R3^2 + 6*R3^3) + 
            L3^2*R3*(145*mu2^3 - 123*mu2^2*R3 - 49*mu2*R3^2 + 27*R3^3))/
           (2*(L3 - mu2)^2*(L3 - R3)^4*(mu2 - R3)*R3) + 
          ((12*L3*mu2^5 + 2*L3^2*R3^2*(-2*L3^2 - 2*L3*R3 + R3^2) - 
             mu2^4*(25*L3^2 + 28*L3*R3 + R3^2) + 6*L3*mu2^3*(3*L3^2 + 10*L3*
                R3 + 3*R3^2) - L3*mu2^2*(3*L3^3 + 44*L3^2*R3 + 41*L3*R3^2 - 4*
                R3^3) + 4*L3*mu2*R3*(2*L3^3 + 7*L3^2*R3 + L3*R3^2 - R3^3))*
            Log[mu2])/((L3 - mu2)^2*(L3 - R3)^4*(mu2 - R3)^2) + 
          ((13*L3^2 - 7*R3^2 + 6*L3*(-2*mu2 + R3))*Log[q2])/(L3 - R3)^4 + 
          ((10*L3)/(L3 - R3)^3 - 2/(L3 - R3)^2 + (22*L3*R3)/(L3 - R3)^4 + 
            (6*R3)/(-L3 + R3)^3 - (2*L3*(L3 - 2*mu2)*(L3 + R3))/
             ((L3 - mu2)^2*(L3 - R3)^3) + ((2*mu2 - R3)*R3*(L3 + R3))/
             ((mu2 - R3)^2*(-L3 + R3)^3) + (4*(L3 + R3)^2)/(L3 - R3)^4)*
           Log[R3]) + Log[mu2]*((-2*mu2*(-(mu2^2*(mu2 - R3)^2*R3^2) + 
             L3*mu2^2*R3*(2*mu2^2 - 5*mu2*R3 + 2*R3^2) + 
             L3^4*(2*mu2^2 - 3*mu2*R3 + 2*R3^2) + L3^3*(-4*mu2^3 + 8*mu2^2*
                R3 - 9*mu2*R3^2 + 2*R3^3) + L3^2*(2*mu2^4 - 5*mu2^3*R3 + 7*
                mu2^2*R3^2 - R3^4)))/(L3*(L3 - mu2)^2*(L3 - R3)^3*
            (mu2 - R3)^2*R3) + ((-12*L3*mu2^5 + 2*L3^2*R3^2*(2*L3^2 + 2*L3*
                R3 - R3^2) + mu2^4*(25*L3^2 + 28*L3*R3 + R3^2) - 
             6*L3*mu2^3*(3*L3^2 + 10*L3*R3 + 3*R3^2) + L3*mu2^2*
              (3*L3^3 + 44*L3^2*R3 + 41*L3*R3^2 - 4*R3^3) - 
             4*L3*mu2*R3*(2*L3^3 + 7*L3^2*R3 + L3*R3^2 - R3^3))*Log[R3])/
           ((L3 - mu2)^2*(L3 - R3)^4*(mu2 - R3)^2)) + 
        ((4*L3^2 + 6*mu2^2 - 2*R3^2 + 4*L3*(-3*mu2 + R3))*
          PolyLog[2, (-L3 + mu2)/mu2])/(L3 - R3)^4 + 
        (2*(-2*L3^2 + 6*L3*mu2 - 3*mu2^2 - 2*L3*R3 + R3^2)*
          PolyLog[2, 1 - R3/mu2])/(L3 - R3)^4))/cbe^2) + 
  (gtau^4*Nc*(Xb*Xtau*(Log[L3]*(4/(-L3 + R3) + (4*Log[q2])/(-L3 + R3) + 
         (4*Q3*Log[Q3])/((-D3 + Q3)*(L3 - R3))) + (4*Log[R3])/(L3 - R3) + 
       (4*Log[q2]*Log[R3])/(L3 - R3) + (4*Q3*Log[Q3]*Log[R3])/
        ((D3 - Q3)*(L3 - R3)) + Log[D3]*
        ((4*D3*Log[L3])/((D3 - Q3)*(L3 - R3)) + (4*D3*Log[R3])/
          ((D3 - Q3)*(-L3 + R3)))) + Xb*Xtau^3*(-8/(L3 - R3)^2 + 
       Log[L3]*((4*(L3 + R3))/(L3 - R3)^3 + (4*(L3 + R3)*Log[q2])/
          (L3 - R3)^3 + (4*Q3*(L3 + R3)*Log[Q3])/((D3 - Q3)*(L3 - R3)^3)) - 
       (4*(L3 + R3)*Log[R3])/(L3 - R3)^3 + Log[q2]*(-8/(L3 - R3)^2 - 
         (4*(L3 + R3)*Log[R3])/(L3 - R3)^3) + 
       Log[D3]*((8*D3)/((D3 - Q3)*(L3 - R3)^2) - (4*D3*(L3 + R3)*Log[L3])/
          ((D3 - Q3)*(L3 - R3)^3) + (4*D3*(L3 + R3)*Log[R3])/
          ((D3 - Q3)*(L3 - R3)^3)) + Log[Q3]*
        ((8*Q3)/((-D3 + Q3)*(L3 - R3)^2) - (4*Q3*(L3 + R3)*Log[R3])/
          ((D3 - Q3)*(L3 - R3)^3))))*ybMSSM[mQ3, mU3, mD3, M3, Mu, TanBeta, 
      Xt, Xb]^2)/cbe^2 + 
  (gtau^2*(48*k^2*(Nc*Xtau^2*(Log[L3]*((-4*L3*R3*Log[q2])/(L3 - R3)^3 + 
           (2*L3*R3*Log[Q3])/(L3 - R3)^3) + Log[Q3]*
          (-((L3 + R3)/(L3 - R3)^2) - (2*L3*R3*Log[R3])/(L3 - R3)^3) + 
         Log[D3]*(-((L3 + R3)/(L3 - R3)^2) + (2*L3*R3*Log[L3])/(L3 - R3)^3 - 
           (2*L3*R3*Log[R3])/(L3 - R3)^3) + Log[q2]*
          ((2*(L3 + R3))/(L3 - R3)^2 + (4*L3*R3*Log[R3])/(L3 - R3)^3)) + 
       Nc*Xb^2*Xtau^2*((-4*L3*R3*Log[L3]*Log[Q3])/((D3 - Q3)*(L3 - R3)^3) + 
         Log[D3]*((-2*(L3 + R3))/((D3 - Q3)*(L3 - R3)^2) + 
           (4*L3*R3*Log[L3])/((D3 - Q3)*(L3 - R3)^3) - (4*L3*R3*Log[R3])/
            ((D3 - Q3)*(L3 - R3)^3)) + Log[Q3]*
          ((2*(L3 + R3))/((D3 - Q3)*(L3 - R3)^2) + (4*L3*R3*Log[R3])/
            ((D3 - Q3)*(L3 - R3)^3))) + Nc*Xb^4*Xtau^2*
        ((-2*(L3 + R3))/((D3 - Q3)^2*(L3 - R3)^2) + 
         Log[L3]*((4*L3*R3)/((D3 - Q3)^2*(L3 - R3)^3) + 
           (2*L3*(D3 + Q3)*R3*Log[Q3])/((D3 - Q3)^3*(L3 - R3)^3)) - 
         (4*L3*R3*Log[R3])/((D3 - Q3)^2*(L3 - R3)^3) + 
         Log[Q3]*(-(((D3 + Q3)*(L3 + R3))/((D3 - Q3)^3*(L3 - R3)^2)) - 
           (2*L3*(D3 + Q3)*R3*Log[R3])/((D3 - Q3)^3*(L3 - R3)^3)) + 
         Log[D3]*(((D3 + Q3)*(L3 + R3))/((D3 - Q3)^3*(L3 - R3)^2) - 
           (2*L3*(D3 + Q3)*R3*Log[L3])/((D3 - Q3)^3*(L3 - R3)^3) + 
           (2*L3*(D3 + Q3)*R3*Log[R3])/((D3 - Q3)^3*(L3 - R3)^3)))) + 
     (48*k^2*Nc*Xtau*(Xb^3*(-8/(D3 - Q3)^2 + Log[q2]*(-8/(D3 - Q3)^2 - 
            (4*(D3 + Q3)*Log[Q3])/(D3 - Q3)^3) + 
          Log[L3]*((8*L3)/((D3 - Q3)^2*(L3 - R3)) + (4*L3*(D3 + Q3)*Log[Q3])/
             ((D3 - Q3)^3*(L3 - R3))) + (8*R3*Log[R3])/((D3 - Q3)^2*
            (-L3 + R3)) + Log[Q3]*((-4*(D3 + Q3))/(D3 - Q3)^3 - 
            (4*(D3 + Q3)*R3*Log[R3])/((D3 - Q3)^3*(L3 - R3))) + 
          Log[D3]*((4*(D3 + Q3))/(D3 - Q3)^3 + (4*L3*(D3 + Q3)*Log[L3])/
             ((D3 - Q3)^3*(-L3 + R3)) + (4*(D3 + Q3)*Log[q2])/(D3 - Q3)^3 + 
            (4*(D3 + Q3)*R3*Log[R3])/((D3 - Q3)^3*(L3 - R3)))) + 
        Xb*((4*L3*Log[L3]*Log[Q3])/((D3 - Q3)*(-L3 + R3)) + 
          (4*Log[q2]*Log[Q3])/(D3 - Q3) + Log[Q3]*(4/(D3 - Q3) + 
            (4*R3*Log[R3])/((D3 - Q3)*(L3 - R3))) + 
          Log[D3]*(4/(-D3 + Q3) + (4*L3*Log[L3])/((D3 - Q3)*(L3 - R3)) + 
            (4*Log[q2])/(-D3 + Q3) + (4*R3*Log[R3])/((D3 - Q3)*
              (-L3 + R3))))))/cbe^2)*ybMSSM[mQ3, mU3, mD3, M3, Mu, TanBeta, 
      Xt, Xb]^4)/(48*k^2)) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
  ];
(* Two loop alphas-alphab corrections, computation for general masses *)
lambda2LPhiHSSAlphaBAlphaSFull = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    Q  = msq2[3,3]*(1+0.02),
    U  = msu2[3,3]*(1-0.02),
    mD  = msd2[3,3]*(1-0.015),
    mu2 = MuInput^2,
    q2 = SCALE^2, (* renormalization/matching scale *)
    g = M3Input^2,
    sgn = M3Input/Abs[M3Input],
    gt = Yu[3,3], (* SM Yukawa coupling *)
    gb = Yd[3,3],
    gtau = Ye[3,3],
    Xt = xt,
    Yt = yt,
    Xb = xb,
    Yb = yb,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    Q2 = SCALE^2,
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
    },
    ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];

    (
       g3^2*ybMSSM[mQ3, mU3, mD3, M3, Mu, TanBeta, Xt, Xb]^4*k^2*((8*(-3*mD^2*Q^2 + 2*g^3*(mD + Q) + 6*g*mD*Q*(mD + Q) - 
            g^2*(2*mD^2 + 9*mD*Q + 2*Q^2)))/((g - mD)*mD*(g - Q)*Q) + 
         (24*(2*g^2 - 2*g*mD + mD^2)*Log[mD])/(g - mD)^2 - 
         (8*(3*g^2 - 2*g*mD + mD^2)*Log[mD]^2)/(g - mD)^2 - 
         (8*(3*g^2 - 2*g*Q + Q^2)*Log[Q]^2)/(g - Q)^2 + 
         ((16*(-3*mD^2*Q^2 + g^3*(mD + Q) + 3*g*mD*Q*(mD + Q) - 
              g^2*(mD^2 + 3*mD*Q + Q^2)))/((g - mD)*mD*(g - Q)*Q) + 
           (16*(2*g^2 - 2*g*mD + mD^2)*Log[mD])/(g - mD)^2)*Log[q2] - 16*Log[q2]^2 + 
         Log[Q]*((24*(2*g^2 - 2*g*Q + Q^2))/(g - Q)^2 + 
           (16*(2*g^2 - 2*g*Q + Q^2)*Log[q2])/(g - Q)^2) + 
         Log[g]*((-8*g^2*(2*g^3*(mD + Q) + 2*g*(-mD + Q)^2*(mD + Q) + 
              g^2*(-4*mD^2 + 2*mD*Q - 4*Q^2) + mD*Q*(mD^2 + Q^2)))/
            ((g - mD)^2*mD*(g - Q)^2*Q) + (16*g^2*Log[mD])/(g - mD)^2 + 
           (16*g^2*Log[Q])/(g - Q)^2 - (16*g^2*(2*g^2 + mD^2 + Q^2 - 2*g*(mD + Q))*
             Log[q2])/((g - mD)^2*(g - Q)^2)) + 
         Xb^2*((-32*g)/(mD*Q) + (32*g^2*(g - mD - Q)*Log[g])/
            ((g - mD)*mD*(g - Q)*Q) + (32*(3*g - 2*mD)*Log[mD])/
            ((g - mD)*(mD - Q)) + (16*(-3*mD + Q)*Log[mD]^2)/(-mD + Q)^2 + 
           (16*(mD - 3*Q)*Log[Q]^2)/(-mD + Q)^2 + 
           ((-32*g)/(mD*Q) - (64*Log[mD])/(-mD + Q))*Log[q2] + 
           Log[Q]*((32*(3*g - 2*Q))/((g - Q)*(-mD + Q)) + (32*(mD + Q)*Log[mD])/
              (-mD + Q)^2 + (64*Log[q2])/(-mD + Q))) - 
         (32*g^2*PolyLog[2, 1 - g/mD])/(g - mD)^2 - (32*g^2*PolyLog[2, 1 - g/Q])/
          (g - Q)^2 + Xb^4*((16*(6*mD*Q + g*(mD + Q)))/(mD*Q*(-mD + Q)^2) - 
           (16*g*(mD + Q)*Log[g])/(mD*Q*(-mD + Q)^2) + 
           (16*(2*g + 7*mD + 3*Q)*Log[mD])/(-mD + Q)^3 + 
           (16*(g*(mD - Q) + 2*mD*(mD + Q))*Log[mD]^2)/(-mD + Q)^4 + 
           (16*(g*(-mD + Q) + 2*Q*(mD + Q))*Log[Q]^2)/(-mD + Q)^4 + 
           ((16*(4*mD*Q + g*(mD + Q)))/(mD*Q*(-mD + Q)^2) + 
             (32*(g + mD + Q)*Log[mD])/(-mD + Q)^3)*Log[q2] + 
           Log[Q]*((-16*(2*g + 3*mD + 7*Q))/(-mD + Q)^3 - (32*(mD + Q)^2*Log[mD])/
              (-mD + Q)^4 - (32*(g + mD + Q)*Log[q2])/(-mD + Q)^3) + 
           (16*(-2*g + mD + Q)*PolyLog[2, 1 - g/mD])/(-mD + Q)^3 + 
           (16*(2*g - mD - Q)*PolyLog[2, 1 - g/Q])/(-mD + Q)^3) + 
         Xb^3*((-256*Sqrt[g]*sgn)/(-mD + Q)^2 - (64*Sqrt[g]*(3*mD + Q)*sgn*Log[mD])/
            (-mD + Q)^3 + (32*Sqrt[g]*(2*g - mD - Q)*sgn*Log[mD]^2)/(mD - Q)^3 + 
           (32*Sqrt[g]*(2*g - mD - Q)*sgn*Log[Q]^2)/(-mD + Q)^3 + 
           Log[g]*((128*g^(3/2)*sgn*Log[mD])/(-mD + Q)^3 - (128*g^(3/2)*sgn*Log[Q])/
              (-mD + Q)^3) + ((-128*Sqrt[g]*sgn)/(-mD + Q)^2 - 
             (64*Sqrt[g]*(mD + Q)*sgn*Log[mD])/(-mD + Q)^3)*Log[q2] + 
           Log[Q]*((64*Sqrt[g]*(mD + 3*Q)*sgn)/(-mD + Q)^3 + 
             (64*Sqrt[g]*(mD + Q)*sgn*Log[q2])/(-mD + Q)^3) + 
           (64*Sqrt[g]*(2*g - mD - Q)*sgn*PolyLog[2, 1 - g/mD])/(mD - Q)^3 + 
           (64*Sqrt[g]*(2*g - mD - Q)*sgn*PolyLog[2, 1 - g/Q])/(-mD + Q)^3) + 
         Xb*((64*Sqrt[g]*sgn*Log[mD])/(-mD + Q) + (64*g^(3/2)*sgn*Log[mD]^2)/
            ((g - mD)*(mD - Q)) + (64*g^(3/2)*sgn*Log[Q]^2)/((g - Q)*(-mD + Q)) - 
           (64*g^(3/2)*sgn*Log[mD]*Log[q2])/((g - mD)*(mD - Q)) + 
           Log[g]*((-64*g^(3/2)*sgn*Log[mD])/((g - mD)*(mD - Q)) - 
             (64*g^(3/2)*sgn*Log[Q])/((g - Q)*(-mD + Q)) + (64*g^(3/2)*sgn*Log[q2])/
              ((g - mD)*(g - Q))) + Log[Q]*((-64*Sqrt[g]*sgn)/(-mD + Q) - 
             (64*g^(3/2)*sgn*Log[q2])/((g - Q)*(-mD + Q))) + 
           (128*g^(3/2)*sgn*PolyLog[2, 1 - g/mD])/((g - mD)*(mD - Q)) + 
           (128*g^(3/2)*sgn*PolyLog[2, 1 - g/Q])/((g - Q)*(-mD + Q))))
      ) /. { delta[x_,y_,z_] :> TDelta[x,y,z], phi[x_,y_,z_] :> TPhi[x,y,z] }
      ];

(* Two loop alphas-alphab corrections, computation for degenerate masses *)
lambda2LPhiHSSAlphaBAlphaSDegenerate = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    q2 = SCALE^2, (* renormalization/matching scale *)
    MS2 = msq2[3,3],
    gb = Yd[3,3],
    gtau = Ye[3,3],
    Xb = xb,
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    Mu = MuInput,
    Q2 = SCALE^2,
    gt = Yu[3,3],
    Xt = xt,
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
  },
  ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
    ];
  g3^2*(ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb])^4*k^2*((-32*Xb)/Sqrt[MS2] + (16*Xb^2)/MS2 + (16*Xb^3)/(3*MS2^(3/2)) -
  (4*Xb^4)/(3*MS2^2) + 32*Log[MS2] + (32*Xb*Log[MS2])/Sqrt[MS2] -
  (32*Xb^2*Log[MS2])/MS2 - (32*Xb^3*Log[MS2])/(3*MS2^(3/2)) - 16*Log[MS2]^2 -
  32*Log[q2] - (32*Xb*Log[q2])/Sqrt[MS2] + (32*Xb^2*Log[q2])/MS2 +
  (32*Xb^3*Log[q2])/(3*MS2^(3/2)) + 32*Log[MS2]*Log[q2] - 16*Log[q2]^2)
];

(* Two loop alphas-alphab corrections, computation for degenerate squark *)
lambda2LPhiHSSAlphaBAlphaSDegenerateSquark = With[{
    k = 1/(4*Pi)^2,
    sbe = Sqrt[TanBeta^2/(1+TanBeta^2)],
    cbe = Sqrt[1/(1+TanBeta^2)],
    Nc = 3, (* number of colors *)
    q2 = SCALE^2, (* renormalization/matching scale *)
    MS2 = msq2[3,3],
    gb = Yd[3,3],
    gtau = Ye[3,3],
    Xb = xb,      
    mQ3 = Sqrt[msq2[3,3]],
    mU3 = Sqrt[msu2[3,3]],
    mD3 = Sqrt[msd2[3,3]],
    M3 = M3Input,
    g = M3Input^2,
    sgn = M3Input/Abs[M3Input],  
    Mu = MuInput,
    Q2 = SCALE^2,
    gt = Yu[3,3],
    Xt = xt,
    Xtildet = xtt,
    Xtildeb = xbb,
    Xtildetau = xtaut,
    cosb = Cos[ArcTan[TanBeta]],
    sinb = Sin[ArcTan[TanBeta]],
    MA = mAInput,
    CF = 4/3
  },
  ybMSSM[mQ3_,mU3_,mD3_,M3_,Mu_,TanBeta_,Xt_,Xb_] := Module[{deltagsb, deltagbyL1, deltagbyL2, deltagbyL3, deltagbyL4},
          deltagsb = - g3^2*CF*k*(1+Log[M3^2/Q2]+TCF[6][mQ3/M3]+TCF[6][mD3/M3]-Xb/M3*TCF[9][mQ3/M3,mD3/M3]);
          deltagbyL1 = - gb^2/cosb^2*k*(3/4*Log[Mu^2/Q2]+3/8*sinb^2*(2*Log[MA^2/Q2]-1)+TCF[6][mQ3/Mu]+1/2*TCF[6][mD3/Mu]);
          deltagbyL2 = - gt^2/sinb^2*k*(1/4*Log[Mu^2/Q2]+1/8*cosb^2*(2*Log[MA^2/Q2]-1)+sinb^2*(Log[MA^2/Q2]-1));
          deltagbyL3 = - gt^2/sinb^2*k*(1/2*TCF[6][mU3/Mu]+(Xt*TanBeta)/(2*Mu)*TCF[9][mQ3/Mu,mU3/Mu] );
          deltagbyL4 = - 1/2*(-gt^2*Nc*k*Xtildet/6*TCF[5][xQU]-gb^2*Nc*k*Xtildeb/6*TCF[5][xQD]-gtau^2*k*Xtildetau/6*TCF[5][xLE]);
          (gb/(1-deltagsb-(deltagbyL1+deltagbyL2+deltagbyL3+deltagbyL4)))
      ];
(-8*g3^2*ybMSSM[mQ3,mU3,mD3,M3,Mu,TanBeta,Xt,Xb]^4*k^2*(6*MS2^3*(3*g^2 - 2*g*MS2 + MS2^2 - 4*g^(3/2)*sgn*Xb)*
    Log[MS2]^2 + g^(3/2)*Log[g]*(6*Sqrt[g]*MS2^2*(2*g + MS2) - 
     24*MS2^2*(g + MS2)*sgn*Xb - 12*Sqrt[g]*(g - 2*MS2)*MS2*Xb^2 + 
     4*(g - 2*MS2)*MS2*sgn*Xb^3 + Sqrt[g]*(2*g - 3*MS2)*Xb^4 - 
     12*MS2^3*(Sqrt[g] - 2*sgn*Xb)*(Log[MS2] - Log[q2])) + 
   MS2*Log[MS2]*(-18*MS2^2*(2*g^2 - 2*g*MS2 + MS2^2) + 
     48*g^(3/2)*MS2^2*sgn*Xb + 12*(g - 2*MS2)*(2*g - MS2)*MS2*Xb^2 + 
     4*Sqrt[g]*MS2^2*sgn*Xb^3 + (-3*g^2 + 6*g*MS2 - 2*MS2^2)*Xb^4 - 
     12*MS2^2*(2*g^2 - 2*g*MS2 + MS2^2 - 2*g^(3/2)*sgn*Xb)*Log[q2]) - 
   (g - MS2)*(3*MS2^2*(4*g^2 - 9*g*MS2 + 3*MS2^2) + 
     24*Sqrt[g]*MS2^2*(-g + MS2)*sgn*Xb - 12*MS2*(g^2 - 3*g*MS2 + MS2^2)*
      Xb^2 + 4*Sqrt[g]*(g - 2*MS2)*MS2*sgn*Xb^3 + (2*g^2 - 4*g*MS2 + MS2^2)*
      Xb^4 + 2*Log[q2]*(3*MS2^2*(2*g^2 - 3*g*MS2 + 3*MS2^2) - 
       12*g^(3/2)*MS2^2*sgn*Xb - 6*(g - 2*MS2)*(g - MS2)*MS2*Xb^2 + 
       2*Sqrt[g]*(g - MS2)*MS2*sgn*Xb^3 + (g - MS2)^2*Xb^4 + 
       3*MS2^3*(-g + MS2)*Log[q2])) + 24*g^(3/2)*MS2^3*(Sqrt[g] - 2*sgn*Xb)*
    PolyLog[2, 1 - g/MS2]))/(3*(g - MS2)^2*MS2^3)
];

(* convert Piecewise to Which statement *)
PiecewiseToWhich[conds_List] := 
    Sequence @@ (Sequence @@@ (Reverse /@ conds));

PiecewiseToWhich[HoldPattern[Piecewise[conds_List, default_]]] :=
    Which[Evaluate@PiecewiseToWhich[conds], True, default];

PiecewiseToWhich[expr_] := expr;

(* shift DR -> OS at O(at*as), from SUSYHD 1.1 *)
lambda2LDRtoOSAtAs = With[{
    gt = Yu[3,3],
    mtilde = SCALE,
    mQ3 = Sqrt[Abs[msq2[3,3]]],
    mU3 = Sqrt[Abs[msu2[3,3]]],
    Xt = xt,
    M3 = M3Input
    },
    (
    g3^2 gt^4 With[{\[Delta] = 
    8.5 10^-2, \[Delta]2 = 1 10^-2}, 
  Piecewise[{{shiftDROS\[Alpha]t\[Alpha]Sdegenerate[mtilde, mQ3, Xt], 
     Abs[mQ3/mU3 - 1] < \[Delta]2 && Abs[M3/mU3 - 1] < \[Delta]2 && 
      Abs[M3/mQ3 - 1] < \[Delta]2}, {Re[
        shiftDROS\[Alpha]t\[Alpha]Sgeneral[1. mtilde, 
         mU3 (1 + \[Delta]), mU3, Xt , M3]] (mQ3 - 
         mU3 (1 - \[Delta]))/(2 mU3 \[Delta]) + 
      Re[shiftDROS\[Alpha]t\[Alpha]Sgeneral[1. mtilde, 
         mU3 (1 - \[Delta]), mU3, Xt , M3]] (1 - (
         mQ3 - mU3 (1 - \[Delta]))/(2 mU3 \[Delta])), 
     Abs[mQ3/mU3 - 1] < \[Delta]},
    {Re[shiftDROS\[Alpha]t\[Alpha]Sgeneral[1. mtilde, mQ3, mU3, Xt , 
         mQ3 (1 + \[Delta])]] (M3 - mQ3 (1 - \[Delta]))/(
       2 mQ3 \[Delta]) + 
      Re[shiftDROS\[Alpha]t\[Alpha]Sgeneral[1. mtilde, mQ3, mU3, Xt , 
         mQ3 (1 - \[Delta])]] (1 - (M3 - mQ3 (1 - \[Delta]))/(
         2 mQ3 \[Delta])), Abs[mQ3/M3 - 1] < \[Delta]},
    {Re[shiftDROS\[Alpha]t\[Alpha]Sgeneral[1. mtilde, mQ3, mU3, Xt , 
         mU3 (1 + \[Delta])]] (M3 - mU3 (1 - \[Delta]))/(
       2 mU3 \[Delta]) + 
      Re[shiftDROS\[Alpha]t\[Alpha]Sgeneral[1. mtilde, mQ3, mU3, Xt , 
         mU3 (1 - \[Delta])]] (1 - (M3 - mU3 (1 - \[Delta]))/(
         2 mU3 \[Delta])), Abs[mU3/M3 - 1] < \[Delta]}}, 
   Re[shiftDROS\[Alpha]t\[Alpha]Sgeneral[mtilde, 1. mQ3, mU3, Xt , 
     M3]]]]
    ) /. a:Piecewise[__] :> PiecewiseToWhich[a]
    ];

(* shift DR -> OS at O(at*as), general expression, from SUSYHD 1.1 *)
shiftDROS\[Alpha]t\[Alpha]Sgeneral[mtilde_, mQ3_, mU3_, QXt_, M3_] := 
 Re@ Block[{x1 = (mQ3)^2/M3^2, x2 = (mU3)^2/M3^2, xt = QXt/M3, 
            xs = mtilde^2/M3^2}, 
    (
   1/(64 \[Pi]^4) (-((
       2 ((-1 + x1)^2 ComplexLog[1 - x1] + 
          x1 (3 + x1 - (-2 + x1) Log[x1] + 2 Log[xs/x1])))/x1^2) - (
      2 ((-1 + x2)^2 ComplexLog[1 - x2] + 
         x2 (3 + x2 - (-2 + x2) Log[x2] + 2 Log[xs/x2])))/x2^2 + 
      1/(x1^2 x2^2)
        xt (1/(x1 - 
            x2)^3 x1^2 x2^2 (2 (x1 - x2) xt^2 + (2 x1^2 + 
               x2 (2 x2 - xt^2) - x1 (4 x2 + xt^2)) Log[x1/x2]) (-((
             4 xt^2 (x1 (-1 + x2) Log[x1] - (-1 + x1) x2 Log[
                  x2]))/((-1 + x1) (x1 - x2) (-1 + x2))) + 
            1/(x1 - x2)
              xt (-7 x1 + x1/(-1 + x1) + x1/(-1 + x2) + 7 x2 + x2/(
               1 - x1) + x2/(1 - x2) - (4 (-1 + x1)^2 ComplexLog[1 - x1])/
               x1 + ((-x1 + x2) Log[x1])/(-1 + x1)^2 - 
               8 ComplexLog[1 - x2] + (4 ComplexLog[1 - x2])/x2 + 
               4 x2 ComplexLog[1 - x2] - (x1 Log[x2])/(-1 + x2)^2 + (
               x2 Log[x2])/(-1 + x2)^2 + 9 x1 Log[x1/xs] - 
               x2 Log[x1/xs] + x1 Log[x2/xs] - 9 x2 Log[x2/xs] + 
               4 x1 Log[xs] - 4 x2 Log[xs]) - (
            4 ((-1 + x1) x2 ComplexLog[1 - x1] + 
               x1 ((-1 + x2) ComplexLog[1 - x2] - 2 x2 (2 + Log[xs]))))/(
            x1 x2) + (
            xt ((-1 + x1)^2 ComplexLog[1 - x1] + 
               x1 (3 + x1 - (-2 + x1) Log[x1] + 2 Log[xs/x1])))/
            x1^2 + (
            xt ((-1 + x2)^2 ComplexLog[1 - x2] + 
               x2 (3 + x2 - (-2 + x2) Log[x2] + 2 Log[xs/x2])))/
            x2^2) + 
         xt (1/(x1 - x2)^2 2 (-2 x1 + 
               2 x2 + (x1 + x2) Log[x1/x2]) ((-1 + x1)^2 x2^2 ComplexLog[
                 1 - x1] - 
               x1 ((-2 + x1) x2^2 Log[x1] + 
                  x1 (-1 + x2)^2 ComplexLog[1 - x2] + 
                  x2 (3 x1 - 3 x2 - x1 (-2 + x2) Log[x2] - 
                    2 x2 Log[xs/x1] + 2 x1 Log[xs/x2]))) + 
            x1^2 xt (1/(x1^2 (x1 - x2)^3)
                 2 xt (-x1 + x2 + 
                  x1 Log[x1/x2]) ((-1 + x1)^2 x2^2 ComplexLog[1 - x1] - 
                  x1 ((-2 + x1) x2^2 Log[x1] + 
                    x1 (-1 + x2)^2 ComplexLog[1 - x2] + 
                    x2 (3 x1 - 3 x2 - x1 (-2 + x2) Log[x2] - 
                    2 x2 Log[xs/x1] + 2 x1 Log[xs/x2]))) - 
               1/((-1 + x1/x2)^4 x2) (-2 x1 + 
                  2 x2 + (x1 + x2) Log[x1/
                    x2]) ((-1 + x1/
                    x2) (-((
                    4 xt^2 (x1 (-1 + x2) Log[x1] - (-1 + x1) x2 Log[
                    x2]))/((-1 + x1) (x1 - x2) (-1 + x2))) + 
                    1/(x1 - x2)
                     xt (-7 x1 + x1/(-1 + x1) + x1/(-1 + x2) + 7 x2 + 
                    x2/(1 - x1) + x2/(1 - x2) - (
                    4 (-1 + x1)^2 ComplexLog[1 - x1])/
                    x1 + ((-x1 + x2) Log[x1])/(-1 + x1)^2 - 
                    8 ComplexLog[1 - x2] + (4 ComplexLog[1 - x2])/x2 + 
                    4 x2 ComplexLog[1 - x2] - (x1 Log[x2])/(-1 + x2)^2 + (
                    x2 Log[x2])/(-1 + x2)^2 + 9 x1 Log[x1/xs] - 
                    x2 Log[x1/xs] + x1 Log[x2/xs] - 9 x2 Log[x2/xs] + 
                    4 x1 Log[xs] - 4 x2 Log[xs]) - (
                    4 ((-1 + x1) x2 ComplexLog[1 - x1] + 
                    x1 ((-1 + x2) ComplexLog[1 - x2] - 
                    2 x2 (2 + Log[xs]))))/(x1 x2) - (
                    xt ((-1 + x1)^2 ComplexLog[1 - x1] + 
                    x1 (3 + x1 - (-2 + x1) Log[x1] + 2 Log[xs/x1])))/
                    x1^2 + (
                    3 xt ((-1 + x2)^2 ComplexLog[1 - x2] + 
                    x2 (3 + x2 - (-2 + x2) Log[x2] + 2 Log[xs/x2])))/
                    x2^2) + 
                  1/(x1 x2^3)
                    6 xt ((-1 + x1)^2 x2^2 ComplexLog[1 - x1] - 
                    x1 ((-2 + x1) x2^2 Log[x1] + 
                    x1 (-1 + x2)^2 ComplexLog[1 - x2] + 
                    x2 (3 x1 - 3 x2 - x1 (-2 + x2) Log[x2] - 
                    2 x2 Log[xs/x1] + 2 x1 Log[xs/x2]))))))))
    )
];

(* shift DR -> OS at O(at*as), degenerate masses, from SUSYHD 1.1 *)
shiftDROS\[Alpha]t\[Alpha]Sdegenerate[mtilde_, M3_, Xt_] := 
 Re@ Block[{xt = Xt/M3, xs = mtilde^2/M3^2},
           ((2 + xt)^2 (-6 + 18 xt - 9 xt^2 + xt^3) + (-12 + 24 xt - 
           6 xt^2 - 4 xt^3 + xt^4) Log[xs])/(96 \[Pi]^4)];

(* shift DR -> OS at O(at*at), general expression, from SUSYHD 1.1 *)
lambda2LDRtoOSAtAt = With[{
    gt = Yu[3,3],
    Q = SCALE,
    mst = Sqrt[Sqrt[msq2[3,3] msu2[3,3]]],
    \[Mu] = MuInput,
    tan\[Beta] = TanBeta,
    Xt = xt,
    \[Delta] = 10^-5
    },
    (
    gt^6 Block[{\[Beta] = ArcTan[tan\[Beta]], xt = Xt/mst, 
    yt = Xt/mst + (2 \[Mu])/(mst Sin[2 ArcTan[tan\[Beta]]]), 
    x\[Mu] = (\[Mu])^2/mst^2, xq = Q^2/mst^2},
   Sin[\[Beta]]^-2 (1/(
       1024 \[Pi]^4) (-36 - 36 xt^2 + 6 xt^4 + 72 x\[Mu] - 
         72 xt^2 x\[Mu] + 12 xt^4 x\[Mu] - 36 x\[Mu]^2 + 
         108 xt^2 x\[Mu]^2 - 18 xt^4 x\[Mu]^2) log[1 - x\[Mu]] + 
      1/(1024 \[Pi]^4) (108 - 54 xt^2 + 9 xt^4 - 108 x\[Mu] + 
         180 xt^2 x\[Mu] - 30 xt^4 x\[Mu] - 
         6 xt^2 (-6 + xt^2) x\[Mu] f2[Sqrt[x\[Mu]]] + 
         3 (-24 (-1 + x\[Mu]) + 12 xt^2 (3 + 2 x\[Mu]) - 
            2 xt^4 (3 + 2 x\[Mu])) Log[xq] + 
         Cos[\[Beta]]^2 (36 - 48 (-6 + Sqrt[3] \[Pi]) xt yt + 
            8 (-6 + Sqrt[3] \[Pi]) xt^3 yt + 72 yt^2 - 
            12 Sqrt[3] \[Pi] yt^2 + 
            xt^4 (9 - 2 (-6 + Sqrt[3] \[Pi]) yt^2) + 
            6 xt^2 (-9 + 2 (-6 + Sqrt[3] \[Pi]) yt^2) + 
            6 (24 xt yt - 4 xt^3 yt + 6 (1 + yt^2) - 
               6 xt^2 (2 + yt^2) + xt^4 (2 + yt^2)) Log[xq]) + 
         36 x\[Mu]^2 Log[x\[Mu]] - 108 xt^2 x\[Mu]^2 Log[x\[Mu]] + 
         18 xt^4 x\[Mu]^2 Log[x\[Mu]] + 
         6 xt^2 (30 - 10 xt^2 + xt^4) (2 + Log[xq]) Sin[\[Beta]]^2))]
    ) /. {
        a:log[_] :> FiniteLog[a],
        f2[\[Mu]_] :> Piecewise[{{1/2, Abs[\[Mu]^2 - 1] < \[Delta]}}, 
                                1/(1 - \[Mu]^2) (1 + \[Mu]^2/(1 - \[Mu]^2) Log[\[Mu]^2])]
    } /. Piecewise[{{val_, cond_}}, default_] :> If[cond, val, default]
    ];

replaceThresholdLoopFunctions = {
    ftilde1 -> TCF[1],
    ftilde2 -> TCF[2],
    ftilde3 -> TCF[3],
    ftilde4 -> TCF[4],
    ftilde5 -> TCF[5],
    ftilde6 -> TCF[6],
    ftilde7 -> TCF[7],
    ftilde8 -> TCF[8],
    ftilde9 -> TCF[9]
};

(* 2-loop shift O(at*as) when re-parametrizing 1-loop contribution in terms of yt^MSSM *)
lambda2LAtAsYtShift = With[{
    gt = Yu[3,3],
    Q2 = SCALE^2,
    MQ32 = msq2[3,3],
    MU32 = msu2[3,3],
    mu2 = MuInput^2,
    M32 = M3Input^2,
    Xt = xt,
    cosbeta = Cos[ArcTan[TanBeta]],
    sinbeta = Sin[ArcTan[TanBeta]]
    },
    computeStopShift[] := Module[{oneLoopStop,deltagtrengs,oneLoopStopShift,oneloopStopShiftExp,deltagtphigs},
	oneLoopStop = gt^4/(4*Pi)^2*(3*Log[MQ3^2/Q2]+3*Log[MU3^2/Q2]+6*Xtildet*(ftilde1[MQ3/MU3]-Xtildet/12*ftilde2[MQ3/MU3])) /. {Xtildet-> Xt^2/(MQ3*MU3)};
	deltagtphigs = 1/(4*Pi)^2*(-4/3)*g3^2*(Log[M3^2/Q2]+ftilde6[MQ3/M3]+ftilde6[MU3/M3]-Xt/M3*ftilde9[MQ3/M3,MU3/M3]);
	deltagtrengs = -4/3*g3^2*1/(4*Pi)^2;
	oneLoopStopShift = oneLoopStop /. {gt -> gt*(1+deltagtphigs+deltagtrengs)} /. {MU3->Sqrt[MU32], MQ3-> Sqrt[MQ32], M3-> Sqrt[M32]};
	oneloopStopShiftExp = g3^2*Collect[Normal[SeriesCoefficient[Series[oneLoopStopShift,{g3,0,2}],2]], gt, (Collect[#,g3, Collect[#,Pi,Collect[#,Log[_],Simplify]&]&])&];
	oneloopStopShiftExp
    ];
    -computeStopShift[] /. replaceThresholdLoopFunctions
    ];

(* 2-loop shift O(at^2) when re-parametrizing 1-loop contribution in terms of yt^MSSM *)
lambda2LAtAtYtShift = With[{
    gt = Yu[3,3],
    Q2 = SCALE^2,
    MQ32 = msq2[3,3],
    MU32 = msu2[3,3],
    mu2 = MuInput^2,
    M32 = M3Input^2,
    ma2 = mAInput^2,
    Xt = xt,
    cosbeta = Cos[ArcTan[TanBeta]],
    sinbeta = Sin[ArcTan[TanBeta]],
    k = 1/(4 Pi)^2
    },
    computegtStopShift[] := Module[{oneLoopStop,oneLoopStopShift,oneloopStopShiftExp,deltagtphigt},
	oneLoopStop = gt^4*k*(3*Log[MQ3^2/Q2]+3*Log[MU3^2/Q2]+6*Xtildet*(ftilde1[MQ3/MU3]-Xtildet/12*ftilde2[MQ3/MU3])) /. {Xtildet-> Xt^2/(MQ3*MU3)};
	deltagtphigt = k*(-1)*gt^2*( 3/(4*(sinbeta)^2)*Log[mu2/Q2]+3/8*(cosbeta/sinbeta)^2*(2*Log[ma2/Q2]-1)-Xtildet/4*ftilde5[MQ3/MU3] 
										 +1/(sinbeta)^2*ftilde6[MQ3/mu]+1/(2*(sinbeta)^2)*ftilde6[MU3/mu]) /. {Xtildet-> Xt^2/(MQ3*MU3)};
	oneLoopStopShift = oneLoopStop /. {gt -> gt*(1+deltagtphigt)} /. {MU3->Sqrt[MU32], MQ3-> Sqrt[MQ32], M3-> Sqrt[M32],mu->Sqrt[mu2]};
	oneloopStopShiftExp = gt^6*Collect[Normal[SeriesCoefficient[Series[oneLoopStopShift,{gt,0,6}],6]], gt, (Collect[#,gt, Collect[#,Pi,Collect[#,Log[_],Simplify]&]&])&];
	oneloopStopShiftExp
    ];
    -computegtStopShift[] /. replaceThresholdLoopFunctions
    ];

(* 3-loop shift O(at*as^2) when re-parametrizing 1- and 2-loop contribution 
   in terms of yt^MSSM and g3^MSSM
   [degenerate mass case]
 *)
lambda3LAtAsAsYtShiftDeg = With[{
    MR = SCALE,
    MS = Sqrt[msq2[3,3]],
    Xt = xt
    },
    1/2 g3^4 Yu[3, 3]^4/(4 Pi)^6 (
    -2*((1/2 + 2*(-4/3 + (4*Xt)/(3*MS) - (4*Log[MS^2/MR^2])/3) - 
    2*Log[MS^2/MR^2])*((-64*Xt)/MS - (32*Xt^2)/MS^2 + (224*Xt^3)/(3*MS^3) + 
    (8*Xt^4)/(3*MS^4) - (16*Xt^5)/(3*MS^5) + 
    ((128*Xt)/MS - (128*Xt^2)/MS^2 - (64*Xt^3)/(3*MS^3) + (16*Xt^4)/(3*MS^4))*
     Log[MS^2/MR^2] - 96*Log[MS^2/MR^2]^2) + 
  ((12*Xt^2)/MS^2 - Xt^4/MS^4 + 12*Log[MS^2/MR^2])*
   (3*(-4/3 + (4*Xt)/(3*MS) - (4*Log[MS^2/MR^2])/3)^2 + 
    2*(-2075/54 + (356*Xt)/(27*MS) + (16*Xt^2)/(9*MS^2) + 
      (-670/27 - (208*Xt)/(27*MS))*Log[MS^2/MR^2] + (2*Log[MS^2/MR^2]^2)/9)))
    )
    ];

(* 3-loop shift O(at*as^2) when re-parametrizing 1- and 2-loop contribution
   in terms of yt^MSSM
   [general mass case]
 *)
lambda3LAtAsAsYtShift = With[{
    MR = SCALE,
    mQ3 = Sqrt[msq2[3,3] (1 + 0.02)],
    mU3 = Sqrt[msu2[3,3] (1 - 0.02)],
    m3  = M3Input,
    msq = (msq2[1,1] msu2[1,1] msd2[1,1]
           msq2[2,2] msu2[2,2] msd2[2,2])^(1/12),
    Xt = xt
    },
    1/2 g3^4 Yu[3, 3]^4/(4 Pi)^6 (
    -2*((1/2 - Log[m3^2/MR^2] + (-Log[mQ3^2/MR^2] - 10*Log[msq^2/MR^2] - 
      Log[mU3^2/MR^2])/12 + 2*((2*m3^4)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) - 
      (2*mQ3^2*mU3^2)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) + 
      ((-2*(2*m3^8 - 2*m3^6*mQ3^2 + m3^4*mQ3^4 - 2*m3^6*mU3^2 + m3^4*mU3^4))/
         (3*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) + (8*m3^3*Xt)/
         (3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)))*Log[m3^2/MR^2] + 
      ((4*m3^2*mQ3^2)/(3*(m3^2 - mQ3^2)^2) - (2*mQ3^4)/(3*(m3^2 - mQ3^2)^2) - 
        (8*m3*mQ3^2*Xt)/(3*(m3^2 - mQ3^2)*(mQ3^2 - mU3^2)))*Log[mQ3^2/MR^2] + 
      ((4*m3^2*mU3^2)/(3*(m3^2 - mU3^2)^2) - (2*mU3^4)/(3*(m3^2 - mU3^2)^2) - 
        (8*m3*mU3^2*Xt)/(3*(m3^2 - mU3^2)*(-mQ3^2 + mU3^2)))*
       Log[mU3^2/MR^2]))*
   ((16*(-2*m3^6*mQ3^2 + 2*m3^4*mQ3^4 - 2*m3^6*mU3^2 + 9*m3^4*mQ3^2*mU3^2 - 
       6*m3^2*mQ3^4*mU3^2 + 2*m3^4*mU3^4 - 6*m3^2*mQ3^2*mU3^4 + 
       3*mQ3^4*mU3^4))/(mQ3^2*(-m3^2 + mQ3^2)*mU3^2*(m3^2 - mU3^2)) - 
    (64*m3^2*Xt^2)/(mQ3^2*mU3^2) - (512*m3*Xt^3)/(mQ3^2 - mU3^2)^2 + 
    (32*(-(m3^6*mQ3^2) + m3^4*mQ3^4 - m3^6*mU3^2 - 5*m3^4*mQ3^2*mU3^2 + 
       5*m3^2*mQ3^4*mU3^2 + m3^4*mU3^4 + 5*m3^2*mQ3^2*mU3^4 - 5*mQ3^4*mU3^4)*
      Xt^4)/(mQ3^2*(-m3^2 + mQ3^2)*mU3^2*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^2) + 
    ((-16*(3*m3^4 - 4*m3^2*mQ3^2 + 2*mQ3^4))/(m3^2 - mQ3^2)^2 + 
      (64*(2*m3^3 - m3*mQ3^2)*Xt)/((m3^2 - mQ3^2)*(mQ3^2 - mU3^2)) - 
      (32*(4*mQ3^6 - 2*mQ3^4*mU3^2 + m3^4*(3*mQ3^2 - mU3^2) + 
         m3^2*(-8*mQ3^4 + 4*mQ3^2*mU3^2))*Xt^2)/((m3^2 - mQ3^2)^2*
        (mQ3^2 - mU3^2)^2) + (64*m3*(2*m3^4 - mQ3^4 + 3*mQ3^2*mU3^2 - 
         m3^2*(3*mQ3^2 + mU3^2))*Xt^3)/((m3^2 - mQ3^2)*(mQ3^2 - mU3^2)^3) + 
      (16*(5*mQ3^8 + 8*m3^4*mQ3^2*mU3^2 + 4*mQ3^6*mU3^2 - mQ3^4*mU3^4 + 
         2*m3^6*(mQ3^2 - mU3^2) - 2*m3^2*(4*mQ3^6 + 5*mQ3^4*mU3^2 - 
           mQ3^2*mU3^4))*Xt^4)/((m3^2 - mQ3^2)^2*(mQ3^2 - mU3^2)^4) + 
      (64*m3*mQ3^2*(mQ3^2 + mU3^2)*Xt^5)/((m3^2 - mQ3^2)*(mQ3^2 - mU3^2)^4))*
     Log[mQ3^2/MR^2]^2 + 
    ((16*(7*m3^6 - 2*mQ3^2*mU3^4 - m3^4*(6*mQ3^2 + 7*mU3^2) + 
         m3^2*(5*mQ3^2*mU3^2 + 3*mU3^4)))/((m3^2 - mQ3^2)*(m3^2 - mU3^2)^2) + 
      (128*m3*Xt)/(mQ3^2 - mU3^2) - 
      (32*(7*m3^4 + 3*mQ3^2*mU3^2 - 2*m3^2*(3*mQ3^2 + 2*mU3^2))*Xt^2)/
       ((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)) - 
      (128*m3*(mQ3^2 + 3*mU3^2)*Xt^3)/(mQ3^2 - mU3^2)^3 + 
      (16*(4*m3^8 - 3*mQ3^2*mU3^4*(mQ3^2 + 5*mU3^2) + 
         m3^6*(3*mQ3^2 + 7*mU3^2) - m3^4*(6*mQ3^4 + 15*mQ3^2*mU3^2 + 
           29*mU3^4) + m3^2*(7*mQ3^4*mU3^2 + 31*mQ3^2*mU3^4 + 16*mU3^6))*
        Xt^4)/((m3^2 - mQ3^2)*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)^3) - 
      (128*m3*mU3^2*Xt^5)/((m3^2 - mU3^2)*(-mQ3^2 + mU3^2)^3))*
     Log[mU3^2/MR^2] + ((-16*(3*m3^4 - 4*m3^2*mU3^2 + 2*mU3^4))/
       (m3^2 - mU3^2)^2 + (64*(2*m3^3 - m3*mU3^2)*Xt)/
       ((m3^2 - mU3^2)*(-mQ3^2 + mU3^2)) + 
      (32*(m3^4*(mQ3^2 - 3*mU3^2) + 2*mU3^4*(mQ3^2 - 2*mU3^2) + 
         m3^2*(-4*mQ3^2*mU3^2 + 8*mU3^4))*Xt^2)/((m3^2 - mU3^2)^2*
        (mQ3^2 - mU3^2)^2) + (64*m3*(2*m3^4 + 3*mQ3^2*mU3^2 - mU3^4 - 
         m3^2*(mQ3^2 + 3*mU3^2))*Xt^3)/((m3^2 - mU3^2)*(-mQ3^2 + mU3^2)^3) + 
      (16*(8*m3^4*mQ3^2*mU3^2 - mQ3^4*mU3^4 + 4*mQ3^2*mU3^6 + 5*mU3^8 - 
         2*m3^6*(mQ3^2 - mU3^2) + 2*m3^2*mU3^2*(mQ3^4 - 5*mQ3^2*mU3^2 - 
           4*mU3^4))*Xt^4)/((m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)^4) + 
      (64*m3*mU3^2*(mQ3^2 + mU3^2)*Xt^5)/((m3^2 - mU3^2)*(mQ3^2 - mU3^2)^4))*
     Log[mU3^2/MR^2]^2 - 4*((2*m3^4)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) - 
      (2*mQ3^2*mU3^2)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) + 
      ((-2*(2*m3^8 - 2*m3^6*mQ3^2 + m3^4*mQ3^4 - 2*m3^6*mU3^2 + m3^4*mU3^4))/
         (3*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) + (8*m3^3*Xt)/
         (3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)))*Log[m3^2/MR^2] + 
      ((4*m3^2*mQ3^2)/(3*(m3^2 - mQ3^2)^2) - (2*mQ3^4)/(3*(m3^2 - mQ3^2)^2) - 
        (8*m3*mQ3^2*Xt)/(3*(m3^2 - mQ3^2)*(mQ3^2 - mU3^2)))*Log[mQ3^2/MR^2] + 
      ((4*m3^2*mU3^2)/(3*(m3^2 - mU3^2)^2) - (2*mU3^4)/(3*(m3^2 - mU3^2)^2) - 
        (8*m3*mU3^2*Xt)/(3*(m3^2 - mU3^2)*(-mQ3^2 + mU3^2)))*Log[mU3^2/MR^2])*
     ((12*Xt^4)/(mQ3^2 - mU3^2)^2 + (6 + (12*Xt^2)/(mQ3^2 - mU3^2) + 
        ((-6*mQ3^2)/(mQ3^2 - mU3^2)^3 - (6*mU3^2)/(mQ3^2 - mU3^2)^3)*Xt^4)*
       Log[mQ3^2/MR^2] + (6 - (12*Xt^2)/(mQ3^2 - mU3^2) + 
        ((6*mQ3^2)/(mQ3^2 - mU3^2)^3 + (6*mU3^2)/(mQ3^2 - mU3^2)^3)*Xt^4)*
       Log[mU3^2/MR^2]) + Log[m3^2/MR^2]*
     ((-16*m3^4*(2*m3^6*(mQ3^2 + mU3^2) + 2*m3^2*(mQ3^2 - mU3^2)^2*
          (mQ3^2 + mU3^2) + m3^4*(-4*mQ3^4 + 2*mQ3^2*mU3^2 - 4*mU3^4) + 
         mQ3^2*mU3^2*(mQ3^4 + mU3^4)))/(mQ3^2*(m3^2 - mQ3^2)^2*mU3^2*
        (m3^2 - mU3^2)^2) + (64*m3^4*(-m3^2 + mQ3^2 + mU3^2)*Xt^2)/
       (mQ3^2*(-m3^2 + mQ3^2)*mU3^2*(m3^2 - mU3^2)) - 
      (32*m3^2*(m3^8*(mQ3^2 + mU3^2) + mQ3^4*mU3^4*(mQ3^2 + mU3^2) + 
         m3^4*(mQ3^2 + mU3^2)^3 - 2*m3^6*(mQ3^4 + mQ3^2*mU3^2 + mU3^4) - 
         m3^2*mQ3^2*mU3^2*(mQ3^4 + 4*mQ3^2*mU3^2 + mU3^4))*Xt^4)/
       (mQ3^2*(m3^2 - mQ3^2)^2*mU3^2*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)^2) + 
      (128*m3^3*Xt^5)/((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^2) + 
      ((16*m3^4*(2*m3^2 - mQ3^2 - mU3^2)*(mQ3^2 - mU3^2))/
         ((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) + 
        (64*m3^3*(-2*m3^2 + mQ3^2 + mU3^2)*Xt)/((m3^2 - mQ3^2)*(m3^2 - mU3^2)*
          (mQ3^2 - mU3^2)) - (32*m3^4*(2*m3^4 + mQ3^4 + mU3^4 - 
           2*m3^2*(mQ3^2 + mU3^2))*Xt^2)/((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2*
          (mQ3^2 - mU3^2)) - (128*m3^3*(2*m3^4 - mQ3^4 + 4*mQ3^2*mU3^2 - 
           mU3^4 - 2*m3^2*(mQ3^2 + mU3^2))*Xt^3)/((m3^2 - mQ3^2)*
          (m3^2 - mU3^2)*(mQ3^2 - mU3^2)^3) + 
        (16*m3^4*(mQ3^2 + mU3^2)*(2*m3^4 + mQ3^4 + mU3^4 - 
           2*m3^2*(mQ3^2 + mU3^2))*Xt^4)/((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2*
          (mQ3^2 - mU3^2)^3) - (64*m3^3*(mQ3^2 + mU3^2)*Xt^5)/
         ((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^3))*Log[mQ3^2/MR^2] + 
      ((-16*m3^4*(2*m3^2 - mQ3^2 - mU3^2)*(mQ3^2 - mU3^2))/
         ((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) + 
        (64*m3^3*(2*m3^2 - mQ3^2 - mU3^2)*Xt)/((m3^2 - mQ3^2)*(m3^2 - mU3^2)*
          (mQ3^2 - mU3^2)) + (32*m3^4*(2*m3^4 + mQ3^4 + mU3^4 - 
           2*m3^2*(mQ3^2 + mU3^2))*Xt^2)/((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2*
          (mQ3^2 - mU3^2)) + (128*m3^3*(2*m3^4 - mQ3^4 + 4*mQ3^2*mU3^2 - 
           mU3^4 - 2*m3^2*(mQ3^2 + mU3^2))*Xt^3)/((m3^2 - mQ3^2)*
          (m3^2 - mU3^2)*(mQ3^2 - mU3^2)^3) - 
        (16*m3^4*(mQ3^2 + mU3^2)*(2*m3^4 + mQ3^4 + mU3^4 - 
           2*m3^2*(mQ3^2 + mU3^2))*Xt^4)/((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2*
          (mQ3^2 - mU3^2)^3) + (64*m3^3*(mQ3^2 + mU3^2)*Xt^5)/
         ((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^3))*
       Log[mU3^2/MR^2]) + Log[mQ3^2/MR^2]*
     ((16*(7*m3^6 - 2*mQ3^4*mU3^2 - m3^4*(7*mQ3^2 + 6*mU3^2) + 
         m3^2*(3*mQ3^4 + 5*mQ3^2*mU3^2)))/((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)) - 
      (128*m3*Xt)/(mQ3^2 - mU3^2) + 
      (32*(7*m3^4 + 3*mQ3^2*mU3^2 - 2*m3^2*(2*mQ3^2 + 3*mU3^2))*Xt^2)/
       ((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)) + 
      (128*m3*(3*mQ3^2 + mU3^2)*Xt^3)/(mQ3^2 - mU3^2)^3 - 
      (16*(4*m3^8 - 3*mQ3^4*mU3^2*(5*mQ3^2 + mU3^2) + 
         m3^6*(7*mQ3^2 + 3*mU3^2) - m3^4*(29*mQ3^4 + 15*mQ3^2*mU3^2 + 
           6*mU3^4) + m3^2*(16*mQ3^6 + 31*mQ3^4*mU3^2 + 7*mQ3^2*mU3^4))*Xt^4)/
       ((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^3) - 
      (128*m3*mQ3^2*Xt^5)/((m3^2 - mQ3^2)*(mQ3^2 - mU3^2)^3) + 
      ((16*(-2*mQ3^4*mU3^4 + 2*m3^6*(mQ3^2 + mU3^2) + 4*m3^2*mQ3^2*mU3^2*
            (mQ3^2 + mU3^2) - m3^4*(mQ3^4 + 8*mQ3^2*mU3^2 + mU3^4)))/
         ((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) - (64*m3^3*Xt)/
         ((m3^2 - mQ3^2)*(m3^2 - mU3^2)) + 
        (32*(2*m3^8*(mQ3^2 + mU3^2) + 2*mQ3^4*mU3^4*(mQ3^2 + mU3^2) - 
           4*m3^2*mQ3^2*mU3^2*(mQ3^2 + mU3^2)^2 + 3*m3^4*(mQ3^2 + mU3^2)^3 - 
           2*m3^6*(3*mQ3^4 + 2*mQ3^2*mU3^2 + 3*mU3^4))*Xt^2)/
         ((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)^2) + 
        (128*(-2*m3*mQ3^2*mU3^2 + m3^3*(mQ3^2 + mU3^2))*Xt^3)/
         ((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^2) - 
        (16*(mQ3^2 + mU3^2)*(4*m3^8*(mQ3^2 + mU3^2) + 4*mQ3^4*mU3^4*
            (mQ3^2 + mU3^2) - 8*m3^2*mQ3^2*mU3^2*(mQ3^2 + mU3^2)^2 - 
           2*m3^6*(5*mQ3^4 + 6*mQ3^2*mU3^2 + 5*mU3^4) + 
           m3^4*(5*mQ3^6 + 19*mQ3^4*mU3^2 + 19*mQ3^2*mU3^4 + 5*mU3^6))*Xt^4)/
         ((m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)^4) - 
        (64*m3*(mQ3^2 + mU3^2)*(-2*mQ3^2*mU3^2 + m3^2*(mQ3^2 + mU3^2))*Xt^5)/
         ((m3^2 - mQ3^2)*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)^4))*
       Log[mU3^2/MR^2]) + ((-64*m3^4)/(-m3^2 + mQ3^2)^2 + 
      (256*m3^3*Xt)/((m3^2 - mQ3^2)*(mQ3^2 - mU3^2)) + 
      (128*m3*(2*m3^2 - mQ3^2 - mU3^2)*Xt^3)/(mQ3^2 - mU3^2)^3 - 
      (32*(-2*m3^2 + mQ3^2 + mU3^2)*Xt^4)/(mQ3^2 - mU3^2)^3)*
     PolyLog[2, 1 - m3^2/mQ3^2] + ((-64*m3^4)/(m3^2 - mU3^2)^2 + 
      (256*m3^3*Xt)/((m3^2 - mU3^2)*(-mQ3^2 + mU3^2)) + 
      (128*m3*(-2*m3^2 + mQ3^2 + mU3^2)*Xt^3)/(mQ3^2 - mU3^2)^3 + 
      (32*(-2*m3^2 + mQ3^2 + mU3^2)*Xt^4)/(mQ3^2 - mU3^2)^3)*
     PolyLog[2, 1 - m3^2/mU3^2]) + 
  ((12*Xt^4)/(mQ3^2 - mU3^2)^2 + (6 + (12*Xt^2)/(mQ3^2 - mU3^2) + 
      ((-6*mQ3^2)/(mQ3^2 - mU3^2)^3 - (6*mU3^2)/(mQ3^2 - mU3^2)^3)*Xt^4)*
     Log[mQ3^2/MR^2] + (6 - (12*Xt^2)/(mQ3^2 - mU3^2) + 
      ((6*mQ3^2)/(mQ3^2 - mU3^2)^3 + (6*mU3^2)/(mQ3^2 - mU3^2)^3)*Xt^4)*
     Log[mU3^2/MR^2])*(3*((2*m3^4)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) - 
       (2*mQ3^2*mU3^2)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) + 
       ((-2*(2*m3^8 - 2*m3^6*mQ3^2 + m3^4*mQ3^4 - 2*m3^6*mU3^2 + m3^4*mU3^4))/
          (3*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) + (8*m3^3*Xt)/
          (3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)))*Log[m3^2/MR^2] + 
       ((4*m3^2*mQ3^2)/(3*(m3^2 - mQ3^2)^2) - (2*mQ3^4)/
          (3*(m3^2 - mQ3^2)^2) - (8*m3*mQ3^2*Xt)/(3*(m3^2 - mQ3^2)*
           (mQ3^2 - mU3^2)))*Log[mQ3^2/MR^2] + 
       ((4*m3^2*mU3^2)/(3*(m3^2 - mU3^2)^2) - (2*mU3^4)/
          (3*(m3^2 - mU3^2)^2) - (8*m3*mU3^2*Xt)/(3*(m3^2 - mU3^2)*
           (-mQ3^2 + mU3^2)))*Log[mU3^2/MR^2])^2 + 
    2*((291*m3^8 - 282*m3^6*mQ3^2 + 93*m3^4*mQ3^4 - 18*m3^2*mQ3^6 - 
        360*m3^6*msq^2 + 300*m3^4*mQ3^2*msq^2 - 180*m3^2*mQ3^4*msq^2 - 
        282*m3^6*mU3^2 - 112*m3^4*mQ3^2*mU3^2 + 208*m3^2*mQ3^4*mU3^2 - 
        6*mQ3^6*mU3^2 + 300*m3^4*msq^2*mU3^2 + 240*m3^2*mQ3^2*msq^2*mU3^2 - 
        60*mQ3^4*msq^2*mU3^2 + 93*m3^4*mU3^4 + 208*m3^2*mQ3^2*mU3^4 - 
        169*mQ3^4*mU3^4 - 180*m3^2*msq^2*mU3^4 - 60*mQ3^2*msq^2*mU3^4 - 
        18*m3^2*mU3^6 - 6*mQ3^2*mU3^6)/(18*(m3^2 - mQ3^2)^2*
        (m3^2 - mU3^2)^2) - (8*(14*m3^3 - 3*m3*mQ3^2 - 30*m3*msq^2 - 
         3*m3*mU3^2)*Xt)/(9*(m3^2 - mQ3^2)*(m3^2 - mU3^2)) + 
      ((262*m3^16 - 878*m3^14*mQ3^2 + 885*m3^12*mQ3^4 - 320*m3^10*mQ3^6 + 
          49*m3^8*mQ3^8 - 878*m3^14*mU3^2 + 2972*m3^12*mQ3^2*mU3^2 - 
          2932*m3^10*mQ3^4*mU3^2 + 1004*m3^8*mQ3^6*mU3^2 - 
          158*m3^6*mQ3^8*mU3^2 + 885*m3^12*mU3^4 - 2932*m3^10*mQ3^2*mU3^4 + 
          2404*m3^8*mQ3^4*mU3^4 - 368*m3^6*mQ3^6*mU3^4 - m3^4*mQ3^8*mU3^4 - 
          320*m3^10*mU3^6 + 1004*m3^8*mQ3^2*mU3^6 - 368*m3^6*mQ3^4*mU3^6 - 
          452*m3^4*mQ3^6*mU3^6 + 144*m3^2*mQ3^8*mU3^6 + 49*m3^8*mU3^8 - 
          158*m3^6*mQ3^2*mU3^8 - m3^4*mQ3^4*mU3^8 + 144*m3^2*mQ3^6*mU3^8 - 
          36*mQ3^8*mU3^8)/(9*(m3^2 - mQ3^2)^4*(m3^2 - mU3^2)^4) - 
        (8*(27*m3^11 - 82*m3^9*mQ3^2 + 57*m3^7*mQ3^4 - 82*m3^9*mU3^2 + 
           220*m3^7*mQ3^2*mU3^2 - 142*m3^5*mQ3^4*mU3^2 + 57*m3^7*mU3^4 - 
           142*m3^5*mQ3^2*mU3^4 + 87*m3^3*mQ3^4*mU3^4)*Xt)/
         (9*(m3^2 - mQ3^2)^3*(m3^2 - mU3^2)^3) + (64*m3^6*Xt^2)/
         (9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2))*Log[m3^2/MR^2]^2 + 
      ((128*m3^16 - 556*m3^14*mQ3^2 + 714*m3^12*mQ3^4 - 373*m3^10*mQ3^6 + 
          96*m3^8*mQ3^8 - 42*m3^6*mQ3^10 + 38*m3^4*mQ3^12 - 9*m3^2*mQ3^14 + 
          60*m3^12*mQ3^2*msq^2 + 120*m3^10*mQ3^4*msq^2 - 
          180*m3^8*mQ3^6*msq^2 - 90*m3^10*mQ3^2*msq^4 + 60*m3^8*mQ3^4*msq^4 + 
          30*m3^6*mQ3^6*msq^4 - 340*m3^14*mU3^2 + 1712*m3^12*mQ3^2*mU3^2 - 
          2209*m3^10*mQ3^4*mU3^2 + 933*m3^8*mQ3^6*mU3^2 + 
          106*m3^6*mQ3^8*mU3^2 - 230*m3^4*mQ3^10*mU3^2 + 
          47*m3^2*mQ3^12*mU3^2 - 3*mQ3^14*mU3^2 - 60*m3^12*msq^2*mU3^2 - 
          300*m3^10*mQ3^2*msq^2*mU3^2 - 180*m3^8*mQ3^4*msq^2*mU3^2 + 
          540*m3^6*mQ3^6*msq^2*mU3^2 + 90*m3^10*msq^4*mU3^2 + 
          210*m3^8*mQ3^2*msq^4*mU3^2 - 210*m3^6*mQ3^4*msq^4*mU3^2 - 
          90*m3^4*mQ3^6*msq^4*mU3^2 + 262*m3^12*mU3^4 - 1903*m3^10*mQ3^2*
           mU3^4 + 2771*m3^8*mQ3^4*mU3^4 - 1506*m3^6*mQ3^6*mU3^4 + 
          316*m3^4*mQ3^8*mU3^4 + 37*m3^2*mQ3^10*mU3^4 - mQ3^12*mU3^4 + 
          180*m3^10*msq^2*mU3^4 + 540*m3^8*mQ3^2*msq^2*mU3^4 - 
          180*m3^6*mQ3^4*msq^2*mU3^4 - 540*m3^4*mQ3^6*msq^2*mU3^4 - 
          270*m3^8*msq^4*mU3^4 - 90*m3^6*mQ3^2*msq^4*mU3^4 + 
          270*m3^4*mQ3^4*msq^4*mU3^4 + 90*m3^2*mQ3^6*msq^4*mU3^4 + 
          5*m3^10*mU3^6 + 723*m3^8*mQ3^2*mU3^6 - 1186*m3^6*mQ3^4*mU3^6 + 
          610*m3^4*mQ3^6*mU3^6 - 119*m3^2*mQ3^8*mU3^6 - 17*mQ3^10*mU3^6 - 
          180*m3^8*msq^2*mU3^6 - 420*m3^6*mQ3^2*msq^2*mU3^6 + 
          420*m3^4*mQ3^4*msq^2*mU3^6 + 180*m3^2*mQ3^6*msq^2*mU3^6 + 
          270*m3^6*msq^4*mU3^6 - 90*m3^4*mQ3^2*msq^4*mU3^6 - 
          150*m3^2*mQ3^4*msq^4*mU3^6 - 30*mQ3^6*msq^4*mU3^6 - 43*m3^8*mU3^8 - 
          60*m3^6*mQ3^2*mU3^8 + 162*m3^4*mQ3^4*mU3^8 - 84*m3^2*mQ3^6*mU3^8 + 
          21*mQ3^8*mU3^8 + 60*m3^6*msq^2*mU3^8 + 120*m3^4*mQ3^2*msq^2*mU3^8 - 
          180*m3^2*mQ3^4*msq^2*mU3^8 - 90*m3^4*msq^4*mU3^8 + 
          60*m3^2*mQ3^2*msq^4*mU3^8 + 30*mQ3^4*msq^4*mU3^8)/
         (18*(m3^2 - mQ3^2)^4*(m3^2 - mU3^2)^3*(mQ3^2 - mU3^2)) + 
        ((40*m3*(mQ3^2 - msq^2)^2)/(3*(m3^2 - mQ3^2)^2*(mQ3^2 - mU3^2)) + 
          (4*(18*m3^5 + m3^3*mQ3^2 + 3*m3*mQ3^4 - 37*m3^3*mU3^2 - 
             7*m3*mQ3^2*mU3^2 + 22*m3*mU3^4))/(9*(m3^2 - mU3^2)^2*
            (-mQ3^2 + mU3^2)) - (4*(5*m3^5*mQ3^2 - 5*m3^3*mQ3^4 + 
             3*m3*mQ3^6 - 5*m3^5*mU3^2 - 4*m3*mQ3^4*mU3^2 + 5*m3^3*mU3^4 + 
             4*m3*mQ3^2*mU3^4 - 3*m3*mU3^6))/(9*(m3^2 - mQ3^2)^2*
            (m3^2 - mU3^2)^2) + (4*(33*m3^9*mQ3^4 - 75*m3^7*mQ3^6 + 
             33*m3^5*mQ3^8 + 8*m3^3*mQ3^10 - 3*m3*mQ3^12 - 17*m3^9*mQ3^2*
              mU3^2 - 33*m3^7*mQ3^4*mU3^2 + 158*m3^5*mQ3^6*mU3^2 - 
             106*m3^3*mQ3^8*mU3^2 + 10*m3*mQ3^10*mU3^2 + 49*m3^7*mQ3^2*
              mU3^4 - 78*m3^5*mQ3^4*mU3^4 - 11*m3^3*mQ3^6*mU3^4 + 
             28*m3*mQ3^8*mU3^4 - 5*m3^7*mU3^6 - 22*m3^5*mQ3^2*mU3^6 + 
             43*m3^3*mQ3^4*mU3^6 - 12*m3*mQ3^6*mU3^6 + 5*m3^5*mU3^8 + 
             5*m3^3*mQ3^2*mU3^8 - 10*m3*mQ3^4*mU3^8 - 3*m3^3*mU3^10 + 
             3*m3*mQ3^2*mU3^10))/(9*(m3^2 - mQ3^2)^3*(m3^2 - mU3^2)^2*
            (mQ3^2 - mU3^2)^2))*Xt + (64*m3^2*mQ3^4*Xt^2)/
         (9*(m3^2 - mQ3^2)^2*(mQ3^2 - mU3^2)^2))*Log[mQ3^2/MR^2]^2 + 
      ((-10*(12*m3^8 - 23*m3^6*mQ3^2 + 11*m3^4*mQ3^4 + 18*m3^6*msq^2 - 
           15*m3^4*mQ3^2*msq^2 + 9*m3^2*mQ3^4*msq^2 - 23*m3^6*mU3^2 + 
           44*m3^4*mQ3^2*mU3^2 - 21*m3^2*mQ3^4*mU3^2 - 15*m3^4*msq^2*mU3^2 - 
           12*m3^2*mQ3^2*msq^2*mU3^2 + 3*mQ3^4*msq^2*mU3^2 + 11*m3^4*mU3^4 - 
           21*m3^2*mQ3^2*mU3^4 + 10*mQ3^4*mU3^4 + 9*m3^2*msq^2*mU3^4 + 
           3*mQ3^2*msq^2*mU3^4))/(9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2) + 
        (80*m3*msq^2*Xt)/(3*(m3^2 - mQ3^2)*(m3^2 - mU3^2)))*Log[msq^2/MR^2] + 
      ((-5*(2*mQ3^2 - 2*msq^2)*(-2*m3^4 - 3*m3^2*mQ3^2 + mQ3^4 + 
           3*m3^2*msq^2 + mQ3^2*msq^2))/(6*(-m3^2 + mQ3^2)^3) - 
        (5*(-2*msq^2 + 2*mU3^2)*(-2*m3^4 + 3*m3^2*msq^2 - 3*m3^2*mU3^2 + 
           msq^2*mU3^2 + mU3^4))/(6*(-m3^2 + mU3^2)^3) + 
        ((40*m3*(mQ3^2 - msq^2)^2)/(3*(m3^2 - mQ3^2)^2*(mQ3^2 - mU3^2)) + 
          (40*m3*(-msq^2 + mU3^2)^2)/(3*(m3^2 - mU3^2)^2*(-mQ3^2 + mU3^2)))*
         Xt)*Log[msq^2/MR^2]^2 + 
      ((52*m3^10*mQ3^2 - 105*m3^8*mQ3^4 + 53*m3^6*mQ3^6 - 308*m3^10*mU3^2 + 
          764*m3^8*mQ3^2*mU3^2 - 622*m3^6*mQ3^4*mU3^2 + 163*m3^4*mQ3^6*
           mU3^2 - 9*m3^2*mQ3^8*mU3^2 - 90*m3^6*mQ3^2*msq^2*mU3^2 + 
          180*m3^4*mQ3^4*msq^2*mU3^2 - 90*m3^2*mQ3^6*msq^2*mU3^2 + 
          365*m3^8*mU3^4 - 985*m3^6*mQ3^2*mU3^4 + 942*m3^4*mQ3^4*mU3^4 - 
          271*m3^2*mQ3^6*mU3^4 - 3*mQ3^8*mU3^4 + 90*m3^6*msq^2*mU3^4 - 
          210*m3^4*mQ3^2*msq^2*mU3^4 + 150*m3^2*mQ3^4*msq^2*mU3^4 - 
          30*mQ3^6*msq^2*mU3^4 + 18*m3^6*mU3^6 + 25*m3^4*mQ3^2*mU3^6 - 
          246*m3^2*mQ3^4*mU3^6 + 131*mQ3^6*mU3^6 + 30*m3^4*msq^2*mU3^6 - 
          60*m3^2*mQ3^2*msq^2*mU3^6 + 30*mQ3^4*msq^2*mU3^6 - 106*m3^4*mU3^8 + 
          279*m3^2*mQ3^2*mU3^8 - 125*mQ3^4*mU3^8 - 9*m3^2*mU3^10 - 
          3*mQ3^2*mU3^10)/(9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^3*
          (mQ3^2 - mU3^2)) + (8*(16*m3^7 - 16*m3^5*mQ3^2 + 55*m3^5*mU3^2 - 
           53*m3^3*mQ3^2*mU3^2 + 3*m3*mQ3^4*mU3^2 - 30*m3^3*msq^2*mU3^2 + 
           30*m3*mQ3^2*msq^2*mU3^2 - 61*m3^3*mU3^4 + 53*m3*mQ3^2*mU3^4 + 
           3*m3*mU3^6)*Xt)/(9*(m3^2 - mQ3^2)*(m3^2 - mU3^2)^2*
          (mQ3^2 - mU3^2)) + ((-10*(6*m3^4*msq^2 - 9*m3^2*msq^4 + 
             2*m3^4*mU3^2 + 18*m3^2*msq^2*mU3^2 - 3*msq^4*mU3^2 - 
             3*m3^2*mU3^4 + mU3^6))/(9*(m3^2 - mU3^2)^3) + 
          (40*(-6*m3*msq^4 + m3^3*mU3^2 + 12*m3*msq^2*mU3^2 - m3*mU3^4)*Xt)/
           (9*(m3^2 - mU3^2)^2*(-mQ3^2 + mU3^2)))*Log[msq^2/MR^2])*
       Log[mU3^2/MR^2] + ((128*m3^10 + 44*m3^8*mQ3^2 - 2*m3^6*mQ3^4 + 
          9*m3^4*mQ3^6 - 60*m3^6*mQ3^2*msq^2 + 90*m3^4*mQ3^2*msq^4 - 
          556*m3^8*mU3^2 + 68*m3^6*mQ3^2*mU3^2 - 33*m3^4*mQ3^4*mU3^2 - 
          6*m3^2*mQ3^6*mU3^2 + 60*m3^6*msq^2*mU3^2 - 120*m3^4*mQ3^2*msq^2*
           mU3^2 - 90*m3^4*msq^4*mU3^2 - 60*m3^2*mQ3^2*msq^4*mU3^2 + 
          702*m3^6*mU3^4 - 141*m3^4*mQ3^2*mU3^4 + 36*m3^2*mQ3^4*mU3^4 - 
          3*mQ3^6*mU3^4 + 120*m3^4*msq^2*mU3^4 + 180*m3^2*mQ3^2*msq^2*mU3^4 + 
          60*m3^2*msq^4*mU3^4 - 30*mQ3^2*msq^4*mU3^4 - 347*m3^4*mU3^6 + 
          50*m3^2*mQ3^2*mU3^6 - mQ3^4*mU3^6 - 180*m3^2*msq^2*mU3^6 + 
          30*msq^4*mU3^6 + 48*m3^2*mU3^8 - 17*mQ3^2*mU3^8 + 21*mU3^10)/
         (18*(m3^2 - mU3^2)^4*(-mQ3^2 + mU3^2)) - 
        (4*(-18*m3^7*mQ3^2 + m3^5*mQ3^4 + 3*m3^3*mQ3^6 + 30*m3^3*mQ3^2*
            msq^4 + 18*m3^7*mU3^2 + 69*m3^5*mQ3^2*mU3^2 - 
           11*m3^3*mQ3^4*mU3^2 - 3*m3*mQ3^6*mU3^2 - 60*m3^3*mQ3^2*msq^2*
            mU3^2 - 30*m3^3*msq^4*mU3^2 - 30*m3*mQ3^2*msq^4*mU3^2 - 
           86*m3^5*mU3^4 - 59*m3^3*mQ3^2*mU3^4 + 10*m3*mQ3^4*mU3^4 + 
           60*m3^3*msq^2*mU3^4 + 60*m3*mQ3^2*msq^2*mU3^4 + 
           30*m3*msq^4*mU3^4 + 99*m3^3*mU3^6 + 4*m3*mQ3^2*mU3^6 - 
           60*m3*msq^2*mU3^6 - 27*m3*mU3^8)*Xt)/(9*(m3^2 - mU3^2)^3*
          (-mQ3^2 + mU3^2)^2) + (64*m3^2*mU3^4*Xt^2)/(9*(m3^2 - mU3^2)^2*
          (-mQ3^2 + mU3^2)^2))*Log[mU3^2/MR^2]^2 + 
      Log[m3^2/MR^2]*((-864*m3^12 + 1942*m3^10*mQ3^2 - 1460*m3^8*mQ3^4 + 
          351*m3^6*mQ3^6 - 9*m3^4*mQ3^8 + 180*m3^10*msq^2 - 
          240*m3^8*mQ3^2*msq^2 + 270*m3^6*mQ3^4*msq^2 - 90*m3^4*mQ3^6*msq^2 + 
          1942*m3^10*mU3^2 - 4088*m3^8*mQ3^2*mU3^2 + 2853*m3^6*mQ3^4*mU3^2 - 
          572*m3^4*mQ3^6*mU3^2 - 3*m3^2*mQ3^8*mU3^2 - 240*m3^8*msq^2*mU3^2 - 
          180*m3^6*mQ3^2*msq^2*mU3^2 + 90*m3^4*mQ3^4*msq^2*mU3^2 - 
          30*m3^2*mQ3^6*msq^2*mU3^2 - 1460*m3^8*mU3^4 + 2853*m3^6*mQ3^2*
           mU3^4 - 1894*m3^4*mQ3^4*mU3^4 + 345*m3^2*mQ3^6*mU3^4 + 
          270*m3^6*msq^2*mU3^4 + 90*m3^4*mQ3^2*msq^2*mU3^4 + 351*m3^6*mU3^6 - 
          572*m3^4*mQ3^2*mU3^6 + 345*m3^2*mQ3^4*mU3^6 - 48*mQ3^6*mU3^6 - 
          90*m3^4*msq^2*mU3^6 - 30*m3^2*mQ3^2*msq^2*mU3^6 - 9*m3^4*mU3^8 - 
          3*m3^2*mQ3^2*mU3^8)/(9*(m3^2 - mQ3^2)^3*(m3^2 - mU3^2)^3) + 
        (8*(85*m3^7 - 75*m3^5*mQ3^2 + 3*m3^3*mQ3^4 - 60*m3^5*msq^2 + 
           30*m3^3*mQ3^2*msq^2 - 75*m3^5*mU3^2 + 59*m3^3*mQ3^2*mU3^2 + 
           30*m3^3*msq^2*mU3^2 + 3*m3^3*mU3^4)*Xt)/(9*(m3^2 - mQ3^2)^2*
          (m3^2 - mU3^2)^2) + 
        ((-2*(64*m3^16 - 180*m3^14*mQ3^2 + 122*m3^12*mQ3^4 - 31*m3^10*mQ3^6 + 
             30*m3^8*mQ3^8 - 7*m3^6*mQ3^10 - 204*m3^14*mU3^2 + 
             604*m3^12*mQ3^2*mU3^2 - 369*m3^10*mQ3^4*mU3^2 + 
             13*m3^8*mQ3^6*mU3^2 - 49*m3^6*mQ3^8*mU3^2 + 13*m3^4*mQ3^10*
              mU3^2 + 234*m3^12*mU3^4 - 779*m3^10*mQ3^2*mU3^4 + 
             503*m3^8*mQ3^4*mU3^4 - 13*m3^6*mQ3^6*mU3^4 + 61*m3^4*mQ3^8*
              mU3^4 - 18*m3^2*mQ3^10*mU3^4 - 101*m3^10*mU3^6 + 
             401*m3^8*mQ3^2*mU3^6 - 239*m3^6*mQ3^4*mU3^6 - 53*m3^4*mQ3^6*
              mU3^6 - 6*m3^2*mQ3^8*mU3^6 + 6*mQ3^10*mU3^6 + 13*m3^8*mU3^8 - 
             76*m3^6*mQ3^2*mU3^8 + 43*m3^4*mQ3^4*mU3^8 + 24*m3^2*mQ3^6*
              mU3^8 - 6*mQ3^8*mU3^8))/(9*(m3^2 - mQ3^2)^4*(m3^2 - mU3^2)^3*
            (mQ3^2 - mU3^2)) + (4*(68*m3^11 - 183*m3^9*mQ3^2 + 
             120*m3^7*mQ3^4 + 3*m3^5*mQ3^6 - 137*m3^9*mU3^2 + 
             353*m3^7*mQ3^2*mU3^2 - 215*m3^5*mQ3^4*mU3^2 - 17*m3^3*mQ3^6*
              mU3^2 + 75*m3^7*mU3^4 - 184*m3^5*mQ3^2*mU3^4 + 
             105*m3^3*mQ3^4*mU3^4 + 12*m3*mQ3^6*mU3^4)*Xt)/
           (9*(m3^2 - mQ3^2)^3*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)) - 
          (128*m3^4*mQ3^2*Xt^2)/(9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)*
            (mQ3^2 - mU3^2)))*Log[mQ3^2/MR^2] + 
        ((-40*m3^12 + 200*m3^10*mQ3^2 - 60*m3^8*mQ3^4 + 20*m3^6*mQ3^6 + 
            200*m3^10*mU3^2 - 840*m3^8*mQ3^2*mU3^2 + 420*m3^6*mQ3^4*mU3^2 - 
            140*m3^4*mQ3^6*mU3^2 - 60*m3^8*mU3^4 + 420*m3^6*mQ3^2*mU3^4 + 
            20*m3^6*mU3^6 - 140*m3^4*mQ3^2*mU3^6)/(9*(m3^2 - mQ3^2)^3*
            (m3^2 - mU3^2)^3) - (40*(m3^7 + 5*m3^5*mQ3^2 + 5*m3^5*mU3^2 - 
             11*m3^3*mQ3^2*mU3^2)*Xt)/(9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2))*
         Log[msq^2/MR^2] + ((2*(64*m3^16 - 204*m3^14*mQ3^2 + 
             234*m3^12*mQ3^4 - 101*m3^10*mQ3^6 + 13*m3^8*mQ3^8 - 
             180*m3^14*mU3^2 + 604*m3^12*mQ3^2*mU3^2 - 779*m3^10*mQ3^4*
              mU3^2 + 401*m3^8*mQ3^6*mU3^2 - 76*m3^6*mQ3^8*mU3^2 + 
             122*m3^12*mU3^4 - 369*m3^10*mQ3^2*mU3^4 + 503*m3^8*mQ3^4*mU3^4 - 
             239*m3^6*mQ3^6*mU3^4 + 43*m3^4*mQ3^8*mU3^4 - 31*m3^10*mU3^6 + 
             13*m3^8*mQ3^2*mU3^6 - 13*m3^6*mQ3^4*mU3^6 - 53*m3^4*mQ3^6*
              mU3^6 + 24*m3^2*mQ3^8*mU3^6 + 30*m3^8*mU3^8 - 
             49*m3^6*mQ3^2*mU3^8 + 61*m3^4*mQ3^4*mU3^8 - 6*m3^2*mQ3^6*mU3^8 - 
             6*mQ3^8*mU3^8 - 7*m3^6*mU3^10 + 13*m3^4*mQ3^2*mU3^10 - 
             18*m3^2*mQ3^4*mU3^10 + 6*mQ3^6*mU3^10))/(9*(m3^2 - mQ3^2)^3*
            (m3^2 - mU3^2)^4*(mQ3^2 - mU3^2)) - 
          (4*(68*m3^11 - 137*m3^9*mQ3^2 + 75*m3^7*mQ3^4 - 183*m3^9*mU3^2 + 
             353*m3^7*mQ3^2*mU3^2 - 184*m3^5*mQ3^4*mU3^2 + 120*m3^7*mU3^4 - 
             215*m3^5*mQ3^2*mU3^4 + 105*m3^3*mQ3^4*mU3^4 + 3*m3^5*mU3^6 - 
             17*m3^3*mQ3^2*mU3^6 + 12*m3*mQ3^4*mU3^6)*Xt)/(9*(m3^2 - mQ3^2)^2*
            (m3^2 - mU3^2)^3*(mQ3^2 - mU3^2)) + (128*m3^4*mU3^2*Xt^2)/
           (9*(m3^2 - mQ3^2)*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)))*
         Log[mU3^2/MR^2]) + Log[mQ3^2/MR^2]*
       ((308*m3^10*mQ3^2 - 365*m3^8*mQ3^4 - 18*m3^6*mQ3^6 + 106*m3^4*mQ3^8 + 
          9*m3^2*mQ3^10 - 90*m3^6*mQ3^4*msq^2 - 30*m3^4*mQ3^6*msq^2 - 
          52*m3^10*mU3^2 - 764*m3^8*mQ3^2*mU3^2 + 985*m3^6*mQ3^4*mU3^2 - 
          25*m3^4*mQ3^6*mU3^2 - 279*m3^2*mQ3^8*mU3^2 + 3*mQ3^10*mU3^2 + 
          90*m3^6*mQ3^2*msq^2*mU3^2 + 210*m3^4*mQ3^4*msq^2*mU3^2 + 
          60*m3^2*mQ3^6*msq^2*mU3^2 + 105*m3^8*mU3^4 + 622*m3^6*mQ3^2*mU3^4 - 
          942*m3^4*mQ3^4*mU3^4 + 246*m3^2*mQ3^6*mU3^4 + 125*mQ3^8*mU3^4 - 
          180*m3^4*mQ3^2*msq^2*mU3^4 - 150*m3^2*mQ3^4*msq^2*mU3^4 - 
          30*mQ3^6*msq^2*mU3^4 - 53*m3^6*mU3^6 - 163*m3^4*mQ3^2*mU3^6 + 
          271*m3^2*mQ3^4*mU3^6 - 131*mQ3^6*mU3^6 + 90*m3^2*mQ3^2*msq^2*
           mU3^6 + 30*mQ3^4*msq^2*mU3^6 + 9*m3^2*mQ3^2*mU3^8 + 3*mQ3^4*mU3^8)/
         (9*(m3^2 - mQ3^2)^3*(m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)) - 
        (8*(16*m3^7 + 55*m3^5*mQ3^2 - 61*m3^3*mQ3^4 + 3*m3*mQ3^6 - 
           30*m3^3*mQ3^2*msq^2 - 16*m3^5*mU3^2 - 53*m3^3*mQ3^2*mU3^2 + 
           53*m3*mQ3^4*mU3^2 + 30*m3*mQ3^2*msq^2*mU3^2 + 3*m3*mQ3^2*mU3^4)*
          Xt)/(9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)*(mQ3^2 - mU3^2)) + 
        ((-10*(2*m3^4*mQ3^2 - 3*m3^2*mQ3^4 + mQ3^6 + 6*m3^4*msq^2 + 
             18*m3^2*mQ3^2*msq^2 - 9*m3^2*msq^4 - 3*mQ3^2*msq^4))/
           (9*(m3^2 - mQ3^2)^3) + (40*(m3^3*mQ3^2 - m3*mQ3^4 + 
             12*m3*mQ3^2*msq^2 - 6*m3*msq^4)*Xt)/(9*(m3^2 - mQ3^2)^2*
            (mQ3^2 - mU3^2)))*Log[msq^2/MR^2] + 
        ((-14*m3^10*mQ3^2 + 17*m3^8*mQ3^4 - 35*m3^6*mQ3^6 + 29*m3^4*mQ3^8 - 
            9*m3^2*mQ3^10 - 2*m3^10*mU3^2 + 28*m3^8*mQ3^2*mU3^2 + 
            53*m3^6*mQ3^4*mU3^2 - 75*m3^4*mQ3^6*mU3^2 + 35*m3^2*mQ3^8*mU3^2 - 
            3*mQ3^10*mU3^2 + 3*m3^8*mU3^4 - 57*m3^6*mQ3^2*mU3^4 + 
            23*m3^4*mQ3^4*mU3^4 - m3^2*mQ3^6*mU3^4 - 4*mQ3^8*mU3^4 - 
            m3^6*mU3^6 + 19*m3^4*mQ3^2*mU3^6 - 9*m3^2*mQ3^4*mU3^6 + 
            3*mQ3^6*mU3^6)/(9*(m3^2 - mQ3^2)^3*(m3^2 - mU3^2)^3) + 
          (4*(m3^7*mQ3^4 + 9*m3^5*mQ3^6 - 10*m3^3*mQ3^8 + 6*m3*mQ3^10 - 
             34*m3^7*mQ3^2*mU3^2 + 29*m3^5*mQ3^4*mU3^2 + 11*m3^3*mQ3^6*
              mU3^2 - 20*m3*mQ3^8*mU3^2 + m3^7*mU3^4 + 59*m3^5*mQ3^2*mU3^4 - 
             88*m3^3*mQ3^4*mU3^4 + 38*m3*mQ3^6*mU3^4 - m3^5*mU3^6 - 
             9*m3^3*mQ3^2*mU3^6 + 8*m3*mQ3^4*mU3^6)*Xt)/(9*(m3^2 - mQ3^2)^2*
            (m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)^2) - (128*m3^2*mQ3^2*mU3^2*Xt^2)/
           (9*(m3^2 - mQ3^2)*(m3^2 - mU3^2)*(-mQ3^2 + mU3^2)^2))*
         Log[mU3^2/MR^2]) + 
      (-(-128*m3^12 + 294*m3^10*mQ3^2 - 98*m3^8*mQ3^4 - 9*m3^6*mQ3^6 - 
           20*m3^4*mQ3^8 + 9*m3^2*mQ3^10 + 346*m3^10*mU3^2 - 
           902*m3^8*mQ3^2*mU3^2 + 251*m3^6*mQ3^4*mU3^2 + 151*m3^4*mQ3^6*
            mU3^2 - 41*m3^2*mQ3^8*mU3^2 + 3*mQ3^10*mU3^2 - 280*m3^8*mU3^4 + 
           1025*m3^6*mQ3^2*mU3^4 - 417*m3^4*mQ3^4*mU3^4 - 
           41*m3^2*mQ3^6*mU3^4 + mQ3^8*mU3^4 + 13*m3^6*mU3^6 - 
           391*m3^4*mQ3^2*mU3^6 + 167*m3^2*mQ3^4*mU3^6 + 19*mQ3^6*mU3^6 + 
           37*m3^4*mU3^8 + 34*m3^2*mQ3^2*mU3^8 - 23*mQ3^4*mU3^8)/
         (9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^3*(mQ3^2 - mU3^2)) + 
        (8*(18*m3^5 + m3^3*mQ3^2 + 3*m3*mQ3^4 - 37*m3^3*mU3^2 - 
           7*m3*mQ3^2*mU3^2 + 22*m3*mU3^4)*Xt)/(9*(m3^2 - mU3^2)^2*
          (-mQ3^2 + mU3^2)))*PolyLog[2, 1 - m3^2/mQ3^2] + 
      (((2*m3^2 - 2*msq^2)*((40*mQ3^4)/(3*(-m3^2 + mQ3^2)^3) - 
           (10*mQ3^2)/(-m3^2 + mQ3^2)^2 - 10/(3*(-m3^2 + mQ3^2)) - 
           (40*mQ3^2*msq^2)/(3*(-m3^2 + mQ3^2)^3) + (10*msq^2)/
            (-m3^2 + mQ3^2)^2 - (40*msq^2*mU3^2)/(3*(-m3^2 + mU3^2)^3) + 
           (40*mU3^4)/(3*(-m3^2 + mU3^2)^3) + (10*msq^2)/(-m3^2 + mU3^2)^2 - 
           (10*mU3^2)/(-m3^2 + mU3^2)^2 - 10/(3*(-m3^2 + mU3^2))))/2 + 
        (80*(m3^2 - msq^2)*(m3^3*mQ3^2 - 2*m3^3*msq^2 + m3*mQ3^2*msq^2 + 
           m3^3*mU3^2 - 2*m3*mQ3^2*mU3^2 + m3*msq^2*mU3^2)*Xt)/
         (3*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2))*PolyLog[2, 1 - m3^2/msq^2] + 
      ((-5*(2*mQ3^2 - 2*msq^2)*(-2*m3^4 - 3*m3^2*mQ3^2 + mQ3^4 + 
           3*m3^2*msq^2 + mQ3^2*msq^2))/(3*(-m3^2 + mQ3^2)^3) + 
        (80*m3*(mQ3^2 - msq^2)^2*Xt)/(3*(m3^2 - mQ3^2)^2*(mQ3^2 - mU3^2)))*
       PolyLog[2, 1 - msq^2/mQ3^2] + 
      ((-128*m3^12 + 346*m3^10*mQ3^2 - 280*m3^8*mQ3^4 + 13*m3^6*mQ3^6 + 
          37*m3^4*mQ3^8 + 294*m3^10*mU3^2 - 902*m3^8*mQ3^2*mU3^2 + 
          1025*m3^6*mQ3^4*mU3^2 - 391*m3^4*mQ3^6*mU3^2 + 
          34*m3^2*mQ3^8*mU3^2 - 98*m3^8*mU3^4 + 251*m3^6*mQ3^2*mU3^4 - 
          417*m3^4*mQ3^4*mU3^4 + 167*m3^2*mQ3^6*mU3^4 - 23*mQ3^8*mU3^4 - 
          9*m3^6*mU3^6 + 151*m3^4*mQ3^2*mU3^6 - 41*m3^2*mQ3^4*mU3^6 + 
          19*mQ3^6*mU3^6 - 20*m3^4*mU3^8 - 41*m3^2*mQ3^2*mU3^8 + 
          mQ3^4*mU3^8 + 9*m3^2*mU3^10 + 3*mQ3^2*mU3^10)/(9*(m3^2 - mQ3^2)^3*
          (m3^2 - mU3^2)^2*(mQ3^2 - mU3^2)) + 
        (8*(18*m3^5 - 37*m3^3*mQ3^2 + 22*m3*mQ3^4 + m3^3*mU3^2 - 
           7*m3*mQ3^2*mU3^2 + 3*m3*mU3^4)*Xt)/(9*(m3^2 - mQ3^2)^2*
          (mQ3^2 - mU3^2)))*PolyLog[2, 1 - m3^2/mU3^2] + 
      ((-24*m3^10*mQ3^2 + 28*m3^8*mQ3^4 - 68*m3^6*mQ3^6 + 58*m3^4*mQ3^8 - 
          18*m3^2*mQ3^10 + 24*m3^10*mU3^2 + 220*m3^6*mQ3^4*mU3^2 - 
          188*m3^4*mQ3^6*mU3^2 + 70*m3^2*mQ3^8*mU3^2 - 6*mQ3^10*mU3^2 - 
          28*m3^8*mU3^4 - 220*m3^6*mQ3^2*mU3^4 + 16*m3^2*mQ3^6*mU3^4 - 
          8*mQ3^8*mU3^4 + 68*m3^6*mU3^6 + 188*m3^4*mQ3^2*mU3^6 - 
          16*m3^2*mQ3^4*mU3^6 - 58*m3^4*mU3^8 - 70*m3^2*mQ3^2*mU3^8 + 
          8*mQ3^4*mU3^8 + 18*m3^2*mU3^10 + 6*mQ3^2*mU3^10)/
         (18*(-m3^2 + mQ3^2)^3*(m3^2 - mU3^2)^3) - 
        (4*(10*m3^5*mQ3^2 - 10*m3^3*mQ3^4 + 6*m3*mQ3^6 - 10*m3^5*mU3^2 - 
           8*m3*mQ3^4*mU3^2 + 10*m3^3*mU3^4 + 8*m3*mQ3^2*mU3^4 - 6*m3*mU3^6)*
          Xt)/(9*(m3^2 - mQ3^2)^2*(m3^2 - mU3^2)^2))*
       PolyLog[2, 1 - mQ3^2/mU3^2] + 
      ((-5*(-2*msq^2 + 2*mU3^2)*(-2*m3^4 + 3*m3^2*msq^2 - 3*m3^2*mU3^2 + 
           msq^2*mU3^2 + mU3^4))/(3*(-m3^2 + mU3^2)^3) + 
        (80*m3*(-msq^2 + mU3^2)^2*Xt)/(3*(m3^2 - mU3^2)^2*(-mQ3^2 + mU3^2)))*
       PolyLog[2, 1 - msq^2/mU3^2])))
    )
    ];

lambda3LATASASDegenerate = With[{
    k = 1/(4*Pi)^2,
    gt = Yu[3,3],
    MR = SCALE,
    MS = Sqrt[Sqrt[Abs[msq2[3,3] msu2[3,3]]]],
    Xtt = xt/Sqrt[Sqrt[Abs[msq2[3,3] msu2[3,3]]]]
    },
    (1/2)*k^3*(gt^4 g3^4) (
    (* 16/27 Xtt^3 (3568 - 1664 Log[MS^2/MR^2] + 444 Log[MS^2/MR^2]^2 - 2259 Zeta[3]) +  *)
     8/9 Xtt^2 (149 - 448 Log[MS^2/MR^2] + 912 Log[MS^2/MR^2]^2 - 990 Zeta[3]) - 
     64/27 Xtt (823 - 200 Log[MS^2/MR^2] + 954 Log[MS^2/MR^2]^2 - 477 Zeta[3]) - 
     4/135 (176 \[Pi]^4 + 960 \[Pi]^2 Log[2]^2 + 9900 Log[MS^2/MR^2]^2 - 24840 Log[MS^2/MR^2]^3 - 
        5 (4577 + 192 Log[2]^4 + 4608 PolyLog[4, 1/2] - 9864 Zeta[3]) + 
        10 Log[MS^2/MR^2] (-877 + 432 Zeta[3]))
    ) /. { a:PolyLog[_,_] :> N[a], a:Zeta[__] :> N[a] }
];
