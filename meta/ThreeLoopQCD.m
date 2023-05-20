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

BeginPackage["ThreeLoopQCD`", {"SARAH`"}];
EndPackage[];

(* parameters *)
{Q};

BeginPackage["ThreeLoopQCD`", {"SARAH`"}];

GetMTopMSbarOverMTopPole::usage = "Returns the ratio of the MS-bar top
 mass over the top pole mass, Eq. (10) of arxiv:hep-ph/9912391 and at
 an arbitrary scale Q The logarithmic scale dependence is taken from
 arxiv:hep-ph/9911434.

Example call:

   mTopMSbar = MtPole GetMTopMSbarOverMTopPole[{1,1,1,1}]

The first function argument is a list of numbers, which are multiplied
by the tree-level, 1-loop, 2-loop and 3-loop contributions,
respectively.  I.e. by setting factors to 0, different loop orders can
be disabled.

The second function argument is the symbol for the quark under
consideration.  The quark symbol is needed to obtain the Casimir
operators which apper in the expressions of arxiv:hep-ph/9912391.  The
default is SARAH`TopQuark.

The third function argument is the symbol to be used for the
renormalization scale.  The default is Q.

The last two parameters are NH and NL, the number of heavy and light
quarks, respectively.  The default is NH = 1, NL = 5.

The full signature of function reads

   GetMTopMSbarOverMTopPole[loopOrder_List:{1,1,1,1},
                            quark_:SARAH`TopQuark, Q_:Q, NH_:1, NL_:5]
";

GetMTopPoleOverMTopMSbar::usage = "Returns the ratio of the pole top
 mass over the top MS-bar, being the inverse of
 GetMTopMSbarOverMTopPole[]."

Begin["`Private`"];

GetMTopMSbarOverMTopPole0L[] := 1;

GetMTopMSbarOverMTopPole1L[Q_, quark_] :=
    Module[{CF, colorPosition, as},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           as = SARAH`strongCoupling^2 / (4 Pi);
           CF as / Pi (-1 + 3/4 Log[FlexibleSUSY`Pole[FlexibleSUSY`M[quark]]^2/Q^2])
          ];

(* 2-loop coefficients *)
d12 = 7/128 - 3/4 Zeta[3] + 1/2 Pi^2 Log[2] - 5/16 Pi^2;
d22 = -1111/384 + 3/8 Zeta[3] - 1/4 Pi^2 Log[2] + Pi^2/12;
d32 = 71/96 + Pi^2/12;
d42 = 143/96 - 1/6 Pi^2;

(* 3-loop coefficients *)
a4 = PolyLog[4, 1/2];
d13 = (
    -2969/768 - 1/16 Pi^2 Zeta[3] - 81/16 Zeta[3] + 5/8 Zeta[5] + 29/4 Pi^2 Log[2]
    +1/2 Pi^2 Log[2]^2 - 613/192 Pi^2 - 1/48 Pi^4 - 1/2 Log[2]^4 - 12 a4
);
d23 = (
    13189/4608 - 19/16 Pi^2 Zeta[3] - 773/96 Zeta[3] + 45/16 Zeta[5] - 31/72 Pi^2 Log[2]
    -31/36 Pi^2 Log[2]^2 + 509/576 Pi^2 + 65/432 Pi^4 - 1/18 Log[2]^4 - 4/3 a4
);
d33 = (
    -1322545/124416 + 51/64 Pi^2 Zeta[3] + 1343/288 Zeta[3] - 65/32 Zeta[5]
    -115/72 Pi^2 Log[2] + 11/36 Pi^2 Log[2]^2 - 1955/3456 Pi^2 - 179/3456 Pi^4
    +11/72 Log[2]^4 + 11/3 a4
);
d43 = (
    1283/576 + 55/24 Zeta[3] - 11/9 Pi^2 Log[2] + 2/9 Pi^2 Log[2]^2 + 13/18 Pi^2
    -119/2160 Pi^4 + 1/9 Log[2]^4 + 8/3 a4
);
d53 = (
    1067/576 - 53/24 Zeta[3] + 8/9 Pi^2 Log[2] - 1/9 Pi^2 Log[2]^2 - 85/108 Pi^2
    +91/2160 Pi^4 + 1/9 Log[2]^4 + 8/3 a4
);
d63 = (
    70763/15552 + 89/144 Zeta[3] + 11/18 Pi^2 Log[2] - 1/9 Pi^2 Log[2]^2
    +175/432 Pi^2 + 19/2160 Pi^4 - 1/18 Log[2]^4 - 4/3 a4
);
d73 = (
    144959/15552 + 1/8 Pi^2 Zeta[3] - 109/144 Zeta[3] - 5/8 Zeta[5] + 32/9 Pi^2 Log[2]
    +1/18 Pi^2 Log[2]^2 - 449/144 Pi^2 - 43/1080 Pi^4 - 1/18 Log[2]^4 - 4/3 a4
);
d83 = -5917/3888 + 2/9 Zeta[3] + 13/108 Pi^2;
d93 = -9481/7776 + 11/18 Zeta[3] + 4/135 Pi^2;
d103 = -2353/7776 - 7/18 Zeta[3] - 13/108 Pi^2;

Get2LLogs[M_, Q_, CF_, CA_] := (
    (4*(-156 + 185*CA + 81*CF)*Log[M^2/Q^2] - 12*(-12 + 11*CA +
     9*CF)*Log[M^2/Q^2]^2 - 576*CF*Log[M^2/Q^2] +
     216*CF*Log[M^2/Q^2]^2)/384
);

GetMTopMSbarOverMTopPole2L[Q_, quark_, NH_, NL_] :=
    Module[{CF, CA, TR, colorPosition, as},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           TR = 1/2;
           as = SARAH`strongCoupling^2 / (4 Pi);
           CF (as/Pi)^2 (
               CF d12 + CA d22 + TR NL d32 + TR NH d42 +
               Get2LLogs[FlexibleSUSY`Pole[FlexibleSUSY`M[quark]], Q, CF, CA]
           )
          ];

(* Q-dependend terms arxiv:hep-ph/9911434 Eq. (29) *)
Get3LLogs[M_, Q_, CF_, CA_, TR_, NL_] := (
  CF^3 (
    - 9/128 Log[Q^2/M^2]^3
    - 27/128 Log[Q^2/M^2]^2
    + (
      -489/512
      + 45/32 Zeta[2]
      - 9/4 Zeta[2] Log[2]
      + 9/16 Zeta[3]
      ) Log[Q^2/M^2]
  )
  + CF^2 CA (
    33/128 Log[Q^2/M^2]^3
    +109/64 Log[Q^2/M^2]^2
    +(
      5813/1536
      - 61/16 Zeta[2]
      + 53/8 Zeta[2] Log[2]
      - 53/32 Zeta[3]
    ) Log[Q^2/M^2]
  )
  +CF^2 TR NL (
    - 3/32 Log[Q^2/M^2]^3
    - 13/32 Log[Q^2/M^2]^2
    + (
      65/384
      +7/8 Zeta[2]
      - 2 Zeta[2] Log[2]
      - 1/4 Zeta[3]
    ) Log[Q^2/M^2]
  )
  +CF^2 TR (
    - 3/32 Log[Q^2/M^2]^3
    - 13/32 Log[Q^2/M^2]^2
    + (
      -151/384
      + 2 Zeta[2]
      - 2 Zeta[2]  Log[2]
      - 1/4 Zeta[3]
    ) Log[Q^2/M^2]
  )
  +CF CA^2 (
    - 121/576 Log[Q^2/M^2]^3
    - 2341/1152 Log[Q^2/M^2]^2
    + (
      - 13243/1728
      + 11/12 Zeta[2]
      - 11/4 Zeta[2] Log[2]
      + 11/16 Zeta[3]
    ) Log[Q^2/M^2]
  )
  +CF CA TR NL (
    11/72 Log[Q^2/M^2]^3
    +373/288 Log[Q^2/M^2]^2
    + (
      869/216
      + 7/12 Zeta[2]
      + Zeta[2] Log[2]
      + 1/2 Zeta[3]
    ) Log[Q^2/M^2]
  )
  +CF CA TR (
    + 11/72 Log[Q^2/M^2]^3
    + 373/288 Log[Q^2/M^2]^2
    + (
      583/108
      - 13/6 Zeta[2]
      + Zeta[2] Log[2]
      + 1/2 Zeta[3]
    )  Log[Q^2/M^2]
  )
  +CF TR^2 NL^2 (
    - 1/36 Log[Q^2/M^2]^3
    - 13/72 Log[Q^2/M^2]^2
    + (
      -89/216
      -1/3 Zeta[2]
    ) Log[Q^2/M^2]
  )
  +CF TR^2 NL (
    - 1/18 Log[Q^2/M^2]^3
    - 13/36 Log[Q^2/M^2]^2
    + (
      -143/108
      + 1/3 Zeta[2]
    ) Log[Q^2/M^2]
  )
  +CF TR^2 (
    - 1/36 Log[Q^2/M^2]^3
    - 13/72 Log[Q^2/M^2]^2
    + (
      -197/216
      +2/3 Zeta[2]
    ) Log[Q^2/M^2]
  )
);

GetMTopMSbarOverMTopPole3L[Q_, quark_, NH_, NL_] :=
    Module[{CF, CA, TR, colorPosition, as},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           TR = 1/2;
           as = SARAH`strongCoupling^2 / (4 Pi);
           (as/Pi)^3 (
               CF (
                   CF^2 d13 + CF CA d23 + CA^2 d33 +
                   CF TR NL d43 + CF TR NH d53 + CA TR NL d63 +
                   CA TR NH d73 + TR^2 NL NH d83 + TR^2 NH^2 d93 +
                   TR^2 NL^2 d103
                  ) +
               Get3LLogs[FlexibleSUSY`M[quark], Q, CF, CA, TR, NL]
           )
          ];

GetMTopMSbarOverMTopPole[loopOrder_List:{1,1,1,1}, quark_:SARAH`TopQuark, Q_:Q, NH_:1, NL_:5] :=
    Module[{result, Mpole, h},
           Assert[Length[loopOrder] == 4];
           result = (
               h^0 GetMTopMSbarOverMTopPole0L[] +
               h^1 GetMTopMSbarOverMTopPole1L[Q, quark] +
               h^2 GetMTopMSbarOverMTopPole2L[Q, quark, NH, NL] +
               h^3 GetMTopMSbarOverMTopPole3L[Q, quark, NH, NL]
           );
           (* rewrite in terms of the running mass *)
           Mpole = FlexibleSUSY`M[quark] / result;
           result = Normal @ Series[
               result /. FlexibleSUSY`Pole[FlexibleSUSY`M[quark]] -> Mpole /.
                         FlexibleSUSY`Pole[FlexibleSUSY`M[quark]] -> Mpole /.
                         FlexibleSUSY`Pole[FlexibleSUSY`M[quark]] -> FlexibleSUSY`M[quark],
               {h,0,3}];
           (
               loopOrder[[1]] Coefficient[result, h, 0] +
               loopOrder[[2]] Coefficient[result, h, 1] +
               loopOrder[[3]] Coefficient[result, h, 2] +
               loopOrder[[4]] Coefficient[result, h, 3]
           )
          ];

GetMTopPoleOverMTopMSbar[loopOrder_List:{1,1,1,1}, quark_:SARAH`TopQuark, Q_:Q, NH_:1, NL_:5] :=
    Module[{h, inv},
           inv = Normal @ Series[1 / GetMTopMSbarOverMTopPole[{1, h, h^2, h^3}, quark, Q, NH, NL], {h,0,3}];
           (
               loopOrder[[1]] Coefficient[inv, h, 0] +
               loopOrder[[2]] Coefficient[inv, h, 1] +
               loopOrder[[3]] Coefficient[inv, h, 2] +
               loopOrder[[4]] Coefficient[inv, h, 3]
           )
          ];

End[];

EndPackage[];
