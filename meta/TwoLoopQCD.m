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

BeginPackage["TwoLoopQCD`", {"SARAH`"}];

GetDeltaMPoleOverMRunningQCDOneLoop::usage="Returns one-loop QCD
contributions to Delta M_f/m_f in the given renormalization scheme.
(M_f = pole mass, m_f = running mass)";

GetDeltaMPoleOverMRunningQCDOneLoopDRbar::usage="Returns one-loop QCD
contributions to Delta M_f/m_f in the DR-bar scheme.  Taken from
hep-ph/0210258 Eq. (58).  (M_f = pole mass, m_f = DR-bar mass)";

GetDeltaMPoleOverMRunningQCDOneLoopMSbar::usage="Returns one-loop QCD
contributions to Delta M_f/m_f in the MS-bar scheme.  Taken from
hep-ph/9803493 Eq. (16).  (M_f = pole mass, m_f = MS-bar mass)";

GetDeltaMPoleOverMRunningQCDTwoLoop::usage="Returns two-loop QCD
contributions to Delta M_f/m_f in the given renormalization scheme.
(M_f = pole mass, m_f = running mass)";

GetDeltaMPoleOverMRunningQCDTwoLoopDRbar::usage="Returns two-loop QCD
contributions to Delta M_f/m_f in the DR-bar scheme.  Taken from
hep-ph/0210258 Eq. (60)-(61).  (M_f = pole mass, m_f = DR-bar mass)";

GetDeltaMPoleOverMRunningQCDTwoLoopMSbar::usage="Returns two-loop QCD
contributions to Delta M_f/m_f in the MS-bar scheme.  Taken from
hep-ph/9803493 Eq. (17).  (M_f = pole mass, m_f = MS-bar mass)

Note: Eq. (17) of [hep-ph/9803493] is expressed in terms of
Log[Q^2/M_f^2].  This function, however, returns Delta M_f/m_f as a
function of Log[Q^2/m_f^2].  We have accounted for the difference.";

Begin["`Private`"];

SelectRenormalizationScheme::UnknownRenormalizationScheme = "Unknown\
 renormalization scheme `1`.";

(* ******* two-loop ******* *)

GetDeltaMPoleOverMRunningQCDTwoLoop[quark_, renScale_, renScheme_] :=
    Switch[renScheme,
           FlexibleSUSY`DRbar, GetDeltaMPoleOverMRunningQCDTwoLoopDRbar[quark,renScale],
           FlexibleSUSY`MSbar, GetDeltaMPoleOverMRunningQCDTwoLoopMSbar[quark,renScale],
           _, Message[SelectRenormalizationScheme::UnknownRenormalizationScheme, renScheme];
              Quit[1];
          ];

GetDeltaMPoleOverMRunningQCDTwoLoopDRbar[quark_ /; quark === TreeMasses`GetSMTopQuarkMultiplet[], renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, log, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           log = Log[(mf / renScale)^2];
           result = CF (alphaStrong / (4 Pi))^2 (
               -43 - 12 Zeta[2] + 26 log - 6 log^2
               + CF (-59/8 + 30 Zeta[2] - 48 Log[2] Zeta[2]
                     + 12 Zeta[3] + (3/2) log + (9/2) log^2)
               + CA (1093/24 - 8 Zeta[2] + 24 Log[2] Zeta[2]
                     - 6 Zeta[3] - (179/6) log + 11/2 log^2));
           result
          ];

GetDeltaMPoleOverMRunningQCDTwoLoopDRbar[quark_ /; quark === TreeMasses`GetSMBottomQuarkMultiplet[], renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, mt, log, logMt, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           mt = FlexibleSUSY`M[TreeMasses`GetSMTopQuarkMultiplet[]];
           log = Log[(mf / renScale)^2];
           logMt = Log[(mf / mt)^2];
           result = CF (alphaStrong / (4 Pi))^2 (
               -623/18 - 8 Zeta[2] + 26 log - 6 log^2
               + CF (-59/8 + 30 Zeta[2] - 48 Log[2] Zeta[2]
                     + 12 Zeta[3] + (3/2) log + (9/2) log^2)
               + CA (1093/24 - 8 Zeta[2] + 24 Log[2] Zeta[2]
                     - 6 Zeta[3] - (179/6) log + 11/2 log^2)
               - 13/3 logMt + logMt^2);
           result
          ];

GetDeltaMPoleOverMRunningQCDTwoLoopDRbar[_, _] := 0;

GetDeltaMPoleOverMRunningQCDTwoLoopMSbar[quark_, renScale_] :=
    Module[{CF, CA, colorPosition, alphaStrong, mf, log,
            t = 1/2, Nf = 6, I31, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           CA = SA`Casimir[SARAH`VectorG, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           log = Log[(renScale / mf)^2];
           I31 = 3/2 Zeta[3] - 6 Zeta[2] Log[2];
           result = (alphaStrong / (4 Pi))^2 (
               CF CA (1111/24 - 8 Zeta[2] - 4 I31 + 185/6 log + 11/2 log^2)
               - CF t Nf (71/6 + 8 Zeta[2] + 26/3 log + 2 log^2)
               + CF^2 (-71/8 + 30 Zeta[2] + 8 I31 - 9/2 log + 9/2 log^2)
               - 12 CF t (1 - 2 Zeta[2]));
           result
          ];

(* ******* one-loop ******* *)

GetDeltaMPoleOverMRunningQCDOneLoop[quark_, renScale_, renScheme_] :=
    Switch[renScheme,
           FlexibleSUSY`DRbar, GetDeltaMPoleOverMRunningQCDOneLoopDRbar[quark,renScale],
           FlexibleSUSY`MSbar, GetDeltaMPoleOverMRunningQCDOneLoopMSbar[quark,renScale],
           _, Message[SelectRenormalizationScheme::UnknownRenormalizationScheme, renScheme];
              Quit[1];
          ];

GetDeltaMPoleOverMRunningQCDOneLoopDRbar[quark_, renScale_] :=
    Module[{CF, colorPosition, alphaStrong, mf, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           result = CF alphaStrong / (4 Pi) (5 - 3 Log[(mf / renScale)^2]);
           result
          ];

GetDeltaMPoleOverMRunningQCDOneLoopMSbar[quark_, renScale_] :=
    Module[{CF, colorPosition, alphaStrong, mf, result},
           colorPosition = Position[SARAH`Gauge, SARAH`color][[1,1]];
           CF = SA`Casimir[quark, colorPosition];
           alphaStrong = SARAH`strongCoupling^2 / (4 Pi);
           mf = FlexibleSUSY`M[quark];
           result = CF alphaStrong / (4 Pi) (4 - 3 Log[(mf / renScale)^2]);
           result
          ];

End[];

EndPackage[];
