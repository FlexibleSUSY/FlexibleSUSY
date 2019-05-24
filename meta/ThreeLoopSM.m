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

BeginPackage["ThreeLoopSM`", {"SARAH`", "Parameters`"}];

BetaSM::usage = "3-loop beta functions of the SM from SUSYHD v1.0.2
 and 4-loop beta functions from [arxiv:1508.00912]";

Begin["`Private`"];

subDir = FileNameJoin[{FlexibleSUSY`$flexiblesusyMetaDir, "SM"}];

IsDefinedAndEqual[descr_String, c_] :=
    Parameters`GetParameterFromDescription[descr] =!= Null &&
    Parameters`GetParameterFromDescription[descr] === c

G1GUTNormalization[] :=
    Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling] / Sqrt[3/5];

BetaSM[gc_] :=
    Switch[gc,
           SARAH`hyperchargeCoupling, Get[FileNameJoin[{subDir, "beta_g1.m"}]] / G1GUTNormalization[],
           SARAH`leftCoupling       , Get[FileNameJoin[{subDir, "beta_g2.m"}]],
           SARAH`strongCoupling     , Get[FileNameJoin[{subDir, "beta_g3.m"}]],
           SARAH`UpYukawa           , Get[FileNameJoin[{subDir, "beta_gt.m"  }]],
           SARAH`DownYukawa         , Get[FileNameJoin[{subDir, "beta_gb.m"  }]],
           SARAH`ElectronYukawa     , Get[FileNameJoin[{subDir, "beta_gtau.m"}]],
           \[Lambda]                , Get[FileNameJoin[{subDir, "beta_lambda.m"}]],
           m2                       , Get[FileNameJoin[{subDir, "beta_m2.m"}]],
           _, Which[IsDefinedAndEqual["SM Mu Parameter", gc],
                    Get[FileNameJoin[{subDir, "beta_m2.m"}]],
                    IsDefinedAndEqual["SM Higgs Selfcouplings", gc],
                    Get[FileNameJoin[{subDir, "beta_lambda.m"}]],
                    True, Print["Error: unknown coupling: ", gc]; {0,0,0}
                    ]
          ] /. ThreeLoopSM`Private`ToSARAHNamingConvention[] /. Zeta[n_] :> N[Zeta[n]];

(* Note:
   g1, g2, g3, gb are global variables in SARAH.
   Since version 4.14.0 the symbol m2 is defined in Susyno`LieGroups`.
 *)
ToSARAHNamingConvention[] := {
    g1 -> SARAH`hyperchargeCoupling G1GUTNormalization[],
    g2 -> SARAH`leftCoupling,
    g3 -> SARAH`strongCoupling,
    Global`gt -> SARAH`UpYukawa[3,3],
    gb -> SARAH`DownYukawa[3,3],
    Global`g\[Tau] -> SARAH`ElectronYukawa[3,3],
    \[Lambda] -> Parameters`GetParameterFromDescription["SM Higgs Selfcouplings"],
    m2        -> Parameters`GetParameterFromDescription["SM Mu Parameter"],
    Susyno`LieGroups`m2 -> Parameters`GetParameterFromDescription["SM Mu Parameter"]
};

End[];

EndPackage[];
