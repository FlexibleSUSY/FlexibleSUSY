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

BeginPackage["ThreeLoopMSSM`", {"SARAH`", "Parameters`"}];

BetaMSSM::usage = "beta functions of the MSSM from
 http://www.liv.ac.uk/~dij/betas/allgennb.log";

Begin["`Private`"];

subDir = FileNameJoin[{FlexibleSUSY`$flexiblesusyMetaDir, "MSSM"}];

IsDefinedAndEqual[descr_String, c_] :=
    Parameters`GetParameterFromDescription[descr] =!= Null &&
    Parameters`GetParameterFromDescription[descr] === c

G1GUTNormalization[] :=
    Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling] / Sqrt[3/5];

BetaMSSM[gc_] :=
    Switch[gc,
           SARAH`hyperchargeCoupling, Get[FileNameJoin[{subDir, "beta_g1.m"}]] / G1GUTNormalization[],
           SARAH`leftCoupling       , Get[FileNameJoin[{subDir, "beta_g2.m"}]],
           SARAH`strongCoupling     , Get[FileNameJoin[{subDir, "beta_g3.m"}]],
           SARAH`UpYukawa           , Get[FileNameJoin[{subDir, "beta_Yu.m"}]],
           SARAH`DownYukawa         , Get[FileNameJoin[{subDir, "beta_Yd.m"}]],
           SARAH`ElectronYukawa     , Get[FileNameJoin[{subDir, "beta_Ye.m"}]],
           SARAH`TrilinearUp        , Get[FileNameJoin[{subDir, "beta_TYu.m"}]],
           SARAH`TrilinearDown      , Get[FileNameJoin[{subDir, "beta_TYd.m"}]],
           SARAH`TrilinearLepton    , Get[FileNameJoin[{subDir, "beta_TYe.m"}]],
           SARAH`SoftSquark         , Get[FileNameJoin[{subDir, "beta_mq2.m"}]],
           SARAH`SoftUp             , Get[FileNameJoin[{subDir, "beta_mu2.m"}]],
           SARAH`SoftDown           , Get[FileNameJoin[{subDir, "beta_md2.m"}]],
           SARAH`SoftLeftLepton     , Get[FileNameJoin[{subDir, "beta_ml2.m"}]],
           SARAH`SoftRightLepton    , Get[FileNameJoin[{subDir, "beta_me2.m"}]],
           _, Which[IsDefinedAndEqual["Bino Mass parameter", gc]  , Get[FileNameJoin[{subDir, "beta_M1.m"}]],
                    IsDefinedAndEqual["Wino Mass parameter", gc]  , Get[FileNameJoin[{subDir, "beta_M2.m"}]],
                    IsDefinedAndEqual["Gluino Mass parameter", gc], Get[FileNameJoin[{subDir, "beta_M3.m"}]],
                    IsDefinedAndEqual["Mu-parameter"         , gc], Get[FileNameJoin[{subDir, "beta_Mu.m"}]],
                    IsDefinedAndEqual["Bmu-parameter"        , gc], Get[FileNameJoin[{subDir, "beta_BMu.m"}]],
                    IsDefinedAndEqual["Softbreaking Down-Higgs Mass", gc], Get[FileNameJoin[{subDir, "beta_mHd2.m"}]],
                    IsDefinedAndEqual["Softbreaking Up-Higgs Mass"  , gc], Get[FileNameJoin[{subDir, "beta_mHu2.m"}]],
                    True, Print["Error: unknown coupling: ", gc]; {0,0,0}
                    ]
          ] /. ThreeLoopMSSM`Private`ToSARAHNamingConvention[] /. Zeta[s_] :> N[Zeta[s]];

(* Note:
   g1, g2, g3, Ye, M1, M2, M3, mb are in SARAH` context
 *)
ToSARAHNamingConvention[] := {
    g1 -> SARAH`hyperchargeCoupling G1GUTNormalization[],
    g2 -> SARAH`leftCoupling,
    g3 -> SARAH`strongCoupling,
    Global`Yt -> SARAH`UpYukawa,
    Global`Yb -> SARAH`DownYukawa,
    Ye -> SARAH`ElectronYukawa,
    M1 -> Parameters`GetParameterFromDescription["Bino Mass parameter"],
    M2 -> Parameters`GetParameterFromDescription["Wino Mass parameter"],
    M3 -> Parameters`GetParameterFromDescription["Gluino Mass parameter"],
    Global`Mu -> Parameters`GetParameterFromDescription["Mu-parameter"],
    Global`BMu-> Parameters`GetParameterFromDescription["Bmu-parameter"],
    Global`ht -> SARAH`TrilinearUp,
    Global`hb -> SARAH`TrilinearDown,
    Global`he -> SARAH`TrilinearLepton,
    Global`mq -> SARAH`SoftSquark,
    Global`mt -> SARAH`SoftUp,
    mb -> SARAH`SoftDown,
    Global`ml -> SARAH`SoftLeftLepton,
    Global`me -> SARAH`SoftRightLepton,
    Global`mh1 -> Parameters`GetParameterFromDescription["Softbreaking Down-Higgs Mass"],
    Global`mh2 -> Parameters`GetParameterFromDescription["Softbreaking Up-Higgs Mass"]
};

End[];

EndPackage[];
