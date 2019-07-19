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

BeginPackage["AMuon`", {"SARAH`", "CXXDiagrams`", "TextFormatting`", "TreeMasses`", "LoopMasses`"}];

AMuonGetMuon::usage="";
AMuonGetMSUSY::usage="";

Begin["`Private`"];

AMuonGetMuon[] := If[TreeMasses`GetDimension[TreeMasses`GetSMMuonLeptonMultiplet[]] =!= 1,
                TreeMasses`GetSMMuonLeptonMultiplet[],
                Cases[SARAH`ParticleDefinitions[FlexibleSUSY`FSEigenstates],
                      {p_, {Description -> "Muon", ___}} -> p, 1][[1]]
               ]

GetCXXMuonIndex[] := If[TreeMasses`GetDimension[TreeMasses`GetSMMuonLeptonMultiplet[]] =!= 1,
                        1,
                        Null]

GetMinMass[particle_] :=
    Module[{dim = TreeMasses`GetDimension[particle],
            mStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]],
            tail},
           If[dim == 1,
              "model.get_" <> mStr <> "()",
              tail = ToString[GetDimension[particle] - GetDimensionStartSkippingGoldstones[particle] + 1];
              "model.get_" <> mStr <> "().tail<" <> tail <> ">().minCoeff()"
             ]
          ];

AMuonGetMSUSY[] :=
    Module[{susyParticles},
           susyParticles = Select[TreeMasses`GetSusyParticles[], IsElectricallyCharged];
           If[susyParticles === {},
              "return 0.;",
              "return Min(" <>
                 StringJoin[Riffle[GetMinMass /@ susyParticles, ", "]] <>
              ");"
             ]
          ];

End[];
EndPackage[];

