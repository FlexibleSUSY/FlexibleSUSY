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

(*
   Converts SPheno.m to FlexibleSUSY.m.in file.
   Author: Alexander Voigt

   Note: The conversion from a SPheno.m to a FlexibleSUSY.m.in file is
         not unique, because in FlexibleSUSY additional settings
         can/should be made.

   Usage:

      Read SPheno.m in the current directory:

         $ math -noprompt -run "<< convert_SPheno_to_FlexibleSUSY.m; Quit[]"

      Read a specific SPheno.m file:

         $ math -noprompt -run "sphenoFile = \"/path/to/SPheno.m\"; << convert_SPheno_to_FlexibleSUSY.m; Quit[]"

 *)

DIAGONAL = UNITMATRIX[unknownMatrixSize];

GetSPhenoFile[] :=
    Module[{},
           If[!ValueQ[sphenoFile],
              sphenoFile = "SPheno.m";
             ];
           If[FileExistsQ[sphenoFile],
              Get[sphenoFile];,
              Print["Error: cannot find SPheno file ", sphenoFile];
              Quit[1];
             ];
          ];

ConvertToString[expr_] := ToString[expr, InputForm];

SetAttributes[WriteVariable, HoldFirst];

WriteVariable[var_] :=
    WriteVariable[Unevaluated[var], ToString[Unevaluated[var]]];

WriteVariable[var_, name_String] :=
    WriteVariable[Unevaluated[var], name, Null];

WriteVariable[var_, name_String, defaultValue_, F_:Identity] :=
    name <> " = " <>
    If[ValueQ[Unevaluated[var]],
       ConvertToString[F[var]],
       ConvertToString[defaultValue]
      ] <> ";\n\n";

ListDepth[{}]     := 1;
ListDepth[l_List] := 1 + Max[ListDepth /@ l];
ListDepth[_]      := 0;

GetSPhenoFile[];

output = "\
FSModelName = \"@CLASSNAME@\";
FSDefaultSARAHModel = UnknownSARAHModel;
FSRGELoopOrder = 2;

";

(*** OnlyLowEnergySPheno ***)

If[OnlyLowEnergySPheno === True,
   output = output <> "OnlyLowEnergyFlexibleSUSY = True;\n\n";
  ];
If[OnlyLowEnergySPheno === False,
   output = output <> "OnlyLowEnergyFlexibleSUSY = False;\n\n";
  ];

(*** MINPAR ***)

(* if multiple input parameter sets are given, take 1st one *)
If[ListDepth[MINPAR] === 3,
   MINPAR = MINPAR[[1]];
  ];
output = output <> WriteVariable[MINPAR, "MINPAR", {}];

(*** EXTPAR ***)

output = output <> WriteVariable[EXTPAR, "EXTPAR", {}];

output = output <> "\
FSAuxiliaryParameterInfo = {};

";

(*** ParametersToSolveTadpoles ***)

output = output <> WriteVariable[ParametersToSolveTadpoles, "EWSBOutputParameters", {}];

(*** BoundaryHighScale ***)

If[OnlyLowEnergySPheno =!= True,
   output = output <> WriteVariable[ConditionGUTscale, "HighScale"];
   output = output <> "HighScaleFirstGuess =. ; (* please set to a reasonable value! *) \n\n";
   (* if multiple BCs are given, take 1st one *)
   If[ListDepth[BoundaryHighScale] === 3,
      BoundaryHighScale = BoundaryHighScale[[1]];
     ];
   output = output <> WriteVariable[BoundaryHighScale, "HighScaleInput"];
  ];

(*** BoundarySUSYScale ***)

output = output <> "\
(* SUSYScale is the EWSB scale by default *)
" <>
WriteVariable[RenormalizationScale, "SUSYScale", LowScale, Sqrt];
output = output <> WriteVariable[RenormalizationScaleFirstGuess, "SUSYScaleFirstGuess", LowScaleFirstGuess, Sqrt];
output = output <> WriteVariable[BoundarySUSYScale, "SUSYScaleInput", {}];

(*** BoundaryLowScaleInput ***)

output = output <> "\
LowScale = LowEnergyConstant[MZ]; (* or LowEnergyConstant[MT] *)

LowScaleFirstGuess = LowScale;

(* Note: In FlexibleSUSY the strong, left and hypercharge gauge
   couplings are initialized automatically at the `LowScale' and don't
   need to be set in `LowScaleInput'.  The Standard Model Yukawa
   couplings can be set to `Automatic'.  The VEV(s) can be set from
   the running Z and/or W masses `MZDRbar' and/or `MWDRbar',
   respectively.  Example (in the MSSM):

   LowScaleInput = {
      {vd, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Cos[ArcTan[TanBeta]]},
      {vu, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Sin[ArcTan[TanBeta]]},
      {SARAH`UpYukawa, Automatic},
      {SARAH`DownYukawa, Automatic},
      {SARAH`ElectronYukawa, Automatic}
   };
*)
";
output = output <> WriteVariable[BoundaryLowScaleInput, "LowScaleInput", {}];

(*** initial guesses ***)

If[OnlyLowEnergySPheno =!= True,
   output = output <> WriteVariable[InitializationValues, "InitialGuessAtHighScale", {}];
  ];

output = output <> "\
InitialGuessAtLowScale = {
   (* Important: need to initialize the VEVs
      for the fermions to non-zero value! *)
   (*
   {SARAH`VEVSM1, LowEnergyConstant[vev] Cos[ArcTan[TanBeta]]},
   {SARAH`VEVSM2, LowEnergyConstant[vev] Sin[ArcTan[TanBeta]]},
   {SARAH`VEVSM , LowEnergyConstant[vev]},
   *)
   {SARAH`UpYukawa, Automatic},
   {SARAH`DownYukawa, Automatic},
   {SARAH`ElectronYukawa, Automatic}
};
"

(*** further settings ***)

output = output <> "
(* set to True to enable Pietro's 2L Higgs mass corrections
   in MSSM-like models (2 CP-even Higges, 1 CP-odd Higgs) *)
UseHiggs2LoopMSSM = False;
EffectiveMu =. ;

(* set to True to enable Pietro's 2L Higgs mass corrections
   in NMSSM-like models (3 CP-even Higges, 1 or 2 CP-odd Higgess) *)
UseHiggs2LoopNMSSM = False;
EffectiveMu =. ;
EffectiveMASqr =. ;

PotentialLSPParticles = {};

ExtraSLHAOutputBlocks = {
   {FlexibleSUSYOutput,
           {{1, Hold[SUSYScale]},
            {2, Hold[LowScale]} } }
};
";

WriteString["stdout", output];
