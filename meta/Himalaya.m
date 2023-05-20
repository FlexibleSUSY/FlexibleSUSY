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

BeginPackage["Himalaya`", {"SARAH`", "TextFormatting`", "CConversion`",
                           "Parameters`", "TreeMasses`", "Utils`"}];

FillHimalayaInput::usage = "Fills parameters into Himalaya input
 struct.";

{ RenormalizationScheme, Lambda3L, Lambda3LUncertainty,
  \[Mu], SARAH`g1, Susyno`LieGroups`g2, SARAH`g3, SARAH`vd, SARAH`vu,
  MSQ2, MSD2, MSU2, MSL2, MSE2,
  Au, Ad, Ae, Yu, Yd, Ye, M1, M2, M3, mA };

Begin["`Private`"];

AppendAtEnd[str_, char_] :=
    StringReplace[str, EndOfLine|EndOfString -> char];

HimalayaIsScalar[par_] :=
    MatchQ[Parameters`GetType[par], CConversion`ScalarType[_]];

FillHimalayaInput[inputPars_, struct_:"pars"] :=
    Module[{result = "",
            inPars = Parameters`DecreaseIndexLiterals[inputPars],
            pars = First /@ Reverse /@ inputPars,
            isScalarAu, isScalarAd, isScalarAe
           },
           isScalarAu = HimalayaIsScalar[Au /. inPars];
           isScalarAd = HimalayaIsScalar[Ad /. inPars];
           isScalarAe = HimalayaIsScalar[Ae /. inPars];
           result = Parameters`CreateLocalConstRefs[pars] <> "

pars.scale = MODELPARAMETER(scale);
pars.mu = Re(" <> CConversion`RValueToCFormString[\[Mu] /. inPars] <> ");
pars.g1 = Re(" <> CConversion`RValueToCFormString[Sqrt[5/3] SARAH`g1 FlexibleSUSY`GUTNormalization[SARAH`g1] /. inPars] <> ");
pars.g2 = Re(" <> CConversion`RValueToCFormString[Susyno`LieGroups`g2 /. inPars] <> ");
pars.g3 = Re(" <> CConversion`RValueToCFormString[g3 /. inPars] <> ");
pars.vd = Re(" <> CConversion`RValueToCFormString[vd /. inPars] <> ");
pars.vu = Re(" <> CConversion`RValueToCFormString[vu /. inPars] <> ");
pars.mq2 = Re(" <> CConversion`RValueToCFormString[MSQ2 /. inPars] <> ");
pars.md2 = Re(" <> CConversion`RValueToCFormString[MSD2 /. inPars] <> ");
pars.mu2 = Re(" <> CConversion`RValueToCFormString[MSU2 /. inPars] <> ");
pars.ml2 = Re(" <> CConversion`RValueToCFormString[MSL2 /. inPars] <> ");
pars.me2 = Re(" <> CConversion`RValueToCFormString[MSE2 /. inPars] <> ");
" <> If[isScalarAu,
        "pars.Au(2,2) = Re(" <> CConversion`RValueToCFormString[Au /. inPars] <> ");",
        "pars.Au = Re(" <> CConversion`RValueToCFormString[Au /. inPars] <> ");"
       ] <> "
" <> If[isScalarAd,
        "pars.Ad(2,2) = Re(" <> CConversion`RValueToCFormString[Ad /. inPars] <> ");",
        "pars.Ad = Re(" <> CConversion`RValueToCFormString[Ad /. inPars] <> ");"
       ] <> "
" <> If[isScalarAe,
        "pars.Ae(2,2) = Re(" <> CConversion`RValueToCFormString[Ae /. inPars] <> ");",
        "pars.Ae = Re(" <> CConversion`RValueToCFormString[Ae /. inPars] <> ");"
       ] <> "
pars.Yu = Re(" <> CConversion`RValueToCFormString[Yu /. inPars] <> ");
pars.Yd = Re(" <> CConversion`RValueToCFormString[Yd /. inPars] <> ");
pars.Ye = Re(" <> CConversion`RValueToCFormString[Ye /. inPars] <> ");
pars.M1 = " <> CConversion`RValueToCFormString[M1 /. inPars] <> ";
pars.M2 = " <> CConversion`RValueToCFormString[M2 /. inPars] <> ";
pars.MG = " <> CConversion`RValueToCFormString[M3 /. inPars] <> ";
pars.MA = " <> CConversion`RValueToCFormString[mA /. inPars] <> ";

const double msbar_scheme = " <> ToString[If[(RenormalizationScheme /. inPars) === DRbar, 0, 1]] <> ";
const double lambda_3L_eft = " <> CConversion`RValueToCFormString[Lambda3L /. inPars] <> ";
const double lambda_3L_uncertainty = " <> CConversion`RValueToCFormString[Lambda3LUncertainty /. inPars] <> ";
";
           AppendAtEnd[result, " \\"]
          ];

End[];

EndPackage[];
