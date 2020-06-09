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

BeginPackage["TerminalFormatting`", {"SARAH`", "Utils`"}];
Begin["`internal`"];
If[!$Notebooks,

showTerms[info:_Integer, add:_Integer:3, pre:_String:""] :=
   "\033["<>ToString[Length@IntegerDigits@#+add]<>"D"<>
   pre<>"("<>ToString@#<>" terms"&[info];
showTerms // Utils`MakeUnknownInputDefinition;
showTerms ~ SetAttributes ~ {Protected, Locked};

showProgress[ind:_String:"      "] := "\r\033[K"<>ind<>"(in progress";
showProgress // Utils`MakeUnknownInputDefinition;
showProgress ~ SetAttributes ~ {Protected, Locked};

define[s:_Symbol] := TagSetDelayed[s, Dynamic@s, ""];
define[(s:_Symbol)[p:_]] := TagSetDelayed[s, Dynamic@s@p, ""];
define[RuleDelayed[s:_Symbol, what:_ ]] := TagSetDelayed[s, Dynamic@s, what];
define[RuleDelayed[(s:_Symbol)[p:_], what:_ ]] := TagSetDelayed[s, Dynamic@s@p, what];
define // Utils`MakeUnknownInputDefinition;
define ~ SetAttributes ~ {Protected, Locked};

nameMassQ = True;

define /@ {
   DynamicTermSuperpotentialNr, DynamicFTermNr, DynamicDTermsNr, DynamicOneLoopNameMM,
   DynamicMatterNr, DynamicDGnr, DynamicKineticScalarNr, DynamicKineticFermionNr,
   DynamicGauginoMatter, DynamicGauginoVector, DynamicVectorNr, DynamicGaugeTNr,
   DynamicOneLoopNrMM, DynamicOneLoopTadName, DynamicOneLoopTadNr, DynamicOneLoopNrNM,
   DynamicOneLoopNameNM,
   DynamicCheckModelFiles :> "done", DynamicInitGaugeG :> "done", DynamicInitFields :> "done",
   DynamicInitMisc :> "done", DynamicCheckAnomalies :> "done", DynamicCheckingCCSup :> "done",
   DynamicSoftTermsCurrent :> "done", DynamicSpectrumFileInput :> "done",
   DynamicCalcTreeMasses :> "for all eigenstates",
   DynamicStatusAddTerms[_], DynamicGFnr[_], DynamicUGT[_], DynamicStatusAddTerms[_],
   DynamicMMgaugeNr[_], DynamicNrMass[_], DynamicProgressRGE[_], progressNrGV[_],
   DynamicSaveInfo[_] :> "done", DynamicRotateLag[_] :> "14",
   DynamicTermSuperpotential :> showTerms[Length@SuperPotential, 2],
   DynamicFTermName :> showTerms@Length@SFieldList,
   DynamicMatterName :> showTerms[Length[SFieldList]^2],
   DynamicDGname :> showTerms@Length@Gauge,
   DynamicKineticScalarName :> showTerms@AnzahlChiral,
   DynamicKineticFermionName :> showTerms@AnzahlChiral,
   DynamicDTermsName :> showTerms[AnzahlGauge*AnzahlChiral],
   DynamicGauginoMatterName :> showTerms[AnzahlGauge*AnzahlChiral],
   DynamicGauginoVectorName :> showTerms@AnzahlGauge,
   DynamicVectorName :> showTerms@AnzahlGauge,
   DynamicGaugeTName :> showTerms[AnzahlChiral+Length@Gauge],
   DynamicOneLoopNrMM :> Length@basis,
   DynamicOneLoopTadNrAll :> showProgress[],
   DynamicOneLoopNrNM :> Length@listNotMixedMasses,
   DynamicGFname[_] :> showTerms@Length@gb,
   DynamicUGTname[_] :> showTerms@Length@Particles@Current,
   DynamicCoupProgess@trace :> showProgress["   "]<>")",
   DynamicCoupProgess[_] :> showProgress["   "],
   progressCurrentGV[_] :> showProgress[],
   DynamicMMgaugeName[_] :> showTerms[ Length[DEFINITION[NameOfStates[[rotNr]]][GaugeSector]], 2],
   DynamicNameMass[_] :> showTerms[Length@If[nameMassQ,nameMassQ=False;mixBasis,nameMassQ=True;mixBasisNoFV], 5, ": "]
}

]; (* !$Notebooks *)
End[];
EndPackage[];
