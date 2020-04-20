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

Off[FileByteCount::fdnfnd];

onSuppressedMessages::usage =
"Turns on corresponding warning messages";
onFileByteCount[] := On[FileByteCount::fdnfnd];
onFileByteCount // Utils`MakeUnknownInputDefinition;
onFileByteCount ~ SetAttributes ~ {Protected, Locked};

Begin["`internal`"];

showSingle[] := "done";
showSingle // Utils`MakeUnknownInputDefinition;
showSingle ~ SetAttributes ~ {Protected, Locked};

showTerms[info:_Integer, add:_Integer:3, pre:_String:""] :=
   "\033["<>ToString[Length@IntegerDigits@#+add]<>"D"<>
   pre<>"("<>ToString@#<>" terms"&[info];
showTerms // Utils`MakeUnknownInputDefinition;
showTerms ~ SetAttributes ~ {Protected, Locked};

showProgress[ind:_String:"      "] := "\r\033[K"<>ind<>"(in progress";
showProgress // Utils`MakeUnknownInputDefinition;
showProgress ~ SetAttributes ~ {Protected, Locked};

If[!$Notebooks,

DynamicCheckModelFiles /: Dynamic[DynamicCheckModelFiles] := showSingle[];
DynamicInitGaugeG /: Dynamic[DynamicInitGaugeG] := showSingle[];
DynamicInitFields /: Dynamic[DynamicInitFields] := showSingle[];
DynamicInitMisc /: Dynamic[DynamicInitMisc] := showSingle[];
DynamicCheckAnomalies /: Dynamic[DynamicCheckAnomalies] := showSingle[];

DynamicTermSuperpotentialNr /: Dynamic[DynamicTermSuperpotentialNr] := "";
DynamicTermSuperpotential /: Dynamic[DynamicTermSuperpotential] :=
   showTerms[Length@SuperPotential, 2];

DynamicCheckingCCSup /: Dynamic[DynamicCheckingCCSup] := showSingle[];

DynamicFTermNr /: Dynamic[DynamicFTermNr] := "";
DynamicFTermName /: Dynamic[DynamicFTermName] := showTerms@Length@SFieldList;

DynamicMatterNr /: Dynamic[DynamicMatterNr] := "";
DynamicMatterName /: Dynamic[DynamicMatterName] := showTerms[Length[SFieldList]^2];

DynamicSoftTermsCurrent /: Dynamic[DynamicSoftTermsCurrent] := showSingle[];

DynamicDGnr /: Dynamic[DynamicDGnr] := "";
DynamicDGname /: Dynamic[DynamicDGname] := showTerms@Length@Gauge;

DynamicKineticScalarNr /: Dynamic[DynamicKineticScalarNr] := "";
DynamicKineticScalarName /: Dynamic[DynamicKineticScalarName] := showTerms@AnzahlChiral;

DynamicKineticFermionNr /: Dynamic[DynamicKineticFermionNr] := "";
DynamicKineticFermionName /: Dynamic[DynamicKineticFermionName] := showTerms@AnzahlChiral;

DynamicDTermsNr /: Dynamic[DynamicDTermsNr] := "";
DynamicDTermsName /: Dynamic[DynamicDTermsName] := showTerms[AnzahlGauge*AnzahlChiral];

DynamicGauginoMatter /: Dynamic[DynamicGauginoMatter] := "";
DynamicGauginoMatterName /: Dynamic[DynamicGauginoMatterName] := showTerms[AnzahlGauge*AnzahlChiral];

DynamicGauginoVector /: Dynamic[DynamicGauginoVector] := "";
DynamicGauginoVectorName /: Dynamic[DynamicGauginoVectorName] := showTerms@AnzahlGauge;

DynamicVectorNr /: Dynamic[DynamicVectorNr] := "";
DynamicVectorName /: Dynamic[DynamicVectorName] := showTerms@AnzahlGauge;

DynamicGaugeTNr /: Dynamic[DynamicGaugeTNr] := "";
DynamicGaugeTName /: Dynamic[DynamicGaugeTName] := showTerms[AnzahlChiral+Length@Gauge];

DynamicStatusAddTerms /: Dynamic[DynamicStatusAddTerms[__]] := "";

DynamicRotateLag /: Dynamic[DynamicRotateLag[_]] := "14";

DynamicGFnr /: Dynamic[DynamicGFnr[_]] := "";
DynamicGFname /: Dynamic[DynamicGFname[_]] := showTerms@Length@gb;

DynamicSaveInfo /: Dynamic[DynamicSaveInfo[_]] := showSingle[];

DynamicUGT /: Dynamic[DynamicUGT[_]] := "";
DynamicUGTname /: Dynamic[DynamicUGTname[_]] := showTerms[Length@Particles@Current];

DynamicMMgaugeNr /: Dynamic[DynamicMMgaugeNr[_]] := "";
DynamicMMgaugeName /: Dynamic[DynamicMMgaugeName[_]] :=
   showTerms[ Length[DEFINITION[NameOfStates[[rotNr]]][GaugeSector]], 2];

nameMassQ = True;
DynamicNrMass /: Dynamic[DynamicNrMass[_]] := "";
DynamicNameMass /: Dynamic[DynamicNameMass[_]] :=
   showTerms[Length@If[nameMassQ,nameMassQ=False;mixBasis,nameMassQ=True;mixBasisNoFV], 5, ": "];

DynamicCalcTreeMasses /: Dynamic[DynamicCalcTreeMasses] := "for all eigenstates";
DynamicSpectrumFileInput /: Dynamic[DynamicSpectrumFileInput] := showSingle[];

DynamicProgressRGE /: Dynamic[DynamicProgressRGE[_]] = "";
DynamicCoupProgess /: Dynamic[DynamicCoupProgess[trace]] = showProgress["   "]<>")";
DynamicCoupProgess /: Dynamic[DynamicCoupProgess[_]] = showProgress["   "];

progressNrGV /: Dynamic[progressNrGV[_]] := "";
progressCurrentGV /: Dynamic[progressCurrentGV[_]] := showProgress[];

DynamicOneLoopNameMM /: Dynamic[DynamicOneLoopNameMM] := "";
DynamicOneLoopNrMM /: Dynamic[DynamicOneLoopNrMM] := Length@basis;

DynamicOneLoopTadName /: Dynamic[DynamicOneLoopTadName] := "";
DynamicOneLoopTadNr /: Dynamic[DynamicOneLoopTadNr] := "";
DynamicOneLoopTadNrAll /: Dynamic[DynamicOneLoopTadNrAll] := showProgress[];

DynamicOneLoopNameNM /: Dynamic[DynamicOneLoopNameNM] := "";
DynamicOneLoopNrNM /: Dynamic[DynamicOneLoopNrNM] := Length@listNotMixedMasses;
]; (* !$Notebooks *)

End[];
EndPackage[];
