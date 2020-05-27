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

BeginPackage["TreeMasses`", {"SARAH`", "TextFormatting`", "CConversion`", "Parameters`", "Utils`"}];

FSMassMatrix::usage="Head of a mass matrix";

ConvertSarahMassMatrices::usage="creates list of mass matrices using
SARAH's MassMatrix[] function";

GetUnmixedParticleMasses::usage="returns list of masses of unmixed
 particles";

CreateMassGetter::usage="creates a C function for
the mass getter";

CreateSLHAPoleMassGetter::usage="creates a C function for the pole mass
getter";

CreateParticleLaTeXNames::usage="creates a list of the particle's
LaTeX names";

CreateParticleNames::usage="creates a list of the particle's names";

CreateParticleMixingNames::usage="creates a list of the mixing matrix element names";

CreateParticleEnum::usage="creates an enum of the particles";

CreateParticleMassEnum::usage="creates an enum of the particle masses";

CreateParticleMixingEnum::usage="creates enum with mixing matrices";

CreateParticleMultiplicity::usage="creates array of the particle
multiplicities";

prefixNamespace::usage="";
CreateFieldClassName::usage="creates the names of the class associated with
the given field";

FillSpectrumVector::usage="";

CreateMixingMatrixGetter::usage="creates a getter for the mixing
matrix";

CreateSLHAPoleMixingMatrixGetter::usage="creates a getter for the pole
mixing matrix";

CreateMassCalculationPrototype::usage="creates a C function prototype
from a mass matrix";

CreateMassCalculationFunction::usage="creates a C function that
calculates the mass eigenstates by diagonalizing the mass matrix";

CallMassCalculationFunctions::usage="creates C function calls of all
mass matrix calcualtion functions";

CreatePhysicalMassDefinition::usage="creates definition of physical
mass.";

CreateMixingMatrixDefinition::usage="creates definition of mixing
matrix";

ClearOutputParameters::usage="clears masses and mixing matrices";

CopyDRBarMassesToPoleMasses::usage="copies DRbar mass to pole mass";
CopyRunningMassesFromTo::usage="copies masses between two objects";

CreateMassArrayGetter::usage="";
CreateMassArraySetter::usage="";
CreateMixingArrayGetter::usage="";
CreateMixingArraySetter::usage="";

GetSMVEVExpr::usage = "Returns expression for SM-like VEV";

GetParticles::usage="returns list of particles";

GetSusyParticles::usage="returns list of susy particles";

GetSMParticles::usage="returns list of Standard Model particles";

GetGoldstoneBosons::usage="returns list of all goldstone bosons";

GetSMGoldstoneBosons::usage="returns list of all Standard Model
goldstone bosons";

GetVectorBosons::usage="returns list of all vector bosons";

GetDimension::usage="returns the size of the particle multiplet";

GetDimensionWithoutGoldstones::usage = "returns the size of the
particle multiplet, ignoring Goldstone bosons";

GetDimensionStartSkippingGoldstones::usage="return first index,
skipping goldstone bosons";

GetDimensionStartSkippingSMGoldstones::usage="return first index,
 skipping Standard Model goldstone bosons";

GetParticleIndices::usage = "returns list of particle indices with
 names";

FindMixingMatrixSymbolFor::usage="returns the mixing matrix symbol for
a given field";

SetUnrotatedParticles::usage="set list of unrotated particles in SARAH
format (result of ListUnmixed[EWSB] or content of
UnrotatedParticles.m)";

GetMassEigenstate::usage="get mass eigenstates symbol from mass
matrix";

GetMassMatrix::usage="get mass matrix from FSMassMatrix object";

GetMixingMatrixSymbol::usage="get mixing matrix symbol from mass
matrix";

MakeESSymbol::usage="Combines a list of particles to a symbol";

GetMassOfUnmixedParticle::usage="returns mass of unmixed particle";

GetMassType::usage="returns mass array type of particle";
GetMassMatrixType::usage="returns mass matrix type of particle";
GetMixingMatrixType::usage="returns mixing matrix type of particle";

ReplaceDependencies::usage="returs expression with dependencies
(ThetaW etc.) replaced by the user-defined expressions (";

ReplaceDependenciesReverse::usage="returs expression with dependencies
(ThetaW etc.) replaced back";

FindColorGaugeGroup::usage="returns triplet of color gauge coupling,
group and SARAH name";

FindColorGaugeCoupling::usage="returns symbol of color gauge coupling";
FindLeftGaugeCoupling::usage="returns symbol of weak (left) gauge coupling";
FindHyperchargeGaugeCoupling::usage="returns symbol of hypercharge gauge coupling";

CreateDependencePrototypes::usage="";
CreateDependenceFunctions::usage="";

ColorChargedQ::usage="";

FieldInfo::usage="";
includeLorentzIndices::usage="";
includeColourIndices::usage="";

IsParticle::usage = "returns True if argument is a particle or anti-particle, False otherwise"
IsScalar::usage="";
IsFermion::usage="";
IsVector::usage="";
IsGhost::usage="";
IsGoldstone::usage="";
IsSMGoldstone::usage="";
IsSMHiggs::usage="";
IsAuxiliary::usage="";
IsMajoranaFermion::usage="";
IsDiracFermion::usage="";
IsComplexScalar::usage="";
IsRealScalar::usage="";
IsComplexVector::usage="";
IsRealVector::usage="";
IsMassless::usage="";
IsUnmixed::usage="";
IsQuark::usage="";
IsLepton::usage="";
IsSMChargedLepton::usage="";
IsSMNeutralLepton::usage="";
IsSMLepton::usage="";
IsSMUpQuark::usage="";
IsSMDownQuark::usage="";
IsSMQuark::usage="";
IsSMParticle::usage="";
IsSMParticleElementwise::usage=
"For a given multiplet this function returns a list indicating whether
the element is a SM-like field or not.  The function assumes that BSM
fields are always heavier than the SM fields.";

IsElectricallyCharged::usage="";
ContainsGoldstone::usage="";

FSAntiField::usage = "Returns the anti-field of a given field";

GetSMChargedLeptons::usage="";
GetSMNeutralLeptons::usage="";
GetSMLeptons::usage="";
GetSMUpQuarks::usage="";
GetSMDownQuarks::usage="";
GetSMQuarks::usage="";
GetColoredParticles::usage="";
GetColorRepresentation::usage="";

GetUpQuark::usage="";
GetDownQuark::usage="";
GetUpLepton::usage="";
GetDownLepton::usage="";

GetSMUpQuark::usage        = "returns SM up quark, Fu[1] or Fu";
GetSMCharmQuark::usage     = "returns SM charm quark, Fu[3] or Ft";
GetSMTopQuark::usage       = "returns SM top quark, Fu[3] or Ft";
GetSMDownQuark::usage      = "returns SM down quark, Fd[1] or Fd";
GetSMStrangeQuark::usage   = "returns SM strange quark, Fd[2] or Fs";
GetSMBottomQuark::usage    = "returns SM bottom quark, Fd[3] or Fb";
GetSMElectronLepton::usage = "returns SM electon, Fe[1] or Fe";
GetSMMuonLepton::usage     = "returns SM muon, Fe[2] or Fm";
GetSMTauLepton::usage      = "returns SM tau, Fe[3] or Ftau";
GetSMNeutrino1::usage      = "returns SM neutrino 1, Fv[1] or FveL";
GetSMNeutrino2::usage      = "returns SM neutrino 2, Fv[2] or FvmL";
GetSMNeutrino3::usage      = "returns SM neutrino 3, Fv[3] or FvmL";
GetPhoton::usage           = "returns the photon";
GetGluon::usage            = "returns the gluon";
GetZBoson::usage           = "returns the Z boson";
GetWBoson::usage           = "returns the W boson";
GetHiggsBoson::usage       = "return Higgs boson(s)";
GetPseudoscalarHiggsBoson::usage = "";
GetChargedHiggsBoson::usage = "";

GetSMTopQuarkMultiplet::usage    = "Returns multiplet containing the top quark, Fu or Ft";
GetSMStrangeQuarkMultiplet::usage = "Returns multiplet containing the bottom quark, Fd or Fs";
GetSMBottomQuarkMultiplet::usage = "Returns multiplet containing the bottom quark, Fd or Fb";
GetSMElectronLeptonMultiplet::usage = "Returns multiplet containing the electron, Fe";
GetSMMuonLeptonMultiplet::usage     = "Returns multiplet containing the muon, Fe or Fm";
GetSMTauLeptonMultiplet::usage      = "Returns multiplet containing the tau lepton, Fe or Ftau";

GetMass::usage="wraps M[] head around particle";

GetElectricCharge::usage="Returns electric charge defined in the SARAH particle.m file";

StripGenerators::usage="removes all generators Lam, Sig, fSU2, fSU3
and removes Delta with the given indices";

CreateGenerationHelpers::usage="";
CallGenerationHelperFunctionName::usage="";

GetThirdGenerationMass::usage;
GetLightestMass::usage;

ReorderGoldstoneBosons::usage="";
CheckPoleMassesForTachyons::usage="";
CreateHiggsMassGetters::usage="";
CallPseudoscalarHiggsMassGetterFunction::usage="";

GetCorrespondingVectorBosons::usage="returns list of vector bosons
corresponding to a given goldstone boson";

(* exported for the use in LoopMasses.m *)
CallSVDFunction::usage="";
CallDiagonalizeSymmetricFunction::usage="";
CallDiagonalizeHermitianFunction::usage="";

FlagPoleTachyon::usage = "";
FlagRunningTachyon::usage = "";

GetStrongCoupling::usage = "Returns the name of the QCD coupling, e.g. g3. If not present returns Null."

Begin["`Private`"];

unrotatedParticles = {};

SetUnrotatedParticles[list_List] :=
    unrotatedParticles = ({#[[1]], #[[4]]})& /@ list;

GetVectorBosons[states_:FlexibleSUSY`FSEigenstates] :=
    #[[1]]& /@ Cases[SARAH`Particles[states], {_,_,_,V,__}] /.
       SARAH`diracSubBack1[SARAH`ALL] /.
       SARAH`diracSubBack2[SARAH`ALL];

(* Create list of mass eigenstate particles *)
GetParticles[states_:FlexibleSUSY`FSEigenstates] :=
    Module[{particles = {}},
           particles = Cases[SARAH`Masses[states],
                             HoldPattern[SARAH`Mass[p_] -> _] | HoldPattern[_ -> SARAH`MassGiven[p_]]
                             :> CConversion`GetHead[p]];
           particles = particles /.
                       SARAH`diracSubBack1[SARAH`ALL] /.
                       SARAH`diracSubBack2[SARAH`ALL];
           particles = Join[particles, GetVectorBosons[states]];
           particles = Select[particles, (!IsOfType[#,SARAH`NoField])&];
           Return[DeleteDuplicates[particles]];
          ];

GetSusyParticles[states_:FlexibleSUSY`FSEigenstates] :=
    Select[GetParticles[states], (!IsSMParticle[#] && !IsGhost[#])&];

GetSMParticles[states_:FlexibleSUSY`FSEigenstates] :=
    Select[GetParticles[states], IsSMParticle];

IsParticle[p_, states_:FlexibleSUSY`FSEigenstates] :=
    MemberQ[GetParticles[states], p] || MemberQ[GetParticles[states], FSAntiField[p]];

FieldInfo[field_, OptionsPattern[{includeLorentzIndices -> False,
	includeColourIndices -> False}]] := 
	Module[{fieldInfo = Cases[SARAH`Particles[FlexibleSUSY`FSEigenstates],
		{SARAH`getParticleName @ field, ___}][[1]]},
		fieldInfo = DeleteCases[fieldInfo, {SARAH`generation, 1}, {2}];
		
		fieldInfo = If[!OptionValue[includeLorentzIndices],
			DeleteCases[fieldInfo, {SARAH`lorentz, _}, {2}],
			fieldInfo];
		If[!OptionValue[includeColourIndices],
			DeleteCases[fieldInfo, {SARAH`color, _}, {2}],
			fieldInfo]
	]

IsOfType[sym_Symbol, type_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    SARAH`getType[sym, False, states] === type;

IsOfType[sym_[__], type_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    IsOfType[sym, type, states];

IsSMParticle[sym_List] := And @@ (IsSMParticle /@ sym);
IsSMParticle[Susyno`LieGroups`conj[sym_]] := IsSMParticle[sym];
IsSMParticle[SARAH`bar[sym_]] := IsSMParticle[sym];
IsSMParticle[sym_[__]] := IsSMParticle[sym];
IsSMParticle[sym_] := SARAH`SMQ[sym, Higgs -> True];

MakeTrueFalse[n_, t_] :=
    Join[Array[True&, n]] /; t >= n;

MakeTrueFalse[n_, t_] :=
    Join[Array[True&, t], Array[False&, n - t]];

Options[IsSMParticleElementwise] :=
    {
        IncludeHiggs -> False
    };

IsSMParticleElementwise[sym_, OptionsPattern[]] :=
    Which[
        IsSMLepton[sym], MakeTrueFalse[GetDimension[sym], 3],
        IsSMQuark[sym], MakeTrueFalse[GetDimension[sym], 3],
        True, (IsSMParticle[#] || IsSMGoldstone[#] ||
               (OptionValue[IncludeHiggs] && IsSMHiggs[#]))& /@
              Table[sym[i], {i, GetDimension[sym]}]
    ];

IsScalar[Susyno`LieGroups`conj[sym_]] := IsScalar[sym];
IsScalar[SARAH`bar[sym_]] := IsScalar[sym];
IsScalar[sym_] := IsOfType[sym, S];
IsScalar[sym_List] := And @@ (IsScalar /@ sym);

IsFermion[Susyno`LieGroups`conj[sym_]] := IsFermion[sym];
IsFermion[SARAH`bar[sym_]] := IsFermion[sym];
IsFermion[sym_] := IsOfType[sym, F];
IsFermion[sym_List] := And @@ (IsFermion /@ sym);

IsVector[Susyno`LieGroups`conj[sym_]] := IsVector[sym];
IsVector[SARAH`bar[sym_]] := IsVector[sym];
IsVector[sym_] := IsOfType[sym, V];
IsVector[sym_List] := And @@ (IsVector /@ sym);

IsGhost[Susyno`LieGroups`conj[sym_]] := IsGhost[sym];
IsGhost[SARAH`bar[sym_]] := IsGhost[sym];
IsGhost[sym_] := IsOfType[sym, G];
IsGhost[sym_List] := And @@ (IsGhost /@ sym);

IsGoldstone[Susyno`LieGroups`conj[sym_]] := IsGoldstone[sym];
IsGoldstone[SARAH`bar[sym_]] := IsGoldstone[sym];
IsGoldstone[sym_] := MemberQ[
    Join[GetGoldstoneBosons[],
         GetGoldstoneBosons[] /. a_[{idx__}] :> a[idx]],
    sym
];
IsGoldstone[sym_List] := And @@ (IsGoldstone /@ sym);

GetSMGoldstones[] :=
    Cases[SARAH`GoldstoneGhost /. a_[{idx__}] :> a[idx], {v_?IsSMParticle, goldstone_} :> goldstone];

IsSMGoldstone[Susyno`LieGroups`conj[sym_]] := IsSMGoldstone[sym];
IsSMGoldstone[SARAH`bar[sym_]] := IsSMGoldstone[sym];
IsSMGoldstone[sym_] :=
    MemberQ[GetSMGoldstones[], sym];

IsSMHiggs[Susyno`LieGroups`conj[sym_]] := IsSMHiggs[sym];
IsSMHiggs[SARAH`bar[sym_]] := IsSMHiggs[sym];
IsSMHiggs[sym_] :=
    Module[{higgs = Parameters`GetParticleFromDescription["Higgs"]},
           If[GetDimension[sym] == 1,
              SameQ[sym, higgs],
              SameQ[sym, higgs[GetDimensionStartSkippingGoldstones[higgs]]]
           ]
    ];

IsChargino[Susyno`LieGroups`conj[p_]] := IsChargino[p];
IsChargino[SARAH`bar[p_]] := IsChargino[p];
IsChargino[p_] :=
    p === Parameters`GetParticleFromDescription["Charginos"];

IsElectricallyCharged[par_] := GetElectricCharge[par] != 0;

ContainsGoldstone[sym_] := MemberQ[GetGoldstoneBosons[] /. a_[{idx__}] :> a, sym];

ContainsGoldstone[sym_[__]] := MemberQ[GetGoldstoneBosons[] /. a_[{idx__}] :> a, sym];

IsAuxiliary[Susyno`LieGroups`conj[sym_]] := IsAuxiliary[sym];
IsAuxiliary[SARAH`bar[sym_]] := IsAuxiliary[sym];
IsAuxiliary[sym_Symbol] := IsOfType[sym, A];

IsMajoranaFermion[Susyno`LieGroups`conj[sym_]] := IsMajoranaFermion[sym];
IsMajoranaFermion[SARAH`bar[sym_]] := IsMajoranaFermion[sym];
IsMajoranaFermion[sym_Symbol] :=
    And[IsFermion[sym], MemberQ[SARAH`MajoranaPart, sym]];

IsMajoranaFermion[sym_List] :=
    And @@ (IsMajoranaFermion /@ sym);

IsDiracFermion[Susyno`LieGroups`conj[sym_]] := IsDiracFermion[sym];
IsDiracFermion[SARAH`bar[sym_]] := IsDiracFermion[sym];
IsDiracFermion[sym_Symbol] :=
    And[IsFermion[sym], !MemberQ[SARAH`MajoranaPart, sym]];

IsDiracFermion[sym_List] :=
    And @@ (IsDiracFermion /@ sym);

IsComplexScalar[Susyno`LieGroups`conj[sym_]] := IsComplexScalar[sym];
IsComplexScalar[SARAH`bar[sym_]] := IsComplexScalar[sym];
IsComplexScalar[sym_Symbol] :=
    And[IsScalar[sym], Parameters`IsComplexParameter[sym]];

IsComplexScalar[sym_List] :=
    And @@ (IsComplexScalar /@ sym);

IsRealScalar[Susyno`LieGroups`conj[sym_]] := IsRealScalar[sym];
IsRealScalar[SARAH`bar[sym_]] := IsRealScalar[sym];
IsRealScalar[sym_Symbol] :=
    And[IsScalar[sym], Parameters`IsRealParameter[sym]];

IsRealScalar[sym_List] :=
    And[IsScalar[sym], And @@ (Parameters`IsRealParameter /@ sym)];

IsRealVector[p_] := IsVector[p] && Parameters`IsRealParameter[p];

IsComplexVector[p_] := IsVector[p] && Parameters`IsComplexParameter[p];

IsMassless[Susyno`LieGroups`conj[sym_], states_:FlexibleSUSY`FSEigenstates] := IsMassless[sym, states];

IsMassless[SARAH`bar[sym_], states_:FlexibleSUSY`FSEigenstates] := IsMassless[sym, states];

(* Massless ghosts are not stored in SARAH`Massless[FSEigenstates],
   so use mass of the corresponding vector boson. *)
IsMassless[sym_?IsGhost] :=
    Module[{v = Symbol["V" <> StringDrop[ToString[sym],1]]},
           Switch[RXi[v],
                  0, True,
                  _, IsMassless[v]
                 ]
         ];

IsMassless[sym_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    MemberQ[SARAH`Massless[states], sym];

IsMassless[sym_List, states_:FlexibleSUSY`FSEigenstates] :=
    And @@ (IsMassless /@ sym);

ContainsMassless[sym_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    IsMassless[sym, states];

ContainsMassless[sym_List, states_:FlexibleSUSY`FSEigenstates] :=
    Or @@ (IsMassless[#,states]& /@ sym);

ColorChargedQ[field_] :=
    !FreeQ[FieldInfo[field, includeColourIndices -> True], SARAH`color];

GetColoredParticles[] :=
    Select[GetParticles[], ColorChargedQ];

GetColorRepresentation[SARAH`bar[particle_]] :=
    Module[{rep = GetColorRepresentation[particle]},
           If[rep =!= S && rep =!= O,
              rep = -rep;
             ];
           rep
          ];

GetColorRepresentation[Susyno`LieGroups`conj[particle_]] :=
    Module[{rep = GetColorRepresentation[particle]},
           If[rep =!= S && rep =!= O,
              rep = -rep;
             ];
           rep
          ];

GetColorRepresentation[particle_] :=
    SARAH`getColorRep[particle];

IsQuark[Susyno`LieGroups`conj[sym_]] := IsQuark[sym];
IsQuark[SARAH`bar[sym_]] := IsQuark[sym];
IsQuark[sym_[___]] := IsQuark[sym];
IsQuark[sym_Symbol] := MemberQ[GetColoredParticles[], sym];

IsLepton[Susyno`LieGroups`conj[sym_]] := IsLepton[sym];
IsLepton[SARAH`bar[sym_]] := IsLepton[sym];
IsLepton[sym_[___]] := IsLepton[sym];
IsLepton[sym_Symbol] :=
    MemberQ[Complement[GetParticles[], GetColoredParticles[]], sym] && IsFermion[sym] && IsSMParticle[sym];

IsSMChargedLepton[Susyno`LieGroups`conj[sym_]] := IsSMChargedLepton[sym];
IsSMChargedLepton[SARAH`bar[sym_]] := IsSMChargedLepton[sym];
IsSMChargedLepton[sym_[__]] := IsSMChargedLepton[sym];
IsSMChargedLepton[sym_]     := MemberQ[GetSMChargedLeptons[], sym];

IsSMNeutralLepton[Susyno`LieGroups`conj[sym_]] := IsSMNeutralLepton[sym];
IsSMNeutralLepton[SARAH`bar[sym_]] := IsSMNeutralLepton[sym];
IsSMNeutralLepton[sym_[__]] := IsSMNeutralLepton[sym];
IsSMNeutralLepton[sym_]     := MemberQ[GetSMNeutralLeptons[], sym];

IsSMLepton[Susyno`LieGroups`conj[sym_]] := IsSMLepton[sym];
IsSMLepton[SARAH`bar[sym_]] := IsSMLepton[sym];
IsSMLepton[sym_[__]]        := IsSMLepton[sym];
IsSMLepton[sym_]            := MemberQ[GetSMLeptons[], sym];

IsSMUpQuark[Susyno`LieGroups`conj[sym_]] := IsSMUpQuark[sym];
IsSMUpQuark[SARAH`bar[sym_]] := IsSMUpQuark[sym];
IsSMUpQuark[sym_[__]]       := IsSMUpQuark[sym];
IsSMUpQuark[sym_]           := MemberQ[GetSMUpQuarks[], sym];

IsSMDownQuark[Susyno`LieGroups`conj[sym_]] := IsSMDownQuark[sym];
IsSMDownQuark[SARAH`bar[sym_]] := IsSMDownQuark[sym];
IsSMDownQuark[sym_[__]]     := IsSMDownQuark[sym];
IsSMDownQuark[sym_]         := MemberQ[GetSMDownQuarks[], sym];

IsSMQuark[Susyno`LieGroups`conj[sym_]] := IsSMQuark[sym];
IsSMQuark[SARAH`bar[sym_]] := IsSMQuark[sym];
IsSMQuark[sym_[__]]         := IsSMQuark[sym];
IsSMQuark[sym_]             := MemberQ[GetSMQuarks[], sym];

FSAntiField[p_?IsRealScalar] := p;
FSAntiField[p_?IsComplexScalar] := Susyno`LieGroups`conj[p];
FSAntiField[p_?IsMajoranaFermion] := p;
FSAntiField[p_?IsDiracFermion] := SARAH`bar[p];
FSAntiField[p_?IsVector] := Susyno`LieGroups`conj[p];
FSAntiField[p_?IsGhost] := Susyno`LieGroups`conj[p];

GetSMChargedLeptons[] :=
    Parameters`GetParticleFromDescription["Leptons", {"Electron","Muon","Tau"}];

GetSMNeutralLeptons[] :=
    Parameters`GetParticleFromDescription["Neutrinos", {"Electron Neutrino","Muon Neutrino","Tau Neutrino"}];

GetSMLeptons[] :=
    Join[GetSMNeutralLeptons[], GetSMChargedLeptons[]];

GetSMUpQuarks[] :=
    Parameters`GetParticleFromDescription["Up-Quarks", {"Up Quark","Charmed Quark","Top Quark"}];

GetSMDownQuarks[] :=
    Parameters`GetParticleFromDescription["Down-Quarks", {"Down Quark","Strange Quark","Bottom Quark"}];

GetSMQuarks[] :=
    Join[GetSMDownQuarks[], GetSMUpQuarks[]];

GetUpQuark[gen_, cConvention_:False] :=
    Module[{fields = GetSMUpQuarks[]},
           Switch[Length[fields],
                  1, fields[[1]][gen - If[cConvention === True, 1, 0]],
                  3, fields[[gen]],
                  _, Print["Error: Number of up quarks != 1 and != 3"]; Null
                 ]
          ];

GetDownQuark[gen_, cConvention_:False] :=
    Module[{fields = GetSMDownQuarks[]},
           Switch[Length[fields],
                  1, fields[[1]][gen - If[cConvention === True, 1, 0]],
                  3, fields[[gen]],
                  _, Print["Error: Number of down quarks != 1 and != 3"]; Null
                 ]
          ];

GetUpLepton[gen_, cConvention_:False] :=
    Module[{fields = GetSMNeutralLeptons[]},
           Switch[Length[fields],
                  1, fields[[1]][gen - If[cConvention === True, 1, 0]],
                  3, fields[[gen]],
                  _, Print["Error: Number of up leptons != 1 and != 3"]; Null
                 ]
          ];

GetDownLepton[gen_, cConvention_:False] :=
    Module[{fields = GetSMChargedLeptons[]},
           Switch[Length[fields],
                  1, fields[[1]][gen - If[cConvention === True, 1, 0]],
                  3, fields[[gen]],
                  _, Print["Error: Number of down leptons != 1 and != 3"]; Null
                 ]
          ];

GetSMUpQuark[]        := GetUpQuark[1];
GetSMCharmQuark[]     := GetUpQuark[2];
GetSMTopQuark[]       := GetUpQuark[3];
GetSMDownQuark[]      := GetDownQuark[1];
GetSMStrangeQuark[]   := GetDownQuark[2];
GetSMBottomQuark[]    := GetDownQuark[3];
GetSMElectronLepton[] := GetDownLepton[1];
GetSMMuonLepton[]     := GetDownLepton[2];
GetSMTauLepton[]      := GetDownLepton[3];
GetSMNeutrino1[]      := GetUpLepton[1];
GetSMNeutrino2[]      := GetUpLepton[2];
GetSMNeutrino3[]      := GetUpLepton[3];

GetSMTopQuarkMultiplet[]       := GetUpQuark[3]    /. head_[_] :> head;
GetSMStrangeQuarkMultiplet[]    := GetDownQuark[2]  /. head_[_] :> head;
GetSMBottomQuarkMultiplet[]    := GetDownQuark[3]  /. head_[_] :> head;
GetSMElectronLeptonMultiplet[] := GetDownLepton[1] /. head_[_] :> head;
GetSMMuonLeptonMultiplet[]     := GetDownLepton[2] /. head_[_] :> head;
GetSMTauLeptonMultiplet[]      := GetDownLepton[3] /. head_[_] :> head;

GetMass[particle_[idx__]] := GetMass[particle][idx];
GetMass[particle_Symbol] := FlexibleSUSY`M[particle];

GetElectricCharge[p_] :=
    Module[{charge},
           If[p === SARAH`AntiField[p],
              charge = 0;,
              charge = SARAH`getElectricCharge[p];
              If[!NumericQ[charge],
                 charge = Cases[-I SARAH`Vertex[{SARAH`AntiField[p], p, SARAH`VectorP},
                                                UseDependences -> True][[2,1]], _?NumberQ];
                 If[charge === {},
                    charge = 0;,
                    charge = First[charge];
                   ];
                ];
             ];
           charge
          ];

(* Returns list of pairs {p,v}, where p is the given golstone
   boson and v is the corresponding vector boson.

   Example (MSSM):
     GetCorrespondingVectorBosons[Ah]      ->  {{Ah[1], VZ}}
     GetCorrespondingVectorBosons[Ah[1]]   ->  {{Ah[1], VZ}}
     GetCorrespondingVectorBosons[Ah[{1}]] ->  {{Ah[1], VZ}}
*)
GetCorrespondingVectorBosons[goldstone_[idx_Integer]] :=
    GetCorrespondingVectorBosons[goldstone[{idx}]];

GetCorrespondingVectorBosons[goldstone_[idx_Symbol]] :=
    GetCorrespondingVectorBosons[goldstone[{idx}]];

GetCorrespondingVectorBosons[goldstone_] :=
    Module[{vector, idx, sym, association},
           association = Cases[SARAH`GoldstoneGhost, {vector_, goldstone | goldstone[{idx_}]}];
           Reverse /@ association /. sym_[{idx_}] :> sym[idx]
          ];

GetGoldstoneBosons[] :=
    Transpose[SARAH`GoldstoneGhost][[2]];

GetSMGoldstoneBosons[] :=
    Cases[SARAH`GoldstoneGhost, {vector_?IsSMParticle, goldstone_} :> goldstone];

GetDimension[sym_List, states_:FlexibleSUSY`FSEigenstates] :=
    Plus @@ (GetDimension[#, states]& /@ sym);

GetDimension[(SARAH`bar|Susyno`LieGroups`conj)[sym_], states_:FlexibleSUSY`FSEigenstates] :=
    GetDimension[sym, states];

GetDimension[sym_[__], states_:FlexibleSUSY`FSEigenstates] := GetDimension[sym, states];

GetDimension[sym_Symbol, states_:FlexibleSUSY`FSEigenstates] :=
    SARAH`getGen[sym, states];

GetDimensionStartSkippingGoldstones[sym_[__]] :=
    GetDimensionStartSkippingGoldstones[sym];

GetDimensionStartSkippingGoldstones[(SARAH`bar|Susyno`LieGroups`conj)[sym_]] :=
    GetDimensionStartSkippingGoldstones[sym];

GetDimensionStartSkippingGoldstones[sym_, goldstoneGhost_] :=
    Module[{goldstones, max = 1},
           goldstones = Transpose[goldstoneGhost][[2]];
           If[FreeQ[goldstones, sym],
              Return[1];,
              If[GetDimension[sym] === 1,
                 Return[2];,
                 While[!FreeQ[goldstones, sym[{max}]],
                       max++];
                 Return[max];
                ];
             ];
          ];

GetDimensionStartSkippingGoldstones[sym_] :=
    GetDimensionStartSkippingGoldstones[sym, SARAH`GoldstoneGhost];

GetDimensionStartSkippingSMGoldstones[sym_] :=
    GetDimensionStartSkippingGoldstones[sym, Cases[SARAH`GoldstoneGhost, {_?IsSMParticle, _}]];

GetDimensionWithoutGoldstones[sym_[__], states_:FlexibleSUSY`FSEigenstates] :=
    GetDimensionWithoutGoldstones[sym, states];

GetDimensionWithoutGoldstones[(SARAH`bar|Susyno`LieGroups`conj)[sym_], states_:FlexibleSUSY`FSEigenstates] :=
    GetDimensionWithoutGoldstones[sym, states];

GetDimensionWithoutGoldstones[sym_, states_:FlexibleSUSY`FSEigenstates] :=
    Module[{dim, numberOfGoldstones},
           numberOfGoldstones = GetDimensionStartSkippingGoldstones[sym] - 1;
           dim = GetDimension[sym] - numberOfGoldstones;
           If[dim <= 0, 0, dim]
          ];

GetParticleIndices[sym_[__]] :=
    GetParticleIndices[sym];

GetParticleIndices[sym_] :=
    Module[{result},
           result = Cases[SARAH`Particles[SARAH`EWSB], {sym, __, indexList_} :> indexList];
           If[Length[result] > 0,
              result = result[[1]]
             ];
           result
          ];

DimOf[CConversion`ScalarType[CConversion`realScalarCType]] := 1;
DimOf[CConversion`VectorType[CConversion`realScalarCType, dim]] := dim;
DimOf[t_] := (Print["Unknown type: ", t]; Quit[1]);

GetMassType[FlexibleSUSY`M[particles_List]] :=
    CConversion`VectorType[CConversion`realScalarCType,
                           Plus @@ (DimOf /@ (GetMassType /@ particles))];

GetMassType[FlexibleSUSY`M[particle_]] := GetMassType[particle];

GetMassType[particle_] :=
    Module[{dim = GetDimension[particle]},
           If[dim == 1,
              CConversion`ScalarType[CConversion`realScalarCType],
              CConversion`VectorType[CConversion`realScalarCType, dim]
             ]
          ];

GetMassMatrixType[particle_] :=
    Module[{dim = GetDimension[particle]},
           If[dim == 1,
              If[Parameters`AllModelParametersAreReal[],
                 CConversion`ScalarType[CConversion`realScalarCType]
                 ,
                 Which[IsFermion[particle],
                       CConversion`ScalarType[CConversion`complexScalarCType],
                       True,
                       CConversion`ScalarType[CConversion`realScalarCType]
                      ]
                ]
              ,
              If[Parameters`AllModelParametersAreReal[],
                 CConversion`MatrixType[CConversion`realScalarCType, dim, dim]
                 ,
                 Which[IsFermion[particle],
                       CConversion`MatrixType[CConversion`complexScalarCType, dim, dim],
                       IsRealScalar[particle],
                       CConversion`MatrixType[CConversion`realScalarCType, dim, dim],
                       IsComplexScalar[particle],
                       CConversion`MatrixType[CConversion`complexScalarCType, dim, dim],
                       IsVector[particle],
                       CConversion`MatrixType[CConversion`realScalarCType, dim, dim],
                       True,
                       Print["Error: GetMassMatrixType: unknown particle type: ", particle];
                       Quit[1];
                      ]
                ]
             ]
          ];

GetMixingMatrixType[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{type, eigenstate, mixingMatrixSymbol, dim},
           eigenstate = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = Length[GetMassMatrix[massMatrix]];
           Which[IsRealScalar[eigenstate],
                 type = CConversion`realScalarCType;,
                 IsScalar[eigenstate] && Parameters`AllModelParametersAreReal[],
                 type = CConversion`realScalarCType;,
                 IsVector[eigenstate],
                 type = CConversion`realScalarCType;,
                 True,
                 type = CConversion`complexScalarCType;
                ];
           Return[CConversion`MatrixType[type, dim, dim]];
          ];

(* Removes generators and Delta with the given indices.
 * Especially the following replacements are done:
 *
 * SARAH`Lam[__] -> 2    corresponds to SU(3) generator T[__] -> 1
 * SARAH`Sig[__] -> 2    corresponds to SU(2) generator T[__] -> 1
 * SARAH`fSU2[__] -> 1   SU(2) structure function
 * SARAH`fSU3[__] -> 1   SU(3) structure function
 * SARAH`Delta[a,b] -> 1 where a and b are from the `indices' list
 *)
StripGenerators[expr_, indices_List] :=
    Module[{headers = {SARAH`Delta}, h, indexCombinations, removeSymbols = {}},
           indexCombinations = DeleteCases[Subsets[indices],{}];
           For[h = 1, h <= Length[headers], h++,
               AppendTo[removeSymbols, (headers[[h]] @@ #)& /@ indexCombinations];
              ];
           removeSymbols = (Rule[#, 1])& /@ Flatten[removeSymbols];
           removeSymbols = Join[removeSymbols,
                                { Rule[SARAH`Lam[__], 2],
                                  Rule[SARAH`Sig[__], 2],
                                  Rule[SARAH`fSU2[__], 1],
                                  Rule[SARAH`fSU3[__], 1] }];
           expr /. removeSymbols
          ];

FindMixingMatrixSymbolFor[particle_Symbol] :=
    Module[{k, i, l, mixingMatrixSymbol = Null, mixingList = {}, mixingScheme = {},
            currentName},
           For[k = 1, k <= Length[NameOfStates], k++,
               If[Head[DEFINITION[NameOfStates[[k]]][MatterSector]] === List,
                  mixingList = DEFINITION[NameOfStates[[k]]][MatterSector];
                  For[i = 1, i <= Length[mixingList], i++,
                      If[Length[mixingList[[i]]] != 2, Continue[]];
                      mixingScheme = mixingList[[i,2]];
                      If[Length[mixingScheme] != 2, Continue[]];
                      If[Head[mixingScheme[[1]]] === Symbol,
                         currentName = mixingScheme[[1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                         If[currentName === particle,
                            mixingMatrixSymbol = mixingScheme[[2]];
                            Return[mixingMatrixSymbol];
                           ];
                         ,
                         currentName = mixingScheme[[1,1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                         If[currentName === particle,
                            mixingMatrixSymbol = {mixingScheme[[1,2]], mixingScheme[[2,2]]};
                            Return[mixingMatrixSymbol];
                           ];
                        ];
                     ];
                 ];
               If[Head[DEFINITION[NameOfStates[[k]]][GaugeSector]] === List,
                  mixingList = DEFINITION[NameOfStates[[k]]][GaugeSector];
                  For[i = 1, i <= Length[mixingList], i++,
                      mixingScheme = mixingList[[i]];
                      If[Length[mixingScheme] != 3, Continue[]];
                      If[MemberQ[mixingScheme[[2]], particle],
                         mixingMatrixSymbol = mixingScheme[[3]];
                         Return[mixingMatrixSymbol];
                        ];
                     ];
                 ];
              ];
           Return[mixingMatrixSymbol];
          ];

IsUnmixed[particle_Symbol] :=
    MemberQ[(#[[1]])& /@ unrotatedParticles, particle];

GetMassOfUnmixedParticle[particle_Symbol] :=
    Cases[unrotatedParticles, {particle, _}][[1,2]];

FindMassEigenstateForMixingMatrix[mixingMatrixSymbol_Symbol] :=
    Module[{k, i, l, particle, mixingList = {}, mixingScheme = {},
            currentName},
           For[k = 1, k <= Length[NameOfStates], k++,
               If[Head[DEFINITION[NameOfStates[[k]]][MatterSector]] === List,
                  mixingList = DEFINITION[NameOfStates[[k]]][MatterSector];
                  For[i = 1, i <= Length[mixingList], i++,
                      If[Length[mixingList[[i]]] != 2, Continue[]];
                      mixingScheme = mixingList[[i,2]];
                      If[Length[mixingScheme] != 2, Continue[]];
                      If[mixingScheme[[2]] === mixingMatrixSymbol,
                         particle = mixingScheme[[1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                         Return[particle];
                        ];
                      If[Head[mixingScheme[[1]]] === List,
                         If[mixingScheme[[1,2]] === mixingMatrixSymbol,
                            particle = mixingScheme[[1,1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                            Return[particle];
                           ];
                        ];
                      If[Head[mixingScheme[[2]]] === List,
                         If[mixingScheme[[2,2]] === mixingMatrixSymbol,
                            particle = mixingScheme[[2,1]] /.
                            SARAH`diracSubBack1[NameOfStates[[k]]] /.
                            SARAH`diracSubBack2[NameOfStates[[k]]];
                            Return[particle];
                           ];
                        ];
                     ];
                 ];
               If[Head[DEFINITION[NameOfStates[[k]]][GaugeSector]] === List,
                  mixingList = DEFINITION[NameOfStates[[k]]][GaugeSector];
                  For[i = 1, i <= Length[mixingList], i++,
                      mixingScheme = mixingList[[i]];
                      If[Length[mixingScheme] != 3, Continue[]];
                      If[mixingScheme[[3]] === mixingMatrixSymbol,
                         Return[mixingScheme[[2]]];
                        ];
                     ];
                 ];
              ];
           Null
          ];

DeleteDuplicateSinglets[massMatrices_List] :=
    Module[{result = massMatrices, i, m, me, mm, p, other},
           (* delete duplicates of the form

              FSMassMatrix[_, {me}, _]
              FSMassMatrix[_,  me , _]   <--  delete

              but keep mixing matrix from one or the other
            *)
           For[i = 1, i <= Length[massMatrices], i++,
               me = GetMassEigenstate[massMatrices[[i]]];
               mm = GetMixingMatrixSymbol[massMatrices[[i]]];
               If[Head[me] =!= List,
                  other = Cases[result, p:TreeMasses`FSMassMatrix[_, {me}, _] :> p];
                  If[other === {}, Continue[]];
                  other = other[[1]];
                  result = DeleteCases[result, TreeMasses`FSMassMatrix[_, {me}, _]];
                  If[Head[GetMassEigenstate[other]] === List &&
                     Length[GetMassEigenstate[other]] > 0 &&
                     GetMixingMatrixSymbol[other] =!= Null,
                     result = result /.
                         TreeMasses`FSMassMatrix[m_, me, Null] :> TreeMasses`FSMassMatrix[m, me, GetMixingMatrixSymbol[other]];
                    ];
                 ];
              ];
           result
          ];

GetIntermediateMassMatrices[massMatrices_List] :=
    Module[{intermediatePars, massEigenstates},
           CreateMMs[{massEigenstate_, mixingMatrix_}] :=
               Module[{massMatrix},
                      If[Head[massEigenstate] === List,
                         massMatrix = ReplaceDependencies[SARAH`MassMatrix[massEigenstate[[1]]]];,
                         massMatrix = ReplaceDependencies[SARAH`MassMatrix[massEigenstate]];
                        ];
                      If[Head[massMatrix] === List,
                         TreeMasses`FSMassMatrix[massMatrix, massEigenstate, mixingMatrix],
                         Null
                        ]
                     ];
           intermediatePars = DeleteCases[
               Parameters`GetIntermediateOutputParameterDependencies[GetMassMatrix /@ massMatrices], _?NumberQ];
           massEigenstates = (FindMassEigenstateForMixingMatrix /@ intermediatePars) /. Susyno`LieGroups`conj[a_] :> a;
           CreateMMs /@ DeleteCases[Utils`Zip[massEigenstates, intermediatePars], {Null, _}]
          ];

ConvertSarahMassMatrices[] :=
    Module[{particles = {}, result = {}, eigenstateName, massMatrix,
            k, mixingMatrixSymbol},
           (* the ghost masses will later be replaced explicitely by the
              corresponding vector boson masses *)
           particles = Select[GetParticles[], (!IsGhost[#])&];
           For[k = 1, k <= Length[particles], k++,
               eigenstateName = particles[[k]];
               If[IsUnmixed[eigenstateName],
                  massMatrix = { ReplaceDependencies[GetMassOfUnmixedParticle[eigenstateName]] };
                  AppendTo[result, TreeMasses`FSMassMatrix[massMatrix, eigenstateName, Null]];
                  ,
                  massMatrix = ReplaceDependencies[SARAH`MassMatrix[eigenstateName]];
                  mixingMatrixSymbol = FindMixingMatrixSymbolFor[eigenstateName];
                  If[Head[massMatrix] === List,
                     AppendTo[result, TreeMasses`FSMassMatrix[massMatrix, eigenstateName, mixingMatrixSymbol]];
                    ];
                 ];
              ];
           (* append mass matrix for intermediate output parameters *)
           result = DeleteDuplicateSinglets[Join[result, GetIntermediateMassMatrices[result]]];
           result = DeleteDuplicateVectors[result];
           Return[result];
          ];

(* returns masses of unmixed particles *)
GetUnmixedParticleMasses[] :=
    Module[{particles = {}, result = {}, eigenstateName, massMatrix, k},
           particles = Select[GetParticles[], (!IsGhost[#] && IsUnmixed[#])&];
           For[k = 1, k <= Length[particles], k++,
               eigenstateName = particles[[k]];
               massMatrix = { ReplaceDependencies[GetMassOfUnmixedParticle[eigenstateName]] };
               AppendTo[result, TreeMasses`FSMassMatrix[massMatrix, eigenstateName, Null]];
              ];
           result
          ];

GetMixingMatrixSymbol[massMatrix_TreeMasses`FSMassMatrix] := massMatrix[[3]];

GetMassEigenstate[massMatrix_TreeMasses`FSMassMatrix] := massMatrix[[2]];

GetMassMatrix[massMatrix_TreeMasses`FSMassMatrix] := massMatrix[[1]];

MakeESSymbol[p_List] := Symbol[StringJoin[ToString /@ p]];
MakeESSymbol[FlexibleSUSY`M[p_List]] := FlexibleSUSY`M[MakeESSymbol[p]];
MakeESSymbol[p_] := p;

CreateMassGetter[p:TreeMasses`FSMassMatrix[_,massESSymbols_List,_], postFix_String:"", wrapper_String:""] :=
    Module[{massMatrices},
           massMatrices = DeleteDuplicates[TreeMasses`FSMassMatrix[0, #, Null]& /@ massESSymbols];
           StringJoin[CreateMassGetter[#,postFix,wrapper]& /@ massMatrices]
          ];

CreateMassGetter[massMatrix_TreeMasses`FSMassMatrix, postFix_String:"", wrapper_String:""] :=
    Module[{massESSymbol, returnType, dim, dimStr, massESSymbolStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           massESSymbolStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[MakeESSymbol[massESSymbol]]];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           If[dim == 1,
              returnType = CConversion`ScalarType[CConversion`realScalarCType];,
              returnType = CConversion`ArrayType[CConversion`realScalarCType, dim];
             ];
           CConversion`CreateInlineGetters[massESSymbolStr, massESSymbolStr, returnType, postFix, wrapper]
          ];

CreateSLHAPoleMassGetter[massMatrix_TreeMasses`FSMassMatrix] :=
    CreateMassGetter[massMatrix, "_pole_slha", "PHYSICAL_SLHA"];

CreateParticleEnum[particles_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CConversion`ToValidCSymbolString /@ particles, ", "];
           If[Length[particles] > 0, result = result <> ", ";];
           "enum Particles : int { " <> result <> "NUMBER_OF_PARTICLES };\n"
          ];

DecomposeParticle[particle_] :=
    If[GetDimension[particle] == 1,
       { particle },
       Array[particle, GetDimension[particle]]
      ];

CreateParticleMassEnumName[particle_[idx_]] :=
    CreateParticleMassEnumName[particle] <> "_" <> ToString[idx];

CreateParticleMassEnumName[particle_] :=
    "M" <> CConversion`ToValidCSymbolString[particle];

CreateParticleMassEnum[particles_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[CreateParticleMassEnum /@ particles, ", "];
           If[Length[particles] > 0, result = result <> ", ";];
           "enum Masses : int { " <> result <> "NUMBER_OF_MASSES };\n"
          ];

CreateParticleMassEnum[p_] :=
    Utils`StringJoinWithSeparator[CreateParticleMassEnumName /@ DecomposeParticle[p], ", "];

DecomposeMixingMatrix[mm_List, type_] := { #, type }& /@ mm;
DecomposeMixingMatrix[mm_    , type_] := {{ mm, type }};

GetMixingMatrixAndTypeFrom[mixing_] :=
    DecomposeMixingMatrix[GetMixingMatrixSymbol[mixing], GetMixingMatrixType[mixing]];

GetMixingMatricesAndTypesFrom[mixings_List] :=
    Join @@ (GetMixingMatrixAndTypeFrom /@ mixings);

CreateParticleMixingEnum[mixings_List] :=
    Module[{nonNullMixings, result},
           nonNullMixings = Select[mixings, (GetMixingMatrixSymbol[#] =!= Null)&];
           result = Utils`StringJoinWithSeparator[
               Parameters`CreateParameterEnums[#[[1]], #[[2]]]& /@ GetMixingMatricesAndTypesFrom[nonNullMixings], ", "];
           If[Length[nonNullMixings] > 0, result = result <> ", ";];
           "enum Mixings : int { " <> result <> "NUMBER_OF_MIXINGS };\n"
          ];

SARAHNameStr[p_] :=
    "\"" <> CConversion`ToValidCSymbolString[p] <> "\"";

CreateParticleNames[particles_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[SARAHNameStr /@ particles, ", "];
           "const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {" <>
           result <> "};\n"
          ];

CreateParticleMultiplicity[particles_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[GetDimension /@ particles, ", "];
           "const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {" <>
           result <> "};\n"
          ];

TeXNameStr[p_] :=
    "\"" <> StringReplace[SARAH`getLaTeXField[p], "\\" -> "\\\\"] <> "\"";

CreateParticleLaTeXNames[particles_List] :=
    Module[{result},
           result = Utils`StringJoinWithSeparator[TeXNameStr /@ particles, ", "];
           "const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {" <>
           IndentText[result] <> "};\n"
          ];

CreateParticleMixingNames[mixings_List] :=
    Module[{nonNullMixings, result},
           nonNullMixings = Select[mixings, (GetMixingMatrixSymbol[#] =!= Null)&];
           result = Utils`StringJoinWithSeparator[
               Parameters`CreateParameterSARAHNames[#[[1]], #[[2]]]& /@ GetMixingMatricesAndTypesFrom[nonNullMixings], ", "];
           "const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {" <>
           IndentText[result] <> "};\n"
          ];

(* Note that "bar" and "conj" get turned into bar<...>::type and
 conj<...>::type respectively! *)
CreateFieldClassName[p_, OptionsPattern[{prefixNamespace -> False}]] :=
    If[StringQ[OptionValue[prefixNamespace]],
       OptionValue[prefixNamespace] <> "::",
       ""] <> SymbolName[p];

CreateFieldClassName[SARAH`bar[p_], OptionsPattern[{prefixNamespace -> False}]] :=
    "typename field_traits::bar<" <> CreateFieldClassName[p, prefixNamespace -> OptionValue[prefixNamespace]] <> ">::type";

CreateFieldClassName[Susyno`LieGroups`conj[p_], OptionsPattern[{prefixNamespace -> False}]] :=
    "typename field_traits::conj<" <> CreateFieldClassName[p, prefixNamespace -> OptionValue[prefixNamespace]] <> ">::type";

FillSpectrumVector[particles_List] :=
    Module[{par, parStr, massStr, latexName, result = ""},
           For[i = 1, i <= Length[particles], i++,
               par = particles[[i]];
               parStr = CConversion`ToValidCSymbolString[par];
               massStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[par]];
               latexName = StringReplace[SARAH`getLaTeXField[par], "\\" -> "\\\\"];
               result = result <> "spectrum.emplace_back(TParticle(\"" <> parStr <>
                        "\", \"" <> latexName <> "\", to_valarray(PHYSICAL(" <>
                        massStr <> "))));\n";
              ];
           Return[result];
          ];

CreateMixingMatrixGetter[massMatrix_TreeMasses`FSMassMatrix, postFix_String:"", wrapper_String:""] :=
    Module[{mixingMatrixSymbol, returnType},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           returnType = GetMixingMatrixType[massMatrix];
           CreateMixingMatrixGetter[mixingMatrixSymbol, returnType, postFix, wrapper]
          ];

CreateMixingMatrixGetter[mixingMatrixSymbol_List, returnType_, postFix_String:"", wrapper_String:""] :=
    Module[{result = ""},
           (result = result <> CreateMixingMatrixGetter[#,returnType, postFix, wrapper])& /@ mixingMatrixSymbol;
           Return[result];
          ];

CreateMixingMatrixGetter[Null, returnType_, postFix_String:"", wrapper_String:""] := "";

CreateMixingMatrixGetter[mixingMatrixSymbol_Symbol, returnType_, postFix_String:"", wrapper_String:""] :=
    CConversion`CreateInlineGetters[CConversion`ToValidCSymbolString[mixingMatrixSymbol],
                                    CConversion`ToValidCSymbolString[mixingMatrixSymbol],
                                    returnType, postFix, wrapper];

CreateSLHAPoleMixingMatrixGetter[massMatrix_TreeMasses`FSMassMatrix /; GetMixingMatrixSymbol[massMatrix] === Null] := "";

CreateSLHAPoleMixingMatrixGetter[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{mixingMatrixSymbol, particle, dim, returnType},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           particle = GetMassEigenstate[massMatrix];
           dim = Length[GetMassMatrix[massMatrix]];
           returnType = GetMixingMatrixType[massMatrix];
           (* mixing matrices for majorana Fermions and
              2-component charginos will be made real *)
           If[IsMajoranaFermion[particle] || (IsChargino[particle] && dim <= 2),
              Utils`ApplyAndConcatenate[
                  Function[m,
                           CConversion`CreateInlineGetter[
                               CConversion`ToValidCSymbolString[m],
                               CConversion`ToValidCSymbolString[m],
                               returnType, "_pole_slha", "PHYSICAL_SLHA"] <>
                           CConversion`CreateInlineElementGetter[
                               CConversion`ToValidCSymbolString[m],
                               CConversion`ToValidCSymbolString[m],
                               CConversion`ToRealType[returnType], "_pole_slha", "PHYSICAL_SLHA_REAL"]],
                  mixingMatrixSymbol
              ]
              ,
              CreateMixingMatrixGetter[massMatrix, "_pole_slha", "PHYSICAL_SLHA"]
             ]
          ];

CreateFSMassMatrixForUnmixedParticle[TreeMasses`FSMassMatrix[expr_, massESSymbol_, Null]] :=
    Module[{matrix, dim},
           dim = GetDimension[massESSymbol];
           If[dim == 1,
              matrix = expr;
              ,
              matrix = Table[expr /. List -> Identity,
                             {SARAH`gt1, 1, dim}, {SARAH`gt2, 1, dim}];
             ];
           TreeMasses`FSMassMatrix[matrix, massESSymbol, Null]
          ];

CreateMassCalculationPrototype[m:TreeMasses`FSMassMatrix[expr_, massESSymbol_, Null]] :=
    Module[{result, ev = CConversion`ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]],
            massMatrix},
           result = "void calculate_" <> ev <> "();\n";
           massMatrix = CreateFSMassMatrixForUnmixedParticle[m];
           result = CreateMassMatrixGetterPrototype[massMatrix] <>
                    result;
           Return[result];
          ];

CreateMassCalculationPrototype[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           result = CreateMassMatrixGetterPrototype[massMatrix] <>
                    "void calculate_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[MakeESSymbol[massESSymbol]]] <>
                    "();\n";
           Return[result];
          ];

BubbleSort[l_List, Pred_] :=
    Module[{i,k,t,ltemp = l},
           For[i = 1, i <= Length[ltemp], i++,
               For[k = 1, k <= Length[ltemp], k++,
                   If[!Pred[ltemp[[i]], ltemp[[k]]],
                      t = ltemp[[i]];
                      ltemp[[i]] = ltemp[[k]];
                      ltemp[[k]] = t;
                     ];
                  ];
              ];
           ltemp
          ];

(* delete duplicates of the form

   FSMassMatrix[_, {m1, m2}, _]
   FSMassMatrix[_,  m1     , _]   <--  delete
   FSMassMatrix[_,  m2     , _]   <--  delete
 *)
DeleteDuplicateVectors[massMatrices_List] :=
    Module[{result = massMatrices, i, m, me},
           For[i = 1, i <= Length[massMatrices], i++,
               me = GetMassEigenstate[massMatrices[[i]]];
               If[Head[me] === List,
                  result = DeleteCases[result, TreeMasses`FSMassMatrix[_, m_ /; MemberQ[me, m], _]];
                 ];
              ];
           result
          ];

CallMassCalculationFunctions[massMatrices_List] :=
    Module[{result = "", k, sortedMassMatrices, matrix, PredVectorsFirst},
           (* Predicate function which returns false if the mass matrix
              of m1 depends on the mixing matrix of m2.  True otherwise. *)
           PredVectorsFirst[m1_TreeMasses`FSMassMatrix, m2_TreeMasses`FSMassMatrix] :=
               Module[{mm1, z2},
                      mm1 = GetMassMatrix[m1] /. Parameters`GetDependenceSPhenoRules[];
                      z2  = GetMixingMatrixSymbol[m2];
                      If[Head[z2] === List && Length[z2] == 2,
                         FreeQ[mm1, z2[[1]]] && FreeQ[mm1, z2[[2]]],
                         FreeQ[mm1, z2]
                        ]
                     ];
           (* Sort mass matrices such that vector boson masses get
              calculated first.  This is necessary because the later
              calculated masses might depend on some SM mixing angles,
              as ThetaW. *)
           (* Note: Due to the chosen predicate, a stable bubble sort must be used *)
           sortedMassMatrices = Reverse @ BubbleSort[massMatrices, PredVectorsFirst];
           For[k = 1, k <= Length[sortedMassMatrices], k++,
               matrix = sortedMassMatrices[[k]];
               result = result <> CallMassCalculationFunction[matrix];
              ];
           Return[result];
          ];

CallMassCalculationFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", k, massESSymbol},
           massESSymbol = GetMassEigenstate[massMatrix];
           result = "calculate_" <> CConversion`ToValidCSymbolString[FlexibleSUSY`M[MakeESSymbol[massESSymbol]]]
                    <> "();\n";
           Return[result];
          ];

IsSymmetric[matrix_List] := IsHermitian[matrix, Identity];

IsHermitian[matrix_List, op_:Susyno`LieGroups`conj] :=
    Module[{rows, cols, i, k, difference},
           rows = Length[matrix];
           For[i = 1, i <= rows, i++,
               cols = Length[matrix[[i]]];
               If[rows =!= cols, Return[False];];
               For[k = 1, k <= i, k++,
                   difference = matrix[[i,k]] - op[matrix[[k,i]]] //. {
                       Susyno`LieGroups`conj[SARAH`sum[ind_,a_,b_,expr_]] :>
                           SARAH`sum[ind,a,b,Susyno`LieGroups`conj[expr]],
                       Susyno`LieGroups`conj[m_[a_,b_]] :>
                           m[b,a] /; MemberQ[SARAH`ListSoftBreakingScalarMasses, m]
                   };
                   If[!PossibleZeroQ[difference], Return[False];];
                  ];
              ];
           Return[True];
          ];

MatrixToCFormString[matrix_List /; MatrixQ[matrix], symbol_String] :=
    Module[{dim = Length[matrix], result = "", i, k},
           For[i = 1, i <= dim, i++,
               For[k = 1, k <= dim, k++,
                   result = result <> symbol <> "(" <> ToString[i-1] <>
                            "," <> ToString[k-1] <> ") = " <>
                            RValueToCFormString[matrix[[i,k]]] <> ";\n";
                  ];
              ];
           Return[result];
          ];

(* fill upper triangle of matrix *)
UpperTriangleMatrixToCFormString[matrix_List /; MatrixQ[matrix], symbol_String] :=
    Module[{dim = Length[matrix], result = "", i, k},
           For[i = 1, i <= dim, i++,
               For[k = i, k <= dim, k++,
                   result = result <> symbol <> "(" <> ToString[i-1] <>
                            "," <> ToString[k-1] <> ") = " <>
                            RValueToCFormString[matrix[[i,k]]] <> ";\n";
                  ];
              ];
           Return[result];
          ];

CastMatrixElementToReal[el_] :=
    If[Or @@ (Parameters`IsComplexParameter /@ Parameters`FindAllParameters[el]),
       Re[el],
       el
      ];

MatrixToCFormString[matrix_List /; MatrixQ[matrix], symbol_String, matrixElementType_] :=
    Module[{dim, result = "", i, k, isSymmetric = IsSymmetric[matrix],
            isHermitian = IsHermitian[matrix], matrixType, dimStr, castedMatrix},
           dim = Length[matrix];
           dimStr = ToString[dim];
           matrixType = CreateCType[CConversion`MatrixType[matrixElementType, dim, dim]];
           castedMatrix = If[matrixElementType === CConversion`realScalarCType,
                             CastMatrixElementToReal /@ matrix, matrix];
           result = matrixType <> " " <> symbol <> ";\n\n"; (* not initialized *)
           Which[isSymmetric,
                 result = result <> UpperTriangleMatrixToCFormString[castedMatrix, symbol] <> "\n" <>
                          "Symmetrize(" <> symbol <> ");\n";,
                 isHermitian,
                 result = result <> UpperTriangleMatrixToCFormString[castedMatrix, symbol] <> "\n" <>
                          "Hermitianize(" <> symbol <> ");\n";,
                 True,
                 result = result <> MatrixToCFormString[castedMatrix, symbol];
                ];
           Return[result];
          ];

MatrixToCFormString[matrix_List, symbol_String, matrixElementType_] :=
    Module[{result = "", type, ctype},
           If[Length[matrix] != 1,
              Print["Error: Expression is not a 1-element list: ", matrix];
              Return["return 0.;"];
             ];
           type = CConversion`ScalarType[matrixElementType];
           ctype = CreateCType[type];
           result = "const " <> ctype <> " " <> symbol <> " = " <>
                    CastTo[RValueToCFormString[matrix[[1]]], type] <>
                    ";\n"; (* not initialized *)
           Return[result];
          ];

CreateMassMatrixGetterFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, body, ev, matrixSymbol, matrix, massESSymbol,
            inputParsDecl, matrixType, matrixElementType, dim, dimStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           ev = CConversion`ToValidCSymbolString[CConversion`GetHead[MakeESSymbol[massESSymbol]]];
           matrixSymbol = "mass_matrix_" <> ev;
           matrix = GetMassMatrix[massMatrix];
           (* Remove color SU(3) generators, structure functions and
              Kronecker delta with color indices.
              Note: ct1 ... ct4 are reserved SU(3) color indices of
              the fundamental representation of SU(3) in SARAH.
           *)
           matrix = StripGenerators[matrix,
                                    {SARAH`ct1, SARAH`ct2, SARAH`ct3, SARAH`ct4}];
           dim = Length[matrix];
           dimStr = ToString[dim];
           (* convert 1-dimensional matrix to scalar *)
           If[dim == 1 && MatrixQ[matrix],
              matrix = matrix[[1]];
             ];
           matrixType = GetMassMatrixType[massESSymbol];
           matrixElementType = CConversion`GetElementType[matrixType];
           matrixType = CreateCType[matrixType];
           inputParsDecl = Parameters`CreateLocalConstRefsForInputParameters[matrix, "LOCALINPUT"];
           body = inputParsDecl <> "\n" <> MatrixToCFormString[matrix, matrixSymbol, matrixElementType] <> "\n";
           result = matrixType <> " CLASSNAME::get_" <> matrixSymbol <> "() const\n{\n" <>
                    IndentText[body] <>
                    "return " <> matrixSymbol <> ";\n}\n";
           Return[result];
          ];

CreateMassMatrixGetterPrototype[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, ev, matrixSymbol, matrix, massESSymbol, matrixType,
            dim, dimStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           ev = CConversion`ToValidCSymbolString[CConversion`GetHead[MakeESSymbol[massESSymbol]]];
           matrix = GetMassMatrix[massMatrix];
           dim = Length[matrix];
           dimStr = ToString[dim];
           matrixType = GetMassMatrixType[massESSymbol];
           matrixType = CreateCType[matrixType];
           matrixSymbol = "mass_matrix_" <> ev;
           result = matrixType <> " get_" <> matrixSymbol <> "() const;\n";
           Return[result];
          ];

WrapMacro[sym_String, ""] := sym;
WrapMacro[sym_String, macro_String] := macro <> "(" <> sym <> ")";

CreateReorderingFunctionCalls[{idx_Integer, vector_, higgs_, mixingMatrix_}, macro_String:""] :=
    "move_goldstone_to(" <> ToString[idx-1] <> ", " <>
    WrapMacro[CConversion`ToValidCSymbolString[FlexibleSUSY`M[vector]], ""] <> ", " <>
    WrapMacro[CConversion`ToValidCSymbolString[FlexibleSUSY`M[higgs]], macro] <> ", " <>
    WrapMacro[CConversion`ToValidCSymbolString[mixingMatrix], macro] <> ");\n";

CreateReorderingFunctionCalls[___] := "";

ReorderGoldstoneBosons[particle_, mixingMatrix_, macro_String] :=
    Module[{goldstoneList},
           goldstoneList = Cases[SARAH`GoldstoneGhost,
                                 {vector_, particle[{idx_}]} :> {idx, vector, particle, mixingMatrix}];
           StringJoin[CreateReorderingFunctionCalls[#,macro]& /@ goldstoneList]
          ];

ReorderGoldstoneBosons[particles_List, macro_String] :=
    Module[{result = ""},
           (result = result <> ReorderGoldstoneBosons[#,macro])& /@ particles;
           Return[result];
          ];

ReorderGoldstoneBosons[particle_Symbol, macro_String] :=
    ReorderGoldstoneBosons[particle, FindMixingMatrixSymbolFor[particle], macro];

ReorderGoldstoneBosons[particle_[___], macro_String] :=
    ReorderGoldstoneBosons[particle, macro];

ReorderGoldstoneBosons[macro_String] :=
    ReorderGoldstoneBosons[GetParticles[], macro];

CheckPoleMassesForTachyons[particles_List, macro_String] :=
    StringJoinWithSeparator[CheckPoleMassesForTachyons[#, macro]& /@ particles, "\n"];

CheckPoleMassesForTachyons[particle_, macro_String] :=
    Module[{dimStart, dimEnd, particleName},
           dimStart = GetDimensionStartSkippingGoldstones[particle];
           dimEnd   = GetDimension[particle];
           particleName = CConversion`ToValidCSymbolString[particle];
           (* skip check if all particles in the multiplet are tachyons *)
           If[(dimEnd == 1 && IsGoldstone[particle]) ||
              (dimStart >= dimEnd + 1),
              Return[""];
             ];
           "if (" <>
           WrapMacro[CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]],macro] <>
           If[dimEnd > 1, ".tail<" <> ToString[dimEnd - dimStart + 1] <> ">().minCoeff()", ""]<>
           " < 0.) { " <> FlagPoleTachyon[particleName] <> " }"
          ];

CheckPoleMassesForTachyons[macro_String] :=
    CheckPoleMassesForTachyons[Select[GetParticles[], (IsScalar[#] && !IsGoldstone[#])&], macro];

GetHiggsName[sym_] :=
    Switch[sym,
           SARAH`ChargedHiggs, "ChargedHiggs",
           SARAH`PseudoScalar, "PseudoscalarHiggs",
           SARAH`HiggsBoson  , "Higgs",
           _                 , ""
          ];

CallHiggsMassGetterFunction[name_String] :=
    "get_M" <> name <> "()";

CallPseudoscalarHiggsMassGetterFunction[] :=
    CallHiggsMassGetterFunction[GetHiggsName[SARAH`PseudoScalar]];

(* function that fills array with vector boson masses *)
FillGoldstoneMassVector[targetVector_String, vectorList_List] :=
    Module[{i, result = ""},
           For[i = 0, i < Length[vectorList], i++,
               result = result <> targetVector <> "(" <> ToString[i] <> ") = " <>
                        CConversion`ToValidCSymbolString[FlexibleSUSY`M[vectorList[[i+1]]]] <>
                        ";\n";
              ];
           result
          ];

CreateHiggsMassGetters[particle_[___], macro_String] :=
    CreateHiggsMassGetters[particle, macro];

CreateHiggsMassGetters[particle_, macro_String] :=
    Module[{vectorList, prototype, def, particleStr,
            particleHiggsStr, particleGoldstoneStr,
            typeHiggs, typeGoldstone, body, name,
            dim, dimGoldstone, dimHiggs},
           name                 = GetHiggsName[particle];
           particleStr          = CConversion`ToValidCSymbolString[FlexibleSUSY`M[particle]];
           particleHiggsStr     = particleStr <> "_" <> name;
           particleGoldstoneStr = particleStr <> "_goldstone";
           vectorList = Cases[SARAH`GoldstoneGhost,
                                 {vector_, particle[{_}]} :> vector];
           dim          = GetDimension[particle];
           dimGoldstone = Length[vectorList];
           (* number of physical (non-goldstone) particles *)
           dimHiggs     = dim - dimGoldstone;
           (* If dimHiggs == 0, all particles in the particle
              multiplet are Goldstone bosons and no one is a Higgs.

              If dimGoldstone == 0, all particles in the particle
              multiplet are Higgs bosons and there is no point in
              generating this function.
            *)
           If[dimHiggs <= 0 || dimGoldstone == 0
              || !MemberQ[GetParticles[], particle]
              ,
              If[dimHiggs < 0,
                 Print["Error: CreateHiggsMassGetters: There are more",
                       " Goldstone bosons than Higgs bosons."];
                 Print["   Dimension of ", particle, " = ", dim];
                 Print["   Number of Goldstones = ", dimGoldstone];
                 Return[{"",""}];
                ];
              Return[{"",""}];
             ];
           typeHiggs     = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dimHiggs]];
           typeGoldstone = CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dimGoldstone]];
           prototype = typeHiggs <> " get_M" <> name <> "() const;\n";
           body =
               typeGoldstone <> " " <> particleGoldstoneStr <> ";\n" <>
               FillGoldstoneMassVector[particleGoldstoneStr, vectorList] <>
               "\n" <>
               "return remove_if_equal(" <> particleStr <> ", " <>
                                       particleGoldstoneStr <> ");\n";
           def = typeHiggs <> " CLASSNAME::get_M" <> name <> "() const\n{\n" <>
               IndentText[body] <>
               "}\n";
           {prototype, def}
          ];

FlagPoleTachyon[particle_String, problems_String:"problems."] :=
    problems <> "flag_pole_tachyon(" <>
    FlexibleSUSY`FSModelName <> "_info::" <> particle <>
    ");";

FlagRunningTachyon[particle_String, problems_String:"problems."] :=
    problems <> "flag_running_tachyon(" <>
    FlexibleSUSY`FSModelName <> "_info::" <> particle <>
    ");";

FlagRunningTachyon[particles_List] :=
    StringJoinWithSeparator[FlagRunningTachyon /@ particles, "\n"];

FlagRunningTachyon[particle_] :=
    FlagRunningTachyon[CConversion`ToValidCSymbolString[CConversion`GetHead[particle]]];

CheckRunningTachyon[particle_, eigenvector_String] :=
    "if (" <> eigenvector <> If[GetDimension[particle] > 1, ".minCoeff()", ""] <> " < 0.) {\n" <>
    IndentText[FlagRunningTachyon[particle]] <>
    "\n}\n";

FlagBadMass[particle_String, eigenvalue_String] :=
    "problems.flag_bad_mass(" <> FlexibleSUSY`FSModelName <> "_info::" <> particle <>
    ", eigenvalue_error > precision * Abs(" <> eigenvalue <> "(0)));\n";

FlagBadMass[particles_List, eigenvalue_String] :=
    StringJoin[FlagBadMass[#,eigenvalue]& /@ particles];

FlagBadMass[particle_, eigenvalue_String] :=
    FlagBadMass[CConversion`ToValidCSymbolString[CConversion`GetHead[particle]], eigenvalue];

CallSVDFunction[particle_, matrix_String, eigenvalue_String, U_String, V_String] :=
    "\
#ifdef CHECK_EIGENVALUE_ERROR
" <> IndentText[
"double eigenvalue_error;
fs_svd(" <> matrix <> ", " <> eigenvalue <> ", " <> U <> ", " <> V <> ", eigenvalue_error);\n" <>
    If[ContainsMassless[particle],"",FlagBadMass[particle, eigenvalue]]
] <> "#else
" <> IndentText["\
fs_svd(" <> matrix <> ", " <> eigenvalue <> ", " <> U <> ", " <> V <> ");"] <> "
#endif
";

CallDiagonalizationFunction[particle_, matrix_String, eigenvalue_String, U_String, function_String] :=
"#ifdef CHECK_EIGENVALUE_ERROR
" <> IndentText[
"double eigenvalue_error;
" <> function <> "(" <> matrix <> ", " <> eigenvalue <> ", " <> U <> ", eigenvalue_error);\n" <>
    If[ContainsMassless[particle],"",FlagBadMass[particle, eigenvalue]]
] <> "#else
" <> IndentText["\n" <> function <> "(" <> matrix <> ", " <> eigenvalue <> ", " <> U <> ");\n"
] <> "#endif\n" <>
IndentText[
    If[IsVector[particle], U <> ".transposeInPlace();\n", ""] <>
"normalize_to_interval(" <> U <> ");"
] <> "\n";

CallDiagonalizeSymmetricFunction[particle_, matrix_String, eigenvector_String, U_String] :=
    CallDiagonalizationFunction[particle, matrix, eigenvector, U, "fs_diagonalize_symmetric"];

CallDiagonalizeHermitianFunction[particle_, matrix_String, eigenvector_String, U_String] :=
    CallDiagonalizationFunction[particle, matrix, eigenvector, U, "fs_diagonalize_hermitian"];

AssignVectorBosonMassesFrom[eigenVector_List] :=
    Module[{i, ev, result = ""},
           ev = CConversion`ToValidCSymbolString[FlexibleSUSY`M[CConversion`GetHead[MakeESSymbol[eigenVector]]]];
           For[i = 1, i <= Length[eigenVector], i++,
               result = result <>
                        CConversion`ToValidCSymbolString[FlexibleSUSY`M[eigenVector[[i]]]] <>
                        " = " <>
                        If[IsMassless[eigenVector[[i]]],
                           "0.;\n",
                           ev <> "(" <> ToString[i-1] <> ");\n"
                          ];
             ];
           result
          ];

AssignVectorBosonMassesFrom[eigenVector_] := "";

CreateDiagonalizationFunction[matrix_List, eigenVector_, mixingMatrixSymbol_] :=
    Module[{dim, body, result, U = "", V = "", dimStr = "", ev, evMap, particle, k,
            diagFunctionStr, matrixType, matrixElementType, vectorType,
            OneDimMappingPre = "", OneDimMappingPost = ""},
           dim = Length[matrix];
           dimStr = ToString[dim];
           particle = CConversion`ToValidCSymbolString[CConversion`GetHead[MakeESSymbol[eigenVector]]];
           matrixSymbol = "mass_matrix_" <> particle;
           ev = CConversion`ToValidCSymbolString[FlexibleSUSY`M[CConversion`GetHead[MakeESSymbol[eigenVector]]]];
           evMap = ev;
           result = "void CLASSNAME::calculate_" <> ev <> "()\n{\n";
           body = IndentText["const auto " <> matrixSymbol <> "(get_" <> matrixSymbol <> "());\n"];
           (* declare vector boson multiplet *)
           If[IsVector[eigenVector] && Head[eigenVector] === List,
              body = body <>
                     IndentText[
                         CreateCType[CConversion`ArrayType[CConversion`realScalarCType, dim]] <>
                         " " <> ev <> ";\n"
                     ];
             ];
           (* map scalar matrix and eigenvector to matrix and array types *)
           If[dim == 1,
              matrixElementType = CConversion`GetElementType[GetMassMatrixType[CConversion`GetHead[eigenvector]]];
              matrixType = CreateCType[CConversion`MatrixType[CConversion`realScalarCType,1,1]];
              vectorType = CreateCType[CConversion`ArrayType[CConversion`realScalarCType,1]];
              OneDimMappingPre = IndentText[
                  matrixType <> " " <> matrixSymbol <> "_map;\n" <>
                  matrixSymbol <> "_map(0,0) = " <> matrixSymbol <> ";\n" <>
                  vectorType <> " " <> ev <> "_map;"] <> "\n";
              OneDimMappingPost = IndentText[
                  ev <> " = " <> ev <> "_map(0);"] <> "\n";
              matrixSymbol = matrixSymbol <> "_map";
              evMap = ev <> "_map";
             ];
           If[Head[mixingMatrixSymbol] === List && Length[mixingMatrixSymbol] == 2,
              (* use SVD *)
              U = CConversion`ToValidCSymbolString[mixingMatrixSymbol[[1]]];
              V = CConversion`ToValidCSymbolString[mixingMatrixSymbol[[2]]];
              body = body <> "\n" <> OneDimMappingPre <> "\n" <>
                     CallSVDFunction[eigenVector, matrixSymbol, evMap, U, V] <> "\n" <>
                     OneDimMappingPost;
              ,
              (* use conventional diagonalization *)
              U = CConversion`ToValidCSymbolString[mixingMatrixSymbol];
              If[IsSymmetric[matrix] && Head[eigenVector] =!= List && IsFermion[CConversion`GetHead[eigenVector]],
                 body = body <> "\n" <> OneDimMappingPre <> "\n" <>
                        CallDiagonalizeSymmetricFunction[eigenVector, matrixSymbol, evMap, U] <> "\n" <>
                        OneDimMappingPost;,
                 body = body <> "\n" <> OneDimMappingPre <> "\n" <>
                        CallDiagonalizeHermitianFunction[eigenVector, matrixSymbol, evMap, U] <> "\n" <>
                        OneDimMappingPost;
                ];
             ];
           If[IsScalar[eigenVector] || IsVector[eigenVector],
              (* check for tachyons *)
              body = body <> "\n" <>
                     IndentText[
                         If[ContainsMassless[eigenVector], "",
                            CheckRunningTachyon[eigenVector, ev] <> "\n"] <>
                         ev <> " = AbsSqrt(" <> ev <> ");\n"
                     ];
             ];
           If[IsVector[eigenVector],
              (* assign individual vector boson masses to the
                 calculated mass eigenvalues *)
              body = body <> "\n" <>
                     IndentText[AssignVectorBosonMassesFrom[eigenVector]];
             ];
           Return[result <> body <> "}\n"];
          ];

CreateMassCalculationFunction[m:TreeMasses`FSMassMatrix[mass_, massESSymbol_, Null]] :=
    Module[{result, ev = CConversion`ToValidCSymbolString[FlexibleSUSY`M[massESSymbol]], body,
            inputParsDecl, expr, particle, dim, dimStr, phase, massMatrix,
            massMatrixStr},
           result = "void CLASSNAME::calculate_" <> ev <> "()\n{\n";
           (* Remove color SU(3) generators, structure functions and
              Kronecker delta with color indices.
              Note: ct1 ... ct4 are reserved SU(3) color indices of
              the fundamental representation of SU(3) in SARAH.
           *)
           expr = StripGenerators[mass[[1]],
                                  {SARAH`ct1, SARAH`ct2, SARAH`ct3, SARAH`ct4}];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           inputParsDecl = Parameters`CreateLocalConstRefsForInputParameters[expr, "LOCALINPUT"];
           particle = CConversion`ToValidCSymbolString[massESSymbol];
           massMatrixStr = "mass_matrix_" <> particle;
           If[dim == 1,
              body = inputParsDecl <> "\n" <>
                     "const auto " <> massMatrixStr <> " = " <>
                     "get_" <> massMatrixStr <> "();\n";
              (* adapt phases of massive fermions *)
              phase = Parameters`GetPhase[massESSymbol];
              Which[IsFermion[massESSymbol] && phase =!= Null && !IsMassless[massESSymbol],
                    body = body <> ev <> " = calculate_" <>
                           If[IsMajoranaFermion[massESSymbol], "majorana", "dirac"] <>
                           "_singlet_mass(" <> massMatrixStr <> ", " <>
                           CConversion`ToValidCSymbolString[phase] <> ");\n";
                    ,
                    IsFermion[massESSymbol],
                    body = body <> ev <> " = calculate_singlet_mass(" <> massMatrixStr <> ");\n";
                    ,
                    IsVector[massESSymbol] || IsScalar[massESSymbol],
                    body = body <> ev <> " = " <> massMatrixStr <> ";\n";
                    ,
                    True,
                    Print["Error: unknown particle type of ", massESSymbol];
                    Quit[1];
                   ];
              ,
              If[FreeQ[expr, SARAH`gt1] && FreeQ[expr, SARAH`gt2],
                 body = inputParsDecl <> "\n" <> ev <>
                        ".setConstant(" <> RValueToCFormString[expr] <> ");\n";
                 ,
                 body = inputParsDecl <> "\n" <>
                        "for (int gt1 = 1; gt1 <= " <> dimStr <> "; gt1++) {\n" <>
                        IndentText[ev <> "(gt1 - 1) = " <> RValueToCFormString[expr /. SARAH`gt2 -> SARAH`gt1] <> ";"] <>
                        "\n}\n";
                ];
             ];
           (* check for tachyons *)
           If[(IsVector[massESSymbol] || IsScalar[massESSymbol]) &&
              !IsMassless[massESSymbol],
              body = body <> "\n" <>
                     If[ContainsMassless[eigenVector], "",
                        CheckRunningTachyon[massESSymbol, ev] <> "\n"] <>
                     ev <> " = AbsSqrt(" <> ev <> ");\n";
             ];
           body = IndentText[body];
           result = result <> body <> "}\n\n";
           massMatrix = CreateFSMassMatrixForUnmixedParticle[m];
           result = CreateMassMatrixGetterFunction[massMatrix] <>
                    "\n" <> result;
           Return[result];
          ];

CreateMassCalculationFunction[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol, mixingMatrixSymbol, matrix},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           matrix = GetMassMatrix[massMatrix];
           result = result <>
                    CreateMassMatrixGetterFunction[massMatrix] <> "\n" <>
                    CreateDiagonalizationFunction[matrix, massESSymbol,
                                                  mixingMatrixSymbol]
                    <> "\n";
           Return[result];
          ];

CreatePhysicalMassDefinition[p:TreeMasses`FSMassMatrix[_, massESSymbols_List, _]] :=
    Module[{massMatrices},
           massMatrices = DeleteDuplicates[TreeMasses`FSMassMatrix[0, #, Null]& /@ massESSymbols];
           StringJoin[CreatePhysicalMassDefinition /@ massMatrices]
          ];

CreatePhysicalMassDefinition[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result = "", massESSymbol, dim, dimStr, type},
           massESSymbol = GetMassEigenstate[massMatrix];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`realScalarCType];,
              type = CConversion`ArrayType[CConversion`realScalarCType, dim];
             ];
           Parameters`CreateParameterDefinitionAndDefaultInitialize[
               {CConversion`ToValidCSymbolString[FlexibleSUSY`M[MakeESSymbol[massESSymbol]]], type}
           ]
          ];

DefineMatrix[Null, _] := "";

DefineMatrix[matrix_Symbol, type_] :=
    Parameters`CreateParameterDefinitionAndDefaultInitialize[{CConversion`ToValidCSymbolString[matrix], type}];

DefineMatrix[matrix_List, type_] :=
    StringJoin[DefineMatrix[#,type]& /@ matrix];

CreateMixingMatrixDefinition[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{mixingMatrixSymbol, type},
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           type = GetMixingMatrixType[massMatrix];
           DefineMatrix[mixingMatrixSymbol, type]
          ];

ClearOutputParameters[p:TreeMasses`FSMassMatrix[_,massESSymbols_List,_]] :=
    Module[{massMatrices},
           massMatrices = DeleteDuplicates[TreeMasses`FSMassMatrix[0, #, Null]& /@ massESSymbols];
           StringJoin[ClearOutputParameters /@ massMatrices]
          ];

ClearOutputParameters[massMatrix_TreeMasses`FSMassMatrix] :=
    Module[{result, massESSymbol, mixingMatrixSymbol, matrixType,
            dim, massESType},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = GetDimension[massESSymbol];
           massESType = Parameters`GetRealTypeFromDimension[{dim}];
           result = CConversion`SetToDefault[CConversion`ToValidCSymbolString[FlexibleSUSY`M[MakeESSymbol[massESSymbol]]],
                                             massESType];
           If[mixingMatrixSymbol =!= Null,
              matrixType = GetMixingMatrixType[massMatrix];
              If[Head[mixingMatrixSymbol] === List,
                 (result = result <>
                           CConversion`SetToDefault[CConversion`ToValidCSymbolString[#], matrixType])& /@ mixingMatrixSymbol;
                 ,
                 result = result <>
                          CConversion`SetToDefault[CConversion`ToValidCSymbolString[mixingMatrixSymbol], matrixType];
                ];
             ];
           Return[result];
          ];

CopyRunningMassesFromTo[p:TreeMasses`FSMassMatrix[_, massESSymbols_List, mix_], from_String, to_String] :=
    Module[{massMatrices},
           massMatrices = Join[
                   DeleteDuplicates[TreeMasses`FSMassMatrix[0, #, Null]& /@ massESSymbols],
                   { TreeMasses`FSMassMatrix[0, Null, mix] }
           ];
           StringJoin[CopyRunningMassesFromTo[#, from, to]& /@ massMatrices]
          ];

CopyRunningMassesFromTo[massMatrix_TreeMasses`FSMassMatrix, from_String, to_String] :=
    Module[{result = "", massESSymbol, mixingMatrixSymbol, dim, dimStr,
            i, massStr, mixStr},
           massESSymbol = GetMassEigenstate[massMatrix];
           mixingMatrixSymbol = GetMixingMatrixSymbol[massMatrix];
           dim = GetDimension[massESSymbol];
           dimStr = ToString[dim];
           (* copy mass *)
           If[massESSymbol =!= Null,
              massStr = CConversion`ToValidCSymbolString[FlexibleSUSY`M[MakeESSymbol[massESSymbol]]];
              result = WrapMacro[massStr, to] <> " = " <> WrapMacro[massStr, from] <> ";\n";
           ];
           (* copy mixings *)
           If[mixingMatrixSymbol =!= Null,
              If[Head[mixingMatrixSymbol] === List,
                 For[i = 1, i <= Length[mixingMatrixSymbol], i++,
                     mixStr = CConversion`ToValidCSymbolString[mixingMatrixSymbol[[i]]];
                     result = result <> WrapMacro[mixStr, to] <> " = " <> WrapMacro[mixStr, from] <> ";\n";
                    ];
                 ,
                 mixStr = CConversion`ToValidCSymbolString[mixingMatrixSymbol];
                 result = result <> WrapMacro[mixStr, to] <> " = " <> WrapMacro[mixStr, from] <> ";\n";
                ];
             ];
           result
          ];

CopyDRBarMassesToPoleMasses[p_] :=
        CopyRunningMassesFromTo[p, "", "PHYSICAL"];

FindColorGaugeGroup[] :=
    Module[{coupling, gaugeGroup, result},
           coupling = FindColorGaugeCoupling[];
           gaugeGroup = Cases[SARAH`Gauge, {_, group_, name_, coupling, ___}];
           If[gaugeGroup === {},
              Print["Error: could not find color gauge group"];
              result = Null;
              ,
              result = {coupling, gaugeGroup[[1,3]], gaugeGroup[[1,2]]};
             ];
           result
          ];

FindLeftGaugeGroup[] :=
    Module[{coupling, gaugeGroup, result},
           coupling = FindLeftGaugeCoupling[];
           gaugeGroup = Cases[SARAH`Gauge, {_, group_, name_, coupling, ___}];
           If[gaugeGroup === {},
              Print["Error: could not weak gauge group"];
              result = Null;
              ,
              result = {coupling, gaugeGroup[[1,3]], gaugeGroup[[1,2]]};
             ];
           result
          ];

FindHyperchargeGaugeGroup[] :=
    Module[{coupling, gaugeGroup, result},
           coupling = FindHyperchargeGaugeGroup[];
           gaugeGroup = Cases[SARAH`Gauge, {_, group_, name_, coupling, ___}];
           If[gaugeGroup === {},
              Print["Error: could not find Hypercharge gauge group"];
              result = Null;
              ,
              result = {coupling, gaugeGroup[[1,3]], gaugeGroup[[1,2]]};
             ];
           result
          ];

FindColorGaugeCoupling[] := SARAH`strongCoupling;

FindLeftGaugeCoupling[] := SARAH`leftCoupling;

FindHyperchargeGaugeCoupling[] := SARAH`hyperchargeCoupling;

GetSMVEVExpr[symbIfUndefined_:Undefined] :=
    Module[{vexp},
           If[ValueQ[SARAH`VEVSM],
              vexp = Cases[Parameters`GetAllDependenceSPhenoRules[],
                           RuleDelayed[SARAH`VEVSM, expr_] :> expr];
              If[vexp === {},
                 If[Parameters`IsParameter[SARAH`VEVSM],
                    SARAH`VEVSM,
                    DebugPrint["Warning: SM-like Higgs vev is not define in the SARAH model file!"];
                    symbIfUndefined]
                 ,
                 First[vexp]
              ]
              ,
              DebugPrint["Warning: SM-like Higgs vev is not define in the SARAH model file!"];
              symbIfUndefined
           ]
    ];

PrivateGetDependenceSPhenoRules[] :=
    DecreaseIndexLiterals @ Join[
        Parameters`GetDependenceSPhenoRules[],
        { FlexibleSUSY`VEV -> GetSMVEVExpr[0] }
    ];

CreateDependencePrototype[(Rule | RuleDelayed)[parameter_, _]] :=
    CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]] <>
    " " <> CConversion`ToValidCSymbolString[parameter] <> "() const;\n";

CreateDependencePrototypes[] :=
    StringJoin[CreateDependencePrototype /@ PrivateGetDependenceSPhenoRules[]]

CreateDependenceFunction[(Rule | RuleDelayed)[parameter_, value_]] :=
    Module[{result, body, parStr},
           parStr = CConversion`ToValidCSymbolString[parameter];
           body = Parameters`CreateLocalConstRefsForInputParameters[value, "LOCALINPUT"] <> "\n" <>
                  "return " <> RValueToCFormString[Simplify[value]] <> ";\n";
           result = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]] <>
                    " CLASSNAME::" <> parStr <> "() const\n{\n" <>
                    IndentText[body] <> "}\n\n";
           Return[result];
          ];

CreateDependenceFunctions[] :=
    StringJoin[CreateDependenceFunction /@ PrivateGetDependenceSPhenoRules[]]

CreateDependencyFunctionSymbols[] :=
    RuleDelayed[#,#[]]& /@ Parameters`GetDependenceSPhenoSymbols[];

ReplaceDependencies[expr_] :=
    expr /. CreateDependencyFunctionSymbols[];

ReplaceDependenciesReverse[expr_] :=
    expr /. (Reverse /@ CreateDependencyFunctionSymbols[]);

CallGenerationHelperFunctionName[gen_, fermion_, msf1_String, msf2_String, theta_] :=
    "calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_" <> NiceGenNum[gen] <> "_generation(" <> msf1 <> ", " <> msf2 <> ", " <> theta <> ")";

CreateGenerationHelperPrototype[gen_, fermion_] :=
    "void calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_" <> NiceGenNum[gen] <> "_generation(double&, double&, double&) const;\n";

CreateGenerationHelperFunction[gen_, fermion_ /; fermion === SARAH`TopSquark] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_" <> NiceGenNum[gen] <> "_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftSquark[gen-1,gen-1]] <> ");
   sf_data.mr2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftUp[gen-1,gen-1]] <> ");
   sf_data.yf  = Re(" <> CConversion`RValueToCFormString[SARAH`UpYukawa[gen-1,gen-1]] <> ");
   sf_data.vd  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ");
   sf_data.vu  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ");
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = Re(" <> CConversion`RValueToCFormString[SARAH`TrilinearUp[gen-1,gen-1]] <> ");
   sf_data.mu  = Re(" <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ");
   sf_data.T3  = sfermions::Isospin[sfermions::up];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::up];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::up];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG(\"diagonalize_sfermions_2x2: stop tachyon\");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}
";

CreateGenerationHelperFunction[gen_, fermion_ /; fermion === SARAH`BottomSquark] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_" <> NiceGenNum[gen] <> "_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftSquark[gen-1,gen-1]] <> ");
   sf_data.mr2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftDown[gen-1,gen-1]] <> ");
   sf_data.yf  = Re(" <> CConversion`RValueToCFormString[SARAH`DownYukawa[gen-1,gen-1]] <> ");
   sf_data.vd  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ");
   sf_data.vu  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ");
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = Re(" <> CConversion`RValueToCFormString[SARAH`TrilinearDown[gen-1,gen-1]] <> ");
   sf_data.mu  = Re(" <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ");
   sf_data.T3  = sfermions::Isospin[sfermions::down];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::down];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::down];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG(\"diagonalize_sfermions_2x2: sbottom tachyon\");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}
";

CreateGenerationHelperFunction[gen_, fermion_ /; fermion === SARAH`Selectron] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_" <> NiceGenNum[gen] <> "_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftLeftLepton[gen-1,gen-1]] <> ");
   sf_data.mr2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftRightLepton[gen-1,gen-1]] <> ");
   sf_data.yf  = Re(" <> CConversion`RValueToCFormString[SARAH`ElectronYukawa[gen-1,gen-1]] <> ");
   sf_data.vd  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ");
   sf_data.vu  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ");
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = Re(" <> CConversion`RValueToCFormString[SARAH`TrilinearLepton[gen-1,gen-1]] <> ");
   sf_data.mu  = Re(" <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ");
   sf_data.T3  = sfermions::Isospin[sfermions::electron];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::electron];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::electron];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG(\"diagonalize_sfermions_2x2: selecton tachyon\");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}
";

CreateGenerationHelperFunction[gen_, fermion_ /; fermion === SARAH`Sneutrino] :=
    "void CLASSNAME::calculate_" <>
    CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
    "_" <> NiceGenNum[gen] <> "_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(" <> CConversion`RValueToCFormString[SARAH`SoftLeftLepton[gen-1,gen-1]] <> ");
   sf_data.mr2 = 0.;
   sf_data.yf  = 0.;
   sf_data.vd  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> ");
   sf_data.vu  = Re(" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> ");
   sf_data.gY  = " <> CConversion`RValueToCFormString[SARAH`hyperchargeCoupling /.
                                                      Parameters`ApplyGUTNormalization[]] <> ";
   sf_data.g2  = " <> CConversion`RValueToCFormString[SARAH`leftCoupling] <> ";
   sf_data.Tyf = 0.;
   sf_data.mu  = Re(" <> CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]] <> ");
   sf_data.T3  = sfermions::Isospin[sfermions::neutrino];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::neutrino];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::neutrino];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG(\"diagonalize_sfermions_2x2: sneutrino tachyon\");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}
";

NiceGenNum[gen_] :=
    Switch[gen,
           1, "1st",
           2, "2nd",
           3, "3rd",
           _, ToString[gen] <> "th"
          ];

CreateGenerationHelperFunction[gen_, fermion_] :=
    Module[{},
           Print["Error: ", fermion, " does not seem to be a", gen, "th"
                 " generation SM fermion."];
           "void CLASSNAME::calculate_" <>
           CConversion`ToValidCSymbolString[FlexibleSUSY`M[fermion]] <>
           "_" <> NiceGenNum[gen] <> "_generation(double& msf1, double& msf2, double& theta) const {}"
          ];

CreateGenerationHelpers[gen_] :=
    Module[{prototypes, functions},
           functions = CreateGenerationHelperFunction[gen, SARAH`TopSquark] <> "\n" <>
                       CreateGenerationHelperFunction[gen, SARAH`BottomSquark] <> "\n" <>
                       CreateGenerationHelperFunction[gen, SARAH`Sneutrino] <> "\n" <>
                       CreateGenerationHelperFunction[gen, SARAH`Selectron];
           prototypes = CreateGenerationHelperPrototype[gen, SARAH`TopSquark] <>
                        CreateGenerationHelperPrototype[gen, SARAH`BottomSquark] <>
                        CreateGenerationHelperPrototype[gen, SARAH`Sneutrino] <>
                        CreateGenerationHelperPrototype[gen, SARAH`Selectron];
           {prototypes, functions}
          ];

GetThirdGenerationMass[fermion_, cConvention_:True, bracket_:False] :=
    If[GetDimension[fermion] == 1,
       If[bracket,
          FlexibleSUSY`M[fermion][],
          FlexibleSUSY`M[fermion]
         ]
       ,
       FlexibleSUSY`M[fermion][3 - If[cConvention === True, 1, 0]]
      ];

GetLightestMass[par_] :=
    Module[{dim, mass},
           dim = GetDimension[par];
           If[dim == 1,
              mass = FlexibleSUSY`M[par];,
              mass = FlexibleSUSY`M[par][0];
             ];
           Return[mass];
          ];

CountNumberOfMasses[masses_List] :=
    Plus @@ (CountNumberOfMasses /@ masses);

CountNumberOfMasses[mass_FSMassMatrix] :=
    CountNumberOfMasses[GetMassEigenstate[mass]];

CountNumberOfMasses[mass_] :=
    GetDimension[mass];

CreateMixingMatrixArrayGetterSetter[masses_List, offset_Integer, func_] :=
    Module[{display = "", paramCount = 0, name = "", mass,
            i, assignment = "", nAssignments = 0, CreateMixingAssignment},
           CreateMixingAssignment[mix_ /; mix === Null, _, _] := {"", 0};
           CreateMixingAssignment[{mix1_, mix2_}, paramCount_, type_] :=
               Module[{assignment, nAssignments, a, n},
                      {a, n} = CreateMixingAssignment[mix1,paramCount,type];
                      assignment = a;
                      nAssignments = n;
                      {a, n} = CreateMixingAssignment[mix2,paramCount + n,type];
                      assignment = assignment <> a;
                      nAssignments += n;
                      {assignment, nAssignments}
                     ];
           CreateMixingAssignment[mix_, paramCount_, type_] :=
               func[CConversion`ToValidCSymbolString[mix], paramCount, type];
           CreateMixingAssignment[mix_, paramCount_] :=
               CreateMixingAssignment[GetMixingMatrixSymbol[mix], paramCount, GetMixingMatrixType[mix]];
           For[i = 1, i <= Length[masses], i++,
               mix = masses[[i]];
               {assignment, nAssignments} = CreateMixingAssignment[mix, paramCount + offset];
               display = display <> assignment;
               paramCount += nAssignments;
              ];
           {display, paramCount}
          ];

CreateMassArrayGetter[masses_List] :=
    Module[{display = "", paramCount = 0, mass,
            type, i, v, assignment = "", nAssignments = 0, mes},
           For[i = 1, i <= Length[masses], i++,
               If[Head[GetMassEigenstate[masses[[i]]]] === List,
                  mes = DeleteDuplicates[GetMassEigenstate[masses[[i]]]];
                  For[v = 1, v <= Length[mes], v++,
                      mass = FlexibleSUSY`M[mes[[v]]];
                      type = GetMassType[mass];
                      {assignment, nAssignments} = Parameters`CreateDisplayAssignment[mass, paramCount, type];
                      display = display <> assignment;
                      paramCount += nAssignments;
                     ];
                  ,
                  mass = FlexibleSUSY`M[GetMassEigenstate[masses[[i]]]];
                  type = GetMassType[mass];
                  {assignment, nAssignments} = Parameters`CreateDisplayAssignment[mass, paramCount, type];
                  display = display <> assignment;
                  paramCount += nAssignments;
                 ];
              ];
           display = "Eigen::ArrayXd pars(" <> ToString[paramCount] <> ");\n\n" <>
                     display <> "\n" <>
                     "return pars;";
           Return[display];
          ];

CreateMixingArrayGetter[masses_List] :=
    Module[{display = "", paramCount, assignment = "", nAssignments = 0},
           paramCount = CountNumberOfMasses[masses];
           {assignment, nAssignments} = CreateMixingMatrixArrayGetterSetter[masses, paramCount, Parameters`CreateDisplayAssignment];
           display = display <> assignment;
           paramCount += nAssignments;
           display = "pars.conservativeResize(" <> ToString[paramCount] <> ");\n\n" <>
                     display;
           Return[display];
          ];

CreateMassArraySetter[masses_List, array_String] :=
    Module[{set = "", paramCount = 0, name = "", mass,
            type, i, v, assignment = "", nAssignments = 0, mes},
           For[i = 1, i <= Length[masses], i++,
               If[Head[GetMassEigenstate[masses[[i]]]] === List,
                  mes = DeleteDuplicates[GetMassEigenstate[masses[[i]]]];
                  For[v = 1, v <= Length[mes], v++,
                      mass = FlexibleSUSY`M[mes[[v]]];
                      type = GetMassType[mass];
                      name = CConversion`ToValidCSymbolString[mass];
                      {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
                      set = set <> assignment;
                      paramCount += nAssignments;
                     ];
                  ,
                  mass = FlexibleSUSY`M[GetMassEigenstate[masses[[i]]]];
                  type = GetMassType[mass];
                  name = CConversion`ToValidCSymbolString[mass];
                  {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type];
                  set = set <> assignment;
                  paramCount += nAssignments;
                 ];
              ];
           Return[set];
          ];

CreateMixingArraySetter[masses_List, array_String] :=
    Module[{set = "", paramCount, assignment = "", nAssignments = 0},
           paramCount = CountNumberOfMasses[masses];
           {assignment, nAssignments} = CreateMixingMatrixArrayGetterSetter[masses, paramCount, Parameters`CreateSetAssignment];
           set = set <> assignment;
           paramCount += nAssignments;
           Return[set];
          ];

(*
   1. Once Dominik wanted to have functions identifying SM particles.
      This might be non-trivial in some models.
      For now, we just have wrappers that return SM particles using SARAH symbols
   2. If a particle does not exist in the model, we don't stop.
      Instead, we return Null and let the calling code decide how to handle it.
*)

GetPhoton[] :=
   If[ValueQ[SARAH`Photon],
      SARAH`Photon
   ];

GetGluon[] :=
   If[ValueQ[SARAH`Gluon],
      SARAH`Gluon
   ];

GetZBoson[] :=
   If[ValueQ[SARAH`Zboson],
      SARAH`Zboson
   ];

GetWBoson[] :=
   Module[{temp = Select[Unevaluated[{SARAH`Wboson, SARAH`VectorW}], ValueQ]},
      If[Length @ DeleteDuplicates[temp] === 1,
         temp[[1]]
      ]
   ];

GetHiggsBoson[] :=
   If[ValueQ[SARAH`HiggsBoson],
      SARAH`HiggsBoson
   ];

GetChargedHiggsBoson[] :=
   If[ValueQ[SARAH`ChargedHiggs],
      SARAH`ChargedHiggs
   ];

GetPseudoscalarHiggsBoson[] :=
   If[ValueQ[SARAH`PseudoScalarBoson], SARAH`PseudoScalarBoson];

GetStrongCoupling[] :=
   If[ValueQ[SARAH`strongCoupling], SARAH`strongCoupling];

End[];

EndPackage[];
