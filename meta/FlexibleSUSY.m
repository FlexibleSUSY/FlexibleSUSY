(* ::Package:: *)

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

BeginPackage["FlexibleSUSY`",
             {"SARAH`",
              "AnomalousDimension`",
              "BetaFunction`",
              "Parameters`",
              "TextFormatting`",
              "CConversion`",
              "TreeMasses`",
              "EWSB`",
              "Traces`",
              "SelfEnergies`",
              "Vertices`",
              "Phases`",
              "LatticeUtils`",
              "LoopMasses`",
              "WriteOut`",
              "Constraint`",
              "ThresholdCorrections`",
              "ConvergenceTester`",
              "FunctionModifiers`",
              "Utils`",
              "References`",
              "SemiAnalytic`",
              "ThreeLoopSM`",
              "TerminalFormatting`",
              "ThreeLoopMSSM`",
              "CXXDiagrams`",
              "AMM`",
              "Decays`",
              "EDM`",
              "FFVFormFactors`",
              "BtoSGamma`",
              "FlexibleEFTHiggsMatching`",
              "FSMathLink`",
              "FlexibleTower`",
              "WeinbergAngle`",
              "Wrappers`",
              "Himalaya`",
              "GM2Calc`",
              "Unitarity`"
}];

$flexiblesusyMetaDir     = DirectoryName[FindFile[$Input]];
$flexiblesusyConfigDir   = FileNameJoin[{ParentDirectory[$flexiblesusyMetaDir], "config"}];
$flexiblesusyTemplateDir = FileNameJoin[{ParentDirectory[$flexiblesusyMetaDir], "templates"}];
$observablesWildcard[file:_:""] := FileNameJoin@{$flexiblesusyMetaDir, "Observables", "*", file};

FS`Version = StringTrim[FSImportString[FileNameJoin[{$flexiblesusyConfigDir,"flexiblesusy-version"}]]];
FS`GitCommit = StringTrim[FSImportString[FileNameJoin[{$flexiblesusyConfigDir,"git_commit"}]]];
FS`Authors = {"P. Athron", "M. Bach", "D. Harries", "U. Khasianevich",
              "W. Kotlarski", "T. Kwasnitza", "J.-h. Park", "T. Steudtner",
              "D. St\[ODoubleDot]ckinger", "A. Voigt", "J. Ziebell"};
FS`Contributors = {};
FS`Years   = "2013-2025";
FS`References = Get[FileNameJoin[{$flexiblesusyConfigDir,"references"}]];

Print[""];
Utils`FSFancyLine["="];
Utils`FSFancyPrint["FlexibleSUSY " <> FS`Version, 0];
Print["  by " <> StringJoin[Riffle[Riffle[FS`Authors, ", "], "\n  ", 11]] <>
      "\n  " <> FS`Years];
If[FS`Contributors =!= {},
   Print["  contributions by " <> Utils`StringJoinWithSeparator[FS`Contributors, ", "]];
  ];
Print[""];
Utils`FSFancyPrint["References:"];
Print["  " <> #]& /@ FS`References;
Print[""];
Utils`FSFancyPrint["Download and Documentation:"];
Print["  https://flexiblesusy.hepforge.org"];
Utils`FSFancyLine["="];
Print[""];

LoadModelFile::usage="";
ReadSARAHBetaFunctions::usage="";
SetupModelParameters::usage="";
SetupMassMatrices::usage="";
SetupOutputParameters::usage="";
WriteClass::usage="Writes a C++ class for an observable";

MakeFlexibleSUSY::usage="Creates a spectrum generator given a
 FlexibleSUSY model file (FlexibleSUSY.m).

Example:

  MakeFlexibleSUSY[
      InputFile -> \"models/<model>/FlexibleSUSY.m\",
      OutputDirectory -> \"models/<model>/\",
      DebugOutput -> False];

Options:

  InputFile: The name of the model file.

  OutputDirectory: The output directory for the generated code.

  DebugOutput (True|False): Enable/Disable debug output while running
    the Mathematica meta code.
";

LowPrecision::usage="";
MediumPrecision::usage="";
HighPrecision::usage="";
GUTNormalization::usage="Returns GUT normalization of a given coupling";

BETA::usage = "Head for beta functions";
FSModelName;
FSOutputDir = ""; (* directory for generated code *)
FSLesHouchesList;
FSUnfixedParameters;
EWSBOutputParameters = {};
EWSBInitialGuess = {};
EWSBSubstitutions = {};
SUSYScale;
SUSYScaleFirstGuess;
SUSYScaleInput = {};
SUSYScaleMinimum;
SUSYScaleMaximum;
HighScale;
HighScaleFirstGuess;
HighScaleInput = {};
HighScaleMinimum;
HighScaleMaximum;
LowScale;
LowScaleFirstGuess;
LowScaleInput = {};
LowScaleMinimum;
LowScaleMaximum;
InitialGuessAtLowScale = {};
InitialGuessAtSUSYScale = {};
InitialGuessAtHighScale = {};
OnlyLowEnergyFlexibleSUSY = False;
FlexibleEFTHiggs = False;
MatchingScaleInput={};
AutomaticInputAtMSUSY = True; (* input unfixed parameters at MSUSY *)
TreeLevelEWSBSolution = {};
Pole;
LowEnergyConstant;
LowEnergyGaugeCoupling;
FSMinimize;
FSFindRoot;
FSFlagProblem;       (* flag a problem *)
FSFlagWarning;       (* flag a problem *)
(* potential problems *)
{ FSNoProblem, FSNonPerturbativeParameter, FSInvalidInputParameter };
FSRestrictParameter; (* restrict parameter to interval *)
FSInitialSetting;    (* set parameter before calculating masses *)
FSSolveEWSBFor;
FSSolveEWSBTreeLevelFor = {};
FSTemporary;
MZ;
MT;
MZDRbar;
MWDRbar;
MZMSbar;
MWMSbar;
EDRbar;
EMSbar;
ThetaWDRbar;
SCALE;
THRESHOLD;
VEV::usage = "running SM-like VEV in the full model";
UseHiggs2LoopNMSSM = False;
UseHiggs3LoopMSSM = False;
UseHiggs3LoopNMSSM = False;
EffectiveMu;
EffectiveMASqr;
UseSM3LoopRGEs = False;
UseSM4LoopRGEs = False;
UseSM5LoopRGEs = False;
UseSMAlphaS3Loop = False;
UseSMAlphaS4Loop = False;
UseSMYukawa2Loop = False;
UseMSSM3LoopRGEs = False;
UseMSSMYukawa2Loop = False;
UseMSSMAlphaS2Loop = False;
UseHiggs2LoopSM = False;
UseHiggs3LoopSM = False;
UseHiggs4LoopSM = False;
UseHiggs3LoopSplit = False;
UseYukawa3LoopQCD = Automatic;
UseYukawa4LoopQCD = Automatic;
FSRGELoopOrder = 2; (* RGE loop order (0, 1 or 2) *)
PotentialLSPParticles = {};
ExtraSLHAOutputBlocks;
FWCOEF ~ SetAttributes ~ {Protected, Locked}; (*reserved for a Block name*)
IMFWCOEF ~ SetAttributes ~ {Protected, Locked}; (*reserved for a Block name*)

FSAuxiliaryParameterInfo = {};
IMMINPAR = {};
IMEXTPAR = {};
FSCalculateDecays = False;
FSDecayParticles = Automatic;
FSEnableParallelism = True;
FSUnitarityConstraints = True;
FSEnableCompile;
FSGaugeLess::usage = "Symbol to represent a small number to impose the gauge-less limit.";
FSGaugeLessLimit = {
   {SARAH`hyperchargeCoupling, FSGaugeLess/Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]},
   {SARAH`leftCoupling, FSGaugeLess/Parameters`GetGUTNormalization[SARAH`leftCoupling]}
};
FSGaugeLessLimit::usage = "List of 2-component lists {parameter, value} to set parameter to obtain the gauge-less limit.";
FSYukawaLessLimit = {
    {SARAH`UpYukawa[1, 1], 0},
    {SARAH`UpYukawa[1, 2], 0},
    {SARAH`UpYukawa[1, 3], 0},
    {SARAH`UpYukawa[2, 1], 0},
    {SARAH`UpYukawa[2, 2], 0},
    {SARAH`UpYukawa[2, 3], 0},
    {SARAH`UpYukawa[3, 1], 0},
    {SARAH`UpYukawa[3, 2], 0},
    (* *)
    {SARAH`DownYukawa[1, 1], 0},
    {SARAH`DownYukawa[1, 2], 0},
    {SARAH`DownYukawa[1, 3], 0},
    {SARAH`DownYukawa[2, 1], 0},
    {SARAH`DownYukawa[2, 2], 0},
    {SARAH`DownYukawa[2, 3], 0},
    {SARAH`DownYukawa[3, 1], 0},
    {SARAH`DownYukawa[3, 2], 0},
    (* *)
    {SARAH`ElectronYukawa[1, 1], 0},
    {SARAH`ElectronYukawa[1, 2], 0},
    {SARAH`ElectronYukawa[1, 3], 0},
    {SARAH`ElectronYukawa[2, 1], 0},
    {SARAH`ElectronYukawa[2, 2], 0},
    {SARAH`ElectronYukawa[2, 3], 0},
    {SARAH`ElectronYukawa[3, 1], 0},
    {SARAH`ElectronYukawa[3, 2], 0}
};
FSYukawaLessLimit::usage = "List of 2-component lists {parameter, value} to set parameter to obtain the yukawa-less limit (leaving only the 3rd generation Yukawa couplings non-zero).";
FSMSSMLimit = {};
FSMSSMLimit::usage = "List of 2-component lists {parameter, value} to set parameter to obtain the MSSM-limit.";

(* Standard Model input parameters (SLHA input parameters) *)
(* {parameter, {"block", entry}, type}                     *)
SMINPUTS = {
    {AlphaEMInvInput    , {"SMINPUTS",  1}, CConversion`ScalarType[CConversion`realScalarCType]},
    {GFermiInput        , {"SMINPUTS",  2}, CConversion`ScalarType[CConversion`realScalarCType]},
    {AlphaSInput        , {"SMINPUTS",  3}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MZPoleInput        , {"SMINPUTS",  4}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MBottomMbottomInput, {"SMINPUTS",  5}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MTopPoleInput      , {"SMINPUTS",  6}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MTauPoleInput      , {"SMINPUTS",  7}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MNeutrino3PoleInput, {"SMINPUTS",  8}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MWPoleInput        , {"SMINPUTS", 10}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MElectronPoleInput , {"SMINPUTS", 11}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MNeutrino1PoleInput, {"SMINPUTS", 12}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MMuonPoleInput     , {"SMINPUTS", 13}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MNeutrino2PoleInput, {"SMINPUTS", 14}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MDown2GeVInput     , {"SMINPUTS", 21}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MUp2GeVInput       , {"SMINPUTS", 22}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MStrange2GeVInput  , {"SMINPUTS", 23}, CConversion`ScalarType[CConversion`realScalarCType]},
    {MCharmMCharm       , {"SMINPUTS", 24}, CConversion`ScalarType[CConversion`realScalarCType]}
};

(* renormalization schemes *)
DRbar;
MSbar;
FSRenormalizationScheme = DRbar;

(* all model parameters are real by default *)
SARAH`RealParameters = { All };

(* precision of pole mass calculation *)
DefaultPoleMassPrecision = MediumPrecision;
HighPoleMassPrecision    = {};
MediumPoleMassPrecision  = {};
LowPoleMassPrecision     = {};

FSEigenstates = SARAH`EWSB;
FSSolveEWSBTimeConstraint = 120;
FSSimplifyBetaFunctionsTimeConstraint = 120;
FSSolveWeinbergAngleTimeConstraint = 120;
FSCheckPerturbativityOfDimensionlessParameters = True;
FSPerturbativityThreshold = N[Sqrt[4 Pi]];
FSMaximumExpressionSize = 100;

(* list of masses and parameters to check for convergence

   Example:

   FSConvergenceCheck = {
      M[hh], g3, Yu, Yd[3,3], Ye, B[\[Mu]]
   };
*)
FSConvergenceCheck = Automatic;

(* list of parameters used by the semi-analytic solver to
   check for convergence of the inner two-scale iteration *)
SemiAnalyticSolverInnerConvergenceCheck = Automatic;

(* EWSB solvers *)
GSLHybrid;   (* hybrid method *)
GSLHybridS;  (* hybrid method with dynamic step size *)
GSLBroyden;  (* Broyden method *)
GSLNewton;   (* Newton method *)
FPIRelative; (* Fixed point iteration, convergence crit. relative step size *)
FPIAbsolute; (* Fixed point iteration, convergence crit. absolute step size *)
FPITadpole;  (* Fixed point iteration, convergence crit. relative step size + tadpoles *)
FSEWSBSolvers = { FPIRelative, GSLHybridS, GSLBroyden };

(* BVP solvers *)
TwoScaleSolver;      (* two-scale algorithm *)
LatticeSolver;       (* lattice algorithm *)
SemiAnalyticSolver;  (* semi-analytic algorithm *)
ShootingSolver;      (* shooting method, can be used for FlexibleEFTHiggs *)
FSBVPSolvers = { TwoScaleSolver };

(* macros *)
IF;
SUM;

(* methods for the calculation of the weak mixing angle *)
{ FSFermiConstant, FSMassW };

(* rules to apply on expressions *)
FSSelfEnergyRules = {{}, {}}; (* 1L, 2L *)
FSVertexRules = {};
FSBetaFunctionRules = {{}, {}, {}}; (* 1L, 2L, 3L *)

FSWeakMixingAngleInput = Automatic;
FSWeakMixingAngleInput::usage = "Method to determine the weak mixing
 angle. Possible values = { Automatic, FSFermiConstant, FSMassW }.
 Default value: Automatic";

FSWeakMixingAngleExpr = ArcSin[Sqrt[1 - Mass[SARAH`VectorW]^2/Mass[SARAH`VectorZ]^2]];
FSWeakMixingAngleExpr::usage = "Expression to fix weak mixing angle.
 Default value: ArcSin[Sqrt[1 -
    Mass[SARAH`VectorW]^2/Mass[SARAH`VectorZ]^2]]";

ReadPoleMassPrecisions::ImpreciseHiggs="Warning: Calculating the Higgs pole mass M[`1`] with `2` will lead to an inaccurate result!  Please select MediumPrecision or HighPrecision (recommended) for `1`.";

tadpole::usage="symbolic expression for a tadpole contribution in the
EWSB eqs.  The index corresponds to the ordering of the tadpole
equations in SARAH`TadpoleEquations[] .";

NoScale::usage="placeholder indicating an SLHA block should not
have a scale associated with it.";
CurrentScale::usage="placeholder indicating the current renormalization
scale of the model.";

(* input parameters for GM2Calc, see [arxiv:2110.13238] *)
FSGM2CalcInput = {
    yukawaType -> 0,
    lambda1 -> 0,
    lambda2 -> 0,
    lambda3 -> 0,
    lambda4 -> 0,
    lambda5 -> 0,
    lambda6 -> 0,
    lambda7 -> 0,
    tanBeta -> SARAH`VEVSM2 / SARAH`VEVSM1,
    m122 -> 0,
    zetau -> 0,
    zetad -> 0,
    zetal -> 0,
    deltau -> ZEROMATRIX[3,3],
    deltad -> ZEROMATRIX[3,3],
    deltal -> ZEROMATRIX[3,3],
    piu -> ZEROMATRIX[3,3],
    pid -> ZEROMATRIX[3,3],
    pil -> ZEROMATRIX[3,3]
};

(* input parameters for Himalaya *)
FSHimalayaInput = {
    RenormalizationScheme -> DRbar,
    Lambda3L -> 1,
    Lambda3LUncertainty -> 0,
    \[Mu] -> \[Mu],
    SARAH`g1 -> SARAH`hyperchargeCoupling,
    Susyno`LieGroups`g2 -> SARAH`leftCoupling,
    g3 -> SARAH`strongCoupling,
    vu -> SARAH`VEVSM2,
    vd -> SARAH`VEVSM1,
    MSQ2 -> SARAH`UpMatrixL,
    MSD2 -> SARAH`DownMatrixR,
    MSU2 -> SARAH`UpMatrixR,
    MSL2 -> SARAH`ElectronMatrixL,
    MSE2 -> SARAH`ElectronMatrixR,
    Au -> SARAH`TrilinearUp,
    Ad -> SARAH`TrilinearDown,
    Ae -> SARAH`TrilinearLepton,
    Yu -> SARAH`UpYukawa,
    Yd -> SARAH`DownYukawa,
    Ye -> SARAH`ElectronYukawa,
    M1 -> 0,
    M2 -> 0,
    M3 -> FlexibleSUSY`M[SARAH`Gluino],
    mA -> FlexibleSUSY`M[SARAH`PseudoScalar]
};

FSDebugOutput = False;

FSLoopLibraries::usage = "Contains a List of enabled loop libraries.";
FSLoopTools;
FSFFlite;
FSCOLLIER;
FSLoopLibraries = { FSSOFTSUSY };

FSFeynArtsAvailable = False;
FSFormCalcAvailable = False;
FSEnableColors = False;

Begin["`Private`"];

allIndexReplacementRules = {};

GetIndexReplacementRules[] := allIndexReplacementRules;

allBetaFunctions = {};

GetBetaFunctions[] := allBetaFunctions;

numberOfModelParameters = 0;

allEWSBSolvers = { GSLHybrid, GSLHybridS, GSLBroyden, GSLNewton,
                   FPIRelative, FPIAbsolute, FPITadpole };

allBVPSolvers = { TwoScaleSolver, LatticeSolver, SemiAnalyticSolver, ShootingSolver };

HaveEWSBSolver[solver_] := MemberQ[FlexibleSUSY`FSEWSBSolvers, solver];

HaveBVPSolver[solver_] := MemberQ[FlexibleSUSY`FSBVPSolvers, solver];

ToVersionString[{major_Integer, minor_Integer, patch_Integer}] :=
    ToString[major] <> "." <> ToString[minor] <> "." <> ToString[patch];

DebugPrint[msg___] :=
    If[FlexibleSUSY`FSDebugOutput,
       Print["Debug<FlexibleSUSY>: ", Sequence @@ InputFormOfNonStrings /@ {msg}]];

ReplaceSymbolsInUserInput[rules_] :=
    Module[{},
           (* input/output blocks *)
           SARAH`MINPAR                          = SARAH`MINPAR                          /. rules;
           SARAH`EXTPAR                          = SARAH`EXTPAR                          /. rules;
           IMMINPAR                              = IMMINPAR                              /. rules;
           IMEXTPAR                              = IMEXTPAR                              /. rules;
           FlexibleSUSY`ExtraSLHAOutputBlocks    = FlexibleSUSY`ExtraSLHAOutputBlocks    /. rules;
           FlexibleSUSY`FSCalculateDecays        = FlexibleSUSY`FSCalculateDecays        /. rules;
           FlexibleSUSY`FSUnitarityConstraints   = FlexibleSUSY`FSUnitarityConstraints   /. rules;
           FlexibleSUSY`FSDecayParticles         = FlexibleSUSY`FSDecayParticles         /. rules;
           FlexibleSUSY`FSExtraInputParameters   = FlexibleSUSY`FSExtraInputParameters   /. rules;
           FlexibleSUSY`FSAuxiliaryParameterInfo = FlexibleSUSY`FSAuxiliaryParameterInfo /. rules;
           FlexibleSUSY`FSLesHouchesList         = FlexibleSUSY`FSLesHouchesList         /. rules;
           FlexibleSUSY`FSHimalayaInput          = FlexibleSUSY`FSHimalayaInput          /. rules;
           (* boundary conditions *)
           FlexibleSUSY`InitialGuessAtLowScale   = FlexibleSUSY`InitialGuessAtLowScale   /. rules;
           FlexibleSUSY`InitialGuessAtSUSYScale  = FlexibleSUSY`InitialGuessAtSUSYScale  /. rules;
           FlexibleSUSY`InitialGuessAtHighScale  = FlexibleSUSY`InitialGuessAtHighScale  /. rules;
           FlexibleSUSY`LowScale                 = FlexibleSUSY`LowScale                 /. rules;
           FlexibleSUSY`LowScaleFirstGuess       = FlexibleSUSY`LowScaleFirstGuess       /. rules;
           FlexibleSUSY`LowScaleInput            = FlexibleSUSY`LowScaleInput            /. rules;
           FlexibleSUSY`LowScaleMinimum          = FlexibleSUSY`LowScaleMinimum          /. rules;
           FlexibleSUSY`LowScaleMaximum          = FlexibleSUSY`LowScaleMaximum          /. rules;
           FlexibleSUSY`SUSYScale                = FlexibleSUSY`SUSYScale                /. rules;
           FlexibleSUSY`SUSYScaleFirstGuess      = FlexibleSUSY`SUSYScaleFirstGuess      /. rules;
           FlexibleSUSY`SUSYScaleInput           = FlexibleSUSY`SUSYScaleInput           /. rules;
           FlexibleSUSY`SUSYScaleMinimum         = FlexibleSUSY`SUSYScaleMinimum         /. rules;
           FlexibleSUSY`SUSYScaleMaximum         = FlexibleSUSY`SUSYScaleMaximum         /. rules;
           FlexibleSUSY`HighScale                = FlexibleSUSY`HighScale                /. rules;
           FlexibleSUSY`HighScaleFirstGuess      = FlexibleSUSY`HighScaleFirstGuess      /. rules;
           FlexibleSUSY`HighScaleInput           = FlexibleSUSY`HighScaleInput           /. rules;
           FlexibleSUSY`HighScaleMinimum         = FlexibleSUSY`HighScaleMinimum         /. rules;
           FlexibleSUSY`HighScaleMaximum         = FlexibleSUSY`HighScaleMaximum         /. rules;
           FlexibleSUSY`SUSYScaleMatching        = FlexibleSUSY`SUSYScaleMatching        /. rules;
           FlexibleSUSY`MatchingScaleInput       = FlexibleSUSY`MatchingScaleInput       /. rules;
           (* EWSB *)
           FlexibleSUSY`EWSBInitialGuess         = FlexibleSUSY`EWSBInitialGuess         /. rules;
           FlexibleSUSY`EWSBOutputParameters     = FlexibleSUSY`EWSBOutputParameters     /. rules;
           FlexibleSUSY`EWSBSubstitutions        = FlexibleSUSY`EWSBSubstitutions        /. rules;
           FlexibleSUSY`TreeLevelEWSBSolution    = FlexibleSUSY`TreeLevelEWSBSolution    /. rules;
           (* mass calculation *)
           FlexibleSUSY`LowPrecision             = FlexibleSUSY`LowPrecision             /. rules;
           FlexibleSUSY`MediumPrecision          = FlexibleSUSY`MediumPrecision          /. rules;
           FlexibleSUSY`HighPrecision            = FlexibleSUSY`HighPrecision            /. rules;
           FlexibleSUSY`DefaultPoleMassPrecision = FlexibleSUSY`DefaultPoleMassPrecision /. rules;
           FlexibleSUSY`HighPoleMassPrecision    = FlexibleSUSY`HighPoleMassPrecision    /. rules;
           FlexibleSUSY`MediumPoleMassPrecision  = FlexibleSUSY`MediumPoleMassPrecision  /. rules;
           FlexibleSUSY`LowPoleMassPrecision     = FlexibleSUSY`LowPoleMassPrecision     /. rules;
           (* RG running *)
           FlexibleSUSY`FSPerturbativityThreshold = FlexibleSUSY`FSPerturbativityThreshold /. rules;
           FlexibleSUSY`FSConvergenceCheck        = FlexibleSUSY`FSConvergenceCheck        /. rules;
           FlexibleSUSY`SemiAnalyticSolverInnerConvergenceCheck = FlexibleSUSY`SemiAnalyticSolverInnerConvergenceCheck /. rules;
           (* model-specific corrections *)
           FlexibleSUSY`EffectiveMu             = FlexibleSUSY`EffectiveMu               /. rules;
           FlexibleSUSY`EffectiveMASqr          = FlexibleSUSY`EffectiveMASqr            /. rules;
           (* simplifications *)
           FlexibleSUSY`FSSelfEnergyRules       = FlexibleSUSY`FSSelfEnergyRules         /. rules;
           FlexibleSUSY`FSVertexRules           = FlexibleSUSY`FSVertexRules             /. rules;
           FlexibleSUSY`FSBetaFunctionRules     = FlexibleSUSY`FSBetaFunctionRules       /. rules;

           SetAttributes[CheckParticleInPrecision, HoldFirst];
           CheckParticleInPrecision[precision_] :=
              If[!SubsetQ[TreeMasses`GetParticles[], precision],
                 Utils`FSFancyWarning["Particle(s) ",
                    Complement[precision, TreeMasses`GetParticles[]],
                    " cannot be used in ", Unevaluated@precision,
                    " as it is not part of ", ToString@FSModelName,
                    ". Removing it/them."
                 ];
                 precision = Intersection[precision, TreeMasses`GetParticles[]]
              ];
           CheckParticleInPrecision /@
              {Unevaluated@FlexibleSUSY`HighPoleMassPrecision, Unevaluated@FlexibleSUSY`MediumPoleMassPrecision, Unevaluated@FlexibleSUSY`LowPoleMassPrecision};

           (* decay calculation require 3- and 4-point loop functions *)
           If[FlexibleSUSY`FSCalculateDecays && DisjointQ[FSLoopLibraries, {FSLoopTools, FSCOLLIER}] && FSEnableCompile,
              Utils`FSFancyWarning[
                 "Decay calculation requires a dedicated loop library.",
                 " Currently it's either LoopTools or Collier but",
                 " FlexibleSUSY was only configured with internal libraries.",
                 " Disabling decays."
              ];
              FlexibleSUSY`FSCalculateDecays = False
           ];

           With[{sarahVersion = DecomposeVersionString[SA`Version],
                 minimRequired = {4, 13, 0}},
              If[!TrueQ@Utils`VersionOrderGtEqThan[sarahVersion, minimRequired],
                 Print["Warning: using SARAH version ", SA`Version, " but unitarity calculation requires ", StringRiffle[ToString/@minimRequired, "."]];
                 Print["         Disabling unitarity calculation."];
                 FlexibleSUSY`FSUnitarityConstraints = False;
              ];
           ];
          ];

CheckSARAHVersion[] :=
    Module[{minimRequired, minimRequiredVersionFile, sarahVersion},
           Print["Checking SARAH version ..."];
           minimRequiredVersionFile = FileNameJoin[{$flexiblesusyConfigDir,
                                                    "required_sarah_version.m"}];
           (* reading minimum required SARAH version from config file *)
           minimRequired = Get[minimRequiredVersionFile];
           If[minimRequired === $Failed,
              Print["Error: Cannot read required SARAH version from file ",
                    minimRequiredVersionFile];
              Print["   Did you run configure?"];
              Quit[1];
             ];
           sarahVersion = Utils`DecomposeVersionString[SA`Version];
           If[!TrueQ@Utils`VersionOrderGtEqThan[sarahVersion, minimRequired],
              Print["Error: SARAH version ", SA`Version, " no longer supported!"];
              Print["Please use version ", ToVersionString[minimRequired],
                    " or higher"];
              Quit[1];
             ];
          ];

CheckEWSBSolvers[solvers_List] :=
    Module[{invalidSolvers},
           invalidSolvers = Complement[solvers, allEWSBSolvers];
           If[invalidSolvers =!= {},
              Print["Error: invalid EWSB solvers requested: ", invalidSolvers];
              Quit[1];
             ];
          ];

CheckBVPSolvers[solvers_List] :=
    Module[{invalidSolvers},
           invalidSolvers = Complement[solvers, allBVPSolvers];
           If[invalidSolvers =!= {},
              Print["Error: invalid BVP solvers requested: ", invalidSolvers];
              Quit[1];
             ];
          ];

CheckDecaysOptions[] :=
   If[FlexibleSUSY`FSDecayParticles =!= Automatic && FlexibleSUSY`FSDecayParticles =!= All && !ListQ[FlexibleSUSY`FSDecayParticles],
      Utils`FSFancyWarning[
         "Allowed values for FSDecayParticles are Automatic, All or",
         " list of particles. Got ", FlexibleSUSY`FSDecayParticles,
         ". Disabling decays."
      ];
      FlexibleSUSY`FSCalculateDecays = False,
      If[FlexibleSUSY`FSDecayParticles === Automatic,
         FlexibleSUSY`FSDecayParticles =
            DeleteCases[{TreeMasses`GetHiggsBoson[], TreeMasses`GetChargedHiggsBoson[], TreeMasses`GetPseudoscalarHiggsBoson[]}, Null],
         If[FlexibleSUSY`FSDecayParticles === All,
            FlexibleSUSY`FSDecayParticles = TreeMasses`GetParticles[],
            FlexibleSUSY`FSDecayParticles = FlexibleSUSY`FSDecayParticles /. SARAH`bar|Susyno`LieGroups`conj -> Identity;
            If[!SubsetQ[TreeMasses`GetParticles[], FlexibleSUSY`FSDecayParticles],
               Utils`FSFancyWarning[
                  "Requested decay of particles ",
                  Complement[FlexibleSUSY`FSDecayParticles, TreeMasses`GetParticles[]],
                  " which are not part of the model. Removing them."];
               FlexibleSUSY`FSDecayParticles = Intersection[TreeMasses`GetParticles[], FlexibleSUSY`FSDecayParticles]
            ]
         ]
      ];
   ];

(* sets model file variables to default values, after SARAH`Start[] has been called *)
FSSetDefaultModelFileSettings[] :=
    Module[{},
           FlexibleSUSY`HighPoleMassPrecision = DeleteCases[{TreeMasses`GetHiggsBoson[], TreeMasses`GetChargedHiggsBoson[], TreeMasses`GetPseudoscalarHiggsBoson[]}, Null];
    ];

CheckModelFileSettings[] :=
    Module[{},
           (* FlexibleSUSY model name *)
           If[!ValueQ[FlexibleSUSY`FSModelName] || Head[FlexibleSUSY`FSModelName] =!= String,
              Utils`FSFancyWarning[
                 "FlexibleSUSY`FSModelName not defined!",
                 " I'm using Model`Name from SARAH: ", Model`Name
              ];
              FlexibleSUSY`FSModelName = Model`Name;
             ];
           (* Set OnlyLowEnergyFlexibleSUSY to False by default *)
           If[!ValueQ[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY] ||
              (FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= True &&
               FlexibleSUSY`OnlyLowEnergyFlexibleSUSY =!= False),
              FlexibleSUSY`OnlyLowEnergyFlexibleSUSY = False;
             ];
           If[Head[FlexibleSUSY`InitialGuessAtLowScale] =!= List,
              FlexibleSUSY`InitialGuessAtLowScale = {};
             ];
           If[Head[FlexibleSUSY`InitialGuessAtSUSYScale] =!= List,
              FlexibleSUSY`InitialGuessAtSUSYScale = {};
             ];
           If[Head[FlexibleSUSY`InitialGuessAtHighScale] =!= List,
              FlexibleSUSY`InitialGuessAtHighScale = {};
             ];
           (* HighScale *)
           If[!ValueQ[FlexibleSUSY`HighScale],
              If[!FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 Utils`FSFancyWarning[
                    "FlexibleSUSY`HighScale should be set in the model file!"
                 ];
              ];
              FlexibleSUSY`HighScale := 2 10^16;
             ];
           If[!ValueQ[FlexibleSUSY`HighScaleFirstGuess],
              If[!FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 Utils`FSFancyWarning[
                    "FlexibleSUSY`HighScaleFirstGuess should be set in the model file!"
                 ];
              ];
              FlexibleSUSY`HighScaleFirstGuess = 2.0 10^16;
             ];
           If[Head[FlexibleSUSY`HighScaleInput] =!= List,
              FlexibleSUSY`HighScaleInput = {};
             ];
           (* LowScale *)
           If[!ValueQ[FlexibleSUSY`LowScale],
              Utils`FSFancyWarning[
                 "FlexibleSUSY`LowScale should be set in the model file!"
              ];
              FlexibleSUSY`LowScale := LowEnergyConstant[MZ];
             ];
           If[!ValueQ[FlexibleSUSY`LowScaleFirstGuess],
              Utils`FSFancyWarning[
                 "FlexibleSUSY`LowScaleFirstGuess should be set in the model file!"
              ];
              FlexibleSUSY`LowScaleFirstGuess = LowEnergyConstant[MZ];
             ];
           If[Head[FlexibleSUSY`LowScaleInput] =!= List,
              FlexibleSUSY`LowScaleInput = {};
             ];
           (* SUSYScale *)
           If[!ValueQ[FlexibleSUSY`SUSYScale],
              Utils`FSFancyWarning[
                 "FlexibleSUSY`SUSYScale should be set in the model file!"
              ];
              FlexibleSUSY`SUSYScale := 1000;
             ];
           If[!ValueQ[FlexibleSUSY`SUSYScaleFirstGuess],
              Utils`FSFancyWarning[
                 "FlexibleSUSY`SUSYScaleFirstGuess should be set in the model file!"
              ];
              FlexibleSUSY`SUSYScaleFirstGuess = 1000;
             ];
           If[Head[FlexibleSUSY`SUSYScaleInput] =!= List,
              FlexibleSUSY`SUSYScaleInput = {};
             ];
           If[Head[FlexibleSUSY`SUSYScaleMatching] === List &&
              Head[FlexibleSUSY`MatchingScaleInput] =!= List,
              Utils`FSFancyWarning[
                 "SUSYScaleMatching is deprecated. Please use MatchingScaleInput instead!"
              ];
              FlexibleSUSY`MatchingScaleInput = FlexibleSUSY`SUSYScaleMatching;
             ];

           If[Head[SARAH`MINPAR] =!= List,
              SARAH`MINPAR = {};
             ];
           If[Head[SARAH`EXTPAR] =!= List,
              SARAH`EXTPAR = {};
             ];
           If[Head[IMMINPAR] =!= List,
              IMMINPAR = {};
             ];
           If[Head[IMEXTPAR] =!= List,
              IMEXTPAR = {};
             ];
           If[Head[FlexibleSUSY`TreeLevelEWSBSolution] =!= List,
              FlexibleSUSY`TreeLevelEWSBSolution = {};
             ];
           If[Head[FlexibleSUSY`ExtraSLHAOutputBlocks] =!= List,
              FlexibleSUSY`ExtraSLHAOutputBlocks = {};
             ];
           If[MemberQ[FlexibleSUSY`ExtraSLHAOutputBlocks, {FlexibleSUSY`EFFHIGGSCOUPLINGS, __}],
              Utils`FSFancyWarning[
                 "Effective coupling module has been disabled since v2.6.0.",
                 " Please use FlexibleDecay instead."
              ];
              FlexibleSUSY`ExtraSLHAOutputBlocks =
                 DeleteCases[FlexibleSUSY`ExtraSLHAOutputBlocks, {FlexibleSUSY`EFFHIGGSCOUPLINGS, __}];
           ];
           If[Head[FlexibleSUSY`EWSBOutputParameters] =!= List,
              Print["Error: EWSBOutputParameters has to be set to a list",
                    " of model parameters chosen to be output of the EWSB eqs."];
              Quit[1];
             ];
           If[Head[FlexibleSUSY`EWSBInitialGuess] =!= List,
              FlexibleSUSY`EWSBInitialGuess = {};
             ];
           If[Head[FlexibleSUSY`EWSBSubstitutions] =!= List,
              FlexibleSUSY`EWSBSubstitutions = {};
             ];
           If[ValueQ[FlexibleSUSY`FSExtraInputParameters],
              Print["Error: the use of FSExtraInputParameters is no longer supported!"];
              Print["   Please add the entries in FSExtraInputParameters to"];
              Print["   the parameters defined in FSAuxiliaryParameterInfo."];
              Quit[1];
             ];
           If[Head[FlexibleSUSY`FSAuxiliaryParameterInfo] =!= List,
              Print["Error: FSAuxiliaryParameterInfo has to be set to a list!"];
              Quit[1];
              ,
              If[!(And @@ (MatchQ[#,{_, {__}}]& /@ FlexibleSUSY`FSAuxiliaryParameterInfo)),
                 Print["Error: FSAuxiliaryParameterInfo must be of the form",
                       " {{par, {property -> value, ...}}, ... }"];
                ];
             ];
           If[FlexibleSUSY`FlexibleEFTHiggs === True && HaveBVPSolver[FlexibleSUSY`SemiAnalyticSolver],
              Print["Error: the use of FlexibleEFTHiggs with the semi-analytic solver"];
              Print["   is not yet supported.  Please either set FlexibleEFTHiggs = False"];
              Print["   or remove the entry SemiAnalyticSolver from FSBVPSolvers in the"];
              Print["   model file."];
              Quit[1];
             ];
           CheckEWSBSolvers[FlexibleSUSY`FSEWSBSolvers];
           CheckBVPSolvers[FlexibleSUSY`FSBVPSolvers];
           CheckDecaysOptions[];
           ReplaceSymbolsInUserInput[{Susyno`LieGroups`M -> FlexibleSUSY`M}];
          ];

CheckExtraParametersUsage[parameters_List, boundaryConditions_List] :=
    Module[{usedCases, multiplyUsedPars},
           usedCases = Function[par, !FreeQ[#, par]& /@ boundaryConditions] /@ parameters;
           multiplyUsedPars = Position[Count[#, True]& /@ usedCases, n_ /; n > 1];
           If[multiplyUsedPars =!= {},
              Utils`FSFancyWarning[
                 "The following auxiliary parameters appear at multiple",
                 " scales, but do not run: ", Extract[parameters, multiplyUsedPars]
              ];
             ];
          ];

ReplaceIndicesInUserInput[rules_] :=
    Block[{},
          FlexibleSUSY`InitialGuessAtLowScale  = FlexibleSUSY`InitialGuessAtLowScale  /. rules;
          FlexibleSUSY`InitialGuessAtSUSYScale = FlexibleSUSY`InitialGuessAtSUSYScale /. rules;
          FlexibleSUSY`InitialGuessAtHighScale = FlexibleSUSY`InitialGuessAtHighScale /. rules;
          FlexibleSUSY`HighScale               = FlexibleSUSY`HighScale               /. rules;
          FlexibleSUSY`HighScaleFirstGuess     = FlexibleSUSY`HighScaleFirstGuess     /. rules;
          FlexibleSUSY`HighScaleInput          = FlexibleSUSY`HighScaleInput          /. rules;
          FlexibleSUSY`LowScale                = FlexibleSUSY`LowScale                /. rules;
          FlexibleSUSY`LowScaleFirstGuess      = FlexibleSUSY`LowScaleFirstGuess      /. rules;
          FlexibleSUSY`LowScaleInput           = FlexibleSUSY`LowScaleInput           /. rules;
          FlexibleSUSY`SUSYScale               = FlexibleSUSY`SUSYScale               /. rules;
          FlexibleSUSY`SUSYScaleFirstGuess     = FlexibleSUSY`SUSYScaleFirstGuess     /. rules;
          FlexibleSUSY`SUSYScaleInput          = FlexibleSUSY`SUSYScaleInput          /. rules;
         ];

EvaluateUserInput[] :=
    Block[{},
          FlexibleSUSY`HighScaleInput          = Map[Evaluate, FlexibleSUSY`HighScaleInput, {0,Infinity}];
          FlexibleSUSY`LowScaleInput           = Map[Evaluate, FlexibleSUSY`LowScaleInput , {0,Infinity}];
          FlexibleSUSY`SUSYScaleInput          = Map[Evaluate, FlexibleSUSY`SUSYScaleInput, {0,Infinity}];
         ];

GUTNormalization[coupling_] :=
    Parameters`GetGUTNormalization[coupling];

ParticleIndexRule[par_, name_String] := {
    "@" <> name <> "@" -> CConversion`ToValidCSymbolString[par],
    "@" <> name <> "_" ~~ num___ ~~ "@" /; StringFreeQ[num, "@"] :>
    CConversion`ToValidCSymbolString[par] <> If[TreeMasses`GetDimension[par] > 1, "(" <> num <> ")", ""],
    "@" <> name <> "(" ~~ num___ ~~ ")@" /; StringFreeQ[num, "@"] :>
    CConversion`ToValidCSymbolString[par] <> If[TreeMasses`GetDimension[par] > 1, "(" <> num <> ")", "()"]
};

GenerationIndexRule[par_, name_String] :=
    "@Generations(" ~~ name ~~ ")@" :>
    ToString[TreeMasses`GetDimension[par]];

GeneralReplacementRules[] :=
    Join[
    ParticleIndexRule[SARAH`VectorZ, "VectorZ"],
    ParticleIndexRule[SARAH`VectorW, "VectorW"],
    ParticleIndexRule[SARAH`VectorP, "VectorP"],
    ParticleIndexRule[SARAH`VectorG, "VectorG"],
    ParticleIndexRule[SARAH`TopQuark, "TopQuark"],
    ParticleIndexRule[SARAH`BottomQuark, "BottomQuark"],
    ParticleIndexRule[SARAH`Electron, "Electron"],
    ParticleIndexRule[SARAH`Neutrino, "Neutrino"],
    ParticleIndexRule[SARAH`HiggsBoson, "HiggsBoson"],
    ParticleIndexRule[SARAH`PseudoScalarBoson, "PseudoScalarBoson"],
    ParticleIndexRule[SARAH`ChargedHiggs, "ChargedHiggs"],
    ParticleIndexRule[SARAH`TopSquark, "TopSquark"],
    ParticleIndexRule[SARAH`BottomSquark, "BottomSquark"],
    ParticleIndexRule[SARAH`Sneutrino, "Sneutrino"],
    ParticleIndexRule[SARAH`Selectron, "Selectron"],
    ParticleIndexRule[SARAH`Gluino, "Gluino"],
    {
        GenerationIndexRule[SARAH`VectorZ, "VectorZ"],
        GenerationIndexRule[SARAH`VectorW, "VectorW"],
        GenerationIndexRule[SARAH`VectorP, "VectorP"],
        GenerationIndexRule[SARAH`VectorG, "VectorG"],
        GenerationIndexRule[SARAH`TopQuark, "TopQuark"],
        GenerationIndexRule[SARAH`BottomQuark, "BottomQuark"],
        GenerationIndexRule[SARAH`Electron, "Electron"],
        GenerationIndexRule[SARAH`Neutrino, "Neutrino"],
        GenerationIndexRule[SARAH`HiggsBoson, "HiggsBoson"],
        GenerationIndexRule[SARAH`PseudoScalarBoson, "PseudoScalarBoson"],
        GenerationIndexRule[SARAH`ChargedHiggs, "ChargedHiggs"],
        GenerationIndexRule[SARAH`TopSquark, "TopSquark"],
        GenerationIndexRule[SARAH`BottomSquark, "BottomSquark"],
        GenerationIndexRule[SARAH`Sneutrino, "Sneutrino"],
        GenerationIndexRule[SARAH`Selectron, "Selectron"],
        GenerationIndexRule[SARAH`Gluino, "Gluino"]
    },
    { "@UpYukawa@"       -> CConversion`ToValidCSymbolString[SARAH`UpYukawa],
      "@DownYukawa@"     -> CConversion`ToValidCSymbolString[SARAH`DownYukawa],
      "@ElectronYukawa@" -> CConversion`ToValidCSymbolString[SARAH`ElectronYukawa],
      "@TrilinearUp@"          -> CConversion`ToValidCSymbolString[SARAH`TrilinearUp],
      "@TrilinearDown@"        -> CConversion`ToValidCSymbolString[SARAH`TrilinearDown],
      "@TrilinearLepton@"      -> CConversion`ToValidCSymbolString[SARAH`TrilinearLepton],
      "@LeftUpMixingMatrix@"   -> CConversion`ToValidCSymbolString[SARAH`UpMatrixL],
      "@LeftDownMixingMatrix@" -> CConversion`ToValidCSymbolString[SARAH`DownMatrixL],
      "@RightUpMixingMatrix@"  -> CConversion`ToValidCSymbolString[SARAH`UpMatrixR],
      "@RightDownMixingMatrix@"-> CConversion`ToValidCSymbolString[SARAH`DownMatrixR],
      "@hyperchargeCoupling@" -> CConversion`ToValidCSymbolString[SARAH`hyperchargeCoupling],
      "@leftCoupling@"        -> CConversion`ToValidCSymbolString[SARAH`leftCoupling],
      "@strongCoupling@"      -> CConversion`ToValidCSymbolString[SARAH`strongCoupling],
      "@hyperchargeCouplingGutNormalization@"  -> CConversion`RValueToCFormString[Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]],
      "@leftCouplingGutNormalization@"  -> CConversion`RValueToCFormString[Parameters`GetGUTNormalization[SARAH`leftCoupling]],
      "@hyperchargeCouplingInverseGutNormalization@" -> CConversion`RValueToCFormString[1/Parameters`GetGUTNormalization[SARAH`hyperchargeCoupling]],
      "@leftCouplingInverseGutNormalization@" -> CConversion`RValueToCFormString[1/Parameters`GetGUTNormalization[SARAH`leftCoupling]],
      "@perturbativityThreshold@" -> ToString[N[FlexibleSUSY`FSPerturbativityThreshold]],
      "@ModelName@"           -> FlexibleSUSY`FSModelName,
      "@numberOfModelParameters@" -> ToString[numberOfModelParameters],
      "@numberOfParticles@"    -> ToString[Length @ GetLoopCorrectedParticles[FlexibleSUSY`FSEigenstates]],
      "@numberOfSMParticles@"  -> ToString[Length @ Select[GetLoopCorrectedParticles[FlexibleSUSY`FSEigenstates], TreeMasses`IsSMParticle]],
      "@numberOfBSMParticles@" -> ToString[Length @ Complement[GetLoopCorrectedParticles[FlexibleSUSY`FSEigenstates],
                                                               Select[GetLoopCorrectedParticles[FlexibleSUSY`FSEigenstates], TreeMasses`IsSMParticle]]],
      "@InputParameter_" ~~ num_ ~~ "@" /; IntegerQ[ToExpression[num]] :> CConversion`ToValidCSymbolString[
          If[Parameters`GetInputParameters[] === {},
             "",
             Parameters`GetInputParameters[][[ToExpression[num]]]
            ]
      ],
      "@setInputParameterTo[" ~~ num_ ~~ "," ~~ value__ ~~ "]@" /; IntegerQ[ToExpression[num]] :>
          If[Parameters`GetInputParameters[] === {},
             "",
             IndentText[IndentText[
                 Parameters`SetInputParameter[
                     Parameters`GetInputParameters[][[ToExpression[num]]],
                     value,
                     "INPUTPARAMETER"
                 ]
             ]]
            ],
      "@RenScheme@"           -> ToString[FlexibleSUSY`FSRenormalizationScheme],
      "@ModelTypes@"          -> FlexibleTower`GetModelTypes[],
      "@DateAndTime@"         -> DateString[],
      "@SARAHVersion@"        -> SA`Version,
      "@FlexibleSUSYVersion@" -> FS`Version,
      "@FlexibleSUSYGitCommit@" -> FS`GitCommit
    }
    ];


WriteRGEClass[betaFun_List, anomDim_List, files_List,
              templateFile_String, makefileModuleTemplates_List,
              additionalTraces_List:{}, numberOfBaseClassParameters_:0] :=
   Module[{beta, setter, getter, parameterDef, set,
           display,
           cCtorParameterList, parameterCopyInit, betaParameterList,
           anomDimPrototypes, anomDimFunctions, printParameters, parameters,
           numberOfParameters, clearParameters,
           singleBetaFunctionsDecls, singleBetaFunctionsDefsFiles,
           traceDefs, calcTraces, sarahTraces},
          (* extract list of parameters from the beta functions *)
          parameters = BetaFunction`GetName[#]& /@ betaFun;
          (* count number of parameters *)
          numberOfParameters = BetaFunction`CountNumberOfParameters[betaFun] + numberOfBaseClassParameters;
          (* create C++ functions and parameter declarations *)
          sarahTraces          = Traces`ConvertSARAHTraces[additionalTraces];
          beta                 = BetaFunction`CreateBetaFunction[betaFun];
          setter               = BetaFunction`CreateSetters[betaFun];
          getter               = BetaFunction`CreateGetters[betaFun];
          parameterDef         = BetaFunction`CreateParameterDefinitions[betaFun];
          set                  = BetaFunction`CreateSetFunction[betaFun, numberOfBaseClassParameters];
          display              = BetaFunction`CreateDisplayFunction[betaFun, numberOfBaseClassParameters];
          cCtorParameterList   = BetaFunction`CreateCCtorParameterList[betaFun];
          parameterCopyInit    = BetaFunction`CreateCCtorInitialization[betaFun];
          betaParameterList    = BetaFunction`FSCreateParameterList[betaFun, "beta_"];
          clearParameters      = BetaFunction`ClearParameters[betaFun];
          anomDimPrototypes    = AnomalousDimension`CreateAnomDimPrototypes[anomDim];
          anomDimFunctions     = AnomalousDimension`CreateAnomDimFunctions[anomDim];
          printParameters      = WriteOut`PrintParameters[parameters, "ostr"];
          singleBetaFunctionsDecls = BetaFunction`CreateSingleBetaFunctionDecl[betaFun];
          traceDefs            = Traces`CreateTraceDefs[betaFun];
          traceDefs            = traceDefs <> Traces`CreateSARAHTraceDefs[sarahTraces];
          calcTraces           = {Traces`CreateSARAHTraceCalculation[sarahTraces, "TRACE_STRUCT"],
                                  Sequence @@ Traces`CreateTraceCalculation[betaFun, "TRACE_STRUCT"] };
          WriteOut`ReplaceInFiles[files,
                 { "@beta@"                 -> IndentText[WrapLines[beta]],
                   "@clearParameters@"      -> IndentText[WrapLines[clearParameters]],
                   "@display@"              -> IndentText[display],
                   "@set@"                  -> IndentText[set],
                   "@cCtorParameterList@"   -> WrapLines[cCtorParameterList],
                   "@parameterCopyInit@"    -> WrapLines[parameterCopyInit],
                   "@betaParameterList@"    -> betaParameterList,
                   "@parameterDef@"         -> IndentText[parameterDef],
                   "@cCtorParameterList@"   -> WrapLines[cCtorParameterList],
                   "@setter@"               -> IndentText[setter],
                   "@getter@"               -> IndentText[getter],
                   "@anomDimPrototypes@"    -> IndentText[anomDimPrototypes],
                   "@anomDimFunctions@"     -> WrapLines[anomDimFunctions],
                   "@numberOfParameters@"   -> RValueToCFormString[numberOfParameters],
                   "@printParameters@"      -> IndentText[printParameters],
                   "@singleBetaFunctionsDecls@" -> IndentText[singleBetaFunctionsDecls],
                   "@traceDefs@"            -> IndentText[IndentText[traceDefs]],
                   "@calc1LTraces@"         -> IndentText @ IndentText[WrapLines[calcTraces[[1]] <> "\n" <> calcTraces[[2]]]],
                   "@calc2LTraces@"         -> IndentText @ IndentText[WrapLines[calcTraces[[3]]]],
                   "@calc3LTraces@"         -> IndentText @ IndentText[WrapLines[calcTraces[[4]]]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          singleBetaFunctionsDefsFiles = BetaFunction`CreateSingleBetaFunctionDefs[betaFun, templateFile, sarahTraces];
          Print["Creating makefile module for the beta functions ..."];
          WriteMakefileModule[singleBetaFunctionsDefsFiles,
                              makefileModuleTemplates];
         ];

WriteInputParameterClass[inputParameters_List, files_List] :=
   Module[{defineInputParameters, printInputParameters, get, set, inputPars},
          inputPars = {First[#], #[[3]]}& /@ inputParameters;
          defineInputParameters = Constraint`DefineInputParameters[inputPars];
          printInputParameters = WriteOut`PrintInputParameters[inputPars,"ostr"];
          get = Parameters`CreateInputParameterArrayGetter[inputPars];
          set = Parameters`CreateInputParameterArraySetter[inputPars];
          WriteOut`ReplaceInFiles[files,
                         { "@defineInputParameters@" -> IndentText[defineInputParameters],
                           "@printInputParameters@"       -> IndentText[printInputParameters],
                           "@get@"                        -> IndentText[get],
                           "@set@"                        -> IndentText[set],
                        Sequence @@ GeneralReplacementRules[]
                      } ];
       ];

WriteConstraintClass[condition_, settings_List, scaleFirstGuess_,
                  {minimumScale_, maximumScale_}, files_List] :=
   Module[{applyConstraint = "", calculateScale, scaleGuess,
           restrictScale,
           initialSetting,
           setDRbarYukawaCouplings,
           calculateDRbarMasses,
           calculateDeltaAlphaEm, calculateDeltaAlphaS,
           calculateGaugeCouplings,
           calculateThetaW, calculateSMHiggsPoleMass,
           fillHimalayaInput,
           checkPerturbativityForDimensionlessParameters = "",
           twoLoopThresholdHeaders = "" },
          Constraint`SetBetaFunctions[GetBetaFunctions[]];
          applyConstraint = Constraint`ApplyConstraints[settings];
          calculateScale  = Constraint`CalculateScale[condition, "scale"];
          scaleGuess      = Constraint`CalculateScale[scaleFirstGuess, "initial_scale_guess"];
          restrictScale   = Constraint`RestrictScale[{minimumScale, maximumScale}];
          initialSetting  = Constraint`InitialApplyConstraint[settings];
          calculateDeltaAlphaEm   = ThresholdCorrections`CalculateDeltaAlphaEm[FlexibleSUSY`FSRenormalizationScheme];
          calculateDeltaAlphaS    = ThresholdCorrections`CalculateDeltaAlphaS[FlexibleSUSY`FSRenormalizationScheme];
          calculateGaugeCouplings = ThresholdCorrections`CalculateGaugeCouplings[];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[settings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[settings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[settings]
          };
          calculateDRbarMasses = {
              LoopMasses`CallCalculateDRbarMass["Up Quark"         , "Up-Quarks"  , 1, "upQuarksDRbar", "qedqcd.displayMass(softsusy::mUp)"      ],
              LoopMasses`CallCalculateDRbarMass["Charmed Quark"    , "Up-Quarks"  , 2, "upQuarksDRbar", "qedqcd.displayMass(softsusy::mCharm)"   ],
              LoopMasses`CallCalculateDRbarMass["Top Quark"        , "Up-Quarks"  , 3, "upQuarksDRbar", "qedqcd.displayPoleMt()"                 ],
              LoopMasses`CallCalculateDRbarMass["Down Quark"       , "Down-Quarks", 1, "downQuarksDRbar", "qedqcd.displayMass(softsusy::mDown)"    ],
              LoopMasses`CallCalculateDRbarMass["Strange Quark"    , "Down-Quarks", 2, "downQuarksDRbar", "qedqcd.displayMass(softsusy::mStrange)" ],
              LoopMasses`CallCalculateDRbarMass["Bottom Quark"     , "Down-Quarks", 3, "downQuarksDRbar", "qedqcd.displayMass(softsusy::mBottom)"  ],
              LoopMasses`CallCalculateDRbarMass["Electron"         , "Leptons"    , 1, "downLeptonsDRbar", "qedqcd.displayMass(softsusy::mElectron)"],
              LoopMasses`CallCalculateDRbarMass["Muon"             , "Leptons"    , 2, "downLeptonsDRbar", "qedqcd.displayMass(softsusy::mMuon)"    ],
              LoopMasses`CallCalculateDRbarMass["Tau"              , "Leptons"    , 3, "downLeptonsDRbar", "qedqcd.displayMass(softsusy::mTau)"     ],
              LoopMasses`CallCalculateDRbarMass["Electron Neutrino", "Neutrinos"  , 1, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(1)"      ],
              LoopMasses`CallCalculateDRbarMass["Muon Neutrino"    , "Neutrinos"  , 2, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(2)"      ],
              LoopMasses`CallCalculateDRbarMass["Tau Neutrino"     , "Neutrinos"  , 3, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(3)"      ]
          };
          If[FSCheckPerturbativityOfDimensionlessParameters,
             checkPerturbativityForDimensionlessParameters =
                 Constraint`CheckPerturbativityForParameters[
                     Parameters`GetModelParametersWithMassDimension[0],
                     FlexibleSUSY`FSPerturbativityThreshold
                 ];
            ];
          If[FlexibleSUSY`FSWeakMixingAngleInput === FlexibleSUSY`FSFermiConstant &&
             !WeinbergAngle`CheckMuonDecayRunning[],
             Print["Error: Cannot use Fermi constant to determine weak mixing angle!"];
             Print["   Please use FSWeakMixingAngleInput = Automatic or FSMassW"];
             Quit[1];
            ];
          If[FlexibleSUSY`FSWeakMixingAngleInput === Automatic,
             If[WeinbergAngle`CheckMuonDecayRunning[],
                FlexibleSUSY`FSWeakMixingAngleInput = FlexibleSUSY`FSFermiConstant,
                FlexibleSUSY`FSWeakMixingAngleInput = FlexibleSUSY`FSMassW
               ];
            ];
          calculateThetaW   = ThresholdCorrections`CalculateThetaW[FlexibleSUSY`FSWeakMixingAngleInput];
          calculateSMHiggsPoleMass = LoopMasses`CalculateSMHiggsPoleMass[FlexibleSUSY`FSWeakMixingAngleInput];
          fillHimalayaInput = Himalaya`FillHimalayaInput[FSHimalayaInput];
          twoLoopThresholdHeaders = ThresholdCorrections`GetTwoLoopThresholdHeaders[];
          WriteOut`ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@calculateScale@"       -> IndentText[WrapLines[calculateScale]],
                   "@scaleGuess@"           -> IndentText[WrapLines[scaleGuess]],
                   "@restrictScale@"        -> IndentText[WrapLines[restrictScale]],
                   "@initialSetting@"       -> IndentText[WrapLines[initialSetting]],
                   "@calculateGaugeCouplings@" -> IndentText[WrapLines[calculateGaugeCouplings]],
                   "@calculateDeltaAlphaEm@" -> IndentText[WrapLines[calculateDeltaAlphaEm]],
                   "@calculateDeltaAlphaS@"  -> IndentText[WrapLines[calculateDeltaAlphaS]],
                   "@calculateThetaW@"       -> IndentText[WrapLines[calculateThetaW]],
                   "@calculateSMHiggsPoleMass@" -> IndentText[WrapLines[calculateSMHiggsPoleMass]],
                   "@calculateDRbarMassUp@"      -> IndentText[IndentText[calculateDRbarMasses[[1]]]],
                   "@calculateDRbarMassCharm@"   -> IndentText[IndentText[calculateDRbarMasses[[2]]]],
                   "@calculateDRbarMassTop@"     -> IndentText[IndentText[calculateDRbarMasses[[3]]]],
                   "@calculateDRbarMassDown@"    -> IndentText[IndentText[calculateDRbarMasses[[4]]]],
                   "@calculateDRbarMassStrange@" -> IndentText[IndentText[calculateDRbarMasses[[5]]]],
                   "@calculateDRbarMassBottom@"  -> IndentText[IndentText[calculateDRbarMasses[[6]]]],
                   "@calculateDRbarMassElectron@"-> IndentText[IndentText[calculateDRbarMasses[[7]]]],
                   "@calculateDRbarMassMuon@"    -> IndentText[IndentText[calculateDRbarMasses[[8]]]],
                   "@calculateDRbarMassTau@"     -> IndentText[IndentText[calculateDRbarMasses[[9]]]],
                   "@calculateDRbarMassElectronNeutrino@"-> IndentText[IndentText[calculateDRbarMasses[[10]]]],
                   "@calculateDRbarMassMuonNeutrino@"    -> IndentText[IndentText[calculateDRbarMasses[[11]]]],
                   "@calculateDRbarMassTauNeutrino@"     -> IndentText[IndentText[calculateDRbarMasses[[12]]]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   "@fillHimalayaInput@"                -> IndentText[IndentText[fillHimalayaInput]],
                   "@checkPerturbativityForDimensionlessParameters@" -> IndentText[checkPerturbativityForDimensionlessParameters],
                   "@twoLoopThresholdHeaders@" -> twoLoopThresholdHeaders,
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteSemiAnalyticConstraintClass[condition_, settings_List, initialGuessSettings_List,
                                 scaleFirstGuess_, {minimumScale_, maximumScale_},
                                 isBoundaryConstraint_, isSemiAnalyticConstraint_, semiAnalyticSolns_List, files_List] :=
   Module[{innerSettings = {}, outerSettings = {}, innerInitialGuessSettings = {},
           applyConstraint = "", applyOuterConstraint = "", calculateScale, scaleGuess,
           restrictScale, calculateOuterScale, outerScaleGuess, restrictOuterScale,
           initialSetting, temporaryResetting = "",
           setDRbarYukawaCouplings,
           calculateDRbarMasses,
           calculateDeltaAlphaEm, calculateDeltaAlphaS,
           calculateGaugeCouplings,
           calculateThetaW, calculateSMHiggsPoleMass,
           checkPerturbativityForDimensionlessParameters = "",
           saveBoundaryValueParameters = "",
           usingSemiAnalyticScaleGetter = "",
           setSemiAnalyticScaleGetter = "",
           semiAnalyticScaleGetter = "",
           getConstraintScale = "return scale;",
           twoLoopThresholdHeaders = "" },
          Constraint`SetBetaFunctions[GetBetaFunctions[]];
          innerSettings = Select[settings, (!SemiAnalytic`IsSemiAnalyticSetting[#]
                                            && !SemiAnalytic`IsBasisParameterSetting[#, semiAnalyticSolns])&];
          If[isSemiAnalyticConstraint,
             innerSettings = DeleteCases[innerSettings, FlexibleSUSY`FSSolveEWSBFor[___]];
             outerSettings = Cases[settings, FlexibleSUSY`FSSolveEWSBFor[___]];
             applyOuterConstraint = Constraint`ApplyConstraints[outerSettings];
            ];
          applyConstraint = Constraint`ApplyConstraints[innerSettings];
          If[initialGuessSettings =!= {},
             innerInitialGuessSettings = Select[initialGuessSettings, (!SemiAnalytic`IsSemiAnalyticSetting[#]
                                                                       && !SemiAnalytic`IsBasisParameterSetting[#, semiAnalyticSolns])&];
             applyConstraint = applyConstraint <> "if (is_initial_guess) {\n";
             applyConstraint = applyConstraint <> IndentText[Constraint`ApplyConstraints[innerInitialGuessSettings]]
                               <> "}\n";
            ];
          calculateScale  = Constraint`CalculateScale[condition, "scale"];
          scaleGuess      = Constraint`CalculateScale[scaleFirstGuess, "initial_scale_guess"];
          restrictScale   = Constraint`RestrictScale[{minimumScale, maximumScale}];
          initialSetting  = Constraint`InitialApplyConstraint[innerSettings];
          calculateDeltaAlphaEm   = ThresholdCorrections`CalculateDeltaAlphaEm[FlexibleSUSY`FSRenormalizationScheme];
          calculateDeltaAlphaS    = ThresholdCorrections`CalculateDeltaAlphaS[FlexibleSUSY`FSRenormalizationScheme];
          calculateGaugeCouplings = ThresholdCorrections`CalculateGaugeCouplings[];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[innerSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[innerSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[innerSettings]
          };
          calculateDRbarMasses = {
              LoopMasses`CallCalculateDRbarMass["Up Quark"         , "Up-Quarks"  , 1, "upQuarksDRbar", "qedqcd.displayMass(softsusy::mUp)"      ],
              LoopMasses`CallCalculateDRbarMass["Charmed Quark"    , "Up-Quarks"  , 2, "upQuarksDRbar", "qedqcd.displayMass(softsusy::mCharm)"   ],
              LoopMasses`CallCalculateDRbarMass["Top Quark"        , "Up-Quarks"  , 3, "upQuarksDRbar", "qedqcd.displayPoleMt()"                 ],
              LoopMasses`CallCalculateDRbarMass["Down Quark"       , "Down-Quarks", 1, "downQuarksDRbar", "qedqcd.displayMass(softsusy::mDown)"    ],
              LoopMasses`CallCalculateDRbarMass["Strange Quark"    , "Down-Quarks", 2, "downQuarksDRbar", "qedqcd.displayMass(softsusy::mStrange)" ],
              LoopMasses`CallCalculateDRbarMass["Bottom Quark"     , "Down-Quarks", 3, "downQuarksDRbar", "qedqcd.displayMass(softsusy::mBottom)"  ],
              LoopMasses`CallCalculateDRbarMass["Electron"         , "Leptons"    , 1, "downLeptonsDRbar", "qedqcd.displayMass(softsusy::mElectron)"],
              LoopMasses`CallCalculateDRbarMass["Muon"             , "Leptons"    , 2, "downLeptonsDRbar", "qedqcd.displayMass(softsusy::mMuon)"    ],
              LoopMasses`CallCalculateDRbarMass["Tau"              , "Leptons"    , 3, "downLeptonsDRbar", "qedqcd.displayMass(softsusy::mTau)"     ],
              LoopMasses`CallCalculateDRbarMass["Electron Neutrino", "Neutrinos"  , 1, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(1)"      ],
              LoopMasses`CallCalculateDRbarMass["Muon Neutrino"    , "Neutrinos"  , 2, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(2)"      ],
              LoopMasses`CallCalculateDRbarMass["Tau Neutrino"     , "Neutrinos"  , 3, "neutrinoDRbar", "qedqcd.displayNeutrinoPoleMass(3)"      ]
          };
          If[FSCheckPerturbativityOfDimensionlessParameters,
             checkPerturbativityForDimensionlessParameters =
                 Constraint`CheckPerturbativityForParameters[
                     Parameters`GetModelParametersWithMassDimension[0],
                     FlexibleSUSY`FSPerturbativityThreshold
                 ];
            ];
          If[FlexibleSUSY`FSWeakMixingAngleInput === FlexibleSUSY`FSFermiConstant &&
             !WeinbergAngle`CheckMuonDecayRunning[],
             Print["Error: Cannot use Fermi constant to determine weak mixing angle!"];
             Print["   Please use FSWeakMixingAngleInput = Automatic or FSMassW"];
             Quit[1];
            ];
          If[FlexibleSUSY`FSWeakMixingAngleInput === Automatic,
             If[WeinbergAngle`CheckMuonDecayRunning[],
                FlexibleSUSY`FSWeakMixingAngleInput = FlexibleSUSY`FSFermiConstant,
                FlexibleSUSY`FSWeakMixingAngleInput = FlexibleSUSY`FSMassW
               ];
            ];
          calculateThetaW   = ThresholdCorrections`CalculateThetaW[FlexibleSUSY`FSWeakMixingAngleInput];
          calculateSMHiggsPoleMass = LoopMasses`CalculateSMHiggsPoleMass[FlexibleSUSY`FSWeakMixingAngleInput];
          If[isBoundaryConstraint,
             saveBoundaryValueParameters = SemiAnalytic`SaveBoundaryValueParameters[semiAnalyticSolns];
            ];
          If[isSemiAnalyticConstraint,
             usingSemiAnalyticScaleGetter = "using Scale_getter = std::function<double()>;\n";
             setSemiAnalyticScaleGetter = "void set_scale(const Scale_getter& s) { scale_getter = s; }\n"
                                          <> "void set_scale(Scale_getter&& s) { scale_getter = std::move(s); }\n";
             semiAnalyticScaleGetter = "Scale_getter scale_getter{};\n";
             getConstraintScale = "if (!scale_getter) {\n" <> IndentText[getConstraintScale] <> "\n}\nreturn scale_getter();\n";
             calculateOuterScale  = calculateScale;
             outerScaleGuess      = scaleGuess;
             restrictOuterScale   = restrictScale;
             calculateScale = "";
             restrictScale = "";
            ];
          twoLoopThresholdHeaders = ThresholdCorrections`GetTwoLoopThresholdHeaders[];
          WriteOut`ReplaceInFiles[files,
                 { "@applyConstraint@"      -> IndentText[WrapLines[applyConstraint]],
                   "@applyOuterConstraint@" -> IndentText[WrapLines[applyOuterConstraint]],
                   "@calculateScale@"       -> IndentText[WrapLines[calculateScale]],
                   "@scaleGuess@"           -> IndentText[WrapLines[scaleGuess]],
                   "@restrictScale@"        -> IndentText[WrapLines[restrictScale]],
                   "@calculateOuterScale@"       -> IndentText[WrapLines[calculateOuterScale]],
                   "@outerScaleGuess@"           -> IndentText[WrapLines[outerScaleGuess]],
                   "@restrictOuterScale@"        -> IndentText[WrapLines[restrictOuterScale]],
                   "@initialSetting@"       -> IndentText[WrapLines[initialSetting]],
                   "@temporaryResetting@"   -> IndentText[WrapLines[temporaryResetting]],
                   "@calculateGaugeCouplings@" -> IndentText[WrapLines[calculateGaugeCouplings]],
                   "@calculateDeltaAlphaEm@" -> IndentText[WrapLines[calculateDeltaAlphaEm]],
                   "@calculateDeltaAlphaS@"  -> IndentText[WrapLines[calculateDeltaAlphaS]],
                   "@calculateThetaW@"       -> IndentText[WrapLines[calculateThetaW]],
                   "@calculateSMHiggsPoleMass@" -> IndentText[WrapLines[calculateSMHiggsPoleMass]],
                   "@calculateDRbarMassUp@"      -> IndentText[IndentText[calculateDRbarMasses[[1]]]],
                   "@calculateDRbarMassCharm@"   -> IndentText[IndentText[calculateDRbarMasses[[2]]]],
                   "@calculateDRbarMassTop@"     -> IndentText[IndentText[calculateDRbarMasses[[3]]]],
                   "@calculateDRbarMassDown@"    -> IndentText[IndentText[calculateDRbarMasses[[4]]]],
                   "@calculateDRbarMassStrange@" -> IndentText[IndentText[calculateDRbarMasses[[5]]]],
                   "@calculateDRbarMassBottom@"  -> IndentText[IndentText[calculateDRbarMasses[[6]]]],
                   "@calculateDRbarMassElectron@"-> IndentText[IndentText[calculateDRbarMasses[[7]]]],
                   "@calculateDRbarMassMuon@"    -> IndentText[IndentText[calculateDRbarMasses[[8]]]],
                   "@calculateDRbarMassTau@"     -> IndentText[IndentText[calculateDRbarMasses[[9]]]],
                   "@calculateDRbarMassElectronNeutrino@"-> IndentText[IndentText[calculateDRbarMasses[[10]]]],
                   "@calculateDRbarMassMuonNeutrino@"    -> IndentText[IndentText[calculateDRbarMasses[[11]]]],
                   "@calculateDRbarMassTauNeutrino@"     -> IndentText[IndentText[calculateDRbarMasses[[12]]]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   "@checkPerturbativityForDimensionlessParameters@" -> IndentText[checkPerturbativityForDimensionlessParameters],
                   "@saveBoundaryValueParameters@" -> IndentText[WrapLines[saveBoundaryValueParameters]],
                   "@usingSemiAnalyticScaleGetter@" -> IndentText[usingSemiAnalyticScaleGetter],
                   "@setSemiAnalyticScaleGetter@" -> IndentText[setSemiAnalyticScaleGetter],
                   "@semiAnalyticScaleGetter@" -> IndentText[semiAnalyticScaleGetter],
                   "@getConstraintScale@" -> IndentText[getConstraintScale],
                   "@twoLoopThresholdHeaders@" -> twoLoopThresholdHeaders,
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteInitialGuesserClass[lowScaleGuess_List, susyScaleGuess_List, highScaleGuess_List, files_List] :=
   Module[{initialGuessAtLowScale, initialGuessAtLowScaleGaugeCouplings = "",
           initialGuessAtHighScale, setDRbarYukawaCouplings,
           allSettings},
          initialGuessAtLowScale  = Constraint`ApplyConstraints[lowScaleGuess];
          initialGuessAtLowScaleGaugeCouplings = Constraint`InitialGuessAtLowScaleGaugeCouplings[];
          initialGuessAtSUSYScale = Constraint`ApplyConstraints[susyScaleGuess];
          initialGuessAtHighScale = Constraint`ApplyConstraints[highScaleGuess];
          allSettings             = Join[lowScaleGuess, highScaleGuess];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[allSettings]
          };
          WriteOut`ReplaceInFiles[files,
                 { "@initialGuessAtLowScale@"  -> IndentText[WrapLines[initialGuessAtLowScale]],
                   "@initialGuessAtLowScaleGaugeCouplings@" -> IndentText[WrapLines[initialGuessAtLowScaleGaugeCouplings]],
                   "@initialGuessAtSUSYScale@" -> IndentText[WrapLines[initialGuessAtSUSYScale]],
                   "@initialGuessAtHighScale@" -> IndentText[WrapLines[initialGuessAtHighScale]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteSemiAnalyticInitialGuesserClass[lowScaleGuess_List, susyScaleGuess_List, highScaleGuess_List,
                                     solutionsScale_String, files_List] :=
   Module[{initialGuessAtLowScale, initialGuessAtLowScaleGaugeCouplings = "",
           initialGuessAtHighScale, setDRbarYukawaCouplings,
           allSettings},
          initialGuessAtLowScale  = Constraint`ApplyConstraints[lowScaleGuess];
          initialGuessAtLowScaleGaugeCouplings = Constraint`InitialGuessAtLowScaleGaugeCouplings[];
          initialGuessAtSUSYScale = Constraint`ApplyConstraints[susyScaleGuess];
          initialGuessAtHighScale = Constraint`ApplyConstraints[highScaleGuess];
          allSettings             = Join[lowScaleGuess, highScaleGuess];
          setDRbarYukawaCouplings = {
              ThresholdCorrections`SetDRbarYukawaCouplingTop[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingBottom[allSettings],
              ThresholdCorrections`SetDRbarYukawaCouplingElectron[allSettings]
          };
          WriteOut`ReplaceInFiles[files,
                 { "@initialGuessAtLowScale@"  -> IndentText[WrapLines[initialGuessAtLowScale]],
                   "@initialGuessAtLowScaleGaugeCouplings@" -> IndentText[WrapLines[initialGuessAtLowScaleGaugeCouplings]],
                   "@initialGuessAtSUSYScale@" -> IndentText[WrapLines[initialGuessAtSUSYScale]],
                   "@initialGuessAtHighScale@" -> IndentText[WrapLines[initialGuessAtHighScale]],
                   "@setDRbarUpQuarkYukawaCouplings@"   -> IndentText[WrapLines[setDRbarYukawaCouplings[[1]]]],
                   "@setDRbarDownQuarkYukawaCouplings@" -> IndentText[WrapLines[setDRbarYukawaCouplings[[2]]]],
                   "@setDRbarElectronYukawaCouplings@"  -> IndentText[WrapLines[setDRbarYukawaCouplings[[3]]]],
                   "@inputScaleGuess@" -> solutionsScale,
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteConvergenceTesterClass[parameters_, files_List] :=
   Module[{compareFunction},
          compareFunction = ConvergenceTester`CreateCompareFunction[parameters];
          WriteOut`ReplaceInFiles[files,
                 { "@compareFunction@"      -> IndentText[WrapLines[compareFunction]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

WriteWeinbergAngleClass[deltaVBcontributions_List, vertexRules_List, files_List] :=
   Module[{deltaVBprototypes = "", deltaVBfunctions = "", deltaVBcalculation = "",
           scheme = GetRenormalizationScheme[]},
          {deltaVBprototypes, deltaVBfunctions} =
             WeinbergAngle`CreateDeltaVBContributions[deltaVBcontributions, vertexRules];
          deltaVBcalculation = WeinbergAngle`CreateDeltaVBCalculation[deltaVBcontributions];
          WriteOut`ReplaceInFiles[files,
                 { "@DefSMhyperCoupling@" -> WeinbergAngle`DefSMhyperCoupling[],
                   "@DefSMleftCoupling@"  -> WeinbergAngle`DefSMleftCoupling[],
                   "@DeltaRhoHat2LoopSM@" -> IndentText[IndentText[WrapLines[WeinbergAngle`DeltaRhoHat2LoopSM[]]]],
                   "@DeltaRHat2LoopSM@"   -> IndentText[IndentText[WrapLines[WeinbergAngle`DeltaRHat2LoopSM[]]]],
                   "@RhoHatTree@"         -> IndentText[WrapLines[WeinbergAngle`RhoHatTree[]]],
                   "@GetBottomMass@"      -> WeinbergAngle`GetBottomMass[],
                   "@GetTopMass@"         -> WeinbergAngle`GetTopMass[],
                   "@DefVZSelfEnergy@"    -> WeinbergAngle`DefVZSelfEnergy[],
                   "@DefVWSelfEnergy@"    -> WeinbergAngle`DefVWSelfEnergy[],
                   "@GetNeutrinoIndex@"   -> IndentText[WeinbergAngle`GetNeutrinoIndex[]],
                   "@DeltaVBprototypes@"  -> IndentText[deltaVBprototypes],
                   "@DeltaVBfunctions@"   -> deltaVBfunctions,
                   "@DeltaVBcalculation@" -> IndentText[deltaVBcalculation],
                   "@YukawaMatching@"     -> IndentText[WeinbergAngle`YukawaMatching[]],
                   "@DefVZVWSelfEnergies@" -> IndentText[WeinbergAngle`DefVZVWSelfEnergies[]],
                   "@DeltaAlphaHatBSM@"   -> IndentText[WrapLines[WeinbergAngle`DeltaAlphaHatBSM[scheme]]],
                   Sequence @@ GeneralReplacementRules[]
                 } ];
          ];

FindVEV[gauge_] :=
    Module[{vev},
           vev = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                       {_,{v_,_},{gauge,_},{p_,_},___} | {_,{v_,_},{s_,_},{gauge,_},___} :> v];
           If[vev === {},
              Print["Error: could not find VEV for gauge eigenstate ", gauge];
              Quit[1];
             ];
           First[vev]
          ];

(* returns VEV normalization w.r.t. the corresponding gauge eigenstate *)
FindVEVNormalization[gauge_] :=
    Module[{vev},
           vev = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                       {_,{v_,n_},{gauge,m_},{p_,_},___} | {_,{v_,n_},{s_,_},{gauge,m_},___} :> Abs[n/m]];
           If[vev === {},
              Print["Error: could not find VEV for gauge eigenstate ", gauge];
              Quit[1];
             ];
           First[vev]
          ];

GetDimOfVEV[vev_] :=
    Switch[SARAH`getDimParameters[vev],
           {}                         , 1,
           {0}                        , 1,
           {1}                        , 1,
           {idx_}                     , SARAH`getDimParameters[vev][[1]]
          ];

ExpandIndices[sym_, 1] := sym;

ExpandIndices[sym_, number_] :=
    Table[sym[i], {i,1,number}];

ExpandGaugeIndices[gauge_List] :=
    Flatten[ExpandIndices[#, GetDimOfVEV[FindVEV[#]]]& /@ gauge];

ExpandVEVIndices[vev_] :=
    ExpandIndices[vev, GetDimOfVEV[vev]];

(* Returns a list of three-component lists where the information is
   stored which Higgs corresponds to which EWSB eq. and whether the
   corresponding tadpole is real or imaginary (only in models with CP
   violation).

   Example: MRSSM
   In[] := CreateHiggsToEWSBEqAssociation[]
   Out[] = {{hh, 1, Re}, {hh, 2, Re}, {hh, 4, Re}, {hh, 3, Re}}

   This result means:

   EWSB eq. 1 corresponds to hh[1], the 1L tadpole[1] is real
   EWSB eq. 2 corresponds to hh[2], the 1L tadpole[2] is real
   EWSB eq. 3 corresponds to hh[4], the 1L tadpole[3] is real
   EWSB eq. 4 corresponds to hh[3], the 1L tadpole[4] is real
 *)
CreateHiggsToEWSBEqAssociation[] :=
    Module[{vevs},
           vevs = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                        {_,{v_,_},{s_,_},{p_,_},___} :> {v,s,p}];
           If[Length[vevs] == 1,
              Return[{{SARAH`HiggsBoson, 1, Re}}];
             ];
           FindPositions[es_] :=
               Module[{gaugeES, higgsGaugeES},
                      gaugeES = ExpandGaugeIndices[es];
                      (* list of gauge eigenstate fields, ordered according to Higgs mixing *)
                      higgsGaugeES = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`MatterSector],
                                           {gauge_List, {SARAH`HiggsBoson, _}} :> gauge][[1]];
                      higgsGaugeES = ExpandGaugeIndices[higgsGaugeES];
                      (* find positions of gaugeES in higgsGaugeES *)
                      {SARAH`HiggsBoson,#}& /@ (Flatten[Position[higgsGaugeES, #]& /@ gaugeES])
                     ];
           Join[Append[#,Re]& /@ FindPositions[Transpose[vevs][[3]]],
                Append[#,Re]& /@ FindPositions[Transpose[vevs][[2]]]]
          ];

WriteModelSLHAClass[massMatrices_List, files_List] :=
    Module[{k,
            slhaYukawaDef = "",
            slhaYukawaGetter = "",
            convertYukawaCouplingsToSLHA = "",
            slhaTrilinearCouplingsDef = "",
            slhaTrilinearCouplingsGetter = "",
            convertTrilinearCouplingsToSLHA = "",
            slhaSoftSquaredMassesDef = "",
            slhaSoftSquaredMassesGetter = "",
            convertSoftSquaredMassesToSLHA = "",
            slhaFerimonMixingMatricesDef = "",
            slhaFerimonMixingMatricesGetters = "",
            slhaPoleMassGetters = "",
            slhaPoleMixingMatrixGetters = "",
            calculateCKMMatrix = "",
            calculatePMNSMatrix = ""
           },
           slhaYukawaDef        = WriteOut`CreateSLHAYukawaDefinition[];
           slhaYukawaGetter     = WriteOut`CreateSLHAYukawaGetters[];
           convertYukawaCouplingsToSLHA = WriteOut`ConvertYukawaCouplingsToSLHA[];
           slhaTrilinearCouplingsDef    = WriteOut`CreateSLHATrilinearCouplingDefinition[];
           slhaTrilinearCouplingsGetter = WriteOut`CreateSLHATrilinearCouplingGetters[];
           convertTrilinearCouplingsToSLHA = WriteOut`ConvertTrilinearCouplingsToSLHA[];
           slhaSoftSquaredMassesDef    = WriteOut`CreateSLHASoftSquaredMassesDefinition[];
           slhaSoftSquaredMassesGetter = WriteOut`CreateSLHASoftSquaredMassesGetters[];
           convertSoftSquaredMassesToSLHA = WriteOut`ConvertSoftSquaredMassesToSLHA[];
           slhaFerimonMixingMatricesDef = WriteOut`CreateSLHAFermionMixingMatricesDef[];
           slhaFerimonMixingMatricesGetters = WriteOut`CreateSLHAFermionMixingMatricesGetters[];
           calculateCKMMatrix = WriteOut`CalculateCKMMatrix[];
           calculatePMNSMatrix = WriteOut`CalculatePMNSMatrix[];
           For[k = 1, k <= Length[massMatrices], k++,
               slhaPoleMassGetters         = slhaPoleMassGetters <> TreeMasses`CreateSLHAPoleMassGetter[massMatrices[[k]]];
               slhaPoleMixingMatrixGetters = slhaPoleMixingMatrixGetters <> TreeMasses`CreateSLHAPoleMixingMatrixGetter[massMatrices[[k]]];
              ];
           WriteOut`ReplaceInFiles[files,
                          { "@slhaYukawaDef@"                  -> IndentText[slhaYukawaDef],
                            "@slhaYukawaGetter@"               -> IndentText[slhaYukawaGetter],
                            "@convertYukawaCouplingsToSLHA@"   -> IndentText[convertYukawaCouplingsToSLHA],
                            "@slhaFerimonMixingMatricesDef@"   -> IndentText[slhaFerimonMixingMatricesDef],
                            "@slhaFerimonMixingMatricesGetters@" -> IndentText[slhaFerimonMixingMatricesGetters],
                            "@slhaTrilinearCouplingsDef@"      -> IndentText[slhaTrilinearCouplingsDef],
                            "@slhaTrilinearCouplingsGetter@"   -> IndentText[slhaTrilinearCouplingsGetter],
                            "@convertTrilinearCouplingsToSLHA@"-> IndentText[convertTrilinearCouplingsToSLHA],
                            "@slhaSoftSquaredMassesDef@"       -> IndentText[slhaSoftSquaredMassesDef],
                            "@slhaSoftSquaredMassesGetter@"    -> IndentText[slhaSoftSquaredMassesGetter],
                            "@convertSoftSquaredMassesToSLHA@" -> IndentText[convertSoftSquaredMassesToSLHA],
                            "@slhaPoleMassGetters@"            -> IndentText[slhaPoleMassGetters],
                            "@slhaPoleMixingMatrixGetters@"    -> IndentText[slhaPoleMixingMatrixGetters],
                            "@calculateCKMMatrix@"             -> IndentText[calculateCKMMatrix],
                            "@calculatePMNSMatrix@"             -> IndentText[calculatePMNSMatrix],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

(* Returns a list of four-component lists where the information is
   stored which VEV corresponds to which Tadpole eq. and the
   normalization w.r.t. the corresponding scalar field (last element).

   Example: MRSSM
   In[] := CreateVEVToTadpoleAssociation[]
   Out[] = {{hh, 1, vd, 1}, {hh, 2, vu, 1}, {hh, 4, vS, 1}, {hh, 3, vT, 1}}
 *)
CreateVEVToTadpoleAssociation[] :=
    Module[{association, vs, vevs, norms, i},
           vs = Cases[SARAH`DEFINITION[FlexibleSUSY`FSEigenstates][SARAH`VEVs],
                      {_,{v_,_},{s_,_},{p_,_},___} :> {v,s,p}];
           (* find VEVs associtated to the scalar/pseudoscalar gauge eigenstates *)
           vevs = Flatten @
                  Join[ExpandVEVIndices[FindVEV[#]]& /@ Transpose[vs][[3]],
                       ExpandVEVIndices[FindVEV[#]]& /@ Transpose[vs][[2]]];
           (* find corresponding norms *)
           norms = Flatten @
                  Join[Table[FindVEVNormalization[#], {i,GetDimOfVEV[FindVEV[#]]}]& /@ Transpose[vs][[3]],
                       Table[FindVEVNormalization[#], {i,GetDimOfVEV[FindVEV[#]]}]& /@ Transpose[vs][[2]]];
           association = CreateHiggsToEWSBEqAssociation[];
           {#[[1]], #[[2]], vevs[[#[[2]]]], norms[[#[[2]]]]}& /@ association
          ];

GetRenormalizationScheme[] :=
    If[SARAH`SupersymmetricModel, FlexibleSUSY`DRbar, FlexibleSUSY`MSbar];

WriteFlexibleEFTHiggsMakefileModule[files_List] :=
    Module[{source = "", header = ""},
           If[FlexibleSUSY`FlexibleEFTHiggs === True,
              source = "\t\t" <> FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_matching.cpp"}];
              header = "\t\t" <> FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_matching.hpp"}];
             ];
           WriteOut`ReplaceInFiles[files,
                  { "@FlexibleEFTHiggsSource@" -> source,
                    "@FlexibleEFTHiggsHeader@" -> header,
                    Sequence @@ GeneralReplacementRules[]
                  } ];
          ];

WriteMatchingClass[susyScaleMatching_List, massMatrices_List, files_List] :=
    Module[{scheme = GetRenormalizationScheme[], userMatching = "",
            alphaS1Lmatching = "", alphaEM1Lmatching = "",
            setBSMParameters = "", higgsMassMatrix,
            setRunningUpQuarkMasses = "", setRunningDownQuarkMasses = "",
            setRunningDownLeptonMasses = "", setYukawas = "",
            calculateMUpQuarkPole1L = "", calculateMDownQuarkPole1L = "",
            calculateMDownLeptonPole1L = "",
            calculateMHiggs2LoopShift = "throw SetupError(\"2-loop Higgs self-energy not enabled.\");",
            calculateMHiggs3LoopShift = "throw SetupError(\"3-loop Higgs self-energy not enabled.\");",
            threeLoopLambdaMatching = "throw SetupError(\"3-loop matching not enabled.\");",
            twoLoopLambdaMatching = "throw SetupError(\"2-loop matching not enabled.\");",
            setGaugeLessLimit = "",
            setYukawaLessLimit = "",
            setMSSMLimit = "",
            includeMSSMTwoLoopTopMassHeader = "",
            createSMMt2LoopFunction = ""},
           If[FlexibleSUSY`FlexibleEFTHiggs === True,
              If[Head[susyScaleMatching] === List,
                 userMatching = Constraint`ApplyConstraints[susyScaleMatching];
                ];
              alphaS1Lmatching = Parameters`CreateLocalConstRefs[ThresholdCorrections`CalculateColorCoupling[scheme]] <> "\n" <>
                                 "delta_alpha_s += alpha_s/(2.*Pi)*(" <>
                                 CConversion`RValueToCFormString[ThresholdCorrections`CalculateColorCoupling[scheme]] <> ");\n";
              alphaEM1Lmatching = Parameters`CreateLocalConstRefs[ThresholdCorrections`CalculateElectromagneticCoupling[scheme]] <> "\n" <>
                                  "delta_alpha_em += alpha_em/(2.*Pi)*(" <>
                                  CConversion`RValueToCFormString[ThresholdCorrections`CalculateElectromagneticCoupling[scheme]] <> ");\n";
              higgsMassMatrix = Select[massMatrices, (TreeMasses`GetMassEigenstate[#] === SARAH`HiggsBoson)&];
              If[higgsMassMatrix === {},
                 Print["Error: Could not find mass matrix for ", SARAH`HiggsBoson];
                 Quit[1];
                ];
              setBSMParameters                  = FlexibleEFTHiggsMatching`SetBSMParameters[susyScaleMatching, GetMassMatrix[higgsMassMatrix[[1]]], "model."];
              setRunningUpQuarkMasses           = FlexibleEFTHiggsMatching`CalculateRunningUpQuarkMasses[];
              setRunningDownQuarkMasses         = FlexibleEFTHiggsMatching`CalculateRunningDownQuarkMasses[];
              setRunningDownLeptonMasses        = FlexibleEFTHiggsMatching`CalculateRunningDownLeptonMasses[];
              setYukawas                        = ThresholdCorrections`SetDRbarYukawaCouplings[];
              setGaugeLessLimit                 = FlexibleEFTHiggsMatching`SetLimit["model.", Parameters`DecreaseIndexLiterals @ FlexibleSUSY`FSGaugeLessLimit];
              setYukawaLessLimit                = FlexibleEFTHiggsMatching`SetLimit["model.", Parameters`DecreaseIndexLiterals @ FlexibleSUSY`FSYukawaLessLimit];
	      setMSSMLimit                      = FlexibleEFTHiggsMatching`SetLimit["model.", Parameters`DecreaseIndexLiterals @ FlexibleSUSY`FSMSSMLimit];
              If[SARAH`UseHiggs2LoopMSSM === True || FlexibleSUSY`UseHiggs2LoopNMSSM === True, 
 		 twoLoopLambdaMatching = FlexibleEFTHiggsMatching`Create2LoopMatching["model_input", "sm", SARAH`HiggsBoson, "idx"];
		 calculateMHiggs2LoopShift = FlexibleEFTHiggsMatching`CalculateMHiggs2LoopShift["model", SARAH`HiggsBoson, "idx"];
	      ];
              If[FlexibleSUSY`UseHiggs3LoopMSSM === True || FlexibleSUSY`UseHiggs3LoopNMSSM === True, 
 		 threeLoopLambdaMatching = FlexibleEFTHiggsMatching`Create3LoopMatching["model_input", "sm", SARAH`HiggsBoson, "idx"];
		 calculateMHiggs3LoopShift = FlexibleEFTHiggsMatching`CalculateMHiggs3LoopShift["model", "sm", SARAH`HiggsBoson, "idx"];
		 createSMMt2LoopFunction = FlexibleEFTHiggsMatching`CreateSMMtop2LoopFunction[];
		 includeMSSMTwoLoopTopMassHeader = "#include \"mssm_twoloop_mt.hpp\"";
              ];

              calculateMUpQuarkPole1L    = FlexibleEFTHiggsMatching`CalculateMUpQuarkPole1L[];
              calculateMDownQuarkPole1L  = FlexibleEFTHiggsMatching`CalculateMDownQuarkPole1L[];
              calculateMDownLeptonPole1L = FlexibleEFTHiggsMatching`CalculateMDownLeptonPole1L[];
           ];
           WriteOut`ReplaceInFiles[files,
                       { "@alphaS1Lmatching@"        -> IndentText[WrapLines[alphaS1Lmatching]],
                         "@alphaEM1Lmatching@"       -> IndentText[WrapLines[alphaEM1Lmatching]],
                         "@setBSMParameters@"        -> IndentText[setBSMParameters],
                         "@setRunningUpQuarkMasses@" -> IndentText[setRunningUpQuarkMasses],
                         "@setRunningDownQuarkMasses@" -> IndentText[setRunningDownQuarkMasses],
                         "@setRunningDownLeptonMasses@" -> IndentText[setRunningDownLeptonMasses],
                         "@calculateMUpQuarkPole1L@"    -> IndentText[calculateMUpQuarkPole1L],
                         "@calculateMDownQuarkPole1L@"  -> IndentText[calculateMDownQuarkPole1L],
                         "@calculateMDownLeptonPole1L@" -> IndentText[calculateMDownLeptonPole1L],
                         "@setYukawas@"              -> IndentText[WrapLines[setYukawas]],
                         "@applyUserMatching@"       -> IndentText[WrapLines[userMatching]],
			 "@calculateMHiggs2LoopShift@" -> IndentText[calculateMHiggs2LoopShift],
			 "@calculateMHiggs3LoopShift@" -> IndentText[calculateMHiggs3LoopShift],
                         "@numberOfEWSBEquations@" -> ToString[TreeMasses`GetDimension[SARAH`HiggsBoson]],
                         "@threeLoopLambdaMatching@" -> IndentText[threeLoopLambdaMatching],
                         "@twoLoopLambdaMatching@" -> IndentText[twoLoopLambdaMatching],
                         "@createSMMt2LoopFunction@" -> createSMMt2LoopFunction,
                         "@setGaugeLessLimit@" -> IndentText[setGaugeLessLimit],
                         "@setYukawaLessLimit@" -> IndentText[setYukawaLessLimit],
                         "@setMSSMLimit@" -> IndentText[setMSSMLimit],
                         "@includeMSSMTwoLoopTopMassHeader@" -> includeMSSMTwoLoopTopMassHeader,
                         Sequence @@ GeneralReplacementRules[]
                       } ];
        ];

WriteEWSBSolverClass[ewsbEquations_List, parametersFixedByEWSB_List, ewsbInitialGuessValues_List,
                     ewsbSubstitutions_List, ewsbSolution_List, freePhases_List, allowedEwsbSolvers_List,
                     files_List] :=
    Module[{numberOfIndependentEWSBEquations,
            ewsbEquationsTreeLevel, independentEwsbEquationsTreeLevel,
            independentEwsbEquations, higgsToEWSBEqAssociation,
            calculateOneLoopTadpolesNoStruct = "", calculateTwoLoopTadpolesNoStruct = "",
            ewsbInitialGuess = "", solveEwsbTreeLevel = "", setTreeLevelSolution = "", EWSBSolvers = "",
            setEWSBSolution = "", fillArrayWithEWSBParameters = "",
            solveEwsbWithTadpoles = "", getEWSBParametersFromVector = "",
            setEWSBParametersFromLocalCopies = "", applyEWSBSubstitutions = "",
            setModelParametersFromEWSB = ""},
           independentEwsbEquations = EWSB`GetLinearlyIndependentEqs[ewsbEquations, parametersFixedByEWSB,
                                                                     ewsbSubstitutions];
           numberOfIndependentEWSBEquations = Length[independentEwsbEquations];
           ewsbEquationsTreeLevel = ewsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
           independentEwsbEquationsTreeLevel = independentEwsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
           If[ewsbEquations =!= Table[0, {Length[ewsbEquations]}] &&
              Length[parametersFixedByEWSB] != numberOfIndependentEWSBEquations,
              Print["Error: There are ", numberOfIndependentEWSBEquations, " independent EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
	      Quit[1];
             ];
           higgsToEWSBEqAssociation     = CreateHiggsToEWSBEqAssociation[];
           calculateOneLoopTadpolesNoStruct = SelfEnergies`FillArrayWithLoopTadpoles[1, higgsToEWSBEqAssociation, "tadpole", "+", "model."];
           If[SARAH`UseHiggs2LoopMSSM === True || FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              calculateTwoLoopTadpolesNoStruct = SelfEnergies`FillArrayWithTwoLoopTadpoles[SARAH`HiggsBoson, "tadpole", "+", "model."];
             ];
           ewsbInitialGuess             = EWSB`FillInitialGuessArray[parametersFixedByEWSB, ewsbInitialGuessValues];
           solveEwsbTreeLevel           = EWSB`CreateTreeLevelEwsbSolver[ewsbSolution /. FlexibleSUSY`tadpole[_] -> 0];
           setTreeLevelSolution         = EWSB`SetTreeLevelSolution[ewsbSolution, ewsbSubstitutions];
           EWSBSolvers                  = EWSB`CreateEWSBRootFinders[allowedEwsbSolvers];
           setEWSBSolution              = EWSB`SetEWSBSolution[parametersFixedByEWSB, freePhases, "solution", "model."];
           If[ewsbSolution =!= {},
              fillArrayWithEWSBParameters  = EWSB`FillArrayWithParameters["ewsb_parameters", parametersFixedByEWSB];
             ];
           solveEwsbWithTadpoles        = EWSB`CreateEwsbSolverWithTadpoles[ewsbSolution];
           getEWSBParametersFromVector  = EWSB`GetEWSBParametersFromVector[parametersFixedByEWSB, freePhases, "ewsb_pars"];
           setEWSBParametersFromLocalCopies = EWSB`SetEWSBParametersFromLocalCopies[parametersFixedByEWSB, "model."];
           setModelParametersFromEWSB   = EWSB`SetModelParametersFromEWSB[parametersFixedByEWSB, ewsbSubstitutions, "model."];
           applyEWSBSubstitutions       = EWSB`ApplyEWSBSubstitutions[parametersFixedByEWSB, ewsbSubstitutions];
           WriteOut`ReplaceInFiles[files,
                          { "@calculateOneLoopTadpolesNoStruct@" -> IndentText[calculateOneLoopTadpolesNoStruct],
                            "@calculateTwoLoopTadpolesNoStruct@" -> IndentText[calculateTwoLoopTadpolesNoStruct],
                            "@numberOfEWSBEquations@"-> ToString[TreeMasses`GetDimension[SARAH`HiggsBoson]],
                            "@ewsbInitialGuess@"       -> IndentText[ewsbInitialGuess],
                            "@solveEwsbTreeLevel@"           -> IndentText[WrapLines[solveEwsbTreeLevel]],
                            "@setTreeLevelSolution@"         -> IndentText[WrapLines[setTreeLevelSolution]],
                            "@EWSBSolvers@"                  -> IndentText[IndentText[EWSBSolvers]],
                            "@fillArrayWithEWSBParameters@"  -> IndentText[fillArrayWithEWSBParameters],
                            "@solveEwsbWithTadpoles@"        -> IndentText[WrapLines[solveEwsbWithTadpoles]],
                            "@getEWSBParametersFromVector@"  -> IndentText[IndentText[getEWSBParametersFromVector]],
                            "@setEWSBParametersFromLocalCopies@" -> IndentText[IndentText[setEWSBParametersFromLocalCopies]],
                            "@setEWSBSolution@"              -> IndentText[setEWSBSolution],
                            "@applyEWSBSubstitutions@"       -> IndentText[IndentText[WrapLines[applyEWSBSubstitutions]]],
                            "@setModelParametersFromEWSB@"   -> IndentText[WrapLines[setModelParametersFromEWSB]],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteSemiAnalyticEWSBSolverClass[ewsbEquations_List, parametersFixedByEWSB_List, ewsbInitialGuessValues_List,
                                 ewsbSubstitutions_List, ewsbSolution_List, freePhases_List, allowedEwsbSolvers_List,
                                 solutions_List, files_List] :=
    Module[{semiAnalyticSubs, additionalEwsbSubs, numberOfIndependentEWSBEquations,
            ewsbEquationsTreeLevel,
            independentEwsbEquations, higgsToEWSBEqAssociation,
            calculateOneLoopTadpolesNoStruct = "", calculateTwoLoopTadpolesNoStruct = "",
            ewsbInitialGuess = "", solveEwsbTreeLevel = "", setTreeLevelSolution = "", EWSBSolvers = "",
            setEWSBSolution = "", fillArrayWithEWSBParameters = "",
            solveEwsbWithTadpoles = "", getEWSBParametersFromVector = "",
            setEWSBParametersFromLocalCopies = "", applyEWSBSubstitutions = "",
            setModelParametersFromEWSB = ""},
           semiAnalyticSubs = SemiAnalytic`GetSemiAnalyticEWSBSubstitutions[solutions];
           additionalEwsbSubs = Complement[ewsbSubstitutions, semiAnalyticSubs];
           independentEwsbEquations = EWSB`GetLinearlyIndependentEqs[ewsbEquations, parametersFixedByEWSB, ewsbSubstitutions];
           numberOfIndependentEWSBEquations = Length[independentEwsbEquations];
           ewsbEquationsTreeLevel = ewsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
           If[ewsbEquations =!= Table[0, {Length[ewsbEquations]}] &&
              Length[parametersFixedByEWSB] != numberOfIndependentEWSBEquations,
              Print["Error: There are ", numberOfIndependentEWSBEquations, " independent EWSB ",
                    "equations, but you want to fix ", Length[parametersFixedByEWSB],
                    " parameters: ", parametersFixedByEWSB];
	      Quit[1];
             ];
           higgsToEWSBEqAssociation     = CreateHiggsToEWSBEqAssociation[];
           calculateOneLoopTadpolesNoStruct = SelfEnergies`FillArrayWithLoopTadpoles[1, higgsToEWSBEqAssociation, "tadpole", "+", "model."];
           If[SARAH`UseHiggs2LoopMSSM === True || FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              calculateTwoLoopTadpolesNoStruct = SelfEnergies`FillArrayWithTwoLoopTadpoles[SARAH`HiggsBoson, "tadpole", "+", "model."];
             ];
           ewsbInitialGuess             = EWSB`FillInitialGuessArray[parametersFixedByEWSB, ewsbInitialGuessValues];
           solveEwsbTreeLevel = EWSB`CreateTreeLevelEwsbSolver[ewsbSolution /. FlexibleSUSY`tadpole[_] -> 0];
           solveEwsbTreeLevel = SemiAnalytic`ReplacePreprocessorMacros[solveEwsbTreeLevel, solutions];
           setTreeLevelSolution = SemiAnalytic`SetTreeLevelEWSBSolution[ewsbSolution, solutions, additionalEwsbSubs];
           solveEwsbWithTadpoles        = EWSB`CreateEwsbSolverWithTadpoles[ewsbSolution];
           solveEwsbWithTadpoles        = SemiAnalytic`ReplacePreprocessorMacros[solveEwsbWithTadpoles, solutions];
           EWSBSolvers                  = EWSB`CreateEWSBRootFinders[allowedEwsbSolvers];
           setEWSBSolution              = EWSB`SetEWSBSolution[parametersFixedByEWSB, freePhases, "solution", "model."];
           If[ewsbSolution =!= {},
              fillArrayWithEWSBParameters  = EWSB`FillArrayWithParameters["ewsb_parameters", parametersFixedByEWSB];
             ];
           getEWSBParametersFromVector  = EWSB`GetEWSBParametersFromVector[parametersFixedByEWSB, freePhases, "ewsb_pars"];
           setEWSBParametersFromLocalCopies = EWSB`SetEWSBParametersFromLocalCopies[parametersFixedByEWSB, "model."];
           setModelParametersFromEWSB   = EWSB`SetModelParametersFromEWSB[parametersFixedByEWSB, additionalEwsbSubs, "model."];
           applyEWSBSubstitutions       = EWSB`ApplyEWSBSubstitutions[parametersFixedByEWSB, additionalEwsbSubs];
           WriteOut`ReplaceInFiles[files,
                          { "@calculateOneLoopTadpolesNoStruct@" -> IndentText[calculateOneLoopTadpolesNoStruct],
                            "@calculateTwoLoopTadpolesNoStruct@" -> IndentText[calculateTwoLoopTadpolesNoStruct],
                            "@numberOfEWSBEquations@"-> ToString[TreeMasses`GetDimension[SARAH`HiggsBoson]],
                            "@ewsbInitialGuess@"       -> IndentText[ewsbInitialGuess],
                            "@solveEwsbTreeLevel@"           -> IndentText[WrapLines[solveEwsbTreeLevel]],
                            "@setTreeLevelSolution@"         -> IndentText[WrapLines[setTreeLevelSolution]],
                            "@EWSBSolvers@"                  -> IndentText[IndentText[EWSBSolvers]],
                            "@fillArrayWithEWSBParameters@"  -> IndentText[fillArrayWithEWSBParameters],
                            "@solveEwsbWithTadpoles@"        -> IndentText[WrapLines[solveEwsbWithTadpoles]],
                            "@getEWSBParametersFromVector@"  -> IndentText[IndentText[getEWSBParametersFromVector]],
                            "@setEWSBParametersFromLocalCopies@" -> IndentText[IndentText[setEWSBParametersFromLocalCopies]],
                            "@setEWSBSolution@"              -> IndentText[setEWSBSolution],
                            "@applyEWSBSubstitutions@"       -> IndentText[IndentText[WrapLines[applyEWSBSubstitutions]]],
                            "@setModelParametersFromEWSB@"   -> IndentText[WrapLines[setModelParametersFromEWSB]],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

CreateDefaultEWSBSolverConstructor[solvers_List] :=
    Module[{i, makeSharedPointer, solver, body, init = ""},
           If[Length[solvers] > 0,
              makeSharedPointer[s_] := IndentText[WrapLines[", ewsb_solver(new " <> FlexibleSUSY`FSModelName
                                                            <> "_ewsb_solver<" <> GetBVPSolverTemplateParameter[s]
                                                            <> ">())\n"]];
              solver = solvers[[1]];
              body = makeSharedPointer[solver];
              init = init <> "#if defined(" <> GetBVPSolverEnabledMacro[solver] <> ")\n" <> body;
              For[i = 2, i <= Length[solvers], i++,
                  solver = solvers[[i]];
                  body = makeSharedPointer[solver];
                  init = init <> "#elif defined(" <> GetBVPSolverEnabledMacro[solver] <> ")\n" <> body;
                 ];
              init = init <> "#endif";
             ];
           init
          ];

ParameterAppearsExactlyOnceIn[eqs_List, par_] :=
    Length[Select[eqs, (!FreeQ[#, par])&]] === 1;

WriteModelClass[massMatrices_List, ewsbEquations_List,
                parametersFixedByEWSB_List, ewsbSubstitutions_List,
                nPointFunctions_List, vertexRules_List, phases_List,
                files_List, diagonalizationPrecision_List] :=
    Module[{ewsbEquationsTreeLevel, independentEwsbEquationsTreeLevel,
            independentEwsbEquations,
            massGetters = "", k,
            mixingMatrixGetters = "",
            slhaPoleMassGetters = "", slhaPoleMixingMatrixGetters = "",
            higgsMassGetters = "", higgsToEWSBEqAssociation,
            tadpoleEqPrototypes = "", tadpoleEqFunctions = "",
            calculateTreeLevelTadpoles = "", divideTadpoleByVEV = "",
            calculateOneLoopTadpoles = "", calculateTwoLoopTadpoles = "",
            physicalMassesDef = "", mixingMatricesDef = "",
            massCalculationPrototypes = "", massCalculationFunctions = "",
            calculateAllMasses = "",
            selfEnergyPrototypes = "", selfEnergyFunctions = "", 
            selfEnergyDerivPrototypes = "",   selfEnergyDerivFunctions = "",
            twoLoopTadpolePrototypes = "", twoLoopTadpoleFunctions = "",
            twoLoopSelfEnergyPrototypes = "", twoLoopSelfEnergyFunctions = "",
            threeLoopSelfEnergyPrototypes = "", threeLoopSelfEnergyFunctions = "",
            fourLoopSelfEnergyPrototypes = "", fourLoopSelfEnergyFunctions = "",
            secondGenerationHelperPrototypes = "", secondGenerationHelperFunctions = "",
            thirdGenerationHelperPrototypes = "", thirdGenerationHelperFunctions = "",
            phasesDefinition = "", phasesGetterSetters = "",
            extraParameterDefs = "",
            extraParameterSetters = "", extraParameterGetters = "",
            loopMassesPrototypes = "", loopMassesFunctions = "",
            runningDRbarMassesPrototypes = "", runningDRbarMassesFunctions = "",
            callAllLoopMassFunctions = "",
            callAllLoopMassFunctionsInThreads = "",
            printMasses = "", printMixingMatrices = "",
            getMixings = "", setMixings = "",
            getMasses = "", setMasses = "",
            getExtraParameters = "", setExtraParameters = "",
            masses, mixingMatrices,
            dependencePrototypes, dependenceFunctions,
            clearOutputParameters = "",
            clearPhases = "", clearExtraParameters = "",
            softScalarMasses, treeLevelEWSBOutputParameters,
            parametersToSave, saveEWSBOutputParameters,
            solveTreeLevelEWSBviaSoftHiggsMasses,
            solveEWSBTemporarily,
            copyDRbarMassesToPoleMasses = "",
            copyRunningBSMMassesToDecouplingMasses = "",
            reorderDRbarMasses = "", reorderPoleMasses = "",
            checkPoleMassesForTachyons = "",
            twoLoopHiggsHeaders = "", threeLoopHiggsHeaders = "", fourLoopHiggsHeaders = "",
            twoLoopThresholdHeaders = "",
            lspGetters = "", lspFunctions = "",
            convertMixingsToSLHAConvention = "",
            convertMixingsToHKConvention = "",
            enablePoleMassThreads = True,
            ewsbSolverHeaders = "", defaultEWSBSolverCctor = "",
            smDecouplingParameters = { SARAH`hyperchargeCoupling,
                                       SARAH`leftCoupling,
                                       SARAH`strongCoupling,
                                       SARAH`UpYukawa,
                                       SARAH`DownYukawa,
                                       SARAH`ElectronYukawa, MZDRbar,
                                       MZMSbar, MWDRbar, MWMSbar,
                                       EDRbar, EMSbar, ThetaWDRbar,
                                       VEV }
           },
           convertMixingsToSLHAConvention = WriteOut`ConvertMixingsToSLHAConvention[massMatrices];
           convertMixingsToHKConvention   = WriteOut`ConvertMixingsToHKConvention[massMatrices];
           independentEwsbEquations = EWSB`GetLinearlyIndependentEqs[ewsbEquations, parametersFixedByEWSB,
                                                                     ewsbSubstitutions];
           ewsbEquationsTreeLevel = ewsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
           independentEwsbEquationsTreeLevel = independentEwsbEquations /. FlexibleSUSY`tadpole[_] -> 0;
           For[k = 1, k <= Length[massMatrices], k++,
               massGetters          = massGetters <> TreeMasses`CreateMassGetter[massMatrices[[k]]];
               mixingMatrixGetters  = mixingMatrixGetters <> TreeMasses`CreateMixingMatrixGetter[massMatrices[[k]]];
               physicalMassesDef    = physicalMassesDef <> TreeMasses`CreatePhysicalMassDefinition[massMatrices[[k]]];
               mixingMatricesDef    = mixingMatricesDef <> TreeMasses`CreateMixingMatrixDefinition[massMatrices[[k]]];
               clearOutputParameters = clearOutputParameters <> TreeMasses`ClearOutputParameters[massMatrices[[k]]];
               copyDRbarMassesToPoleMasses = copyDRbarMassesToPoleMasses <> TreeMasses`CopyDRBarMassesToPoleMasses[massMatrices[[k]]];
               copyRunningBSMMassesToDecouplingMasses = copyRunningBSMMassesToDecouplingMasses <> TreeMasses`CopyRunningMassesFromTo[massMatrices[[k]], "OTHER", ""];
               massCalculationPrototypes = massCalculationPrototypes <> TreeMasses`CreateMassCalculationPrototype[massMatrices[[k]]];
               massCalculationFunctions  = massCalculationFunctions  <> TreeMasses`CreateMassCalculationFunction[massMatrices[[k]]];
              ];
           higgsMassGetters =
               Utils`StringZipWithSeparator[
                   TreeMasses`CreateHiggsMassGetters[SARAH`HiggsBoson,""],
                   TreeMasses`CreateHiggsMassGetters[SARAH`ChargedHiggs,""],
                   TreeMasses`CreateHiggsMassGetters[SARAH`PseudoScalar,""],
                   "\n"
               ];
           clearPhases = Phases`ClearPhases[phases];
           calculateAllMasses = TreeMasses`CallMassCalculationFunctions[massMatrices];
           tadpoleEqPrototypes = EWSB`CreateEWSBEqPrototype[SARAH`HiggsBoson];
           tadpoleEqFunctions  = EWSB`CreateEWSBEqFunction[SARAH`HiggsBoson, ewsbEquationsTreeLevel];
           higgsToEWSBEqAssociation     = CreateHiggsToEWSBEqAssociation[];
           calculateTreeLevelTadpoles   = EWSB`FillArrayWithEWSBEqs[SARAH`HiggsBoson, "tadpole"];
           calculateOneLoopTadpoles     = SelfEnergies`FillArrayWithLoopTadpoles[1, higgsToEWSBEqAssociation, "tadpole", "-"];
           divideTadpoleByVEV           = SelfEnergies`DivideTadpoleByVEV[Parameters`DecreaseIndexLiterals @ CreateVEVToTadpoleAssociation[], "tadpole"];
           If[SARAH`UseHiggs2LoopMSSM === True || FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              calculateTwoLoopTadpoles  = SelfEnergies`FillArrayWithTwoLoopTadpoles[SARAH`HiggsBoson, "tadpole", "-"];
             ];
           If[FlexibleSUSY`UseHiggs2LoopSM === True,
              {twoLoopSelfEnergyPrototypes, twoLoopSelfEnergyFunctions} = SelfEnergies`CreateTwoLoopSelfEnergiesSM[{SARAH`HiggsBoson}];
              twoLoopHiggsHeaders = "#include \"sm_twoloophiggs.hpp\"\n";
             ];
           If[FlexibleSUSY`UseHiggs3LoopSM === True,
              {threeLoopSelfEnergyPrototypes, threeLoopSelfEnergyFunctions} = SelfEnergies`CreateThreeLoopSelfEnergiesSM[{SARAH`HiggsBoson}];
              threeLoopHiggsHeaders = "#include \"sm_threeloophiggs.hpp\"\n";
             ];
           If[FlexibleSUSY`UseHiggs4LoopSM === True,
              {fourLoopSelfEnergyPrototypes, fourLoopSelfEnergyFunctions} = SelfEnergies`CreateFourLoopSelfEnergiesSM[{SARAH`HiggsBoson}];
              fourLoopHiggsHeaders = "#include \"sm_fourloophiggs.hpp\"\n";
             ];
           If[FlexibleSUSY`UseHiggs3LoopSplit === True,
              {threeLoopSelfEnergyPrototypes, threeLoopSelfEnergyFunctions} = SelfEnergies`CreateThreeLoopSelfEnergiesSplit[{SARAH`HiggsBoson}];
              threeLoopHiggsHeaders = "#include \"splitmssm_threeloophiggs.hpp\"\n";
             ];
           If[SARAH`UseHiggs2LoopMSSM === True,
              {twoLoopTadpolePrototypes, twoLoopTadpoleFunctions} = SelfEnergies`CreateTwoLoopTadpolesMSSM[SARAH`HiggsBoson];
              {twoLoopSelfEnergyPrototypes, twoLoopSelfEnergyFunctions} = SelfEnergies`CreateTwoLoopSelfEnergiesMSSM[{SARAH`HiggsBoson, SARAH`PseudoScalar}];
              twoLoopHiggsHeaders = "#include \"sfermions.hpp\"\n#include \"mssm_twoloophiggs.hpp\"\n";
             ];
           If[FlexibleSUSY`UseHiggs3LoopMSSM === True,
              {threeLoopSelfEnergyPrototypes, threeLoopSelfEnergyFunctions} = SelfEnergies`CreateThreeLoopSelfEnergiesMSSM[{SARAH`HiggsBoson}];
              threeLoopHiggsHeaders = threeLoopHiggsHeaders <> "\
#ifdef ENABLE_HIMALAYA
#include \"himalaya/HierarchyCalculator.hpp\"
#include \"himalaya/version.hpp\"
#endif
";
             ];
           If[FlexibleSUSY`UseHiggs3LoopNMSSM === True,
              {threeLoopSelfEnergyPrototypes, threeLoopSelfEnergyFunctions} = SelfEnergies`CreateThreeLoopSelfEnergiesNMSSM[{SARAH`HiggsBoson}];
              threeLoopHiggsHeaders = threeLoopHiggsHeaders <> "\
#ifdef ENABLE_HIMALAYA
#include \"himalaya/HierarchyCalculator.hpp\"
#include \"himalaya/version.hpp\"
#endif
";
             ];
           If[FlexibleSUSY`UseHiggs2LoopNMSSM === True,
              {twoLoopTadpolePrototypes, twoLoopTadpoleFunctions} = SelfEnergies`CreateTwoLoopTadpolesNMSSM[SARAH`HiggsBoson];
              {twoLoopSelfEnergyPrototypes, twoLoopSelfEnergyFunctions} = SelfEnergies`CreateTwoLoopSelfEnergiesNMSSM[{SARAH`HiggsBoson, SARAH`PseudoScalar}];
              twoLoopHiggsHeaders = "#include \"sfermions.hpp\"\n#include \"mssm_twoloophiggs.hpp\"\n#include \"nmssm_twoloophiggs.hpp\"\n";
             ];
           twoLoopThresholdHeaders = ThresholdCorrections`GetTwoLoopThresholdHeaders[];
           If[SARAH`UseHiggs2LoopMSSM === True ||
              FlexibleSUSY`UseHiggs2LoopNMSSM === True ||
              FlexibleSUSY`UseMSSMYukawa2Loop === True ||
              FlexibleSUSY`UseMSSMAlphaS2Loop === True ||
              FlexibleSUSY`UseHiggs3LoopMSSM === True ||
              FlexibleSUSY`UseHiggs3LoopNMSSM === True,
              {secondGenerationHelperPrototypes, secondGenerationHelperFunctions} = TreeMasses`CreateGenerationHelpers[2];
              {thirdGenerationHelperPrototypes, thirdGenerationHelperFunctions} = TreeMasses`CreateGenerationHelpers[3];
             ];
           {selfEnergyPrototypes, selfEnergyFunctions} = SelfEnergies`CreateNPointFunctions[nPointFunctions, vertexRules];
           phasesDefinition             = Phases`CreatePhasesDefinition[phases];
           phasesGetterSetters          = Phases`CreatePhasesGetterSetters[phases];
           If[Parameters`GetExtraParameters[] =!= {},
              extraParameterDefs           = StringJoin[Parameters`CreateParameterDefinitionAndDefaultInitialize
                                                        /@ Parameters`GetExtraParameters[]];
              extraParameterGetters        = StringJoin[CConversion`CreateInlineGetters[CConversion`ToValidCSymbolString[#],
                                                                                        CConversion`ToValidCSymbolString[#],
                                                                                        Parameters`GetType[#]]& /@
                                                        Parameters`GetExtraParameters[]];
              extraParameterSetters        = StringJoin[CConversion`CreateInlineSetters[CConversion`ToValidCSymbolString[#],
                                                                                        Parameters`GetType[#]]& /@
                                                        Parameters`GetExtraParameters[]];
              clearExtraParameters         = StringJoin[CConversion`SetToDefault[CConversion`ToValidCSymbolString[#],
                                                                                 Parameters`GetType[#]]& /@
                                                        Parameters`GetExtraParameters[]];
             ];
           loopMassesPrototypes         = LoopMasses`CreateOneLoopPoleMassPrototypes[];
           (* If you want to add tadpoles, call the following routine like this:
              CreateOneLoopPoleMassFunctions[diagonalizationPrecision, Cases[nPointFunctions, SelfEnergies`Tadpole[___]], vevs];
              *)
           loopMassesFunctions          = LoopMasses`CreateOneLoopPoleMassFunctions[diagonalizationPrecision, {}, {}];
           runningDRbarMassesPrototypes = LoopMasses`CreateRunningDRbarMassPrototypes[];
           runningDRbarMassesFunctions  = LoopMasses`CreateRunningDRbarMassFunctions[FlexibleSUSY`FSRenormalizationScheme];
           enablePoleMassThreads = False;
           callAllLoopMassFunctions     = LoopMasses`CallAllPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           enablePoleMassThreads = True;
           callAllLoopMassFunctionsInThreads = LoopMasses`CallAllPoleMassFunctions[FlexibleSUSY`FSEigenstates, enablePoleMassThreads];
           masses                       = Flatten[(FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices) /.
                                                  FlexibleSUSY`M[p_List] :> (FlexibleSUSY`M /@ p)];
           {lspGetters, lspFunctions}   = LoopMasses`CreateLSPFunctions[FlexibleSUSY`PotentialLSPParticles];
           printMasses                  = WriteOut`PrintParameters[masses, "ostr"];
           getMixings                   = TreeMasses`CreateMixingArrayGetter[massMatrices];
           setMixings                   = TreeMasses`CreateMixingArraySetter[massMatrices, "pars"];
           getMasses                    = TreeMasses`CreateMassArrayGetter[massMatrices];
           setMasses                    = TreeMasses`CreateMassArraySetter[massMatrices, "pars"];
           getExtraParameters           = Parameters`CreateExtraParameterArrayGetter[Parameters`GetExtraParameters[]];
           setExtraParameters           = Parameters`CreateExtraParameterArraySetter[Parameters`GetExtraParameters[]];
           mixingMatrices               = Flatten[TreeMasses`GetMixingMatrixSymbol[#]& /@ massMatrices];
           printMixingMatrices          = WriteOut`PrintParameters[mixingMatrices, "ostr"];
           dependencePrototypes         = TreeMasses`CreateDependencePrototypes[];
           dependenceFunctions          = TreeMasses`CreateDependenceFunctions[];
           softScalarMasses =
               If[SARAH`SupersymmetricModel,
                  DeleteDuplicates[SARAH`ListSoftBreakingScalarMasses],
                  Select[Parameters`GetModelParametersWithMassDimension[2], Parameters`IsRealParameter]
                 ];
           (* find soft Higgs masses that appear in tree-level EWSB eqs. *)
           treeLevelEWSBOutputParameters =
               Parameters`DecreaseIndexLiterals @
               Parameters`ExpandExpressions @
               Parameters`AppendGenerationIndices @
               If[MatchQ[FlexibleSUSY`FSSolveEWSBTreeLevelFor, {__}],
                  FlexibleSUSY`FSSolveEWSBTreeLevelFor,
                  (* each softScalarMasses should appear only in exactly 1 EWSB eq.! *)
                  Select[softScalarMasses, ParameterAppearsExactlyOnceIn[ewsbEquations, #]&]
                 ];
           If[MatchQ[treeLevelEWSBOutputParameters, {__}],
              parametersToSave = treeLevelEWSBOutputParameters;
              solveTreeLevelEWSBviaSoftHiggsMasses = First @ EWSB`FindSolutionAndFreePhases[independentEwsbEquationsTreeLevel,
                                                                                            treeLevelEWSBOutputParameters];
              If[solveTreeLevelEWSBviaSoftHiggsMasses === {},
                 Print["Error: could not find an analytic solution to the tree-level EWSB eqs."];
                 Print["   for the parameters ", treeLevelEWSBOutputParameters];
                 Quit[1];
                ];
              solveTreeLevelEWSBviaSoftHiggsMasses = EWSB`CreateMemberTreeLevelEwsbSolver[solveTreeLevelEWSBviaSoftHiggsMasses];
              solveEWSBTemporarily = IndentText["solve_ewsb_tree_level_custom();"];
              ,
              parametersToSave = parametersFixedByEWSB;
              solveTreeLevelEWSBviaSoftHiggsMasses = "";
              solveEWSBTemporarily = EWSB`SolveEWSBIgnoringFailures[0];
             ];
           parametersToSave = Join[parametersToSave,
                                   #[[1]]& /@ (Select[ewsbSubstitutions,
                                                      Function[sub, Or @@ (!FreeQ[sub[[2]], #]& /@ parametersToSave)]])];
           parametersToSave = DeleteDuplicates[parametersToSave /.
                                               {
                                                   Susyno`LieGroups`conj -> Identity,
                                                   SARAH`Conj            -> Identity,
                                                   SARAH`Tp              -> Identity,
                                                   SARAH`Adj             -> Identity,
                                                   SARAH`bar             -> Identity,
                                                   Conjugate             -> Identity,
                                                   Transpose             -> Identity,
                                                   Re                    -> Identity,
                                                   Im                    -> Identity
                                               }];
           saveEWSBOutputParameters = Parameters`SaveParameterLocally[parametersToSave];
           (ewsbSolverHeaders = ewsbSolverHeaders
                                <> EnableForBVPSolver[#, ("#include \"" <> FlexibleSUSY`FSModelName
                                                          <> "_" <> GetBVPSolverHeaderName[#] <> "_ewsb_solver.hpp\"\n")] <> "\n")&
                                /@ FlexibleSUSY`FSBVPSolvers;
           defaultEWSBSolverCctor = CreateDefaultEWSBSolverConstructor[FlexibleSUSY`FSBVPSolvers];
           reorderDRbarMasses           = TreeMasses`ReorderGoldstoneBosons[""];
           reorderPoleMasses            = TreeMasses`ReorderGoldstoneBosons["PHYSICAL"];
           checkPoleMassesForTachyons   = TreeMasses`CheckPoleMassesForTachyons["PHYSICAL"];
           WriteOut`ReplaceInFiles[files,
                          { "@[abstract]lspGetters@"           -> IndentText[FunctionModifiers`MakeAbstract[lspGetters]],
                            "@[override]lspGetters@"           -> IndentText[FunctionModifiers`MakeOverride[lspGetters]],
                            "@lspFunctions@"         -> lspFunctions,
                            "@[abstract]parameterGetters@" -> IndentText[FunctionModifiers`MakeAbstract[StringJoin[Parameters`CreateModelParameterGetter /@ Parameters`GetModelParameters[]]]],
                            "@[abstract]parameterSetters@" -> IndentText[FunctionModifiers`MakeAbstract[StringJoin[Parameters`CreateModelParameterSetter /@ Parameters`GetModelParameters[]]]],
                            "@[override]delegateParameterGetters@" -> IndentText[FunctionModifiers`MakeOverride[StringJoin[Parameters`CreateDelegateModelParameterGetter /@ Parameters`GetModelParameters[]]]],
                            "@[override]delegateParameterSetters@" -> IndentText[FunctionModifiers`MakeOverride[StringJoin[Parameters`CreateModelParameterSetter /@ Parameters`GetModelParameters[]]]],
                            "@massGetters@"          -> IndentText[massGetters],
                            "@[abstract]massGetters@" -> IndentText[FunctionModifiers`MakeAbstract[massGetters]],
                            "@[override]massGetters@" -> IndentText[FunctionModifiers`MakeOverride[massGetters]],
                            "@mixingMatrixGetters@"  -> IndentText[mixingMatrixGetters],
                            "@[abstract]mixingMatrixGetters@" -> IndentText[FunctionModifiers`MakeAbstract[mixingMatrixGetters]],
                            "@[override]mixingMatrixGetters@" -> IndentText[FunctionModifiers`MakeOverride[mixingMatrixGetters]],
                            "@slhaPoleMassGetters@"  -> IndentText[slhaPoleMassGetters],
                            "@slhaPoleMixingMatrixGetters@" -> IndentText[slhaPoleMixingMatrixGetters],
                            "@higgsMassGetterPrototypes@"   -> IndentText[higgsMassGetters[[1]]],
                            "@[abstract]higgsMassGetterPrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[higgsMassGetters[[1]]]],
                            "@[override]higgsMassGetterPrototypes@" -> IndentText[FunctionModifiers`MakeOverride[higgsMassGetters[[1]]]],
                            "@higgsMassGetters@"     -> higgsMassGetters[[2]],
                            "@tadpoleEqPrototypes@"  -> IndentText[tadpoleEqPrototypes],
                            "@[abstract]tadpoleEqPrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[tadpoleEqPrototypes]],
                            "@[override]tadpoleEqPrototypes@" -> IndentText[FunctionModifiers`MakeOverride[tadpoleEqPrototypes]],
                            "@tadpoleEqFunctions@"   -> tadpoleEqFunctions,
                            "@numberOfEWSBEquations@"-> ToString[TreeMasses`GetDimension[SARAH`HiggsBoson]],
                            "@calculateTreeLevelTadpoles@" -> IndentText[calculateTreeLevelTadpoles],
                            "@calculateOneLoopTadpoles@"   -> IndentText @ IndentText[calculateOneLoopTadpoles],
                            "@calculateTwoLoopTadpoles@"   -> IndentText @ IndentText @ IndentText[calculateTwoLoopTadpoles],
                            "@divideTadpoleByVEV@"     -> IndentText[divideTadpoleByVEV],
                            "@clearOutputParameters@"  -> IndentText[clearOutputParameters],
                            "@clearPhases@"            -> IndentText[clearPhases],
                            "@copyDRbarMassesToPoleMasses@" -> IndentText[copyDRbarMassesToPoleMasses],
                            "@calculateDecouplingGaugeCouplings@" -> IndentText @ IndentText[ThresholdCorrections`CalculateGaugeCouplings[]],
                            "@setDecouplingGaugeCouplings@" -> IndentText @ IndentText[
                                Constraint`ApplyConstraints[
                                    Cases[FlexibleSUSY`LowScaleInput,
                                          {SARAH`hyperchargeCoupling, __} |
                                          {SARAH`leftCoupling, __} |
                                          {SARAH`strongCoupling, __} |
                                          (* any gauge coupling that is fixed in terms of SM parameters *)
                                          {_?Parameters`IsGaugeCoupling | (_?Parameters`IsGaugeCoupling)[__],
                                           _?(!FreeQ[smDecouplingParameters,#])&}
                                    ]
                                ]
                                                                            ],
                            "@setDecouplingVEV@" -> IndentText @ IndentText @ WrapLines[
                                Constraint`ApplyConstraints[
                                    Cases[FlexibleSUSY`LowScaleInput,
                                          {SARAH`VEVSM, __} |
                                          {SARAH`VEVSM1, __} |
                                          {SARAH`VEVSM2, __} |
                                          (* any VEV that is fixed in terms of SM parameters *)
                                          {_?Parameters`IsVEV | (_?Parameters`IsVEV)[__],
                                           _?(!FreeQ[smDecouplingParameters,#])&}
                                    ]
                                ]
                                                                              ],
                            "@setDecouplingYukawaUpQuarks@"  -> IndentText @ IndentText[
                                ThresholdCorrections`SetDRbarYukawaCouplingTop[FlexibleSUSY`LowScaleInput]
                                                                             ],
                            "@setDecouplingYukawaDownQuarks@"  -> IndentText @ IndentText[
                                ThresholdCorrections`SetDRbarYukawaCouplingBottom[FlexibleSUSY`LowScaleInput]
                                                                             ],
                            "@setDecouplingYukawaDownLeptons@"  -> IndentText @ IndentText[
                                ThresholdCorrections`SetDRbarYukawaCouplingElectron[FlexibleSUSY`LowScaleInput]
                                                                             ],
                            "@overrideTreeHiggsMixings@" -> IndentText@IndentText@StringRiffle[
                               With[{mixingMatrix = FindMixingMatrixSymbolFor[#]},
                                  If[mixingMatrix =!= Null,
                                     ToString@mixingMatrix <> " = this->get_physical()." <> ToString@mixingMatrix <> ";",
                                     ""
                                  ]
                               ]& /@ DeleteCases[{TreeMasses`GetHiggsBoson[], TreeMasses`GetPseudoscalarHiggsBoson[], TreeMasses`GetChargedHiggsBoson[]}, Null],
                               "\n"
                            ],
                            "@copyRunningBSMMassesToDecouplingMasses@" -> IndentText[copyRunningBSMMassesToDecouplingMasses],
                            "@reorderDRbarMasses@"     -> IndentText[reorderDRbarMasses],
                            "@reorderPoleMasses@"      -> IndentText[reorderPoleMasses],
                            "@checkPoleMassesForTachyons@" -> IndentText[checkPoleMassesForTachyons],
                            "@physicalMassesDef@"      -> IndentText[physicalMassesDef],
                            "@mixingMatricesDef@"      -> IndentText[mixingMatricesDef],
                            "@massCalculationPrototypes@" -> IndentText[massCalculationPrototypes],
                            "@[abstract]massCalculationPrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[massCalculationPrototypes]],
                            "@[override]massCalculationPrototypes@" -> IndentText[FunctionModifiers`MakeOverride[massCalculationPrototypes]],
                            "@massCalculationFunctions@"  -> WrapLines[massCalculationFunctions],
                            "@calculateAllMasses@"        -> IndentText[calculateAllMasses],
                            "@selfEnergyPrototypes@"      -> IndentText[selfEnergyPrototypes],
                            "@[abstract]selfEnergyPrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[selfEnergyPrototypes]],
                            "@[override]selfEnergyPrototypes@" -> IndentText[FunctionModifiers`MakeOverride[selfEnergyPrototypes]],
                            "@selfEnergyFunctions@"       -> selfEnergyFunctions,
                            "@selfEnergyDerivPrototypes@" -> selfEnergyDerivPrototypes,
                            "@[abstract]selfEnergyDerivPrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[selfEnergyDerivPrototypes]],
                            "@[override]selfEnergyDerivPrototypes@" -> IndentText[FunctionModifiers`MakeOverride[selfEnergyDerivPrototypes]],
                            "@selfEnergyDerivFunctions@"  -> selfEnergyDerivFunctions,
                            "@twoLoopTadpolePrototypes@"  -> IndentText[twoLoopTadpolePrototypes],
                            "@twoLoopTadpoleFunctions@"   -> twoLoopTadpoleFunctions,
                            "@twoLoopSelfEnergyPrototypes@" -> IndentText[twoLoopSelfEnergyPrototypes],
                            "@twoLoopSelfEnergyFunctions@"  -> twoLoopSelfEnergyFunctions,
                            "@twoLoopHiggsHeaders@"       -> twoLoopHiggsHeaders,
                            "@twoLoopThresholdHeaders@"   -> twoLoopThresholdHeaders,
                            "@threeLoopSelfEnergyPrototypes@" -> IndentText[threeLoopSelfEnergyPrototypes],
                            "@threeLoopSelfEnergyFunctions@"  -> threeLoopSelfEnergyFunctions,
                            "@threeLoopHiggsHeaders@"         -> threeLoopHiggsHeaders,
                            "@fourLoopSelfEnergyPrototypes@"  -> IndentText[fourLoopSelfEnergyPrototypes],
                            "@fourLoopSelfEnergyFunctions@"   -> fourLoopSelfEnergyFunctions,
                            "@fourLoopHiggsHeaders@"          -> fourLoopHiggsHeaders,
                            "@secondGenerationHelperPrototypes@"-> IndentText[secondGenerationHelperPrototypes],
                            "@secondGenerationHelperFunctions@" -> secondGenerationHelperFunctions,
                            "@thirdGenerationHelperPrototypes@" -> IndentText[thirdGenerationHelperPrototypes],
                            "@thirdGenerationHelperFunctions@"  -> thirdGenerationHelperFunctions,
                            "@phasesDefinition@"          -> IndentText[phasesDefinition],
                            "@phasesGetterSetters@"          -> IndentText[phasesGetterSetters],
                            "@[abstract]phasesGetterSetters@"-> IndentText[FunctionModifiers`MakeAbstract[phasesGetterSetters]],
                            "@[override]phasesGetterSetters@"-> IndentText[FunctionModifiers`MakeOverride[phasesGetterSetters]],
                            "@extraParameterDefs@"           -> IndentText[extraParameterDefs],
                            "@extraParameterGetters@"        -> IndentText[extraParameterGetters],
                            "@extraParameterSetters@"        -> IndentText[extraParameterSetters],
                            "@[abstract]extraParameterGetters@" -> IndentText[FunctionModifiers`MakeAbstract[extraParameterGetters]],
                            "@[abstract]extraParameterSetters@" -> IndentText[FunctionModifiers`MakeAbstract[extraParameterSetters]],
                            "@[override]extraParameterGetters@" -> IndentText[FunctionModifiers`MakeOverride[extraParameterGetters]],
                            "@[override]extraParameterSetters@" -> IndentText[FunctionModifiers`MakeOverride[extraParameterSetters]],
                            "@clearExtraParameters@"         -> IndentText[clearExtraParameters],
                            "@loopMassesPrototypes@"         -> IndentText[WrapLines[loopMassesPrototypes]],
                            "@[abstract]loopMassesPrototypes@"-> IndentText[WrapLines[FunctionModifiers`MakeAbstract[loopMassesPrototypes]]],
                            "@[override]loopMassesPrototypes@"-> IndentText[WrapLines[FunctionModifiers`MakeOverride[loopMassesPrototypes]]],
                            "@loopMassesFunctions@"          -> WrapLines[loopMassesFunctions],
                            "@runningDRbarMassesPrototypes@" -> IndentText[runningDRbarMassesPrototypes],
                            "@[abstract]runningDRbarMassesPrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[runningDRbarMassesPrototypes]],
                            "@[override]runningDRbarMassesPrototypes@" -> IndentText[FunctionModifiers`MakeOverride[runningDRbarMassesPrototypes]],
                            "@runningDRbarMassesFunctions@"  -> WrapLines[runningDRbarMassesFunctions],
                            "@callAllLoopMassFunctions@"     -> IndentText[callAllLoopMassFunctions],
                            "@callAllLoopMassFunctionsInThreads@" -> IndentText[callAllLoopMassFunctionsInThreads],
                            "@printMasses@"                  -> IndentText[printMasses],
                            "@getMixings@"                   -> IndentText[getMixings],
                            "@setMixings@"                   -> IndentText[setMixings],
                            "@getMasses@"                    -> IndentText[getMasses],
                            "@setMasses@"                    -> IndentText[setMasses],
                            "@getExtraParameters@"           -> IndentText[getExtraParameters],
                            "@setExtraParameters@"           -> IndentText[setExtraParameters],
                            "@printMixingMatrices@"          -> IndentText[printMixingMatrices],
                            "@dependencePrototypes@"         -> IndentText[dependencePrototypes],
                            "@[abstract]dependencePrototypes@" -> IndentText[FunctionModifiers`MakeAbstract[dependencePrototypes]],
                            "@[override]dependencePrototypes@" -> IndentText[FunctionModifiers`MakeOverride[dependencePrototypes]],
                            "@dependenceFunctions@"          -> WrapLines[dependenceFunctions],
                            "@saveEWSBOutputParameters@"     -> IndentText[saveEWSBOutputParameters],
                            "@solveTreeLevelEWSBviaSoftHiggsMasses@" -> IndentText[WrapLines[solveTreeLevelEWSBviaSoftHiggsMasses]],
                            "@solveEWSBTemporarily@"         -> solveEWSBTemporarily,
                            "@convertMixingsToSLHAConvention@" -> IndentText[convertMixingsToSLHAConvention],
                            "@convertMixingsToHKConvention@"   -> IndentText[convertMixingsToHKConvention],
                            "@ewsbSolverHeaders@"            -> ewsbSolverHeaders,
                            "@defaultEWSBSolverCctor@"       -> defaultEWSBSolverCctor,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteDecaysClass[decayParticles_List, finalStateParticles_List, files_List] :=
    Module[{maxFinalStateParticles = 2, decaysLists = {}, decaysVertices, numberOfDecayParticles = 0,
            enableDecaysCalculationThreads,
            callAllDecaysFunctions = "", callAllDecaysFunctionsInThreads = "",
            decaysListGettersPrototypes = "", decaysListGettersFunctions = "",
            decaysGetters = "", initDecayTable = "",
            decaysCalculationPrototypes = "", decaysCalculationFunctions = "",
            partialWidthCalculationPrototypes = "", partialWidthCalculationFunctions = "",
            calcAmplitudeSpecializationDecls = "", calcAmplitudeSpecializationDefs = "",
            partialWidthSpecializationDecls = "", partialWidthSpecializationDefs = "",
            solverIncludes = "", solver = "", contentOfPath = $Path, modelName, bsmParticleAliasList},

            modelName =
               If[SARAH`submodeldir =!= False,
                  SARAH`modelDir <> "/" <> SARAH`submodeldir,
                  SARAH`modelDir
               ];

           (solverIncludes = solverIncludes <> EnableSpectrumGenerator[#])& /@ FlexibleSUSY`FSBVPSolvers;

           numberOfDecayParticles = Plus @@ (TreeMasses`GetDimensionWithoutGoldstones /@ decayParticles);

           (* in parallel model decay of every particle is computed in a separate theread
              but if we have only one particle then turning parallelism on only introduces
              an overhead with no speed up *)
           Block[{FlexibleSUSY`FSEnableParallelism = If[Length[decayParticles] < 2, False, FlexibleSUSY`FSEnableParallelism]},
              (* create list containing elements {field, {FSParticleDecay 'objects'}} *)
              If[FlexibleSUSY`FSEnableParallelism,
                 If[Head[SARAH`VertexList3] === Symbol || Length[SARAH`VertexList3] === 0,
                    SA`CurrentStates = FlexibleSUSY`FSEigenstates;
                    SARAH`InitVertexCalculation[FlexibleSUSY`FSEigenstates, False];
                    SARAH`partDefinition = ParticleDefinitions[FlexibleSUSY`FSEigenstates];
                    SARAH`Particles[SARAH`Current] = SARAH`Particles[FlexibleSUSY`FSEigenstates];
                    SARAH`ReadVertexList[FlexibleSUSY`FSEigenstates, False, False, True];
                    SARAH`MakeCouplingLists;
                 ];
                 LaunchKernels[];
                 DistributeDefinitions[contentOfPath, modelName];
                 ParallelEvaluate[
                    (* subkernels have different $Path variable than the main kernel
                       https://mathematica.stackexchange.com/questions/11595/package-found-with-needs-but-not-with-parallelneeds *)
                    $Path = contentOfPath;
                    (* don't pollute terminal with SARAH initialization message *)
                    Block[{Print,workingDirectory = Directory[]},
                       << SARAH`;
                       SARAH`SARAH[OutputDirectory] = FileNameJoin[{workingDirectory, "Output"}];
                       SARAH`SARAH[InputDirectories] = {
                          FileNameJoin[{workingDirectory, "sarah"}],
                          ToFileName[{$sarahDir, "Models"}]
                       };
                       (* SARAH looks for models in dirs in SARAH`SARAH[InputDirectories] list.
                          If the model doesn't exist in any particular location mathematica prints a
                          FileByteCount::fdnfnd warning. This is harmless because at this point in
                          the code we know that the model does exist at least somewhere. *)
                       Off[FileByteCount::fdnfnd];
                       Off[Superpotential::ViolationGlobal];
                       Start@modelName;
                       On[FileByteCount::fdnfnd];
                       On[Superpotential::ViolationGlobal];
                    ];,
                    DistributedContexts -> None
                 ];
                 DistributeDefinitions[SARAH`VertexList3, SARAH`VertexList4];
                 decaysLists =
                    AbsoluteTiming @ ParallelMap[
                       {#, Decays`GetDecaysForParticle[#, maxFinalStateParticles, finalStateParticles]}&,
                       decayParticles,
                       DistributedContexts -> All,
                       Method -> "FinestGrained"
                    ];
                 Needs["Parallel`Developer`"];
                 Parallel`Developer`ClearDistributedDefinitions[];
                 Parallel`Developer`ClearKernels[];
                 CloseKernels[];
                 ,
                 decaysLists =
                    AbsoluteTiming @ Map[
                       {#, Decays`GetDecaysForParticle[#, maxFinalStateParticles, finalStateParticles]}&,
                       decayParticles
                    ];
              ]
           ];
           Print[""];
           Print["Creation of decay amplitudes took", FSRound[First@decaysLists, 1], "s"];
           decaysLists = Last@decaysLists;

           (* get from generated FSParticleDecay 'objects' vertices needed in decay calculation *)
           decaysVertices = DeleteDuplicates[Flatten[Decays`GetVerticesForDecays[Last[#]]& /@ decaysLists, 1]];
           decaysVertices = SortFieldsInCp /@ decaysVertices;
           With[{Wm = If[GetElectricCharge[TreeMasses`GetWBoson[]] < 0, Susyno`LieGroups`conj[TreeMasses`GetWBoson[]], Susyno`LieGroups`conj[TreeMasses`GetWBoson[]]]},
              AppendTo[decaysVertices, {TreeMasses`GetHiggsBoson[], Susyno`LieGroups`conj[Wm], Wm}]
           ];

           enableDecaysCalculationThreads = False;
           callAllDecaysFunctions = Decays`CallDecaysCalculationFunctions[decayParticles, enableDecaysCalculationThreads];
           enableDecaysCalculationThreads = True;
           callAllDecaysFunctionsInThreads = Decays`CallDecaysCalculationFunctions[decayParticles, enableDecaysCalculationThreads];
           decaysCalculationPrototypes = Decays`CreateDecaysCalculationPrototypes[decaysLists];
           decaysCalculationFunctions = Decays`CreateDecaysCalculationFunctions[decaysLists];
           partialWidthCalculationPrototypes = Decays`CreatePartialWidthCalculationPrototypes[decaysLists];
           partialWidthCalculationFunctions = Decays`CreatePartialWidthCalculationFunctions[decaysLists, FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"];
           decaysGetters = Decays`CreateDecaysGetterFunctions[decayParticles];
           decaysListGettersPrototypes = Decays`CreateDecayTableGetterPrototypes[decayParticles];
           decaysListGettersFunctions = Decays`CreateDecayTableGetterFunctions[decayParticles, FlexibleSUSY`FSModelName <> "_decay_table"];
           initDecayTable = Decays`CreateDecayTableInitialization[decayParticles];

           {calcAmplitudeSpecializationDecls, calcAmplitudeSpecializationDefs}
               = Decays`CreateTotalAmplitudeSpecializations[decaysLists, FlexibleSUSY`FSModelName];
           With[{temp = Decays`CreatePartialWidthSpecializations[decaysLists, FlexibleSUSY`FSModelName]},
              If[temp =!= {},
                 {partialWidthSpecializationDecls, partialWidthSpecializationDefs}
                    = temp;
              ]
           ];

           bsmParticleAliasList = Decays`CreateBSMParticleAliasList["fields"];

           solver = GetBVPSolverTemplateParameter[First@FlexibleSUSY`FSBVPSolvers];

           WriteOut`ReplaceInFiles[files,
                          { "@callAllDecaysFunctions@" -> IndentText[callAllDecaysFunctions],
                            "@callAllDecaysFunctionsInThreads@" -> IndentText[callAllDecaysFunctionsInThreads],
                            "@decaysGetters@" -> IndentText[WrapLines[decaysGetters]],
                            "@decaysCalculationPrototypes@" -> IndentText[decaysCalculationPrototypes],
                            "@decaysCalculationFunctions@" -> WrapLines[decaysCalculationFunctions],
                            "@partialWidthCalculationPrototypes@" -> TextFormatting`IndentText[partialWidthCalculationPrototypes],
                            "@partialWidthCalculationFunctions@" -> partialWidthCalculationFunctions,
                            "@calcAmplitudeSpecializationDecls@" -> calcAmplitudeSpecializationDecls,
                            "@calcAmplitudeSpecializationDefs@" -> calcAmplitudeSpecializationDefs,
                            "@partialWidthSpecializationDecls@" -> partialWidthSpecializationDecls,
                            "@partialWidthSpecializationDefs@" -> partialWidthSpecializationDefs,
                            "@decaysListGettersPrototypes@" -> IndentText[decaysListGettersPrototypes],
                            "@decaysListGettersFunctions@" -> decaysListGettersFunctions,
                            "@initDecayTable@" -> IndentText[WrapLines[initDecayTable]],
                            "@numberOfDecayParticles@" -> ToString[numberOfDecayParticles],
                            "@create_BSM_particle_list@" -> Last@bsmParticleAliasList,
                            "@gs_name@" -> ToString[TreeMasses`GetStrongCoupling[]],
                            "@solver@" -> solver,
                            "@solverIncludes@" -> solverIncludes,
                            "@isCPodd@" -> If[TreeMasses`GetPseudoscalarHiggsBoson[] =!= Null && GetDimensionWithoutGoldstones[TreeMasses`GetPseudoscalarHiggsBoson[]] > 0, " || std::is_same_v<FieldIn, " <> FSModelName <> "_cxx_diagrams::fields::PseudoscalarHiggs>", ""],
                            Sequence @@ GeneralReplacementRules[]
                          } ];

           DeleteDuplicates@Join[decaysVertices, First@bsmParticleAliasList]
          ];

WriteBVPSolverTemplates[files_List] :=
    WriteOut`ReplaceInFiles[files, { Sequence @@ GeneralReplacementRules[] }];

WriteTwoScaleModelClass[files_List] :=
    WriteOut`ReplaceInFiles[files, { Sequence @@ GeneralReplacementRules[] }];

WriteSemiAnalyticSolutionsClass[semiAnalyticBCs_List, semiAnalyticSolns_List, files_List] :=
    Module[{semiAnalyticSolutionsDefs = "", boundaryValueStructDefs = "",
            coefficientGetters = "",
            createBasisEvaluators = "", applyBoundaryConditions = "",
            datasets, numberOfTrialPoints, initializeTrialBoundaryValues = "",
            createLinearSystemSolvers = "", calculateCoefficients = "",
            evaluateSemiAnalyticSolns = "",
            calculateCoefficientsPrototypes = "", calculateCoefficientsFunctions = ""},
           semiAnalyticSolutionsDefs = SemiAnalytic`CreateSemiAnalyticSolutionsDefinitions[semiAnalyticSolns];
           boundaryValueStructDefs = SemiAnalytic`CreateLocalBoundaryValuesDefinitions[semiAnalyticSolns];
           coefficientGetters = SemiAnalytic`CreateSemiAnalyticCoefficientGetters[semiAnalyticSolns];
           createBasisEvaluators = SemiAnalytic`CreateBasisEvaluators[semiAnalyticSolns];
           datasets = SemiAnalytic`ConstructTrialDatasets[semiAnalyticSolns];
           createLinearSystemSolvers = SemiAnalytic`CreateLinearSystemSolvers[datasets, semiAnalyticSolns];
           {numberOfTrialPoints, initializeTrialBoundaryValues} = SemiAnalytic`InitializeTrialInputValues[datasets];
           applyBoundaryConditions = SemiAnalytic`ApplySemiAnalyticBoundaryConditions[semiAnalyticBCs, semiAnalyticSolns];
           evaluateSemiAnalyticSolns = SemiAnalytic`EvaluateSemiAnalyticSolutions[semiAnalyticSolns];
           calculateCoefficients = SemiAnalytic`CalculateCoefficients[datasets];
           {calculateCoefficientsPrototypes, calculateCoefficientsFunctions} = SemiAnalytic`CreateCoefficientsCalculations[semiAnalyticSolns];
           WriteOut`ReplaceInFiles[files, { "@semiAnalyticSolutionsDefs@" -> IndentText[WrapLines[semiAnalyticSolutionsDefs]],
                                            "@boundaryValuesDefs@" -> IndentText[WrapLines[boundaryValuesDefs]],
                                            "@boundaryValueStructDefs@" -> IndentText[IndentText[WrapLines[boundaryValueStructDefs]]],
                                            "@coefficientGetters@" -> IndentText[WrapLines[coefficientGetters]],
                                            "@numberOfTrialPoints@" -> ToString[numberOfTrialPoints],
                                            "@initializeTrialBoundaryValues@" -> IndentText[WrapLines[initializeTrialBoundaryValues]],
                                            "@createBasisEvaluators@" -> IndentText[WrapLines[createBasisEvaluators]],
                                            "@createLinearSystemSolvers@" -> IndentText[WrapLines[createLinearSystemSolvers]],
                                            "@calculateCoefficients@" -> IndentText[calculateCoefficients],
                                            "@applyBoundaryConditions@" -> IndentText[WrapLines[applyBoundaryConditions]],
                                            "@evaluateSemiAnalyticSolns@" -> IndentText[WrapLines[evaluateSemiAnalyticSolns]],
                                            "@calculateCoefficientsPrototypes@" -> IndentText[calculateCoefficientsPrototypes],
                                            "@calculateCoefficientsFunctions@" -> calculateCoefficientsFunctions,
                                            Sequence @@ GeneralReplacementRules[] }];
          ];

WriteSemiAnalyticModelClass[semiAnalyticBCs_List, semiAnalyticSolns_List, files_List] :=
    Module[{getSemiAnalyticCoefficients = "", printSemiAnalyticCoefficients = ""},
           getSemiAnalyticCoefficients = SemiAnalytic`GetModelCoefficients[semiAnalyticSolns];
           printSemiAnalyticCoefficients = SemiAnalytic`PrintModelCoefficients[semiAnalyticSolns, "out"];
           WriteOut`ReplaceInFiles[files, { "@getSemiAnalyticCoefficients@"  -> IndentText[WrapLines[getSemiAnalyticCoefficients]],
                                            "@printSemiAnalyticCoefficients@" -> IndentText[WrapLines[printSemiAnalyticCoefficients]],
                                            Sequence @@ GeneralReplacementRules[] }];
          ];

WriteTwoScaleSpectrumGeneratorClass[files_List] :=
    Module[{fillSMFermionPoleMasses = ""},
           fillSMFermionPoleMasses = FlexibleEFTHiggsMatching`FillSMFermionPoleMasses[];
           WriteOut`ReplaceInFiles[files,
                          { "@fillSMFermionPoleMasses@" -> IndentText[fillSMFermionPoleMasses],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteShootingSpectrumGeneratorClass[files_List] :=
    Module[{fillSMFermionPoleMasses = ""},
           fillSMFermionPoleMasses = FlexibleEFTHiggsMatching`FillSMFermionPoleMasses[];
           WriteOut`ReplaceInFiles[files,
                          {
                            "@fillSMFermionPoleMasses@" -> IndentText[fillSMFermionPoleMasses],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteSemiAnalyticSpectrumGeneratorClass[files_List] :=
    Module[{boundaryConstraint = "", semiAnalyticConstraint = "", getBoundaryScale = ""},
           Which[SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`HighScaleInput],
                 boundaryConstraint = "high_scale_constraint";
                 getBoundaryScale = "get_high_scale()";,
                 SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`SUSYScaleInput],
                 boundaryConstraint = "susy_scale_constraint";
                 getBoundaryScale = "get_susy_scale()";,
                 SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`LowScaleInput],
                 boundaryConstraint = "low_scale_constraint";
                 getBoundaryScale = "get_low_scale();",
                 True,
                 boundaryConstraint = "high_scale_constraint";
                 getBoundaryScale = "get_high_scale()";
                ];
           Which[SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`HighScaleInput],
                 semiAnalyticConstraint = "high_scale_constraint";,
                 SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`SUSYScaleInput],
                 semiAnalyticConstraint = "susy_scale_constraint";,
                 SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`LowScaleInput],
                 semiAnalyticConstraint = "low_scale_constraint";,
                 True,
                 semiAnalyticConstraint = "susy_scale_constraint";
                ];
           WriteOut`ReplaceInFiles[files,
                          { "@boundaryConstraint@" -> boundaryConstraint,
                            "@getBoundaryScale@" -> getBoundaryScale,
                            "@semiAnalyticConstraint@" -> semiAnalyticConstraint,
                            "@fillSMFermionPoleMasses@" -> IndentText[fillSMFermionPoleMasses],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

(* Write the observables files *)
WriteObservables[extraSLHAOutputBlocks_, files_List] :=
    Module[{requestedObservables, numberOfObservables, observablesDef,
            observablesInit, getObservables, getObservablesNames,
            clearObservables, setObservables, calculateObservables},
           requestedObservables = Observables`GetRequestedObservables[extraSLHAOutputBlocks];
           numberOfObservables = Observables`CountNumberOfObservables[requestedObservables];
           observablesDef = Observables`CreateObservablesDefinitions[requestedObservables];
           observablesInit = Observables`CreateObservablesInitialization[requestedObservables];
           {getObservables, getObservablesNames, setObservables} =
               Observables`CreateSetAndDisplayObservablesFunctions[requestedObservables];
           clearObservables = Observables`CreateClearObservablesFunction[requestedObservables];
           calculateObservables = Observables`CalculateObservables[requestedObservables, "observables"];


           WriteOut`ReplaceInFiles[files,
                                   {   "@numberOfObservables@" -> ToString[numberOfObservables],
                                       "@observablesDef@" -> IndentText[observablesDef],
                                       "@observablesInit@" -> IndentText[WrapLines[observablesInit]],
                                       "@getObservables@" -> IndentText[getObservables],
                                       "@getObservablesNames@" -> IndentText[getObservablesNames],
                                       "@clearObservables@" -> IndentText[clearObservables],
                                       "@setObservables@" -> IndentText[setObservables],
                                       "@calculateObservables@" -> IndentText @ IndentText[calculateObservables],
                                       "@observablesHeaders@" -> Observables`GetObservablesHeaders[],
                                       Sequence @@ GeneralReplacementRules[]
                                   } ];
           ];

(* Write the CXXDiagrams c++ files *)
WriteCXXDiagramClass[vertices_List, files_List,
    cxxQFTVerticesTemplate_, cxxQFTVerticesOutputDirectory_,
    cxxQFTVerticesMakefileTemplates_] :=
    Module[{fields = "", cxxVerticesParts = {},
            massFunctions, physicalMassFunctions,
            unitCharge,
            defineFieldTraits,
            sarahOutputDir = SARAH`$sarahCurrentOutputMainDir,
            outputDir, cxxDiagramsDir, createdVerticesFile, fileHandle,
            cxxQFTVerticesFiles, realFieldsconjtraits},

        massFunctions = CXXDiagrams`CreateMassFunctions[];
        physicalMassFunctions = CXXDiagrams`CreatePhysicalMassFunctions[];
        {fields, realFieldsconjtraits} = CXXDiagrams`CreateFields[PotentialLSPParticles];
        defineFieldTraits =
           CXXDiagrams`CreateFieldTraitsDefinitions[
              TreeMasses`GetParticles[], "flexiblesusy::" <> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"
           ];

        If[vertices =!= {},

            cxxVerticesParts = CXXDiagrams`CreateVertices[vertices];

            (* Document which vertices are created. This is mainly useful for
               unit testing. See e.g test/test_MSSM_npointfunctions.m *)
            outputDir = FileNameJoin[{sarahOutputDir, ToString[FlexibleSUSY`FSEigenstates]}];
            cxxDiagramsDir = FileNameJoin[{outputDir, "CXXDiagrams"}];
            createdVerticesFile = FileNameJoin[{cxxDiagramsDir, "CreatedVertices.m"}];

            If[DirectoryQ[cxxDiagramsDir] === False,
                CreateDirectory[cxxDiagramsDir]];

            (* There is a bug in WriteString[] in older Mathematica versions
               that causes the files to be left open. *)
            fileHandle = OpenWrite[createdVerticesFile];
            Write[fileHandle, vertices];
            Close[fileHandle];
        ];

        unitCharge = CXXDiagrams`CreateUnitCharge[];
        AppendTo[cxxVerticesParts, {"", unitCharge}];

        WriteOut`ReplaceInFiles[files,
                            {"@CXXDiagrams_Fields@"                -> fields,
                             "@CXXDiagrams_realFields@"            -> realFieldsconjtraits,
                             "@CXXDiagrams_MassFunctions@"         -> massFunctions,
                             "@CXXDiagrams_PhysicalMassFunctions@" -> physicalMassFunctions,
                             "@defineFieldTraits@"                 -> defineFieldTraits,
                             "@CXXDiagrams_VertexPrototypes@"  ->
                                StringRiffle[cxxVerticesParts[[All, 1]], "\n\n"],
                             Sequence @@ GeneralReplacementRules[]
                            }];

        cxxQFTVerticesFiles = Table[
            {cxxQFTVerticesTemplate,
		      FileNameJoin[{cxxQFTVerticesOutputDirectory,
			   FSModelName <> "_" <> FileNameTake[StringReplace[cxxQFTVerticesTemplate,
			     {".cpp.in" -> ToString[k] <> ".cpp"}]]}]
			   },
		      {k, Length[cxxVerticesParts]}];

        WriteOut`ReplaceInFiles[{#[[1]]},
            {"@CXXDiagrams_VertexDefinitions@" -> #[[2, 2]],
            Sequence @@ GeneralReplacementRules[]
            }] & /@ Transpose[{cxxQFTVerticesFiles, cxxVerticesParts}];
        WriteOut`ReplaceInFiles[cxxQFTVerticesMakefileTemplates,
            {"@generatedCXXVerticesFiles@" ->
                "\t" <> StringRiffle[cxxQFTVerticesFiles[[All, 2]], " \\\n\t"],
            Sequence @@ GeneralReplacementRules[]
            }];
    ];

(* Write the EDM c++ files *)
WriteEDMClass[fields_List, files_List] :=
  Module[{(* in models without flavour violation (no FV models) lepton does not have an index *)
          leptonIndex = If[Length[fields] > 0, If[TreeMasses`GetDimension[First@fields] =!= 1, "idx", ""], ""],
          (* we want to calculate an offset of g-2 compared to the SM *)
          discardSMcontributions = CConversion`CreateCBoolValue[False],
          calculation, calculateForwadDeclaration},

    calculation =
       If[Length[fields] =!= 0,
            "const auto form_factors = " <>
            FSModelName <> "_FFV_form_factors::calculate_form_factors<Lepton,Lepton," <>
            CXXDiagrams`CXXNameOfField[TreeMasses`GetPhoton[]] <> ">(" <>
            leptonIndex <> If[leptonIndex === "", "", ", "] <>
            leptonIndex <> If[leptonIndex === "", "", ", "] <>
            "model, " <> discardSMcontributions <> ");",
            "const std::valarray<std::complex<double>> form_factors {0., 0., 0., 0.};"
       ];

    calculateForwadDeclaration = StringRiffle[EDM`EDMForwardDeclaration[#, "calculate_edm"]& /@ fields, "\n"];

    WriteOut`ReplaceInFiles[files,
                            {"@EDMCalculation@"       -> TextFormatting`IndentText[calculation],
                             "@extraIdxDecl@" -> If[GetParticleFromDescription["Leptons"] =!= Null, ", int idx", ""],
                             "@calculateForwadDeclaration@" -> calculateForwadDeclaration,
                             "@extraIdxUsageNoComma@" -> If[GetParticleFromDescription["Leptons"] =!= Null, "idx", ""],
                             Sequence @@ GeneralReplacementRules[]
                            }];

    vertices
  ];

(* Write the FFV c++ files *)
WriteFFVFormFactorsClass[extParticles_List, files_List] :=
   Module[{
         interfacePrototypes = "", interfaceDefinitions = "", templateWrapperDecl = "", templateWrapperDef = "",
         graphs, diagrams, vertices = {}
      },

      If[extParticles =!= {},

      	graphs = FFVFormFactors`FFVGraphs[];
	      diagrams =
            Outer[FFVFormFactors`FFVContributingDiagramsForGraph, graphs, extParticles, 1];

         (* group things not according to graphs but according to external states *)
         (* diagrams[[i,j]] will be: i is for a given external state, j for a topology *)
         diagrams = Transpose @ diagrams;

      	vertices = Flatten[CXXDiagrams`VerticesForDiagram /@ Flatten[diagrams,2], 1];

         {interfacePrototypes, interfaceDefinitions, templateWrapperDecl, templateWrapperDef} =
            StringJoin /@ (Riffle[#, "\n\n"]& /@ Transpose[
               FFVFormFactors`FFVFormFactorsCreateInterfaceFunction[#1, graphs, #2]& @@@
                  Transpose[{extParticles, diagrams}]
            ]);
         ];

      WriteOut`ReplaceInFiles[files,
         {"@FFVFormFactors_InterfacePrototypes@"   -> interfacePrototypes,
          "@FFVFormFactors_InterfaceDefinitions@"  -> interfaceDefinitions,
          "@TemplateWrapper_InterfacePrototypes@" -> templateWrapperDecl,
          "@TemplateWrapper_InterfaceDefinitions@" -> templateWrapperDef,
          Sequence @@ GeneralReplacementRules[]}
      ];

      vertices
   ];

WriteBToSGammaClass[decays_List, files_List] :=
   Module[{createInterface = False,
          btosgammaInterfaceDefinitions},
       If[decays =!= {},
         createInterface = True];
       btosgammaInterfaceDefinitions =
         BtoSGamma`CreateInterfaceBtoSGamma[createInterface];

       WriteOut`ReplaceInFiles[files, {
            "@BtoSGammaInterface@" -> btosgammaInterfaceDefinitions,
            Sequence @@ GeneralReplacementRules[]}];
   ];

WriteUnitarityClass[files_List] :=
    Module[{res},
      res = If[FSUnitarityConstraints, Unitarity`GetScatteringMatrix[], {0, ""}];
      WriteOut`ReplaceInFiles[files,
        {"@scatteringPairsLength@" -> ToString@res[[1]],
         "@infiniteSMatrix@" -> res[[2]],
         Sequence @@ GeneralReplacementRules[]
        }];
    ];

(* Write the AMM c++ files *)
WriteAMMClass[fields_List, files_List] :=
    Module[{calculation, getMLCP,
            (* in models without flavour violation (no FV models) lepton does not have an index *)
            leptonIndex = If[Length[fields] > 0, If[TreeMasses`GetDimension[First@fields] =!= 1, "idx", ""], ""],
            (* we want to calculate an offset of g-2 compared to the SM *)
            discardSMcontributions = CConversion`CreateCBoolValue[True],
            graphs, diagrams, vertices, barZee = "", calculateForwadDeclaration, uncertaintyForwadDeclaration, leptonPoleMass,
            ammWrapperDecl, ammWrapperDef, ammUncWrapperDecl, ammUncWrapperDef},

      calculation =
         If[Length[fields] =!= 0,
            "const auto form_factors = " <>
            FSModelName <> "_FFV_form_factors::calculate_form_factors<Lepton,Lepton," <>
            CXXDiagrams`CXXNameOfField[TreeMasses`GetPhoton[]] <> ">(" <>
            leptonIndex <> If[leptonIndex === "", "", ", "] <>
            leptonIndex <> If[leptonIndex === "", "", ", "] <>
            "model, " <> discardSMcontributions <> ");",
            "const std::valarray<std::complex<double>> form_factors {0., 0., 0., 0.};"
         ];

      getMLCP = AMM`AMMGetMLCP[];

      (* only Barr-Zee graphs
         1-loop diagrams will be provided by the FFMasslessV module and are taken care of elsewhere *)
      graphs = If[Length[fields] > 0, AMM`AMMContributingGraphs[], {}];
      diagrams = If[Length[fields] > 0, Flatten[Outer[AMM`AMMContributingDiagramsForGraph, graphs, fields, 1], 1], {}];

      vertices = Flatten[CXXDiagrams`VerticesForDiagram /@ Flatten[diagrams, 1], 1];
      calculateForwadDeclaration = StringRiffle[AMM`AMMForwardDeclaration[#, "calculate_amm"]& /@ fields, "\n"];
      uncertaintyForwadDeclaration = StringRiffle[AMM`AMMForwardDeclaration[#, "calculate_amm_uncertainty"]& /@ fields, "\n"];

      For[i = 1, i <= Length[graphs], i++,
         For[j = 1, j <= Length[diagrams[[i]]], j++,
            barZee = barZee <>
               "val += " <> ToString @ N[
                  With[{colFac = CXXDiagrams`ColourFactorForIndexedDiagramFromGraph[
                                    CXXDiagrams`IndexDiagramFromGraph[diagrams[[i,j]], graphs[[i]]],
                                    graphs[[i]]
                                 ]},
                     If[Im[colFac] == 0, colFac, Print["Error: Colour prefactor of a Barr-Zee diagram should be real!"]; Quit[1]]
                  ], 16] <> " * " <>
                ToString @ AMM`CXXEvaluatorForDiagramFromGraph[diagrams[[i,j]], graphs[[i]]] <>
                "::value({" <> leptonIndex <> "}, context, ml_pole);\n"
         ];
      ];

      leptonPoleMass =
         If[GetParticleFromDescription["Leptons"] =!= Null,
"template <typename Lepton>
double lepton_pole_mass(const softsusy::QedQcd& qedqcd, int idx)
{
   return qedqcd.displayLeptonPoleMass(idx);
}",
StringRiffle[
(
"template <typename Lepton>
std::enable_if_t<std::is_same_v<Lepton, " <> CXXDiagrams`CXXNameOfField[#, prefixNamespace-> FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"] <> ">, double>
lepton_pole_mass(const softsusy::QedQcd& qedqcd)
{\n" <>
TextFormatting`IndentText[
   Switch[#,
      GetParticleFromDescription["Electron"], "return qedqcd.displayPoleMel()",
      GetParticleFromDescription["Muon"],     "return qedqcd.displayPoleMmuon()",
      GetParticleFromDescription["Tau"],      "return qedqcd.displayPoleMtau()"
   ] <> ";\n"
] <>
"}")& /@ If[fields =!= {}, fields, {GetParticleFromDescription["Electron"], GetParticleFromDescription["Muon"], GetParticleFromDescription["Tau"]}],
"\n\n"
]
         ];

      ammWrapperDecl = StringRiffle[("double calculate_" <> CXXDiagrams`CXXNameOfField[#] <> "_amm(const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings" <> If[GetParticleFromDescription["Leptons"] =!= Null, ", int idx", ""] <> ");")& /@ fields, "\n"];
      ammUncWrapperDecl = StringRiffle[("double calculate_" <> CXXDiagrams`CXXNameOfField[#] <> "_amm_uncertainty(const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings" <> If[GetParticleFromDescription["Leptons"] =!= Null, ", int idx", ""] <> ");")& /@ fields, "\n"];
      ammWrapperDef = StringRiffle[
("double calculate_" <> CXXDiagrams`CXXNameOfField[#] <> "_amm(const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings" <> If[GetParticleFromDescription["Leptons"] =!= Null, ", int idx", ""] <> ") {
   return " <> FlexibleSUSY`FSModelName <> "_amm::calculate_amm<fields::" <> CXXDiagrams`CXXNameOfField[#] <> ">(model, qedqcd, settings" <> If[GetParticleFromDescription["Leptons"] =!= Null, ", idx", ""] <> ");
}")& /@ fields, "\n"
      ];
      ammUncWrapperDef = StringRiffle[
("double calculate_" <> CXXDiagrams`CXXNameOfField[#] <> "_amm_uncertainty(const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings& settings" <> If[GetParticleFromDescription["Leptons"] =!= Null, ", int idx", ""] <> ") {
   return " <> FlexibleSUSY`FSModelName <> "_amm::calculate_amm_uncertainty<fields::" <> CXXDiagrams`CXXNameOfField[#] <> ">(model, qedqcd, settings" <> If[GetParticleFromDescription["Leptons"] =!= Null, ", idx", ""] <> ");
}")& /@ fields, "\n"
      ];

      WriteOut`ReplaceInFiles[files,
        {"@AMMZBosonField@"       -> CXXDiagrams`CXXNameOfField[TreeMasses`GetZBoson[]],
         "@AMMCalculation@"       -> Nest[TextFormatting`IndentText, calculation, 1],
         "@AMMGetMLCP@"           -> TextFormatting`IndentText[WrapLines[getMLCP]],
         "@calculateAForwardDeclaration@" -> calculateForwadDeclaration,
         "@calculateAUncertaintyForwardDeclaration@" -> uncertaintyForwadDeclaration,
         "@extraIdxDecl@" -> If[GetParticleFromDescription["Leptons"] =!= Null, ", int idx", ""],
         "@extraIdxUsage@" -> If[GetParticleFromDescription["Leptons"] =!= Null, ", idx", ""],
         "@extraIdxUsageNoComma@" -> If[GetParticleFromDescription["Leptons"] =!= Null, "idx", ""],
         "@leptonPoleMass@" -> leptonPoleMass,
         "@AMMBarrZeeCalculation@" -> Nest[TextFormatting`IndentText, barZee, 1],
         "@ammWrapperDecl@" -> ammWrapperDecl,
         "@ammWrapperDef@" -> ammWrapperDef,
         "@ammUncWrapperDecl@" -> ammUncWrapperDecl,
         "@ammUncWrapperDef@" -> ammUncWrapperDef,
         Sequence @@ GeneralReplacementRules[]
        }];

        vertices
      ];

GetBVPSolverHeaderName[solver_] :=
    Switch[solver,
           FlexibleSUSY`TwoScaleSolver, "two_scale",
           FlexibleSUSY`SemiAnalyticSolver, "semi_analytic",
           FlexibleSUSY`ShootingSolver, "shooting",
           FlexibleSUSY`LatticeSolver, "lattice",
           _, Print["Error: invalid BVP solver requested: ", solver];
              Quit[1];
          ];

GetBVPSolverEnabledMacro[solver_] :=
    Switch[solver,
           FlexibleSUSY`TwoScaleSolver, "ENABLE_TWO_SCALE_SOLVER",
           FlexibleSUSY`SemiAnalyticSolver, "ENABLE_SEMI_ANALYTIC_SOLVER",
           FlexibleSUSY`ShootingSolver, "ENABLE_SHOOTING_SOLVER",
           FlexibleSUSY`LatticeSolver, "ENABLE_LATTICE_SOLVER",
           _, Print["Error: invalid BVP solver requested: ", solver];
              Quit[1];
          ];

GetBVPSolverSLHAOptionKey[solver_] :=
    Switch[solver,
           FlexibleSUSY`TwoScaleSolver, "1",
           FlexibleSUSY`SemiAnalyticSolver, "2",
           FlexibleSUSY`LatticeSolver, "3",
           FlexibleSUSY`ShootingSolver, "4",
           _, Print["Error: invalid BVP solver requested: ", solver];
              Quit[1];
          ];

GetBVPSolverTemplateParameter[solver_] :=
    Switch[solver,
           FlexibleSUSY`TwoScaleSolver, "Two_scale",
           FlexibleSUSY`SemiAnalyticSolver, "Semi_analytic",
           FlexibleSUSY`ShootingSolver, "Shooting",
           FlexibleSUSY`LatticeSolver, "Lattice",
           _, Print["Error: invalid BVP solver requested: ", solver];
              Quit[1];
          ];

EnableForBVPSolver[solver_, statements_String] :=
    "#ifdef " <> GetBVPSolverEnabledMacro[solver] <> "\n" <> statements <> "#endif";

EnableSpectrumGenerator[solver_] :=
    Module[{header = "#include \"" <> FlexibleSUSY`FSModelName},
           header = header <> "_" <> GetBVPSolverHeaderName[solver];
           header = header <> "_spectrum_generator.hpp\"\n";
           EnableForBVPSolver[solver, header] <> "\n"
          ];

RunEnabledSpectrumGenerator[solver_] :=
    Module[{key = "", class = "", body = "", result = ""},
           key = GetBVPSolverSLHAOptionKey[solver];
           class = GetBVPSolverTemplateParameter[solver];
           body = "exit_code = run_solver<" <> class <> ">(\n"
                  <> IndentText[
                        "slha_io, spectrum_generator_settings, " <>
                        If[FSCalculateDecays, "flexibledecay_settings, ", ""] <>
                        "slha_output_file,\n"]
                  <> IndentText["database_output_file, spectrum_file, rgflow_file, higgsbounds_dataset, higgssignals_dataset, lilith_db);\n"]
                  <> "if (!exit_code || solver_type != 0) break;\n"
                  <> "[[fallthrough]];\n";
           result = "case " <> key <> ":\n" <> IndentText[body];
           EnableForBVPSolver[solver, IndentText[result]] <> "\n"
          ];

ScanEnabledSpectrumGenerator[solver_] :=
    Module[{key = "", class = "", body = "", result = ""},
           key = GetBVPSolverSLHAOptionKey[solver];
           class = GetBVPSolverTemplateParameter[solver];
           body = "result = run_parameter_point<" <> class <> ">(loop_library, qedqcd, input);\n"
                  <> "if (!result.problems.have_problem() || solver_type != 0) break;\n"
                  <> "[[fallthrough]];\n";
           result = "case " <> key <> ":\n" <> IndentText[body];
           EnableForBVPSolver[solver, IndentText[IndentText[result]]] <> "\n"
          ];

RunCmdLineEnabledSpectrumGenerator[solver_] :=
    Module[{key = "", class = "", body = "", result = ""},
           key = GetBVPSolverSLHAOptionKey[solver];
           class = GetBVPSolverTemplateParameter[solver];
           body = "exit_code = run_solver<" <> class <> ">(loop_library,input);\n"
                  <> "if (!exit_code || solver_type != 0) break;\n"
                  <> "[[fallthrough]];\n";
           result = "case " <> key <> ":\n" <> IndentText[body];
           EnableForBVPSolver[solver, IndentText[result]] <> "\n"
          ];

ExampleDecaysIncludes[] :=
    Utils`StringJoinWithSeparator[
       ("#include \"" <> # <> "\"")& /@ {
         "decays/" <> FlexibleSUSY`FSModelName <> "_decays.hpp",
         "decays/flexibledecay_problems.hpp",
         FlexibleSUSY`FSModelName <> "_mass_eigenstates_decoupling_scheme.hpp",
         "loop_libraries/loop_library.hpp"},
       "\n"
    ];

ExampleCalculateDecaysForModel[] :=
IndentText[
"if (flexibledecay_settings.get(FlexibleDecay_settings::calculate_decays) &&
     (spectrum_generator_settings.get(Spectrum_generator_settings::force_output) ||
      !problems.have_problem())) {
   if (loop_library_for_decays) {
      decays.calculate_decays();\n"
] <>
IndentText@IndentText@IndentText[
"effc =
   get_normalized_effective_couplings(decays.get_neutral_higgs_effc(), physical_input, qedqcd, spectrum_generator_settings, flexibledecay_settings);\n" <>
IndentText[
   "try {\n" <>
   IndentText[
      "// structured bindings creates new variables - need to use std::tie
#ifdef ENABLE_HIGGSTOOLS
if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools)) {
         std::tie(hs, higgsbounds_v) =
            call_higgstools(effc, physical_input, higgsbounds_dataset, higgssignals_dataset);
}
#endif
#ifdef ENABLE_LILITH
if (flexibledecay_settings.get(FlexibleDecay_settings::call_lilith)) {
   lilith = call_lilith(effc, physical_input, lilith_db);
}
#endif\n"
   ] <>
   "}\n" <>
   "catch (const std::exception& error) {\n" <>
      IndentText["ERROR(error.what());\n"] <>
   "}\n"
]
] <>
IndentText[IndentText[
"}
else if (!loop_library_for_decays) {
   WARNING(\"Decay module requires a dedicated loop library. Configure FlexibleSUSY with Collier or LoopTools and set appropriately flag 31 in Block FlexibleSUSY of the LesHouches input.\");
}"
] <> "\n}"
];

ExampleSetDecaysSLHAOutput[] := "\
const bool show_decays = !decays.get_problems().have_problem() ||
   spectrum_generator_settings.get(Spectrum_generator_settings::force_output);

if (show_decays && flexibledecay_settings.get(FlexibleDecay_settings::calculate_decays) && loop_library_for_decays) {
   slha_io.set_dcinfo(decays.get_problems());
   slha_io.set_decays(decays.get_decay_table(), flexibledecay_settings);
   if (flexibledecay_settings.get(FlexibleDecay_settings::print_effc_block)) {
      slha_io.set_effectivecouplings_block(decays.get_effhiggscouplings_block_input());
   }
   if (flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc)) {
      slha_io.set_normalized_effectivecouplings_block(effc);\n" <>
      If[SA`CPViolationHiggsSector || GetDimensionWithoutGoldstones@TreeMasses`GetPseudoscalarHiggsBoson[] != 0, IndentText@IndentText["slha_io.set_imnormalized_effectivecouplings_block(effc);\n"], ""] <>
"   }
#ifdef ENABLE_HIGGSTOOLS
   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools)) {
      if (hs.has_value()) {
         SignalResult hs_ = hs.value();
         slha_io.set_hs_or_lilith(\"HIGGSSIGNALS\", hs_.ndof, hs_.chi2BSM, hs_.chi2SM, hs_.mhRef, hs_.pval);
      }
      if (higgsbounds_v.size() > 0) {
         slha_io.set_higgsbounds(higgsbounds_v);
      }
   }
#endif
#ifdef ENABLE_LILITH
   if (flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) && lilith.has_value()) {
      slha_io.set_hs_or_lilith(\"LILITH\", lilith.value().ndof, lilith.value().chi2BSM, lilith.value().chi2SM, lilith.value().mhRef, lilith.value().pval);
   }
#endif
}";

ExampleCalculateCmdLineDecays[] :=
FlexibleSUSY`FSModelName <> "_decays decays " <>
"= " <> FlexibleSUSY`FSModelName <> "_decays(std::get<0>(models), qedqcd, physical_input, flexibledecay_settings);
const bool loop_library_for_decays =
   (Loop_library::get_type() == Loop_library::Library::Collier) ||
   (Loop_library::get_type() == Loop_library::Library::Looptools);
if (spectrum_generator.get_exit_code() == 0 && loop_library_for_decays) {
   decays.calculate_decays();
}";

WriteExampleCmdLineOutput[enableDecays_] :=
    If[enableDecays,
       "SLHAea::Coll slhaea(" <> FlexibleSUSY`FSModelName <> "_slha_io::fill_slhaea(\n" <>
       "                       model, qedqcd, scales, observables, decays));",
       "SLHAea::Coll slhaea(" <> FlexibleSUSY`FSModelName <> "_slha_io::fill_slhaea(\n" <>
       "                       model, qedqcd, scales, observables));"
      ];

WriteUserExample[inputParameters_List, files_List] :=
    Module[{parseCmdLineOptions, printCommandLineOptions, inputPars,
            solverIncludes = "", runEnabledSolvers = "", scanEnabledSolvers = "",
            runEnabledCmdLineSolvers = "", defaultSolverType,
            decaysIncludes = "", calculateDecaysForModel = "",
            decaysObject = "", decaySettingsObj = "", fillDecaySettings = "",
            setDecaysSLHAOutput = "", decaySetttingsOverride = "", flexibleDecaySettingsVarInDecl = "", flexibleDecaySettingsVarInDef = "",
            calculateCmdLineDecays = "", writeCmdLineOutput = "", fillSLHAIO = ""},
           inputPars = {First[#], #[[3]]}& /@ inputParameters;
           parseCmdLineOptions = WriteOut`ParseCmdLineOptions[inputPars];
           printCommandLineOptions = WriteOut`PrintCmdLineOptions[inputPars];
           (solverIncludes = solverIncludes <> EnableSpectrumGenerator[#])& /@ FlexibleSUSY`FSBVPSolvers;
           (runEnabledSolvers = runEnabledSolvers <> RunEnabledSpectrumGenerator[#])& /@ FlexibleSUSY`FSBVPSolvers;
           (scanEnabledSolvers = scanEnabledSolvers <> ScanEnabledSpectrumGenerator[#])& /@ FlexibleSUSY`FSBVPSolvers;
           (runEnabledCmdLineSolvers = runEnabledCmdLineSolvers <> RunCmdLineEnabledSpectrumGenerator[#])& /@ FlexibleSUSY`FSBVPSolvers;
           If[Length[FlexibleSUSY`FSBVPSolvers] == 0,
              defaultSolverType = "-1",
              defaultSolverType = GetBVPSolverSLHAOptionKey[FlexibleSUSY`FSBVPSolvers[[1]]]
             ];
           If[FlexibleSUSY`FSCalculateDecays,
              decaysIncludes = ExampleDecaysIncludes[];
              decaysObject =
                  IndentText[
                     FlexibleSUSY`FSModelName <> "_decays decays(std::get<0>(models), qedqcd, physical_input, flexibledecay_settings);\n" <>
                     "const bool loop_library_for_decays =\n" <>
                     IndentText[
                        "(Loop_library::get_type() == Loop_library::Library::Collier) ||\n" <>
                        "(Loop_library::get_type() == Loop_library::Library::Looptools);\n"
                     ]
                  ];
              decaySettingsObj = "FlexibleDecay_settings flexibledecay_settings;\n";
              fillDecaySettings = "slha_io.fill(flexibledecay_settings);\n";
              flexibleDecaySettingsVarInDecl = "flexibledecay_settings, ";
              flexibleDecaySettingsVarInDef = "const flexiblesusy::FlexibleDecay_settings& flexibledecay_settings,";
              calculateDecaysForModel = ExampleCalculateDecaysForModel[];
              setDecaysSLHAOutput = ExampleSetDecaysSLHAOutput[];
              calculateCmdLineDecays = ExampleCalculateCmdLineDecays[];
              fillSLHAIO = "slha_io.fill(models, qedqcd, scales, observables, settings, flexibledecay_settings, &decays);";
              decaySetttingsOverride =
"if (flexibledecay_settings.get(FlexibleDecay_settings::calculate_decays)) {
   if (!spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses)) {
      WARNING(\"Decay module requires SM pole masses. Setting FlexibleSUSY[3] = 1.\");
      spectrum_generator_settings.set(
         Spectrum_generator_settings::calculate_sm_masses, 1.0);
   }
   if (!spectrum_generator_settings.get(Spectrum_generator_settings::calculate_bsm_masses)) {
      WARNING(\"Decay module requires BSM pole masses. Setting FlexibleSUSY[23] = 1.\");
      spectrum_generator_settings.set(
         Spectrum_generator_settings::calculate_bsm_masses, 1.0);
   }
   if (flexibledecay_settings.get(FlexibleDecay_settings::print_effc_block) && " <> FlexibleSUSY`FSModelName <> "_info::is_CP_violating_Higgs_sector) {
      WARNING(\"Printing of EFFHIGGSCOUPLINGS block is disabled in models with CP-violating Higgs sector\");
      flexibledecay_settings.set(FlexibleDecay_settings::print_effc_block, 0.);
   }
}",
              fillSLHAIO = "slha_io.fill(models, qedqcd, scales, observables, settings, flexibledecay_settings);"
             ];
           writeCmdLineOutput = WriteExampleCmdLineOutput[FlexibleSUSY`FSCalculateDecays];
           WriteOut`ReplaceInFiles[files,
                          { "@parseCmdLineOptions@" -> IndentText[IndentText[parseCmdLineOptions]],
                            "@printCommandLineOptions@" -> IndentText[IndentText[printCommandLineOptions]],
                            "@solverIncludes@" -> solverIncludes,
                            "@runEnabledSolvers@" -> runEnabledSolvers,
                            "@scanEnabledSolvers@" -> scanEnabledSolvers,
                            "@runEnabledCmdLineSolvers@" -> runEnabledCmdLineSolvers,
                            "@defaultSolverType@" -> defaultSolverType,
                            "@decaysIncludes@" -> decaysIncludes,
                            "@decaysObject@" -> decaysObject,
                            "@decaySettingsObj@" -> IndentText@decaySettingsObj,
                            "@fillDecaySettings@" -> IndentText@IndentText@fillDecaySettings,
                            "@flexibleDecaySettingsVarInDef@" -> flexibleDecaySettingsVarInDef,
                            "@flexibleDecaySettingsVarInDecl@" -> flexibleDecaySettingsVarInDecl,
                            "@calculateDecaysForModel@" -> calculateDecaysForModel,
                            "@setDecaysSLHAOutput@" -> IndentText[IndentText[setDecaysSLHAOutput]],
                            "@calculateCmdLineDecays@" -> IndentText[calculateCmdLineDecays],
                            "@writeCmdLineOutput@" -> IndentText[writeCmdLineOutput],
                            "@fillSLHAIO@" -> fillSLHAIO,
                            "@decaySetttingsOverride@" -> IndentText[decaySetttingsOverride],
                            "@calculateUnitarity@" -> If[FSUnitarityConstraints,
                                 "const auto unitarityStruct = " <> FSModelName <> "_unitarity::max_scattering_eigenvalue_infinite_s(std::get<0>(models));\n" <>
                                  "slha_io.set_unitarity_infinite_s(spectrum_generator_settings, unitarityStruct);",
                                  ""
                               ],
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

EnableMathlinkSpectrumGenerator[solver_] :=
    Module[{type, headers = ""},
           type = GetBVPSolverHeaderName[solver];
           headers = "#include \"" <> FlexibleSUSY`FSModelName <> "_" <> type <> "_model.hpp\"\n";
           headers = headers <> "#include \"" <> FlexibleSUSY`FSModelName <> "_" <> type <> "_spectrum_generator.hpp\"\n";
           EnableForBVPSolver[solver, headers] <> "\n"
          ];

RunEnabledModelType[solver_] :=
    Module[{key, class, body, result},
           key = GetBVPSolverSLHAOptionKey[solver];
           class = GetBVPSolverTemplateParameter[solver];
           body = "spectrum.reset(new " <> FlexibleSUSY`FSModelName <>
                  "_spectrum_impl<" <> class <> ">());\n"
                  <> "spectrum->calculate_spectrum(settings, modsel, qedqcd, input);\n"
                  <> "if (!spectrum->get_problems().have_problem() || solver_type != 0) break;\n"
                  <> "[[fallthrough]];\n";
           result = "case " <> key <> ":\n" <> IndentText[body];
           EnableForBVPSolver[solver, IndentText[result]] <> "\n"
          ];

WriteMathLink[inputParameters_List, extraSLHAOutputBlocks_List, files_List] :=
    Module[{numberOfInputParameters, numberOfInputParameterRules,
            putInputParameters,
            setInputParameterDefaultArguments,
            setInputParameterArguments,
            numberOfSpectrumEntries, putSpectrum, setInputParameters,
            numberOfObservables, putObservables,
            inputPars, outPars, requestedObservables, defaultSolverType,
            solverIncludes = "", runEnabledSolvers = "",
            decaysData = "", calculateDecaysVirtualFunc = "", calculateDecaysEffCVirtualFunc = "", calculateSpectrumDecaysPrototype = "", calculateSpectrumDecaysEffCPrototype = "",
            calculateSpectrumDecaysFunction = "", calculateModelDecaysPrototype = "", calculateModelDecaysEffCFunction = "",
            calculateSpectrumDecaysEffCFunction = "", calculateModelDecaysEffCPrototype = "",
            calculateModelDecaysFunction = "", fillDecaysSLHA = "", getDecaysVirtualFunc = "",
            getSpectrumDecays = "", putDecaysPrototype = "", putDecaysFunction = "", putEffCPrototype = "", putEffCFunction = "", callLilithMessages = "",
            mathlinkDecaysCalculationFunction = "", loadCalculateDecaysFunction = "",
            mathlinkCalcNormalizedEffC = "", loadCalculateEffCFunction = "", loadCallHiggsToolsFunction = "", loadCallLilithFunction = "",
            setUnitarity = "", loadCalculateUnitarityFunction = "", calculateUnitarityMessages = "",
            calculateDecaysMessages = "", calculateDecaysExample = "", decaysIncludes = "", fdDefaultSettings = "",
            addFDOptions1 = "", addFDOptions2 = "", setFDOptions = "", setDecayOptions = "", fillFDSettings = "",
            decayIndex = "const Index_t n_fd_settings = 0;"},
           inputPars = {#[[1]], #[[3]]}& /@ inputParameters;
           numberOfInputParameters = Total[CConversion`CountNumberOfEntries[#[[2]]]& /@ inputPars];
           numberOfInputParameterRules = FSMathLink`GetNumberOfInputParameterRules[inputPars];
           putInputParameters = FSMathLink`PutInputParameters[inputPars, "link"];
           setInputParameters = FSMathLink`SetInputParametersFromArguments[inputPars];
           setInputParameterDefaultArguments = FSMathLink`SetInputParameterDefaultArguments[inputPars];
           setInputParameterArguments = FSMathLink`SetInputParameterArguments[inputPars];
           outPars = Parameters`GetOutputParameters[] /. FlexibleSUSY`M[p_List] :> Sequence @@ (FlexibleSUSY`M /@ p);
           outPars = Join[outPars, FlexibleSUSY`Pole /@ outPars, Parameters`GetModelParameters[],
                          Parameters`GetExtraParameters[], {FlexibleSUSY`SCALE}];
           numberOfSpectrumEntries = FSMathLink`GetNumberOfSpectrumEntries[outPars];
           putSpectrum = FSMathLink`PutSpectrum[outPars, "link"];
           (* get observables *)
           requestedObservables = Observables`GetRequestedObservables[extraSLHAOutputBlocks];
           numberOfObservables = Length[requestedObservables];
           putObservables = FSMathLink`PutObservables[requestedObservables, "link"];
           (solverIncludes = solverIncludes <> EnableMathlinkSpectrumGenerator[#])& /@ FlexibleSUSY`FSBVPSolvers;
           (runEnabledSolvers = runEnabledSolvers <> RunEnabledModelType[#])& /@ FlexibleSUSY`FSBVPSolvers;
           If[Length[FlexibleSUSY`FSBVPSolvers] == 0,
              defaultSolverType = "-1",
              defaultSolverType = GetBVPSolverSLHAOptionKey[FlexibleSUSY`FSBVPSolvers[[1]]];
             ];
           If[FSUnitarityConstraints,
              loadCalculateUnitarityFunction = "FS" <> FlexibleSUSY`FSModelName <> "CalculateUnitarity = LibraryFunctionLoad[lib" <>
                                            FlexibleSUSY`FSModelName <> ", \"FS" <> FlexibleSUSY`FSModelName <>
                                            "CalculateUnitarity\", LinkObject, LinkObject];\n";
              calculateUnitarityMessages = "\n" <> "FS" <> FlexibleSUSY`FSModelName <> "CalculateUnitarity::error = \"`1`\";\n" <>
                                        "FS" <> FlexibleSUSY`FSModelName <> "CalculateUnitarity::warning = \"`1`\";\n";
              setUnitarity = "slha_io.set_unitarity_infinite_s(settings, unitarityData);";
           ];
           If[FlexibleSUSY`FSCalculateDecays,
              decaysData = FlexibleSUSY`FSModelName <> "_decays decays{};              ///< decays";
              getDecaysVirtualFunc = FSMathLink`CreateSpectrumDecaysGetterInterface[FlexibleSUSY`FSModelName];
              getSpectrumDecays = CreateSpectrumDecaysGetter[FlexibleSUSY`FSModelName];
              calculateDecaysVirtualFunc = FSMathLink`CreateSpectrumDecaysInterface[FlexibleSUSY`FSModelName];
              calculateDecaysEffCVirtualFunc = FSMathLink`CreateSpectrumDecaysEffCInterface[FlexibleSUSY`FSModelName];
              {calculateSpectrumDecaysPrototype, calculateSpectrumDecaysFunction} =
                  FSMathLink`CreateSpectrumDecaysCalculation[FlexibleSUSY`FSModelName];
              {calculateSpectrumDecaysEffCPrototype, calculateSpectrumDecaysEffCFunction} =
                  FSMathLink`CreateSpectrumDecaysEffCCalculation[FlexibleSUSY`FSModelName];
              {calculateModelDecaysPrototype, calculateModelDecaysFunction} =
                  FSMathLink`CreateModelDecaysCalculation[FlexibleSUSY`FSModelName];
              {calculateModelDecaysEffCPrototype, calculateModelDecaysEffCFunction} =
                  FSMathLink`CreateModelDecaysEffCCalculation[FlexibleSUSY`FSModelName];
              fillDecaysSLHA = FSMathLink`FillDecaysSLHAData[];
              {putDecaysPrototype, putDecaysFunction} = FSMathLink`PutDecays[FlexibleSUSY`FSModelName];
              {putEffCPrototype, putEffCFunction} = FSMathLink`PutEffCouplings[FlexibleSUSY`FSModelName];
              mathlinkDecaysCalculationFunction = FSMathLink`CreateMathLinkDecaysCalculation[FlexibleSUSY`FSModelName];
              loadCalculateDecaysFunction = "FS" <> FlexibleSUSY`FSModelName <> "CalculateDecays = LibraryFunctionLoad[lib" <>
                                            FlexibleSUSY`FSModelName <> ", \"FS" <> FlexibleSUSY`FSModelName <>
                                            "CalculateDecays\", LinkObject, LinkObject];\n";
              loadCalculateEffCFunction = "FS" <> FlexibleSUSY`FSModelName <> "CalculateNormalizedEffectiveCouplings = LibraryFunctionLoad[lib" <>
                                            FlexibleSUSY`FSModelName <> ", \"FS" <> FlexibleSUSY`FSModelName <>
                                            "CalculateNormalizedEffectiveCouplings\", LinkObject, LinkObject];";
              loadCallHiggsToolsFunction = "FS" <> FlexibleSUSY`FSModelName <> "CallHiggsTools = LibraryFunctionLoad[lib" <>
                                            FlexibleSUSY`FSModelName <> ", \"FS" <> FlexibleSUSY`FSModelName <>
                                            "CallHiggsTools\", LinkObject, LinkObject];";
              loadCallLilithFunction = "FS" <> FlexibleSUSY`FSModelName <> "CallLilith = LibraryFunctionLoad[lib" <>
                                            FlexibleSUSY`FSModelName <> ", \"FS" <> FlexibleSUSY`FSModelName <>
                                            "CallLilith\", LinkObject, LinkObject];";
              calculateDecaysMessages = "\n" <> "FS" <> FlexibleSUSY`FSModelName <> "CalculateDecays::error = \"`1`\";\n" <>
                                        "FS" <> FlexibleSUSY`FSModelName <> "CalculateDecays::warning = \"`1`\";\n";
              callLilithMessages = "\n" <> "FS" <> FlexibleSUSY`FSModelName <> "CallLilith::error = \"`1`\";\n" <>
                                        "FS" <> FlexibleSUSY`FSModelName <> "CallLilith::warning = \"`1`\";\n";
              calculateDecaysExample = "decays      = FS" <> FlexibleSUSY`FSModelName <> "CalculateDecays[handle];\n";
              decaysIncludes = "#include \"loop_libraries/loop_library.hpp\"";
              mathlinkCalcNormalizedEffC = FSMathLink`CalculateNormalizedEffectiveCouplings[FlexibleSUSY`FSModelName];
              fdDefaultSettings =
"\nfdDefaultSettings = {
   minBRtoPrint -> 1*^-5,
   maxHigherOrderCorrections -> 4,
   alphaThomson -> 1,
   offShellVV -> 2,
   printEffCBlock -> 1,
   calcNormalizedEffC -> 0,
   callHiggsTools -> 0,
   callLilith -> 0,
   usePoleHiggsMixings -> 1
};\n";
               addFDOptions1 = ", Sequence @@ fdDefaultSettings";
               addFDOptions2 = "| fdSettings";
setFDOptions =
",
OptionValue[minBRtoPrint],
OptionValue[maxHigherOrderCorrections],
OptionValue[alphaThomson],
OptionValue[offShellVV],
OptionValue[printEffCBlock],
OptionValue[calcNormalizedEffC],
OptionValue[callHiggsTools],
OptionValue[callLilith],
OptionValue[usePoleHiggsMixings]";
setDecayOptions =
"FlexibleDecay_settings flexibledecay_settings;
flexibledecay_settings.set(FlexibleDecay_settings::min_br_to_print, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::offshell_VV_decays, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::print_effc_block, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::calculate_normalized_effc, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::call_higgstools, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::call_lilith, pars[c++]);
flexibledecay_settings.set(FlexibleDecay_settings::use_pole_higgs_mixings, pars[c++]);
";
decayIndex = "const Index_t n_fd_settings = 9;";
fillFDSettings = "data.set_fd_settings(flexibledecay_settings);\n"
             ];
           WriteOut`ReplaceInFiles[files,
                          { "@numberOfInputParameters@" -> ToString[numberOfInputParameters],
                            "@numberOfInputParameterRules@" -> ToString[numberOfInputParameterRules],
                            "@putInputParameters@" -> IndentText[putInputParameters],
                            "@setInputParameters@" -> IndentText[setInputParameters],
                            "@setInputParameterArguments@" -> IndentText[setInputParameterArguments, 12],
                            "@setInputParameterDefaultArguments@" -> IndentText[setInputParameterDefaultArguments],
                            "@setDefaultInputParameters@" -> IndentText[setInputParameterDefaultArguments,8],
                            "@numberOfSpectrumEntries@" -> ToString[numberOfSpectrumEntries],
                            "@putSpectrum@" -> IndentText[putSpectrum],
                            "@numberOfObservables@" -> ToString[numberOfObservables],
                            "@putObservables@" -> IndentText[putObservables],
                            "@solverIncludes@" -> solverIncludes,
                            "@runEnabledSolvers@" -> runEnabledSolvers,
                            "@defaultSolverType@" -> defaultSolverType,
                            "@calculateDecaysVirtualFunc@" -> IndentText[calculateDecaysVirtualFunc],
                            "@calculateDecaysEffCVirtualFunc@" -> IndentText[calculateDecaysEffCVirtualFunc],
                            "@calculateSpectrumDecaysPrototype@" -> IndentText[calculateSpectrumDecaysPrototype],
                            "@calculateSpectrumDecaysEffCPrototype@" -> IndentText[calculateSpectrumDecaysEffCPrototype],
                            "@calculateSpectrumDecaysFunction@" -> calculateSpectrumDecaysFunction,
                            "@calculateSpectrumDecaysEffCFunction@" -> calculateSpectrumDecaysEffCFunction,
                            "@calculateModelDecaysPrototype@" -> IndentText[calculateModelDecaysPrototype],
                            "@calculateModelDecaysFunction@" -> calculateModelDecaysFunction,
                            "@calculateModelDecaysEffCPrototype@" -> IndentText[calculateModelDecaysEffCPrototype],
                            "@calculateModelDecaysEffCFunction@" -> calculateModelDecaysEffCFunction,
                            "@decaysData@" -> IndentText[decaysData],
                            "@fillDecaysSLHA@" -> IndentText[fillDecaysSLHA],
                            "@getDecaysVirtualFunc@" -> IndentText[getDecaysVirtualFunc],
                            "@getSpectrumDecays@" -> IndentText[getSpectrumDecays],
                            "@putDecaysPrototype@" -> IndentText[putDecaysPrototype],
                            "@putDecaysFunction@" -> putDecaysFunction,
                            "@putEffCPrototype@" -> IndentText[putEffCPrototype],
                            "@putEffCFunction@" -> putEffCFunction,
                            "@mathlinkDecaysCalculationFunction@" -> mathlinkDecaysCalculationFunction,
                            "@loadCalculateDecaysFunction@" -> loadCalculateDecaysFunction,
                            "@loadCalculateEffCFunction@" -> loadCalculateEffCFunction,
                            "@loadCallHiggsToolsFunction@" -> loadCallHiggsToolsFunction,
                            "@callLilithMessages@" -> callLilithMessages,
                            "@loadCallLilithFunction@" -> loadCallLilithFunction,
                            "@calculateDecaysMessages@" -> calculateDecaysMessages,
                            "@calculateDecaysExample@" -> calculateDecaysExample,
                            "@decaysIncludes@" -> decaysIncludes,
                            "@fdDefaultSettings@" -> fdDefaultSettings,
                            "@addFDOptions1@" -> IndentText[addFDOptions1],
                            "@addFDOptions2@" -> addFDOptions2,
                            "@setFDOptions@" -> IndentText @ IndentText @ IndentText @ IndentText @ setFDOptions,
                            "@setDecayOptions@" -> IndentText @ setDecayOptions,
                            "@fillFDSettings@" -> fillFDSettings,
                            "@decayIndex@" -> decayIndex,
                            "@mathlinkCalcNormalizedEffC@" -> mathlinkCalcNormalizedEffC,
                            "@loadCalculateUnitarityFunction@" -> loadCalculateUnitarityFunction,
                            "@calculateUnitarityMessages@" -> calculateUnitarityMessages,
                            "@setUnitarity@" -> setUnitarity,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WritePlotScripts[files_List] :=
    Module[{},
           WriteOut`ReplaceInFiles[files,
                          { Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteSLHAInputFile[inputParameters_List, files_List] :=
    Module[{formattedSLHAInputBlocks},
           formattedSLHAInputBlocks = CreateFormattedSLHABlocks[inputParameters];
           WriteOut`ReplaceInFiles[files,
                          { "@formattedSLHAInputBlocks@" -> formattedSLHAInputBlocks,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteMakefileModule[rgeFile_List, files_List] :=
    Module[{concatenatedFileList},
           concatenatedFileList = "\t" <> Utils`StringJoinWithSeparator[rgeFile, " \\\n\t"];
           WriteOut`ReplaceInFiles[files,
                          { "@generatedBetaFunctionModules@" -> concatenatedFileList,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteBVPSolverMakefile[files_List] :=
    Module[{twoScaleSource = "", twoScaleHeader = ""},
           If[FlexibleSUSY`FlexibleEFTHiggs === True,
              twoScaleSource = "\t\t" <> FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_two_scale_matching_interface.cpp"}];
              twoScaleHeader = "\t\t" <> FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_two_scale_matching_interface.hpp"}];
             ];
           WriteOut`ReplaceInFiles[files,
                   { "@FlexibleEFTHiggsTwoScaleSource@" -> twoScaleSource,
                     "@FlexibleEFTHiggsTwoScaleHeader@" -> twoScaleHeader,
                     "@FlexibleEFTHiggsShootingSource@" -> "",
                     "@FlexibleEFTHiggsShootingHeader@" -> "",
                     Sequence @@ GeneralReplacementRules[]
                   } ];
          ];

WriteUtilitiesClass[massMatrices_List, betaFun_List, inputParameters_List, extraParameters_List,
                    lesHouchesParameters_List, extraSLHAOutputBlocks_List,
                    decaysSLHAIncludeFiles_List, files_List] :=
    Module[{k, particles, susyParticles, smParticles,
            minpar, extpar, imminpar, imextpar, extraSLHAInputParameters,
            fillSpectrumVectorWithSusyParticles = "",
            fillSpectrumVectorWithSMParticles = "",
            particleLaTeXNames = "",
            particleNames = "", particleEnum = "", particleMassEnum, particleMultiplicity = "",
            particleMixingEnum = "", particleMixingNames = "",
            parameterNames = "", parameterEnum = "", numberOfParameters = 0,
            inputParameterEnum = "", inputParameterNames = "",
            extraParameterEnum = "", extraParameterNames = "",
            isLowEnergyModel = "false",
            isSupersymmetricModel = "false",
            isFlexibleEFTHiggs = "false",
            getPDGCodeFromParticleEnumNoIndex = "", getPDGCodeFromParticleEnumIndex = "",
            setParticleMultipletNameAndIndexFromPDG = "",
            fillInputParametersFromMINPAR = "", fillInputParametersFromEXTPAR = "",
            fillInputParametersFromIMMINPAR = "",
            fillInputParametersFromIMEXTPAR = "",
            writeSLHAMassBlock = "", writeSLHAMixingMatricesBlocks = "",
            writeSLHAModelParametersBlocks = "", writeSLHAPhasesBlocks = "",
            writeSLHAMinparBlock = "", writeSLHAExtparBlock = "",
            writeSLHAImMinparBlock = "", writeSLHAImExtparBlock = "",
            writeSLHAInputParameterBlocks = "",
            readLesHouchesInputParameters, writeExtraSLHAOutputBlock = "",
            readLesHouchesOutputParameters, readLesHouchesPhysicalParameters,
            gaugeCouplingNormalizationDecls = "",
            gaugeCouplingNormalizationDefs = "",
            numberOfDRbarBlocks, drBarBlockNames,
            setDecaysPrototypes = "", setDecaysFunctions = "",
            fillDecaysDataPrototypes = "", fillDecaysDataFunctions = "",
            decaysHeaderIncludes = "", useDecaysData = "",
            unitarityIncludes = ""
           },
           particles = DeleteDuplicates @ Flatten[TreeMasses`GetMassEigenstate /@ massMatrices];
           susyParticles = Select[particles, (!TreeMasses`IsSMParticle[#])&];
           smParticles   = Complement[particles, susyParticles];
           minpar = Cases[inputParameters, {p_, {"MINPAR", idx_}, ___} :> {idx, p}];
           extpar = Cases[inputParameters, {p_, {"EXTPAR", idx_}, ___} :> {idx, p}];
           imminpar = Cases[inputParameters, {p_, {"IMMINPAR", idx_}, ___} :> {idx, p}];
           imextpar = Cases[inputParameters, {p_, {"IMEXTPAR", idx_}, ___} :> {idx, p}];
           extraSLHAInputParameters = Complement[
               inputParameters,
               Cases[inputParameters, {_, {"MINPAR", _}, ___}],
               Cases[inputParameters, {_, {"EXTPAR", _}, ___}],
               Cases[inputParameters, {_, {"IMMINPAR", _}, ___}],
               Cases[inputParameters, {_, {"IMEXTPAR", _}, ___}]
           ];
           particleEnum       = TreeMasses`CreateParticleEnum[particles];
           particleMassEnum   = TreeMasses`CreateParticleMassEnum[particles];
           particleMixingEnum = TreeMasses`CreateParticleMixingEnum[massMatrices];
           particleMultiplicity = TreeMasses`CreateParticleMultiplicity[particles];
           particleNames      = TreeMasses`CreateParticleNames[particles];
           particleLaTeXNames = TreeMasses`CreateParticleLaTeXNames[particles];
           particleMixingNames= TreeMasses`CreateParticleMixingNames[massMatrices];
           inputParameterEnum  = Parameters`CreateInputParameterEnum[inputParameters];
           inputParameterNames = Parameters`CreateInputParameterNames[inputParameters];
           extraParameterEnum  = Parameters`CreateExtraParameterEnum[extraParameters];
           extraParameterNames = Parameters`CreateExtraParameterNames[extraParameters];
           fillSpectrumVectorWithSusyParticles = TreeMasses`FillSpectrumVector[susyParticles];
           fillSpectrumVectorWithSMParticles   = TreeMasses`FillSpectrumVector[smParticles];
           numberOfParameters = BetaFunction`CountNumberOfParameters[betaFun];
           parameterEnum      = BetaFunction`CreateParameterEnum[betaFun];
           parameterNames     = BetaFunction`CreateParameterNames[betaFun];
           isLowEnergyModel = If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True, "true", "false"];
           isSupersymmetricModel = If[SARAH`SupersymmetricModel === True, "true", "false"];
           isFlexibleEFTHiggs = If[FlexibleSUSY`FlexibleEFTHiggs === True, "true", "false"];
           fillInputParametersFromMINPAR = Parameters`FillInputParametersFromTuples[minpar, "MINPAR"];
           fillInputParametersFromEXTPAR = Parameters`FillInputParametersFromTuples[extpar, "EXTPAR"];
           fillInputParametersFromIMMINPAR = Parameters`FillInputParametersFromTuples[imminpar, "IMMINPAR"];
           fillInputParametersFromIMEXTPAR = Parameters`FillInputParametersFromTuples[imextpar, "IMEXTPAR"];
           readLesHouchesInputParameters = WriteOut`ReadLesHouchesInputParameters[{First[#], #[[2]]}& /@ extraSLHAInputParameters];
           readLesHouchesOutputParameters = WriteOut`ReadLesHouchesOutputParameters[lesHouchesParameters];
           readLesHouchesPhysicalParameters = WriteOut`ReadLesHouchesPhysicalParameters[lesHouchesParameters, "LOCALPHYSICAL",
                                                                                        "DEFINE_PHYSICAL_PARAMETER"];
           writeSLHAMassBlock = WriteOut`WriteSLHAMassBlock[massMatrices];
           writeSLHAMixingMatricesBlocks  = WriteOut`WriteSLHAMixingMatricesBlocks[lesHouchesParameters];
           writeSLHAModelParametersBlocks = WriteOut`WriteSLHAModelParametersBlocks[lesHouchesParameters];
           writeSLHAPhasesBlocks = WriteOut`WriteSLHAPhasesBlocks[lesHouchesParameters];
           writeSLHAMinparBlock = WriteOut`WriteSLHAMinparBlock[minpar];
           writeSLHAExtparBlock = WriteOut`WriteSLHAExtparBlock[extpar];
           writeSLHAImMinparBlock = WriteOut`WriteSLHAImMinparBlock[imminpar];
           writeSLHAImExtparBlock = WriteOut`WriteSLHAImExtparBlock[imextpar];
           writeSLHAInputParameterBlocks = WriteSLHAInputParameterBlocks[extraSLHAInputParameters];
           writeExtraSLHAOutputBlock = WriteOut`WriteExtraSLHAOutputBlock[extraSLHAOutputBlocks];
           numberOfDRbarBlocks  = WriteOut`GetNumberOfDRbarBlocks[lesHouchesParameters];
           drBarBlockNames      = WriteOut`GetDRbarBlockNames[lesHouchesParameters];
           gaugeCouplingNormalizationDecls = WriteOut`GetGaugeCouplingNormalizationsDecls[SARAH`Gauge];
           gaugeCouplingNormalizationDefs  = WriteOut`GetGaugeCouplingNormalizationsDefs[SARAH`Gauge];
           If[FSUnitarityConstraints,
              unitarityIncludes = "#include \"" <> FlexibleSUSY`FSModelName <> "_unitarity.hpp\""
           ];
           If[FlexibleSUSY`FSCalculateDecays,
              setDecaysPrototypes = WriteOut`CreateSetDecaysPrototypes[FlexibleSUSY`FSModelName];
              setDecaysFunctions = WriteOut`CreateSetDecaysFunctions[FlexibleSUSY`FSModelName];
              fillDecaysDataPrototypes = WriteOut`CreateFillDecaysDataPrototypes[FlexibleSUSY`FSModelName];
              fillDecaysDataFunctions = WriteOut`CreateFillDecaysDataFunctions[FlexibleSUSY`FSModelName];
              decaysHeaderIncludes = Utils`StringJoinWithSeparator[("#include \"decays/" <> # <> "\"")& /@ decaysSLHAIncludeFiles, "\n"];
              useDecaysData = "fill_decays_data(*decays, flexibledecay_settings);"
             ];
           getPDGCodeFromParticleEnumNoIndex = Parameters`CreatePDGCodeFromParticleCases[particles];
           getPDGCodeFromParticleEnumIndex = Parameters`CreatePDGCodeFromParticleIndexedCases[particles];
           setParticleMultipletNameAndIndexFromPDG = Parameters`CreateParticleMultipletNameAndIndexFromPDGCases[DeleteDuplicates[Join[particles, SARAH`AntiField /@ particles]]];
           WriteOut`ReplaceInFiles[files,
                          { "@fillSpectrumVectorWithSusyParticles@" -> IndentText[fillSpectrumVectorWithSusyParticles],
                            "@fillSpectrumVectorWithSMParticles@"   -> IndentText[IndentText[fillSpectrumVectorWithSMParticles]],
                            "@particleEnum@"       -> IndentText[WrapLines[particleEnum]],
                            "@particleMassEnum@"   -> IndentText[WrapLines[particleMassEnum]],
                            "@particleMixingEnum@" -> IndentText[WrapLines[particleMixingEnum]],
                            "@particleMultiplicity@" -> IndentText[WrapLines[particleMultiplicity]],
                            "@particleNames@"      -> IndentText[WrapLines[particleNames]],
                            "@particleLaTeXNames@" -> IndentText[WrapLines[particleLaTeXNames]],
                            "@parameterEnum@"     -> IndentText[WrapLines[parameterEnum]],
                            "@parameterNames@"     -> IndentText[WrapLines[parameterNames]],
                            "@particleMixingNames@"-> IndentText[WrapLines[particleMixingNames]],
                            "@inputParameterEnum@" -> IndentText[WrapLines[inputParameterEnum]],
                            "@inputParameterNames@"-> IndentText[WrapLines[inputParameterNames]],
                            "@extraParameterEnum@" -> IndentText[WrapLines[extraParameterEnum]],
                            "@extraParameterNames@"-> IndentText[WrapLines[extraParameterNames]],
                            "@isLowEnergyModel@"   -> isLowEnergyModel,
                            "@isSupersymmetricModel@" -> isSupersymmetricModel,
                            "@isFlexibleEFTHiggs@" -> isFlexibleEFTHiggs,
                            "@fillInputParametersFromMINPAR@" -> IndentText[fillInputParametersFromMINPAR],
                            "@fillInputParametersFromEXTPAR@" -> IndentText[fillInputParametersFromEXTPAR],
                            "@fillInputParametersFromIMMINPAR@" -> IndentText[fillInputParametersFromIMMINPAR],
                            "@fillInputParametersFromIMEXTPAR@" -> IndentText[fillInputParametersFromIMEXTPAR],
                            "@readLesHouchesInputParameters@" -> IndentText[readLesHouchesInputParameters],
                            "@readLesHouchesOutputParameters@" -> IndentText[readLesHouchesOutputParameters],
                            "@readLesHouchesPhysicalParameters@" -> IndentText[readLesHouchesPhysicalParameters],
                            "@writeSLHAMassBlock@" -> IndentText[writeSLHAMassBlock],
                            "@writeSLHAMixingMatricesBlocks@"  -> IndentText[writeSLHAMixingMatricesBlocks],
                            "@writeSLHAModelParametersBlocks@" -> IndentText[writeSLHAModelParametersBlocks],
                            "@writeSLHAPhasesBlocks@"          -> IndentText[writeSLHAPhasesBlocks],
                            "@writeSLHAMinparBlock@"           -> IndentText[writeSLHAMinparBlock],
                            "@writeSLHAExtparBlock@"           -> IndentText[writeSLHAExtparBlock],
                            "@writeSLHAImMinparBlock@"         -> IndentText[writeSLHAImMinparBlock],
                            "@writeSLHAImExtparBlock@"         -> IndentText[writeSLHAImExtparBlock],
                            "@writeSLHAInputParameterBlocks@"  -> IndentText[writeSLHAInputParameterBlocks],
                            "@writeExtraSLHAOutputBlock@"      -> IndentText[writeExtraSLHAOutputBlock],
                            "@gaugeCouplingNormalizationDecls@"-> IndentText[gaugeCouplingNormalizationDecls],
                            "@gaugeCouplingNormalizationDefs@" -> IndentText[gaugeCouplingNormalizationDefs],
                            "@numberOfDRbarBlocks@"            -> ToString[numberOfDRbarBlocks],
                            "@decaysHeaderIncludes@"           -> decaysHeaderIncludes,
                            "@setDecaysPrototypes@"            -> IndentText[setDecaysPrototypes],
                            "@setDecaysFunctions@"             -> setDecaysFunctions,
                            "@fillDecaysDataPrototypes@"       -> IndentText[fillDecaysDataPrototypes],
                            "@fillDecaysDataFunctions@"        -> fillDecaysDataFunctions,
                            "@drBarBlockNames@"                -> WrapLines[drBarBlockNames],
                            "@getPDGCodeFromParticleEnumNoIndex@" -> IndentText[getPDGCodeFromParticleEnumNoIndex],
                            "@getPDGCodeFromParticleEnumIndex@" -> IndentText[getPDGCodeFromParticleEnumIndex],
                            "@setParticleMultipletNameAndIndexFromPDG@" -> IndentText[setParticleMultipletNameAndIndexFromPDG],
                            "@isCPViolatingHiggsSector@"       -> CreateCBoolValue @ SA`CPViolationHiggsSector,
                            "@useDecaysData@"                   -> useDecaysData,
                            "@numberOfNeutralGoldstones@"       -> IndentText["static constexpr int number_of_neutral_goldstones = " <> ToString[TreeMasses`GetDimensionStartSkippingGoldstones[TreeMasses`GetPseudoscalarHiggsBoson[]]-1] <> ";"],
                            "@numberOfChargedGoldstones@"       -> IndentText["static constexpr int number_of_charged_goldstones = " <> ToString[TreeMasses`GetDimensionStartSkippingGoldstones[TreeMasses`GetChargedHiggsBoson[]]-1] <> ";"],
                            "@unitarityIncludes@" -> unitarityIncludes,
                            Sequence @@ GeneralReplacementRules[]
                          } ];
          ];

WriteReferences[files_List] :=
    Module[{},
           WriteOut`ReplaceInFiles[files,
               { "@referencesForComponents@" -> References`CreateCitation[],
                 Sequence @@ GeneralReplacementRules[]
               } ];
          ];

WriteSMParticlesAliases[files_List] := Module[{},
   SimplifiedName[particle_ /; TreeMasses`IsSMChargedLepton[particle] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMChargedLeptons[] === 1] := "ChargedLepton";
   SimplifiedName[particle_ /; TreeMasses`GetSMElectronLeptonMultiplet[] === particle && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMChargedLeptons[] > 1] := "Electron";
   SimplifiedName[particle_ /; TreeMasses`GetSMMuonLeptonMultiplet[] === particle && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMChargedLeptons[] > 1] := "Muon";
   SimplifiedName[particle_ /; TreeMasses`GetSMTauLeptonMultiplet[] === particle && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMChargedLeptons[] > 1] := "Tauon";
   SimplifiedName[particle_ /; TreeMasses`IsSMNeutralLepton[particle] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMNeutralLeptons[] === 1] := "Neutrino";
   SimplifiedName[particle_ /; particle === TreeMasses`GetSMNeutrino1[] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMNeutralLeptons[] > 1] := "ElectronNeutrino";
   SimplifiedName[particle_ /; particle === TreeMasses`GetSMNeutrino2[] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMNeutralLeptons[] > 1] := "MuonNeutrino";
   SimplifiedName[particle_ /; particle === TreeMasses`GetSMNeutrino3[] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMNeutralLeptons[] > 1] := "TauNeutrino";
   SimplifiedName[particle_ /; TreeMasses`IsSMDownQuark[particle] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMDownQuarks[] === 1] := "DownTypeQuark";
   SimplifiedName[particle_ /; TreeMasses`IsSMDownQuark[particle] && Head[particle] === SARAH`bar] := "AntiDownQuark";
   SimplifiedName[particle_ /; TreeMasses`IsSMUpQuark[particle] && Head[particle] =!= SARAH`bar && Length@TreeMasses`GetSMUpQuarks[] === 1] := "UpTypeQuark";
   SimplifiedName[particle_ /; TreeMasses`IsSMUpQuark[particle] && Head[particle] === SARAH`bar] := "AntiUpQuark";
   SimplifiedName[particle_ /; TreeMasses`GetHiggsBoson[] =!= Null && particle === TreeMasses`GetHiggsBoson[]] := "Higgs";
   SimplifiedName[particle_ /; TreeMasses`GetPseudoscalarHiggsBoson[] =!= Null && particle === TreeMasses`GetPseudoscalarHiggsBoson[]] := "PseudoscalarHiggs";
   SimplifiedName[particle_ /; TreeMasses`GetWBoson[] =!= Null && particle === If[GetElectricCharge[TreeMasses`GetWBoson[]] < 0, TreeMasses`GetWBoson[], Susyno`LieGroups`conj[TreeMasses`GetWBoson[]]]] := "WmBoson";
   SimplifiedName[particle_ /; TreeMasses`GetWBoson[] =!= Null && particle === If[GetElectricCharge[TreeMasses`GetWBoson[]] < 0, Susyno`LieGroups`conj[TreeMasses`GetWBoson[]], TreeMasses`GetWBoson[]]] := "WpBoson";
   SimplifiedName[particle_ /; TreeMasses`GetZBoson[] =!= Null && particle === TreeMasses`GetZBoson[]] := "ZBoson";
   SimplifiedName[particle_ /; TreeMasses`GetPhoton[] =!= Null && particle === TreeMasses`GetPhoton[]] := "Photon";
   SimplifiedName[particle_ /; TreeMasses`GetGluon[] =!= Null && particle === TreeMasses`GetGluon[]] := "Gluon";
   SimplifiedName[particle_ /; TreeMasses`GetChargedHiggsBoson[] =!= Null && particle === If[GetElectricCharge[TreeMasses`GetChargedHiggsBoson[]] < 0, TreeMasses`GetChargedHiggsBoson[], Susyno`LieGroups`conj[TreeMasses`GetChargedHiggsBoson[]]]] := "Hm";
   SimplifiedName[particle_ /; TreeMasses`GetChargedHiggsBoson[] =!= Null && particle === If[GetElectricCharge[TreeMasses`GetChargedHiggsBoson[]] < 0, Susyno`LieGroups`conj[TreeMasses`GetChargedHiggsBoson[]], TreeMasses`GetChargedHiggsBoson[]]] := "Hp";
   SimplifiedName[particle_] := particle;

   CreateParticleAlias[particle_, namespace_String] :=
      "using " <> SimplifiedName[particle] <> " = " <>
      CXXDiagrams`CXXNameOfField[particle] <> ";";

   CreateParticleAliases[particles_, namespace_:""] :=
      Utils`StringJoinWithSeparator[CreateParticleAlias[#, namespace]& /@ particles, "\n"];

   CreateSMParticleAliases[namespace_:""] :=
      Module[{smParticlesToAlias},
           smParticlesToAlias = Select[Flatten[{
                                        (* neutral Higgs bosons *)
                                        TreeMasses`GetHiggsBoson[],
                                        If[GetDimensionWithoutGoldstones[TreeMasses`GetPseudoscalarHiggsBoson[]] > 0, TreeMasses`GetPseudoscalarHiggsBoson[]],
                                        (* charged Higgs bosons *)
                                        If[GetDimensionWithoutGoldstones[TreeMasses`GetChargedHiggsBoson[]] > 0,
                                           {TreeMasses`GetChargedHiggsBoson[], Susyno`LieGroups`conj[TreeMasses`GetChargedHiggsBoson[]]}
                                        ],
                                        (* W bosons *)
                                        TreeMasses`GetWBoson[], Susyno`LieGroups`conj[TreeMasses`GetWBoson[]],
                                        (* neutral gauge bosons *)
                                        TreeMasses`GetPhoton[], TreeMasses`GetZBoson[], TreeMasses`GetGluon[],
                                        (* leptons *)
                                        TreeMasses`GetSMChargedLeptons[], TreeMasses`GetSMNeutralLeptons[],
                                        If[Length@TreeMasses`GetSMUpQuarks[] === 1,
                                           TreeMasses`GetSMUpQuarks[]
                                        ],
                                        If[Length@TreeMasses`GetSMDownQuarks[] === 1,
                                           TreeMasses`GetSMDownQuarks[]
                                        ]
                                       }, 1], (# =!= Null)&];
           CreateParticleAliases[smParticlesToAlias, namespace]
      ];

      WriteOut`ReplaceInFiles[files,
              {
                 "@ModelName@"          -> FlexibleSUSY`FSModelName,
                 "@SMParticlesAliases@" -> CreateSMParticleAliases[FlexibleSUSY`FSModelName <> "_cxx_diagrams::fields"]
              }
      ];
   ];

FilesExist[fileNames_List] :=
    And @@ (FileExistsQ /@ fileNames);

LatestModificationTimeInSeconds[file_String] :=
    If[FileExistsQ[file],
       AbsoluteTime[FileDate[file, "Modification"]], 0];

LatestModificationTimeInSeconds[files_List] :=
    Max[LatestModificationTimeInSeconds /@ files];

SARAHModelFileModificationTimeInSeconds[] :=
    LatestModificationTimeInSeconds @ \
    Join[{SARAH`ModelFile},
         FileNameJoin[{$sarahCurrentModelDir, #}]& /@ {"parameters.m", "particles.m"}];

GetRGEFileNames[outputDir_String] :=
    Module[{rgeDir, fileNames},
           rgeDir = FileNameJoin[{outputDir, "RGEs"}];
           fileNames = { "BetaYijk.m", "BetaGauge.m", "BetaMuij.m",
                         "BetaTijk.m", "BetaBij.m", "BetaVEV.m" };
           If[SARAH`AddDiracGauginos === True,
              AppendTo[fileNames, "BetaDGi.m"];
             ];
           If[SARAH`SupersymmetricModel === False,
              AppendTo[fileNames, "BetaLijkl.m"];
             ];
           If[SARAH`SupersymmetricModel === True,
              fileNames = Join[fileNames,
                               { "BetaWijkl.m", "BetaQijkl.m", "BetaLSi.m",
                                 "BetaLi.m", "Betam2ij.m", "BetaMi.m" }];
             ];
           FileNameJoin[{rgeDir, #}]& /@ fileNames
          ];

GetSelfEnergyFileNames[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "SelfEnergy.m"}];

NeedToCalculateSelfEnergies[eigenstates_] :=
    NeedToUpdateTarget[
        "self-energy",
        GetSelfEnergyFileNames[$sarahCurrentOutputMainDir, eigenstates]];

GetTadpoleFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "Tadpoles1Loop.m"}];

NeedToCalculateTadpoles[eigenstates_] :=
    NeedToUpdateTarget[
        "tadpole",
        GetTadpoleFileName[$sarahCurrentOutputMainDir, eigenstates]];

GetUnrotatedParticlesFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates],
                  "One-Loop", "UnrotatedParticles.m"}];

NeedToCalculateUnrotatedParticles[eigenstates_] :=
    NeedToUpdateTarget[
        "unrotated particle",
        GetUnrotatedParticlesFileName[$sarahCurrentOutputMainDir,eigenstates]];

NeedToCalculateRGEs[] :=
    NeedToUpdateTarget["RGE", GetRGEFileNames[$sarahCurrentOutputMainDir]];

GetVertexRuleFileName[outputDir_String, eigenstates_] :=
    FileNameJoin[{outputDir, ToString[eigenstates], "Vertices",
                  "FSVertexRules.m"}];

NeedToCalculateVertices[eigenstates_] :=
    NeedToUpdateTarget[
        "vertex",
        { GetVertexRuleFileName[$sarahCurrentOutputMainDir, eigenstates] }];

NeedToUpdateTarget[name_String, targets_List] := Module[{
        targetsExist = FilesExist[targets],
        targetTimeStamp = LatestModificationTimeInSeconds[targets],
        sarahModelFileTimeStamp = SARAHModelFileModificationTimeInSeconds[],
        files = If[Length[targets] === 1, "file", "files"],
        them = If[Length[targets] === 1, "it", "them"]
    },
    If[targetsExist,
       If[sarahModelFileTimeStamp > targetTimeStamp,
          Print["SARAH model files are newer than ", name,
                " ", files, ", updating ", them, " ..."];
          True,
          Print["Found up-to-date ", name, " ", files, "."];
          False
       ],
       Print[name, " ", files, " not found, producing ", them, " ..."];
       True
    ]
];

NeedToUpdateTarget[name_String, target_] :=
    NeedToUpdateTarget[name, {target}];

PrepareFSRules[] :=
    Block[{},
          If[Head[FlexibleSUSY`FSSelfEnergyRules]   =!= List, FlexibleSUSY`FSSelfEnergyRules   = {}];
          If[Head[FlexibleSUSY`FSVertexRules]       =!= List, FlexibleSUSY`FSVertexRules       = {}];
          If[Head[FlexibleSUSY`FSBetaFunctionRules] =!= List, FlexibleSUSY`FSBetaFunctionRules = {}];
         ];

ApplyRulesAtParts[orig_List, rules_List] :=
    (#1 /. #2)& @@@ Utils`Zip[orig, PadRight[rules, Length[orig], {{}}]];

ApplyFSBetaFunctionRules[betas_List] :=
    ({ First[#],
       Sequence @@ ApplyRulesAtParts[Drop[#,1], FlexibleSUSY`FSBetaFunctionRules] }& /@ betas) //.
    {
        SARAH`Adj[0] -> 0,
        SARAH`Conj[0] -> 0,
        SARAH`Tp[0] -> 0,
        SARAH`trace[a__]  /; MemberQ[{a}, 0] -> 0,
        SARAH`MatMul[a__] /; MemberQ[{a}, 0] -> 0
    };

ApplyFSSelfEnergyRules[se_List] :=
    { First[#], Sequence @@ ApplyRulesAtParts[Drop[#,1], FlexibleSUSY`FSSelfEnergyRules] }& /@ se;

FSPrepareRGEs[loopOrder_] :=
    Module[{needToCalculateRGEs, betas},
           If[loopOrder > 0,
              needToCalculateRGEs = NeedToCalculateRGEs[];
              SARAH`CalcRGEs[ReadLists -> !needToCalculateRGEs,
                             TwoLoop -> If[loopOrder < 2, False, True],
                             NoMatrixMultiplication -> False];
              ,
              (* create Beta* symbols with beta functions set to 0 *)
              SARAH`MakeDummyListRGEs[];
             ];
           (* check if the beta functions were calculated correctly *)
           betas = { SARAH`BetaWijkl, SARAH`BetaYijk, SARAH`BetaMuij,
                     SARAH`BetaLi, SARAH`BetaGauge, SARAH`BetaVEV,
                     SARAH`BetaQijkl, SARAH`BetaTijk, SARAH`BetaBij,
                     SARAH`BetaLSi, SARAH`Betam2ij, SARAH`BetaMi,
                     SARAH`BetaDGi, SARAH`BetaLijkl };
           If[Head[#] === Symbol && !ValueQ[#], Set[#,{}]]& /@ betas;
           If[!ValueQ[SARAH`Gij] || Head[SARAH`Gij] =!= List,
              SARAH`Gij = {};
             ];
          ];

FSCheckLoopCorrections[eigenstates_] :=
    Module[{needToCalculateLoopCorrections},
           needToCalculateLoopCorrections = Or[
               NeedToCalculateSelfEnergies[eigenstates],
               NeedToCalculateTadpoles[eigenstates],
               NeedToCalculateUnrotatedParticles[eigenstates]
                                              ];
           If[needToCalculateLoopCorrections,
              SARAH`CalcLoopCorrections[eigenstates];
             ];
          ];

FSCheckFlags[] :=
    Module[{},
           If[FlexibleSUSY`UseHiggs3LoopMSSM === True,
              SARAH`UseHiggs2LoopMSSM = True;
              FlexibleSUSY`UseMSSMYukawa2Loop = True;
              FlexibleSUSY`UseMSSMAlphaS2Loop = True;
              FlexibleSUSY`UseMSSM3LoopRGEs = True;
             ];

           If[FlexibleSUSY`UseHiggs3LoopNMSSM === True,
              FlexibleSUSY`UseHiggs2LoopNMSSM === True;
              FlexibleSUSY`UseMSSMYukawa2Loop = True;
              FlexibleSUSY`UseMSSMAlphaS2Loop = True;
              FlexibleSUSY`UseMSSM3LoopRGEs = True;
             ];

           If[FlexibleSUSY`UseHiggs3LoopSM === True,
              FlexibleSUSY`UseHiggs2LoopSM = True;
              FlexibleSUSY`UseSMAlphaS3Loop = True;
              (* FlexibleSUSY`UseSMYukawa2Loop = True; *)
              FlexibleSUSY`UseYukawa3LoopQCD = True;
              FlexibleSUSY`UseSM3LoopRGEs = True;
             ];

           If[FlexibleSUSY`UseHiggs4LoopSM === True,
              FlexibleSUSY`UseHiggs2LoopSM = True;
              FlexibleSUSY`UseHiggs3LoopSM = True;
              FlexibleSUSY`UseSMAlphaS3Loop = True;
              FlexibleSUSY`UseSMAlphaS4Loop = True;
              (* FlexibleSUSY`UseSMYukawa2Loop = True; *)
              FlexibleSUSY`UseYukawa3LoopQCD = True;
              FlexibleSUSY`UseYukawa4LoopQCD = True;
              FlexibleSUSY`UseSM3LoopRGEs = True;
              FlexibleSUSY`UseSM4LoopRGEs = True;
              FlexibleSUSY`UseSM5LoopRGEs = True;
             ];

           If[FlexibleSUSY`UseYukawa4LoopQCD === True,
              FlexibleSUSY`UseYukawa3LoopQCD = True;
              FlexibleSUSY`UseSMAlphaS3Loop = True;
              FlexibleSUSY`UseSMAlphaS4Loop = True;
              FlexibleSUSY`UseSM3LoopRGEs = True;
              FlexibleSUSY`UseSM4LoopRGEs = True;
             ];

           If[FlexibleSUSY`FlexibleEFTHiggs,
              References`AddReference["Athron:2016fuq"];
             ];

           If[FlexibleSUSY`FSCalculateDecays,
              References`AddReference["Athron:2021kve"];
              References`AddReference["Sjodahl:2012nk"];
             ];

           If[FlexibleSUSY`UseYukawa3LoopQCD || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 3-loop SM QCD corrections to yt from ",
                    "[arxiv:hep-ph/9911434, arxiv:hep-ph/9912391]"];
              References`AddReference["Chetyrkin:1999qi"];
              References`AddReference["Melnikov:2000qh"];
             ];

           If[FlexibleSUSY`UseSMYukawa2Loop || FlexibleSUSY`UseYukawa4LoopQCD ||
              FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 2-loop SM O(as^2,at*as,at^2) corrections to yt from ",
                    "[arXiv:1604.01134]"];
              Print["Adding 4-loop SM QCD corrections to yt from ",
                    "[arxiv:1502.01030, arxiv:1604.01134, arxiv:1606.06754]"];
              References`AddReference["Martin:2016xsp"];
             ];

           If[FlexibleSUSY`UseSMAlphaS3Loop || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 3-loop SM QCD threshold corrections to alpha_s ",
                    "[arxiv:hep-ph/9708255]"];
              References`AddReference["Chetyrkin:1997un"];
             ];

           If[FlexibleSUSY`UseSMAlphaS4Loop || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 4-loop SM QCD threshold corrections to alpha_s ",
                    "[arxiv:hep-ph/0512060]"];
              References`AddReference["Chetyrkin:2005ia"];
             ];

           If[FlexibleSUSY`UseMSSMYukawa2Loop,
              Print["Adding 2-loop MSSM SQCD threshold corrections for yt and yb from ",
                    "[arxiv:hep-ph/0210258, arxiv:hep-ph/0507139, arxiv:hep-ph/0707.0650]"];
              References`AddReference["Bednyakov:2002sf"];
              References`AddReference["Bednyakov:2005kt"];
              References`AddReference["Bednyakov:2007vm"];
             ];

           If[FlexibleSUSY`UseMSSMAlphaS2Loop,
              Print["Adding 2-loop MSSM SQCD threshold corrections for \[Alpha]_s from ",
                    "[arxiv:hep-ph/0509048, arxiv:0810.5101, arxiv:1009.5455]"];
              References`AddReference["Harlander:2005wm"];
              References`AddReference["Bauer:2008bj"];
              References`AddReference["Bednyakov:2010ni"];
             ];

           If[FlexibleSUSY`UseHiggs2LoopSM || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 2-loop SM Higgs mass contributions from ",
                    "[arxiv:1205.6497, arxiv:1407.4336]"];
              References`AddReference["Degrassi:2012ry"];
              References`AddReference["Martin:2014cxa"];
             ];

           If[SARAH`UseHiggs2LoopMSSM,
              Print["Adding 2-loop MSSM Higgs mass contributions from ",
                    "[arxiv:hep-ph/0105096, arxiv:hep-ph/0112177, arxiv:hep-ph/0212132,",
                    " arxiv:hep-ph/0206101, arxiv:hep-ph/0305127]"];
              References`AddReference["Degrassi:2001yf"];
              References`AddReference["Brignole:2001jy"];
              References`AddReference["Dedes:2002dy"];
              References`AddReference["Brignole:2002bz"];
              References`AddReference["Dedes:2003km"];
             ];

           If[FlexibleSUSY`UseHiggs2LoopNMSSM,
              Print["Adding 2-loop NMSSM Higgs mass contributions from ",
                    "[arxiv:hep-ph/0105096, arxiv:hep-ph/0112177, arxiv:hep-ph/0212132,",
                    " arxiv:hep-ph/0206101, arxiv:hep-ph/0305127, arxiv:0907.4682]"];
              References`AddReference["Degrassi:2001yf"];
              References`AddReference["Brignole:2001jy"];
              References`AddReference["Dedes:2002dy"];
              References`AddReference["Brignole:2002bz"];
              References`AddReference["Dedes:2003km"];
              References`AddReference["Degrassi:2009yq"];
             ];

           If[FlexibleSUSY`UseHiggs3LoopSM || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 3-loop SM Higgs mass contributions from ",
                    "[arxiv:1407.4336]"];
              References`AddReference["Martin:2014cxa"];
             ];

           If[FlexibleSUSY`UseHiggs4LoopSM || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 4-loop SM Higgs mass contributions from ",
                    "[arxiv:1508.00912]"];
              References`AddReference["Martin:2015eia"];
             ];

           If[FlexibleSUSY`UseHiggs3LoopSplit,
              Print["Adding 3-loop split-MSSM Higgs mass contributions from ",
                    "[arxiv:1312.5220]"];
              References`AddReference["Benakli:2013msa"];
             ];

           If[FlexibleSUSY`UseHiggs3LoopMSSM || FlexibleSUSY`UseHiggs3LoopNMSSM,
              Print["Adding 3-loop MSSM Higgs mass contributions from ",
                    "[arxiv:hep-ph/0803.0672, arxiv:hep-ph/1005.5709,",
                    " arxiv:1409.2297, arxiv:1708.05720]"];
              References`AddReference["Harlander:2008ju"];
              References`AddReference["Kant:2010tf"];
              References`AddReference["Kunz:2014gya"];
              References`AddReference["Harlander:2017kuc"];
             ];

           If[FlexibleSUSY`UseSM3LoopRGEs || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 3-loop SM beta-functions from ",
                    "[arxiv:1303.4364, arxiv:1307.3536,",
                    " arxiv:1504.05200 (SUSYHD v1.0.1)]"];
              References`AddReference["Bednyakov:2013eba"];
              References`AddReference["Buttazzo:2013uya"];
              References`AddReference["Vega:2015fna"];
             ];

           If[FlexibleSUSY`UseSM4LoopRGEs || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 4-loop SM beta-function from ",
                    "[arxiv:1508.00912, arXiv:1604.00853, arxiv:1508.02680]"];
              References`AddReference["Martin:2015eia"];
              References`AddReference["Chetyrkin:2016ruf"];
              References`AddReference["Bednyakov:2015ooa"];
             ];

           If[FlexibleSUSY`UseSM5LoopRGEs || FlexibleSUSY`FlexibleEFTHiggs,
              Print["Adding 5-loop SM beta-function from ",
                    "[arxiv:1606.08659]"];
              References`AddReference["Baikov:2016tgj"];
             ];

           If[FlexibleSUSY`UseMSSM3LoopRGEs,
              Print["Adding 3-loop MSSM beta-functions from ",
                    "[arxiv:hep-ph/0308231, arxiv:hep-ph/0408128]"];
              References`AddReference["Jack:2003sx"];
              References`AddReference["Jack:2004ch"];
             ];
          ];

PrepareSelfEnergies[eigenstates_] :=
    Module[{selfEnergies = {}, selfEnergiesFile},
           selfEnergiesFile = GetSelfEnergyFileNames[$sarahCurrentOutputMainDir, eigenstates];
           If[!FileExistsQ[selfEnergiesFile],
              Print["Error: self-energy files not found: ", selfEnergiesFile];
              Quit[1];
             ];
           Print["Reading self-energies from file ", selfEnergiesFile, " ..."];
           selfEnergies = Get[selfEnergiesFile];
           If[selfEnergies === Null,
              Print["Error: Could not read self-energies from ", selfEnergiesFile];
              Quit[1];
             ];
           Print["Converting self-energies ..."];
           ConvertSarahSelfEnergies[ApplyFSSelfEnergyRules @ selfEnergies]
          ];

PrepareTadpoles[eigenstates_] :=
    Module[{tadpoles = {}, tadpolesFile},
           tadpolesFile = GetTadpoleFileName[$sarahCurrentOutputMainDir, eigenstates];
           If[!FilesExist[tadpolesFile],
              Print["Error: tadpole file not found: ", tadpolesFile];
              Quit[1];
             ];
           Print["Reading tadpoles from file ", tadpolesFile, " ..."];
           tadpoles = Get[tadpolesFile];
           If[tadpoles === Null,
              Print["Error: Could not read tadpoles from ", tadpolesFile];
              Quit[1];
             ];
           Print["Converting tadpoles ..."];
           ConvertSarahTadpoles[ApplyFSSelfEnergyRules @ tadpoles]
          ];

PrepareUnrotatedParticles[eigenstates_] :=
    Module[{nonMixedParticles = {}, nonMixedParticlesFile},
           nonMixedParticlesFile = GetUnrotatedParticlesFileName[$sarahCurrentOutputMainDir, eigenstates];
           If[!FilesExist[nonMixedParticlesFile],
              Print["Error: file with unrotated fields not found: ", nonMixedParticlesFile];
              Quit[1];
             ];
           Print["Reading unrotated particles from file ", nonMixedParticlesFile, " ..."];
           nonMixedParticles = Get[nonMixedParticlesFile];
           DebugPrint["unrotated particles: ", nonMixedParticles];
           TreeMasses`SetUnrotatedParticles[nonMixedParticles];
          ];

PrepareEWSBEquations[indexReplacementRules_] :=
    Module[{ewsbEquations},
           ewsbEquations = SARAH`TadpoleEquations[FSEigenstates] /.
                           Parameters`ApplyGUTNormalization[] /.
                           indexReplacementRules /.
                           SARAH`sum[idx_, start_, stop_, expr_] :> Sum[expr, {idx,start,stop}];
           If[Head[ewsbEquations] =!= List,
              Print["Error: Could not find EWSB equations for eigenstates ",
                    FSEigenstates];
              Quit[1];
             ];
           (* filter out trivial EWSB eqs. *)
           ewsbEquations = Select[ewsbEquations, (#=!=0)&];
           ewsbEquations = Parameters`ExpandExpressions[ewsbEquations];
           (* add tadpoles to the EWSB eqs. *)
           MapIndexed[#1 - tadpole[First[#2]]&, ewsbEquations]
          ];


AddEWSBSubstitutionsForSolver[solver_, currentSubs_, extraSubs_] :=
    Module[{pos, oldSubs, newSubs},
           pos = Position[currentSubs, solver -> subs_];
           If[pos === {},
              newSubs = Append[currentSubs, Rule[solver, extraSubs]];,
              oldSubs = Extract[currentSubs, pos][[1,-1]];
              newSubs = ReplacePart[currentSubs, pos -> Rule[solver, Join[oldSubs, extraSubs]]];
             ];
           newSubs
          ];

SolveEWSBEquations[ewsbEquations_, ewsbOutputParameters_, ewsbSubstitutions_, treeLevelSolution_, treeLevelEwsbSolutionOutputFile_] :=
    Module[{i, independentEwsbEquations, ewsbSolution, freePhases},
           Print["Searching for independent EWSB equations ..."];
           independentEwsbEquations = EWSB`GetLinearlyIndependentEqs[ewsbEquations, ewsbOutputParameters,
                                                                     ewsbSubstitutions];
           If[treeLevelSolution === {},
              (* trying to find an analytic solution for the EWSB eqs. *)
              Print["Solving ", Length[independentEwsbEquations],
                    " independent EWSB equations for ",
                    ewsbOutputParameters, " ..."];
              {ewsbSolution, freePhases} = EWSB`FindSolutionAndFreePhases[independentEwsbEquations,
                                                                          ewsbOutputParameters,
                                                                          treeLevelEwsbSolutionOutputFile,
                                                                          ewsbSubstitutions];
              Print["   The EWSB solution was written to the file:"];
              Print["      ", treeLevelEwsbSolutionOutputFile];
             ,
              If[Length[treeLevelSolution] != Length[independentEwsbEquations],
                 Print["Error: not enough EWSB solutions given!"];
                 Print["   You provided solutions for ", Length[treeLevelSolution],
                       " parameters."];
                 Print["   However, there are ", Length[independentEwsbEquations],
                       " independent EWSB eqs."];
                 Quit[1];
                ];
              If[Sort[#[[1]]& /@ treeLevelSolution] =!= Sort[ewsbOutputParameters],
                 Print["Error: Parameters given in TreeLevelEWSBSolution, do not match"];
                 Print["   the Parameters given in FlexibleSUSY`EWSBOutputParameters!"];
                 Quit[1];
                ];
              Print["Using user-defined EWSB eqs. solution"];
              freePhases = {};
              ewsbSolution = Rule[#[[1]], #[[2]]]& /@ treeLevelSolution;
             ];
           {ewsbSolution, freePhases}
          ];

SolveEWSBEquationsForSolvers[solvers_List, ewsbEquations_, ewsbOutputParameters_,
                             solverSubstitutions_, treeLevelSolution_, solutionOutputFiles_] :=
    Module[{i, solver, substitutions, outputFile, solution, freePhases,
            allSolutions = {}, allFreePhases = {}},
           For[i = 1, i <= Length[solvers], i++,
               solver = solvers[[i]];
               substitutions = solver /. solverSubstitutions;
               outputFile = solver /. solutionOutputFiles;
               {solution, freePhases} = SolveEWSBEquations[ewsbEquations, ewsbOutputParameters,
                                                           substitutions, treeLevelSolution, outputFile];
               allSolutions = Append[allSolutions, solver -> solution];
               allFreePhases = Append[allFreePhases, solver -> freePhases];
              ];
           {allSolutions, allFreePhases}
          ];

SelectValidEWSBSolvers[solverSolutions_, ewsbSolvers_] :=
    Module[{i, solver, solution, validSolvers, solverEwsbSolvers = {}},
           For[i = 1, i <= Length[solverSolutions], i++,
               solver = First[solverSolutions[[i]]];
               solution = Last[solverSolutions[[i]]];
               validSolvers = ewsbSolvers;
               If[solution === {},
                  (* Fixed-point iteration can only be used if an analytic EWSB solution exists *)
                  If[MemberQ[validSolvers, FlexibleSUSY`FPIRelative],
                     Utils`FSFancyWarning[
                        "FPIRelative was selected, but no analytic solution",
                        " to the EWSB eqs. is provided. FPIRelative will be",
                        " removed from the list of EWSB solvers."
                     ];
                     validSolvers = Cases[validSolvers, Except[FlexibleSUSY`FPIRelative]];
                    ];
                  If[MemberQ[validSolvers, FlexibleSUSY`FPIAbsolute],
                     Utils`FSFancyWarning[
                        "FPIAbsolute was selected, but no analytic solution",
                        " to the EWSB eqs. is provided. FPIAbsolute will be",
                        " removed from the list of EWSB solvers."
                     ];
                     validSolvers = Cases[validSolvers, Except[FlexibleSUSY`FPIAbsolute]];
                    ];
                 ];
               solverEwsbSolvers = Append[solverEwsbSolvers, solver -> validSolvers];
              ];
           solverEwsbSolvers
          ];

GetAllFreePhases[solverFreePhases_List] := DeleteDuplicates[Flatten[#[[2]]& /@ solverFreePhases]];

ReadPoleMassPrecisions[defaultPrecision_Symbol, highPrecisionList_List,
                       mediumPrecisionList_List, lowPrecisionList_List, eigenstates_] :=
    Module[{particles, particle, i, precisionList = {}, higgs},
           If[!MemberQ[{LowPrecision, MediumPrecision, HighPrecision}, defaultPrecision],
              Print["Error: ", defaultPrecision, " is not a valid",
                    " diagonalization precision!"];
              Print["   Available are: LowPrecision, MediumPrecision, HighPrecision"];
              Quit[1];
             ];
           particles = LoopMasses`GetLoopCorrectedParticles[eigenstates];
           For[i = 1, i <= Length[particles], i++,
               particle = particles[[i]];
               Which[MemberQ[highPrecisionList  , particle], AppendTo[precisionList, {particle, HighPrecision}],
                     MemberQ[mediumPrecisionList, particle], AppendTo[precisionList, {particle, MediumPrecision}],
                     MemberQ[lowPrecisionList   , particle], AppendTo[precisionList, {particle, LowPrecision}],
                     True, AppendTo[precisionList, {particle, defaultPrecision}]
                    ];
              ];
           higgs = Cases[precisionList, {SARAH`HiggsBoson | SARAH`PseudoScalar | SARAH`ChargedHiggs, LowPrecision}];
           Message[ReadPoleMassPrecisions::ImpreciseHiggs, #[[1]], #[[2]]]& /@ higgs;
           Return[precisionList];
          ];

LoadModelFile[file_String] :=
    Module[{},
           Utils`PrintHeadline["Loading FlexibleSUSY model file"];
           If[FileExistsQ[file],
              Get[file];
              CheckModelFileSettings[];
              ,
              Print["Error: model file not found: ", file];
              Quit[1];
             ];
          ];

StripSARAHIndices[expr_, numToStrip_Integer:4] :=
    Module[{i, rules, result = expr},
           rules = Table[Parameters`StripSARAHIndicesRules[i], {i, 1, numToStrip}];
           For[i = 1, i <= numToStrip, i++,
               result = result /. rules[[i]];
              ];
           result
          ];

EnsureSMGaugeCouplingsSet[] :=
    Block[{},
          If[ValueQ[SARAH`hyperchargeCoupling] &&
             !Constraint`IsFixed[SARAH`hyperchargeCoupling,
                                 Join[FlexibleSUSY`LowScaleInput,
                                      FlexibleSUSY`SUSYScaleInput,
                                      FlexibleSUSY`HighScaleInput]],
             AppendTo[FlexibleSUSY`LowScaleInput,
                      {SARAH`hyperchargeCoupling, "new_g1"}];
            ];
          If[ValueQ[SARAH`leftCoupling] &&
             !Constraint`IsFixed[SARAH`leftCoupling,
                                 Join[FlexibleSUSY`LowScaleInput,
                                      FlexibleSUSY`SUSYScaleInput,
                                      FlexibleSUSY`HighScaleInput]],
             AppendTo[FlexibleSUSY`LowScaleInput,
                      {SARAH`leftCoupling, "new_g2"}];
            ];
          If[ValueQ[SARAH`strongCoupling] &&
             !Constraint`IsFixed[SARAH`strongCoupling,
                                 Join[FlexibleSUSY`LowScaleInput,
                                      FlexibleSUSY`SUSYScaleInput,
                                      FlexibleSUSY`HighScaleInput]],
             AppendTo[FlexibleSUSY`LowScaleInput,
                      {SARAH`strongCoupling, "new_g3"}];
            ];
         ];

EnsureEWSBConstraintApplied[] :=
    Block[{},
          If[FlexibleSUSY`FlexibleEFTHiggs === True,
             If[FreeQ[Join[FlexibleSUSY`SUSYScaleInput, FlexibleSUSY`HighScaleInput],
                      FlexibleSUSY`FSSolveEWSBFor[___]],
                AppendTo[FlexibleSUSY`SUSYScaleInput,
                         FlexibleSUSY`FSSolveEWSBFor[FlexibleSUSY`EWSBOutputParameters]];
               ];
             ,
             If[FreeQ[Join[FlexibleSUSY`LowScaleInput, FlexibleSUSY`SUSYScaleInput, FlexibleSUSY`HighScaleInput],
                      FlexibleSUSY`FSSolveEWSBFor[___]],
                AppendTo[FlexibleSUSY`SUSYScaleInput,
                         FlexibleSUSY`FSSolveEWSBFor[FlexibleSUSY`EWSBOutputParameters]];
               ];
            ];
         ];

EnforceSLHA1Compliance[{parameter_, properties_List}] :=
    Module[{dims},
           dims = Select[properties, (First[#] === Parameters`ParameterDimensions)&];
           dims = Select[dims, And[Last[#] =!= {}, Last[#] =!= {0}, Last[#] =!= {1}]&];
           If[dims =!= {},
              Print["Error: the SLHA1 input parameter ", parameter, " must be a scalar!"];
              Print["   Please define ", parameter, " in a different block,"];
              Print["   or define it to be a scalar."];
              Quit[1];
             ];
          ];

AddSLHA1InputParameterInfo[parameter_, blockName_String, blockEntry_] :=
    Module[{i, definedPars, defaultInfo, oldInfo, pos, property},
           definedPars = First /@ FlexibleSUSY`FSAuxiliaryParameterInfo;
           defaultInfo = {parameter, { Parameters`ParameterDimensions -> {1},
                                       SARAH`LesHouches -> {blockName, blockEntry},
                                       Parameters`InputParameter -> True } };
           If[!MemberQ[definedPars, parameter],
              PrependTo[FlexibleSUSY`FSAuxiliaryParameterInfo, defaultInfo];,
              pos = Position[FlexibleSUSY`FSAuxiliaryParameterInfo, {parameter, {__}}];
              oldInfo = Extract[FlexibleSUSY`FSAuxiliaryParameterInfo, pos];
              EnforceSLHA1Compliance /@ oldInfo;
              For[i = 1, i <= Length[defaultInfo[[2]]], i++,
                  property = defaultInfo[[2,i]];
                  oldInfo = {#[[1]], Utils`AppendOrReplaceInList[#[[2]], property, (First[#1] === First[#2])&]}& /@ oldInfo;
                 ];
              FlexibleSUSY`FSAuxiliaryParameterInfo = ReplacePart[FlexibleSUSY`FSAuxiliaryParameterInfo,
                                                                  MapThread[Rule, {pos, oldInfo}]];
             ];
          ];

AddSLHA1InputBlockInfo[blockName_String, inputParameters_List] :=
    AddSLHA1InputParameterInfo[#[[2]], blockName, #[[1]]]& /@ inputParameters;

AddUnfixedParameterInfo[par_, inputPar_, blockList_] :=
    Module[{definedPars, outputBlock, inputBlock, info},
           definedPars = First /@ FlexibleSUSY`FSAuxiliaryParameterInfo;
           If[!MemberQ[definedPars, inputPar],
              outputBlock = Parameters`FindSLHABlock[blockList, par];
              inputBlock = WriteOut`CreateInputBlockName[outputBlock];
              info = {inputPar, { Parameters`InputParameter -> True,
                                  SARAH`LesHouches -> inputBlock,
                                  Parameters`ParameterDimensions -> Parameters`GetParameterDimensions[par],
                                  Parameters`MassDimension -> Parameters`GetModelParameterMassDimension[par]
                                } };
              AppendTo[FlexibleSUSY`FSAuxiliaryParameterInfo, info];
             ];
          ];

AddUnfixedParameterBlockInfo[unfixedParameters_List, blockList_List] :=
    AddUnfixedParameterInfo[#[[1]], #[[2]], blockList]& /@ unfixedParameters;

FindFixedParameters[] :=
    If[FlexibleSUSY`FlexibleEFTHiggs === True,
       Join[Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`SUSYScaleInput],
            Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`HighScaleInput],
            Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`MatchingScaleInput],
            {SARAH`hyperchargeCoupling, SARAH`leftCoupling, SARAH`strongCoupling,
             SARAH`UpYukawa, SARAH`DownYukawa, SARAH`ElectronYukawa}]
       ,
       Join[Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`LowScaleInput],
            Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`SUSYScaleInput],
            Constraint`FindFixedParametersFromConstraint[FlexibleSUSY`HighScaleInput]]
      ];

FindUnfixedParameters[parameters_List, fixed_List] :=
    Complement[parameters, DeleteDuplicates[Flatten[fixed]]];

AddLesHouchesInputParameterInfo[par_, inputPar_, blockList_List] :=
    Module[{definedPars, outputBlock, inputBlock, info},
           definedPars = First /@ FlexibleSUSY`FSAuxiliaryParameterInfo;
           If[!MemberQ[definedPars, inputPar],
              outputBlock = Parameters`FindSLHABlock[blockList, par];
              inputBlock = WriteOut`CreateInputBlockName[outputBlock];
              info = {inputPar, { Parameters`InputParameter -> True,
                                  SARAH`LesHouches -> inputBlock,
                                  Parameters`ParameterDimensions -> Parameters`GetParameterDimensions[par],
                                  Parameters`MassDimension -> Parameters`GetModelParameterMassDimension[par]
                                } };
              AppendTo[FlexibleSUSY`FSAuxiliaryParameterInfo, info];
             ];
          ];

AddLesHouchesInputParameterBlockInfo[inputPars_List, blockList_List] :=
    AddLesHouchesInputParameterInfo[#[[1]], #[[2]], blockList]& /@ inputPars;

AppendLesHouchesInfo[lesHouchesList_List, auxiliaryInfo_List] :=
    Module[{getSLHAInfo},
           getSLHAInfo[{parameter_, properties_List}] :=
               Module[{block},
                      block = Cases[properties, (SARAH`LesHouches -> value_) :> value];
                      If[block === {},
                         block = None,
                         block = Last[block]
                        ];
                      {parameter, block}
                     ];
           Join[lesHouchesList, getSLHAInfo /@ auxiliaryInfo]
          ];

(* returns beta functions of VEV phases *)
GetVEVPhases[eigenstates_:FlexibleSUSY`FSEigenstates] :=
    Flatten @ Cases[DEFINITION[eigenstates][SARAH`VEVs], {_,_,_,_, p_} :> p];

AddSM3LoopRGE[beta_List, couplings_List] :=
    Module[{rules, MakeRule},
           MakeRule[coupling_] := {
               RuleDelayed[{coupling         , b1_, b2_},
                           {coupling         , b1 , b2, Part[ThreeLoopSM`BetaSM[coupling], 3]}],
               RuleDelayed[{coupling[i1_,i2_], b1_, b2_},
                           {coupling[i1,i2]  , b1 , b2, Part[ThreeLoopSM`BetaSM[coupling], 3] CConversion`PROJECTOR}]
           };
           rules = Flatten[MakeRule /@ couplings];
           beta /. rules
          ];

AddSM3LoopRGEs[] := Module[{
    gauge = { SARAH`hyperchargeCoupling,
              SARAH`leftCoupling,
              SARAH`strongCoupling },
    yuks  = { SARAH`UpYukawa,
              SARAH`DownYukawa,
              SARAH`ElectronYukawa },
    quart = { Parameters`GetParameterFromDescription["SM Higgs Selfcouplings"] },
    bilin = { Parameters`GetParameterFromDescription["SM Mu Parameter"] }
    },
    SARAH`BetaGauge = AddSM3LoopRGE[SARAH`BetaGauge, gauge];
    SARAH`BetaYijk  = AddSM3LoopRGE[SARAH`BetaYijk , yuks];
    SARAH`BetaLijkl = AddSM3LoopRGE[SARAH`BetaLijkl, quart];
    SARAH`BetaBij   = AddSM3LoopRGE[SARAH`BetaBij  , bilin];
    ];

AddSM4LoopRGE[beta_List, couplings_List] :=
    Module[{rules, MakeRule},
           MakeRule[coupling_] := {
               RuleDelayed[{coupling         , b1_, b2_, b3_},
                           {coupling         , b1 , b2 , b3, Part[ThreeLoopSM`BetaSM[coupling], 4]}],
               RuleDelayed[{coupling[i1_,i2_], b1_, b2_, b3_},
                           {coupling[i1,i2]  , b1 , b2 , b3, Part[ThreeLoopSM`BetaSM[coupling], 4] CConversion`PROJECTOR}]
           };
           rules = Flatten[MakeRule /@ couplings];
           beta /. rules
          ];

AddSM5LoopRGE[beta_List, couplings_List] :=
    Module[{rules, MakeRule},
           MakeRule[coupling_] := {
               RuleDelayed[{coupling         , b1_, b2_, b3_, b4_},
                           {coupling         , b1 , b2 , b3 , b4, Part[ThreeLoopSM`BetaSM[coupling], 5]}],
               RuleDelayed[{coupling[i1_,i2_], b1_, b2_, b3_, b4_},
                           {coupling[i1,i2]  , b1 , b2 , b3 , b4, Part[ThreeLoopSM`BetaSM[coupling], 5] CConversion`PROJECTOR}]
           };
           rules = Flatten[MakeRule /@ couplings];
           beta /. rules
          ];

AddSM4LoopRGEs[] := Module[{
    gauge = { SARAH`strongCoupling },
    yuks  = { SARAH`UpYukawa },
    quart = { Parameters`GetParameterFromDescription["SM Higgs Selfcouplings"] }
    },
    SARAH`BetaGauge = AddSM4LoopRGE[SARAH`BetaGauge, gauge];
    SARAH`BetaYijk  = AddSM4LoopRGE[SARAH`BetaYijk , yuks];
    SARAH`BetaLijkl = AddSM4LoopRGE[SARAH`BetaLijkl, quart];
    ];

AddSM5LoopRGEs[] := Module[{
    gauge = { SARAH`strongCoupling }
    },
    SARAH`BetaGauge = AddSM5LoopRGE[SARAH`BetaGauge, gauge];
    ];

AddMSSM3LoopRGE[beta_List, couplings_List] :=
    Module[{rules, MakeRule},
           MakeRule[coupling_] := {
               RuleDelayed[{coupling         , b1_, b2_}, {coupling       , b1, b2, Last[ThreeLoopMSSM`BetaMSSM[coupling]]}],
               RuleDelayed[{coupling[i1_,i2_], b1_, b2_}, {coupling[i1,i2], b1, b2, Last[ThreeLoopMSSM`BetaMSSM[coupling]]}]
           };
           rules = Flatten[MakeRule /@ couplings];
           beta /. rules
          ];

AddMSSM3LoopRGEs[] := Module[{
    gauge = { SARAH`hyperchargeCoupling,
              SARAH`leftCoupling,
              SARAH`strongCoupling },
    yuks  = { SARAH`UpYukawa,
              SARAH`DownYukawa,
              SARAH`ElectronYukawa },
    gaugi = { Parameters`GetParameterFromDescription["Bino Mass parameter"],
              Parameters`GetParameterFromDescription["Wino Mass parameter"],
              Parameters`GetParameterFromDescription["Gluino Mass parameter"] },
    trili = { SARAH`TrilinearUp, SARAH`TrilinearDown, SARAH`TrilinearLepton },
    mass2 = { SARAH`SoftSquark, SARAH`SoftUp, SARAH`SoftDown,
              SARAH`SoftLeftLepton, SARAH`SoftRightLepton,
              Parameters`GetParameterFromDescription["Softbreaking Down-Higgs Mass"],
              Parameters`GetParameterFromDescription["Softbreaking Up-Higgs Mass"] },
    mu    = { Parameters`GetParameterFromDescription["Mu-parameter"] },
    bmu   = { Parameters`GetParameterFromDescription["Bmu-parameter"] }
    },
    SARAH`BetaGauge = AddMSSM3LoopRGE[SARAH`BetaGauge, gauge];
    SARAH`BetaYijk  = AddMSSM3LoopRGE[SARAH`BetaYijk , yuks];
    SARAH`BetaMi    = AddMSSM3LoopRGE[SARAH`BetaMi   , gaugi];
    SARAH`BetaMuij  = AddMSSM3LoopRGE[SARAH`BetaMuij , mu   ];
    SARAH`BetaBij   = AddMSSM3LoopRGE[SARAH`BetaBij  , bmu  ];
    SARAH`BetaTijk  = AddMSSM3LoopRGE[SARAH`BetaTijk , trili];
    SARAH`Betam2ij  = AddMSSM3LoopRGE[SARAH`Betam2ij , mass2];
    ];

SelectRenormalizationScheme::UnknownRenormalizationScheme = "Unknown\
 renormalization scheme `1`.";

SelectRenormalizationScheme[renormalizationScheme_] :=
    Switch[renormalizationScheme,
           FlexibleSUSY`DRbar, 0,
           FlexibleSUSY`MSbar, 1,
           _, Message[SelectRenormalizationScheme::UnknownRenormalizationScheme, renormalizationScheme];
              Quit[1];
          ];

RenameSLHAInputParametersInUserInput[lesHouchesInputParameters_] :=
    Module[{lesHouchesInputParameterReplacementRules},
           lesHouchesInputParameterReplacementRules = Flatten[{
               Rule[SARAH`LHInput[#[[1]]], #[[2]]],
               Rule[SARAH`LHInput[#[[1]][p__]], #[[2]][p]]
           }& /@ lesHouchesInputParameters];

           FlexibleSUSY`LowScaleInput = FlexibleSUSY`LowScaleInput /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`SUSYScaleInput = FlexibleSUSY`SUSYScaleInput /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`HighScaleInput = FlexibleSUSY`HighScaleInput /.
               lesHouchesInputParameterReplacementRules;

           FlexibleSUSY`InitialGuessAtLowScale = FlexibleSUSY`InitialGuessAtLowScale /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`InitialGuessAtSUSYScale = FlexibleSUSY`InitialGuessAtSUSYScale /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`InitialGuessAtHighScale = FlexibleSUSY`InitialGuessAtHighScale /.
               lesHouchesInputParameterReplacementRules;

           FlexibleSUSY`LowScaleFirstGuess = FlexibleSUSY`LowScaleFirstGuess /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`SUSYScaleFirstGuess = FlexibleSUSY`SUSYScaleFirstGuess /.
               lesHouchesInputParameterReplacementRules;
           FlexibleSUSY`HighScaleFirstGuess = FlexibleSUSY`HighScaleFirstGuess /.
               lesHouchesInputParameterReplacementRules;
          ];

ReadSARAHBetaFunctions[] :=
    Module[{susyBetaFunctions, susyBreakingBetaFunctions},
           FSPrepareRGEs[FlexibleSUSY`FSRGELoopOrder];

           FlexibleSUSY`FSRenormalizationScheme = GetRenormalizationScheme[];

           (* adapt SARAH`Conj to our needs *)
           (* Clear[Conj]; *)
           SARAH`Conj[(B_)[b__]] =.;
           SARAH`Conj /: SARAH`Conj[SARAH`Conj[x_]] := x;
           RXi[_] = 1;
           SARAH`Xi = 1;
           SARAH`Xip = 1;
           SARAH`rMS = SelectRenormalizationScheme[FlexibleSUSY`FSRenormalizationScheme];

           If[FlexibleSUSY`UseSM3LoopRGEs,
              AddSM3LoopRGEs[];
             ];

           If[FlexibleSUSY`UseSM4LoopRGEs,
              AddSM4LoopRGEs[];
             ];

           If[FlexibleSUSY`UseSM5LoopRGEs,
              AddSM5LoopRGEs[];
             ];

           If[FlexibleSUSY`UseMSSM3LoopRGEs,
              AddMSSM3LoopRGEs[];
             ];

           If[SARAH`SupersymmetricModel,
              (* pick beta functions of supersymmetric parameters *)
              susyBetaFunctions = { SARAH`BetaLijkl,
                                    SARAH`BetaWijkl,
                                    SARAH`BetaYijk ,
                                    SARAH`BetaMuij ,
                                    SARAH`BetaLi   ,
                                    SARAH`BetaGauge,
                                    SARAH`BetaVEV  };

              (* pick beta functions of non-supersymmetric parameters *)
              susyBreakingBetaFunctions = { SARAH`BetaQijkl,
                                            SARAH`BetaTijk ,
                                            SARAH`BetaBij  ,
                                            SARAH`BetaLSi  ,
                                            SARAH`Betam2ij ,
                                            SARAH`BetaMi   ,
                                            SARAH`BetaDGi  };
              ,
              (* pick beta functions of dimensionless parameters *)
              susyBetaFunctions = { SARAH`BetaGauge,
                                    SARAH`BetaLijkl, (* quartic scalar interactions *)
                                    SARAH`BetaYijk };

              (* pick beta functions of dimensionfull parameters *)
              susyBreakingBetaFunctions = { SARAH`BetaTijk, (* cubic scalar interactions *)
                                            SARAH`BetaMuij, (* bilinear fermion term *)
                                            SARAH`BetaBij , (* bilinear scalar term *)
                                            SARAH`BetaLi  , (* linear scalar term *)
                                            SARAH`BetaVEV };
             ];

           (* filter out buggy and duplicate beta functions *)
           DeleteBuggyBetaFunctions[beta_List] :=
               DeleteDuplicates[Select[beta, (!NumericQ[#[[1]]])&], (#1[[1]] === #2[[1]])&];

           susyBetaFunctions         = DeleteBuggyBetaFunctions @ (Join @@ susyBetaFunctions);
           susyBreakingBetaFunctions = DeleteBuggyBetaFunctions @ (Join @@ susyBreakingBetaFunctions);

           {susyBetaFunctions, susyBreakingBetaFunctions}
    ];

(* disable tensor couplings *)
FSDisableTensorCouplings[parameters_] :=
    Module[{tensorCouplings = Select[parameters, Parameters`IsTensor]},
           If[tensorCouplings =!= {},
              Print["Error: Models with tensor couplings are currently not supported."];
              Print["   The following parameters have a tensor structure:"];
              Print["   ", InputForm[tensorCouplings]];
              Quit[1];
           ];
    ];

SetupModelParameters[susyBetaFunctions_, susyBreakingBetaFunctions_] :=
    Module[{allParameters, phases},
           (* identify real parameters *)
           If[Head[SARAH`RealParameters] === List,
              Parameters`AddRealParameter[SARAH`RealParameters];
             ];

           (* store all model parameters *)
           allParameters = StripSARAHIndices[((#[[1]])& /@ Join[susyBetaFunctions, susyBreakingBetaFunctions])];
           Parameters`SetModelParameters[allParameters];
           DebugPrint["Model parameters: ", allParameters];

           (* collect all phases from SARAH *)
           phases = DeleteDuplicates @ Join[
               ConvertSarahPhases[SARAH`ParticlePhases],
               Exp[I #]& /@ GetVEVPhases[FlexibleSUSY`FSEigenstates]];
           Parameters`SetPhases[phases];

           FSDisableTensorCouplings[allParameters];

           allParameters
    ];

ConvertBetaFunctions[susyBetaFunctionsSARAH_, susyBreakingBetaFunctionsSARAH_] :=
    Module[{susyBetaFunctions, susyBreakingBetaFunctions,
	    numberOfSusyParameters, numberOfSusyBreakingParameters},
	   susyBetaFunctions = BetaFunction`ConvertSarahRGEs[ApplyFSBetaFunctionRules @ susyBetaFunctionsSARAH];
           susyBetaFunctions = Select[susyBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];

           susyBreakingBetaFunctions = BetaFunction`ConvertSarahRGEs[ApplyFSBetaFunctionRules @ susyBreakingBetaFunctionsSARAH];
           susyBreakingBetaFunctions = Select[susyBreakingBetaFunctions, (BetaFunction`GetAllBetaFunctions[#]!={})&];

           allBetaFunctions = Join[susyBetaFunctions, susyBreakingBetaFunctions];

           numberOfSusyParameters = BetaFunction`CountNumberOfParameters[susyBetaFunctions];
           numberOfSusyBreakingParameters = BetaFunction`CountNumberOfParameters[susyBreakingBetaFunctions];
           numberOfModelParameters = numberOfSusyParameters + numberOfSusyBreakingParameters;

           {susyBetaFunctions, susyBreakingBetaFunctions}
    ];

SetupMassMatrices[allParameters_] :=
		Module[{Lat$massMatrices, massMatrices,
		        allIntermediateOutputParameters,
		        allIntermediateOutputParameterIndexReplacementRules},
           allIndexReplacementRules = Join[
             Parameters`CreateIndexReplacementRules[allParameters],
             {Global`upQuarksDRbar[i_,j_] :> Global`upQuarksDRbar[i-1,j-1],
             Global`downQuarksDRbar[i_,j_] :> Global`downQuarksDRbar[i-1,j-1],
             Global`downLeptonsDRbar[i_,j_] :> Global`downLeptonsDRbar[i-1,j-1]}
		       ];

		       Lat$massMatrices = TreeMasses`ConvertSarahMassMatrices[] /.
		         Parameters`ApplyGUTNormalization[] //.
		         { SARAH`sum[j_, start_, end_, expr_] :> (Sum[expr, {j,start,end}]) };

		       massMatrices = Lat$massMatrices /. allIndexReplacementRules;
		       Lat$massMatrices = LatticeUtils`FixDiagonalization[Lat$massMatrices];

		       allIntermediateOutputParameters =
		         Parameters`GetIntermediateOutputParameterDependencies[
		           TreeMasses`GetMassMatrix /@ massMatrices];
		       DebugPrint["intermediate output parameters = ", allIntermediateOutputParameters];

		       (* decrease index literals of intermediate output parameters in mass matrices *)
		       allIntermediateOutputParameterIndexReplacementRules =
		         Parameters`CreateIndexReplacementRules[allIntermediateOutputParameters];

		       massMatrices = massMatrices /. allIntermediateOutputParameterIndexReplacementRules;

		       {massMatrices, Lat$massMatrices}
		];

SetupOutputParameters[massMatrices_] :=
		Module[{allParticles, allOutputParameters},
           allParticles = FlexibleSUSY`M[TreeMasses`GetMassEigenstate[#]]& /@ massMatrices;
           allOutputParameters = DeleteCases[DeleteDuplicates[
               Join[allParticles,
                    Flatten[TreeMasses`GetMixingMatrixSymbol[#]& /@ massMatrices]]], Null];

           Parameters`SetOutputParameters[allOutputParameters];
           DebugPrint["output parameters = ", allOutputParameters];
    ];

CheckObsDependencies[requested_List] :=
Module[{allObs, dir, filtered = requested},
   Needs@"NPointFunctions`";
   If[FlexibleSUSY`FSFeynArtsAvailable && FlexibleSUSY`FSFormCalcAvailable,
      Needs@"WilsonCoeffs`";
      Return[filtered];
   ];

   allObs = Cases[
      requested,
      x_/; Context@Evaluate@Head@x === "FlexibleSUSYObservable`" :> x,
      Infinity,
      Heads -> True
   ];
   allObs = DeleteDuplicates[SymbolName@*Head/@allObs];

   Do[
      dir = FileNameJoin@{$flexiblesusyMetaDir, "Observables", obs};
      If[DirectoryQ@dir,
         If[TextSearch[dir, "NPointFunctions" | "WilsonCoeffs", "Count"] > 0,
            Utils`FSFancyWarning[obs,
               " is requested but FeynArts or FormCalc are disabled. ",
               "Removing ", obs, " from calculated observables."
            ];
            filtered = filtered /. Symbol["FlexibleSUSYObservable`"<>obs][___] :> Null;
         ];
      ];,
      {obs, allObs}
   ];

   filtered /. {_Integer, Null} :> Sequence[] /. {_Symbol, {}} :> Sequence[]
];

Options[MakeFlexibleSUSY] :=
    {
        InputFile -> "FlexibleSUSY.m",
        OutputDirectory -> "",
        DebugOutput -> False
    };

MakeFlexibleSUSY[OptionsPattern[]] :=
    Module[{nPointFunctions, initialGuesserInputFile,
            aMMVertices, edmFields, ammFields,
            QToQGammaFields = {},
            FFMasslessVVertices = {},
            ObservablesExtraOutput, observablesExtraVertices = {},
            cxxQFTTemplateDir, cxxQFTOutputDir, cxxQFTFiles,
            cxxQFTVerticesTemplate, cxxQFTVerticesMakefileTemplates,
            susyBetaFunctions, susyBreakingBetaFunctions,
            anomDim,
            inputParameters (* list of 3-component lists of the form {name, block, type} *),
            massMatrices,
            diagonalizationPrecision,
            allInputParameterIndexReplacementRules = {},
            allExtraParameterIndexReplacementRules = {},
            allParameters,
            ewsbEquations, sharedEwsbSubstitutions = {}, solverEwsbSubstitutions = {},
            freePhases = {}, solverFreePhases = {}, solverEwsbSolutions = {}, missingPhases,
            treeLevelEwsbSolutionOutputFiles = {}, treeLevelEwsbEqsOutputFile,
            solverEwsbSolvers = {}, fixedParameters,
            lesHouchesInputParameters,
            extraSLHAOutputBlocks,
            deltaVBwave, deltaVBvertex, deltaVBbox,
            vertexRules, vertexRuleFileName,
            Lat$massMatrices, spectrumGeneratorInputFile,
            semiAnalyticBCs, semiAnalyticSolns,
            semiAnalyticHighScaleFiles, semiAnalyticSUSYScaleFiles, semiAnalyticLowScaleFiles,
            semiAnalyticSolnsOutputFile, semiAnalyticEWSBSubstitutions = {}, semiAnalyticInputScale = "",
            decaysFinalStateParticles = {}, decaysVertices = {},
            decaysSources = {}, decaysHeaders = {}, decaysSLHAIncludeFiles = {}},

           Utils`PrintHeadline["Starting FlexibleSUSY"];
           FSDebugOutput["meta code directory: ", $flexiblesusyMetaDir];
           FSDebugOutput["config directory   : ", $flexiblesusyConfigDir];
           FSDebugOutput["templates directory: ", $flexiblesusyTemplateDir];

           (* check if SARAH`Start[] was called *)
           If[!ValueQ[Model`Name],
              Print["Error: Model`Name is not defined.  Did you call SARAH`Start[\"Model\"]?"];
              Quit[1];
             ];
           FSDebugOutput = OptionValue[DebugOutput];
           FSOutputDir = OptionValue[OutputDirectory];
           If[!DirectoryQ[FSOutputDir],
              Print["Error: OutputDirectory ", FSOutputDir, " does not exist."];
              Print["   Please run ./createmodel first."];
              Quit[1]];
           CheckSARAHVersion[];
           (* set default values that need SARAH`Start[] to be called first *)
           FSSetDefaultModelFileSettings[];
           (* load model file *)
           LoadModelFile[OptionValue[InputFile]];
           Print["FlexibleSUSY model file loaded"];
           Print["  Model: ", Style[FlexibleSUSY`FSModelName, FSColor]];
           Print["  Model file: ", OptionValue[InputFile]];
           Print["  Model output directory: ", FSOutputDir];

           Utils`PrintHeadline["Reading SARAH output files"];
           PrepareFSRules[];

           {susyBetaFunctionsSARAH, susyBreakingBetaFunctionsSARAH} = ReadSARAHBetaFunctions[];

           FSCheckFlags[];
           FSCheckLoopCorrections[FSEigenstates];
           nPointFunctions = Vertices`EnforceCpColorStructures @
           					 Vertices`SortCps @
             				 Join[PrepareSelfEnergies[FSEigenstates], PrepareTadpoles[FSEigenstates]];
           PrepareUnrotatedParticles[FSEigenstates];

           DebugPrint["particles (mass eigenstates): ", TreeMasses`GetParticles[]];

           allParameters = SetupModelParameters[susyBetaFunctionsSARAH, susyBreakingBetaFunctionsSARAH];

           Needs@"Observables`";
           FlexibleSUSY`ExtraSLHAOutputBlocks = CheckObsDependencies[FlexibleSUSY`ExtraSLHAOutputBlocks];

           Print["Converting SARAH beta functions ..."];
           {susyBetaFunctions, susyBreakingBetaFunctions} =
	       ConvertBetaFunctions[susyBetaFunctionsSARAH, susyBreakingBetaFunctionsSARAH];

           Print["Converting SARAH anomalous dimensions ..."];
           anomDim = AnomalousDimension`ConvertSarahAnomDim[SARAH`Gij];

           FlexibleSUSY`FSLesHouchesList = SA`LHList;

           (* collect input parameters from MINPAR and EXTPAR lists *)
           AddSLHA1InputBlockInfo["IMEXTPAR", Reverse @ IMEXTPAR];
           AddSLHA1InputBlockInfo["IMMINPAR", Reverse @ IMMINPAR];
           AddSLHA1InputBlockInfo["EXTPAR", Reverse @ SARAH`EXTPAR];
           AddSLHA1InputBlockInfo["MINPAR", Reverse @ SARAH`MINPAR];

           (* search for unfixed parameters *)
           Constraint`SanityCheck[Join[If[FlexibleEFTHiggs === True,
                                          FlexibleSUSY`InitialGuessAtSUSYScale,
                                          FlexibleSUSY`InitialGuessAtLowScale],
                                       FlexibleSUSY`InitialGuessAtHighScale],
                                  "initial guess"
                                 ];

           (* add SM gauge couplings to low-scale constraint if not set anywhere *)
           EnsureSMGaugeCouplingsSet[];

           (* add EWSB constraint to SUSY-scale constraint if not set *)
           EnsureEWSBConstraintApplied[];

           fixedParameters = FindFixedParameters[];
           FlexibleSUSY`FSUnfixedParameters = FindUnfixedParameters[allParameters, fixedParameters];

           If[FlexibleSUSY`FSUnfixedParameters =!= {} &&
              FlexibleSUSY`AutomaticInputAtMSUSY =!= True,
              Utils`FSFancyWarning[
                 "The following parameters are not fixed by any constraint: ",
                 FlexibleSUSY`FSUnfixedParameters
              ];
             ];

           (* add the unfixed parameters to the susy scale constraint *)
           If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY === True &&
              FlexibleSUSY`AutomaticInputAtMSUSY,
              (* adding input names for the parameters *)
              FlexibleSUSY`FSUnfixedParameters = Select[
                  StripSARAHIndices[
                      Join[{BetaFunction`GetName[#], Symbol[CConversion`ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"]}& /@ susyBetaFunctions,
                           {BetaFunction`GetName[#], Symbol[CConversion`ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"]}& /@ susyBreakingBetaFunctions]
                  ],
                  MemberQ[FlexibleSUSY`FSUnfixedParameters,#[[1]]]&
              ];
              FlexibleSUSY`SUSYScaleInput = Join[FlexibleSUSY`SUSYScaleInput,
                                                 {#[[1]],#[[2]]}& /@ FlexibleSUSY`FSUnfixedParameters];
              AddUnfixedParameterBlockInfo[FlexibleSUSY`FSUnfixedParameters, FlexibleSUSY`FSLesHouchesList];
             ];

           lesHouchesInputParameters = DeleteDuplicates[
               Flatten[
                   Cases[
                       Join[FlexibleSUSY`LowScaleInput,
                            FlexibleSUSY`SUSYScaleInput,
                            FlexibleSUSY`HighScaleInput,
                            FlexibleSUSY`InitialGuessAtLowScale,
                            FlexibleSUSY`InitialGuessAtHighScale,
                            {FlexibleSUSY`LowScaleFirstGuess,
                             FlexibleSUSY`SUSYScaleFirstGuess,
                             FlexibleSUSY`HighScaleFirstGuess}
                           ],
                       SARAH`LHInput[p_] :> Parameters`StripIndices[p],
                       Infinity
                        ]
                      ]
           ];

           lesHouchesInputParameters = Select[StripSARAHIndices[{BetaFunction`GetName[#],
                                                                 Symbol[CConversion`ToValidCSymbolString[BetaFunction`GetName[#]] <> "Input"]
                                                                }& /@ Join[susyBetaFunctions, susyBreakingBetaFunctions]],
                                              MemberQ[lesHouchesInputParameters,#[[1]]]&];

           AddLesHouchesInputParameterBlockInfo[lesHouchesInputParameters, FlexibleSUSY`FSLesHouchesList];

           (* apply parameter definitions and properties *)
           Parameters`ApplyAuxiliaryParameterInfo[FlexibleSUSY`FSAuxiliaryParameterInfo];
           Parameters`CheckInputParameterDefinitions[];

           FlexibleSUSY`FSLesHouchesList = AppendLesHouchesInfo[FlexibleSUSY`FSLesHouchesList, FlexibleSUSY`FSAuxiliaryParameterInfo];

           inputParameters = Parameters`GetInputParametersAndTypes[];

           DebugPrint["input parameters: ", Parameters`GetInputParameters[]];
           DebugPrint["auxiliary parameters: ", Parameters`GetExtraParameters[]];

           On[Assert];

           {massMatrices, Lat$massMatrices} = SetupMassMatrices[allParameters];

           allInputParameterIndexReplacementRules = Parameters`CreateIndexReplacementRules[
               Parameters`GetInputParameters[]
            ];

           allExtraParameterIndexReplacementRules = Parameters`CreateIndexReplacementRules[
               Parameters`GetExtraParameters[]
            ];

           SetupOutputParameters[massMatrices];

           (* backwards compatibility replacements in constraints *)
           backwardsCompatRules = {
               Global`topDRbar      -> Global`upQuarksDRbar,
               Global`bottomDRbar   -> Global`downQuarksDRbar,
               Global`electronDRbar -> Global`downLeptonsDRbar
           };
           FlexibleSUSY`LowScaleInput  = FlexibleSUSY`LowScaleInput  /. backwardsCompatRules;
           FlexibleSUSY`SUSYScaleInput = FlexibleSUSY`SUSYScaleInput /. backwardsCompatRules;
           FlexibleSUSY`HighScaleInput = FlexibleSUSY`HighScaleInput /. backwardsCompatRules;
           FlexibleSUSY`InitialGuessAtLowScale  = FlexibleSUSY`InitialGuessAtLowScale  /. backwardsCompatRules;
           FlexibleSUSY`InitialGuessAtHighScale = FlexibleSUSY`InitialGuessAtHighScale /. backwardsCompatRules;

           Constraint`CheckConstraint[FlexibleSUSY`LowScaleInput, "LowScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`SUSYScaleInput, "SUSYScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`HighScaleInput, "HighScaleInput"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtLowScale, "InitialGuessAtLowScale", True];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtSUSYScale, "InitialGuessAtSUSYScale"];
           Constraint`CheckConstraint[FlexibleSUSY`InitialGuessAtHighScale, "InitialGuessAtHighScale"];

           (* warn if extra parameters, which do not run, are used at multiple scales *)
           CheckExtraParametersUsage[Parameters`GetExtraParameters[],
                                     {FlexibleSUSY`LowScaleInput, FlexibleSUSY`SUSYScaleInput, FlexibleSUSY`HighScaleInput}];

           (* replace LHInput[p] by pInput in the constraints *)
           EvaluateUserInput[];
           RenameSLHAInputParametersInUserInput[lesHouchesInputParameters];

           If[HaveBVPSolver[FlexibleSUSY`SemiAnalyticSolver],
              SemiAnalytic`SetSemiAnalyticParameters[BetaFunction`GetName[#]& /@ susyBreakingBetaFunctions];

              (* @note currently require all semi-analytic parameters to be set at same scale *)
              If[!SemiAnalytic`CheckSemiAnalyticBoundaryConditions[{FlexibleSUSY`LowScaleInput,
                                                                    FlexibleSUSY`SUSYScaleInput,
                                                                    FlexibleSUSY`HighScaleInput}],
                 Print["Error: the requested boundary conditions are not"];
                 Print["   supported by the semi-analytic solver."];
                 Print["   Please modify the boundary conditions or disable"];
                 Print["   the semi-analytic solver."];
                 Print["   Alternatively, please contact the developers to"];
                 Print["   discuss adding support for these boundary conditions."];
                 Quit[1];
                ];

              semiAnalyticBCs = SemiAnalytic`SelectBoundaryConstraint[{FlexibleSUSY`LowScaleInput,
                                                                       FlexibleSUSY`SUSYScaleInput,
                                                                       FlexibleSUSY`HighScaleInput}];

              semiAnalyticSolns = SemiAnalytic`GetSemiAnalyticSolutions[semiAnalyticBCs];

              (* add boundary values as additional extra parameters if necessary *)
              Parameters`AddExtraParameters[SemiAnalytic`GetExtraBoundaryParametersToSave[semiAnalyticSolns]];
             ];

           (* replace all indices in the user-defined model file variables *)
           ReplaceIndicesInUserInput[allIndexReplacementRules];
           ReplaceIndicesInUserInput[allInputParameterIndexReplacementRules];
           ReplaceIndicesInUserInput[allExtraParameterIndexReplacementRules];

           Utils`PrintHeadline["Creating model parameter classes"];
           Print["Creating class for susy parameters ..."];
           WriteRGEClass[susyBetaFunctions, anomDim,
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "susy_parameters.hpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_susy_parameters.hpp"}]},
                          {FileNameJoin[{$flexiblesusyTemplateDir, "susy_parameters.cpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_susy_parameters.cpp"}]}},
                         "susy_beta_.cpp.in",
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "betas.mk.in"}],
                           FileNameJoin[{FSOutputDir, "susy_betas.mk"}]}}
                        ];

           Print["Creating class for soft parameters ..."];
           WriteRGEClass[susyBreakingBetaFunctions, {},
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "soft_parameters.hpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_soft_parameters.hpp"}]},
                          {FileNameJoin[{$flexiblesusyTemplateDir, "soft_parameters.cpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_soft_parameters.cpp"}]}},
                         "soft_beta_.cpp.in",
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "betas.mk.in"}],
                           FileNameJoin[{FSOutputDir, "soft_betas.mk"}]}},
                         If[Head[SARAH`TraceAbbr] === List, SARAH`TraceAbbr, {}],
                         BetaFunction`CountNumberOfParameters[susyBetaFunctions]];

           (********************* EWSB *********************)
           ewsbEquations = PrepareEWSBEquations[allIndexReplacementRules];

           If[FlexibleSUSY`EWSBInitialGuess =!= {},
              FlexibleSUSY`EWSBInitialGuess = EWSB`GetValidEWSBInitialGuesses[FlexibleSUSY`EWSBInitialGuess];
             ];

           If[FlexibleSUSY`EWSBSubstitutions =!= {},
              sharedEwsbSubstitutions = EWSB`GetValidEWSBSubstitutions[FlexibleSUSY`EWSBSubstitutions];
              sharedEwsbSubstitutions = (Rule @@ #)& /@ (sharedEwsbSubstitutions /. allIndexReplacementRules);
             ];
           solverEwsbSubstitutions = Rule[#, sharedEwsbSubstitutions]& /@ FlexibleSUSY`FSBVPSolvers;

           If[HaveBVPSolver[FlexibleSUSY`SemiAnalyticSolver],
              semiAnalyticEWSBSubstitutions = SemiAnalytic`GetSemiAnalyticEWSBSubstitutions[semiAnalyticSolns];
              solverEwsbSubstitutions = AddEWSBSubstitutionsForSolver[FlexibleSUSY`SemiAnalyticSolver,
                                                                      solverEwsbSubstitutions,
                                                                      semiAnalyticEWSBSubstitutions];
             ];

           FlexibleSUSY`EWSBOutputParameters = Parameters`DecreaseIndexLiterals[FlexibleSUSY`EWSBOutputParameters];
           If[ewsbEquations =!= {},
              treeLevelEwsbEqsOutputFile = FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_EWSB_equations.m"}];
              Print["Writing EWSB equations to ", treeLevelEwsbEqsOutputFile];
              If[sharedEwsbSubstitutions =!= {},
                 Put[Parameters`ReplaceAllRespectingSARAHHeads[ewsbEquations, sharedEwsbSubstitutions], treeLevelEwsbEqsOutputFile],
                 Put[ewsbEquations, treeLevelEwsbEqsOutputFile]
                ];

              treeLevelEwsbSolutionOutputFiles =
                  (Rule[#, FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_"
                                                      <> GetBVPSolverHeaderName[#] <> "_EWSB_solution.m"}]])& /@ FlexibleSUSY`FSBVPSolvers;
              {solverEwsbSolutions, solverFreePhases} = SolveEWSBEquationsForSolvers[FlexibleSUSY`FSBVPSolvers, ewsbEquations,
                                                                                     FlexibleSUSY`EWSBOutputParameters, solverEwsbSubstitutions,
                                                                                     FlexibleSUSY`TreeLevelEWSBSolution,
                                                                                     treeLevelEwsbSolutionOutputFiles];
              ,
              Print["Note: There are no EWSB equations."];
              solverEwsbSolutions = Rule[#, {}]& /@ FlexibleSUSY`FSBVPSolvers;
              solverFreePhases = Rule[#, {}]& /@ FlexibleSUSY`FSBVPSolvers;
             ];
           freePhases = GetAllFreePhases[solverFreePhases];
           If[freePhases =!= {},
              Print["Note: the following phases are free: ", freePhases];
              missingPhases = Select[freePhases, !MemberQ[#[[1]]& /@ inputParameters, #]&];
              If[missingPhases =!= {},
                 Print["Error: the following phases are not defined as input parameters: ", InputForm[missingPhases]];
                 Print["   Please add them to the MINPAR or EXTPAR input parameter lists."];
                 Quit[1];
                ];
             ];

           If[Cases[solverEwsbSolutions, (Rule[solver_, {}]) :> solver] =!= {},
              Print["Warning: an analytic solution to the EWSB eqs. ",
                    " could not be found for the solvers: ",
                    Cases[solverEwsbSolutions, (Rule[solver_ , {}]) :> solver]];
              Print["   An iterative algorithm will be used.  You can try to set"];
              Print["   the solution by hand in the model file like this:"];
              Print[""];
              Print["   TreeLevelEWSBSolution = {"];
              For[i = 1, i <= Length[FlexibleSUSY`EWSBOutputParameters], i++,
              Print["      { ", FlexibleSUSY`EWSBOutputParameters[[i]], ", ... }" <>
                    If[i != Length[FlexibleSUSY`EWSBOutputParameters], ",", ""]];
                   ];
              Print["   };\n"];
             ];
           solverEwsbSolvers = SelectValidEWSBSolvers[solverEwsbSolutions, FlexibleSUSY`FSEWSBSolvers];

           Print["Input parameters: ", InputForm[Parameters`GetInputParameters[]]];

           Print["Creating class for input parameters ..."];
           WriteInputParameterClass[inputParameters,
                                    {{FileNameJoin[{$flexiblesusyTemplateDir, "input_parameters.hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "input_parameters.cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_input_parameters.cpp"}]}
                                    }
                                   ];

           extraSLHAOutputBlocks = Parameters`DecreaseIndexLiterals[
               FlexibleSUSY`ExtraSLHAOutputBlocks,
               Join[Parameters`GetOutputParameters[], Parameters`GetModelParameters[], Parameters`GetExtraParameters[]]
           ];

           (* determine diagonalization precision for each particle *)
           diagonalizationPrecision = ReadPoleMassPrecisions[
               DefaultPoleMassPrecision,
               Flatten[{HighPoleMassPrecision}],
               Flatten[{MediumPoleMassPrecision}],
               Flatten[{LowPoleMassPrecision}],
               FSEigenstates];

           (*prepare Weinberg angle calculation*)
           WeinbergAngle`InitMuonDecay[];
           deltaVBwave = SortCps @ WeinbergAngle`DeltaVBwave[];
           deltaVBvertex = SortCps @ WeinbergAngle`DeltaVBvertex[];
           deltaVBbox = SortCps @ WeinbergAngle`DeltaVBbox[];

           (* prepare decays calculation *)
           If[FlexibleSUSY`FSCalculateDecays,

              FlexibleSUSY`FSDecayParticles = Select[FlexibleSUSY`FSDecayParticles, Decays`IsSupportedDecayParticle];

              If[FlexibleSUSY`FSDecayParticles === {},
                 Utils`FSFancyWarning[
                    "No supported particles to calculate decays for were found.",
                    " Generation of decays code will be skipped."
                 ];
                 FlexibleSUSY`FSCalculateDecays = False;
                ,
                decaysSLHAIncludeFiles = {FlexibleSUSY`FSModelName <> "_decays.hpp"};
                ];
             ]; (* If[FlexibleSUSY`FSCalculateDecays] *)

           vertexRuleFileName =
              GetVertexRuleFileName[$sarahCurrentOutputMainDir, FSEigenstates];
           If[NeedToCalculateVertices[FSEigenstates],
              Put[vertexRules =
                      Vertices`VertexRules[Join[nPointFunctions,
                                                deltaVBwave, deltaVBvertex, deltaVBbox],
                                           Lat$massMatrices],
                                           vertexRuleFileName],
              vertexRules = Get[vertexRuleFileName];
           ];

           (* apply user-defined rules *)
           vertexRules = vertexRules /. FlexibleSUSY`FSVertexRules;

           Utils`PrintHeadline["Creating model"];
           Print["Creating class for model ..."];
           WriteModelClass[massMatrices, ewsbEquations, FlexibleSUSY`EWSBOutputParameters,
                           DeleteDuplicates[Flatten[#[[2]]& /@ solverEwsbSubstitutions]], nPointFunctions,
                           vertexRules, Parameters`GetPhases[],
                           {{FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates.cpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates_interface.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates_interface.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates_decoupling_scheme.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates_decoupling_scheme.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "mass_eigenstates_decoupling_scheme.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_mass_eigenstates_decoupling_scheme.cpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "physical.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_physical.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "physical.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_physical.cpp"}]}
                           },
                           diagonalizationPrecision];

           Utils`PrintHeadline["Creating SLHA model"];
           Print["Creating class for SLHA model ..."];
           WriteModelSLHAClass[massMatrices,
                               {{FileNameJoin[{$flexiblesusyTemplateDir, "model_slha.hpp.in"}],
                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_model_slha.hpp"}]},
                                {FileNameJoin[{$flexiblesusyTemplateDir, "model_slha.cpp.in"}],
                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_model_slha.cpp"}]}
                               }];

           Utils`PrintHeadline["Creating utilities"];
           Print["Creating utilities class ..."];
           WriteUtilitiesClass[massMatrices, Join[susyBetaFunctions, susyBreakingBetaFunctions],
                               inputParameters, Parameters`GetExtraParameters[],
                               FlexibleSUSY`FSLesHouchesList, extraSLHAOutputBlocks,
                               decaysSLHAIncludeFiles,
               {{FileNameJoin[{$flexiblesusyTemplateDir, "info.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_info.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "info.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_info.cpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "utilities.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_utilities.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "utilities.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_utilities.cpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "slha_io.hpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_slha_io.hpp"}]},
                {FileNameJoin[{$flexiblesusyTemplateDir, "slha_io.cpp.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_slha_io.cpp"}]}
               }
                              ];

           Print["Creating FlexibleEFTHiggs.mk ..."];
           WriteFlexibleEFTHiggsMakefileModule[
                              {{FileNameJoin[{$flexiblesusyTemplateDir, "FlexibleEFTHiggs.mk.in"}],
                                FileNameJoin[{FSOutputDir, "FlexibleEFTHiggs.mk"}]}
                              }];

           Print["Creating list of references to be cited ..."];
           WriteReferences[
               {{FileNameJoin[{$flexiblesusyTemplateDir, "references.tex.in"}],
                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_references.tex"}]}}
           ];

           Print["Creating plot scripts ..."];
           WritePlotScripts[{{FileNameJoin[{$flexiblesusyTemplateDir, "plot_spectrum.gnuplot.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_plot_spectrum.gnuplot"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "plot_rgflow.gnuplot.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_plot_rgflow.gnuplot"}]}}
                           ];

           Utils`PrintHeadline["Creating solver framework"];
           Print["Creating generic solver class templates ..."];
           spectrumGeneratorInterfaceInputFile = If[
               FlexibleSUSY`FlexibleEFTHiggs === True,
               "standard_model_spectrum_generator_interface.hpp.in",
               "spectrum_generator_interface.hpp.in"
           ];
           WriteBVPSolverTemplates[{{FileNameJoin[{$flexiblesusyTemplateDir, "convergence_tester.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_convergence_tester.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "ewsb_solver.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_ewsb_solver.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "ewsb_solver_interface.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_ewsb_solver_interface.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "high_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_high_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "initial_guesser.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_initial_guesser.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "low_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_low_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "model.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_model.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "spectrum_generator.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_spectrum_generator.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInterfaceInputFile}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_spectrum_generator_interface.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "susy_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_susy_scale_constraint.hpp"}]}
                                   }];

           Utils`PrintHeadline["Creating Weinberg angle class ..."];
           WriteWeinbergAngleClass[Join[deltaVBwave, deltaVBvertex, deltaVBbox], vertexRules,
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "weinberg_angle.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_weinberg_angle.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "weinberg_angle.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_weinberg_angle.cpp"}]}
                                   }];

           If[HaveBVPSolver[FlexibleSUSY`TwoScaleSolver],
              Utils`PrintHeadline["Creating two-scale solver"];
              Print["Creating class for convergence tester ..."];
              WriteConvergenceTesterClass[FlexibleSUSY`FSConvergenceCheck,
                  {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_convergence_tester.hpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.hpp"}]},
                   {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_convergence_tester.cpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_convergence_tester.cpp"}]}
                  }
                                         ];

              Print["Creating class for high-scale constraint ..."];
              WriteConstraintClass[FlexibleSUSY`HighScale,
                                   FlexibleSUSY`HighScaleInput,
                                   FlexibleSUSY`HighScaleFirstGuess,
                                   {FlexibleSUSY`HighScaleMinimum, FlexibleSUSY`HighScaleMaximum},
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_high_scale_constraint.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_high_scale_constraint.cpp"}]}
                                   }
                                  ];

              Print["Creating class for susy-scale constraint ..."];
              WriteConstraintClass[FlexibleSUSY`SUSYScale,
                                   FlexibleSUSY`SUSYScaleInput,
                                   FlexibleSUSY`SUSYScaleFirstGuess,
                                   {FlexibleSUSY`SUSYScaleMinimum, FlexibleSUSY`SUSYScaleMaximum},
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_susy_scale_constraint.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_susy_scale_constraint.cpp"}]}
                                   }
                                  ];

              Print["Creating class for low-scale constraint ..."];
              WriteConstraintClass[FlexibleSUSY`LowScale,
                                   FlexibleSUSY`LowScaleInput,
                                   FlexibleSUSY`LowScaleFirstGuess,
                                   {FlexibleSUSY`LowScaleMinimum, FlexibleSUSY`LowScaleMaximum},
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_low_scale_constraint.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_low_scale_constraint.cpp"}]}
                                   }
                                  ];

              Print["Creating class for initial guesser ..."];
              If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 initialGuesserInputFile = "two_scale_low_scale_initial_guesser";,
                 initialGuesserInputFile = "two_scale_high_scale_initial_guesser";
                ];
              If[FlexibleSUSY`FlexibleEFTHiggs === True,
                 initialGuesserInputFile = "standard_model_" <> initialGuesserInputFile;
                ];
              WriteInitialGuesserClass[FlexibleSUSY`InitialGuessAtLowScale,
                                       FlexibleSUSY`InitialGuessAtSUSYScale,
                                       FlexibleSUSY`InitialGuessAtHighScale,
                                       {{FileNameJoin[{$flexiblesusyTemplateDir, initialGuesserInputFile <> ".hpp.in"}],
                                         FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.hpp"}]},
                                        {FileNameJoin[{$flexiblesusyTemplateDir, initialGuesserInputFile <> ".cpp.in"}],
                                         FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_initial_guesser.cpp"}]}
                                       }
                                      ];

              Print["Creating class for two-scale EWSB solver ..."];
              WriteEWSBSolverClass[ewsbEquations, FlexibleSUSY`EWSBOutputParameters, FlexibleSUSY`EWSBInitialGuess,
                                   FlexibleSUSY`TwoScaleSolver /. solverEwsbSubstitutions,
                                   FlexibleSUSY`TwoScaleSolver /. solverEwsbSolutions,
                                   FlexibleSUSY`TwoScaleSolver /. solverFreePhases,
                                   FlexibleSUSY`TwoScaleSolver /. solverEwsbSolvers,
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_ewsb_solver.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_ewsb_solver.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_ewsb_solver.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_ewsb_solver.cpp"}]}}];

              Print["Creating class for two-scale model ..."];
              WriteTwoScaleModelClass[{{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_model.hpp.in"}],
                                        FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.hpp"}]},
                                       {FileNameJoin[{$flexiblesusyTemplateDir, "two_scale_model.cpp.in"}],
                                        FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_model.cpp"}]}}];

              If[FlexibleSUSY`FlexibleEFTHiggs === True,
                 Print["Creating two-scale matching class ..."];
                 WriteMatchingClass[FlexibleSUSY`MatchingScaleInput, massMatrices,
                                    {{FileNameJoin[{$flexiblesusyTemplateDir, "standard_model_two_scale_matching.hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_matching.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "standard_model_two_scale_matching.cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_matching.cpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "standard_model_two_scale_matching_interface.hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_two_scale_matching_interface.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "standard_model_two_scale_matching_interface.cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_two_scale_matching_interface.cpp"}]}
                                    }];
                ];

              spectrumGeneratorInputFile = "two_scale_high_scale_spectrum_generator";
              If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 spectrumGeneratorInputFile = "two_scale_low_scale_spectrum_generator";
                ];
              If[FlexibleSUSY`FlexibleEFTHiggs === True,
                 spectrumGeneratorInputFile = "standard_model_" <> spectrumGeneratorInputFile;
                ];
              Print["Creating class for two-scale spectrum generator ..."];
              WriteTwoScaleSpectrumGeneratorClass[{{FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile <> ".hpp.in"}],
                                                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_spectrum_generator.hpp"}]},
                                                   {FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile <> ".cpp.in"}],
                                                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_two_scale_spectrum_generator.cpp"}]}
                                                   }];

              Print["Creating makefile module for two-scale solver ..."];
              WriteBVPSolverMakefile[{{FileNameJoin[{$flexiblesusyTemplateDir, "two_scale.mk.in"}],
                                       FileNameJoin[{FSOutputDir, "two_scale.mk"}]}}];

             ]; (* If[HaveBVPSolver[FlexibleSUSY`TwoScaleSolver] *)

           If[HaveBVPSolver[FlexibleSUSY`SemiAnalyticSolver],
              Utils`PrintHeadline["Creating semi-analytic solver"];

              Parameters`AddExtraParameters[SemiAnalytic`CreateCoefficientParameters[semiAnalyticSolns]];

              semiAnalyticSolnsOutputFile = FileNameJoin[{FSOutputDir,
                                                          FlexibleSUSY`FSModelName <> "_semi_analytic_solutions.m"}];
              Print["Writing semi-analytic solutions to ", semiAnalyticSolnsOutputFile];
              Put[SemiAnalytic`ExpandSemiAnalyticSolutions[semiAnalyticSolns], semiAnalyticSolnsOutputFile];

              Print["Creating classes for convergence testers ..."];
              If[FlexibleSUSY`SemiAnalyticSolverInnerConvergenceCheck === Automatic,
                 FlexibleSUSY`SemiAnalyticSolverInnerConvergenceCheck = Complement[Parameters`GetModelParameters[], SemiAnalytic`GetSemiAnalyticParameters[]];
                ];
              WriteConvergenceTesterClass[FlexibleSUSY`SemiAnalyticSolverInnerConvergenceCheck,
                  {{FileNameJoin[{$flexiblesusyTemplateDir, "susy_convergence_tester.hpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_susy_convergence_tester.hpp"}]},
                   {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_susy_convergence_tester.hpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_susy_convergence_tester.hpp"}]},
                   {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_susy_convergence_tester.cpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_susy_convergence_tester.cpp"}]}
                  }
                                         ];
              WriteConvergenceTesterClass[FlexibleSUSY`FSConvergenceCheck,
                  {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_convergence_tester.hpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_convergence_tester.hpp"}]},
                   {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_convergence_tester.cpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_convergence_tester.cpp"}]}
                  }
                                         ];

              semiAnalyticHighScaleFiles = {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_high_scale_constraint.hpp.in"}],
                                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_high_scale_constraint.hpp"}]},
                                            {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_high_scale_constraint.cpp.in"}],
                                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_high_scale_constraint.cpp"}]}
                                           };
              semiAnalyticSUSYScaleFiles = {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_susy_scale_constraint.hpp.in"}],
                                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_susy_scale_constraint.hpp"}]},
                                            {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_susy_scale_constraint.cpp.in"}],
                                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_susy_scale_constraint.cpp"}]}
                                           };
              semiAnalyticLowScaleFiles = {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_low_scale_constraint.hpp.in"}],
                                            FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_low_scale_constraint.hpp"}]},
                                           {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_low_scale_constraint.cpp.in"}],
                                            FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_low_scale_constraint.cpp"}]}
                                          };
              semiAnalyticConstraintFiles = {{FileNameJoin[{$flexiblesusyTemplateDir, "soft_parameters_constraint.hpp.in"}],
                                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_soft_parameters_constraint.hpp"}]},
                                             {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_soft_parameters_constraint.hpp.in"}],
                                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_soft_parameters_constraint.hpp"}]},
                                             {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_soft_parameters_constraint.cpp.in"}],
                                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_soft_parameters_constraint.cpp"}]}
                                            };

              Which[SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`HighScaleInput],
                    semiAnalyticHighScaleFiles = Join[semiAnalyticHighScaleFiles, semiAnalyticConstraintFiles];,
                    SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`SUSYScaleInput],
                    semiAnalyticSUSYScaleFiles = Join[semiAnalyticSUSYScaleFiles, semiAnalyticConstraintFiles];,
                    SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`LowScaleInput],
                    semiAnalyticLowScaleFiles = Join[semiAnalyticLowScaleFiles, semiAnalyticConstraintFiles];,
                    True,
                    semiAnalyticSUSYScaleFiles = Join[semiAnalyticSUSYScaleFiles, semiAnalyticConstraintFiles];
                   ];

              Print["Creating class for high-scale constraint ..."];
              WriteSemiAnalyticConstraintClass[FlexibleSUSY`HighScale, FlexibleSUSY`HighScaleInput, {},
                                               FlexibleSUSY`HighScaleFirstGuess,
                                               {FlexibleSUSY`HighScaleMinimum, FlexibleSUSY`HighScaleMaximum},
                                               SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`HighScaleInput],
                                               SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`HighScaleInput],
                                               semiAnalyticSolns, semiAnalyticHighScaleFiles];

              Print["Creating class for susy-scale constraint ..."];
              WriteSemiAnalyticConstraintClass[FlexibleSUSY`SUSYScale, FlexibleSUSY`SUSYScaleInput, {},
                                               FlexibleSUSY`SUSYScaleFirstGuess,
                                               {FlexibleSUSY`SUSYScaleMinimum, FlexibleSUSY`SUSYScaleMaximum},
                                               SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`SUSYScaleInput],
                                               SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`SUSYScaleInput],
                                               semiAnalyticSolns, semiAnalyticSUSYScaleFiles];

              Print["Creating class for low-scale constraint ..."];
              WriteSemiAnalyticConstraintClass[FlexibleSUSY`LowScale, FlexibleSUSY`LowScaleInput,
                                               If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY, {},
                                                  FlexibleSUSY`SUSYScaleInput], FlexibleSUSY`LowScaleFirstGuess,
                                               {FlexibleSUSY`LowScaleMinimum, FlexibleSUSY`LowScaleMaximum},
                                               SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`LowScaleInput],
                                               SemiAnalytic`IsSemiAnalyticConstraint[FlexibleSUSY`LowScaleInput],
                                               semiAnalyticSolns, semiAnalyticLowScaleFiles];

              Print["Creating class for initial guesser ..."];
              If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 initialGuesserInputFile = "semi_analytic_low_scale_initial_guesser";,
                 initialGuesserInputFile = "semi_analytic_high_scale_initial_guesser";
                ];
              Which[SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`HighScaleInput],
                 semiAnalyticInputScale = "high_constraint.get_scale()",
                 SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`SUSYScaleInput],
                 semiAnalyticInputScale = "susy_constraint.get_scale()",
                 SemiAnalytic`IsBoundaryConstraint[FlexibleSUSY`LowScaleInput],
                 semiAnalyticInputScale = "low_constraint.get_scale()",
                 True,
                 semiAnalyticInputScale = "high_constraint.get_scale()"
                ];
              WriteSemiAnalyticInitialGuesserClass[FlexibleSUSY`InitialGuessAtLowScale,
                                                   FlexibleSUSY`InitialGuessAtSUSYScale,
                                                   FlexibleSUSY`InitialGuessAtHighScale,
                                                   semiAnalyticInputScale,
                                                   {{FileNameJoin[{$flexiblesusyTemplateDir, initialGuesserInputFile <> ".hpp.in"}],
                                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_initial_guesser.hpp"}]},
                                                    {FileNameJoin[{$flexiblesusyTemplateDir, initialGuesserInputFile <> ".cpp.in"}],
                                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_initial_guesser.cpp"}]}
                                                   }
                                                  ];

              Print["Creating class for semi-analytic solutions ..."];
              WriteSemiAnalyticSolutionsClass[semiAnalyticBCs, semiAnalyticSolns,
                                              {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_solutions.hpp.in"}],
                                                FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_solutions.hpp"}]},
                                               {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_solutions.cpp.in"}],
                                                FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_solutions.cpp"}]}}
                                             ];

              Print["Creating class for semi-analytic EWSB solver ..."];
              WriteSemiAnalyticEWSBSolverClass[ewsbEquations, FlexibleSUSY`EWSBOutputParameters, FlexibleSUSY`EWSBInitialGuess,
                                               FlexibleSUSY`SemiAnalyticSolver /. solverEwsbSubstitutions,
                                               FlexibleSUSY`SemiAnalyticSolver /. solverEwsbSolutions,
                                               FlexibleSUSY`SemiAnalyticSolver /. solverFreePhases,
                                               FlexibleSUSY`SemiAnalyticSolver /. solverEwsbSolvers, semiAnalyticSolns,
                                               {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_ewsb_solver.hpp.in"}],
                                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_ewsb_solver.hpp"}]},
                                                {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_ewsb_solver.cpp.in"}],
                                                 FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_ewsb_solver.cpp"}]}}];

              Print["Creating class for semi-analytic model ..."];
              WriteSemiAnalyticModelClass[semiAnalyticBCs, semiAnalyticSolns,
                                          {{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_model.hpp.in"}],
                                            FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_model.hpp"}]},
                                           {FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic_model.cpp.in"}],
                                            FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_model.cpp"}]}}];

              spectrumGeneratorInputFile = "semi_analytic_high_scale_spectrum_generator";
              If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                 spectrumGeneratorInputFile = "semi_analytic_low_scale_spectrum_generator";
                ];
              Print["Creating class for semi-analytic spectrum generator ..."];
              WriteSemiAnalyticSpectrumGeneratorClass[{{FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile <> ".hpp.in"}],
                                                        FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_spectrum_generator.hpp"}]},
                                                       {FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile <> ".cpp.in"}],
                                                        FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_semi_analytic_spectrum_generator.cpp"}]}
                                                      }];

              Print["Creating makefile module for semi-analytic solver ..."];
              WriteBVPSolverMakefile[{{FileNameJoin[{$flexiblesusyTemplateDir, "semi_analytic.mk.in"}],
                                       FileNameJoin[{FSOutputDir, "semi_analytic.mk"}]}}];

              Parameters`RemoveExtraParameters[SemiAnalytic`CreateCoefficientParameters[semiAnalyticSolns]];
             ]; (* If[HaveBVPSolver[FlexibleSUSY`SemiAnalyticSolver] *)

           If[HaveBVPSolver[FlexibleSUSY`ShootingSolver],
              Print["Creating FlexibleEFTHiggs.mk ..."];
              WriteFlexibleEFTHiggsMakefileModule[
                  {{FileNameJoin[{$flexiblesusyTemplateDir, "FlexibleEFTHiggs.mk.in"}],
                    FileNameJoin[{FSOutputDir, "FlexibleEFTHiggs.mk"}]}
                  }];

              Print["Creating makefile module for shooting solver ..."];
              WriteBVPSolverMakefile[{{FileNameJoin[{$flexiblesusyTemplateDir, "shooting.mk.in"}],
                                       FileNameJoin[{FSOutputDir, "shooting.mk"}]}}];

              Print["Creating class for shooting model ..."];
              WriteTwoScaleModelClass[{{FileNameJoin[{$flexiblesusyTemplateDir, "shooting_model.hpp.in"}],
                                        FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_model.hpp"}]},
                                       {FileNameJoin[{$flexiblesusyTemplateDir, "shooting_model.cpp.in"}],
                                        FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_model.cpp"}]}}];

              Print["Creating class for shooting EWSB solver ..."];
              WriteEWSBSolverClass[ewsbEquations, FlexibleSUSY`EWSBOutputParameters, FlexibleSUSY`EWSBInitialGuess,
                                   FlexibleSUSY`ShootingSolver /. solverEwsbSubstitutions,
                                   FlexibleSUSY`ShootingSolver /. solverEwsbSolutions,
                                   FlexibleSUSY`ShootingSolver /. solverFreePhases,
                                   FlexibleSUSY`ShootingSolver /. solverEwsbSolvers,
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "shooting_ewsb_solver.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_ewsb_solver.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "shooting_ewsb_solver.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_ewsb_solver.cpp"}]}}];

              Print["Creating class for high-scale constraint ..."];
              WriteConstraintClass[FlexibleSUSY`HighScale,
                                   FlexibleSUSY`HighScaleInput,
                                   FlexibleSUSY`HighScaleFirstGuess,
                                   {FlexibleSUSY`HighScaleMinimum, FlexibleSUSY`HighScaleMaximum},
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "shooting_high_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_high_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "shooting_high_scale_constraint.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_high_scale_constraint.cpp"}]}
                                   }
                                  ];

              Print["Creating class for susy-scale constraint ..."];
              WriteConstraintClass[FlexibleSUSY`SUSYScale,
                                   FlexibleSUSY`SUSYScaleInput,
                                   FlexibleSUSY`SUSYScaleFirstGuess,
                                   {FlexibleSUSY`SUSYScaleMinimum, FlexibleSUSY`SUSYScaleMaximum},
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "shooting_susy_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_susy_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "shooting_susy_scale_constraint.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_susy_scale_constraint.cpp"}]}
                                   }
                                  ];

              Print["Creating class for low-scale constraint ..."];
              WriteConstraintClass[FlexibleSUSY`LowScale,
                                   FlexibleSUSY`LowScaleInput,
                                   FlexibleSUSY`LowScaleFirstGuess,
                                   {FlexibleSUSY`LowScaleMinimum, FlexibleSUSY`LowScaleMaximum},
                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "shooting_low_scale_constraint.hpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_low_scale_constraint.hpp"}]},
                                    {FileNameJoin[{$flexiblesusyTemplateDir, "shooting_low_scale_constraint.cpp.in"}],
                                     FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_low_scale_constraint.cpp"}]}
                                   }
                                  ];

              initialGuesserInputFile = "shooting_initial_guesser";
              If[FlexibleSUSY`FlexibleEFTHiggs === True,
                 initialGuesserInputFile = "standard_model_" <> initialGuesserInputFile;
              ];

              WriteInitialGuesserClass[FlexibleSUSY`InitialGuessAtLowScale,
                                       FlexibleSUSY`InitialGuessAtSUSYScale,
                                       FlexibleSUSY`InitialGuessAtHighScale,
                                       {{FileNameJoin[{$flexiblesusyTemplateDir, initialGuesserInputFile <> ".hpp.in"}],
                                         FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_initial_guesser.hpp"}]},
                                        {FileNameJoin[{$flexiblesusyTemplateDir, initialGuesserInputFile <> ".cpp.in"}],
                                         FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_initial_guesser.cpp"}]}
                                       }
              ];

              If[FlexibleSUSY`FlexibleEFTHiggs === True,
                 Print["Creating shooting matching class ..."];
                 WriteMatchingClass[FlexibleSUSY`MatchingScaleInput, massMatrices,
                                    {{FileNameJoin[{$flexiblesusyTemplateDir, "standard_model_shooting_matching.hpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_matching.hpp"}]},
                                     {FileNameJoin[{$flexiblesusyTemplateDir, "standard_model_shooting_matching.cpp.in"}],
                                      FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_standard_model_matching.cpp"}]}
                                    }];
              ];

              spectrumGeneratorInputFile =
                  If[FlexibleSUSY`FlexibleEFTHiggs,
                     "standard_model_",
                     ""
                  ] <>
                  If[FlexibleSUSY`OnlyLowEnergyFlexibleSUSY,
                     "shooting_low_scale_spectrum_generator",
                     "shooting_high_scale_spectrum_generator"];

              Print["Creating class for shooting spectrum generator ..."];
              WriteShootingSpectrumGeneratorClass[
                  {{FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile <> ".hpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_spectrum_generator.hpp"}]},
                   {FileNameJoin[{$flexiblesusyTemplateDir, spectrumGeneratorInputFile <> ".cpp.in"}],
                    FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_shooting_spectrum_generator.cpp"}]}
                  }];

           ]; (* If[HaveBVPSolver[FlexibleSUSY`ShootingSolver] *)

           Utils`PrintHeadline["Creating observables"];

           (* @todo: should all FlexibleDecay tests be moved into FlexibleDecay.mk
              instead of just a variable that controlls wheter we call them or not? *)
           With[{f = FileNameJoin[{"test", "FlexibleDecay.mk"}]},
              If[FileExistsQ[f], DeleteFile[f]];
              WriteString[f, "ENABLE_FLEXIBLEDECAY := no"]
           ];
           With[{f = FileNameJoin[{FSOutputDir, "decays", "FlexibleDecay.mk"}]},
              If[FileExistsQ[f], DeleteFile[f]]
           ];

           If[FSCalculateDecays,
              PrintHeadline["Creating particle decays"];

              decaysFinalStateParticles = Decays`CreateCompleteParticleList[Select[TreeMasses`GetParticles[], !TreeMasses`IsGhost[#]&]];

              If[FlexibleSUSY`FSDecayParticles =!= {},
                 With[{dir=FileNameJoin[{FSOutputDir, "decays"}]}, If[!DirectoryQ[dir], CreateDirectory[dir]]];
                 decaysSources = Join[decaysSources, {FileNameJoin[{"decays", FlexibleSUSY`FSModelName <> "_decay_table.cpp"}],
                                                      FileNameJoin[{"decays", FlexibleSUSY`FSModelName <> "_decays.cpp"}]}
                                                      ];
                 decaysHeaders = Join[decaysHeaders, {FileNameJoin[{"decays", FlexibleSUSY`FSModelName <> "_decay_table.hpp"}],
                                                      FileNameJoin[{"decays", FlexibleSUSY`FSModelName <> "_decays.hpp"}],
                                                      FileNameJoin[{"decays", FlexibleSUSY`FSModelName <> "_decay_amplitudes.hpp"}]}

                                                      ];
                 decaysVertices = WriteDecaysClass[FlexibleSUSY`FSDecayParticles, decaysFinalStateParticles,
                                                   {{FileNameJoin[{$flexiblesusyTemplateDir, "decays", "decay_table.hpp.in"}],
                                                     FileNameJoin[{FSOutputDir, "decays", FlexibleSUSY`FSModelName <> "_decay_table.hpp"}]},
                                                    {FileNameJoin[{$flexiblesusyTemplateDir, "decays", "decay_table.cpp.in"}],
                                                     FileNameJoin[{FSOutputDir, "decays", FlexibleSUSY`FSModelName <> "_decay_table.cpp"}]},
                                                    {FileNameJoin[{$flexiblesusyTemplateDir, "decays", "decays.hpp.in"}],
                                                     FileNameJoin[{FSOutputDir, "decays", FlexibleSUSY`FSModelName <> "_decays.hpp"}]},
                                                    {FileNameJoin[{$flexiblesusyTemplateDir, "decays", "decays.cpp.in"}],
                                                     FileNameJoin[{FSOutputDir, "decays", FlexibleSUSY`FSModelName <> "_decays.cpp"}]},
                                                    {FileNameJoin[{$flexiblesusyTemplateDir, "run_decays.cpp.in"}],
                                                     FileNameJoin[{FSOutputDir,  "run_decays_" <> FlexibleSUSY`FSModelName <> ".cpp"}]}

                                                   }];
                 WriteOut`ReplaceInFiles[{{FileNameJoin[{$flexiblesusyTemplateDir, "decays", "FlexibleDecay.mk.in"}],
                                           FileNameJoin[{FSOutputDir, "decays", "FlexibleDecay.mk"}]}},
                                         {Sequence @@ GeneralReplacementRules[]}
                 ];
                 ,
                 Print["Skipping calculating decays as no particles to calculate decays for were found."];
                ];

              With[{f = FileNameJoin[{"test", "FlexibleDecay.mk"}]},
                 If[FileExistsQ[f], DeleteFile[f]];
                 WriteString[f, "ENABLE_FLEXIBLEDECAY := yes"]
              ];
              ,
              (* create an empty file (release/generate-models.sh requires this file to exist) *)
              CreateFile[FileNameJoin[{FSOutputDir, "decays", "FlexibleDecay.mk"}]];

             ]; (* If[FSCalculateDecays] *)

           Print["Creating lepton EDM class ..."];
           edmFields = DeleteDuplicates @ Cases[Observables`GetRequestedObservables[extraSLHAOutputBlocks],
                                                FlexibleSUSYObservable`EDM[p_[__]|p_] :> p];
           edmVertices =
             WriteEDMClass[edmFields,
                           {{FileNameJoin[{$flexiblesusyTemplateDir, "edm.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_edm.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "edm.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_edm.cpp"}]}}];

           (* b -> s gamma *)
           If[MemberQ[Observables`GetRequestedObservables[extraSLHAOutputBlocks], FlexibleSUSYObservable`bsgamma],
             Print["Creating b->sγ class ..."];
             QToQGammaFields = Join[{BtoSGamma`GetBottomQuark[] -> {BtoSGamma`GetStrangeQuark[], TreeMasses`GetPhoton[]}},
               {BtoSGamma`GetBottomQuark[] -> {BtoSGamma`GetStrangeQuark[], TreeMasses`GetGluon[]}}],
             QToQGammaFields = {}];
           WriteBToSGammaClass[QToQGammaFields,
                           {{FileNameJoin[{$flexiblesusyTemplateDir, "b_to_s_gamma.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_b_to_s_gamma.hpp"}]},
                            {FileNameJoin[{$flexiblesusyTemplateDir, "b_to_s_gamma.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_b_to_s_gamma.cpp"}]}}];

            (* Load and evaluate NPointFunctions write classes for observables *)
            Module[{files, obs, newRules = {}, down, dir, observablesExtraEntity},
               dir = FileNameJoin@{FSOutputDir, "observables"};
               If[!DirectoryQ@dir, CreateDirectory@dir];

               files = Utils`DynamicInclude@$observablesWildcard@"FlexibleSUSY.m";
               obs = StringSplit[files, $PathnameSeparator][[All, -2]];

               files = {
                     FileNameJoin@{$flexiblesusyTemplateDir, "observables", #<>".in"},
                     FileNameJoin@{dir, FSModelName<>"_"<>#}
                  } &/@ {#<>".hpp", #<>".cpp"} &/@ (Observables`GetObservableFileName/@obs);

               Do[
                  ObservablesExtraOutput@obs[[i]] = WriteClass[Symbol["FlexibleSUSYObservable`"<>obs[[i]]], extraSLHAOutputBlocks, files[[i]]];

                  observablesExtraEntity = FilterRules[ObservablesExtraOutput@obs[[i]], "C++ vertices"];
                  Switch[observablesExtraEntity,
                     {_ -> _}, observablesExtraVertices = Join[observablesExtraVertices, observablesExtraEntity[[1, 2]]],
                     _, Null
                  ];

                  observablesExtraEntity = FilterRules[ObservablesExtraOutput@obs[[i]], "C++ replacements"];
                  Switch[observablesExtraEntity,
                     {_ -> _}, newRules = Join[newRules, observablesExtraEntity[[1, 2]]],
                     _, Null
                  ];,
                  {i, Length@obs}
               ];
               ObservablesExtraOutput@_ = {};

               (* Inserting new rules before default ones *)
               down = DownValues@GeneralReplacementRules;
               down = Insert[down, newRules, {1,2,-1}];
               DownValues@GeneralReplacementRules = down;

               observablesExtraVertices = DeleteDuplicates@observablesExtraVertices;
           ];

           Print["Creating lepton AMM class ..."];
           ammFields = DeleteDuplicates @ Cases[Observables`GetRequestedObservables[extraSLHAOutputBlocks],
                                                FlexibleSUSYObservable`AMM[p_[__]|p_] :> p];
           aMMVertices = WriteAMMClass[
              ammFields,
              {{FileNameJoin[{$flexiblesusyTemplateDir, "amm.hpp.in"}],
                               FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_amm.hpp"}]},
               {FileNameJoin[{$flexiblesusyTemplateDir, "lepton_amm_wrapper.hpp.in"}],
                               FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_lepton_amm_wrapper.hpp"}]},
               {FileNameJoin[{$flexiblesusyTemplateDir, "amm.cpp.in"}],
                               FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_amm.cpp"}]},
               {FileNameJoin[{$flexiblesusyTemplateDir, "lepton_amm_wrapper.cpp.in"}],
                               FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_lepton_amm_wrapper.cpp"}]}}];

           Utils`PrintHeadline["Creating other observables"];
           Print["Creating class for observables ..."];
           WriteObservables[extraSLHAOutputBlocks,
                            {{FileNameJoin[{$flexiblesusyTemplateDir, "observables.hpp.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_observables.hpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "observables.cpp.in"}],
                              FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_observables.cpp"}]}}];

           Print["Creating FFMasslessV form factor class for other observables ..."];
           FFMasslessVVertices =
               WriteFFVFormFactorsClass[
                  (* collect external states from observables needing massless triangles *)
                  DeleteDuplicates @ Join[

                     (* lepton g-2 *)
                     (# -> {#, TreeMasses`GetPhoton[]})& /@ ammFields,

                     (* lepton edm *)
                     (# -> {#, TreeMasses`GetPhoton[]})& /@ edmFields,

                     (* Br(L -> L Gamma) *)
                     Module[{fields = Cases[ObservablesExtraOutput@"BrLToLGamma", ("FFV fields" -> x_) :> x]},
                        Switch[fields, {{}}, {}, _, First/@fields]
                     ],

                     Module[{fields = Flatten@Cases[ObservablesExtraOutput@"BrDLToDL", ("FFV fields" -> x_) :> x]},
                        Switch[fields,
                           {__Symbol}, (# -> {#, TreeMasses`GetPhoton[]})& /@ fields,
                           _, {}
                        ]
                     ],

                     (* b -> s gamma *)
                     QToQGammaFields,

                     Module[{npfFields = FilterRules[ObservablesExtraOutput@"LToLConversion", "FFV fields"]},
                        Switch[npfFields,
                           {_ -> _}, (#[[1]] -> {#[[2]], TreeMasses`GetPhoton[]}) &/@ npfFields[[1, 2]],
                           _, {}
                        ]
                     ],

                     Module[{npfFields = FilterRules[ObservablesExtraOutput@"BrLTo3L", "FFV fields"]},
                        Switch[npfFields,
                           {_ -> _}, (#[[1]] -> {#[[2]], TreeMasses`GetPhoton[]}) &/@ npfFields[[1, 2]],
                           _, {}
                        ]
                     ]
                  ],

                  {{FileNameJoin[{$flexiblesusyTemplateDir, "FFV_form_factors.hpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_FFV_form_factors.hpp"}]},
                     {FileNameJoin[{$flexiblesusyTemplateDir, "FFV_form_factors.cpp.in"}],
                             FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_FFV_form_factors.cpp"}]}}
               ];

           Print["Creating unitarity class..."];
           WriteUnitarityClass[{{FileNameJoin[{$flexiblesusyTemplateDir, "unitarity.hpp.in"}],
                               FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_unitarity.hpp"}]},
                           {FileNameJoin[{$flexiblesusyTemplateDir, "unitarity.cpp.in"}],
                               FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_unitarity.cpp"}]}}
           ];

           Print["Creating C++ QFT class..."];
           cxxQFTTemplateDir = FileNameJoin[{$flexiblesusyTemplateDir, "cxx_qft"}];
           cxxQFTOutputDir = FileNameJoin[{FSOutputDir, "cxx_qft"}];
           cxxQFTFiles = {{FileNameJoin[{cxxQFTTemplateDir, "qft.hpp.in"}],
                           FileNameJoin[{cxxQFTOutputDir, FlexibleSUSY`FSModelName <> "_qft.hpp"}]},
                          {FileNameJoin[{cxxQFTTemplateDir, "fields.hpp.in"}],
                           FileNameJoin[{cxxQFTOutputDir, FlexibleSUSY`FSModelName <> "_fields.hpp"}]},
                          {FileNameJoin[{cxxQFTTemplateDir, "vertices.hpp.in"}],
                           FileNameJoin[{cxxQFTOutputDir, FlexibleSUSY`FSModelName <> "_vertices.hpp"}]},
                          {FileNameJoin[{cxxQFTTemplateDir, "context_base.hpp.in"}],
                           FileNameJoin[{cxxQFTOutputDir, FlexibleSUSY`FSModelName <> "_context_base.hpp"}]},
                          {FileNameJoin[{cxxQFTTemplateDir, "npointfunctions_wilsoncoeffs.hpp.in"}],
                           FileNameJoin[{cxxQFTOutputDir, FlexibleSUSY`FSModelName <> "_npointfunctions_wilsoncoeffs.hpp"}]}
                          };
           cxxQFTVerticesTemplate = FileNameJoin[{cxxQFTTemplateDir, "vertices_.cpp.in"}];
           cxxQFTVerticesMakefileTemplates = {{FileNameJoin[{cxxQFTTemplateDir, "vertices.mk.in"}],
                           FileNameJoin[{cxxQFTOutputDir, "vertices.mk"}]}};

           If[DirectoryQ[cxxQFTOutputDir] === False,
              CreateDirectory[cxxQFTOutputDir]];

           WriteCXXDiagramClass[
              Join[aMMVertices, FFMasslessVVertices, decaysVertices, observablesExtraVertices],
              cxxQFTFiles,
              cxxQFTVerticesTemplate, cxxQFTOutputDir,
              cxxQFTVerticesMakefileTemplates
           ];

           WriteSMParticlesAliases[{{FileNameJoin[{cxxQFTTemplateDir, "particle_aliases.hpp.in"}],
                                     FileNameJoin[{cxxQFTOutputDir, FlexibleSUSY`FSModelName <> "_particle_aliases.hpp"}]}}
                                  ];

           Utils`PrintHeadline["Creating Mathematica interface"];
           Print["Creating LibraryLink ", FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> ".mx"}], " ..."];
           WriteMathLink[inputParameters, extraSLHAOutputBlocks,
                         {{FileNameJoin[{$flexiblesusyTemplateDir, "librarylink.cpp.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_librarylink.cpp"}]},
                          {FileNameJoin[{$flexiblesusyTemplateDir, "librarylink.m.in"}],
                           FileNameJoin[{FSOutputDir, FlexibleSUSY`FSModelName <> "_librarylink.m"}]},
                          {FileNameJoin[{$flexiblesusyTemplateDir, "run.m.in"}],
                           FileNameJoin[{FSOutputDir, "run_" <> FlexibleSUSY`FSModelName <> ".m"}]}
                         }];

           Utils`PrintHeadline["Creating user examples"];
           Print["Creating user example spectrum generator program ..."];
           WriteUserExample[inputParameters,
                            {{FileNameJoin[{$flexiblesusyTemplateDir, "run.cpp.in"}],
                              FileNameJoin[{FSOutputDir, "run_" <> FlexibleSUSY`FSModelName <> ".cpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "run_cmd_line.cpp.in"}],
                              FileNameJoin[{FSOutputDir, "run_cmd_line_" <> FlexibleSUSY`FSModelName <> ".cpp"}]},
                             {FileNameJoin[{$flexiblesusyTemplateDir, "scan.cpp.in"}],
                              FileNameJoin[{FSOutputDir, "scan_" <> FlexibleSUSY`FSModelName <> ".cpp"}]}
                            }
                           ];

           Print["Creating example SLHA input file ..."];
           WriteSLHAInputFile[inputParameters,
                              {{FileNameJoin[{$flexiblesusyTemplateDir, "LesHouches.in"}],
                                FileNameJoin[{FSOutputDir, "LesHouches.in." <> FlexibleSUSY`FSModelName <> "_generated"}]}}
                             ];

           Utils`PrintHeadline["FlexibleSUSY has finished"];
          ];

End[];

EndPackage[];
