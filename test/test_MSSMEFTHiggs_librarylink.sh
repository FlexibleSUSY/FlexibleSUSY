#!/bin/sh

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils
MODELSDIR=${BASEDIR}/../models
MATH=${MATH_CMD:-math}

. "$BASEDIR/test.sh"

slha_input="
Block MODSEL                 # Select model
    6   0                    # flavour violation
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166380000e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.384               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR                 # Input parameters
    4   1                    # SignMu
Block EXTPAR
    0   2000                 # Ms
    1   2000                 # M1(MSUSY)
    2   2000                 # M2(MSUSY)
    3   2000                 # M3(MSUSY)
    4   2000                 # Mu(MSUSY)
    5   2000                 # mA(MSUSY)
   25   5                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     4.00000000E+06   # mq2(1,1)
  2  2     4.00000000E+06   # mq2(2,2)
  3  3     4.00000000E+06   # mq2(3,3)
Block MSE2IN
  1  1     4.00000000E+06   # me2(1,1)
  2  2     4.00000000E+06   # me2(2,2)
  3  3     4.00000000E+06   # me2(3,3)
Block MSL2IN
  1  1     4.00000000E+06   # ml2(1,1)
  2  2     4.00000000E+06   # ml2(2,2)
  3  3     4.00000000E+06   # ml2(3,3)
Block MSU2IN
  1  1     4.00000000E+06   # mu2(1,1)
  2  2     4.00000000E+06   # mu2(2,2)
  3  3     4.00000000E+06   # mu2(3,3)
Block MSD2IN
  1  1     4.00000000E+06   # md2(1,1)
  2  2     4.00000000E+06   # md2(2,2)
  3  3     4.00000000E+06   # md2(3,3)
Block AUIN
  1  1     400   # Au(1,1)
  2  2     400   # Au(2,2)
  3  3     400   # Au(3,3)
Block ADIN
  1  1     1e4   # Ad(1,1)
  2  2     1e4   # Ad(2,2)
  3  3     1e4   # Ad(3,3)
Block AEIN
  1  1     1e4   # Ad(1,1)
  2  2     1e4   # Ad(2,2)
  3  3     1e4   # Ad(3,3)
"

ll_input="
Get[\"${MODELSDIR}/MSSMEFTHiggs/MSSMEFTHiggs_librarylink.m\"];

Off[FSMSSMEFTHiggs::info];
Off[FSMSSMEFTHiggsCalculateSpectrum::warning];
Off[FSMSSMEFTHiggsCalculateSpectrum::error];

handle = FSMSSMEFTHiggsOpenHandle[
    fsSettings -> {
        precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
        maxIterations -> 0,                (* FlexibleSUSY[1] *)
        calculateStandardModelMasses -> 1, (* FlexibleSUSY[3] *)
        poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
        ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
        betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
        thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
        higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
        higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
        higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
        higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
        forceOutput -> 0,                  (* FlexibleSUSY[12] *)
        topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
        betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
        forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
        poleMassScale -> 0,                (* FlexibleSUSY[17] *)
        eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
        eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
        eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
        eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
        eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
        calculateBSMMasses -> 0,           (* FlexibleSUSY[23] *)
        parameterOutputScale -> 0          (* MODSEL[12] *)
    },
    fsSMParameters -> {
        alphaEmMZ -> 1/127.944, (* SMINPUTS[1] *)
        GF -> 1.16637*^-5,      (* SMINPUTS[2] *)
        alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
        MZ -> 91.1876,          (* SMINPUTS[4] *)
        mbmb -> 4.18,           (* SMINPUTS[5] *)
        Mt -> 173.34,           (* SMINPUTS[6] *)
        Mtau -> 1.777,          (* SMINPUTS[7] *)
        Mv3 -> 0,               (* SMINPUTS[8] *)
        MW -> 80.385,           (* SMINPUTS[9] *)
        Me -> 0.000510998902,   (* SMINPUTS[11] *)
        Mv1 -> 0,               (* SMINPUTS[12] *)
        Mm -> 0.1056583715,     (* SMINPUTS[13] *)
        Mv2 -> 0,               (* SMINPUTS[14] *)
        md2GeV -> 0.00475,      (* SMINPUTS[21] *)
        mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
        ms2GeV -> 0.104,        (* SMINPUTS[23] *)
        mcmc -> 1.27,           (* SMINPUTS[24] *)
        CKMTheta12 -> 0,
        CKMTheta13 -> 0,
        CKMTheta23 -> 0,
        CKMDelta -> 0,
        PMNSTheta12 -> 0,
        PMNSTheta13 -> 0,
        PMNSTheta23 -> 0,
        PMNSDelta -> 0,
        PMNSAlpha1 -> 0,
        PMNSAlpha2 -> 0,
        alphaEm0 -> 1/137.035999074,
        Mh -> 125.09
    },
    fsModelParameters -> {
        SignMu -> 0,
        MSUSY -> 2000,
        M1Input -> 2000,
        M2Input -> 2000,
        M3Input -> 2000,
        MuInput -> 2000,
        mAInput -> 2000,
        TanBeta -> 5,
        mq2Input -> 2000^2 IdentityMatrix[3],
        mu2Input -> 2000^2 IdentityMatrix[3],
        md2Input -> 2000^2 IdentityMatrix[3],
        ml2Input -> 2000^2 IdentityMatrix[3],
        me2Input -> 2000^2 IdentityMatrix[3],
        AuInput -> 4 10^2 IdentityMatrix[3],
        AdInput -> 1 10^4 IdentityMatrix[3],
        AeInput -> 1 10^4 IdentityMatrix[3]
    }
];

spectrum    = MSSMEFTHiggs /. FSMSSMEFTHiggsCalculateSpectrum[handle];
observables = FSMSSMEFTHiggsCalculateObservables[handle];
FSMSSMEFTHiggsCloseHandle[handle];

Print[(Pole[M[hh]] /. spectrum)[[1]]];
"

Mh_slha=$(echo "$slha_input" | ${MODELSDIR}/MSSMEFTHiggs/run_MSSMEFTHiggs.x --slha-input-file=- 2>/dev/null | \
                 awk -f ${UTILSDIR}/print_slha_block.awk -v block="MASS" | \
                 awk -f ${UTILSDIR}/print_slha_block_entry.awk -v entries="25" | tail -n 1)

ll_script="${BASEDIR}/run_MSSMEFTHiggs.m"
echo "$ll_input" > "$ll_script"

Mh_ll=$($MATH -noprompt -run "commandLineArg={\"$1\"}; Get[\"$ll_script\"]; Exit[];" | tail -n 1)

echo "Mh SLHA        = $Mh_slha"
echo "Mh LibraryLink = $Mh_ll"

error=0

CHECK_EQUAL_FRACTION "$Mh_slha" "$Mh_ll" "0.00001" || error=$(expr $error + 1)

if [ "x$error" != "x0" ] ; then
    echo "Test FAILED: There were $error errors."
else
    echo "All tests passed."
fi

rm -f "$ll_script"

exit $error
