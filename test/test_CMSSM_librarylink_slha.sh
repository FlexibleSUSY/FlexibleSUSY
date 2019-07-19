#!/bin/sh

BASEDIR="$(dirname $0)"
MODELDIR="${BASEDIR}/../models"
UTILSDIR="${BASEDIR}/../utils"
inputFile1="${BASEDIR}/test_CMSSM_librarylink_slha.in.m";
inputFile2="${BASEDIR}/test_CMSSM_librarylink_slha.in.spc";
outputFile1="${BASEDIR}/test_CMSSM_librarylink_slha.out1.spc";
outputFile2="${BASEDIR}/test_CMSSM_librarylink_slha.out2.spc";
MATH=${MATH_CMD:-math}

rm -f "$inputFile1" "$outputFile1" "$outputFile2"

cat <<EOF > "$inputFile1"
Get["${MODELDIR}/CMSSM/CMSSM_librarylink.m"];

settings = {
    precisionGoal -> 0.0001,
    maxIterations -> 0,
    calculateStandardModelMasses -> 0,
    poleMassLoopOrder -> 2,
    ewsbLoopOrder -> 2,
    betaFunctionLoopOrder -> 3,
    thresholdCorrectionsLoopOrder -> 2,
    higgs2loopCorrectionAtAs -> 1,
    higgs2loopCorrectionAbAs -> 1,
    higgs2loopCorrectionAtAt -> 1,
    higgs2loopCorrectionAtauAtau -> 1,
    forceOutput -> 0,
    topPoleQCDCorrections -> 1,
    betaZeroThreshold -> 1.*10^-11,
    forcePositiveMasses -> 0,
    poleMassScale -> 0.,
    thresholdCorrections -> 123111321,
    parameterOutputScale -> 1000
};

smInputs = {
    alphaEmMZ -> 1./127.934,
    GF -> 0.0000116637,
    alphaSMZ -> 0.1176,
    MZ -> 91.1876,
    mbmb -> 4.2,
    Mt -> 173.3,
    Mtau -> 1.777,
    Mv3 -> 0,
    MW -> 80.404,
    Me -> 0.000510998902,
    Mv1 -> 0.,
    Mm -> 0.105658357,
    Mv2 -> 0.,
    md2GeV -> 0.00475,
    mu2GeV -> 0.0024,
    ms2GeV -> 0.104,
    mcmc -> 1.27,
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
};

handle = FSCMSSMOpenHandle[
    fsSettings -> settings,
    fsSMParameters -> smInputs,
    fsModelParameters -> {
        m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
];

FSCMSSMCalculateSpectrum[handle];
FSCMSSMCalculateObservables[handle];
Export["${outputFile1}", FSCMSSMToSLHA[handle], "String"];
FSCMSSMCloseHandle[handle];
EOF

"$MATH" -run "<< \"$inputFile1\"; Quit[]"

"${MODELDIR}/CMSSM/run_CMSSM.x" \
    --slha-input-file="$inputFile2" \
    --slha-output-file="$outputFile2" > /dev/null

# remove comments
for f in "$outputFile1" "$outputFile2" ; do
    mv "$f" "$f~"
    sed -e 's/ *#.*$//' "$f~" | \
        awk -f "${UTILSDIR}"/remove_slha_block -v block=FlexibleSUSY -v entry=15 \
        > "$f"
done

numdiff --absolute-tolerance=1.0e-12 \
        --relative-tolerance=5.0e-7 \
        "$outputFile1" "$outputFile2"

errors="$?"

if [ $errors = 0 ] ; then
    echo "Test result: OK"
else
    echo "Test result: FAIL"
fi

rm -f "$inputFile1" "$outputFile1" "$outputFile2"

exit $errors
