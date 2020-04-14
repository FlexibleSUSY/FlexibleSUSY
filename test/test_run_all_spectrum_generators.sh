#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
HOMEDIR="${BASEDIR}/.."
FSCONFIG="${HOMEDIR}/flexiblesusy-config"

DEFAULT_CMSSM_INPUT="${HOMEDIR}/model_files/CMSSM/LesHouches.in.CMSSM"
DEFAULT_MSSM_INPUT="${HOMEDIR}/model_files/MSSM/LesHouches.in.MSSM"
DEFAULT_SM_INPUT="${HOMEDIR}/model_files/SM/LesHouches.in.SM"

SGS="\
CMSSM,${DEFAULT_CMSSM_INPUT},0
CMSSMEFTHiggs,_DEFAULT_,0
CMSSMEFTHiggs,${BASEDIR}/test_CMSSMEFTHiggs_no_ewsb.spc,1
CMSSMConvergenceTester,${DEFAULT_CMSSM_INPUT},0
CMSSMFPIAbsolute,${DEFAULT_CMSSM_INPUT},0
CMSSMFPIRelative,${DEFAULT_CMSSM_INPUT},0
CMSSMFPITadpole,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLBroyden,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLHybrid,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLHybridS,${DEFAULT_CMSSM_INPUT},0
CMSSMGSLNewton,${DEFAULT_CMSSM_INPUT},0
CMSSMMassWInput,${DEFAULT_CMSSM_INPUT},0
CMSSMMatchedAtMTDRbar,${DEFAULT_CMSSM_INPUT},0
CMSSMMatchedAtMTPole,${DEFAULT_CMSSM_INPUT},0
CMSSMNoFV,_DEFAULT_,0
CMSSMCKM,_DEFAULT_,0
CMSSMCPV,_DEFAULT_,0
CMSSMCPV,${BASEDIR}/test_CMSSMCPV_wrong_higgs_state.in.spc,0
CMSSMYt2L,${DEFAULT_CMSSM_INPUT},0
MSSMMuBMu,_DEFAULT_,0
NUHMSSMalt,_DEFAULT_,0
NUHMSSMaltEFTHiggs,_DEFAULT_,0
cCMSSM,_DEFAULT_,0
E6SSM,_DEFAULT_,0
E6SSM,${BASEDIR}/test_E6SSM_nan.in.spc,1
CE6SSM,_DEFAULT_,0
E6SSMEFTHiggs,_DEFAULT_,0
InertMSSM,${DEFAULT_CMSSM_INPUT},0
LHInputMSSM,_DEFAULT_,0
LRLR,_DEFAULT_,0
lowMSSM,_DEFAULT_,0
lowNMSSM,_DEFAULT_,0
lowNMSSM,${BASEDIR}/test_lowNMSSM_goldstone_tachyon.in.spc,0
lowNMSSM,${BASEDIR}/test_lowNMSSM_pseudoscalar.in.spc,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.TP1,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.TP2,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.TP3,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.TP4,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.TP5,0
lowNMSSM,${HOMEDIR}/model_files/lowNMSSM/LesHouches.in.TP6,0
lowNMSSMTanBetaAtMZ,_DEFAULT_,0
lowNMSSMTanBetaAtMZ,${HOMEDIR}/model_files/lowNMSSMTanBetaAtMZ/LesHouches.in.TP1,0
lowNMSSMTanBetaAtMZ,${HOMEDIR}/model_files/lowNMSSMTanBetaAtMZ/LesHouches.in.TP2,0
lowNMSSMTanBetaAtMZ,${HOMEDIR}/model_files/lowNMSSMTanBetaAtMZ/LesHouches.in.TP3,0
lowNMSSMTanBetaAtMZ,${HOMEDIR}/model_files/lowNMSSMTanBetaAtMZ/LesHouches.in.TP4,0
lowNMSSMTanBetaAtMZ,${HOMEDIR}/model_files/lowNMSSMTanBetaAtMZ/LesHouches.in.TP5,0
lowNMSSMTanBetaAtMZ,${HOMEDIR}/model_files/lowNMSSMTanBetaAtMZ/LesHouches.in.TP6,0
MDM,_DEFAULT_,0
minMSSM,${DEFAULT_CMSSM_INPUT},1
MRSSM,_DEFAULT_,0
MRSSM2,_DEFAULT_,0
MRSSMEFTHiggs,_DEFAULT_,0
MSSM,_DEFAULT_,0
MSSMCPV,_DEFAULT_,0
MSSMatMGUT,_DEFAULT_,0
MSSMNoFV,_DEFAULT_,0
MSSMNoFVHimalaya,_DEFAULT_,0
MSSMNoFVatMGUT,_DEFAULT_,0
MSSMNoFVatMGUTHimalaya,_DEFAULT_,0
complexMSSM,_DEFAULT_,0
complexMSSM,${DEFAULT_MSSM_INPUT},0
munuSSM,_DEFAULT_,0
NMSSM,_DEFAULT_,0
CNMSSM,_DEFAULT_,0
NMSSMCPV,_DEFAULT_,0
NMSSMEFTHiggs,_DEFAULT_,0
NMSSMEFTHiggs,${HOMEDIR}/model_files/NMSSMEFTHiggs/LesHouches.in.NMSSMEFTHiggs_1507.05093_TP3,0
NoInputParameters,${DEFAULT_SM_INPUT},0
NoYukawaMSSM,${DEFAULT_CMSSM_INPUT},1
NSM,_DEFAULT_,0
NUHMSSM,_DEFAULT_,0
NUHMSSMNoFV,_DEFAULT_,0
NUHMSSMNoFVHimalaya,_DEFAULT_,0
NUTNMSSM,_DEFAULT_,0
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_GTP1,0
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_GTP2,0
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP1,1
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP2,1
NUTNMSSM,${HOMEDIR}/model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_BP3,0
NUTSMSSM,_DEFAULT_,0
rootMSSM,${DEFAULT_CMSSM_INPUT},1
SM,_DEFAULT_,0
SM,${BASEDIR}/test_SM_0L_RGEs.in.spc,0
cSM,${DEFAULT_SM_INPUT},0
SMEWSBAtMZ,${DEFAULT_SM_INPUT},0
SMHighPrecision,${DEFAULT_SM_INPUT},0
SMEFTHiggs,_DEFAULT_,0
SMSSM,_DEFAULT_,0
SplitMSSM,_DEFAULT_,0
SMRules,${DEFAULT_SM_INPUT},0
SSM,_DEFAULT_,0
SMThrow,${DEFAULT_SM_INPUT},0
SMThrow,${HOMEDIR}/model_files/SMThrow/LesHouches.in.SMThrow_large_lambda,1
HSSUSY,_DEFAULT_,0
HSSUSY,${BASEDIR}/test_HSSUSY_SUSYHD_msq_msu_m3_msusy_degenerate.in.spc,0
HSSUSY,${BASEDIR}/test_HSSUSY_SUSYHD_msq_msu_m3_degenerate.in.spc,0
HSSUSY,${BASEDIR}/test_HSSUSY_SUSYHD_msq_m3_degenerate.in.spc,0
HSSUSY,${BASEDIR}/test_HSSUSY_SUSYHD_msu_m3_degenerate.in.spc,0
HSSUSY,${BASEDIR}/test_HSSUSY_SUSYHD_nondegenerate.in.spc,0
THDMII,_DEFAULT_,0
THDMIIMSSMBC,_DEFAULT_,0
HTHDMIIMSSMBC,_DEFAULT_,0
HGTHDMIIMSSMBC,_DEFAULT_,0
THDMIIMSSMBCFull,_DEFAULT_,0
HGTHDMIIMSSMBCFull,_DEFAULT_,0
TMSSM,_DEFAULT_,0
UMSSM,_DEFAULT_,0
YukawaCMSSM,${DEFAULT_CMSSM_INPUT},0
DiracGauginos,_DEFAULT_,0
U1xMSSM3G,_DEFAULT_,0
BLSM,_DEFAULT_,0
BLSMlightZp,_DEFAULT_,0
VCMSSM,_DEFAULT_,0
"

TMP_FILE="${BASEDIR}/test_spectrum_generator.sh.tmp"

rm -f "$TMP_FILE"

errors=0

for setup in ${SGS}
do
    model="`echo ${setup} | tr ',' ' ' | awk '{ print $1 }'`"
    input="`echo ${setup} | tr ',' ' ' | awk '{ print $2 }'`"
    expected_result="`echo ${setup} | tr ',' ' ' | awk '{ print $3 }'`"

    echo "== $model ===================================="

    if [ $("$FSCONFIG" --with-${model}) = no ] ; then
        echo "> skipping, because the model is not configured"
        continue
    fi

    sg="${HOMEDIR}/models/${model}/run_${model}.x"

    if test ! -x "${sg}"; then
        echo "> Error: spectrum generator not built: $sg"
        errors=1
        continue
    fi

    case "$input" in
        _DEFAULT_) input="${HOMEDIR}/model_files/${model}/LesHouches.in.${model}" ;;
    esac

    if test ! -e "${input}"; then
        echo "> Error: input file not found: $input"
        errors=1
        continue
    fi

    cmd="${sg} --slha-input-file=${input} --slha-output-file=${TMP_FILE} > /dev/null 2>&1"

    echo "> running spectrum generator for $model"
    echo "> input file: ${input}"
    echo "> command: ${cmd}"
    eval "${cmd}"

    exit_code="$?"
    echo "> exit code: ${exit_code}"
    echo "> expected result: ${expected_result}"

    if test ${expected_result} -eq ${exit_code} ; then
        echo "> spectrum generator: ok"
    else
        echo "> spectrum generator: FAIL"
        errors=1
    fi
done

rm -f "$TMP_FILE"

echo ""

if test ${errors} -eq 0 ; then
    echo "Test result: OK"
else
    echo "Test result: FAIL"
fi

exit ${errors}
