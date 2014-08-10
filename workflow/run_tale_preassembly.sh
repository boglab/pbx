#!/bin/bash

# Requires
# PBX_SCRIPTS_PATH
# PBX_SMRTANALYSIS_PATH
# PBX_RAW_READS_PATH
# PBX_RESULTS_PATH_WHITELISTING
# PBX_PROTOCOL_TEMPLATES_PATH

PATH="${PBX_SCRIPTS_PATH}/smrtanalysis:${PATH}"
source ${PBX_SMRTANALYSIS_PATH}/current/etc/setup.sh

RESULTS_PATH=${PBX_RESULTS_PATH_PREASSEMBLY}
mkdir -p ${RESULTS_PATH}

find ${PBX_RAW_READS_PATH} -name "*.bax.h5" | sort > ${RESULTS_PATH}/pbx_preassembly_input.fofn

WHITELIST_REPLACEMENT=$(echo "${PBX_RESULTS_PATH_WHITELISTING}" | sed -e 's/[\/&]/\\&/g')

sed "s/__PBX_WHITELISTING_RESULTS__/${WHITELIST_REPLACEMENT}/g" ${PBX_PROTOCOL_TEMPLATES_PATH}/RS_PreAssembler_TALs.1.xml > ${RESULTS_PATH}/pbx_preassembly_protocol.xml

fofnToSmrtpipeInput.py ${RESULTS_PATH}/pbx_preassembly_input.fofn --jobname="TALE_PreAssembly_Job" > ${RESULTS_PATH}/pbx_preassembly_input.xml
smrtpipe.py --params=${RESULTS_PATH}/pbx_preassembly_protocol.xml --output=${RESULTS_PATH} xml:${RESULTS_PATH}/pbx_preassembly_input.xml
python2 ${PBX_SCRIPTS_PATH}/smrtanalysis/trimFastqByQVWindow.py --qvCut=50 --minSeqLen=4000 ${RESULTS_PATH}/data/corrected.fastq > ${RESULTS_PATH}/data/corrected_50_4000.fastq
