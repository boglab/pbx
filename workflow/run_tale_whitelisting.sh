#!/bin/bash

# Requires
# PBX_SMRTANALYSIS_PATH
# PBX_RAW_READS_PATH
# PBX_TALE_SEQS_WHITELISTING
# PBX_SCRIPTS_PATH
# PBX_WHITELISTING_RESULTS_PATH

source ${PBX_SMRTANALYSIS_PATH}/current/etc/setup.sh

RESULTS_PATH=${PBX_WHITELISTING_RESULTS_PATH}
mkdir -p ${RESULTS_PATH}

cd ${PBX_RAW_READS_PATH}
for f in *.bax.h5; do fprefix=`echo "$f" | rev | cut -c 8- | rev`; echo $fprefix; blasr ${f} ${PBX_TALE_SEQS_WHITELISTING} -bestn 1 -header -nproc 8 -out ${RESULTS_PATH}/${fprefix}.align; done;
for f in *.bax.h5; do fprefix=`echo "$f" | rev | cut -c 10- | rev`; echo $fprefix; done | sort | uniq | xargs -n 1 -I {} sh -c "{ cat ${RESULTS_PATH}/{}.* > ${RESULTS_PATH}/{}.align; }"

cd ${RESULTS_PATH}
for f in *_p0.align; do fprefix=`echo "$f" | rev | cut -c 7- | rev`; echo ${fprefix}; python ${PBX_SCRIPTS_PATH}/bash5_alignment_to_filter_lists.py ${fprefix}.align ${fprefix}-talreads.txt ${fprefix}-nontalreads.txt; done;
cat *-talreads.txt > talreads.txt
cat *-nontalreads.txt > nontalreads.txt