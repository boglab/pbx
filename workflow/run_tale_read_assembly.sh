#!/bin/bash

# Requires
# PBX_MUMMER_PATH
# PBX_SCRIPTS_PATH
# PBX_WORKFLOW_PATH
# PBX_PROTOCOL_TEMPLATES_PATH
# PBX_RESULTS_PATH_PREASSEMBLY
# PBX_RESULTS_PATH_ASSEMBLY

RESULTS_PATH=${PBX_RESULTS_PATH_ASSEMBLY}

mkdir -p ${RESULTS_PATH}

python2 ${PBX_SCRIPTS_PATH}/convert_fastq_to_fasta.py ${PBX_RESULTS_PATH_PREASSEMBLY}/data/corrected_50_4000.fastq ${RESULTS_PATH}/Test_4kMinLength_50_4000.fasta y

cd ${RESULTS_PATH}

mkdir banks blast_contigs blast_singletons contigs singletons

parallel --env PBX_PROTOCOL_TEMPLATES_PATH --env PBX_MUMMER_PATH < ${PBX_WORKFLOW_PATH}/run_minimo_nucmer_preassembled.sh