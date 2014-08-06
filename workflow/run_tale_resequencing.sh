#!/bin/bash

# Requires
# PBX_SMRTANALYSIS_PATH
# PBX_RAW_READS_PATH
# PBX_PROTOCOL_TEMPLATES_PATH
# PBX_RESULTS_PATH_WHITELISTING
# PBX_RESULTS_PATH_ASSEMBLY
# PBX_RESULTS_PATH_RESEQUENCING

source ${PBX_SMRTANALYSIS_PATH}/etc/setup.sh

RESULTS_PATH=${PBX_RESULTS_PATH_RESEQUENCING}
mkdir -p ${RESULTS_PATH}
cd ${RESULTS_PATH}

RESULTS_REF_PATH=${RESULTS_PATH}/resequencing_refs
REF_PATH_REPLACEMENT=$(echo "${RESULTS_REF_PATH}" | sed -e 's/[\/&]/\\&/g')

mkdir ${RESULTS_REF_PATH}

WHITELIST_REPLACEMENT=$(echo "${PBX_RESULTS_PATH_WHITELISTING}" | sed -e 's/[\/&]/\\&/g')

find ${PBX_RAW_READS_PATH} -name "*.bax.h5" | sort > ${RESULTS_PATH}/pbx_resequencing_input_files.fofn

SETTINGS_NAMES=(
    "500_91"
    "500_93"
    "500_95"
    "500_97"
    "1000_91"
    "1000_93"
    "1000_95"
    "1000_97"
    "2000_91"
    "2000_93"
    "2000_95"
    "2000_97"
    "3000_91"
    "3000_93"
    "3000_95"
    "3000_97"
)

n=0
while [[ $n -lt ${#SETTINGS_NAMES[@]} ]]; do
    
    settings_name=${SETTINGS_NAMES[$n]}
    minimo_results_file="${PBX_RESULTS_PATH_ASSEMBLY}/contigs/Test_4kMinLength_50_4000_${settings_name}-contigs.fasta"
    
    if [ -f $minimo_results_file ] && [ "`grep -c '^>' ${minimo_results_file}`" -ge 1 ]; then
        
        mkdir resequencing_output_${settings_name}-contigs
        
        # upload contigs as reference
        referenceUploader -p ${RESULTS_REF_PATH} -c -n ref_Test_4kMinLength_50_4000_${settings_name}_contigs -f ${PBX_RESULTS_PATH_ASSEMBLY}/contigs/Test_4kMinLength_50_4000_${settings_name}-contigs.fasta
        
        # generate smrtpipe input file
        fofnToSmrtpipeInput.py pbx_resequencing_input_files.fofn --jobname="Resequencing_Test_4kMinLength_50_4000_${settings_name}-contigs" > resequencing_input_${settings_name}-contigs.xml
        
        #cp ${PBX_PROTOCOL_TEMPLATES_PATH}/resequencing_contigs/resequencing_settings_${settings_name}-contigs.xml ${RESULTS_PATH}/resequencing_settings_${settings_name}-contigs.xml 
        cp ${PBX_PROTOCOL_TEMPLATES_PATH}/RS_Resequencing_TALs.1.xml ${RESULTS_PATH}/resequencing_settings_${settings_name}-contigs.xml
        sed -i "s/__PBX_RESEQUENCING_SETTINGS__/${settings_name}/g" ${RESULTS_PATH}/resequencing_settings_${settings_name}-contigs.xml
        sed -i "s/__PBX_WHITELISTING_RESULTS__/${WHITELIST_REPLACEMENT}/g" ${RESULTS_PATH}/resequencing_settings_${settings_name}-contigs.xml
        sed -i "s/__PBX_RESEQUENCING_REFERENCES__/${REF_PATH_REPLACEMENT}/g" ${RESULTS_PATH}/resequencing_settings_${settings_name}-contigs.xml
        
        # run resequencing
        smrtpipe.py --params=${RESULTS_PATH}/resequencing_settings_${settings_name}-contigs.xml --output=${RESULTS_PATH}/resequencing_output_${settings_name}-contigs xml:${RESULTS_PATH}/resequencing_input_${settings_name}-contigs.xml
        
    fi

    n=$(($n+1))
    
done

for f in resequencing_output_*; do
    gunzip -c $f/data/consensus.fasta.gz > $f/data/consensus.fasta
    gunzip -c $f/data/consensus.fastq.gz > $f/data/consensus.fastq
done;
