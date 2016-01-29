#!/usr/bin/env bash

if [ -f /opt/pbx/run_pbx.sh ]; then
    mv /opt/pbx/run_pbx.sh /opt/pbx/run_pbx_bak.sh
fi

if [ -f /opt/pbx/run_pbx_pre_smrtenv.sh ]; then
    mv /opt/pbx/run_pbx_pre_smrtenv.sh /opt/pbx/run_pbx_pre_smrtenv_bak.sh
fi

if [ -f /opt/pbx/run_pbx_smrtenv.sh ]; then
    mv /opt/pbx/run_pbx_smrtenv.sh /opt/pbx/run_pbx_smrtenv_bak.sh
fi

if [ -f /opt/pbx/run_pbx_post_smrtenv.sh ]; then
    mv /opt/pbx/run_pbx_post_smrtenv.sh /opt/pbx/run_pbx_post_smrtenv_bak.sh
fi

#rm -f /opt/pbx/run_pbx.sh
#rm -rf /opt/pbx/results/*

PATTERNSTRAIN=$(awk -F "=" '/tale_seqs_file_export_boundaries/ {print $2}' /opt/pbx/customizations_xoc/config_template.ini | tr -d ' ' | rev | cut -c 5- | rev)

cd /opt/pbx/results/

CUTOFFS=(16000 12000)

for rawfull in /opt/raw_long_*; do
    
    fprefix=`echo "${rawfull}" | cut -c 15-`
    f="/opt/pbx/results/${fprefix}"
    
    find  ${rawfull} -iname ".DS_Store" -o -iname ".*h5" -o -iname ".*xml" -exec rm {} \;
    
    mkdir -p ${f}/logs/
    mkdir -p ${f}/logs/workflow_scripts
    mkdir -p ${f}/raw_tale_read_dotplots
    
    mkdir -p ${f}/canu_assembly
    mkdir -p ${f}/canu_assembly/dotplots
    mkdir -p ${f}/canu_assembly/exported
    
    mkdir -p ${f}/canu_assembly/resequencing
    find ${rawfull}/ -name "*.bax.h5" | sort > ${f}/canu_assembly/resequencing/input.fofn
    cp /opt/pbx/customizations_xoc/RS_Resequencing_HGAP_Xo.1.xml ${f}/canu_assembly/resequencing/settings.xml
    sed -i "s/HGAP_resequencing/ref_${fprefix/-/_}_Canu/g" ${f}/canu_assembly/resequencing/settings.xml
    
    mkdir -p ${f}/hgap2_assembly
    mkdir -p ${f}/hgap2_assembly/dotplots
    mkdir -p ${f}/hgap2_assembly/exported
    find ${rawfull}/ -name "*.bax.h5" | sort > ${f}/hgap2_assembly/input.fofn
    cp /opt/pbx/customizations_xoc/RS_HGAP2_Xo.1.xml ${f}/hgap2_assembly/settings.xml
    
    mkdir -p ${f}/hgap3_assembly
    mkdir -p ${f}/hgap3_assembly/dotplots
    mkdir -p ${f}/hgap3_assembly/exported
    find ${rawfull}/ -name "*.bax.h5" | sort > ${f}/hgap3_assembly/input.fofn
    cp /opt/pbx/customizations_xoc/RS_HGAP3_Xo.1.xml ${f}/hgap3_assembly/settings.xml
    
    #mkdir -p ${f}/hgap3_resequencing
    #find ${rawfull}/ -name "*.bax.h5" | sort > ${f}/hgap3_resequencing/input.fofn
    #cp /opt/pbx/customizations_xoc/RS_Resequencing_HGAP_Xo.1.xml ${f}/hgap3_resequencing/settings.xml
    #sed -i "s/HGAP_resequencing/ref_${fprefix/-/_}_HGAP3/g" ${f}/hgap3_resequencing/settings.xml
    
    for cutoff in "${CUTOFFS[@]}"; do
        pbxresults=${f}/pbx_${cutoff}
        mkdir -p ${pbxresults}
        cp /opt/pbx/customizations_xoc/config_template.ini ${pbxresults}/config.ini
        sed -i "s/__STRAIN__/${fprefix}/g" ${pbxresults}/config.ini
        sed -i "s/__SEED_CUTOFF__/${cutoff}/g" ${pbxresults}/config.ini
        echo "python2 /opt/pbx/pbx.py ${pbxresults}/config.ini > ${f}/logs/pbx_${cutoff}.txt 2>&1" >> /opt/pbx/run_pbx.sh
    done
    
done

echo "/opt/pbx/run_pbx_pre_smrtenv.sh" >> /opt/pbx/run_pbx.sh
echo "/opt/pbx/run_pbx_smrtenv.sh" >> /opt/pbx/run_pbx.sh
echo "/opt/pbx/run_pbx_post_smrtenv.sh" >> /opt/pbx/run_pbx.sh

echo "source /opt/smrtanalysis/current/etc/setup.sh" >> /opt/pbx/run_pbx_smrtenv.sh

for rawfull in /opt/raw_long_*; do
    
    fprefix=`echo "${rawfull}" | cut -c 15-`
    f="/opt/pbx/results/${fprefix}"
    
    echo "bash -e /opt/pbx/customizations_xoc/run_canu/run_me.sh ${f} ${fprefix} ${rawfull} > ${f}/logs/canu_assembly.txt 2>&1" >> /opt/pbx/run_pbx_pre_smrtenv.sh
    echo "referenceUploader  --samIdx='samtools faidx' --saw='sawriter -blt 8 -welter' -c -n ref_${fprefix/-/_}_Canu -f ${f}/canu_assembly/canu-output/${fprefix}.consensus.fasta" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "cp -R /opt/smrtanalysis/current/common/references/ref_${fprefix/-/_}_Canu/ ${f}/canu_assembly/resequencing/" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "fofnToSmrtpipeInput.py ${f}/canu_assembly/resequencing/input.fofn --jobname='Canu_Resequencing_Job' > ${f}/canu_assembly/resequencing/input.xml" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "smrtpipe.py --params=${f}/canu_assembly/resequencing/settings.xml --output=${f}/canu_assembly/resequencing/ xml:${f}/canu_assembly/resequencing/input.xml > ${f}/logs/canu_resequencing.txt 2>&1" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "gunzip -c ${f}/canu_assembly/resequencing/data/consensus.fasta.gz > ${f}/canu_assembly/resequencing/data/consensus.fasta" >> /opt/pbx/run_pbx_smrtenv.sh
    
    echo "fofnToSmrtpipeInput.py ${f}/hgap2_assembly/input.fofn --jobname='HGAP2_Assembly_Job' > ${f}/hgap2_assembly/input.xml" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "smrtpipe.py --params=${f}/hgap2_assembly/settings.xml --output=${f}/hgap2_assembly/ xml:${f}/hgap2_assembly/input.xml > ${f}/logs/hgap2_assembly.txt 2>&1" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "gunzip -c ${f}/hgap2_assembly/data/polished_assembly.fasta.gz > ${f}/hgap2_assembly/data/polished_assembly.fasta" >> /opt/pbx/run_pbx_smrtenv.sh
    
    echo "fofnToSmrtpipeInput.py ${f}/hgap3_assembly/input.fofn --jobname='HGAP3_Assembly_Job' > ${f}/hgap3_assembly/input.xml" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "smrtpipe.py --params=${f}/hgap3_assembly/settings.xml --output=${f}/hgap3_assembly/ xml:${f}/hgap3_assembly/input.xml > ${f}/logs/hgap3_assembly.txt 2>&1" >> /opt/pbx/run_pbx_smrtenv.sh
    echo "gunzip -c ${f}/hgap3_assembly/data/polished_assembly.fasta.gz > ${f}/hgap3_assembly/data/polished_assembly.fasta" >> /opt/pbx/run_pbx_smrtenv.sh

    #echo "referenceUploader  --samIdx='samtools faidx' --saw='sawriter -blt 8 -welter' -c -n ref_${fprefix/-/_}_HGAP2 -f ${f}/hgap2_assembly/data/polished_assembly.fasta" >> /opt/pbx/run_pbx_smrtenv.sh
    #echo "cp -R /opt/smrtanalysis/current/common/references/ref_${fprefix/-/_}_HGAP2/ ${f}/hgap2_resequencing/"
    #echo "fofnToSmrtpipeInput.py ${f}/hgap2_resequencing/input.fofn --jobname='HGAP2_Resequencing_Job' > ${f}/hgap2_resequencing/input.xml" >> /opt/pbx/run_pbx_smrtenv.sh
    #echo "smrtpipe.py --params=${f}/hgap2_resequencing/settings.xml --output=${f}/hgap2_resequencing/ xml:${f}/hgap2_resequencing/input.xml > ${f}/logs/hgap2_resequencing.txt 2>&1" >> /opt/pbx/run_pbx_smrtenv.sh

    #echo "referenceUploader  --samIdx='samtools faidx' --saw='sawriter -blt 8 -welter' -c -n ref_${fprefix/-/_}_HGAP3 -f ${f}/hgap3_assembly/data/polished_assembly.fasta" >> /opt/pbx/run_pbx_smrtenv.sh
    #echo "cp -R /opt/smrtanalysis/current/common/references/ref_${fprefix/-/_}_HGAP3/ ${f}/hgap3_resequencing/"
    #echo "fofnToSmrtpipeInput.py ${f}/hgap3_resequencing/input.fofn --jobname='HGAP3_Resequencing_Job' > ${f}/hgap3_resequencing/input.xml" >> /opt/pbx/run_pbx_smrtenv.sh
    #echo "smrtpipe.py --params=${f}/hgap3_resequencing/settings.xml --output=${f}/hgap3_resequencing/ xml:${f}/hgap3_resequencing/input.xml > ${f}/logs/hgap3_resequencing.txt 2>&1" >> /opt/pbx/run_pbx_smrtenv.sh
    
    echo "bash -e /opt/pbx/customizations_xoc/run_tale_export_hgap.sh ${f}/hgap2_assembly/data/polished_assembly.fasta ${f}/hgap2_assembly/exported ${PATTERNSTRAIN} > ${f}/logs/hgap2_export.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    echo "bash -e /opt/pbx/customizations_xoc/generate_tiling_dotplots/self_similarity/run_me_nucmer.sh ${f}/hgap2_assembly/data/polished_assembly.fasta ${f}/hgap2_assembly/dotplots > ${f}/logs/hgap2_dotplots.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    echo "python2 /opt/pbx/customizations_xoc/run_overlapper_tests/CA_best_edge_to_GML.py ${f}/hgap2_assembly/data/celera-assembler.gkpStore ${f}/hgap2_assembly/data/celera-assembler.tigStore ${f}/hgap2_assembly/data/celera-assembler.gkpStore.fastqUIDmap ${f}/hgap2_assembly/data/4-unitigger/best.edges ${f}/hgap2_assembly/data/ca_best_edges.gml" >> /opt/pbx/run_pbx_smrtenv.sh
    
    echo "bash -e /opt/pbx/customizations_xoc/run_tale_export_hgap.sh ${f}/hgap3_assembly/data/polished_assembly.fasta ${f}/hgap3_assembly/exported ${PATTERNSTRAIN} > ${f}/logs/hgap3_export.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    echo "bash -e /opt/pbx/customizations_xoc/generate_tiling_dotplots/self_similarity/run_me_nucmer.sh ${f}/hgap3_assembly/data/polished_assembly.fasta ${f}/hgap3_assembly/dotplots > ${f}/logs/hgap3_dotplots.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    echo "python2 /opt/pbx/customizations_xoc/run_overlapper_tests/CA_best_edge_to_GML.py ${f}/hgap3_assembly/data/celera-assembler.gkpStore ${f}/hgap3_assembly/data/celera-assembler.tigStore ${f}/hgap3_assembly/data/celera-assembler.gkpStore.fastqUIDmap ${f}/hgap3_assembly/data/4-unitigger/best.edges ${f}/hgap3_assembly/data/ca_best_edges.gml" >> /opt/pbx/run_pbx_smrtenv.sh
    
    echo "bash -e /opt/pbx/customizations_xoc/run_tale_export_hgap.sh ${f}/canu_assembly/resequencing/data/consensus.fasta ${f}/canu_assembly/exported ${PATTERNSTRAIN} > ${f}/logs/canu_export.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    echo "bash -e /opt/pbx/customizations_xoc/generate_tiling_dotplots/self_similarity/run_me_nucmer.sh ${f}/canu_assembly/resequencing/data/consensus.fasta ${f}/canu_assembly/dotplots > ${f}/logs/canu_dotplots.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    
    echo "bash -e /opt/pbx/customizations_xoc/generate_raw_tale_read_dotplots/run_me.sh ${f} > ${f}/logs/raw_tale_read_dotplots.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    
    echo "bash -e /opt/pbx/customizations_xoc/run_overlapper_tests/run_me.sh ${f} > ${f}/logs/overlapper_tests.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    
    echo "bash -e /opt/pbx/customizations_xoc/generate_qc_report/run_me.sh ${f} > ${f}/qc_report.txt 2>&1" >> /opt/pbx/run_pbx_post_smrtenv.sh
    
done

chmod +x /opt/pbx/run_pbx.sh
chmod +x /opt/pbx/run_pbx_pre_smrtenv.sh
chmod +x /opt/pbx/run_pbx_smrtenv.sh
chmod +x /opt/pbx/run_pbx_post_smrtenv.sh

for rawfull in /opt/raw_long_*; do
    
    fprefix=`echo "${rawfull}" | cut -c 15-`
    f="/opt/pbx/results/${fprefix}"
    
    cp /opt/pbx/run_pbx*.sh ${f}/logs/workflow_scripts

done;

echo "Done writing workflow scripts"
echo "To run, start screen and then do: bash -e /opt/pbx/run_pbx.sh"
