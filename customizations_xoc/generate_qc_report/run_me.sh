#!/usr/bin/env bash

# http://stackoverflow.com/a/2937433/412582
# If it is set, then an unmatched glob is swept away entirely -- 
# replaced with a set of zero words -- 
# instead of remaining in place as a single word.
shopt -s nullglob

if [ "$#" -ne 1 ]; then
    echo "Usage: generate_qc_report.sh /path/to/pbx/results/folder" && exit 1
fi

ORIGDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUTDIR=`readlink -f ${1}`
OUTPUTDIR_PREFIX=`basename ${OUTPUTDIR}`

cd ${OUTPUTDIR}

echo -e "\nQC Report for ${OUTPUTDIR_PREFIX}" >> qc_report.txt

echo -e "\nHGAP2 Assembly Contigs:" >> qc_report.txt
python2 /opt/pbx/customizations_xoc/generate_qc_report/asm_contigs_size.py ${OUTPUTDIR}/hgap2_assembly/data/polished_assembly.fasta >> qc_report.txt

echo -e "\nHGAP3 Assembly Contigs:" >> qc_report.txt
python2 /opt/pbx/customizations_xoc/generate_qc_report/asm_contigs_size.py ${OUTPUTDIR}/hgap3_assembly/data/polished_assembly.fasta >> qc_report.txt

source /opt/smrtanalysis/current/etc/setup.sh

echo -e "\nHGAP2 Assembly CA Overlapper Connected Components:" >> qc_report.txt
python2 /opt/pbx/customizations_xoc/run_overlapper_tests/count_connected_components.py ${OUTPUTDIR}/hgap2_assembly/data/ca_best_edges.gml >> qc_report.txt

echo -e "\nHGAP3 Assembly CA Overlapper Connected Components:" >> qc_report.txt
python2 /opt/pbx/customizations_xoc/run_overlapper_tests/count_connected_components.py ${OUTPUTDIR}/hgap3_assembly/data/ca_best_edges.gml >> qc_report.txt


echo -e "\nCA Overlapper Tests:" >> qc_report.txt

for f in ca_overlapper_tests/gephi_graphs/*.gml; do
    fprefix=`basename ${f}`
    echo ${fprefix}
    python2 /opt/pbx/customizations_xoc/run_overlapper_tests/count_connected_components.py ${f}
done >> qc_report.txt

cd ${ORIGDIR}