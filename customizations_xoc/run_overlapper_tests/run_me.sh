#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: run_me.sh /path/to/pbx/results/folder" && exit 1
fi

ORIGDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUTDIR=`readlink -f ${1}`

mkdir -p ${OUTPUTDIR}/ca_overlapper_tests
cd ${OUTPUTDIR}/ca_overlapper_tests

source /opt/smrtanalysis/current/etc/setup.sh

export PATH=/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin:$PATH

/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_1kMinLength_0p03utgErrorRate -d overlap_1kMinLength_0p03utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_1kMinLength_0p03utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_1kMinLength_0p03utgErrorRate.out 2> overlap_1kMinLength_0p03utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_1kMinLength_0p04utgErrorRate -d overlap_1kMinLength_0p04utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_1kMinLength_0p04utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_1kMinLength_0p04utgErrorRate.out 2> overlap_1kMinLength_0p04utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_1kMinLength_0p05utgErrorRate -d overlap_1kMinLength_0p05utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_1kMinLength_0p05utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_1kMinLength_0p05utgErrorRate.out 2> overlap_1kMinLength_0p05utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_1kMinLength_0p06utgErrorRate -d overlap_1kMinLength_0p06utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_1kMinLength_0p06utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_1kMinLength_0p06utgErrorRate.out 2> overlap_1kMinLength_0p06utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_2kMinLength_0p03utgErrorRate -d overlap_2kMinLength_0p03utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_2kMinLength_0p03utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_2kMinLength_0p03utgErrorRate.out 2> overlap_2kMinLength_0p03utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_2kMinLength_0p04utgErrorRate -d overlap_2kMinLength_0p04utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_2kMinLength_0p04utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_2kMinLength_0p04utgErrorRate.out 2> overlap_2kMinLength_0p04utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_2kMinLength_0p05utgErrorRate -d overlap_2kMinLength_0p05utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_2kMinLength_0p05utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_2kMinLength_0p05utgErrorRate.out 2> overlap_2kMinLength_0p05utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_2kMinLength_0p06utgErrorRate -d overlap_2kMinLength_0p06utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_2kMinLength_0p06utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_2kMinLength_0p06utgErrorRate.out 2> overlap_2kMinLength_0p06utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_3kMinLength_0p03utgErrorRate -d overlap_3kMinLength_0p03utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_3kMinLength_0p03utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_3kMinLength_0p03utgErrorRate.out 2> overlap_3kMinLength_0p03utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_3kMinLength_0p04utgErrorRate -d overlap_3kMinLength_0p04utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_3kMinLength_0p04utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_3kMinLength_0p04utgErrorRate.out 2> overlap_3kMinLength_0p04utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_3kMinLength_0p05utgErrorRate -d overlap_3kMinLength_0p05utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_3kMinLength_0p05utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_3kMinLength_0p05utgErrorRate.out 2> overlap_3kMinLength_0p05utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_3kMinLength_0p06utgErrorRate -d overlap_3kMinLength_0p06utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_3kMinLength_0p06utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_3kMinLength_0p06utgErrorRate.out 2> overlap_3kMinLength_0p06utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_500MinLength_0p03utgErrorRate -d overlap_500MinLength_0p03utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_500MinLength_0p03utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_500MinLength_0p03utgErrorRate.out 2> overlap_500MinLength_0p03utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_500MinLength_0p04utgErrorRate -d overlap_500MinLength_0p04utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_500MinLength_0p04utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_500MinLength_0p04utgErrorRate.out 2> overlap_500MinLength_0p04utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_500MinLength_0p05utgErrorRate -d overlap_500MinLength_0p05utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_500MinLength_0p05utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_500MinLength_0p05utgErrorRate.out 2> overlap_500MinLength_0p05utgErrorRate.err
/opt/smrtanalysis/current/analysis/bin/wgs-8.1/Linux-amd64/bin/runCA -p overlap_500MinLength_0p06utgErrorRate -d overlap_500MinLength_0p06utgErrorRate -s /opt/pbx/customizations_xoc/run_overlapper_tests/ca_specs/ca_500MinLength_0p06utgErrorRate.spec ${OUTPUTDIR}/hgap3_assembly/data/corrected.frg > overlap_500MinLength_0p06utgErrorRate.out 2> overlap_500MinLength_0p06utgErrorRate.err

mkdir gephi_graphs

for rd in overlap_*; do
    
    if [ -d "$rd" ] && [ -d "$rd/4-unitigger" ]; then
        
        rdprint=`basename $rd`
        
        python2 /opt/pbx/customizations_xoc/run_overlapper_tests/CA_best_edge_to_GML.py ${rdprint}/${rdprint}.gkpStore ${rdprint}/${rdprint}.tigStore ${rdprint}/${rdprint}.gkpStore.fastqUIDmap ${rdprint}/4-unitigger/best.edges gephi_graphs/${rdprint}.gml
        
    fi;

done;

cd ${ORIGDIR}