#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Usage: run_canu.sh /path/to/pbx/results/folder /path/to/raw/reads/folder" && exit 1
fi

ORIGDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUTDIR=`readlink -f ${1}`
OUTPUTDIR_PREFIX=`basename ${OUTPUTDIR}`
RAWREADSDIR=`readlink -f ${2}`

mkdir -p ${OUTPUTDIR}/canu_assembly

bash -e ${SCRIPTDIR}/extract_fastq_reads.sh ${RAWREADSDIR}/Analysis_Results ${OUTPUTDIR}/canu_assembly/${OUTPUTDIR_PREFIX}.fastq

cd ${OUTPUTDIR}/canu_assembly

export PATH=/opt/canu/Linux-amd64/bin:$PATH

canu \
 -p ${OUTPUTDIR_PREFIX} -d canu-output \
 genomeSize=5m errorRate=0.03 minReadLength=4000 minOverlapLength=3000 \
 -pacbio-raw ${OUTPUTDIR_PREFIX}.fastq \
 > canu.out 2> canu.err

cd ${ORIGDIR}