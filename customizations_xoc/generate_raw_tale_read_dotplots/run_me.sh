#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Usage: run_me.sh /path/to/pbx/results/folder" && exit 1
fi

ORIGDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUTDIR=`readlink -f ${1}`

mkdir -p ${OUTPUTDIR}/raw_tale_read_dotplots

cd ${SCRIPTDIR}

mkdir -p workdir/sequences
cd workdir

python2 ${SCRIPTDIR}/filter_fasta_by_length.py ${OUTPUTDIR}/pbx_12000/preassembly_job/data/filtered_subreads.fasta 20000 > tale_reads_20k.fasta
python2 ${SCRIPTDIR}/split_reads.py tale_reads_20k.fasta sequences
rm tale_reads_20k.fasta
find sequences/ -name "*.fasta" | xargs -n 1 -P 8 ../make_plot.sh
gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ${OUTPUTDIR}/raw_tale_read_dotplots/raw_tale_read_dotplots_20k.pdf *.ps
rm sequences/*

python2 ${SCRIPTDIR}/filter_fasta_by_length.py ${OUTPUTDIR}/pbx_12000/preassembly_job/data/filtered_subreads.fasta 16000 20000 > tale_reads_16k.fasta
python2 ${SCRIPTDIR}/split_reads.py tale_reads_16k.fasta sequences
rm tale_reads_16k.fasta
find sequences/ -name "*.fasta" | xargs -n 1 -P 8 ../make_plot.sh
gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ${OUTPUTDIR}/raw_tale_read_dotplots/raw_tale_read_dotplots_16k.pdf *.ps
rm sequences/*

python2 ${SCRIPTDIR}/filter_fasta_by_length.py ${OUTPUTDIR}/pbx_12000/preassembly_job/data/filtered_subreads.fasta 10000 16000 > tale_reads_10k.fasta
python2 ${SCRIPTDIR}/split_reads.py tale_reads_10k.fasta sequences
rm tale_reads_10k.fasta
find sequences/ -name "*.fasta" | xargs -n 1 -P 8 ../make_plot.sh
gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ${OUTPUTDIR}/raw_tale_read_dotplots/raw_tale_read_dotplots_10k.pdf *.ps

cd ${SCRIPTDIR}

rm -rf workdir

cd ${ORIGDIR}