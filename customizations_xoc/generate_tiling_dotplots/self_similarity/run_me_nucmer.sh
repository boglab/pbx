#!/usr/bin/env bash

# http://stackoverflow.com/a/2937433/412582
# If it is set, then an unmatched glob is swept away entirely -- 
# replaced with a set of zero words -- 
# instead of remaining in place as a single word.
shopt -s nullglob

if [ "$#" -ne 2 ]; then
    echo "Usage: run_me_nucmer.sh input_seq.fasta /path/to/folder/for/output/pdfs" && exit 1
fi

ORIGDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INPUTSEQ=`readlink -f ${1}`
OUTPUTDIR=`readlink -f ${2}`

make_plot_global_vs_local() {
    global_name=$1
    global_file=$2
    local_file=$3
    local_file_prefix=`basename "${local_file}" | rev | cut -c 7- | rev`;
    /opt/MUMmer3.23/nucmer --maxmatch --nosimplify -p ${global_name}_global_vs_${local_file_prefix} ${local_file} ${global_file}
    /opt/MUMmer3.23/mummerplot -t postscript -p ${global_name}_global_vs_${local_file_prefix} ${global_name}_global_vs_${local_file_prefix}.delta;
}

make_plot_local_vs_local() {
    local_file=$1
    local_file_prefix=`basename "${local_file}" | rev | cut -c 7- | rev`;
    /opt/MUMmer3.23/nucmer --maxmatch --nosimplify -p ${local_file_prefix}_local_vs_local ${local_file} ${local_file}
    /opt/MUMmer3.23/mummerplot -t postscript -p ${local_file_prefix}_local_vs_local ${local_file_prefix}_local_vs_local.delta;
}

export -f make_plot_global_vs_local
export -f make_plot_local_vs_local

cd ${SCRIPTDIR}

python2 split_reads_tile_big.py ${INPUTSEQ} ../small_sequences ../big_sequences

for f in ../big_sequences/*.fasta; do
    
    fprefix=`basename "$f" | rev | cut -c 7- | rev`
    mkdir -p ${fprefix}/sequences
    python2 tile_genome.py ${f} ${fprefix} ${fprefix}/sequences
    
    cd ${fprefix}
    
    # http://stackoverflow.com/questions/11003418/calling-functions-with-xargs-within-a-bash-script
    # _ is a placeholder for bash to use as argv[0]
    # "$@" expands the arguments provided to the bash script into a quoted version of each one
    find sequences/ -name "*.fasta" | xargs -n 1 -P 8 -I{} bash -c 'make_plot_local_vs_local "$@"' _ {}
    ls *_local_vs_local.ps | sort -t _ -k 2n,2 | xargs gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ${OUTPUTDIR}/${fprefix}_local_vs_local.pdf
    
    find sequences/ -name "*.fasta" | xargs -n 1 -P 8 -I{} bash -c 'make_plot_global_vs_local "$@"' _  ${fprefix} ../${f} {}
    ls ${fprefix}_global_vs_*.ps | sort -t _ -k 5n,5 | xargs gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ${OUTPUTDIR}/${fprefix}_global_vs_local.pdf
    
    cd ../
    
    rm -rf ${fprefix}
    
done;

for f in ../small_sequences/*.fasta; do
    
    fprefix=`basename "$f" | rev | cut -c 7- | rev`
    
    mkdir ${fprefix}
    
    cd ${fprefix}
    
    make_plot_local_vs_local ../${f}
    gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ${OUTPUTDIR}/${fprefix}_local_vs_local.pdf *.ps
    
    cd ../
    
    rm -rf ${fprefix}
    
done;

rm ../big_sequences/*.fasta
rm ../small_sequences/*.fasta

cd ${ORIGDIR}