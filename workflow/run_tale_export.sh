#!/bin/bash

# Requires
# PBX_MUMMER_PATH
# PBX_SCRIPTS_PATH
# PBX_WORKFLOW_PATH
# PBX_RESEQUENCING_RESULTS_PATH
# PBX_RESEQUENCED_CONTIG_ASSEMBLY_RESULTS_PATH
# PBX_TALE_SEQS_EXPORT
# PBX_TALE_SEQS_EXPORT_BOUNDARIES

tale_typical_repeat_file=$(mktemp)
python2 ${PBX_SCRIPTS_PATH}/fastx_subset.py ${PBX_TALE_SEQS_EXPORT} typical_repeat > ${tale_typical_repeat_file}

tale_termini_file=$(mktemp)
python2 ${PBX_SCRIPTS_PATH}/fastx_subset.py ${PBX_TALE_SEQS_EXPORT} typical_repeat inv > ${tale_termini_file}

max_terminus_length=$(python2 ${PBX_SCRIPTS_PATH}/get_extreme_sequence_length.py ${tale_termini_file} max)
min_5p_length=$(python2 ${PBX_SCRIPTS_PATH}/get_extreme_sequence_length.py ${tale_termini_file} min start_full)
max_5p_length=$(python2 ${PBX_SCRIPTS_PATH}/get_extreme_sequence_length.py ${tale_termini_file} max start_full)
min_3p_length=$(python2 ${PBX_SCRIPTS_PATH}/get_extreme_sequence_length.py ${tale_termini_file} min end_full)
max_3p_length=$(python2 ${PBX_SCRIPTS_PATH}/get_extreme_sequence_length.py ${tale_termini_file} max end_full)

export_base() {
    fprefix=$1
    assembly=$2
    
    
    blastn -task blastn -strand plus -query ${assembly} -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | \
    sort -k 1,1 -k 7n,7 | \
    python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} y >> tale_ranges/${fprefix}_tale_ranges_temp.txt
    
    blastn -task blastn -strand minus -query ${assembly} -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | \
    sort -k 1,1 -k 7n,7 | \
    python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} y >> tale_ranges/${fprefix}_tale_ranges_temp.txt
    
    cat tale_ranges/${fprefix}_tale_ranges_temp.txt | sort -k 1,1 -k 2n,2 | python2 ${PBX_SCRIPTS_PATH}/name_tal_effectors.py > tale_ranges/${fprefix}_tale_ranges.txt
    rm tale_ranges/${fprefix}_tale_ranges_temp.txt
    
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} ${assembly} tale_ranges/${fprefix}_tale_ranges.txt > tale_seqs_tables/${fprefix}.txt
    
    cat tale_seqs_tables/${fprefix}.txt | cut -f 7 | sort > tale_seqs/${fprefix}.txt
    cat tale_seqs_tables/${fprefix}.txt | column -t -s $'\t' > tale_seqs_tables/${fprefix}_cols.txt
    
    
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_sequences.py ${assembly} tale_ranges/${fprefix}_tale_ranges.txt tale_seqs_dna/${fprefix}.fasta tale_seqs_prot/${fprefix}.fasta
}

export_resequenced_tale_seqs() {
    fprefix=`echo $1 | cut -c 21-`
    assembly=${1}/data/consensus.fasta
    export_base ${fprefix} ${assembly}
}

export_combined_tale_seqs() {
    fprefix=`echo $1 | rev | cut -c 4- | rev`
    fprefix2=${fprefix/\//__}
    assembly=$1
    export_base ${fprefix2} ${assembly}
}

make_mummer_plots() {
    fprefix=`echo $1 | cut -c 21-`;
    assembly=${1}/data/consensus.fasta
    
    mkdir -p contig_dotplots/${fprefix}/sequences
    python2 ${PBX_SCRIPTS_PATH}/split_reads.py ${assembly} contig_dotplots/${fprefix}/sequences
    
    cd contig_dotplots/${fprefix}
    
    for f in sequences/*.fasta; do
        fileprefix=`echo "${f}" | cut -c 11- | rev | cut -c 7- | rev`;
        ${PBX_MUMMER_PATH}/nucmer --maxmatch --nosimplify -p ${fileprefix} ${f} ${f};
        ${PBX_MUMMER_PATH}/mummerplot -t postscript -p ${fileprefix} ${fileprefix}.delta;
    done
    
    ls *.ps | sort -t . -k 1n,1 | xargs gs -sDEVICE=pdfwrite -dBATCH -dPDFSettings=/Screen -o ../${fprefix}.pdf
    
    cd ../
    rm -rf ${fprefix}
    cd ../
}


cd ${PBX_RESEQUENCING_RESULTS_PATH}

# export tales from resequenced contigs

echo "Exporting TALEs from resequenced contigs"

mkdir tale_ranges
mkdir tale_seqs
mkdir tale_seqs_tables
mkdir tale_seqs_dna
mkdir tale_seqs_prot
mkdir contig_dotplots

for f in resequencing_output_*; do

    echo $f >> qc_log.txt 2>&1
    tblastx -query $f/data/consensus.fasta -subject $tale_typical_repeat_file -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen qstart qend sstart send" -max_target_seqs 1000000 | sort -k 1,1 -k 9n,9 | python2 ${PBX_SCRIPTS_PATH}/check_for_boundary_tales.py $max_terminus_length >> qc_log.txt 2>&1
    
    echo $f >> export_log.txt 2>&1
    export_resequenced_tale_seqs $f >> export_log.txt 2>&1
    
    if command -v gnuplot >/dev/null 2>&1; then
        make_mummer_plots $f > /dev/null 2>&1
    fi
    
done

cat tale_seqs/* | sort | uniq -c > unique_tale_seqs.txt

cd ${PBX_RESEQUENCED_CONTIG_ASSEMBLY_RESULTS_PATH}

# export tales from combined resequenced contigs

echo "Exporting TALEs from assembled resequenced contigs"

mkdir tale_ranges
mkdir tale_seqs
mkdir tale_seqs_tables
mkdir tale_seqs_dna
mkdir tale_seqs_prot
mkdir contig_dotplots

for f in */*-combo.fa; do
    
    echo $f >> qc_log.txt 2>&1
    tblastx -query $f -subject $tale_typical_repeat_file -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen qstart qend sstart send" -max_target_seqs 1000000 | sort -k 1,1 -k 9n,9 | python2 ${PBX_SCRIPTS_PATH}/check_for_boundary_tales.py $max_terminus_length >> qc_log.txt 2>&1
    
    echo $f >> export_log.txt 2>&1
    export_combined_tale_seqs $f >> export_log.txt 2>&1
    
done

cat tale_seqs/* | sort | uniq -c > unique_tale_seqs.txt
