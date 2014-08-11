#!/bin/bash

# Requires
# PBX_MUMMER_PATH
# PBX_SCRIPTS_PATH
# PBX_WORKFLOW_PATH
# PBX_RESULTS_PATH_RESEQUENCING
# PBX_RESULTS_PATH_RESEQUENCED_CONTIG_ASSEMBLY
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

export_resequenced_tale_seqs() {
    fprefix=`echo $1 | cut -c 21-`
    blastn -task blastn -strand plus -query $1/data/consensus.fasta -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix}_tale_ranges_plus.txt
    blastn -task blastn -strand minus -query $1/data/consensus.fasta -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix}_tale_ranges_minus.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1/data/consensus.fasta tale_ranges/${fprefix}_tale_ranges_plus.txt plus >> tale_seqs/${fprefix}_temp.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1/data/consensus.fasta tale_ranges/${fprefix}_tale_ranges_minus.txt minus >> tale_seqs/${fprefix}_temp.txt
    cat tale_seqs/${fprefix}_temp.txt | sort > tale_seqs/${fprefix}.txt
    rm tale_seqs/${fprefix}_temp.txt
}

export_resequenced_tale_seqs_tables() {
    fprefix=`echo $1 | cut -c 21-`
    blastn -task blastn -strand plus -query $1/data/consensus.fasta -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix}_tale_ranges_plus.txt
    blastn -task blastn -strand minus -query $1/data/consensus.fasta -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix}_tale_ranges_minus.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1/data/consensus.fasta tale_ranges/${fprefix}_tale_ranges_plus.txt plus y >> tale_seqs_tables/${fprefix}_temp.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1/data/consensus.fasta tale_ranges/${fprefix}_tale_ranges_minus.txt minus y >> tale_seqs_tables/${fprefix}_temp.txt
    cat tale_seqs_tables/${fprefix}_temp.txt | sort -k 1n,1 -k 2n,2 > tale_seqs_tables/${fprefix}.txt
    rm tale_seqs_tables/${fprefix}_temp.txt
    cat tale_seqs_tables/${fprefix}.txt | column -t -s '	' > tale_seqs_tables/${fprefix}_cols.txt
}

export_combined_tale_seqs() {
    fprefix=`echo $1 | rev | cut -c 4- | rev`
    fprefix2=${fprefix/\//__}
    blastn -task blastn -strand plus -query $1 -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix2}_tale_ranges_plus.txt
    blastn -task blastn -strand minus -query $1 -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix2}_tale_ranges_minus.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1 tale_ranges/${fprefix2}_tale_ranges_plus.txt plus >> tale_seqs/${fprefix2}_temp.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1 tale_ranges/${fprefix2}_tale_ranges_minus.txt minus >> tale_seqs/${fprefix2}_temp.txt
    cat tale_seqs/${fprefix2}_temp.txt | sort > tale_seqs/${fprefix2}.txt
    rm tale_seqs/${fprefix2}_temp.txt
}

export_combined_tale_seqs_tables() {
    fprefix=`echo $1 | rev | cut -c 4- | rev`
    fprefix2=${fprefix/\//__}
    blastn -task blastn -strand plus -query $1 -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix2}_tale_ranges_plus.txt
    blastn -task blastn -strand minus -query $1 -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} > tale_ranges/${fprefix2}_tale_ranges_minus.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1 tale_ranges/${fprefix2}_tale_ranges_plus.txt plus y >> tale_seqs_tables/${fprefix2}_temp.txt
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} $1 tale_ranges/${fprefix2}_tale_ranges_minus.txt minus y >> tale_seqs_tables/${fprefix2}_temp.txt
    cat tale_seqs_tables/${fprefix2}_temp.txt | sort -k 1n,1 -k 2n,2 > tale_seqs_tables/${fprefix2}.txt
    rm tale_seqs_tables/${fprefix2}_temp.txt
    cat tale_seqs_tables/${fprefix2}.txt | column -t -s '	' > tale_seqs_tables/${fprefix2}_cols.txt
}

make_mummer_plots() {
    fprefix=`echo $1 | cut -c 21-`;
    
    mkdir -p contig_dotplots/${fprefix}/sequences
    python2 ${PBX_SCRIPTS_PATH}/split_reads.py $1/data/consensus.fasta contig_dotplots/${fprefix}/sequences
    
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


cd ${PBX_RESULTS_PATH_RESEQUENCING}

# export tales from resequenced contigs

echo "Exporting TALEs from resequenced contigs"

mkdir tale_ranges
mkdir tale_seqs
mkdir tale_seqs_tables
mkdir contig_dotplots

for f in resequencing_output_*; do

    echo $f >> qc_log.txt 2>&1
    tblastx -query $f/data/consensus.fasta -subject $tale_typical_repeat_file -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen qstart qend sstart send" -max_target_seqs 1000000 | sort -k 1,1 -k 9n,9 | python2 ${PBX_SCRIPTS_PATH}/check_for_boundary_tales.py $max_terminus_length >> qc_log.txt 2>&1
    
    echo $f >> export_log.txt 2>&1
    export_resequenced_tale_seqs $f >> export_log.txt 2>&1
    
    echo $f >> export_log_tables.txt 2>&1
    export_resequenced_tale_seqs_tables $f >> export_log_tables.txt 2>&1
    
    if command -v gnuplot >/dev/null 2>&1; then
        make_mummer_plots $f > /dev/null 2>&1
    fi
    
done

cat tale_seqs/* | sort | uniq -c > unique_tale_seqs.txt

cd ${PBX_RESULTS_PATH_RESEQUENCED_CONTIG_ASSEMBLY}

 export tales from combined resequenced contigs

echo "Exporting TALEs from assembled resequenced contigs"

mkdir tale_ranges
mkdir tale_seqs
mkdir tale_seqs_tables

for f in */*-combo.fa; do
    
    echo $f >> qc_log.txt 2>&1
    tblastx -query $f -subject $tale_typical_repeat_file -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen qstart qend sstart send" -max_target_seqs 1000000 | sort -k 1,1 -k 9n,9 | python2 ${PBX_SCRIPTS_PATH}/check_for_boundary_tales.py $max_terminus_length >> qc_log.txt 2>&1
    
    echo $f >> export_log.txt 2>&1
    export_combined_tale_seqs $f >> export_log.txt 2>&1
    
    echo $f >> export_log_tables.txt 2>&1
    export_combined_tale_seqs_tables $f >> export_log_tables.txt 2>&1
    
done

cat tale_seqs/* | sort | uniq -c > unique_tale_seqs.txt
