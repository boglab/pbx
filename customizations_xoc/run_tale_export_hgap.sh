#!/usr/bin/env bash

if [ "$#" -ne 3 ]; then
    echo "Usage: run_tale_export_hgap.sh /path/to/assembly/consensus.fasta /path/to/output/folder strain_settings_name" && exit 1
fi

f=`readlink -f ${1}`
ORIGDIR=`pwd`
OUTPUTDIR=`readlink -f ${2}`
STRAIN=${3}

cd ${OUTPUTDIR}

PBX_PATH=/opt/pbx
PBX_TALE_SEQS_EXPORT=${PBX_PATH}/tale_seqs/exporter/${STRAIN}.fasta
PBX_TALE_SEQS_EXPORT_BOUNDARIES=${PBX_PATH}/tale_seqs/exporter/${STRAIN}.ini
PBX_SCRIPTS_PATH=${PBX_PATH}/scripts

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
    
    assembly=$1
    
    if [ -f blast_results_plus.txt ]; then
        
        python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus blast_results_plus.txt ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} y >> tale_ranges_temp.txt
        
    else
        
        blastn -task blastn -strand plus -query ${assembly} -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | \
        sort -k 1,1 -k 7n,7 | \
        python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py plus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} y >> tale_ranges_temp.txt
        
        blastn -task blastn -strand plus -query ${assembly} -subject ${tale_termini_file} -outfmt "6 qseqid sseqid slen length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | \
        awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if ((index($2, "full") != 0 && abs($3-$4) <= 100) || abs($3-$4) <= 5) print $0}' | \
        cut -d $'\t' -f 1,2,4- | \
        sort -k 1,1 -k 7n,7 > blast_results_plus.txt

    fi
    
    if [ -f blast_results_minus.txt ]; then
        
        python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus blast_results_minus.txt ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} y >> tale_ranges_temp.txt
        
    else
        
        blastn -task blastn -strand minus -query ${assembly} -subject ${tale_termini_file} -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | \
        sort -k 1,1 -k 7n,7 | \
        python2 ${PBX_SCRIPTS_PATH}/blast_to_tale_ranges.py minus stdin ${min_5p_length} ${max_5p_length} ${min_3p_length} ${max_3p_length} y >> tale_ranges_temp.txt
        
        blastn -task blastn -strand minus -query ${assembly} -subject ${tale_termini_file} -outfmt "6 qseqid sseqid slen length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | \
        awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if ((index($2, "full") != 0 && abs($3-$4) <= 100) || abs($3-$4) <= 5) print $0}' | \
        cut -d $'\t' -f 1,2,4- | \
        sort -k 1,1 -k 7n,7 > blast_results_minus.txt
        
    fi

    cat tale_ranges_temp.txt | sort -k 1,1 -k 2n,2 | python2 ${PBX_SCRIPTS_PATH}/name_tal_effectors.py > tale_ranges.txt
    rm tale_ranges_temp.txt
    
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_rvds.py ${PBX_TALE_SEQS_EXPORT_BOUNDARIES} ${assembly} tale_ranges.txt > table.txt
    
    cat table.txt | cut -f 7 | sort > rvd_seqs.txt
    cat rvd_seqs.txt | sort | uniq -c > rvd_seqs_unique.txt
    cat table.txt | column -t -s $'\t' > table_cols.txt
    
    python2 ${PBX_SCRIPTS_PATH}/tale_ranges_export_sequences.py ${assembly} tale_ranges.txt seqs_dna.fasta seqs_prot.fasta
    
}

echo $f > qc_log.txt 2>&1
tblastx -query $f -subject $tale_typical_repeat_file -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen qstart qend sstart send" -max_target_seqs 1000000 | sort -k 1,1 -k 9n,9 | python2 ${PBX_SCRIPTS_PATH}/check_for_boundary_tales.py $max_terminus_length >> qc_log.txt 2>&1

echo $f > export_log.txt 2>&1
export_base $f >> export_log.txt 2>&1

cd ${ORIGDIR}