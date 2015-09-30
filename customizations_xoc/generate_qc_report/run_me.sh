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

find_likely_assemblies() {

    assemblytype="${1}"
    
    grep -Ev '(\#|\[|\!$)' unique_tale_seqs.txt | sed 's/^[ ]\+[0-9]\+[ ]\+//g' > unique_tale_seqs_likely.txt
    likelycount=`wc -l unique_tale_seqs_likely.txt | cut -d ' ' -f 1 | tr -d '\n'`
    
    echo "There are ${likelycount} likely TALE RVD sequences."
    
    grepquery=''
    while read line; do line=${line//\*/\\\*}; line=${line//\!/\\\!}; line=${line//\#/\\\#}; line=${line//[]-\>/\\\[\\]->}; grepquery+="|^${line}$"; done < unique_tale_seqs_likely.txt
    grepquery=${grepquery:1}
    
    likelycountmatches=$(for f in tale_seqs/*; do test $(grep -E "${grepquery}" ${f} | sort | uniq | wc -l | cut -d ' ' -f 1 | tr -d '\n') -ge "${likelycount}" && echo ${f}; done || true)
    likelycountmatchcount=`echo -n "${likelycountmatches}" | wc -l | cut -d ' ' -f 1 | tr -d '\n'`
    
    if [ "${likelycountmatchcount}" -eq "0" ]; then
    
        echo "No TALE assemblies contain all likely sequences."
        
        while [ "${likelycount}" -gt "1" ] && [ "${likelycountmatchcount}" -eq "0" ]; do
            likelycount=$[$likelycount-1]
            likelycountmatches=$(for f in tale_seqs/*; do test $(grep -E "${grepquery}" ${f} | sort | uniq | wc -l | cut -d ' ' -f 1 | tr -d '\n') -ge "${likelycount}" && echo ${f}; done || true)
            likelycountmatchcount=`echo -n "${likelycountmatches}" | wc -l | cut -d ' ' -f 1 | tr -d '\n'`
        done;
        
        if [ "${likelycountmatchcount}" -ne "0" ]; then
            echo "Max likely sequences in any assembly is ${likelycount}. Generating metrics for this number instead."
        fi;
        
    fi;
    
    if [ "${likelycountmatchcount}" -ne "0" ]; then
    
        echo "Assemblies with max likely sequences..."
        
        echo "... sorted by by total TALE sequence count:"
        echo "${likelycountmatches}" | xargs -n 1 wc -l | sort -k 1n,1
        
        echo "... sorted by by unique TALE sequence count:"
        for f in tale_seqs/*; do test $(grep -E "${grepquery}" ${f} | sort | uniq | wc -l | cut -d ' ' -f 1 | tr -d '\n') -ge "${likelycount}" && echo ${f}; done | xargs -n 1 -I {} sh -c '{ funiqseqcount=`cat {} | sort | uniq | wc -l`; echo "${funiqseqcount} {}"; }' | sort -k 1n,1
        
        echo "... sorted by by contig count:"
        
        if [ "${assemblytype}" = "resequencing" ]; then
            for f in tale_seqs/*; do test $(grep -E "${grepquery}" ${f} | sort | uniq | wc -l | cut -d ' ' -f 1 | tr -d '\n') -ge "${likelycount}" && echo ${f}; done | xargs -n 1 -I {} sh -c '{ fprefix=`basename "{}" | cut -d '.' -f 1`; fseqcount=`grep -c "^>" resequencing_output_${fprefix}/data/consensus.fasta`; echo "${fseqcount} {}"; }' | sort -k 1n,1
        else
            for f in tale_seqs/*; do test $(grep -E "${grepquery}" ${f} | sort | uniq | wc -l | cut -d ' ' -f 1 | tr -d '\n') -ge "${likelycount}" && echo ${f}; done | xargs -n 1 -I {} sh -c '{ fprefix=`basename "{}" | cut -d '.' -f 1`; fseqcount=`grep -c "^>" ${fprefix/__/\/}.fa`; echo "${fseqcount} {}"; }' | sort -k 1n,1
        fi;
        
    fi;

}




cd ${OUTPUTDIR}

echo -e "\nQC Report for ${OUTPUTDIR_PREFIX}"

echo -e "\nHGAP2 Assembly Contigs:"
python2 /opt/pbx/customizations_xoc/generate_qc_report/asm_contigs_size.py ${OUTPUTDIR}/hgap2_assembly/data/polished_assembly.fasta

echo -e "\nHGAP3 Assembly Contigs:"
python2 /opt/pbx/customizations_xoc/generate_qc_report/asm_contigs_size.py ${OUTPUTDIR}/hgap3_assembly/data/polished_assembly.fasta

source /opt/smrtanalysis/current/etc/setup.sh

echo -e "\nHGAP2 Assembly CA Overlapper Connected Components:"
python2 /opt/pbx/customizations_xoc/run_overlapper_tests/count_connected_components.py ${OUTPUTDIR}/hgap2_assembly/data/ca_best_edges.gml

echo -e "\nHGAP3 Assembly CA Overlapper Connected Components:"
python2 /opt/pbx/customizations_xoc/run_overlapper_tests/count_connected_components.py ${OUTPUTDIR}/hgap3_assembly/data/ca_best_edges.gml


echo -e "\nCA Overlapper Tests:"

for f in ${OUTPUTDIR}/ca_overlapper_tests/gephi_graphs/*.gml; do
    fprefix=`basename ${f}`
    echo ${fprefix}
    python2 /opt/pbx/customizations_xoc/run_overlapper_tests/count_connected_components.py ${f}
done


echo -e "\nLikely good TALE assemblies:"


cd ${OUTPUTDIR}/pbx_16000/resequencing

echo -e "\nPBX 16000 resequencing:"

find_likely_assemblies resequencing

cd ${OUTPUTDIR}/pbx_16000/combine_resequenced_tals

echo -e "\nPBX 16000 combined:"

find_likely_assemblies combined




cd ${OUTPUTDIR}/pbx_12000/resequencing

echo -e "\nPBX 12000 resequencing:"

find_likely_assemblies resequencing

cd ${OUTPUTDIR}/pbx_12000/combine_resequenced_tals

echo -e "\nPBX 12000 combined:"

find_likely_assemblies combined



cd ${ORIGDIR}
