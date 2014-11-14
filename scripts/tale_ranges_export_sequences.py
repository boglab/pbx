import sys
from Bio import SeqIO

ranges = {}

seqs = dict([(entry.id, entry) for entry in SeqIO.parse(sys.argv[1], "fasta")])

output_file = sys.stdout

if len(sys.argv) > 3:
    output_file = open(sys.argv[3], "w")

output_translated = False

if len(sys.argv) > 4:
    output_translated = True
    translated_output_file = open(sys.argv[4], "w")

with open(sys.argv[2], "r") as ranges_file:
    
    for line in ranges_file:
        
        split_line = line.rstrip().split("\t")
        seq_id = split_line[0]
        tal_id = seq_id + "_" + split_line[1]
        
        if seq_id not in seqs:
            continue
        
        start = int(split_line[2])
        end = int(split_line[3])
        
        if split_line[6] == "minus":
            substring = seqs[seq_id][start:end]
            rev_comp_seq = substring.seq.reverse_complement()
            substring.seq = rev_comp_seq
        else:
            substring = seqs[seq_id][start:end]
        
        substring.id = tal_id
        substring.description = ""
        
        SeqIO.write([substring], output_file, "fasta")
        
        if output_translated:
            translated_seq = substring.seq.translate()
            substring.seq = translated_seq
            SeqIO.write([substring], translated_output_file, "fasta")

if len(sys.argv) > 3:
    output_file.close()

if len(sys.argv) > 4:
    translated_output_file.close()