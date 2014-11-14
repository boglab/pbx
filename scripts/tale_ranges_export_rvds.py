import sys
import re
from Bio import SeqIO
import ConfigParser
from collections import defaultdict

config = ConfigParser.RawConfigParser()
config.read(sys.argv[1])

boundary_pattern = re.compile(config.get('General', 'boundary_pattern'))
rvd_start_pattern = re.compile(config.get('General', 'rvd_start_pattern'))
rvd_end_pattern = re.compile(config.get('General', 'rvd_end_pattern'))
typical_repeat_length = config.getint('General', 'typical_repeat_length')
half_repeat_length = config.getint('General', 'half_repeat_length')

seqs = {}

for entry in SeqIO.parse(sys.argv[2], "fasta"):
    seqs[entry.id] = entry

input_table = []

with open(sys.argv[3], "r") as input_file:
    
    line = input_file.readline()
    
    while line:
        
        trimmed_line = line.rstrip().split("\t")
        
        if len(trimmed_line) != 0:
            input_table.append(trimmed_line)
        
        line = input_file.readline()


for row in input_table:
    
    seqid = row[0]
    tale_name = row[1]
    tale_start = int(row[2])
    tale_end = int(row[3])
    repeats_start = int(row[4])
    repeats_end = int(row[5])
    strand = row[6]
    
    tale_len = tale_end - tale_start
    
    #upper makes a copy of the record
    tale_record = seqs[seqid].upper()
    tale_record.seq = tale_record.seq[repeats_start : repeats_end]
    
    if strand == "minus":
        tale_record.seq = tale_record.seq.reverse_complement()
    
    i = 0
    
    #truncate the sequence at end of last full codon
    #test_seq = tale_record.seq[i:]
    #test_seq = test_seq[:len(test_seq) - len(test_seq) % 3]
    tale_seq = str(tale_record.seq.translate())
    
    repeat_list = []
    
    repeat_search_match = boundary_pattern.search(tale_seq)
    
    while repeat_search_match is not None:
        
        repeat_list.append(tale_seq[:repeat_search_match.start() + 2])
        
        tale_seq = tale_seq[repeat_search_match.start() + 2:]
        
        repeat_search_match = boundary_pattern.search(tale_seq)
    
    repeat_list.append(tale_seq)
    
    rvd_list = []
    
    repeat_no = 1
    
    for repeat in repeat_list:
        
        warning_suffix = ""
        
        rvd_start_match = rvd_start_pattern.search(repeat)
        rvd_end_match = rvd_end_pattern.search(repeat)
        
        if len(repeat) not in [half_repeat_length - 1, half_repeat_length, typical_repeat_length - 1, typical_repeat_length]:
            if len(repeat) > (typical_repeat_length + 8) and (repeat_no != len(repeat_list) or (repeat_no == len(repeat_list) and (rvd_start_match is None or rvd_end_match is None))):
                # probable frameshift
                rvd_list.append("[]->")
                repeat_no += 1
                break
            elif len(repeat) > 100:
                # probable frameshift post-RVD
                warning_suffix = " []->"
            else:
                # atypical repeat
                warning_suffix = "!"
        
        if rvd_start_match is None or rvd_end_match is None:
            rvd_list.append("#")
            repeat_no += 1
            break
        
        rvd = repeat[rvd_start_match.end() : rvd_end_match.start()]
        
        # missing 13th residue
        if len(rvd) == 1:
            rvd = rvd + "*"
        
        rvd += warning_suffix
        
        rvd_list.append(rvd)
        repeat_no += 1
        
        if warning_suffix == " []->":
            break
    
    print("\t".join([seqid, tale_name, str(tale_start), str(tale_end), strand, str(tale_len), (" ".join(rvd_list))]))