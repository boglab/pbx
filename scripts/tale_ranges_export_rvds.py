import sys
import re
from Bio import SeqIO
import ConfigParser

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

all_rvd_seqs = []

line_no = 1

with open(sys.argv[3], "r") as input_file:
    for line in input_file:
        
        split_line = line.split('\t')
        seqid = split_line[0]
        tale_start = int(split_line[1])
        tale_end = int(split_line[2])
        tale_len = tale_end - tale_start
        
        #upper makes a copy of the record
        tale_record = seqs[seqid].upper()
        tale_record.seq = tale_record.seq[tale_start : tale_end]
        
        if sys.argv[4] == "minus":
            tale_record.seq = tale_record.seq.reverse_complement()
        
        #truncate the sequence at end of last full codon
        tale_record.seq = tale_record.seq[:tale_len - tale_len % 3]
        
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
        
        all_rvd_seqs.append(" ".join(rvd_list))
        
        if len(sys.argv) > 5:
            print("\t".join([seqid, str(tale_start), str(tale_end), sys.argv[4], str(tale_len)]) + "\t" + " ".join(rvd_list))
        
        line_no += 1

if len(sys.argv) <= 5:
    for rvd_seq in all_rvd_seqs:
        print(rvd_seq)