import sys
from collections import defaultdict

input_table = []

line = sys.stdin.readline()

while line:
    
    trimmed_line = line.rstrip().split("\t")
    
    if len(trimmed_line) != 0:
        input_table.append(trimmed_line)
    
    line = sys.stdin.readline()

if len(input_table) == 0:
    sys.stderr.write("TALE RVD sequence exporter encountered an empty file. Bye!")
    exit(1)

region_counts = defaultdict(int)

current_seq = input_table[0][0]
current_end = int(input_table[0][1]) - 1
current_strand = input_table[0][-1]
current_region = 0

for row in input_table:
    
    seqid = row[0]
    tale_start = int(row[1])
    tale_end = int(row[2])
    repeats_start = int(row[3])
    repeats_end = int(row[4])
    strand = row[5]
    
    if not (seqid == current_seq and strand == current_strand and tale_start - current_end <= 10000):
        current_region += 1
    
    region_counts[current_region] += 1
    current_seq = seqid
    current_end = tale_end
    current_strand = strand

current_seq = input_table[0][0]
current_end = int(input_table[0][1]) - 1
current_strand = input_table[0][5]
current_region = 0
current_tale_in_region = 0

for row in input_table:
    
    seqid = row[0]
    tale_start = int(row[1])
    tale_end = int(row[2])
    repeats_start = int(row[3])
    repeats_end = int(row[4])
    strand = row[5]
    
    tale_len = tale_end - tale_start
    
    if not (seqid == current_seq and strand == current_strand and tale_start - current_end <= 10000):
        current_region += 1
        current_tale_in_region = 0
    
    if strand == "minus":
        tal_name_chr = chr(region_counts[current_region] + 97 - current_tale_in_region - 1) if region_counts[current_region] > 1 else ""
    else:
        tal_name_chr = chr(current_tale_in_region + 97) if region_counts[current_region] > 1 else ""
    
    tal_name = "tal%d%s" % ((current_region + 1), tal_name_chr)
    
    print("\t".join([seqid, tal_name] + row[1:]))
    
    current_seq = seqid
    current_end = tale_end
    current_strand = strand
    current_tale_in_region += 1
