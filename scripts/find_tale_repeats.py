import sys
import math
from Bio import SeqIO

#expects input from:
#tblastx -query Test_4kMinLength_50_16000.fasta -subject bls256_tale_repeat.fasta  -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen qstart qend sstart send frames" -max_target_seqs 1000000 | sort -k 1,1 -k 9n,9 | python2 check_for_boundary_tales.py MAX_TERMINUS_LENGTH

buffer_len = int(sys.argv[1]) + 100

current_read = None
seen_at_start = False
seen_at_end = False

line = sys.stdin.readline()

while line:
    
    line = line.rstrip()
    split_line = line.split('\t')
    
    qname = split_line[0]
    
    qstart = int(split_line[8]) - 1
    qend = int(split_line[9]) - 1
    qlen = int(split_line[2])
    
    slen = int(split_line[3]) - 1
    
    alength = int(split_line[4])
    
    pident = float(split_line[5])
    
    if current_read is None:
        current_read = qname
    
    if qname != current_read:
        seen_at_start = False
        seen_at_end = False
        current_read = qname
    
    if alength > (slen / 3) / 2 and pident >= 60:
        sys.stderr.write("%s\t%d\t%d\t%s" % (qname, qstart, qend, str(pident)) + "\n")
    
    line = sys.stdin.readline()
