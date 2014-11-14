from Bio import SeqIO
import sys

for seq in SeqIO.parse(sys.argv[1], "fasta"):
    print("%s: %d" % (seq.id, len(seq)))