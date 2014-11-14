from Bio import SeqIO
import sys

lower_cutoff = int(sys.argv[2])

upper_cutoff = 0
if len(sys.argv) > 3:
    upper_cutoff = int(sys.argv[3])

for entry in SeqIO.parse(sys.argv[1], "fasta"):
    if upper_cutoff > 0:
        if len(entry.seq) >= lower_cutoff and len(entry.seq) < upper_cutoff:
            SeqIO.write([entry], sys.stdout, "fasta")
    elif len(entry.seq) >= lower_cutoff:
        SeqIO.write([entry], sys.stdout, "fasta")

