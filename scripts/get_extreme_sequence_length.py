import sys
from Bio import SeqIO

if sys.argv[2] == "min":
    if len(sys.argv) > 3:
        print(min([len(entry.seq) for entry in SeqIO.parse(sys.argv[1], "fasta") if entry.id.startswith(sys.argv[3])]))
    else:
        print(min([len(entry.seq) for entry in SeqIO.parse(sys.argv[1], "fasta")]))
else:
    if len(sys.argv) > 3:
        print(max([len(entry.seq) for entry in SeqIO.parse(sys.argv[1], "fasta") if entry.id.startswith(sys.argv[3])]))
    else:
        print(max([len(entry.seq) for entry in SeqIO.parse(sys.argv[1], "fasta")]))
