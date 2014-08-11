from Bio import SeqIO
import sys
for entry in SeqIO.parse(sys.argv[1], "fasta"):
    SeqIO.write([entry], sys.argv[2] + "/" + entry.id.replace("/", "_") + ".fasta", "fasta")