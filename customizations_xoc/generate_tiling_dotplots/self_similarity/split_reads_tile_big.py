from Bio import SeqIO
import sys
for entry in SeqIO.parse(sys.argv[1], "fasta"):
    fixed_id = entry.id.replace("/", "-").replace("|", "-").replace("_", "-").replace(".",  "-")
    if len(entry.seq) < 100000:
        SeqIO.write([entry], sys.argv[2] + "/" + fixed_id + ".fasta", "fasta")
    else:
        SeqIO.write([entry], sys.argv[3] + "/" + fixed_id + ".fasta", "fasta")
