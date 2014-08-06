from Bio import SeqIO
import sys

output_filepath_no_ext = sys.argv[2].split('.')[0]

SeqIO.convert(sys.argv[1], "fastq", sys.argv[2], "fasta")

if len(sys.argv) > 3:
    SeqIO.convert(sys.argv[1], "fastq", output_filepath_no_ext + ".qual", "qual")