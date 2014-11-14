from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

genome = SeqIO.read(sys.argv[1], "fasta")
genome_length = len(genome)

start = 0
end = 50000

while True:
    
    local_end = end
    
    if end > genome_length:
        local_end = genome_length
    
    region_id = "%s_%d_%d" % (sys.argv[2], start, local_end)
    record = SeqRecord(genome.seq[start:local_end], id=region_id, description="")
    SeqIO.write([record], sys.argv[3] + "/" + region_id+".fasta", "fasta")
    
    if end > genome_length:
        break
    
    start += 25000
    end += 25000
