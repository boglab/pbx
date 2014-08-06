from Bio import SeqIO
import sys
import uuid

ext = sys.argv[1][-5:]

records = []

for record in SeqIO.parse(sys.argv[1], ext):
    record.id = str(uuid.uuid1().int>>64)
    records.append(record)

SeqIO.write(records, sys.argv[2] + "." + ext, ext)