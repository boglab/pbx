import sys
from Bio import SeqIO
import os.path

invert = (len(sys.argv) > 3)

keep_me = {}

if len(sys.argv) > 2:
    if os.path.isfile(sys.argv[2]):
        with open(sys.argv[2], "r") as input_file:
            
            line = input_file.readline().strip()
            
            while line:
                keep_me[line] = 1
                line = input_file.readline().strip()
    else:
        keep_me[sys.argv[2]] = 1
    
else:
    
    line = sys.stdin.readline().strip()
    
    while line:
        keep_me[line] = 1
        line = sys.stdin.readline().strip()

ext = sys.argv[1][-5:]

kept_reads = []

if invert:
    for entry in SeqIO.parse(sys.argv[1], ext):
        if not (entry.id in keep_me or ('/' in entry.id and entry.id.rsplit('/', 1)[0] in keep_me) or (',' in entry.id and entry.id.split(',', 1)[0] in keep_me)):
            kept_reads.append(entry)
else:
    kept_reads = [entry for entry in SeqIO.parse(sys.argv[1], ext) if entry.id in keep_me or ('/' in entry.id and entry.id.rsplit('/', 1)[0] in keep_me) or (',' in entry.id and entry.id.split(',', 1)[0] in keep_me)]

SeqIO.write(kept_reads, sys.stdout, ext)
