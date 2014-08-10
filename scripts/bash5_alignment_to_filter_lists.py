# Adapted from https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP-Whitelisting-Tutorial

import sys

tale_like_alignment = sys.argv[1]
num_holes = 163482 # int(sys.argv[2])
tale_like_output_file = sys.argv[2]
non_tale_like_output_file = sys.argv[3]

tale_like_set = set()
movie_prefix = None

with open(tale_like_alignment, "r") as alignment_file:
    
    for line in alignment_file:
        
        if line.startswith('qName'):
            continue
        
        cols = line.split()
        
        alignment_id_parts = cols[0].split('/')
        
        if movie_prefix is None:
            movie_prefix = alignment_id_parts[0]
        
        hole_num = alignment_id_parts[1]
        
        tale_like_set.add(hole_num)

with open(tale_like_output_file, "a") as tale_like_output:
    with open(non_tale_like_output_file, "a") as non_tale_like_output:
        
        for i in range(num_holes):
            
            i = str(i)
            
            if i in tale_like_set:
                tale_like_output.write('/'.join([movie_prefix, i]) + '\n')
            else:
                non_tale_like_output.write('/'.join([movie_prefix, i]) + '\n')
