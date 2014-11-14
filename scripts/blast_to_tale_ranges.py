import sys
import math

# expects input from one of:
#blastn -task blastn -strand plus -query Test_4kMinLength_50_16000.fasta -subject tale_seqs/exporter/xoc_bls256.fasta -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 blast_to_tale_ranges.py plus
#blastn -task blastn -strand minus -query Test_4kMinLength_50_16000.fasta -subject tale_seqs/exporter/xoc_bls256.fasta -outfmt "6 qseqid sseqid length pident mismatch gapopen qstart qend sstart send" -perc_identity 85 | sort -k 1,1 -k 7n,7 | python2 blast_to_tale_ranges.py minus

results = []
notices = []
errors = []

# params
length_cutoff = 100

is_minus_strand = (sys.argv[1] == "minus")

check_oddities = False

min_start = None
max_start = None
min_end = None
max_end = None

if len(sys.argv) > 3:
    check_oddities = True
    min_start = int(sys.argv[3])
    max_start = int(sys.argv[4])
    min_end = int(sys.argv[5])
    max_end = int(sys.argv[6])

#util

def add_result(current_read, start_pos, end_pos):
    
    # uses global variables
    
    #if print_full_range:
    
    #strand = "minus" if is_minus_strand else "plus"
    #
    #if is_minus_strand and start_pos.start > end_pos.start:
    #    repl_tuple = (current_read, end_pos.start, start_pos.end + 1, strand)
    #else:
    #    repl_tuple = (current_read, start_pos.start, end_pos.end + 1, strand)
    #
    #results.append("%s\t%d\t%d\t%s\n" % repl_tuple)
    
    strand = "minus" if is_minus_strand else "plus"
    
    if is_minus_strand and start_pos.start > end_pos.start:
        repl_tuple = (current_read, end_pos.start, start_pos.end + 1, end_pos.end + 1, start_pos.start, strand)
    else:
        repl_tuple = (current_read, start_pos.start, end_pos.end + 1, start_pos.end + 1, end_pos.start, strand)
    
    results.append("%s\t%d\t%d\t%d\t%d\t%s\n" % repl_tuple)
    
    #else:
    #    if is_minus_strand and start_pos.start > end_pos.start:
    #        repl_tuple = (current_read, end_pos.end + 1, start_pos.start)
    #    else:
    #        repl_tuple = (current_read, start_pos.end + 1, end_pos.start)
    #    results.append("%s\t%d\t%d\n" % repl_tuple)
    
    if check_oddities:
        
        notice_name = "%s-%d-%d" % repl_tuple[:3]
        
        if min_start - (start_pos.end + 1 - start_pos.start) >= 5:
            notices.append("Notice: %s has short 5' end (%d)\n" % (notice_name, start_pos.end + 1 - start_pos.start))
        
        if (start_pos.end + 1 - start_pos.start) - max_start >= 5:
            notices.append("Notice: %s has long 5' end(%d)\n" % (notice_name, start_pos.end + 1 - start_pos.start))
        
        if min_end - (end_pos.end + 1 - end_pos.start) >= 5:
            notices.append("Notice: %s has short 3' end (%d)\n" % (notice_name, end_pos.end + 1 - end_pos.start))
        
        if (end_pos.end + 1 - end_pos.start) - max_end >= 5:
            notices.append("Notice: %s has long 3' end (%d)\n" % (notice_name, end_pos.end + 1 - end_pos.start))

class BLASTRange:
    def __init__(self, start, end):
        self.start = start
        self.end = end

# read blast input

input_lines = []

if len(sys.argv) == 2 or sys.argv[2] == "stdin":
    
    line = sys.stdin.readline()
    
    while line:
        input_lines.append(line)
        line = sys.stdin.readline()
    
else:
    
    with open(sys.argv[2], "r") as input_file:
        
        line = input_file.readline()
        
        while line:
            input_lines.append(line)
            line = input_file.readline()

if is_minus_strand:
    input_lines.reverse()

# ping pong time

seen_full_start = False
start_pos = None
seen_full_end = False
end_pos = None
current_read = None

try:
    
    line_iter = iter(input_lines)
    line = next(line_iter)
    
    while True:
        
        line = line.rstrip()
        split_line = line.split('\t')
        
        qname = split_line[0]
        sname = split_line[1]
        alength = int(split_line[2])
        qstart = int(split_line[6]) - 1
        qend = int(split_line[7]) - 1
        
        if current_read is None:
            current_read = qname
        
        if qname != current_read:
            
            # new read, cleanup
            
            if start_pos is not None and end_pos is not None:
                add_result(current_read, start_pos, end_pos)
            
            seen_full_start = False
            start_pos = None
            seen_full_end = False
            end_pos = None
            current_read = qname
        
        i_ate_the_rest = False
        
        if alength >= length_cutoff:
            
            if sname.startswith("start_full") or sname.startswith("start_short"):
                
                if end_pos is not None and start_pos is None:
                    
                    seen_full_end = False
                    end_pos = None
                    
                elif end_pos is not None:
                    
                    #print out previous range
                    add_result(current_read, start_pos, end_pos)
                    #if is_minus_strand and start_pos.start > end_pos.start:
                    #    print("%s\t%d\t%d" % (current_read, end_pos.end + 1, start_pos.start))
                    #else:
                    #    print("%s\t%d\t%d" % (current_read, start_pos.end + 1, end_pos.start))
                    
                    seen_full_start = False
                    start_pos = None
                    seen_full_end = False
                    end_pos = None
                
                if sname.startswith("start_full"):
                    
                    if start_pos is None or (start_pos is not None and qstart == start_pos.start or qend == start_pos.end):
                        seen_full_start = True
                        start_pos = BLASTRange(qstart, qend)
                    else:
                        #shouldn't happen
                        #raise ValueError("Full start match with non-equal start pos or end pos.\n\nLine:\n%s" % line)
                        errors.append("\nFull start match with non-equal start pos or end pos. Skipping remainder of read matches. \n\nLine:\n%s\n" % line)
                        seen_full_start = False
                        start_pos = None
                        seen_full_end = False
                        end_pos = None
                        temp_line = next(line_iter)
                        while temp_line.split('\t')[0] == current_read:
                            temp_line = next(line_iter)
                        line = temp_line
                        i_ate_the_rest = True
                    
                elif sname.startswith("start_short"):
                    
                    if not seen_full_start or (seen_full_start and (qstart == start_pos.start or qend == start_pos.end) and qend - qstart > start_pos.end - start_pos.start):
                        start_pos = BLASTRange(qstart, qend)
                
            elif sname.startswith("end_full") or sname.startswith("end_short"):
                
                if sname.startswith("end_full"):
                    
                    if end_pos is None or (end_pos is not None and qstart == end_pos.start or qend == end_pos.end):
                        seen_full_end = True
                        end_pos = BLASTRange(qstart, qend)
                    else:
                        #shouldn't happen
                        #raise ValueError("Full end match with non-equal start pos or end pos.\n\nLine:\n%s" % line)
                        errors.append("Full end match with non-equal start pos or end pos. Skipping remainder of read matches. \n\nLine:\n%s\n" % line)
                        seen_full_start = False
                        start_pos = None
                        seen_full_end = False
                        end_pos = None
                        temp_line = next(line_iter)
                        while temp_line.split('\t')[0] == current_read:
                            temp_line = next(line_iter)
                        line = temp_line
                        i_ate_the_rest = True
                    
                elif sname.startswith("end_short"):
                    if not seen_full_end or (seen_full_end and (qstart == end_pos.start or qend == end_pos.end) and qend - qstart > end_pos.end - end_pos.start):
                        end_pos = BLASTRange(qstart, qend)
        
        if not i_ate_the_rest:
            line = next(line_iter)
        
except StopIteration:
    pass

# check if there's one remaining result to add
if start_pos is not None and end_pos is not None:
    add_result(current_read, start_pos, end_pos)

#print results and other info

for result in results:
    sys.stdout.write(result)

if len(notices) > 0:
    sys.stderr.write("Notices:\n")
    for notice in notices:
        sys.stderr.write(notice)

if len(errors) > 0:
    sys.stderr.write("Errors:\n")
    for error in errors:
        sys.stderr.write(error)