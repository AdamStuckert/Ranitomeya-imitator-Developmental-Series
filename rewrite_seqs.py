#!/usr/bin/env python3

######
# purpose: 
# input: fasta file, list of positions to change, list of nucleotides to change to, name of output file
# output: fasta file
# return: fasta file with different nucleotides at positions
######


import sys, os, os.path, argparse, re

try:
    os.path.isfile(sys.argv[1])
    input_fasta = sys.argv[1]
except:
    stop_err("Input file does not exist")

#try:
population = str(sys.argv[2])   ### this needs to be a string list 
#except:
#    stop_err("Population name needed")

#try:
positions_file = open(sys.argv[3], "r")   ### this needs to be a string list 
positions_line = str(positions_file.readline())
positions_str = positions_line.split()
positions = list(positions_str)
#except:
#    stop_err("Positions needed")

#try:
nucleotides_file = open(sys.argv[4], "r")   ### this needs to be a string list
nucleotides_line = str(nucleotides_file.readline())
nucleotides_str = nucleotides_line.split()
nucleotides = list(nucleotides_str)
#except:
#    stop_err("Nucleotides to change to needed")

try:
    os.path.isfile(sys.argv[5]) == False
    output_filename = sys.argv[5]
except:
    stop_err("Output file exists, please verify this is what you want to do and delete it")

output_file = open(sys.argv[5], "w+")


print("Input arguments seem ok, starting filtering")
print("Working with fasta file {}".format(input_fasta))

prev_seq = ""
prev_no_linefeed = ""
prev_header = ""
new_seq = ""

my_file = open(input_fasta, "r")
for line in my_file:
    if line.startswith(">"):
        header = line.rstrip()
        
    else:
        prev_header = header
        seq = line.rstrip()
        prev_seq = prev_seq + seq
        # prev_no_linefeed = prev_seq.replace("\n","")
	
new_seq = str(prev_seq)

def replace_str_index(text,index=0,replacement=''):
    return '%s%s%s'%(text[:index],replacement,text[index+1:])

print("Editing this sequence {}".format(new_seq))
count = 1
pos_change_num = 0
tmp = new_seq
for pos in positions:
    adj_pos = int(pos) - 1
    tmp = replace_str_index(tmp, adj_pos, nucleotides[pos_change_num])
    print("Changing position {} from {} to {}".format(pos, new_seq[adj_pos], nucleotides[pos_change_num]))
    print(adj_pos)
    pos_change_num += 1
#    if count == pos:
#        new_seq.append(nucleotides[pos_change_num])
#        pos_change_num += 1
#        else:
#            new_seq += pos

final_seq = tmp

sequence_to_write = "{0}_{1}\n{2}\n".format(prev_header, population, final_seq)
output_file.write(sequence_to_write)
# prev_seq = ""
# prev_header = ""

my_file.close()
output_file.close()
positions_file.close()
nucleotides_file.close()
