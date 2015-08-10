#!/usr/bin/python
from Bio import SeqIO
import sys
import os



#function to calc percent gc from all seqs in a fasta file 
def get_gc(files_dir, filename):
    
    #track number of each base
    bases = {'A':0, 'C':0, 'G':0, 'T':0}
    
    #open the reads file; for each FASTA sequence, track bases in seq
    for seq_record in SeqIO.parse(files_dir + filename, "fasta"):
        for base in seq_record.seq:
            if base == 'A':
                bases['A'] = bases['A'] + 1
            elif base == 'C':
                bases['C'] = bases['C'] + 1
            elif base == 'G':
                bases['G'] = bases['G'] + 1
            else:
                bases['T'] = bases['T'] + 1
    
    #do the stats
    total_bases = float(sum(bases.values()))
    gc = (bases['G'] + bases['C']) / total_bases * 100
    return gc


#input file in fasta
filenames_dir = sys.argv[1]
filenames = os.listdir(filenames_dir)
filenames.sort()

#summary file
results_file = open(filenames_dir + "summary.txt", "w")
results_file.write("filename \t%GC \n")

#process each file
for filename in filenames:
    gc = get_gc(filenames_dir, filename)
    output = filename + "\t" + str(gc) + "\n"
    results_file.write(output)
