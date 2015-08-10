#!/usr/bin/python

from Bio import SeqIO
import sys
import os
import time



#input: directory of files to process; fasta Ec file; fasta vector file
filenames_dir = sys.argv[1]
vector_filename = sys.argv[2]
ec_filename = sys.argv[3]

#get list of filenames into array to process
filenames = os.listdir(filenames_dir)
filenames.sort()
 
#for ec, vector: get the sequence, rev comp of the sequence, in preparation for checking
ec = SeqIO.read(ec_filename, "fasta")
ec_rc = ec.reverse_complement()
vector = SeqIO.read(vector_filename, "fasta")
vector_rc = vector.reverse_complement()


#prep output file
outfile = open(filenames_dir + "results_Ec_or_pJC8.txt", "w")
outfile.write("filename \ttotal \tboth \tEc \tvector \tunaccounted \n")


#process each file 
for filename in filenames:
    
    #check whether each read in the file is from pJC8 or Ec or both; should not be any unaccounted, but track in case
    both_count = 0
    ec_count = 0
    vector_count = 0
    unaccounted = 0
    total = 0
    unaccounted_file = open(filenames_dir + filename + "_unaccounted_reads", "w")
    
    for seq_record in SeqIO.parse(filenames_dir + filename, "fasta"):
        total = total + 1
        
        #if seq in both
        if (seq_record.seq in ec.seq or seq_record.seq in ec_rc.seq):
            ec_count = ec_count + 1
            if (seq_record.seq in vector.seq or seq_record.seq in vector_rc.seq):
                vector_count = vector_count + 1
                both_count = both_count + 1
        
        elif (seq_record.seq in vector.seq or seq_record.seq in vector_rc.seq):
            vector_count = vector_count + 1
        
        #this shouldn't happen
        else:
            unaccounted = unaccounted + 1
            SeqIO.write(seq_record, unaccounted_file, "fasta")


    #write to output file: filename, total num reads, num Ec reads, num pjc8 reads
    output_line = filename + "\t" + str(total) + "\t" + str(both_count) + "\t" + str(ec_count) + "\t" + str(vector_count) + "\t" + str(unaccounted) + "\n"
    outfile.write(output_line)
