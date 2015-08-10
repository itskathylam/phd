#!/usr/bin/python

from Bio import SeqIO
import sys
import os
import time



#FUNCTIONS

#run blat and parse results; return a set of unique read names that are hits to the subject
def run_blat(files_dir, reads_filename, subject_filename):
    
    #run blat in the shell
    results_filename = reads_filename + "_BLAT_" + subject_filename + ".psl"
    os.system("blat " + subject_filename + " " + files_dir + reads_filename + " " + files_dir + results_filename)

    #open results 
    results_file = open(files_dir + results_filename)
    
    #clear the header lines
    for i in range(0,5):
        results_file.readline()

    #track the names of reads that are 100% identical to E. coli (90 base identity)
    match_names = set()
    for line in results_file:
        
        #parse the line
        line = line.split('\t')
        match = line[0]
        mismatch = line[1]
        gaps = line[6]
        query_name = line[9]
        
        #if the match was 100% identical (90 bases), accumulate the name
        if match == '90' and mismatch == '0' and gaps == '0':
            match_names.add(query_name)
    
    #delete psl files
    os.system("rm " + files_dir + "*.psl")
    
    return match_names



#INPUT FILES

filenames_dir = sys.argv[1]
vector_filename = sys.argv[2]
ec_filename = sys.argv[3]

#get list of filenames into array to process
filenames = os.listdir(filenames_dir)
filenames.sort()



#RUN BLAT AND PARSE RESULTS FOR EACH FILE

#write summary file of results
summary_file = open(filenames_dir + "summary.txt", "w")
summary_file.write("filename \ttotal reads \ttotal dirty \tec \tvector \n")

#process files
for filename in filenames:
    
    #get sets of read names that are hits
    ec_hits = run_blat(filenames_dir, filename, ec_filename)
    vector_hits = run_blat(filenames_dir, filename, vector_filename) 
    
    #track for summary file
    total_count = 0
    total_dirty_count = 0
    vector_count = 0
    ec_count = 0
    
    #write clean and dirty reads to  new files; also summary file
    clean_file = open(filenames_dir + filename + "_clean_chked.fa", "w")
    dirty_file = open(filenames_dir + filename + "_dirty_chked.fa", "w")
    
    #open the reads file; for each FASTA sequence read
    for seq_record in SeqIO.parse(filenames_dir + filename, "fasta"):
        total_count = total_count + 1
        
        if (seq_record.id in ec_hits):
            SeqIO.write(seq_record, dirty_file, "fasta")
            ec_hits.remove(seq_record.id) #remove id from set to make following searches faster 
            ec_count = ec_count + 1
            total_dirty_count = total_dirty_count + 1
            
        elif (seq_record.id in vector_hits):
            SeqIO.write(seq_record, dirty_file, "fasta")
            vector_hits.remove(seq_record.id) #remove id from set to make following searches faster 
            vector_count = vector_count + 1
            total_dirty_count = total_dirty_count + 1
        
        #if not in list of read names, it's a clean read
        else:
            
            #write to clean file
            SeqIO.write(seq_record, clean_file, "fasta")
    
    #write to summary
    output = filename + "\t" + str(total_count) + "\t" + str(total_dirty_count) + "\t" + str(ec_count) + "\t" + str(vector_count) + "\n"
    summary_file.write(output)
