#!/usr/bin/python

from Bio import SeqIO
import sys
import os
import re



#FUNCTIONS

#look for consensus sequences 1 promoter; return count
def find_one_consensus(sequence, filename):
    
    #compile regex
    p = re.compile(sequence)
    count = 0
    
    #iterate through each fasta sequence
    for seq_record in SeqIO.parse(filename, "fasta"):
        
        #check the sequence
        for match in p.finditer(str(seq_record.seq)):
            count = count + 1
            
        #check the reverse complement
        for match in p.finditer(str(seq_record.reverse_complement().seq)):
            count = count + 1
    
    return count

#look for consensus sequences for 5 promoters; return a string to be printed to file
def find_all_consensus(files_dir, reads_filename):
    
    #file location
    location = files_dir + reads_filename
    
    #rpoD sigma 70
    rpod_count = find_one_consensus("TTGACA.{15,19}TATAAT", location)
    
    #rpoE sigma 24
    rpoe_count = find_one_consensus("GGAACTT.{15,19}TCAAA", location) 
    
    #rpoH sigma 32
    rpoh_count = find_one_consensus("TTG[AT][AT][AT].{13,14}CCCCAT[AT]T", location) 
    
    #rpoN sigma 54
    rpon_count = find_one_consensus("TGGCA.{7}TGC", location) 
    
    #Bacteroides sigma AB
    bacteroides_count = find_one_consensus("TTTG.{19,21}TA.{2}TTTG", location)  
    
    output = filename + "\t" + str(rpod_count) + "\t" + str(rpoe_count) + "\t" + str(rpoh_count) + "\t" + str(rpon_count)+ "\t" + str(bacteroides_count) + "\n"
    return output



#INPUT FILES

filenames_dir = sys.argv[1]
filenames = os.listdir(filenames_dir)
filenames.sort()


#PROCESS ALL FILES

#write summary file of results
summary_file = open(filenames_dir + "summary.txt", "w")
summary_file.write("filename \trpoD reads \trpoE \trpoH \trpoN \tBacteroides \n")

#process files
for filename in filenames:
    
    #get sets of read names that are hits
    output = find_all_consensus(filenames_dir, filename)
    
    #write to summary
    summary_file.write(output)
