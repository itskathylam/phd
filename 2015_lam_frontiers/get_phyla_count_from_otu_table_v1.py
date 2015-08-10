import sys
import os


#get to working dir and set file name
#os.chdir("/home/kathy/Dropbox/PhD/Data_Analysis/2015_Frontiers/8_otu_bias_analysis/count_phlya_abundance")
otu_filename = sys.argv[1]

#prep outfile
phyla_filename = os.path.splitext(otu_filename)[0] + "_phyla_percent.txt"
phyla_file = open(phyla_filename, "w")

#get otu table
otu_file = open(otu_filename, "r")

#discard first header line
otu_file.readline()

#start dict to keep phyla counts
cosmid = {}
bulk = {}

#process each line, adding to both dicts
for line in otu_file:
    line = line.split(",")
    bulk_count = int(line[1])
    cosmid_count = int(line[2])
    phylum = line[4]
    
    #check if phylum in either dict and add accordingly
    if phylum in cosmid:
        cosmid[phylum] = cosmid[phylum] + cosmid_count
        bulk[phylum] = bulk[phylum] + bulk_count
    else:
        cosmid[phylum] = cosmid_count
        bulk[phylum] = bulk_count

#given a dictionary of phyla counts, return dict of phyla fractions 
def get_phyla_fractions(phyla_dict):
    
    #get total member count 
    total = 0
    for phylum in phyla_dict:
        total = total + phyla_dict[phylum]
    total = float(total)
    
    #make new dict of fractions
    new_dict = {}
    for phylum in phyla_dict:
        new_dict[phylum] = phyla_dict[phylum]/total
        
    return new_dict
    
cosmid_fraction = get_phyla_fractions(cosmid)
bulk_fraction = get_phyla_fractions(bulk)

#write phyla fractions to new file
for item in cosmid_fraction:
    phyla_file.write(item)
    phyla_file.write("\t")
    phyla_file.write(str(format(cosmid_fraction[item], '.9f')))
    phyla_file.write("\t")
    phyla_file.write(str(format(bulk_fraction[item], '.9f')))
    phyla_file.write("\n")

phyla_file.close()

    
    
    
