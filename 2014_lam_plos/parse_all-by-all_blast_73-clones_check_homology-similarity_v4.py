#!/usr/bin/python

from Bio.Blast import NCBIXML
'''
#From command line, execute all-by-all blastn to generate results.xml:
#blastn -query contigs-5.fa -subject contigs-5.fa -evalue .001 -out results.xml -outfmt 5
'''

from interval import Interval, IntervalSet
'''
An Interval is composed of the lower bound, a closed lower bound ^M
    flag, an upper bound, and a closed upper bound flag.  The attributes^M
    are called lower_bound, lower_closed, upper_bound, and upper_closed,^M
    respectively.  For an infinite interval, the bound is set to inf or ^M
    -inf.  IntervalSets are composed of zero to many Intervals.
    
#become familiar with interval usage    
r1 = IntervalSet([Interval(1, 1000)])
r2 = IntervalSet([Interval(30, 50)])
r3 = IntervalSet([Interval(1200, 1300)])
print r1 - r2
print r1 + r2
x = r1 + r3
print x
for interval in x:
    print interval.lower_bound
    print interval.upper_bound
'''

file = open("results.xml")     
blast_records = NCBIXML.parse(file)

##accumulate distance between contig pairs in dictionary
distance = {}

##for each queried sequence
for blast_record in blast_records:
    #print "\n" + blast_record.query
    #print str(blast_record.query_letters)
    
    ##for each subject sequence
    for alignment in blast_record.alignments:
	
        ##accumulate hsp intervals for each subject sequence, by iterating through each hsp
        hsp_interval_list = []
        for hsp in alignment.hsps:

	    ##if alignment was on subject complement, subtract alignment length from start to get interval
	    if hsp.frame == (1,-1):
		hsp_interval = IntervalSet([Interval(hsp.sbjct_start, hsp.sbjct_start - hsp.align_length)])
		hsp_interval_list.append(hsp_interval)
	    
	    ##otherwise, alignment was on subject given strand, add alignment length to start to get interval
	    else:
		hsp_interval = IntervalSet([Interval(hsp.sbjct_start, hsp.sbjct_start + hsp.align_length)])
		hsp_interval_list.append(hsp_interval)

	
	##use interval addition to remove overlapping regions over hsps
	new_intervalset = IntervalSet()
	for interval in hsp_interval_list:        
		new_intervalset = new_intervalset + interval

	##calculate length of the subject sequence that was involved in the alignment = [aligned length]
	range_list =[]
	for interval in new_intervalset:
	    start = interval.lower_bound
	    end = interval.upper_bound
	    for i in range(start, end):
		range_list.append(i)
	    
	##check which of query/subject is shorter; then divide the [aligned length] by length of the shorter one
	##note: blast_record.query_letters = query length; alignment.length = subject length
	##keep track of the fraction and query/subject names for putting in dict
	fraction = 0
	if blast_record.query_letters <= alignment.length:
	    fraction = float(len(range_list))/blast_record.query_letters
	else:
	    fraction = float(len(range_list))/alignment.length
	    
	##save the fraction (distance), which represents the homology between the query and subject
	##put the names into a list to sort; this overwrites duplicate key-value pairs in the dictionary
	name_pair = [str(blast_record.query), str(alignment.hit_def)]
	name_pair = sorted(name_pair)
	new_name_pair = ":".join(name_pair)
	distance[str(new_name_pair)] = fraction    
	    
##write distances to file
out = open("out.txt", "w")
for item in distance:
    #print item + "\t\t\t" + str(distance[item])
    names = item.split(":")
    row = names[0] + "," + names[1] + "," + str(distance[item]) + "\n"
    out.write(row)

    

