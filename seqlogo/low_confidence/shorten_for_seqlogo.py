'''
Basically I only really want a window of ~50nt to display for the weblogo [at most!]. I'm shortening sequence in fasta-format files as a result.
'''

infile = open("rna_subtilis_low_confidence.txt", "r")
outfile = open("rna_subtilis_low_confidence_for_seqlogo.txt", "w")

windowsize = 50
oneitherside = windowsize/2 #so I have at most 25nts on either side of the cleavage site that I keep 
distancefromcleavage = 100 #there are 100nts between the cleavage site and end of line

for line in infile: 
	if ">" in line: #if it's a fasta header
		outfile.write(line)
	else: 
		line = line[:-1] #get rid of \n at the end
		cleavagesite = len(line)-distancefromcleavage-1 #the -1 takes into account 0-indexing

		if cleavagesite < oneitherside: #i.e. if the sequence is too short upstream of the cleavage site for you to get all 25nts
			ontheleft = 0 #start the window at the beginning of the string
		else:
			ontheleft = cleavagesite-oneitherside
	
		ontheright = cleavagesite+oneitherside

		shortenedseq = line[ontheleft:ontheright]
		outfile.write(shortenedseq+"\n")
		

outfile.close()
infile.close()
