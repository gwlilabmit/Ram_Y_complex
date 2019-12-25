#Extracts second-column values from .dat files and prints them out, comma-separated, so they can be used as a colormap in VARNA
#It'll do this for all .dat files you have in your directory. If you don't want this feature just comment out everything with read_files in it
#and unindent as needed. 

#I also plot reads for a given sequence transformed both ways for the sake of comparison. 

import re 
import numpy as np
import glob

read_files = glob.glob("*.dat")

i = 1
for datfile in read_files: 
	infile = open(datfile, "r")
	#comment out this regex stuff if your .dat file isn't named "gene.dat"--with my naming convention this extracts the gene name for me
	regex = r"^[a-zA-Z]+"
	matches = re.findall(regex, datfile) #say the filename is atpi.dat. This extracts "atpi"
	name = matches[0]

	values = [] #array of all second-column values, i.e. the values of interest for the colormap

	for line in infile: 
		reads = line.split("\t")[1] #Each line is tab-separated. We want the value in the second column. 
		reads = reads[:-1] #There's a \n at the end of the "reads" value, which counts as a single character. 
		values.append(reads)
	
	values = np.array(values[:]).astype(float)
	
		
	print name+" reads:\n" + ",".join(values.astype(str))+"\n" #i.e. print "atpI reads: \n" followed by the reads 
		
	infile.close()
