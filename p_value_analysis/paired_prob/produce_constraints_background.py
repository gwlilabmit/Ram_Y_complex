'''
This takes a fasta-format input file
For each fasta header and corresponding sequence produces a separate file, e.g. >atpI becomes atpi_constrained.txt
The output file consists of blocks of the following: 
>header, e.g. >atpI
corresponding sequence
RNAfold constraint [e.g. "...x."] where x refers to an unpaired base

The first block of three lines has the first base marked as unpaired, the second has the second base, etc. 
The idea is that then we can pipe that output file directly into RNAfold using the -C and -p0 options and grab the kTln(Z) values produce when individually constraining
each base in the structure to be unpaired
From there, if we also have kTln(Z) for the Boltzmann ensemble without constraints, we can find the probability that each base is unpaired 
'''

import numpy as np

#use a fasta-format input file 
infile = open("random_sequences.txt", "r")

for line in infile: 
	if ">" in line:
		gene = line[1:-1] #The line will look something like ">atpI\n" and this just gets us "atpI" so we can use that to name the output file
	else: 
		outfile = open(gene+"_constrained_random.txt", "w") 
		sequence = line[:-1] #get rid of the "\n" at the end of the line so now we just have the atpI sequence
		#initialize an array that's the same length as the sequence for the purpose of producing constraints
		#each 
		constraints = np.empty(len(sequence), dtype=object) #dtype=object just lets us fill this array with "."
		constraints.fill(".")
		for i in range(len(sequence)):
			#this is the loop in which we'll write that block of three lines to the output file 
			outfile.write(">"+gene+"\n") #write ">atpI\n"
			outfile.write(line) #write [corresponding sequence]\n
			constraints[np.where(constraints=='x')]='.'	#clear any x's in constraints from the last iteration 
			constraints[i] = 'x'
			outfile.write("".join(constraints)+"\n")
		outfile.close()

infile.close()

