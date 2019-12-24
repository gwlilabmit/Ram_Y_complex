import re
import os
import glob
import random

###########################################################
#sample and fold constrained and unconstrained structures
###########################################################

###open files###

pairdir = "../paired_prob/"

pairedfiles = os.path.join(pairdir, "*_random_paired_folded.txt") #I have to go through this extra step to read files from another directory
paired = glob.glob(pairedfiles) #for constrained structures--using the name of the pairing-constrained structures because there are less of them

###regexes###

nameregex = r"(^[a-zA-Z]\d+)|(^[a-zA-Z]+)"
bsuregex = r"\d+"

###global vars###


#num_samples = 100


#random_sample = random.sample(sequences, num_samples) #randomly sample num_samples files from structures--greatly reduces runtime

#for sequence in random_sample:

for sequence in paired:
	
	seq = sequence.split('/')[2] #just gives you the filename, so ../paired_prob/cwlo_3prime_paired_folded.txt just becomes cwlo_3prime_paired_folded.txt
	nametuple = re.findall(nameregex, seq)[0]
	for i in range(len(nametuple)): 
		if len(nametuple[i])> 0: 
			name = nametuple[i]

	if name == "BSU": #for genes of the name BSU_misc_RNA_3 or something--we want to isolate the 3
		bsunumber = re.findall(bsuregex, seq)[0]
		name = "BSU_misc_RNA_"+bsunumber

	if name == "Random": #random sequences are named something like Random_sequence_3600_random.txt--we want to isolate the 3600
		randomnumber = re.findall(bsuregex, seq)[0]
		name = "Random_sequence_"+randomnumber

	mfedir = "../mfe/"	
	mfe = os.path.join(mfedir, name+"_random.txt")

	centroiddir = "../centroid/"
	centroid = os.path.join(centroiddir, name+"_random.txt")

	meadir = "../mea/"
	mea = os.path.join(meadir, name+"_random.txt")

	###################################################
	#create consensus structures file for consensus.py
	###################################################

	outfile = open(name+"_random_structures.txt", 'w')

	pairedstruct = open(sequence, 'r')
	mfestruct = open(mfe, 'r')
	centroidstruct = open(centroid, 'r')
	meastruct = open(mea, 'r')

	#extract the structures from each file--note that I have to use [0] to get rid of some random free energies etc. following the structures [always separated by a space]

	for line in pairedstruct: 
		if "." in line: 
			pairedline = line.split(' ')[0] 
			outfile.write('#Probability-constrained\n'+pairedline+"\n")

	for line in mfestruct: 
		if "." in line: 
			mfeline = line.split(' ')[0]
			outfile.write('#MFE\n'+mfeline+"\n")

	for line in centroidstruct: 
		if "." in line: 
			centroidline = line.split(' ')[0]
			outfile.write('#Centroid\n'+centroidline+'\n')

	for line in meastruct: 
		if "." in line: 
			mealine = line.split(' ')[0]
			outfile.write('#MEA\n'+mealine+"\n")

	outfile.close()

	

