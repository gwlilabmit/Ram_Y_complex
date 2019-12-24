import re
import os
import glob
import random

###########################################################
#sample and fold constrained and unconstrained structures
###########################################################

###open files###

sequences = glob.glob("*_random_sequence.txt") #for constrained structures

###regexes###

nameregex = r"(^[a-zA-Z]\d+)|(^[a-zA-Z]+)"
bsuregex = r"\d+"

###global vars###
'''

num_samples = 800


random_sample = random.sample(sequences, num_samples) #randomly sample num_samples files from structures--greatly reduces runtime

for sequence in random_sample:
	 
	nametuple = re.findall(nameregex, sequence)[0]
	for i in range(len(nametuple)): 
		if len(nametuple[i])> 0: 
			name = nametuple[i]

	if name == "BSU": #for genes of the name BSU_misc_RNA_3 or something--we want to isolate the 3
		bsunumber = re.findall(bsuregex, sequence)[0]
		name = "BSU_misc_RNA_"+bsunumber

	if name == "Random": #random sequences are named something like Random_sequence_3600_random.txt--we want to isolate the 3600
		randomnumber = re.findall(bsuregex, sequence)[0]
		name = "Random_sequence_"+randomnumber

	
	#check that you haven't already folded this sample--if so, advance to the next loop iter

	constrained = name+"_random_constrained_structures.txt"
	unconstrained = name+"_random_partition.txt"
	paired = name+"random_paired_folded.txt"

	#to get the constrained structures
	
	if os.path.exists(constrained): 
		print constrained
		continue
	else: 
		constrained = "rnafold --noPS -p0 -C --infile="+name+"_constrained_random.txt --outfile="+name+"_random_constrained_structures.txt" 
	
	#to get the unconstrained structures
	
	if os.path.exists(unconstrained):
		print unconstrained 
		continue
	else: 
		unconstrained = "rnafold --noPS -p0 --infile="+name+"_random_sequence.txt --outfile="+name+"_random_partition.txt"

	os.system(constrained)
	os.system(unconstrained)
'''
##############################################################
#fold pairing-prob structures for all sequences with datfiles
##############################################################	
	
###open files###

structures = glob.glob("*_random_constrained_structures.txt") #for constrained structures

for structure in structures:
	 
	nametuple = re.findall(nameregex, structure)[0]
	for i in range(len(nametuple)): 
		if len(nametuple[i])> 0: 
			name = nametuple[i]

	if name == "BSU": #for genes of the name BSU_misc_RNA_3 or something--we want to isolate the 3
		bsunumber = re.findall(bsuregex, structure)[0]
		name = "BSU_misc_RNA_"+bsunumber

	if name == "Random": #random sequences are named something like Random_sequence_3600_random.txt--we want to isolate the 3600
		randomnumber = re.findall(bsuregex, structure)[0]
		name = "Random_sequence_"+randomnumber

	
	#check that you haven't already folded this sample--if so, advance to the next loop iter

	paired = name+"_random_paired_folded.txt"
	
	if os.path.exists(paired): 
		print paired
		continue
	else: 
		pairing = "rnafold --noPS --shape="+name+"_random_paired_probs.dat --infile="+name+"_random_sequence.txt --outfile="+name+"_random_paired_folded.txt"
	
	os.system(pairing)
 
