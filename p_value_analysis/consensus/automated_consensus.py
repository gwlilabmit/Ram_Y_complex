################
# Basically this will create a "consensus structure" from a set of RNAs using the following approach: 
#	-First shorten each sequence to the length of the shortest one provided
#	-So now the indices for one sequence should be the same as that for any other sequence, i.e. a[20] refers to the same position as b[20] for structures a and b
#		-Shortening the sequences prevents issues like "a[20] doesn't exist because len(a)=19."
#	-If the majority of sequences [define as max(num((., (, ))) at that index] have a "(" at index 23, index 23 of the consensus sequence will also have a "(" there. 
#	-I also convert the sequences into a paired [|] or unpaired [.] format--I get a consensus for this by requiring that >=half of the sequences have a nt paired
#	for it to show up as paired on the consensus. 
#
#	-THIS DOES NOT YET CORRECT FOR UNBALANCED PARENS.
#		-For paren matching I use the standard approach--increment a counter +1 for each left paren in the sequence, -1 for each right paren, say the sequence is 
#		unbalanced if the counter!= 0.		  
################

import numpy
import glob
import re 
import os

#For your input structure file, each new structure should be on a separate line. 

structure_files = glob.glob("*random_structures.txt") #file with all the structures you care about
#We'll also write our final paren-matched consensus sequence to a new fasta-format file that has: 
#>sequence name
#sequence
#consensus_matched

#######################################
#get names for each consensus structure
#######################################

nameregex = r"(^[a-zA-Z]\d+)|(^[a-zA-Z]+)"
bsuregex = r"\d+"

for structurefile in structure_files: 
	nametuple = re.findall(nameregex, structurefile)[0]
	for i in range(len(nametuple)): 
		if len(nametuple[i])> 0: 
			name = nametuple[i]

	if name == "BSU": #for genes of the name BSU_misc_RNA_3 or something--we want to isolate the 3
		bsunumber = re.findall(bsuregex, structurefile)[0]
		name = "BSU_misc_RNA_"+bsunumber

	if name == "Random": #random sequences are named something like Random_sequence_3600_random.txt--we want to isolate the 3600
		randomnumber = re.findall(bsuregex, structurefile)[0]
		name = "Random_sequence_"+randomnumber

	############################################
	#get sequences for each consensus structure
	############################################

	sequencedir = "../paired_prob/"
	sequencefiles = os.path.join(sequencedir, name+"_random_sequence.txt")


	sequence_file = open(sequencefiles, "r")
	outfile = open(name+"_random_consensus.txt", "w")

	#Initialize a list of dot brackets based on the structures present in the input file 

	dot_brackets = [] 

	structures = open(structurefile, 'r')
	for line in structures:
		if "#" not in line and "." in line: #I'm commenting my files with "#." Additionally, each structure needs at least one unpaired base so this prevents  
			dot_brackets.append(line[:-2])		#me from adding nonempty lines or lines with nonsense structures--the [:-2] takes off the \n at the end. 

	structures.close()

	#Set the length of every sequence to that of the shortest sequence
	#Why? Rnafold tends to add extra dots/brackets to the end of the regular MFE structure so I need to chop those off 

	minlength = len(min(dot_brackets, key=len))
	for i in range(len(dot_brackets)): 
		structure = dot_brackets[i]
		if len(structure)>minlength:
			dot_brackets[i] = structure[:minlength]
			
	#Initialize counter arrays for the number of (, ), . seen at each index across all sequences
	#So, if sequences a, b and c all start with "(" I'd say that left_paren[0] += 3

	bracketlength = len(min(dot_brackets, key=len))
	dots = numpy.zeros(bracketlength)
	left_paren = numpy.zeros(bracketlength)
	right_paren = numpy.zeros(bracketlength)

	#Initialize arrays for the paired/unpaired sequences [dot-bar annotation]

	paired = dot_brackets #This is an array where, for each structure in dot_brackets, ( and ) is replaced by |, e.g. "(..)"->"|..|"
	pairedvals = numpy.zeros(bracketlength) #This is a counter array for the number of |'s seen at each index across all sequences
	pairedcutoff = 0.5 #Change this as needed--this is the percent of sequences in which a nt must be paired for me to include it in the consensus. 


	#For each location, count how many times we see a ., ( or ) across all sequences considered. 

	for structure in range(len(dot_brackets)):
		rna = dot_brackets[structure]
		for position in range(len(rna)): 
			if rna[position] == ".": 
				dots[position] += 1

			elif rna[position] == "(": 
				left_paren[position] += 1 
				paired[structure]=paired[structure][:position]+"|"+paired[structure][position+1:] #replace ( with |

			elif rna[position] == ")": 
				right_paren[position] += 1
				paired[structure]=paired[structure][:position]+"|"+paired[structure][position+1:] #replace ) with |


	#Normalize the ., ( and ) counts by the number of structures used to calculate them
	#So if all sequences have a "." at an index it'll show up as a 1
	#Sure, this step isn't strictly necessary but it helps if I'd like to use hard numerical cutoffs as a means of determining when to add an element to the consensus sequence
	#E.g. requiring that 75% of sequences contain __ at a given position for it to be included in the consensus. 
	 
	dots = dots/len(dot_brackets)
	left_paren = left_paren/len(dot_brackets)
	right_paren = right_paren/len(dot_brackets)
		
	#Now that we have our dot-bar structures for paired/unpaired nts, let's count the number of times we see a given nt paired across all sequences 

	for structure in range(len(paired)): 
		rna = paired[structure]
		for position in range(len(rna)): #count up the number of paired indices for each structure
			if rna[position] == "|": 
				pairedvals[position] += 1

	pairedvals = pairedvals/len(dot_brackets) #Normalize by the number of structures considered--in this case this step is useful because I'm going to use a numerical cutoff

	consensus_struct = ""
	consensus_paired = ""

	#Now we're going to create the consensus structures. 
	#For the regular consensus structure: if the majority of sequences have, say, ( at index __, index __ of the consensus = (
	#For the paired consensus structure: if a fraction of sequences >= pairedcutoff has index __ paired, it'll be paired in the consensus, else we set it to . 

	for position in range(len(dots)): 
		dotcount = dots[position]
		leftcount = left_paren[position]
		rightcount = right_paren[position]
		pairedcount = pairedvals[position]

		consensus_feature = max(dotcount, leftcount, rightcount) #set the consensus feature to the feature [. / ( / )] that's present in the majority of the sequences 
		if consensus_feature == dotcount: 
			consensus_struct += "."
		elif consensus_feature == leftcount: 
			consensus_struct += "("
		elif consensus_feature == rightcount:
			consensus_struct += ")"

		if pairedcount>=pairedcutoff: #if the majority of sequences have this sequence as paired
			consensus_paired += "|"
		else: 
			consensus_paired += "."

	##############
	#Paren fixing
	##############
	# I parse through the consensus sequences as lists because I can mutate them that way.
	# As I'm iterating through consensus_matched:  
	# I push the indices of all left parens I see in consensus_matched to indices
	# If I find a ) I pop the stack. This ) is the match of the most recent ( I found. 
	# If there's nothing to pop I have an unmatched ) and I set it to .
	# If I get through the whole list and have a nonempty stack I have too many ('s--I go through and set them all to . 
	# I then convert these lists back into strings  

	indices = [] #stack containing all indices corresponding to unpaired left parens 

	consensus_matched = list(consensus_struct) 

	for i in range(len(consensus_matched)): 
		if consensus_matched[i]=="(": #push the index of each ( to indices  
			indices.append(i)
		elif consensus_matched[i]==")":
			if len(indices)==0: #if we have an unpaired ), i.e. if the indices stack is empty, set that ) to .
				consensus_matched[i]="."
			else: 
				indices.pop() #the ) is a match to the last ( in indices so we no longer consider it unpaired

	if len(indices)>0: #if there are unpaired left parens, i.e. if the stack is nonempty after we've gone through the whole sequences
		for i in range(len(indices)): 
			index = indices[i]
			consensus_matched[index] = "."

	paired_matched = consensus_matched[:] 
	#Now that we have a paren-matched consensus_matched list, we can copy it and put it in dot-bar notation to get a paren-matched paired sequences
	#I'm just setting each ( and ) to |
	for i in range(len(paired_matched)): 
		if paired_matched[i] == ")" or paired_matched[i]=="(": 
			paired_matched[i] = "|"

	#convert the lists we've been working with into strings 
	consensus_matched = "".join(consensus_matched)
	paired_matched = "".join(paired_matched)


	####write to output file#####
	#write >sequence_name and sequence to the output file
	for line in sequence_file: 
		outfile.write(line)

	#write consensus_matched to output file
	outfile.write(consensus_matched)

	outfile.close()
	sequence_file.close()

		
