'''
Takes as input different structures for a given site, e.g. cwlo_correct_structures.txt and cwlo_consensus_for_varna.txt, and calculates the hamming distance from each non-consensus structure to the consensus
'''
'''
Gets the hamming distance between two dot-bracket structures--first distance between all pairs of structures in the *structures.txt file used to make the consensus
Then it also gets the hamming distance between each structure in *structures.txt and the consensus
'''

import itertools
import numpy as np
import matplotlib.pyplot as plt

structures = open("atpi_correct_structures.txt", "r")
consensus = open("atpi_correct_consensus_for_varna.txt", "r")

structuredict = {} #key: type of fold [MFE etc.], value: dot-bracket structure for that fold

#########################
#populate structuredict
#########################

for line in structures: 
	if "#" in line: 
		header = line[1:-1] #so "#MFE\n" becomes "MFE"
	else: 
		structuredict[header] = line[:-1] # :-1 gets rid of the \n at the end

for line in consensus: 
	if "." in line: #i.e. if it's the structure
		consensusstructure = line[:-1]

############################################
#distances from each structure to consensus
############################################

for structure in structuredict.keys(): 
	
	distances = {} #hamming distances from structure to every other structure--NOTE: does not include consensus 
	for key in structuredict.keys():
		if key != structure:  
			distances[key] = 0

	hamming = 0	
	
	#calculate hamming distance--iterate position-by-position and see at how many positions the non-consensus structure differs from consensus

	fold = structuredict[structure]

	#make sure the fold and consensus structure are of the same length--since some of the structures have a couple extra "."'s 
	#at the end that accounts for the length difference there's no problem in doing this 	

	if len(fold)>len(consensusstructure): 
		fold = fold[:len(consensusstructure)]

	elif len(consensusstructure)>len(fold): 
		consensusstructure = consensusstructure[:len(fold)]

	for i in range(len(fold)): 
		if fold[i] != consensusstructure[i]: 
			hamming += 1

	consensusdist = hamming #distance to consensus

	#print "The hamming distance between the "+structure+" fold and the consensus structure is "+str(hamming)+"."


	#get hamming distance between structure and other structures in structuredict
	for other in structuredict.keys(): 
		if other != structure: 
			secondfold = structuredict[other]

			otherhamming = 0

			if len(fold)>len(secondfold):
				fold = fold[:len(secondfold)]

			elif len(secondfold)>len(fold): 
				secondfold = secondfold[:len(fold)]

			for i in range(len(fold)):
				if fold[i] != secondfold[i]: 
					otherhamming += 1

			distances[other] = otherhamming

	###plot: x axis, folds, y axis, hamming distances
	#I'm coloring the distance to consensus in a darker color

	xticks = distances.keys()
	xticks.append("consensus")
	xticks = np.array(xticks)

	colors = []
	for tick in xticks: 
		if tick == "consensus": 
			colors.append("darkred")
		else: 	
			colors.append("salmon")
	colors = np.array(colors)

	xaxis = np.array(range(len(xticks)))
	yaxis = distances.values()
	yaxis.append(consensusdist)
	yaxis = np.array(yaxis)

	plt.xticks(xaxis, xticks)
	plt.bar(xaxis, yaxis, color=colors, edgecolor=colors)
	plt.title("Hamming distance between the "+ structure+ " structure and other folds")
	plt.xlabel("Other folds")
	plt.ylabel("Hamming distances")
	plt.show()
	
	
consensus.close()
structures.close()
