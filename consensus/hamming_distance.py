'''
Takes as input different structures for a given site, e.g. cwlo_correct_structures.txt and cwlo_consensus_for_varna.txt, and calculates the hamming distance from each non-consensus structure to the consensus
'''

structures = open("pyrg_structures.txt", "r")
consensus = open("pyrg_consensus_for_varna.txt", "r")

structuredict = {} #key: type of fold [MFE etc.], value: dot-bracket structure for that fold

#populate structuredict

for line in structures: 
	if "#" in line: 
		header = line[1:-1] #so "#MFE\n" becomes "MFE"
	else: 
		structuredict[header] = line[:-1] # :-1 gets rid of the \n at the end

for line in consensus: 
	if "." in line: #i.e. if it's the structure
		consensusstructure = line[:-1]

for structure in structuredict.keys(): 

	hamming = 0	
	
	#calculate hamming distance--iterate position-by-position and see at how many positions the non-consensus structure differs from consensus

	fold = structuredict[structure]

	#make sure the fold and consensus structure are of the same length--since some of the structures have a couple extra "."'s 
	#at the end that accounts for the length difference there's no problem in doing this 	

	if len(fold)>len(consensusstructure): 
		fold = fold[:len(consensusstructure)-1] #-1 needed to account for 0-indexing

	elif len(consensusstructure>len(fold)): 
		consensusstructure = consensusstructure[:len(fold)-1]

	for i in range(len(fold)): 
		if fold[i] != consensusstructure[i]: 
			hamming += 1

	print "The hamming distance between the "+structure+" fold and the consensus structure is "+str(hamming)+"."

consensus.close()
structures.close()
